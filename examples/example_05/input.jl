using Minerva
using Minerva.Lib
using Minerva.Lib.PlasmaFunctions
using Minerva.Lib.Base.Atomic
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkitStandardLibrary.Blocks
using SciCompDSL
using DataFrames, CSV
using Interpolations


# ============================================================================ #
# Simulation parameters
# ============================================================================ #
# Simulation time interval
tspan = (0.0, 60.0*60*24*3);

# INPUT variables
R::Float64 = 9.073; 		      # major radius [m]
a::Float64 = 2.927; 		      # minor radius [m]
Vplasma::Float64 = 2579.0;    # plasma volume [m3]
k::Float64 = Vplasma / (2*pi*R * pi*a^2)             # plasma shape factor
Pfus::Float64 = 2.0;          # fusion power [GW]
xaz::Float64 = 0.796;        # active zone (normalized radius) [-]
S_FW::Float64 = 1462.0;        # First Wall surface [m2]
T_torus_out::Float64 = 273.15 # tempeature at the outlet of the torous [K]
τp::Float64 = 31.63; 		      # particle confinement time [s] - estimated from PROCESS extended output (nD + nT)*Vplasma/Γ_fuelingRate
τp_α::Float64 = 18.20; 		      # alpha particle confinement time [s]

# Residence times
tau_P1::Float64 = 2*60.0 
tau_PI_Q::Float64 = 20*60.0
tau_PI_He::Float64 = 20*60.0
tau_PR::Float64 = 40*60.0
tau_GasPuff::Float64 = 30.0

# Efficiencies / separation factors
eta_DSS::Float64 = 0.80
eta_PI::Float64 = 0.95
eta_PR::Float64 = 0.835
eta_S2::Float64 = 0.2
eta_fuelingEff::Float64 = 0.25 # Fueling efficiency

# Gas puffing
DTpuff_in::Float64 = 1.0e23 # [s-1] target D+T (50/50) puff particle flow rate
Arpuff_in::Float64 = 1.0e19

# Outgassing
OG_H2::Float64 = 9.0e-6        # H2 OutGassing from FW [Pa m3 s-1 m-2]

# Vessels values
T_vessels::Float64 = 300.15 # [K]
V_V1::Float64 = 0.1  # [m3]
V_V2::Float64 = 0.1  # [m3]
V_V3::Float64 = 0.1  # [m3]
p_V1::Float64 = 0.8e5 # [Pa]
p_V2::Float64 = 0.9e5 # [Pa]
p_V3::Float64 = 0.8e5 # [Pa]


# n, T profiles
# references:
# - Siccinio (2022) - https://doi.org/10.1016/j.fusengdes.2022.113047
# - Siccinio (2020) - https://doi.org/10.1016/j.fusengdes.2020.111603
# The cureves are digitalized from the published plots. An accurate analytical fit is hard.
# I use the Interpolations.jl package to query the profile, doing linear interpolations between data points
# and linear extrapolations outside the given range
data = CSV.read(dirname(@__FILE__)*"/2018_DEMO_profiles.csv", DataFrame; delim = ','); # load csv into datafram
dataDict = Dict(name => collect(skipmissing(data[!, name])) for name in names(data)); # convert to Dict to remove "missing" elements
interp_lin_Ti = linear_interpolation(dataDict["rho_Ti"], dataDict["Ti_keV"], extrapolation_bc = Line()); # linear interp of curve
interp_lin_Te = linear_interpolation(dataDict["rho_Te"], dataDict["Te_keV"], extrapolation_bc = Line()); # linear interp of curve
interp_lin_ne = linear_interpolation(dataDict["rho_ne"], dataDict["ne_1e19m-3"], extrapolation_bc = Line()); # linear interp of curve

# Get T [keV] profile
T_profile(x) = interp_lin_Ti(x)
# Get D, T density (n) averaged over the minor radius a
ne_profile(x) = interp_lin_ne(x)*1e19;
n = get_nAveraged_disk_2018(ne_profile) / 2 # nD, nT
# Define scaling factor. Here we pass from [at] to [mol (monoatomic)]
SF::Float64 = Na; # Scaling Factor [-]

# Get reaction rates in [m6 s-1] (they must be multiplied by n^2)
# mol/s = n * n * rR [mol/m3] * [mol/m3] * [m6/(mol s)]
rR_DD_1, rR_DD_2, rR_DT, rR_DHe3, rR_TT, rR_THe3_1, rR_THe3_2, rR_THe3_3 = get_reactionRates(xaz, k, a, R, T_profile)
rR_DD_1 = rR_DD_1 * SF
rR_DD_2 = rR_DD_2 * SF
rR_DT = rR_DT * SF
rR_DHe3 = rR_DHe3 * SF
rR_TT = rR_TT * SF
rR_THe3_1 = rR_THe3_1 * SF
rR_THe3_2 = rR_THe3_2 * SF
rR_THe3_3 = rR_THe3_3 * SF

# Scale density
n = n/SF
# Compute fusion energy
Ev_to_J = 1/1.602176634e−19
Efus = 17.6e+6 / (Ev_to_J * 1.0e+9) # [GW] from MeV to GW
# Compute density averaging corretino factor for the calculations of the reactions
Power_nAverageCorrectionFactor = get_nAvg_disk_powerCorrectionFactor_2018(ne_profile, T_profile, xaz)
# Reference power value and reference flow rate to achieve the desired power
Power_ref = Power_nAverageCorrectionFactor * (n*n*rR_DT) * Efus * SF
GammaTarget = n*Vplasma/τp

# H2 outgassing rate [Pa m3 s-1 m-2] --> H release rate [s-1]
OG_H = (OG_H2 * S_FW) * Na*2 / (SF * 8.314 * T_torus_out) # [s-1]
OG_H_gain = 1.0

fXe = 3.685e-4




# ============================================================================ #
#tag MAIN MODEL
# ============================================================================ #
@mtkmodel TorusTest begin
  @components begin
    sourceD  = Storage2(xH=0.0, xD=1.0, xT=0.0, xHe=0.0, xI=0.0, xXe=0.0)
    sourceT  = Storage2(xH=0.0, xD=0.0, xT=1.0, xHe=0.0, xI=0.0, xXe=0.0)
    source_MFC = MixerTwoWay()

    valve1 = FixedActuatedValve()
    valve2 = FixedActuatedValve()

    BC_PI_b1 = FluidOpenBC()
    BC_torusMixer = Cap()
    BC_V2 = CapOutlet()

    ramp_FM  = Blocks.Ramp(height=2*GammaTarget/eta_fuelingEff, duration=0.0)

    V2 = Vessel3I2O(T=T_vessels, V=V_V2, p_op=p_V2)
    V3 = Vessel1I1O(T=T_vessels, V=V_V3, p_op=p_V3)

    # , k=2000.0, Ti=350.0, Td=50.0)
    PIDsource1 = Blocks.LimPID(u_min=0.0, u_max=40*(DTpuff_in/SF + 2*GammaTarget), k=10.0, Td=0.2)
    PIDsource2 = Blocks.LimPID(u_min=0.0, u_max=40*(DTpuff_in/SF + 2*GammaTarget), k=10.0, Td=0.2)
    PID1 = Blocks.LimPID(u_min=0.0, k=10.0, Td=0.2)

    PI = modPelletInjector(τ_Q=tau_PI_Q, τ_He=tau_PI_He)
    FuelingEfficiency = Splitter(eta = eta_fuelingEff)
    myTorus = Plasma(τp=τp, τp_α=τp_α, Vₚₗₐₛₘₐ=Vplasma, n=n,
      rR_DD_1=rR_DD_1, rR_DD_2=rR_DD_2, rR_DT=rR_DT, rR_TT=rR_TT,
      rR_THe3_1=rR_THe3_1, rR_THe3_2=rR_THe3_2, rR_THe3_3=rR_THe3_3,
      wc=Power_nAverageCorrectionFactor) # plasma model
    myFW = Blocks.Step(height = OG_H, start_time = 0.0) # outgassing
    FWGain = Blocks.Gain(k = OG_H_gain) # outgassing
    myTorusMixer = TorusMixer_withFuelEff() # outlet of the torus
    P1 = Pump(τ=tau_P1)
    DSS = MFP(eta=eta_DSS) # Metal Foil Pump
    BC_DSS = FluidOpenBC()

    S2 = Splitter(eta=eta_S2)
    myPR = IRPR(τ=tau_PR, eta=eta_PR)
    BC_PR = FluidOpenBC()

    # Argon
    sourceXe = Storage2(xH=0.0, xD=0.0, xT=0.0, xHe=0.0, xI=0.0, xXe=1.0)
    ramp_Xe  = Blocks.Ramp(height=fXe*2*GammaTarget/eta_fuelingEff, duration=0.0)
    mixerXe = MixerTwoWay()
  end

  @equations begin
    connect(sourceXe.b, mixerXe.a2)
    connect(ramp_Xe.output, sourceXe.a_G)
    connect(sourceD.b, source_MFC.a1)
    connect(sourceT.b, source_MFC.a2)
    connect(source_MFC.b, V2.a1)
    connect(V2.b1, valve1.a)
    connect(valve1.b, V3.a)
    connect(V3.b, valve2.a)
    connect(valve2.b, mixerXe.a1)
    connect(mixerXe.b, PI.a)
    connect(PI.b1, BC_PI_b1.a)
    connect(PI.b2, FuelingEfficiency.a)
    connect(FuelingEfficiency.b1, myTorus.a)
    connect(FuelingEfficiency.b2, myTorusMixer.a_nf)
    connect(myTorus.b, myTorusMixer.a_plasma)
    connect(myFW.output, FWGain.input)
    connect(FWGain.output, myTorusMixer.a_FW)
    connect(myTorusMixer.b, P1.a)
    connect(P1.b, DSS.a)
    connect(DSS.b1, S2.a)
    connect(S2.b2, V2.a2)
    connect(S2.b1, myPR.a)
    connect(myPR.b2, V2.a3)
    connect(DSS.b2, BC_DSS.a)
    connect(myPR.b1, BC_PR.a)
    #
    connect(BC_torusMixer.b, myTorusMixer.a_DC)
    #
    connect(V2.b2, BC_V2.a)
    #
    # Pressure control in the V2
    connect(V2.ND_ref, PIDsource1.reference)
    connect(V2.ND_meas, PIDsource1.measurement)
    connect(PIDsource1.ctr_output, sourceD.a_G)
    #
    # Composition control in the V2
    connect(V2.NT_ref,  PIDsource2.reference)
    connect(V2.NT_meas, PIDsource2.measurement)
    connect(PIDsource2.ctr_output, sourceT.a_G)
    #
    # Pressure control in the V3
    connect(V3.b_p_ref, PID1.reference)
    connect(V3.b_p_meas, PID1.measurement)
    connect(PID1.ctr_output, valve1.a_ctrl)
    #
    connect(valve2.a_ctrl, ramp_FM.output)
  end
end
println("# Input processed ...")



# Build the model. This marks it as complete and simplifies the equation set symbolically.
print("# Building model ...")
@mtkbuild sys = TorusTest() # mtkbuild
print("built.\n")

#tag ICs: Initial conditions for the variables
u0 = [
  sys.sourceD.M => 1.0e+5
  sys.sourceT.M  => 1.0e+5
  #
  sys.sourceXe.M => 1.0e+5
  #
  sys.V2.NH =>  n/1000
  sys.V2.ND =>  (2.0 * p_V2*V_V2 / (8.314 * T_vessels)) / 2.1
  sys.V2.NT =>  (2.0 * p_V2*V_V2 / (8.314 * T_vessels)) / 2.1
  sys.V2.NHe => n/1000
  sys.V2.NXe => n/1000
  sys.V2.NI =>  n/1000
  #
  sys.V3.NH =>  n/1000
  sys.V3.ND =>  (2.0 * p_V3*V_V3 / (8.314 * T_vessels)) / 2.0
  sys.V3.NT =>  (2.0 * p_V3*V_V3 / (8.314 * T_vessels)) / 2.0
  sys.V3.NHe => n/1000
  sys.V3.NXe => n/1000
  sys.V3.NI =>  n/1000
  #
  sys.PI.NH =>  0.01*GammaTarget*tau_PI_Q/eta_fuelingEff
  sys.PI.ND =>  GammaTarget*tau_PI_Q/eta_fuelingEff
  sys.PI.NT =>  GammaTarget*tau_PI_Q/eta_fuelingEff
  sys.PI.NHe => n/1000
  sys.PI.NXe => n/1000
  sys.PI.NI =>  n/1000
  #
  sys.myTorus.nH  => n/100
  sys.myTorus.nD  => n
  sys.myTorus.nT  => n
  sys.myTorus.nHe => n/50
  sys.myTorus.nXe => n/1000
  sys.myTorus.nI  => n/1000
  #
  sys.P1.NH  => 0.05 #n/1000
  sys.P1.ND  => (GammaTarget + 0.5*DTpuff_in/SF)*tau_P1
  sys.P1.NT  => (GammaTarget + 0.5*DTpuff_in/SF)*tau_P1
  sys.P1.NHe => n/1000
  sys.P1.NXe => n/1000
  sys.P1.NI  => (Arpuff_in/SF)*tau_P1
  #
  sys.myPR.NH  => n/100
  sys.myPR.ND  => n
  sys.myPR.NT  => n
  sys.myPR.NHe => n/1000
  sys.myPR.NXe => n/1000
  sys.myPR.NI  => n/1000
]


# ============================================================================ #
#tag Plot_quantities
# ============================================================================ #
plot_quantities = [
  sys.t,
  #
  sys.myTorus.a.G,
  sys.myTorus.a.xH,
  sys.myTorus.a.xD,
  sys.myTorus.a.xT,
  sys.myTorus.a.xHe,
  sys.myTorus.a.xXe,
  sys.myTorus.a.xI,
  sys.myTorus.b.G,
  sys.myTorus.b.xH,
  sys.myTorus.b.xD,
  sys.myTorus.b.xT,
  sys.myTorus.b.xHe,
  sys.myTorus.b.xXe,
  sys.myTorus.b.xI,
  sys.myTorus.nH,
  sys.myTorus.nD,
  sys.myTorus.nT,
  sys.myTorus.nHe,
  sys.myTorus.nXe,
  sys.myTorus.nI,
  sys.myTorus.DNH,
  sys.myTorus.G_alpha,
  sys.myTorus.G_H,
  #
  sys.myTorusMixer.a_FW.u,
  sys.myTorusMixer.a_plasma.G,
  sys.myTorusMixer.a_plasma.xH,
  sys.myTorusMixer.a_plasma.xD,
  sys.myTorusMixer.a_plasma.xT,
  sys.myTorusMixer.a_plasma.xHe,
  sys.myTorusMixer.a_plasma.xXe,
  sys.myTorusMixer.a_plasma.xI,
  sys.myTorusMixer.a_nf.G,
  sys.myTorusMixer.a_nf.xH,
  sys.myTorusMixer.a_nf.xD,
  sys.myTorusMixer.a_nf.xT,
  sys.myTorusMixer.a_nf.xHe,
  sys.myTorusMixer.a_nf.xXe,
  sys.myTorusMixer.a_nf.xI,
  sys.myTorusMixer.a_DC.G,
  sys.myTorusMixer.a_DC.xH,
  sys.myTorusMixer.a_DC.xD,
  sys.myTorusMixer.a_DC.xT,
  sys.myTorusMixer.a_DC.xHe,
  sys.myTorusMixer.a_DC.xXe,
  sys.myTorusMixer.a_DC.xI,
  sys.myTorusMixer.b.G,
  sys.myTorusMixer.b.xH,
  sys.myTorusMixer.b.xD,
  sys.myTorusMixer.b.xT,
  sys.myTorusMixer.b.xHe,
  sys.myTorusMixer.b.xXe,
  sys.myTorusMixer.b.xI,
  #
  sys.sourceD.b.G,
  sys.sourceD.b.xH,
  sys.sourceD.b.xD,
  sys.sourceD.b.xT,
  sys.sourceD.b.xHe,
  sys.sourceD.b.xXe,
  sys.sourceD.b.xI,
  sys.sourceD.M,
  #
  sys.sourceT.b.G,
  sys.sourceT.b.xH,
  sys.sourceT.b.xD,
  sys.sourceT.b.xT,
  sys.sourceT.b.xHe,
  sys.sourceT.b.xXe,
  sys.sourceT.b.xI,
  sys.sourceT.M,
  #
  sys.sourceXe.b.G,
  sys.sourceXe.b.xH,
  sys.sourceXe.b.xD,
  sys.sourceXe.b.xT,
  sys.sourceXe.b.xHe,
  sys.sourceXe.b.xXe,
  sys.sourceXe.b.xI,
  sys.sourceXe.M,
  #
  sys.source_MFC.a1.G,
  sys.source_MFC.a1.xH,
  sys.source_MFC.a1.xD,
  sys.source_MFC.a1.xT,
  sys.source_MFC.a1.xHe,
  sys.source_MFC.a1.xXe,
  sys.source_MFC.a1.xI,
  sys.source_MFC.a2.G,
  sys.source_MFC.a2.xH,
  sys.source_MFC.a2.xD,
  sys.source_MFC.a2.xT,
  sys.source_MFC.a2.xHe,
  sys.source_MFC.a2.xXe,
  sys.source_MFC.a2.xI,
  sys.source_MFC.b.G,
  sys.source_MFC.b.xH,
  sys.source_MFC.b.xD,
  sys.source_MFC.b.xT,
  sys.source_MFC.b.xHe,
  sys.source_MFC.b.xXe,
  sys.source_MFC.b.xI,
  #
  sys.DSS.a.G,
  sys.DSS.a.xH,
  sys.DSS.a.xD,
  sys.DSS.a.xT,
  sys.DSS.a.xHe,
  sys.DSS.a.xXe,
  sys.DSS.a.xI,
  sys.DSS.b1.G,
  sys.DSS.b1.xH,
  sys.DSS.b1.xD,
  sys.DSS.b1.xT,
  sys.DSS.b1.xHe,
  sys.DSS.b1.xXe,
  sys.DSS.b1.xI,
  sys.DSS.b2.G,
  sys.DSS.b2.xH,
  sys.DSS.b2.xD,
  sys.DSS.b2.xT,
  sys.DSS.b2.xHe,
  sys.DSS.b2.xXe,
  sys.DSS.b2.xI,
  #
  sys.PI.NH,
  sys.PI.ND,
  sys.PI.NT,
  sys.PI.NHe,
  sys.PI.NXe,
  sys.PI.NI,
  sys.PI.a.G,
  sys.PI.a.xH,
  sys.PI.a.xD,
  sys.PI.a.xT,
  sys.PI.a.xHe,
  sys.PI.a.xXe,
  sys.PI.a.xI,
  sys.PI.b1.G,
  sys.PI.b1.xH,
  sys.PI.b1.xD,
  sys.PI.b1.xT,
  sys.PI.b1.xHe,
  sys.PI.b1.xXe,
  sys.PI.b1.xI,
  sys.PI.b2.G,
  sys.PI.b2.xH,
  sys.PI.b2.xD,
  sys.PI.b2.xT,
  sys.PI.b2.xHe,
  sys.PI.b2.xXe,
  sys.PI.b2.xI,
  #
  sys.P1.NH,
  sys.P1.ND,
  sys.P1.NT,
  sys.P1.NHe,
  sys.P1.NXe,
  sys.P1.NI,
  sys.P1.a.G,
  sys.P1.a.xH,
  sys.P1.a.xD,
  sys.P1.a.xT,
  sys.P1.a.xHe,
  sys.P1.a.xXe,
  sys.P1.a.xI,
  sys.P1.b.G,
  sys.P1.b.xH,
  sys.P1.b.xD,
  sys.P1.b.xT,
  sys.P1.b.xHe,
  sys.P1.b.xXe,
  sys.P1.b.xI,
  #
  sys.myPR.NH,
  sys.myPR.ND,
  sys.myPR.NT,
  sys.myPR.NHe,
  sys.myPR.NXe,
  sys.myPR.NI,
  sys.myPR.a.G,
  sys.myPR.a.xH,
  sys.myPR.a.xD,
  sys.myPR.a.xT,
  sys.myPR.a.xHe,
  sys.myPR.a.xXe,
  sys.myPR.a.xI,
  sys.myPR.b1.G,
  sys.myPR.b1.xH,
  sys.myPR.b1.xD,
  sys.myPR.b1.xT,
  sys.myPR.b1.xHe,
  sys.myPR.b1.xXe,
  sys.myPR.b1.xI,
  sys.myPR.b2.G,
  sys.myPR.b2.xH,
  sys.myPR.b2.xD,
  sys.myPR.b2.xT,
  sys.myPR.b2.xHe,
  sys.myPR.b2.xXe,
  sys.myPR.b2.xI,
  #
  sys.V2.NH,
  sys.V2.ND,
  sys.V2.NT,
  sys.V2.NHe,
  sys.V2.NXe,
  sys.V2.NI,
  sys.V2.a1.G,
  sys.V2.a1.xH,
  sys.V2.a1.xD,
  sys.V2.a1.xT,
  sys.V2.a1.xHe,
  sys.V2.a1.xXe,
  sys.V2.a1.xI,
  sys.V2.a2.G,
  sys.V2.a2.xH,
  sys.V2.a2.xD,
  sys.V2.a2.xT,
  sys.V2.a2.xHe,
  sys.V2.a2.xXe,
  sys.V2.a2.xI,
  sys.V2.a3.G,
  sys.V2.a3.xH,
  sys.V2.a3.xD,
  sys.V2.a3.xT,
  sys.V2.a3.xHe,
  sys.V2.a3.xXe,
  sys.V2.a3.xI,
  sys.V2.b1.G,
  sys.V2.b1.xH,
  sys.V2.b1.xD,
  sys.V2.b1.xT,
  sys.V2.b1.xHe,
  sys.V2.b1.xXe,
  sys.V2.b1.xI,
  sys.V2.b2.G,
  sys.V2.b2.xH,
  sys.V2.b2.xD,
  sys.V2.b2.xT,
  sys.V2.b2.xHe,
  sys.V2.b2.xXe,
  sys.V2.b2.xI,
  #
  sys.V3.NH,
  sys.V3.ND,
  sys.V3.NT,
  sys.V3.NHe,
  sys.V3.NXe,
  sys.V3.NI,
  sys.V3.a.G,
  sys.V3.a.xH,
  sys.V3.a.xD,
  sys.V3.a.xT,
  sys.V3.a.xHe,
  sys.V3.a.xXe,
  sys.V3.a.xI,
  sys.V3.b.G,
  sys.V3.b.xH,
  sys.V3.b.xD,
  sys.V3.b.xT,
  sys.V3.b.xHe,
  sys.V3.b.xXe,
  sys.V3.b.xI,
  #
  sys.S2.G,
  sys.S2.xH,
  sys.S2.xD,
  sys.S2.xT,
  sys.S2.xHe,
  sys.S2.xXe,
  sys.S2.xI,
  sys.S2.a.G,
  sys.S2.a.xH,
  sys.S2.a.xD,
  sys.S2.a.xT,
  sys.S2.a.xHe,
  sys.S2.a.xXe,
  sys.S2.a.xI,
  sys.S2.b1.G,
  sys.S2.b1.xH,
  sys.S2.b1.xD,
  sys.S2.b1.xT,
  sys.S2.b1.xHe,
  sys.S2.b1.xXe,
  sys.S2.b1.xI,
  sys.S2.b2.G,
  sys.S2.b2.xH,
  sys.S2.b2.xD,
  sys.S2.b2.xT,
  sys.S2.b2.xHe,
  sys.S2.b2.xXe,
  sys.S2.b2.xI,
  #
  sys.FuelingEfficiency.a.G,
  sys.FuelingEfficiency.a.xH,
  sys.FuelingEfficiency.a.xD,
  sys.FuelingEfficiency.a.xT,
  sys.FuelingEfficiency.a.xHe,
  sys.FuelingEfficiency.a.xXe,
  sys.FuelingEfficiency.a.xI,
  sys.FuelingEfficiency.b1.G,
  sys.FuelingEfficiency.b1.xH,
  sys.FuelingEfficiency.b1.xD,
  sys.FuelingEfficiency.b1.xT,
  sys.FuelingEfficiency.b1.xHe,
  sys.FuelingEfficiency.b1.xXe,
  sys.FuelingEfficiency.b1.xI,
  sys.FuelingEfficiency.b2.G,
  sys.FuelingEfficiency.b2.xH,
  sys.FuelingEfficiency.b2.xD,
  sys.FuelingEfficiency.b2.xT,
  sys.FuelingEfficiency.b2.xHe,
  sys.FuelingEfficiency.b2.xXe,
  sys.FuelingEfficiency.b2.xI,
];

#=
equations(expand_connections(sys))
unknowns(sys)
parameters(sys)
observed(sys)
=#
