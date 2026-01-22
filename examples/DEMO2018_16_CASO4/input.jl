using FuelCycleLibrary
using FuelCycleLibrary.PlasmaFunctions
using FuelCycleLibrary.Base.Atomic
using ModelingToolkit
using ModelingToolkit: t_nounits as t, D_nounits as D
using ModelingToolkitStandardLibrary.Blocks
using SciCompDSL
using DataFrames, CSV
using Interpolations


# ============================================================================ #
#tag Simulation parameters
# ============================================================================ #
tspan = (0.0, 60.0*60*24*365) # Simulation time interval

# INPUT variables
R::Float64  = 9.073; 		      # major radius [m]
a::Float64  = 2.927; 		      # minor radius [m]
Vplasma::Float64 = 2579.0;    # plasma volume [m3]
k::Float64 = Vplasma / (2*pi*R * pi*a^2)             # plasma shape factor
Pfus::Float64 = 2.0;          # fusion power [GW]
xaz::Float64 = 0.796;        # active zone (normalized radius) [-]
S_FW::Float64 = 1462.0;        # First Wall surface [m2]
T_torus_out::Float64 = 273.15 # tempeature at the outlet of the torous [K]
τp::Float64 = 31.63; 		      # particle confinement time [s] - estimated from PROCESS extended output (nD + nT)*Vplasma/Γ_fuelingRate
τp_α::Float64 = 18.20; 		      # alpha particle confinement time [s]




# Residence times
tau_P1::Float64      = 2*60.0 
tau_PI_Q::Float64    = 20*60.0
tau_PI_He::Float64   = 20*60.0
tau_PR::Float64      = 40*60.0
tau_GasPuff::Float64 = 30.0
# tau_total = ?

# Efficiencies / separation factors
eta_DSS::Float64          = 0.99975
eta_PI::Float64           = 0.95
eta_PR::Float64           = 0.835
eta_S1::Float64           = 0.81 #0.78 # 0.78
eta_S2::Float64           = 0.2 #
eta_fuelingEff::Float64   = 0.25 # Fueling efficiency

# Gas puffing
DTpuff_in::Float64        = 1.0e23 # [s-1] target D+T (50/50) puff particle flow rate
Arpuff_in::Float64        = 1.0e19

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
#n = n*sqrt(2.0/1.9588585520327826)
# Define scaling faactor. Here we pass from [at] to [mol (monoatomic)]
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
#Power_nAverageCorrectionFactor = Power_nAverageCorrectionFactor*(2.0/1.9856209192025742)
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
    sourceT  = ControlledDTStorage2_tank(τ=1.0e-5)
    ISSaddInerts = HDTO_to_outlet()
    sourceT2  = ControlledDTStorage2(τ=1.0e-5, xH=0.0, xD=0.0, xT=1.0, xHe=0.0, xI=0.0, xXe=0.0)
    source_MFC = MixerThreeWay()

    valve1 = FixedActuatedValve() #ActuatedValve()
    valve2 = FixedActuatedValve()
    valve3 = FixedActuatedValve() #ActuatedValve()
    valve4 = FixedActuatedValve()

    BC_PI_b1 = FluidOpenBC()
    #BCvalve4 = FluidOpenBC()
    #CapV2_a3 = Cap()
    #CapV1_a1 = Cap()

    # TES Blocks
    ramp_TES = Blocks.Ramp(height=304.0*(1/3)*(1/(60*60*24)*Na/SF), duration=10.0) # Spagnuolo (2021) - 304 g/d
    TES = IdealSource(xH=0.0, xD=0.0, xT=1.0, xHe=0.0, xXe=0.0, xI=0.0)

    ramp_FM  = Blocks.Ramp(height=2*GammaTarget/eta_fuelingEff, duration=0.0)
    ramp_FM2 = Blocks.Ramp(height=DTpuff_in/SF, duration=0.0)
    #ramp_idealS = Blocks.Ramp(height=DTpuff_in/(SF*2), duration=10.0)
    #idealS = IdealSource()


    V2 = Vessel3I2O(T=T_vessels, V=V_V2, p_op=p_V2)
    V3 = Vessel1I1O(T=T_vessels, V=V_V3, p_op=p_V3)
    V1 = Vessel2I1O(T=T_vessels, V=V_V1, p_op=p_V1)

    # , k=2000.0, Ti=350.0, Td=50.0)
    PIDsource1 = Blocks.LimPID(u_min=0.0, u_max=40*(DTpuff_in/SF + 2*GammaTarget), k=10.0, Td=0.2)
    PIDsource2 = Blocks.LimPID(u_min=0.0, u_max=40*(DTpuff_in/SF + 2*GammaTarget), k=10.0, Td=0.2)
    PID1 = Blocks.LimPID(u_min=0.0, k=10.0, Td=0.2)
    PID3 = Blocks.LimPID(u_min=0.0, k=10.0, Td=0.2)

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
    S1 = LimitedSplitter(eta=eta_S1, Gmax=0.9999*DTpuff_in/SF)
    DSS = MFP(eta=eta_DSS) # Metal Foil Pump

    S2 = Splitter(eta=eta_S2)
    myPR = IRPR(τ=tau_PR, eta=eta_PR)
    #BC_PR = FluidOpenBC()

    # Argon
    sourceAr = Storage2(xH=0.0, xD=0.0, xT=0.0, xHe=0.0, xI=1.0, xXe=0.0)
    sourceXe = Storage2(xH=0.0, xD=0.0, xT=0.0, xHe=0.0, xI=0.0, xXe=1.0)
    ramp_Ar  = Blocks.Ramp(height=Arpuff_in/SF, duration=0.0)
    ramp_Xe  = Blocks.Ramp(height=fXe*2*GammaTarget/eta_fuelingEff, duration=0.0)
    #PIDsourceAr = Blocks.LimPID(k=100.0, Td=20.0)
    mixerAr = MixerTwoWay()
    mixerXe = MixerTwoWay()
    #FM_Ar = SpeciesFlowMeter(jH=0, jD=0, jT=0, jHe=0, jI=1)

    # OUTL COMPONENTS -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
    #BC_PI_b2 = FluidOpenBC()
    addO1 = inlet_to_outlet2() # interfaccia myMFP - CuO + ZMS in cui viene aggiunta la specie Ossigeno 'O'
    addO2 = inlet_to_outlet2() # interfaccia myMFP - CuO + ZMS in cui viene aggiunta la specie Ossigeno 'O'
    CuOZMS = CuO_ZMS_bed(eta_ox=0.98, eta_zeo=0.98, τ=15.0*60.0*60.0, fHe_purge=0.001, fO_regen=1.0)
    MR = MembraneReactor(η_r=0.60, η_p=0.98 , τ=0.01*60*60)
    #MR = MembraneReactor(η_r=0.98, η_p=0.98 , τ=15.0*60*60)
    MR_addN = inlet3_to_outletOCN()
    PSS_addCN = inlet2_to_outletOCN()
    removeO1 = RemoveOFromFluidPort()
    removeInert = RemoveInertFromFluidPort()
    removeInert2 = RemoveInertFromFluidPort()
    PSS1 = PSS()

    ##
    Mix1_OUTL = MixerThreeWayOCN()
    Mix2_OUTL_MR_DIRL = MixerTwoWay()
    ##
    #BC_OUTL1 = FluidOpenBC()
    BC_OUTL2 = FluidOpenBC_OCN()
    #BC_OUTL5 = FluidOpenBC_HDTO()
    BC_OUTL6 = FluidOpenBC_HDTO()
    ##

    EDS = WetScrubberColumn(
      externalH2O=43.26/SF,
      τH=0.92*60*60, τD=1.22*60*60, τT=7.94*60*60, τO=0.92*60*60,
      etaH=0.0231, etaD=0.0574, etaT=0.99952, etaO=0.0231 # etaT=0.9973
      #etaH=0.0231, etaD=0.0574, etaT=0.9973, etaO=0.0231
    ) # DATI DA Aspen (Vincenzo) MA etaT scelto invece per rilasciare meno di 0.1 gT/anno

    CPS = WaterDistillation(
      externalH2O=180.0/SF,
      τH=7.86*60*60, τD=70.79*60*60, τT=98.07*60*60, τO=7.87*60*60,
      etaH = 0.992813, etaD=0.584318, etaT=0.278140, etaO=0.992781,
    ) # DATI DA Aspen (Vincenzo)

    WDS = LPCE(
      τH=1.739*60*60, τD=1.720*60*60, τT=0.139*60*60, τO=0.869*60*60,
      etaH=0.75041, etaD=0.75308, etaT=1.0, etaO=0.0,
    ) # ETA E TAU DA RIN + VALUTAZIONI Sulcor Vincenzo - Excel "risultati caso 4", folgio "POSTPROCESSING", tabella "LPCE PARAMETERS"

    ISS = CryogenicDistillation(
      τH=0.284*60*60, τD=4.649*60*60, τT=1.182*60*60,
      etaH=0.0000, etaD=0.50919, etaT=1.0,
    ) # ETA E TAU DA RIN + VALUTAZIONI Sulcor Vincenzo - Excel "risultati caso 4", folgio "POSTPROCESSING", tabella "CD PARAMETERS"

    ##
    ## D natural abundance: 0.0156% of H2 in water
    #capEDS          = CapOCN()
    AirFlow         = FlowMeterOCN(G=(54.529*Na)/SF) # mol/s di N
    sourceAir       = StorageOCN(xH=0.0, xD=0.0, xT=0.0, xHe=0.0, xI=0.0, xXe=0.0, xO=0.0, xC=0.0, xN=1.0)

    sourceH2O_1     = Storage2_withO(xH=0.66656, xD=0.0001038, xT=0.0, xHe=0.0, xI=0.0, xXe=0.0, xO=(1.0 - 0.66656 - 0.0001038))
    sourceH2O_3     = StorageHDTO(xH=0.666614459, xD=5.2e-05, xT=0.000000208, xO=(1.0 - 0.666614459 - 5.2e-05 - 0.000000208))
    sourceCO_MR     = Storage3(xH=0.0, xD=0.0, xT=0.0, xHe=0.0, xI=0.0, xXe=0.0, xO=0.5, xC=0.5)
    sourceHe_CuOZMS = Storage2_withO(xH=0.0, xD=0.0, xT=0.0, xHe=1.0, xI=0.0, xXe=0.0, xO=0.0)
    sourceO_CuOZMS  = Storage2_withO(xH=0.0, xD=0.0, xT=0.0, xHe=0.0, xI=0.0, xXe=0.0, xO=1.0)

    P_ISS = PumpHDTO(τ=20*60.0)

  end

  @equations begin
    connect(sourceXe.b, mixerXe.a2)
    connect(ramp_Xe.output, sourceXe.a_G)

    connect(sourceD.b, source_MFC.a1)
    connect(sourceT.b, source_MFC.a2)
    connect(sourceT2.b, source_MFC.a3)
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
    #connect(valve4.b, myTorusMixer.a_DC)
    connect(myTorusMixer.b, P1.a)
    connect(P1.b, S1.a)
    connect(S1.b1, V1.a2)
    connect(S1.b2, DSS.a)
    connect(DSS.b1, Mix2_OUTL_MR_DIRL.a1)
    connect(Mix2_OUTL_MR_DIRL.b, S2.a)
    connect(S2.b2, V2.a2)
    connect(S2.b1, myPR.a)
    connect(myPR.b2, V2.a3)
    connect(DSS.b2, addO1.a)
    #



    connect(addO1.b, CuOZMS.a1)
    connect(sourceO_CuOZMS.b, CuOZMS.a3)
    connect(sourceHe_CuOZMS.b, CuOZMS.a2)
    connect(CuOZMS.b2, PSS1.a)
    connect(PSS1.b, PSS_addCN.a)
    connect(PSS_addCN.b, Mix1_OUTL.a1)
    connect(CuOZMS.b1, MR.a1)
    connect(sourceCO_MR.b, MR.a2)
    connect(MR.b1, removeO1.a)
    connect(removeO1.b, Mix2_OUTL_MR_DIRL.a2)
    connect(MR.b2, MR_addN.a)
    connect(MR_addN.b, Mix1_OUTL.a2)

    connect(sourceAir.b, AirFlow.a)
    connect(AirFlow.b, Mix1_OUTL.a3)
    #connect(capEDS.b, Mix1_OUTL.a3)


    connect(Mix1_OUTL.b, EDS.a1)
    connect(sourceH2O_1.b, EDS.a2)
    connect(EDS.b2, BC_OUTL2.a)


    connect(EDS.b1, removeInert.a)
    connect(removeInert.b, CPS.a2)
    connect(WDS.b2, EDS.a3)
    connect(WDS.b1, ISS.a1)

    connect(ISS.b2, P_ISS.a)
    connect(P_ISS.b, WDS.a2)

    connect(myPR.b1, addO2.a)
    connect(addO2.b, removeInert2.a)
    connect(removeInert2.b, ISS.a2)
    #connect(ISS.b1, BC_OUTL5.a)
    connect(ISS.b1, ISSaddInerts.a)
    connect(ISSaddInerts.b, sourceT.a1)
    connect(TES.b, sourceT.a2)
    connect(TES.G, ramp_TES.output)


    connect(sourceH2O_3.b, CPS.a1)
    connect(CPS.b1, BC_OUTL6.a)
    connect(CPS.b2, WDS.a1)

    #connect(myPR.b1, BCtmp3.a)

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    connect(valve4.b, mixerAr.a1)
    connect(sourceAr.b, mixerAr.a2)
    connect(mixerAr.b, myTorusMixer.a_DC)
    connect(ramp_Ar.output, sourceAr.a_G)
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    connect(V2.b2, valve3.a)
    connect(valve3.b, V1.a1)
    connect(V1.b, valve4.a)
    connect(valve4.a_ctrl, ramp_FM2.output)
    connect(V1.b_p_ref, PID3.reference)
    connect(V1.b_p_meas, PID3.measurement)
    connect(PID3.ctr_output, valve3.a_ctrl)

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
    #
    connect(sourceT.error_GT,  sourceT2.a_G)
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

  sys.sourceT.NH  => 0.0
  sys.sourceT.ND  => 1.0
  sys.sourceT.NT  => 1.0
  sys.sourceT.NHe => 0.0
  sys.sourceT.NI  => 0.0
  sys.sourceT.NXe => 0.0

  sys.sourceT2.M => 1.0e+5


  sys.sourceAr.M => 1.0e+5
  sys.sourceXe.M => 1.0e+5
  sys.sourceO_CuOZMS.M => 1.0e+5
  sys.sourceHe_CuOZMS.M => 1.0e+5
  sys.sourceCO_MR.M => 1.0e+5
  sys.sourceH2O_1.M => 1.0e+5
  sys.sourceH2O_3.M => 1.0e+5
  sys.sourceAir.M => 1.0e+5
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
  sys.V1.NH =>  n/1000
  sys.V1.ND =>  (2.0 * p_V1*V_V1 / (8.314 * T_vessels)) / 2.0
  sys.V1.NT =>  (2.0 * p_V1*V_V1 / (8.314 * T_vessels)) / 2.0
  sys.V1.NHe => n/1000
  sys.V1.NXe => n/1000
  sys.V1.NI =>  n/1000
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
  #
  sys.CuOZMS.NH  => n/1000
  sys.CuOZMS.ND  => n
  sys.CuOZMS.NT  => n #1.0e-10
  sys.CuOZMS.NHe => 45 #n
  sys.CuOZMS.NXe => n/1000
  sys.CuOZMS.NI  => n/1000
  sys.CuOZMS.NO  => n
  #
  sys.MR.NH  => n/1000
  sys.MR.ND  => n
  sys.MR.NT  => n #1.0e-10
  sys.MR.NHe => n
  sys.MR.NXe => n/1000
  sys.MR.NI  => n/1000
  sys.MR.NO  => n/1000
  sys.MR.NC  => 0.7 #n/1000
  #
  sys.PSS1.NH  => n/1000
  sys.PSS1.ND  => n/1000
  sys.PSS1.NT  => n/1000 #1.0e-10
  sys.PSS1.NHe => 3.0 #n/1000
  sys.PSS1.NXe => n/1000
  sys.PSS1.NI  => n/1000
  sys.PSS1.NO  => n/1000
  #
  sys.EDS.NH  => 4400 #100*n
  sys.EDS.ND  => 1.0 #n/10
  sys.EDS.NT  => n/10 #1.0e-10
  sys.EDS.NHe => n/1000
  sys.EDS.NXe => n/1000
  sys.EDS.NI  => n/1000
  sys.EDS.NO  => 2200 #200*n
  #sys.EDS.NC  => n/1000
  #sys.EDS.NN  => n/1000
  #
  sys.CPS.NH  => 8100 #500*n
  sys.CPS.ND  => 41 #n/1000
  sys.CPS.NT  => n/1000 #1.0e-10
  sys.CPS.NO  => 4000 #1000*n
  #
  sys.WDS.NH  => 500*n
  sys.WDS.ND  => n/1000
  sys.WDS.NT  => n/1000 #1.0e-10
  sys.WDS.NO  => 1000*n
  #
  sys.ISS.NH  => 1500 #n/10 #1.0e-10
  sys.ISS.ND  => 303 #n/10 #1.0e-10
  sys.ISS.NT  => 53 #n/10 #1.0e-10
  sys.ISS.NO  => 0.0
  #
  sys.P_ISS.NH  => 70 #n/10 #1.0e-10
  sys.P_ISS.ND  => 15 #n/10 #1.0e-10
  sys.P_ISS.NT  => n/10 #1.0e-10
  sys.P_ISS.NO  => 0.0

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
  sys.sourceT.NH,
  sys.sourceT.ND,
  sys.sourceT.NT,
  #
  sys.sourceT2.b.G,
  sys.sourceT2.b.xH,
  sys.sourceT2.b.xD,
  sys.sourceT2.b.xT,
  sys.sourceT2.b.xHe,
  sys.sourceT2.b.xXe,
  sys.sourceT2.b.xI,
  sys.sourceT2.M,
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
  sys.sourceAr.b.G,
  sys.sourceAr.b.xI,
  sys.sourceAr.M,
  #
  sys.sourceCO_MR.b.G,
  sys.sourceCO_MR.b.xO,
  sys.sourceCO_MR.b.xC,
  sys.sourceCO_MR.M,
  #
  sys.sourceHe_CuOZMS.b.G,
  sys.sourceHe_CuOZMS.b.xHe,
  sys.sourceHe_CuOZMS.M,
  #
  sys.sourceO_CuOZMS.b.G,
  sys.sourceO_CuOZMS.b.xO,
  sys.sourceO_CuOZMS.M,
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
  sys.V1.NH,
  sys.V1.ND,
  sys.V1.NT,
  sys.V1.NHe,
  sys.V1.NXe,
  sys.V1.NI,
  sys.V1.a1.G,
  sys.V1.a1.xH,
  sys.V1.a1.xD,
  sys.V1.a1.xT,
  sys.V1.a1.xHe,
  sys.V1.a1.xXe,
  sys.V1.a1.xI,
  sys.V1.a2.G,
  sys.V1.a2.xH,
  sys.V1.a2.xD,
  sys.V1.a2.xT,
  sys.V1.a2.xHe,
  sys.V1.a2.xXe,
  sys.V1.a2.xI,
  sys.V1.b.G,
  sys.V1.b.xH,
  sys.V1.b.xD,
  sys.V1.b.xT,
  sys.V1.b.xHe,
  sys.V1.b.xXe,
  sys.V1.b.xI,
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
  sys.S1.G,
  sys.S1.xH,
  sys.S1.xD,
  sys.S1.xT,
  sys.S1.xHe,
  sys.S1.xXe,
  sys.S1.xI,
  sys.S1.a.G,
  sys.S1.a.xH,
  sys.S1.a.xD,
  sys.S1.a.xT,
  sys.S1.a.xHe,
  sys.S1.a.xXe,
  sys.S1.a.xI,
  sys.S1.b1.G,
  sys.S1.b1.xH,
  sys.S1.b1.xD,
  sys.S1.b1.xT,
  sys.S1.b1.xHe,
  sys.S1.b1.xXe,
  sys.S1.b1.xI,
  sys.S1.b2.G,
  sys.S1.b2.xH,
  sys.S1.b2.xD,
  sys.S1.b2.xT,
  sys.S1.b2.xHe,
  sys.S1.b2.xXe,
  sys.S1.b2.xI,
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
  #
  #
  sys.PSS1.NH,
  sys.PSS1.ND,
  sys.PSS1.NT,
  sys.PSS1.NHe,
  sys.PSS1.NXe,
  sys.PSS1.NI,
  sys.PSS1.NO,
  sys.PSS1.a.G,
  sys.PSS1.a.xH,
  sys.PSS1.a.xD,
  sys.PSS1.a.xT,
  sys.PSS1.a.xHe,
  sys.PSS1.a.xXe,
  sys.PSS1.a.xI,
  sys.PSS1.a.xO,
  sys.PSS1.b.G,
  sys.PSS1.b.xH,
  sys.PSS1.b.xD,
  sys.PSS1.b.xT,
  sys.PSS1.b.xHe,
  sys.PSS1.b.xXe,
  sys.PSS1.b.xI,
  sys.PSS1.b.xO,
  #
  sys.CuOZMS.NH,
  sys.CuOZMS.ND,
  sys.CuOZMS.NT,
  sys.CuOZMS.NHe,
  sys.CuOZMS.NXe,
  sys.CuOZMS.NI,
  sys.CuOZMS.NO,
  sys.CuOZMS.a1.G,
  sys.CuOZMS.a1.xH,
  sys.CuOZMS.a1.xD,
  sys.CuOZMS.a1.xT,
  sys.CuOZMS.a1.xHe,
  sys.CuOZMS.a1.xXe,
  sys.CuOZMS.a1.xI,
  sys.CuOZMS.a1.xO,
  sys.CuOZMS.a2.G,
  sys.CuOZMS.a2.xH,
  sys.CuOZMS.a2.xD,
  sys.CuOZMS.a2.xT,
  sys.CuOZMS.a2.xHe,
  sys.CuOZMS.a2.xXe,
  sys.CuOZMS.a2.xI,
  sys.CuOZMS.a2.xO,
  sys.CuOZMS.a3.G,
  sys.CuOZMS.a3.xH,
  sys.CuOZMS.a3.xD,
  sys.CuOZMS.a3.xT,
  sys.CuOZMS.a3.xHe,
  sys.CuOZMS.a3.xXe,
  sys.CuOZMS.a3.xI,
  sys.CuOZMS.a3.xO,
  sys.CuOZMS.b1.G,
  sys.CuOZMS.b1.xH,
  sys.CuOZMS.b1.xD,
  sys.CuOZMS.b1.xT,
  sys.CuOZMS.b1.xHe,
  sys.CuOZMS.b1.xXe,
  sys.CuOZMS.b1.xI,
  sys.CuOZMS.b1.xO,
  sys.CuOZMS.b2.G,
  sys.CuOZMS.b2.xH,
  sys.CuOZMS.b2.xD,
  sys.CuOZMS.b2.xT,
  sys.CuOZMS.b2.xHe,
  sys.CuOZMS.b2.xXe,
  sys.CuOZMS.b2.xI,
  sys.CuOZMS.b2.xO,
  #
  sys.EDS.NH,
  sys.EDS.ND,
  sys.EDS.NT,
  sys.EDS.NHe,
  sys.EDS.NXe,
  sys.EDS.NI,
  sys.EDS.NO,
  sys.EDS.a1.G,
  sys.EDS.a1.xH,
  sys.EDS.a1.xD,
  sys.EDS.a1.xT,
  sys.EDS.a1.xHe,
  sys.EDS.a1.xXe,
  sys.EDS.a1.xI,
  sys.EDS.a1.xO,
  sys.EDS.a1.xC,
  sys.EDS.a1.xN,
  sys.EDS.a2.G,
  sys.EDS.a2.xH,
  sys.EDS.a2.xD,
  sys.EDS.a2.xT,
  sys.EDS.a2.xHe,
  sys.EDS.a2.xXe,
  sys.EDS.a2.xI,
  sys.EDS.a2.xO,
  sys.EDS.a3.G,
  sys.EDS.a3.xH,
  sys.EDS.a3.xD,
  sys.EDS.a3.xT,
  sys.EDS.a3.xO,
  sys.EDS.b1.G,
  sys.EDS.b1.xH,
  sys.EDS.b1.xD,
  sys.EDS.b1.xT,
  sys.EDS.b1.xHe,
  sys.EDS.b1.xXe,
  sys.EDS.b1.xI,
  sys.EDS.b1.xO,
  sys.EDS.b2.G,
  sys.EDS.b2.xH,
  sys.EDS.b2.xD,
  sys.EDS.b2.xT,
  sys.EDS.b2.xHe,
  sys.EDS.b2.xXe,
  sys.EDS.b2.xI,
  sys.EDS.b2.xO,
  sys.EDS.b2.xC,
  sys.EDS.b2.xN,
  #
  sys.MR.NH,
  sys.MR.ND,
  sys.MR.NT,
  sys.MR.NHe,
  sys.MR.NXe,
  sys.MR.NI,
  sys.MR.NO,
  sys.MR.NC,
  sys.MR.a1.G,
  sys.MR.a1.xH,
  sys.MR.a1.xD,
  sys.MR.a1.xT,
  sys.MR.a1.xHe,
  sys.MR.a1.xXe,
  sys.MR.a1.xI,
  sys.MR.a1.xO,
  sys.MR.a2.G,
  sys.MR.a2.xH,
  sys.MR.a2.xD,
  sys.MR.a2.xT,
  sys.MR.a2.xHe,
  sys.MR.a2.xXe,
  sys.MR.a2.xI,
  sys.MR.a2.xO,
  sys.MR.a2.xC,
  sys.MR.b1.G,
  sys.MR.b1.xH,
  sys.MR.b1.xD,
  sys.MR.b1.xT,
  sys.MR.b1.xHe,
  sys.MR.b1.xXe,
  sys.MR.b1.xI,
  sys.MR.b1.xO,
  sys.MR.b2.G,
  sys.MR.b2.xH,
  sys.MR.b2.xD,
  sys.MR.b2.xT,
  sys.MR.b2.xHe,
  sys.MR.b2.xXe,
  sys.MR.b2.xI,
  sys.MR.b2.xO,
  sys.MR.b2.xC,
  #
  sys.CPS.NH,
  sys.CPS.ND,
  sys.CPS.NT,
  sys.CPS.NO,
  sys.CPS.a1.G,
  sys.CPS.a1.xH,
  sys.CPS.a1.xD,
  sys.CPS.a1.xT,
  sys.CPS.a1.xO,
  sys.CPS.a2.G,
  sys.CPS.a2.xH,
  sys.CPS.a2.xD,
  sys.CPS.a2.xT,
  sys.CPS.a2.xO,
  sys.CPS.b1.G,
  sys.CPS.b1.xH,
  sys.CPS.b1.xD,
  sys.CPS.b1.xT,
  sys.CPS.b1.xO,
  sys.CPS.b2.G,
  sys.CPS.b2.xH,
  sys.CPS.b2.xD,
  sys.CPS.b2.xT,
  sys.CPS.b2.xO,
  #
  sys.ISS.NH,
  sys.ISS.ND,
  sys.ISS.NT,
  sys.ISS.NO,
  sys.ISS.a1.G,
  sys.ISS.a1.xH,
  sys.ISS.a1.xD,
  sys.ISS.a1.xT,
  sys.ISS.a1.xO,
  sys.ISS.a2.G,
  sys.ISS.a2.xH,
  sys.ISS.a2.xD,
  sys.ISS.a2.xT,
  sys.ISS.a2.xO,
  sys.ISS.b1.G,
  sys.ISS.b1.xH,
  sys.ISS.b1.xD,
  sys.ISS.b1.xT,
  sys.ISS.b1.xO,
  sys.ISS.b2.G,
  sys.ISS.b2.xH,
  sys.ISS.b2.xD,
  sys.ISS.b2.xT,
  sys.ISS.b2.xO,
  #
  sys.WDS.NH,
  sys.WDS.ND,
  sys.WDS.NT,
  sys.WDS.NO,
  sys.WDS.a1.G,
  sys.WDS.a1.xH,
  sys.WDS.a1.xD,
  sys.WDS.a1.xT,
  sys.WDS.a1.xO,
  sys.WDS.a2.G,
  sys.WDS.a2.xH,
  sys.WDS.a2.xD,
  sys.WDS.a2.xT,
  sys.WDS.a2.xO,
  sys.WDS.b1.G,
  sys.WDS.b1.xH,
  sys.WDS.b1.xD,
  sys.WDS.b1.xT,
  sys.WDS.b1.xO,
  sys.WDS.b2.G,
  sys.WDS.b2.xH,
  sys.WDS.b2.xD,
  sys.WDS.b2.xT,
  sys.WDS.b2.xO,
  #
  sys.Mix1_OUTL.a1.G,
  sys.Mix1_OUTL.a1.xH,
  sys.Mix1_OUTL.a1.xD,
  sys.Mix1_OUTL.a1.xT,
  sys.Mix1_OUTL.a1.xHe,
  sys.Mix1_OUTL.a1.xI,
  sys.Mix1_OUTL.a1.xXe,
  sys.Mix1_OUTL.a1.xO,
  sys.Mix1_OUTL.a1.xC,
  sys.Mix1_OUTL.a1.xN,
  sys.Mix1_OUTL.a3.G,
  sys.Mix1_OUTL.a3.xH,
  sys.Mix1_OUTL.a3.xD,
  sys.Mix1_OUTL.a3.xT,
  sys.Mix1_OUTL.a3.xHe,
  sys.Mix1_OUTL.a3.xI,
  sys.Mix1_OUTL.a3.xXe,
  sys.Mix1_OUTL.a3.xO,
  sys.Mix1_OUTL.a3.xC,
  sys.Mix1_OUTL.a3.xN,
  sys.Mix1_OUTL.a2.G,
  sys.Mix1_OUTL.a2.xH,
  sys.Mix1_OUTL.a2.xD,
  sys.Mix1_OUTL.a2.xT,
  sys.Mix1_OUTL.a2.xHe,
  sys.Mix1_OUTL.a2.xI,
  sys.Mix1_OUTL.a2.xXe,
  sys.Mix1_OUTL.a2.xO,
  sys.Mix1_OUTL.a2.xC,
  sys.Mix1_OUTL.a2.xN,
  sys.Mix1_OUTL.b.G,
  sys.Mix1_OUTL.b.xH,
  sys.Mix1_OUTL.b.xD,
  sys.Mix1_OUTL.b.xT,
  sys.Mix1_OUTL.b.xHe,
  sys.Mix1_OUTL.b.xI,
  sys.Mix1_OUTL.b.xXe,
  sys.Mix1_OUTL.b.xO,
  sys.Mix1_OUTL.b.xC,
  sys.Mix1_OUTL.b.xN,
  #
  sys.TES.b.G,
];

#=
equations(expand_connections(sys))
unknowns(sys)
parameters(sys)
observed(sys)
=#
