Na::Float64 = 6.02214076e+23 # Avogadro number [mol-1]
Rg::Float64 = 8.314 # Ideal gas costant [J mol-1 K-1]
λT::Float64 = 1.782607697460566e-9 # tritium decay rate [s-1]

# Define square root function that avoids the square of negative numbers
regPow(x, a, delta=0.001) = x * (x * x + delta * delta)^((a - 1) / 2);
regRoot(x, delta=0.001) = regPow(x, 0.5, delta)

# ============================================================================ #
#tag AFPInlet, AFPOutlet
# ============================================================================ #
# Connectors for pressure line
@doc raw"""
    AFPInlet(;name)
Inlet connector for species H, D, T, He, Xe and one inert 'I'.
# States:
- `G(t)`: [``\mathrm{s^{-1}}``] flow rate
- `xH(t)`: [-] atomic fraction
- `xD(t)`: [-] atomic fraction
- `xT(t)`: [-] atomic fraction
- `xHe(t)`: [-] atomic fraction
- `xXe(t)`: [-] atomic fraction
- `xI(t)`: [-] atomic fraction
- `p(t)`: [``\mathrm{Pa}``] pressure
- `T(t)`: [``\mathrm{K}``] temperature
"""
@connector AFPInlet begin
  dummy(t)
  G(t), [guess = 0.0, connect = Flow, bounds = (0.0, Inf)]
  xH(t), [guess = 0.0, bounds = (0.0, 1.0), input = true]
  xD(t), [guess = 0.5, bounds = (0.0, 1.0), input = true]
  xT(t), [guess = 0.5, bounds = (0.0, 1.0), input = true]
  xHe(t), [guess = 0.0, bounds = (0.0, 1.0), input = true]
  xXe(t), [guess = 0.0, bounds = (0.0, 1.0), input = true]
  xI(t), [guess = 0.0, bounds = (0.0, 1.0), input = true]
  #p(t), [guess = 1.0e5, bounds = (0.0, Inf), input = true]
  #T(t), [guess = 300.0, bounds = (0.0, Inf), input = true]
end

@doc raw"""
    AFPOutlet(;name)
Outlet connector for species H, D, T, He, Xe and one inert 'I'.
# States:
- `G(t)`: [``\mathrm{s^{-1}}``] flow rate
- `xH(t)`: [-] atomic fraction
- `xD(t)`: [-] atomic fraction
- `xT(t)`: [-] atomic fraction
- `xHe(t)`: [-] atomic fraction
- `xXe(t)`: [-] atomic fraction
- `xI(t)`: [-] atomic fraction
- `p(t)`: [``\mathrm{Pa}``] pressure
- `T(t)`: [``\mathrm{K}``] temperature
"""
@connector AFPOutlet begin
  dummy(t)
  G(t), [guess = 0.0, connect = Flow, bounds = (0.0, Inf)]
  xH(t), [guess = 0.0, bounds = (0.0, 1.0), output = true]
  xD(t), [guess = 0.5, bounds = (0.0, 1.0), output = true]
  xT(t), [guess = 0.5, bounds = (0.0, 1.0), output = true]
  xHe(t), [guess = 0.0, bounds = (0.0, 1.0), output = true]
  xXe(t), [guess = 0.0, bounds = (0.0, 1.0), output = true]
  xI(t), [guess = 0.0, bounds = (0.0, 1.0), output = true]
  #p(t), [guess = 1.0e5, bounds = (0.0, Inf), output = true]
  #T(t), [guess = 300.0, bounds = (0.0, Inf), output = true]
end

@doc raw"""
    AFPInlet_O(;name)
Inlet connector for species H, D, T, He, Xe, O and one inert 'I'.
# States:
- `G(t)`: [``\mathrm{s^{-1}}``] flow rate
- `xH(t)`: [-] atomic fraction
- `xD(t)`: [-] atomic fraction
- `xT(t)`: [-] atomic fraction
- `xHe(t)`: [-] atomic fraction
- `xXe(t)`: [-] atomic fraction
- `xI(t)`: [-] atomic fraction
- `xO(t)`: [-] atomic fraction
- `p(t)`: [``\mathrm{Pa}``] pressure
- `T(t)`: [``\mathrm{K}``] temperature
"""
@connector AFPInlet_O begin
  dummy(t)
  G(t), [guess = 0.0, connect = Flow, bounds = (0.0, Inf)]
  xH(t), [guess = 0.0, bounds = (0.0, 1.0), input = true]
  xD(t), [guess = 0.5, bounds = (0.0, 1.0), input = true]
  xT(t), [guess = 0.5, bounds = (0.0, 1.0), input = true]
  xHe(t), [guess = 0.0, bounds = (0.0, 1.0), input = true]
  xXe(t), [guess = 0.0, bounds = (0.0, 1.0), input = true]
  xI(t), [guess = 0.0, bounds = (0.0, 1.0), input = true]
  xO(t), [guess = 0.0, bounds = (0.0, 1.0), input = true]
  #p(t), [guess = 1.0e5, bounds = (0.0, Inf), input = true]
  #T(t), [guess = 300.0, bounds = (0.0, Inf), input = true]
end

@doc raw"""
    AFPOutlet_O(;name)
Outlet connector for species H, D, T, He, Xe, O and one inert 'I'.
# States:
- `G(t)`: [``\mathrm{s^{-1}}``] flow rate
- `xH(t)`: [-] atomic fraction
- `xD(t)`: [-] atomic fraction
- `xT(t)`: [-] atomic fraction
- `xHe(t)`: [-] atomic fraction
- `xXe(t)`: [-] atomic fraction
- `xI(t)`: [-] atomic fraction
- `xO(t)`: [-] atomic fraction
- `p(t)`: [``\mathrm{Pa}``] pressure
- `T(t)`: [``\mathrm{K}``] temperature
"""
@connector AFPOutlet_O begin
  dummy(t)
  G(t), [guess = 0.0, connect = Flow, bounds = (0.0, Inf)]
  xH(t), [guess = 0.0, bounds = (0.0, 1.0), output = true]
  xD(t), [guess = 0.5, bounds = (0.0, 1.0), output = true]
  xT(t), [guess = 0.5, bounds = (0.0, 1.0), output = true]
  xHe(t), [guess = 0.0, bounds = (0.0, 1.0), output = true]
  xXe(t), [guess = 0.0, bounds = (0.0, 1.0), output = true]
  xI(t), [guess = 0.0, bounds = (0.0, 1.0), output = true]
  xO(t), [guess = 0.0, bounds = (0.0, 1.0), output = true]
  #p(t), [guess = 1.0e5, bounds = (0.0, Inf), output = true]
  #T(t), [guess = 300.0, bounds = (0.0, Inf), output = true]
end

@doc raw"""
    AFPInlet_HDTO(;name)
Inlet connector for species H, D, T, O.
# States:
- `G(t)`: [``\mathrm{s^{-1}}``] flow rate
- `xH(t)`: [-] atomic fraction
- `xD(t)`: [-] atomic fraction
- `xT(t)`: [-] atomic fraction
- `xO(t)`: [-] atomic fraction
- `p(t)`: [``\mathrm{Pa}``] pressure
- `T(t)`: [``\mathrm{K}``] temperature
"""
@connector AFPInlet_HDTO begin
  dummy(t)
  G(t), [guess = 0.0, connect = Flow, bounds = (0.0, Inf)]
  xH(t), [guess = 0.0, bounds = (0.0, 1.0), input = true]
  xD(t), [guess = 0.5, bounds = (0.0, 1.0), input = true]
  xT(t), [guess = 0.5, bounds = (0.0, 1.0), input = true]
  xO(t), [guess = 0.0, bounds = (0.0, 1.0), input = true]
  #p(t), [guess = 1.0e5, bounds = (0.0, Inf), input = true]
  #T(t), [guess = 300.0, bounds = (0.0, Inf), input = true]
end


@doc raw"""
    AFPOutlet_HDTO(;name)
Outlet connector for species H, D, T, O.
# States:
- `G(t)`: [``\mathrm{s^{-1}}``] flow rate
- `xH(t)`: [-] atomic fraction
- `xD(t)`: [-] atomic fraction
- `xT(t)`: [-] atomic fraction
- `xO(t)`: [-] atomic fraction
- `p(t)`: [``\mathrm{Pa}``] pressure
- `T(t)`: [``\mathrm{K}``] temperature
"""
@connector AFPOutlet_HDTO begin
  dummy(t)
  G(t), [guess = 0.0, connect = Flow, bounds = (0.0, Inf)]
  xH(t), [guess = 0.0, bounds = (0.0, 1.0), output = true]
  xD(t), [guess = 0.5, bounds = (0.0, 1.0), output = true]
  xT(t), [guess = 0.5, bounds = (0.0, 1.0), output = true]
  xO(t), [guess = 0.0, bounds = (0.0, 1.0), output = true]
  #p(t), [guess = 1.0e5, bounds = (0.0, Inf), output = true]
  #T(t), [guess = 300.0, bounds = (0.0, Inf), output = true]
end

@doc raw"""
    AFPInlet_OC(;name)
Inlet connector for species H, D, T, O, He, Xe, C, 'I'
# States:
- `G(t)`: [``\mathrm{s^{-1}}``] flow rate
- `xH(t)`:  [-] atomic fraction
- `xD(t)`:  [-] atomic fraction
- `xT(t)`:  [-] atomic fraction
- `xO(t)`:  [-] atomic fraction
- `xHe(t)`: [-] atomic fraction
- `xXe(t)`: [-] atomic fraction
- `xC(t)`:  [-] atomic fraction
- `xI(t)`:  [-] atomic fraction
- `p(t)`: [``\mathrm{Pa}``] pressure
- `T(t)`: [``\mathrm{K}``] temperature
"""
@connector AFPInlet_OC begin
  dummy(t)
  G(t),   [guess = 0.0, connect = Flow, bounds = (0.0, Inf)]
  xH(t),  [guess = 0.0, bounds = (0.0, 1.0), input = true]
  xD(t),  [guess = 0.5, bounds = (0.0, 1.0), input = true]
  xT(t),  [guess = 0.5, bounds = (0.0, 1.0), input = true]
  xHe(t), [guess = 0.0, bounds = (0.0, 1.0), input = true]
  xXe(t), [guess = 0.0, bounds = (0.0, 1.0), input = true]
  xI(t),  [guess = 0.0, bounds = (0.0, 1.0), input = true]
  xO(t),  [guess = 0.0, bounds = (0.0, 1.0), input = true]
  xC(t),  [guess = 0.0, bounds = (0.0, 1.0), input = true]
  #p(t),   [guess = 1.0e5, bounds = (0.0, Inf), input = true]
  #T(t),   [guess = 300.0, bounds = (0.0, Inf), input = true]
end

@doc raw"""
    AFPOutlet_OC(;name)
Outlet connector for species H, D, T, O, He, Xe, C, 'I'
# States:
- `G(t)`: [``\mathrm{s^{-1}}``] flow rate
- `xH(t)`:  [-] atomic fraction
- `xD(t)`:  [-] atomic fraction
- `xT(t)`:  [-] atomic fraction
- `xO(t)`:  [-] atomic fraction
- `xHe(t)`: [-] atomic fraction
- `xXe(t)`: [-] atomic fraction
- `xC(t)`:  [-] atomic fraction
- `xI(t)`:  [-] atomic fraction
- `p(t)`: [``\mathrm{Pa}``] pressure
- `T(t)`: [``\mathrm{K}``] temperature
"""
@connector AFPOutlet_OC begin
  dummy(t)
  G(t),   [guess = 0.0, connect = Flow, bounds = (0.0, Inf)]
  xH(t),  [guess = 0.0, bounds = (0.0, 1.0), output = true]
  xD(t),  [guess = 0.5, bounds = (0.0, 1.0), output = true]
  xT(t),  [guess = 0.5, bounds = (0.0, 1.0), output = true]
  xHe(t), [guess = 0.0, bounds = (0.0, 1.0), output = true]
  xXe(t), [guess = 0.0, bounds = (0.0, 1.0), output = true]
  xI(t),  [guess = 0.0, bounds = (0.0, 1.0), output = true]
  xO(t),  [guess = 0.0, bounds = (0.0, 1.0), output = true]
  xC(t),  [guess = 0.0, bounds = (0.0, 1.0), output = true]
  #p(t),   [guess = 1.0e5, bounds = (0.0, Inf), output = true]
  #T(t),   [guess = 300.0, bounds = (0.0, Inf), output = true]
end

@doc raw"""
    AFPInlet_OCN(;name)
Inlet connector for species H, D, T, 'I', O, C, N.
# States:
- `G(t)`: [``\mathrm{s^{-1}}``] flow rate
- `xH(t)`:  [-] atomic fraction
- `xD(t)`:  [-] atomic fraction
- `xT(t)`:  [-] atomic fraction
- `xI(t)`:  [-] atomic fraction
- `xO(t)`:  [-] atomic fraction
- `xC(t)`:  [-] atomic fraction
- `xN(t)`:  [-] atomic fraction
- `p(t)`: [``\mathrm{Pa}``] pressure
- `T(t)`: [``\mathrm{K}``] temperature
"""
@connector AFPInlet_OCN begin
  dummy(t)
  G(t),   [guess = 0.0, connect = Flow, bounds = (0.0, Inf)]
  xH(t),  [guess = 0.0, bounds = (0.0, 1.0), input = true]
  xD(t),  [guess = 0.5, bounds = (0.0, 1.0), input = true]
  xT(t),  [guess = 0.5, bounds = (0.0, 1.0), input = true]
  xHe(t), [guess = 0.0, bounds = (0.0, 1.0), input = true]
  xXe(t), [guess = 0.0, bounds = (0.0, 1.0), input = true]
  xI(t),  [guess = 0.0, bounds = (0.0, 1.0), input = true]
  xO(t),  [guess = 0.0, bounds = (0.0, 1.0), input = true]
  xC(t),  [guess = 0.0, bounds = (0.0, 1.0), input = true]
  xN(t),  [guess = 0.0, bounds = (0.0, 1.0), input = true]
  #p(t),   [guess = 1.0e5, bounds = (0.0, Inf), input = true]
  #T(t),   [guess = 300.0, bounds = (0.0, Inf), input = true]
end

@doc raw"""
    AFPOutlet_OCN(;name)
Outlet connector for species H, D, T, 'I', O, C, N.
# States:
- `G(t)`: [``\mathrm{s^{-1}}``] flow rate
- `xH(t)`:  [-] atomic fraction
- `xD(t)`:  [-] atomic fraction
- `xT(t)`:  [-] atomic fraction
- `xI(t)`:  [-] atomic fraction
- `xO(t)`:  [-] atomic fraction
- `xC(t)`:  [-] atomic fraction
- `xN(t)`:  [-] atomic fraction
- `p(t)`: [``\mathrm{Pa}``] pressure
- `T(t)`: [``\mathrm{K}``] temperature
"""
@connector AFPOutlet_OCN begin
  dummy(t)
  G(t),   [guess = 0.0, connect = Flow, bounds = (0.0, Inf)]
  xH(t),  [guess = 0.0, bounds = (0.0, 1.0), output = true]
  xD(t),  [guess = 0.5, bounds = (0.0, 1.0), output = true]
  xT(t),  [guess = 0.5, bounds = (0.0, 1.0), output = true]
  xHe(t), [guess = 0.0, bounds = (0.0, 1.0), output = true]
  xXe(t), [guess = 0.0, bounds = (0.0, 1.0), output = true]
  xI(t),  [guess = 0.0, bounds = (0.0, 1.0), output = true]
  xO(t),  [guess = 0.0, bounds = (0.0, 1.0), output = true]
  xC(t),  [guess = 0.0, bounds = (0.0, 1.0), output = true]
  xN(t),  [guess = 0.0, bounds = (0.0, 1.0), output = true]
  #p(t),   [guess = 1.0e5, bounds = (0.0, Inf), output = true]
  #T(t),   [guess = 300.0, bounds = (0.0, Inf), output = true]
end


# ============================================================================ #
#tag inlet_to_outlet2
# ============================================================================ #
@doc raw"""
    inlet_to_outlet2(;name)
Component to pass from AFPOutlet to AFPInlet_O.
Species 'O' is added at the outlet.
# Components:
- `a`: AFPInlet   [inlet]
- `b`: AFPOutlet_O [outlet]
"""
@mtkmodel inlet_to_outlet2 begin
  @components begin
    a = AFPInlet()   # [inlet]
    b = AFPOutlet_O() # [outlet]
  end

  @parameters begin
  end

  @equations begin
    a.xH ~ b.xH
    a.xD ~ b.xD
    a.xT ~ b.xT
    a.xHe ~ b.xHe
    a.xXe ~ b.xXe
    a.xI ~ b.xI
    0.0 ~ b.xO
    0.0 ~ a.G + b.G

    b.dummy ~ 0.0
  end
end


# ============================================================================ #
#tag HDTO_to_outlet
# ============================================================================ #
@doc raw"""
    HDTO_to_outlet(;name)
Component to pass from AFPOutlet_HDTO to AFPInlet.
Species 'O' is removed at the outlet.
Species 'He, Xe, I' are added at the outlet.
# Components:
- `a`: AFPInlet_HDTO [inlet]
- `b`: AFPOutlet    [outlet]
"""
@mtkmodel HDTO_to_outlet begin
  @components begin
    a = AFPInlet_HDTO() # [inlet]
    b = AFPOutlet()    # [outlet]
  end

  @parameters begin
  end

  @equations begin
    a.xH ~ b.xH
    a.xD ~ b.xD
    a.xT ~ b.xT
    0.0  ~ b.xHe
    0.0  ~ b.xXe
    0.0  ~ b.xI 
    0.0 ~ a.G*(a.xH + a.xD + a.xT) + b.G

    b.dummy ~ 0.0
  end
end


# ============================================================================ #
#tag inlet3_to_outletOCN
# ============================================================================ #
@doc raw"""
    inlet3_to_outletOCN(;name)
Component to pass from AFPOutlet_OC to AFPInlet_OCN.
Species 'N' is added at the outlet.
# Components:
- `a`: AFPInlet_OC   [inlet]
- `b`: AFPOutlet_OCN [outlet]
"""
@mtkmodel inlet3_to_outletOCN begin
  @components begin
    a = AFPInlet_OC()   # [inlet]
    b = AFPOutlet_OCN() # [outlet]
  end

  @parameters begin
  end

  @equations begin
    a.xH ~ b.xH
    a.xD ~ b.xD
    a.xT ~ b.xT
    a.xHe ~ b.xHe
    a.xXe ~ b.xXe
    a.xI ~ b.xI
    a.xO ~ b.xO
    a.xC ~ b.xC
    0.0 ~ b.xN
    0.0 ~ a.G + b.G

    b.dummy ~ 0.0
  end
end




# ============================================================================ #
#tag inlet2_to_outletOCN
# ============================================================================ #
# Component to pass from AFPInlet/Outlet to Inlet2/Outlet2 (version 2 also has Oxygen O)
@doc raw"""
    inlet2_to_outletOCN(;name)
Component to pass from AFPOutlet_O to AFPInlet_OCN.
Species 'C, N' are added at the outlet.
# Components:
- `a`: AFPInlet_OC   [inlet]
- `b`: AFPOutlet_OCN [outlet]
"""
@mtkmodel inlet2_to_outletOCN begin
  @components begin
    a = AFPInlet_O()   # [inlet]
    b = AFPOutlet_OCN() # [outlet]
  end

  @parameters begin
  end

  @equations begin
    a.xH ~ b.xH
    a.xD ~ b.xD
    a.xT ~ b.xT
    a.xHe ~ b.xHe
    a.xXe ~ b.xXe
    a.xI ~ b.xI
    a.xO ~ b.xO
    0.0 ~ b.xC
    0.0 ~ b.xN
    0.0 ~ a.G + b.G

    b.dummy ~ 0.0
  end
end




# ============================================================================ #
#tag RemoveOFromFluidPort
# ============================================================================ #
# Component to pass from AFPInlet_O/Outlet2 to Inlet/Outlet (version 2 also has Oxygen O)
@doc raw"""
    RemoveOFromFluidPort(;name)
Component to pass from AFPOutlet_O to AFPInlet.
Species 'O' is removed at the outlet.
# Components:
- `a`: AFPInlet_O  [inlet]
- `b`: AFPOutlet  [outlet]
"""
@mtkmodel RemoveOFromFluidPort begin
  @components begin
    a = AFPInlet_O() # [inlet]
    b = AFPOutlet() # [outlet]
  end

  @parameters begin
  end

  @equations begin
    a.xH *a.G + b.xH *b.G ~ 0.0
    a.xD *a.G + b.xD *b.G ~ 0.0
    a.xT *a.G + b.xT *b.G ~ 0.0
    a.xHe*a.G + b.xHe*b.G ~ 0.0
    a.xXe*a.G + b.xXe*b.G ~ 0.0
    a.xI *a.G + b.xI *b.G ~ 0.0
    0.0 ~ a.G*(1-a.xO) + b.G

    b.dummy ~ 0.0
  end
end




# ============================================================================ #
#tag RemoveInertFromFluidPort
# ============================================================================ #
@doc raw"""
    RemoveInertFromFluidPort(;name)

Connector adapter from `AFPOutlet_O` to `AFPInlet_HDTO`.
Drops the inert species He, Xe, and I from the fluid port, passing only H, D, T, and O.
The outlet flow rate accounts for the removed inert fraction.

# Connectors:
- `a`: `AFPInlet_O` [inlet]
- `b`: `AFPOutlet_HDTO` [outlet]
"""
@mtkmodel RemoveInertFromFluidPort begin
  @components begin
    a = AFPInlet_O() # [inlet]
    b = AFPOutlet_HDTO() # [outlet]
  end

  @parameters begin
  end

  @equations begin
    a.xH*a.G + b.xH*b.G ~ 0.0
    a.xD*a.G + b.xD*b.G ~ 0.0
    a.xT*a.G + b.xT*b.G ~ 0.0
    a.xO*a.G + b.xO*b.G ~ 0.0
    0.0 ~ a.G*(1 - a.xHe - a.xXe - a.xI) + b.G

    b.dummy ~ 0.0
  end
end




# ============================================================================ #
#tag RemoveCFromFluidPort
# ============================================================================ #
@doc raw"""
    RemoveCFromFluidPort(;name)

Connector adapter from `AFPOutlet_OC` to `AFPInlet_O`.
Drops the carbon species C from the fluid port, passing H, D, T, He, Xe, I, and O.
The outlet flow rate accounts for the removed carbon fraction.

# Connectors:
- `a`: `AFPInlet_OC` [inlet]
- `b`: `AFPOutlet_O` [outlet]
"""
@mtkmodel RemoveCFromFluidPort begin
  @components begin
    a = AFPInlet_OC() # [inlet]
    b = AFPOutlet_O() # [outlet]
  end

  @parameters begin
  end

  @equations begin
    a.xH *a.G + b.xH *b.G ~ 0.0
    a.xD *a.G + b.xD *b.G ~ 0.0
    a.xT *a.G + b.xT *b.G ~ 0.0
    a.xHe*a.G + b.xHe*b.G ~ 0.0
    a.xXe*a.G + b.xXe*b.G ~ 0.0
    a.xI *a.G + b.xI *b.G ~ 0.0
    a.xO *a.G + b.xO *b.G ~ 0.0
    0.0 ~ a.G*(1-a.xC) + b.G

    b.dummy ~ 0.0
  end
end


# ============================================================================ #
#tag CuO_ZMS_bed CuOZMS
# ============================================================================ #
@doc raw"""
    CuO_ZMS_bed(;name)

Copper Oxide – Zeolite Molecular Sieve (CuO/ZMS) bed model. Combines an oxidation stage (CuO)
with an adsorption stage (ZMS) to separate hydrogen isotopes from an exhaust gas stream containing
Q species, He, Xe, and inerts. A fraction `eta_ox * eta_zeo` of each Q species and oxygen is
directed to the product outlet `b1`; the remainder, together with inerts, is routed to the reject
outlet `b2`. Helium purge and stoichiometric oxygen regeneration flows are computed internally.

Tritium decay via ``{}^3\mathrm{H} \to {}^3\mathrm{He}`` is included in the `D(NT)` balance.

# States:
- `NH(t)`: [``\mathrm{s^{-1}}``] total number of H atoms in the bed inventory
- `ND(t)`: [``\mathrm{s^{-1}}``] total number of D atoms in the bed inventory
- `NT(t)`: [``\mathrm{s^{-1}}``] total number of T atoms in the bed inventory
- `NHe(t)`: [``\mathrm{s^{-1}}``] total number of He atoms in the bed inventory
- `NXe(t)`: [``\mathrm{s^{-1}}``] total number of Xe atoms in the bed inventory
- `NI(t)`: [``\mathrm{s^{-1}}``] total number of inert atoms in the bed inventory
- `NO(t)`: [``\mathrm{s^{-1}}``] total number of O atoms in the bed inventory

# Connectors:
- `a1`: `AFPInlet_O` [inlet] — exhaust from upstream membranes (retentate), contains Q and He
- `a2`: `AFPInlet_O` [inlet] — O₂ supply for CuO bed regeneration
- `a3`: `AFPInlet_O` [inlet] — He supply for bed purge
- `b1`: `AFPOutlet_O` [outlet] — product stream (Q-enriched)
- `b2`: `AFPOutlet_O` [outlet] — reject stream (inerts and residual Q)

# Parameters:
- `eta_ox`: [-] CuO oxidation efficiency
- `eta_zeo`: [-] ZMS adsorption efficiency
- `τ`: [s] characteristic residence time of the combined bed
- `fHe_purge`: [-] fraction of inlet He flow used for bed purge
- `fO_regen`: [-] stoichiometric ratio of O₂ sent for CuO regeneration (1.0 = stoichiometric)
"""
@mtkmodel CuO_ZMS_bed begin
  @components begin
    a1 = AFPInlet_O()  # [inlet] exhaust from membranes (retentate). Q + PEGs + He
    a2 = AFPInlet_O()  # [inlet] O2 for regeneration
    a3 = AFPInlet_O()  # [inlet] He for regeneration
    b1 = AFPOutlet_O() # [outlet]
    b2 = AFPOutlet_O() # [outlet]
  end

  @parameters begin
    eta_ox = 0.98,     [description = "Oxidation efficiency of CuO [-]", input=true, bounds=(0.0, 1.0)]
    eta_zeo = 0.98,    [description = "Adsorption efficiency of zeolite [-]", input=true, bounds=(0.0, 1.0)]
    τ = 15.0*60*60,    [description = "characteristic time of the whole system [s]", input=true, bounds=(0.0, 100000.0)]
    fHe_purge = 0.001, [description = "fraction of Helium needed for purge as a % of He coming from exhaust [-]", input=true, bounds=(0.0, 1.0)]
    fO_regen = 1.0,    [description = "molar ratio for oxigen sent for regenreation. 1 is stechiometric [-]", input=true, bounds=(0.0, 10000.0)]
  end

  @variables begin
    #
    NH(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    ND(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NT(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NHe(t), [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NXe(t), [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NI(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NO(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
  end

  @equations begin
    D(NH)  ~ a1.xH *a1.G + a2.xH *a2.G + a3.xH *a3.G + b1.xH * b1.G + b2.xH * b2.G
    D(ND)  ~ a1.xD *a1.G + a2.xD *a2.G + a3.xD *a3.G + b1.xD * b1.G + b2.xD * b2.G
    D(NT)  ~ a1.xT *a1.G + a2.xT *a2.G + a3.xT *a3.G + b1.xT * b1.G + b2.xT * b2.G - λT*NT
    D(NHe) ~ a1.xHe*a1.G + a2.xHe*a2.G + a3.xHe*a3.G + b1.xHe* b1.G + b2.xHe* b2.G
    D(NXe) ~ a1.xXe*a1.G + a2.xXe*a2.G + a3.xXe*a3.G + b1.xXe* b1.G + b2.xXe* b2.G
    D(NI)  ~ a1.xI *a1.G + a2.xI *a2.G + a3.xI *a3.G + b1.xI * b1.G + b2.xI * b2.G
    D(NO)  ~ a1.xO *a1.G + a2.xO *a2.G + a3.xO *a3.G + b1.xO * b1.G + b2.xO * b2.G
    #
    # a2.xi ~ ... The cyclinder of type Storage connected here determines the composition. This component defines the flow rate only
    a2.G   ~ fHe_purge*a1.G*a1.xHe
    #
    # a3.xi ~ ... The cyclinder of type Storage connected here determines the composition. This component defines the flow rate only
    a3.G   ~ fO_regen*0.5*a1.G*(a1.xH+a1.xD+a1.xT) # in rapporto stechiometrico
    #
    b1.xH*b1.G  ~ -eta_zeo*eta_ox*NH/τ
    b1.xD*b1.G  ~ -eta_zeo*eta_ox*ND/τ
    b1.xT*b1.G  ~ -eta_zeo*eta_ox*NT/τ
    b1.xO*b1.G  ~ -eta_zeo*eta_ox*NO/τ
    b1.xHe*b1.G ~ -fHe_purge*NHe/τ
    b1.xXe      ~ 0.0
    b1.xI       ~ 0.0
    b1.G        ~ -fHe_purge*NHe/τ -eta_zeo*eta_ox*(NH+ND+NT+NO)/τ
    #
    b2.xH*b2.G  ~ -(1-eta_zeo*eta_ox)*NH/τ
    b2.xD*b2.G  ~ -(1-eta_zeo*eta_ox)*ND/τ
    b2.xT*b2.G  ~ -(1-eta_zeo*eta_ox)*NT/τ
    b2.xO*b2.G  ~ -(1-eta_zeo*eta_ox)*NO/τ
    b2.xHe*b2.G ~ -(1-fHe_purge)*NHe/τ
    b2.xXe*b2.G ~ -NXe/τ
    b2.xI *b2.G ~ -NI/τ
    b2.G        ~ -(1-eta_zeo*eta_ox)*(NH+ND+NT+NO)/τ -(1-fHe_purge)*NHe/τ  -(NXe+NI)/τ

    b1.dummy ~ 0.0
    b2.dummy ~ 0.0
  end
end





# ============================================================================ #
#tag MembraneReactor MR
# ============================================================================ #
@doc raw"""
    MembraneReactor(;name)

Membrane Reactor (MR) model coupling a Water Gas Shift (WGS) reaction stage with a
hydrogen-selective permeation membrane. Inlet `a1` carries a Q₂O-bearing gas stream;
inlet `a2` provides a CO source. The combined efficiency
``\eta_r \cdot \eta_p`` determines the fraction of hydrogen isotopes that react and
permeate through to the permeate outlet `b1`. Unreacted species together with all inerts,
oxygen, and carbon are routed to the retentate outlet `b2`.

Tritium decay via ``{}^3\mathrm{H} \to {}^3\mathrm{He}`` is included in the `D(NT)` balance.

# States:
- `NH(t)`: [``\mathrm{s^{-1}}``] total H inventory
- `ND(t)`: [``\mathrm{s^{-1}}``] total D inventory
- `NT(t)`: [``\mathrm{s^{-1}}``] total T inventory
- `NHe(t)`: [``\mathrm{s^{-1}}``] total He inventory
- `NXe(t)`: [``\mathrm{s^{-1}}``] total Xe inventory
- `NI(t)`: [``\mathrm{s^{-1}}``] total inert inventory
- `NO(t)`: [``\mathrm{s^{-1}}``] total O inventory
- `NC(t)`: [``\mathrm{s^{-1}}``] total C inventory

# Connectors:
- `a1`: `AFPInlet_O` [inlet] — process gas stream (Q₂O)
- `a2`: `AFPInlet_OC` [inlet] — CO source for WGS reaction
- `b1`: `AFPOutlet_O` [outlet] — permeate stream (purified Q)
- `b2`: `AFPOutlet_OC` [outlet] — retentate stream (unreacted species and inerts)

# Parameters:
- `η_r`: [-] WGS reaction efficiency
- `η_p`: [-] membrane permeation efficiency
- `τ`: [s] characteristic residence time
"""
@mtkmodel MembraneReactor begin
  @components begin
    a1 = AFPInlet_O()  # [inlet]
    a2 = AFPInlet_OC()  # [inlet] CO source, connect a cylinder of type Storage (with carbon, Storage3)
    b1 = AFPOutlet_O() # [outlet]
    b2 = AFPOutlet_OC() # [outlet]
  end

  @parameters begin
    η_r = 1.0, [description = "Water Gas Shift (WGS) reaction efficiency [-]", input = true, bounds = (0.0, 1.0)]
    η_p = 1.0, [description = "Permeation efficiency [-]", input = true, bounds = (0.0, 1.0)]
    τ   = 0.01*60*60, [description = "characteristic time [s]", input = true, bounds = (0.0, 100000.0)]
  end

  @variables begin
    #
    NH(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    ND(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NT(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NHe(t), [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NXe(t), [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NI(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NO(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NC(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
  end

  @equations begin
    D(NH)  ~ a1.xH *a1.G + a2.xH *a2.G + b1.xH *b1.G + b2.xH *b2.G
    D(ND)  ~ a1.xD *a1.G + a2.xD *a2.G + b1.xD *b1.G + b2.xD *b2.G
    D(NT)  ~ a1.xT *a1.G + a2.xT *a2.G + b1.xT *b1.G + b2.xT *b2.G - λT*NT
    D(NHe) ~ a1.xHe*a1.G + a2.xHe*a2.G + b1.xHe*b1.G + b2.xHe*b2.G
    D(NXe) ~ a1.xXe*a1.G + a2.xXe*a2.G + b1.xXe*b1.G + b2.xXe*b2.G
    D(NI)  ~ a1.xI *a1.G + a2.xI *a2.G + b1.xI *b1.G + b2.xI *b2.G
    D(NO)  ~ a1.xO *a1.G + a2.xO *a2.G + b1.xO *b1.G + b2.xO *b2.G
    D(NC)  ~ a2.xC *a2.G + b2.xC *b2.G
    #
    a2.G ~ 2*a1.xO*a1.G # At the inlet a1 we get Q2O
    # a2.xi*a2.G ~ ... The cyclinder of type Storage connected here determines the composition. This component defines the flow rate only
    #
    b1.xHe     ~ 0.0
    b1.xXe     ~ 0.0
    b1.xI      ~ 0.0
    b1.xH*b1.G ~ -η_r*η_p*NH/τ
    b1.xD*b1.G ~ -η_r*η_p*ND/τ
    b1.xT*b1.G ~ -η_r*η_p*NT/τ
    b1.xO*b1.G ~ 0.0
    b1.G       ~ -η_r*η_p*(NH+ND+NT)/τ
    #
    b2.xH*b2.G  ~ -(1-η_r*η_p)*NH/τ # tutto ciò che non ha sia reagito che permeato lo trovo al retentato
    b2.xD*b2.G  ~ -(1-η_r*η_p)*ND/τ
    b2.xT*b2.G  ~ -(1-η_r*η_p)*NT/τ
    b2.xHe*b2.G ~ -NHe/τ
    b2.xXe*b2.G ~ -NXe/τ
    b2.xI *b2.G ~ -NI/τ
    b2.xO*b2.G  ~ -NO/τ
    b2.xC*b2.G  ~ -NC/τ
    b2.G        ~ -(1-η_r*η_p)*(NH+ND+NT)/τ -(NHe+NXe+NI+NO+NC)/τ 

    b1.dummy ~ 0.0
    b2.dummy ~ 0.0
  end
end




# ============================================================================ #
#tag WaterDistillation WD
# ============================================================================ #
@doc raw"""
    WaterDistillation(;name)

Water Distillation (WD) column model. Separates a hydrogen-isotope-bearing water stream into a
light-water-enriched top fraction and a tritiated-water-enriched bottom fraction, using
species-specific separation efficiencies and characteristic times. Fresh water input is
specified via the `externalH2O` parameter.

Tritium decay via ``{}^3\mathrm{H} \to {}^3\mathrm{He}`` is included in the `D(NT)` balance.

# States:
- `NH(t)`: [``\mathrm{s^{-1}}``] total H inventory
- `ND(t)`: [``\mathrm{s^{-1}}``] total D inventory
- `NT(t)`: [``\mathrm{s^{-1}}``] total T inventory
- `NO(t)`: [``\mathrm{s^{-1}}``] total O inventory

# Connectors:
- `a1`: `AFPInlet_HDTO` [inlet] — liquid stream from PHTS
- `a2`: `AFPInlet_HDTO` [inlet] — liquid stream from EDS
- `b1`: `AFPOutlet_HDTO` [outlet] — liquid top product returned to PHTS
- `b2`: `AFPOutlet_HDTO` [outlet] — liquid bottom product sent to WDS (electrolyser)

# Parameters:
- `etaH`: [-] fraction of H directed to `b1`
- `etaD`: [-] fraction of D directed to `b1`
- `etaT`: [-] fraction of T directed to `b1`
- `etaO`: [-] fraction of O directed to `b1`
- `τH`: [s] H characteristic residence time
- `τD`: [s] D characteristic residence time
- `τT`: [s] T characteristic residence time
- `τO`: [s] O characteristic residence time
- `externalH2O`: [``\mathrm{kg_{H_2O}\,h^{-1}}``] fresh water feed rate
"""
@mtkmodel WaterDistillation begin
  @components begin
    a1  = AFPInlet_HDTO()  # [inlet][Liq] from PHTS
    a2  = AFPInlet_HDTO()  # [inlet][Liq] from EDS
    b1 = AFPOutlet_HDTO() # [outlet][Liq] to PHTS, top
    b2 = AFPOutlet_HDTO() # [outlet][Liq] to WDS (EL=Electrolyzer), bottom
  end

  @parameters begin
    etaH = 0.992813, [description = "fraction of H inlet flow rate that goes back to PHTS (b1) [-]", input = true, bounds = (0.0, 1.0)]
    etaD = 0.584318, [description = "fraction of D inlet flow rate that goes back to PHTS (b1) [-]", input = true, bounds = (0.0, 1.0)]
    etaT = 0.278140, [description = "fraction of T inlet flow rate that goes back to PHTS (b1) [-]", input = true, bounds = (0.0, 1.0)]
    etaO = 0.992781, [description = "fraction of O inlet flow rate that goes back to PHTS (b1) [-]", input = true, bounds = (0.0, 1.0)]
    τH   = 7.86*60*60,  [description = "H characteristic time [s]", input = true, bounds = (0.0, Inf)]
    τD   = 70.79*60*60, [description = "D characteristic time [s]", input = true, bounds = (0.0, Inf)]
    τT   = 98.07*60*60, [description = "T characteristic time [s]", input = true, bounds = (0.0, Inf)]
    τO   = 7.87*60*60,  [description = "O characteristic time [s]", input = true, bounds = (0.0, Inf)]
    externalH2O = 56.0, [description = "Input of fresh water [kgH2O h-1]", input = true, bounds = (0.0, Inf)]
  end

  @variables begin
    #
    NH(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    ND(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NT(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NO(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
  end

  @equations begin
    D(NH)  ~ a1.xH *a1.G + a2.xH *a2.G + b1.xH *b1.G + b2.xH *b2.G
    D(ND)  ~ a1.xD *a1.G + a2.xD *a2.G + b1.xD *b1.G + b2.xD *b2.G
    D(NT)  ~ a1.xT *a1.G + a2.xT *a2.G + b1.xT *b1.G + b2.xT *b2.G - λT*NT
    D(NO)  ~ a1.xO *a1.G + a2.xO *a2.G + b1.xO *b1.G + b2.xO *b2.G
    #
    a1.G ~ externalH2O * 1000*Na/(60*60*18) * 3 # 3 atoms per molecule of H2O with traces of D according to its natural abundance
    # [moleculesH2O/s] = [kgH2O/h] * (1/3600 h/s) * (1000 g/kg) * (1/18 mol/g) * (Na molecules/mol)
    #
    b1.xT *b1.G ~ -etaT*NT/τT 
    b1.xD *b1.G ~ -etaD*ND/τD
    b1.xH *b1.G ~ -etaH*NH/τH
    b1.xO *b1.G ~ -etaO*NO/τO
    b1.G        ~ -(etaH*NH/τH + etaD*ND/τD + etaT*NT/τT + etaO*NO/τO)
    #
    b2.xT *b2.G ~ -(1-etaT)*NT/τT 
    b2.xD *b2.G ~ -(1-etaD)*ND/τD
    b2.xH *b2.G ~ -(1-etaH)*NH/τH
    b2.xO *b2.G ~ -(1-etaO)*NO/τO
    b2.G        ~ -((1-etaH)*NH/τH + (1-etaD)*ND/τD + (1-etaT)*NT/τT + (1-etaO)*NO/τO)

    b1.dummy ~ 0.0
    b2.dummy ~ 0.0
  end
end
#=
@mtkmodel WaterDistillation begin
  @components begin
    a1 = AFPInlet_O()  # [inlet] from EDS
    a2 = AFPInlet_O()  # [inlet] from CVCS
    b1 = AFPOutlet_O() # [outlet] to LPCE
    b2 = AFPOutlet_O() # [outlet] to 
  end

  @parameters begin
    G_from_CVCS = 100.0, [description = "Water flow rate from the CVCS [kgH2O h-1]", input = true, bounds = (0.0, 10000.0)]
    Gbot  = 1.0, [description = "Water flow rate at the bottom of the column [kgH2O h-1]", input = true, bounds = (0.0, 10000.0)]
    KbotT = 100.0, [description = "Tritium concentration at the bottom of the column [Ci kgH2O-1]", input = true, bounds = (0.0, 10000.0)]
    pHDO = 1.0, [description = "HDO vapor pressure [Pa]", input = true, bounds = (0.0, 1000000.0)]
    pHTO = 1.0, [description = "HTO vapor pressure [Pa]", input = true, bounds = (0.0, 1000000.0)]
    τ  = 15.0*60*60, [description = "characteristic time [s]", input = true, bounds = (0.0, 100000.0)]
  end

  @variables begin
    #
    NH(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    ND(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NT(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NHe(t), [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NXe(t), [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NI(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NO(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    Kbot(t),[guess = 1.0    , description = "H2O fraction that goes to the bottom [-]", bounds = (0.0, Inf)]
    KTbot(t),[guess = 1.0    , description = "T concentration at the bottom [at.T moleculeH2O-1]", bounds = (0.0, Inf)]
  end

  @equations begin
    Kbot ~ Gbot*Na*(1000/(60*60*18)) / (a1.xO*a1.G + a2.xO*a2.G) # fraction of Q2O to the bottom
    KTbot ~ KbotT*(0.001*18/(9.619*3)) # at.T / molecoleH2O
    #
    D(NH)  ~ a1.xH *a1.G + a2.xH *a2.G + b1.xH *b1.G + b2.xH *b2.G
    D(ND)  ~ a1.xD *a1.G + a2.xD *a2.G + b1.xD *b1.G + b2.xD *b2.G
    D(NT)  ~ a1.xT *a1.G + a2.xT *a2.G + b1.xT *b1.G + b2.xT *b2.G
    D(NHe) ~ a1.xHe*a1.G + a2.xHe*a2.G + b1.xHe*b1.G + b2.xHe*b2.G
    D(NXe) ~ a1.xXe*a1.G + a2.xXe*a2.G + b1.xXe*b1.G + b2.xXe*b2.G
    D(NI)  ~ a1.xI *a1.G + a2.xI *a2.G + b1.xI *b1.G + b2.xI *b2.G
    D(NO)  ~ a1.xO *a1.G + a2.xO *a2.G + b1.xO *b1.G + b2.xO *b2.G
    #
    a2.G ~ G_from_CVCS * 1000*Na/(60*60*18) * 3 # 3 atoms per molecule of H2O
    # [moleculesH2O/s] = [kgH2O/h] * (1/3600 h/s) * (1000 g/kg) * (1/18 mol/g) * (Na molecules/mol)
    #
    # b1.xT*b1.G  ~ -KbotT*(0.001*18/(9.619*3)) * Gbot*Na*(1000/(60*60*18))
    # b1.xD*b1.G  ~ -(KbotT / ((a1.xT*a1.G + a2.xT*a2.G)/((0.5*a1.xH+a1.xO)*a1.G + (0.5*a1.xH+a1.xO)*a1.G)))*(0.001*18/(9.619*3))*(pHDO/pHTO) * ND/τ
    # b1.xH*b1.G  ~ -2*Gbot*Na*(1000/(60*60*18)) + (b1.xD + b1.xT)*b1.G # 2 H for every molecule of H2O + 1 H for every HTO + 1 H for every HDO
    # b1.xO*b1.G  ~ -Gbot*Na*(1000/(60*60*18)) + (b1.xD + b1.xT)*b1.G # 1 O for every molecule of H2O + 1 O for every HTO + 1 O for every HDO
    # b1.xHe*b1.G ~ -NHe/τ
    # b1.xXe*b1.G ~ -NXe/τ 
    # b1.xI*b1.G  ~ -NI/τ 
    # b1.G        ~ -KbotT*(0.001*18/(9.619*3)) * Gbot*Na*(1000/(60*60*18)) - 
    #                (KbotT / ((a1.xT*a1.G + a2.xT*a2.G)/((0.5*a1.xH+a1.xO)*a1.G + (0.5*a1.xH+a1.xO)*a1.G)))*(0.001*18/(9.619*3))*(pHDO/pHTO) * ND/τ - 
    #                2*Gbot*Na*(1000/(60*60*18)) + (b1.xD + b1.xT)*b1.G - 
    #                Gbot*Na*(1000/(60*60*18)) + (b1.xD + b1.xT)*b1.G - 
    #                NHe/τ -NXe/τ -NI/τ
    # #
    # b2.xT*b2.G  ~ -NT/τ -(b1.xT*b1.G)
    # b2.xD*b2.G  ~ -ND/τ -(b1.xD*b1.G)
    # b2.xH*b2.G  ~ -NH/τ -(b1.xH*b1.G)
    # b2.xO*b2.G  ~ -NO/τ -(b1.xO*b1.G)
    # b2.xHe*b2.G ~ 0.0
    # b2.xXe*b2.G ~ 0.0
    # b2.xI *b2.G ~ 0.0
    # b2.G        ~ -(NH+NO+NT+ND)/τ -(b1.xH + b1.xO + b1.xT + b1.xD)*b1.G 

    b1.xT*b1.G  ~ -KTbot*(Kbot*NO/τ)
    b1.xD*b1.G  ~ -(pHDO/pHTO)*KTbot*(Kbot*NO/τ)
    b1.xH*b1.G  ~ -(Kbot*NO/τ)*(2 - KTbot*(1 + pHDO/pHTO))
    b1.xO*b1.G  ~ -Kbot*NO/τ
    b1.xHe*b1.G ~ -0.0
    b1.xXe*b1.G ~ -0.0 
    b1.xI*b1.G  ~ -0.0 
    b1.G        ~ -3*Kbot*NO/τ
    #
    b2.xT*b2.G  ~ -NT/τ + KTbot*(Kbot*NO/τ)
    b2.xD*b2.G  ~ -ND/τ + (pHDO/pHTO)*KTbot*(Kbot*NO/τ)
    b2.xH*b2.G  ~ -NH/τ + (Kbot*NO/τ)*(2 - KTbot*(1 + pHDO/pHTO))
    b2.xO*b2.G  ~ -(1-Kbot)*NO/τ
    b2.xHe*b2.G ~ -NHe/τ
    b2.xXe*b2.G ~ -NXe/τ
    b2.xI *b2.G ~ -NI/τ
    b2.G        ~ -(NT+ND+NH)/τ + 2*(Kbot*NO/τ) -(1-Kbot)*NO/τ -NHe/τ -NXe/τ -NI/τ
  end
end
=#



# ============================================================================ #
#tag CECE CECE
# ============================================================================ #
@doc raw"""
    CECE(;name)

Combined Electrolysis and Catalytic Exchange (CECE) system model, integrating a Liquid Phase
Catalytic Exchange (LPCE) column with an electrolyser (EL). Four inlet streams (two liquid
LPCE inlets, one gas LPCE inlet, and one liquid EL inlet) feed the combined system. Three
outlet streams are produced: a gaseous LPCE product, a liquid EL stream to the Isotope
Separation System (ISS), and a separated oxygen stream.

Species-specific efficiencies `etaH`, `etaD`, `etaT`, `etaO` split each species between the
EL outlet (`b1_EL`) and the LPCE outlet (`b_LPCE`). Fresh water input is parameterised by
`externalH2O`.

Tritium decay via ``{}^3\mathrm{H} \to {}^3\mathrm{He}`` is included in the `D(NT)` balance.

# States:
- `NH(t)`: [``\mathrm{s^{-1}}``] total H inventory
- `ND(t)`: [``\mathrm{s^{-1}}``] total D inventory
- `NT(t)`: [``\mathrm{s^{-1}}``] total T inventory
- `NO(t)`: [``\mathrm{s^{-1}}``] total O inventory

# Connectors:
- `a1`: `AFPInlet_HDTO` [inlet] — fresh water to LPCE
- `a2`: `AFPInlet_HDTO` [inlet] — liquid from EDS (WSC) to LPCE
- `a3`: `AFPInlet_HDTO` [inlet] — gas from ISS to LPCE
- `a4`: `AFPInlet_HDTO` [inlet] — liquid from CPS (WD) to EL
- `b_LPCE`: `AFPOutlet_HDTO` [outlet] — gas to H₂ storage
- `b1_EL`: `AFPOutlet_HDTO` [outlet] — liquid to ISS
- `b2_EL`: `AFPOutlet_HDTO` [outlet] — separated oxygen release

# Parameters:
- `etaH`: [-] fraction of H sent to `b1_EL`
- `etaD`: [-] fraction of D sent to `b1_EL`
- `etaT`: [-] fraction of T sent to `b1_EL`
- `etaO`: [-] fraction of O sent to `b1_EL`
- `τH`: [s] H characteristic time
- `τD`: [s] D characteristic time
- `τT`: [s] T characteristic time
- `τO`: [s] O characteristic time
- `externalH2O`: [``\mathrm{kg_{H_2O}\,h^{-1}}``] fresh water feed rate
"""
@mtkmodel CECE begin
  @components begin
    a1     = AFPInlet_HDTO()  # [inlet][Liq] fresh water - LPCE inlet
    a2     = AFPInlet_HDTO()  # [inlet][Liq] from EDS (WSC) - LPCE inlet
    a3     = AFPInlet_HDTO()  # [inlet][Gas] from ISS - LPCE inlet
    a4     = AFPInlet_HDTO()  # [inlet][Liq] from CPS (WD) - EL inlet
    b_LPCE = AFPOutlet_HDTO() # [outlet][Gas] to H2 storage - LPCE outlet
    b1_EL  = AFPOutlet_HDTO() # [outlet][Liq] to ISS - EL outlet
    b2_EL  = AFPOutlet_HDTO() # [outlet][Gas] release of separated oxygen - EL outlet
  end

  @parameters begin
    etaH = 0.0269, [description = "fraction of H inlet flow rate that goes back to PHTS (b1) [-]", input = true, bounds = (0.0, 1.0)]
    etaD = 0.7951, [description = "fraction of D inlet flow rate that goes back to PHTS (b1) [-]", input = true, bounds = (0.0, 1.0)]
    etaT = 1.0   , [description = "fraction of T inlet flow rate that goes back to PHTS (b1) [-]", input = true, bounds = (0.0, 1.0)]
    etaO = 0.0319, [description = "fraction of O inlet flow rate that goes back to PHTS (b1) [-]", input = true, bounds = (0.0, 1.0)]
    τH   = 1.00*60*60, [description = "H characteristic time [s]", input = true, bounds = (0.0, Inf)]
    τD   = 1.00*60*60, [description = "D characteristic time [s]", input = true, bounds = (0.0, Inf)]
    τT   = 1.00*60*60, [description = "T characteristic time [s]", input = true, bounds = (0.0, Inf)]
    τO   = 1.00*60*60, [description = "O characteristic time [s]", input = true, bounds = (0.0, Inf)]
    externalH2O = 56.0, [description = "Input of fresh water [kgH2O h-1]", input = true, bounds = (0.0, Inf)]
  end

  @variables begin
    #
    NH(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    ND(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NT(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NO(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
  end

  @equations begin
    D(NH)  ~ a1.xH *a1.G + a2.xH *a2.G + a3.xH *a3.G + a4.xH *a4.G + b_LPCE.xH *b_LPCE.G + b1_EL.xH *b1_EL.G + b2_EL.xH *b2_EL.G
    D(ND)  ~ a1.xD *a1.G + a2.xD *a2.G + a3.xD *a3.G + a4.xD *a4.G + b_LPCE.xD *b_LPCE.G + b1_EL.xD *b1_EL.G + b2_EL.xD *b2_EL.G
    D(NT)  ~ a1.xT *a1.G + a2.xT *a2.G + a3.xT *a3.G + a4.xT *a4.G + b_LPCE.xT *b_LPCE.G + b1_EL.xT *b1_EL.G + b2_EL.xT *b2_EL.G - λT*NT
    D(NO)  ~ a1.xO *a1.G + a2.xO *a2.G + a3.xO *a3.G + a4.xO *a4.G + b_LPCE.xO *b_LPCE.G + b1_EL.xO *b1_EL.G + b2_EL.xO *b2_EL.G
    #
    a1.G ~ externalH2O * 1000*Na/(60*60*18) * 3 # 3 atoms per molecule of H2O with traces of D according to its natural abundance
    # [moleculesH2O/s] = [kgH2O/h] * (1/3600 h/s) * (1000 g/kg) * (1/18 mol/g) * (Na molecules/mol)
    #
    b2_EL.xH *b2_EL.G ~ -0.0
    b2_EL.xD *b2_EL.G ~ -0.0
    b2_EL.xT *b2_EL.G ~ -0.0 
    b2_EL.xO *b2_EL.G ~ -NO/τO
    b2_EL.G           ~ -NO/τO
    #
    b1_EL.xH *b1_EL.G ~ -etaH*NH/τH
    b1_EL.xD *b1_EL.G ~ -etaD*ND/τD
    b1_EL.xT *b1_EL.G ~ -etaT*NT/τT 
    b1_EL.xO *b1_EL.G ~ -0.0
    b1_EL.G           ~ -etaH*NH/τH - etaD*ND/τD - etaT*NT/τT
    #
    b_LPCE.xH *b_LPCE.G ~ -(1-etaH)*NH/τH
    b_LPCE.xD *b_LPCE.G ~ -(1-etaD)*ND/τD
    b_LPCE.xT *b_LPCE.G ~ -(1-etaT)*NT/τT 
    b_LPCE.xO *b_LPCE.G ~ -0.0
    b_LPCE.G            ~ -(1-etaH)*NH/τH - (1-etaD)*ND/τD - (1-etaT)*NT/τT

    b1_EL.dummy ~ 0.0
    b2_EL.dummy ~ 0.0
    b_LPCE.dummy ~ 0.0
  end
end



# ============================================================================ #
#tag LPCE LPCE
# ============================================================================ #
@doc raw"""
    LPCE(;name)

Liquid Phase Catalytic Exchange (LPCE) column model. Exchanges hydrogen isotopes between a
liquid water stream and a gas stream via a catalytic packed bed. Two inlet streams (liquid
from CPS/WD, and gas ~H₂ from ISS) are processed, yielding a gas product stream to ISS
containing separated Q₂ and Q₂O (`b1`), and a liquid detritiated water stream to EDS (`b2`).

Species-specific efficiencies `etaH`, `etaD`, `etaT`, `etaO` split each species between the
two outlets. Note that `etaO` uses `τH` as the oxygen characteristic time.

Tritium decay via ``{}^3\mathrm{H} \to {}^3\mathrm{He}`` is included in the `D(NT)` balance.

# States:
- `NH(t)`: [``\mathrm{s^{-1}}``] total H inventory
- `ND(t)`: [``\mathrm{s^{-1}}``] total D inventory
- `NT(t)`: [``\mathrm{s^{-1}}``] total T inventory
- `NO(t)`: [``\mathrm{s^{-1}}``] total O inventory

# Connectors:
- `a1`: `AFPInlet_HDTO` [inlet] — liquid stream from CPS (WD)
- `a2`: `AFPInlet_HDTO` [inlet] — gas stream (~H₂) from ISS
- `b1`: `AFPOutlet_HDTO` [outlet] — gas product (Q₂ + Q₂O) to ISS
- `b2`: `AFPOutlet_HDTO` [outlet] — liquid detritiated water to EDS (WSC)

# Parameters:
- `etaH`: [-] fraction of H directed to `b1`
- `etaD`: [-] fraction of D directed to `b1`
- `etaT`: [-] fraction of T directed to `b1`
- `etaO`: [-] fraction of O directed to `b1`
- `τH`: [s] H (and O) characteristic time
- `τD`: [s] D characteristic time
- `τT`: [s] T characteristic time
- `τO`: [s] O characteristic time
"""
@mtkmodel LPCE begin
  @components begin
    a1     = AFPInlet_HDTO()  # [inlet][Liq] from CPS (WD)
    a2     = AFPInlet_HDTO()  # [inlet][Gas] from ISS (~H2)
    b1  = AFPOutlet_HDTO() # [outlet][Gas] to ISS (Q2 + Q2O)
    b2  = AFPOutlet_HDTO() # [outlet][Liq] to EDS (WSC) detritiatited water
  end

  @parameters begin
    etaH = 0.0269, [description = "fraction of H inlet flow rate that goes back to PHTS (b1) [-]", input = true, bounds = (0.0, 1.0)]
    etaD = 0.7951, [description = "fraction of D inlet flow rate that goes back to PHTS (b1) [-]", input = true, bounds = (0.0, 1.0)]
    etaT = 1.0   , [description = "fraction of T inlet flow rate that goes back to PHTS (b1) [-]", input = true, bounds = (0.0, 1.0)]
    etaO = 0.0269, [description = "fraction of O inlet flow rate that goes back to PHTS (b1) [-]", input = true, bounds = (0.0, 1.0)]
    τH   = 1.00*60*60, [description = "H characteristic time [s]", input = true, bounds = (0.0, Inf)]
    τD   = 1.00*60*60, [description = "D characteristic time [s]", input = true, bounds = (0.0, Inf)]
    τT   = 1.00*60*60, [description = "T characteristic time [s]", input = true, bounds = (0.0, Inf)]
    τO   = 1.00*60*60, [description = "O characteristic time [s]", input = true, bounds = (0.0, Inf)]
  end

  @variables begin
    #
    NH(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    ND(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NT(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NO(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
  end

  @equations begin
    D(NH)  ~ a1.xH *a1.G + a2.xH *a2.G + b1.xH *b1.G + b2.xH *b2.G
    D(ND)  ~ a1.xD *a1.G + a2.xD *a2.G + b1.xD *b1.G + b2.xD *b2.G
    D(NT)  ~ a1.xT *a1.G + a2.xT *a2.G + b1.xT *b1.G + b2.xT *b2.G - λT*NT
    D(NO)  ~ a1.xO *a1.G + a2.xO *a2.G + b1.xO *b1.G + b2.xO *b2.G
    #
    # What goes to ISS - separated Q2 + some Q2O
    b1.xH *b1.G ~ -etaH*NH/τH
    b1.xD *b1.G ~ -etaD*ND/τD
    b1.xT *b1.G ~ -etaT*NT/τT 
    b1.xO *b1.G ~ -etaH*NO/τH
    b1.G           ~ -etaH*NH/τH - etaD*ND/τD - etaT*NT/τT - etaH*NO/τH 
    #
    b2.xH *b2.G ~ -(1-etaH)*NH/τH
    b2.xD *b2.G ~ -(1-etaD)*ND/τD
    b2.xT *b2.G ~ -(1-etaT)*NT/τT 
    b2.xO *b2.G ~ -(1-etaH)*NO/τH
    b2.G            ~ -(1-etaH)*NH/τH - (1-etaD)*ND/τD - (1-etaT)*NT/τT - (1-etaH)*NO/τH

    b1.dummy ~ 0.0
    b2.dummy ~ 0.0
  end
end



# ============================================================================ #
#tag PSS Peg Storage System
# ============================================================================ #
@doc raw"""
    PSS(;name)

Permeator–getter Storage System (PSS) model. Acts as a single-port storage volume with a
characteristic emptying time `τ`. All species (H, D, T, He, Xe, I, O) are released through
outlet `b` at a rate proportional to their instantaneous inventory.

Tritium decay via ``{}^3\mathrm{H} \to {}^3\mathrm{He}`` is included in the `D(NT)` balance.

# States:
- `NH(t)`: [``\mathrm{s^{-1}}``] total H inventory
- `ND(t)`: [``\mathrm{s^{-1}}``] total D inventory
- `NT(t)`: [``\mathrm{s^{-1}}``] total T inventory
- `NHe(t)`: [``\mathrm{s^{-1}}``] total He inventory
- `NXe(t)`: [``\mathrm{s^{-1}}``] total Xe inventory
- `NI(t)`: [``\mathrm{s^{-1}}``] total inert inventory
- `NO(t)`: [``\mathrm{s^{-1}}``] total O inventory

# Connectors:
- `a`: `AFPInlet_O` [inlet]
- `b`: `AFPOutlet_O` [outlet]

# Parameters:
- `τ`: [s] characteristic emptying time
"""
@mtkmodel PSS begin
  @components begin
    a = AFPInlet_O()  # [inlet]
    b = AFPOutlet_O() # [outlet]
  end


  @parameters begin
    τ = 60.0*60.0, [description = "characteristic time [s]", input = true, bounds = (0.0, 100000.0)]
  end


  @variables begin
    #
    NH(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    ND(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NT(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NHe(t), [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NXe(t), [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NI(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NO(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
  end


  @equations begin
    D(NH)  ~ a.xH * a.G + b.xH * b.G
    D(ND)  ~ a.xD * a.G + b.xD * b.G
    D(NT)  ~ a.xT * a.G + b.xT * b.G - λT*NT
    D(NHe) ~ a.xHe* a.G + b.xHe* b.G
    D(NXe) ~ a.xXe* a.G + b.xXe* b.G
    D(NI)  ~ a.xI * a.G + b.xI * b.G
    D(NO)  ~ a.xO * a.G + b.xO * b.G
    #
    b.xH * b.G ~ -NH/τ
    b.xD * b.G ~ -ND/τ
    b.xT * b.G ~ -NT/τ
    b.xHe* b.G ~ -NHe/τ
    b.xXe* b.G ~ -NXe/τ
    b.xI * b.G ~ -NI/τ
    b.xO * b.G ~ -NO/τ
    b.G ~ -(NH + ND + NT + NHe + NXe + NI + NO)/τ

    b.dummy ~ 0.0
  end
end



# ============================================================================ #
#tag WetScrubberColumn WSC
# ============================================================================ #
@doc raw"""
    WetScrubberColumn(;name)

Wet Scrubber Column (WSC) model. Absorbs hydrogen isotopes and oxygen from an incoming gas
mixture into a liquid water stream. Three inlet streams are accepted: a gas feed carrying
H, D, T, He, Xe, I, O, C, N from upstream (e.g. Membrane Reactor, room air, glove boxes);
a fresh water makeup stream; and a recycle water stream from WDS. Two outlet streams are
produced: a liquid water product sent to WDS (`b1`) and a gas exhaust to stack (`b2`).

C and N species pass through conservatively to `b2`. The fresh water feed rate is set by
`externalH2O`. Species-specific separation efficiencies and characteristic times govern the
H, D, T, and O splits between `b1` and `b2`.

Tritium decay via ``{}^3\mathrm{H} \to {}^3\mathrm{He}`` is included in the `D(NT)` balance.

# States:
- `NH(t)`: [``\mathrm{s^{-1}}``] total H inventory
- `ND(t)`: [``\mathrm{s^{-1}}``] total D inventory
- `NT(t)`: [``\mathrm{s^{-1}}``] total T inventory
- `NHe(t)`: [``\mathrm{s^{-1}}``] total He inventory
- `NXe(t)`: [``\mathrm{s^{-1}}``] total Xe inventory
- `NI(t)`: [``\mathrm{s^{-1}}``] total inert inventory
- `NO(t)`: [``\mathrm{s^{-1}}``] total O inventory

# Connectors:
- `a1`: `AFPInlet_OCN` [inlet] — gas stream (air, MR exhaust, room/glove box/SVS off-gas)
- `a2`: `AFPInlet_O` [inlet] — fresh water makeup
- `a3`: `AFPInlet_HDTO` [inlet] — recycle water from WDS (LPCE)
- `b1`: `AFPOutlet_O` [outlet] — liquid product to WDS (LPCE)
- `b2`: `AFPOutlet_OCN` [outlet] — gas exhaust to stack

# Parameters:
- `etaH`: [-] fraction of H directed to `b1`
- `etaD`: [-] fraction of D directed to `b1`
- `etaT`: [-] fraction of T directed to `b1`
- `etaO`: [-] fraction of O directed to `b1`
- `τH`: [s] H characteristic time
- `τD`: [s] D characteristic time
- `τT`: [s] T characteristic time
- `τO`: [s] O characteristic time
- `externalH2O`: [``\mathrm{kg_{H_2O}\,h^{-1}}``] fresh water feed rate
"""
@mtkmodel WetScrubberColumn begin
  @components begin
    a1 = AFPInlet_OCN()  # [inlet][Gas] air, from Membrane Reactor (MR) + Rooms + Glove Boxes (GB) + Service Vacuum System (SVS)
    a2 = AFPInlet_O()  # [inlet][Liq] fresh water
    a3 = AFPInlet_HDTO()  # [inlet][Liq] water from WDS (LPCE, caso4)
    b1 = AFPOutlet_O() # [outlet][Liq] to WDS (LPCE)
    b2 = AFPOutlet_OCN() # [outlet][Gas] to stack
  end

  @parameters begin
    etaH = 0.0231, [description = "fraction of H inlet flow rate that goes to WDS (b1) [-]", input = true, bounds = (0.0, 1.0)]
    etaD = 0.0574, [description = "fraction of D inlet flow rate that goes to WDS (b1) [-]", input = true, bounds = (0.0, 1.0)]
    etaT = 0.9528, [description = "fraction of T inlet flow rate that goes to WDS (b1) [-]", input = true, bounds = (0.0, 1.0)]
    etaO = 0.0231, [description = "fraction of O inlet flow rate that goes to WDS (b1) [-]", input = true, bounds = (0.0, 1.0)]
    τH   = 0.92*60*60, [description = "H characteristic time [s]", input = true, bounds = (0.0, Inf)]
    τD   = 1.22*60*60, [description = "D characteristic time [s]", input = true, bounds = (0.0, Inf)]
    τT   = 7.94*60*60, [description = "T characteristic time [s]", input = true, bounds = (0.0, Inf)]
    τO   = 0.92*60*60, [description = "O characteristic time [s]", input = true, bounds = (0.0, Inf)]
    externalH2O = 56.0, [description = "Input of fresh water [kgH2O h-1]", input = true, bounds = (0.0, Inf)]
  end

  @variables begin
    #
    NH(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    ND(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NT(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NHe(t), [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NXe(t), [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NI(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NO(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
  end

  @equations begin
    D(NH)  ~ a1.xH *a1.G + a2.xH *a2.G + a3.xH *a3.G + b1.xH *b1.G + b2.xH *b2.G
    D(ND)  ~ a1.xD *a1.G + a2.xD *a2.G + a3.xD *a3.G + b1.xD *b1.G + b2.xD *b2.G
    D(NT)  ~ a1.xT *a1.G + a2.xT *a2.G + a3.xT *a3.G + b1.xT *b1.G + b2.xT *b2.G - λT*NT
    D(NHe) ~ a1.xHe*a1.G + a2.xHe*a2.G               + b1.xHe*b1.G + b2.xHe*b2.G
    D(NXe) ~ a1.xXe*a1.G + a2.xXe*a2.G               + b1.xXe*b1.G + b2.xXe*b2.G
    D(NI)  ~ a1.xI *a1.G + a2.xI *a2.G               + b1.xI *b1.G + b2.xI *b2.G
    D(NO)  ~ a1.xO *a1.G + a2.xO *a2.G + a3.xO *a3.G + b1.xO *b1.G + b2.xO *b2.G
    0.0  ~ a1.xC *a1.G + b2.xC *b2.G
    0.0  ~ a1.xN *a1.G + b2.xN *b2.G
    #
    a2.G ~ externalH2O * 1000*Na/(60*60*18) * 3 # 3 atoms per molecule of H2O with traces of D according to its natural abundance
    # [moleculesH2O/s] = [kgH2O/h] * (1/3600 h/s) * (1000 g/kg) * (1/18 mol/g) * (Na molecules/mol)
    #
    b1.xT *b1.G ~ -etaT*NT/τT 
    b1.xD *b1.G ~ -etaD*ND/τD
    b1.xH *b1.G ~ -etaH*NH/τH
    b1.xHe*b1.G ~ -0.0
    b1.xXe*b1.G ~ -0.0
    b1.xI *b1.G ~ -0.0
    b1.xO *b1.G ~ -etaO*NO/τO
    b1.G        ~ -(etaH*NH/τH + etaD*ND/τD + etaT*NT/τT + etaO*NO/τO)
    #
    b2.xT *b2.G ~ -(1-etaT)*NT/τT 
    b2.xD *b2.G ~ -(1-etaD)*ND/τD
    b2.xH *b2.G ~ -(1-etaH)*NH/τH
    b2.xHe*b2.G ~ -(a1.xHe*a1.G + a2.xHe*a2.G)
    b2.xXe*b2.G ~ -(a1.xXe*a1.G + a2.xXe*a2.G)
    b2.xI *b2.G ~ -(a1.xI*a1.G + a2.xI*a2.G)
    b2.xO *b2.G ~ -(1-etaO)*NO/τO
    #b2.xC *b2.G ~ -a1.xC*a1.G
    #b2.xN *b2.G ~ -a1.xN*a1.G
    b2.G        ~ -((1-etaH)*NH/τH + (1-etaD)*ND/τD + (1-etaT)*NT/τT + (1-etaO)*NO/τO + a1.G*(a1.xHe+a1.xXe+a1.xI+a1.xC+a1.xN) + a2.G*(a2.xHe+a2.xXe+a2.xI))

    b1.dummy ~ 0.0
    b2.dummy ~ 0.0
  end
end
#=
@mtkmodel WetScrubberColumn begin
  @components begin
    a1 = AFPInlet_O()  # [inlet] from Membrane Reactor (MR) + Rooms + Glove Boxes (GB) + Service Vacuum System (SVS)
    a2 = AFPInlet_O()  # [inlet] from external water source
    a3 = AFPInlet_O()  # [inlet] from LPCE
    b1 = AFPOutlet_O() # [outlet] to Water Distillation (WD)
    b2 = AFPOutlet_O() # [outlet] to stack
  end

  @parameters begin
    DF = 1.0e4, [description = "Detritiation factor [-]", input = true, bounds = (0.0, Inf)]
    τ   = 15.0*60*60, [description = "characteristic time [s]", input = true, bounds = (0.0, Inf)]
    pHDO = 1.0, [description = "HDO vapor pressure [Pa]", input = true, bounds = (0.0, Inf)]
    pHTO = 1.0, [description = "HTO vapor pressure [Pa]", input = true, bounds = (0.0, Inf)]
    pH2O = 1.0, [description = "H2O vapor pressure [Pa]", input = true, bounds = (0.0, Inf)]
    externalH2O = 56.0, [description = "Input of fresh water [kgH2O h-1]", input = true, bounds = (0.0, Inf)]
  end

  @variables begin
    #
    NH(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    ND(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NT(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NHe(t), [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NXe(t), [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NI(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NO(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
  end

  @equations begin
    D(NH)  ~ a1.xH *a1.G + a2.xH *a2.G + a3.xH *a3.G + b1.xH *b1.G + b2.xH *b2.G
    D(ND)  ~ a1.xD *a1.G + a2.xD *a2.G + a3.xD *a3.G + b1.xD *b1.G + b2.xD *b2.G
    D(NT)  ~ a1.xT *a1.G + a2.xT *a2.G + a3.xT *a3.G + b1.xT *b1.G + b2.xT *b2.G
    D(NHe) ~ a1.xHe*a1.G + a2.xHe*a2.G + a3.xHe*a3.G + b1.xHe*b1.G + b2.xHe*b2.G
    D(NXe) ~ a1.xXe*a1.G + a2.xXe*a2.G + a3.xXe*a3.G + b1.xXe*b1.G + b2.xXe*b2.G
    D(NI)  ~ a1.xI *a1.G + a2.xI *a2.G + a3.xI *a3.G + b1.xI *b1.G + b2.xI *b2.G
    D(NO)  ~ a1.xO *a1.G + a2.xO *a2.G + a3.xO *a3.G + b1.xO *b1.G + b2.xO *b2.G
    #
    a2.G ~ externalH2O * 1000*Na/(60*60*18) * 3 # molecules of H2O with traces of D according to its natural abundance
    # [moleculesH2O/s] = [kgH2O/h] * (1/3600 h/s) * (1000 g/kg) * (1/18 mol/g) * (Na molecules/mol)
    #
    b2.xT *b2.G ~ -(1/DF)*NT/τ
    b2.xD *b2.G ~ -(1/((pHDO/pHTO)*DF))*ND/τ
    b2.xH *b2.G ~ -(1/((2*pH2O/pHTO)*DF))*NH/τ -(1/DF)*NT/τ -(1/((pHDO/pHTO)*DF))*ND/τ # H from "DF_H2O" + 1H from HTO + 1H from HDO
    b2.xHe*b2.G ~ -NHe/τ
    b2.xXe*b2.G ~ -NXe/τ
    b2.xI *b2.G ~ -NI/τ
    b2.xO *b2.G ~ -0.5*(1/((2*pH2O/pHTO)*DF))*NH/τ -(1/DF)*NT/τ -(1/((pHDO/pHTO)*DF))*ND/τ # 1 O for every 2 H in H2O + 1 O from HTO + 1O from HDO
    b2.G        ~ -(1/DF)*NT/τ - 
                   (1/((pHDO/pHTO)*DF))*ND/τ - 
                   (1/((2*pH2O/pHTO)*DF))*NH/τ -(1/DF)*NT/τ -(1/((pHDO/pHTO)*DF))*ND/τ - 
                   (NHe+NXe+NI)/τ -
                   0.5*(1/((2*pH2O/pHTO)*DF))*NH/τ -(1/DF)*NT/τ -(1/((pHDO/pHTO)*DF))*ND/τ
    #
    b1.xT *b1.G ~ -(1 - 1/DF)*NT/τ
    b1.xD *b1.G ~ -(1 - (1/((pHDO/pHTO)*DF)))*ND/τ
    b1.xH *b1.G ~ -(1 - (1/((2*pH2O/pHTO)*DF)))*NH/τ +(1/DF)*NT/τ +(1/((pHDO/pHTO)*DF))*ND/τ
    b1.xHe      ~  0.0
    b1.xXe      ~  0.0
    b1.xI       ~  0.0
    b1.xO *b1.G ~ -NO/τ + 0.5*(1/((2*pH2O/pHTO)*DF))*NH/τ +(1/DF)*NT/τ +(1/((pHDO/pHTO)*DF))*ND/τ
    b1.G        ~ -(1 - 1/DF)*NT/τ  -
                   (1 - (1/((pHDO/pHTO)*DF)))*ND/τ   - 
                   (1 - (1/((2*pH2O/pHTO)*DF)))*NH/τ +(1/DF)*NT/τ +(1/((pHDO/pHTO)*DF))*ND/τ - 
                   NO/τ + 0.5*(1/((2*pH2O/pHTO)*DF))*NH/τ +(1/DF)*NT/τ +(1/((pHDO/pHTO)*DF))*ND/τ
  end
end
=#




# ============================================================================ #
#tag CryogenicDistillation CD
# ============================================================================ #
@doc raw"""
    CryogenicDistillation(;name)

Cryogenic Distillation (CD) column model. Separates a hydrogen-isotope mixture into a
tritium-enriched product stream and a deuterium/protium-enriched recycle stream. Two inlet
streams are accepted: a liquid Q₂ feed from WDS (CECE) and a gas stream from the Palladium
Reactor (PR, e.g. TCAP). Species-specific efficiencies `etaH`, `etaD`, `etaT` and
characteristic times `τH`, `τD`, `τT` govern the split between the product outlet `b1`
(to tritium storage) and the recycle outlet `b2` (to WDS/LPCE in CECE). Oxygen is conserved
and exits entirely through `b2`.

Tritium decay via ``{}^3\mathrm{H} \to {}^3\mathrm{He}`` is included in the `D(NT)` balance.

# States:
- `NH(t)`: [``\mathrm{s^{-1}}``] total H inventory
- `ND(t)`: [``\mathrm{s^{-1}}``] total D inventory
- `NT(t)`: [``\mathrm{s^{-1}}``] total T inventory
- `NO(t)`: [``\mathrm{s^{-1}}``] total O inventory

# Connectors:
- `a1`: `AFPInlet_HDTO` [inlet] — liquid feed from WDS (CECE)
- `a2`: `AFPInlet_HDTO` [inlet] — gas feed from PR (TCAP)
- `b1`: `AFPOutlet_HDTO` [outlet] — product stream to tritium storage
- `b2`: `AFPOutlet_HDTO` [outlet] — recycle stream to WDS (LPCE in CECE)

# Parameters:
- `etaH`: [-] fraction of H directed to `b1`
- `etaD`: [-] fraction of D directed to `b1`
- `etaT`: [-] fraction of T directed to `b1`
- `τH`: [s] H characteristic time
- `τD`: [s] D characteristic time
- `τT`: [s] T characteristic time
"""
@mtkmodel CryogenicDistillation begin
  @components begin
    a1 = AFPInlet_HDTO()  # [inlet][Liq] from WDS (CECE)
    a2 = AFPInlet_HDTO()  # [inlet][Gas] from PR (TCAP?)
    b1 = AFPOutlet_HDTO() # [outlet][Gas] to Storage
    b2 = AFPOutlet_HDTO() # [outlet][Gas] to WDS (LPCE in CECE)
  end

  @parameters begin
    etaH = 0.0000, [description = "fraction of H inlet flow rate that goes to Storage (b1) [-]", input = true, bounds = (0.0, 1.0)]
    etaD = 0.0954, [description = "fraction of D inlet flow rate that goes to Storage (b1) [-]", input = true, bounds = (0.0, 1.0)]
    etaT = 0.9999, [description = "fraction of T inlet flow rate that goes to Storage (b1) [-]", input = true, bounds = (0.0, 1.0)]
    τH   = 7.00*60*60, [description = "H characteristic time [s]", input = true, bounds = (0.0, Inf)]
    τD   = 7.00*60*60, [description = "D characteristic time [s]", input = true, bounds = (0.0, Inf)]
    τT   = 7.00*60*60, [description = "T characteristic time [s]", input = true, bounds = (0.0, Inf)]
  end

  @variables begin
    #
    NH(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    ND(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NT(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NO(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
  end

  @equations begin
    D(NH)  ~ a1.xH *a1.G + a2.xH *a2.G + b1.xH *b1.G + b2.xH *b2.G
    D(ND)  ~ a1.xD *a1.G + a2.xD *a2.G + b1.xD *b1.G + b2.xD *b2.G
    D(NT)  ~ a1.xT *a1.G + a2.xT *a2.G + b1.xT *b1.G + b2.xT *b2.G - λT*NT
    D(NO)  ~ a1.xO *a1.G + a2.xO *a2.G + b1.xO *b1.G + b2.xO *b2.G
    #
    b1.xH *b1.G ~ -etaH*NH/τH
    b1.xD *b1.G ~ -etaD*ND/τD
    b1.xT *b1.G ~ -etaT*NT/τT 
    b1.xO *b1.G ~ -0.0
    b1.G        ~ -etaH*NH/τH - etaD*ND/τD - etaT*NT/τT
    #
    b2.xH *b2.G ~ -(1-etaH)*NH/τH
    b2.xD *b2.G ~ -(1-etaD)*ND/τD
    b2.xT *b2.G ~ -(1-etaT)*NT/τT 
    b2.xO *b2.G ~ -NO/τH
    b2.G        ~ -(1-etaH)*NH/τH - (1-etaD)*ND/τD - (1-etaT)*NT/τT - NO/τH

    b1.dummy ~ 0.0
    b2.dummy ~ 0.0
  end
end



#=
# ============================================================================ #
#tag SimpleValve 
# ============================================================================ #
# Simple valve operated by pressure difference
@mtkmodel SimpleValve begin
  @components begin
    a = AFPInlet()  # [inlet]
    b = AFPOutlet() # [outlet]
  end

  @parameters begin
    Cd = 0.1
    A = 0.001
  end

  @equations begin
    #Δp ~ a.p - b.p
    #G  ~ a.G
    a.xH ~ b.xH
    a.xD ~ b.xD
    a.xT ~ b.xT
    a.xHe ~ b.xHe
    a.xXe ~ b.xXe
    a.xI ~ b.xI
    0.0 ~ a.G + b.G
    #
    a.G ~ ifelse(a.p - b.p > 0.0, regRoot(2 * (a.p - b.p) * (a.p / (Rg * a.T))) * Cd * A, 0.0)
  end
end
=#

#=
# ============================================================================ #
#tag Actuatedvalve 
# ============================================================================ #
# Simple valve operated by pressure difference
@mtkmodel ActuatedValve begin
  @components begin
    a = AFPInlet()  # [inlet]
    b = AFPOutlet() # [outlet]
    a_ctrl = RealInput() # [input]
  end

  @variables begin
  end

  @equations begin
    a.xH  ~ b.xH
    a.xD  ~ b.xD
    a.xT  ~ b.xT
    a.xHe ~ b.xHe
    a.xXe ~ b.xXe
    a.xI  ~ b.xI
    0.0 ~ a.G + b.G
    #
    a.G ~ ifelse(a.p - b.p >= 0.0, a_ctrl.u, 0.0)
  end
end
=#


# ============================================================================ #
#tag FixedActuatedValve 
# ============================================================================ #
@doc raw"""
    FixedActuatedValve(;name)

Idealised actuated valve whose flow rate is set directly by an external control signal.
The composition is passed through without modification. No pressure-drop model is included.

# Connectors:
- `a`: `AFPInlet` [inlet]
- `b`: `AFPOutlet` [outlet]
- `a_ctrl`: `RealInput` [input] — prescribed flow rate ``G`` [``\mathrm{s^{-1}}``]
"""
@mtkmodel FixedActuatedValve begin
  @components begin
    a = AFPInlet()  # [inlet]
    b = AFPOutlet() # [outlet]
    a_ctrl = RealInput()  # [input]
  end

  @equations begin
    a.xH ~ b.xH
    a.xD ~ b.xD
    a.xT ~ b.xT
    a.xHe ~ b.xHe
    a.xXe ~ b.xXe
    a.xI ~ b.xI
    0.0 ~ a.G + b.G
    #
    a.G ~ a_ctrl.u

    b.dummy ~ 0.0
  end
end




# ============================================================================ #
#tag Vessel3I2O 
# ============================================================================ #
@doc raw"""
    Vessel3I2O(;name)

Ideal gas storage vessel with three inlet ports and two outlet ports, tracking species
H, D, T, He, Xe, and I. The vessel pressure is computed from the ideal gas law. Both
outlets carry the same composition, set by the instantaneous mole fractions in the vessel.
Control outputs provide the operating-pressure setpoint and the measured deuterium and
tritium inventories for use with external PID controllers.

Tritium decay via ``{}^3\mathrm{H} \to {}^3\mathrm{He}`` is included in the `D(NT)` balance.

# States:
- `NH(t)`: [``\mathrm{mol}``] monoatomic H inventory
- `ND(t)`: [``\mathrm{mol}``] monoatomic D inventory
- `NT(t)`: [``\mathrm{mol}``] monoatomic T inventory
- `NHe(t)`: [``\mathrm{mol}``] monoatomic He inventory
- `NXe(t)`: [``\mathrm{mol}``] monoatomic Xe inventory
- `NI(t)`: [``\mathrm{mol}``] monoatomic inert inventory
- `p(t)`: [``\mathrm{Pa}``] vessel pressure

# Connectors:
- `a1`, `a2`, `a3`: `AFPInlet` [inlet]
- `b1`, `b2`: `AFPOutlet` [outlet]
- `ND_ref`: `RealOutput` [output] — D inventory setpoint at operating pressure
- `NT_ref`: `RealOutput` [output] — T inventory setpoint at operating pressure
- `ND_meas`: `RealOutput` [output] — measured D inventory
- `NT_meas`: `RealOutput` [output] — measured T inventory

# Parameters:
- `T`: [K] vessel temperature
- `V`: [``\mathrm{m^3}``] vessel volume
- `p_op`: [``\mathrm{Pa}``] operating pressure (absolute)
"""
@mtkmodel Vessel3I2O begin
  @components begin
    b1 = AFPOutlet() # [outlet]
    b2 = AFPOutlet() # [outlet]
    a1 = AFPInlet() # [inlet]
    a2 = AFPInlet() # [inlet]
    a3 = AFPInlet() # [inlet]
    ND_ref = RealOutput() # [output] pressure reference
    NT_ref = RealOutput() # [output] pressure reference
    ND_meas = RealOutput() # [output] pressure measurement
    NT_meas = RealOutput() # [output] pressure measurement
    
  end

  @parameters begin
    T = 298.15 # [K] vessel temperature
    V = 1.0 # [m3] vessel volume
    p_op = 0.9e5 # [Pa] operational pressure # WARNING: absolute pressure!
  end

  @variables begin
    #
    NH(t),  [guess = 1.0e-20, description = "moles (monoatomic) [mol]", bounds = (0.0, Inf)]
    ND(t),  [guess = 1.0e-20, description = "moles (monoatomic) [mol]", bounds = (0.0, Inf)]
    NT(t),  [guess = 1.0e-20, description = "moles (monoatomic) [mol]", bounds = (0.0, Inf)]
    NHe(t), [guess = 1.0e-20, description = "moles (monoatomic) [mol]", bounds = (0.0, Inf)]
    NXe(t), [guess = 1.0e-20, description = "moles (monoatomic) [mol]", bounds = (0.0, Inf)]
    NI(t),  [guess = 1.0e-20, description = "moles (monoatomic) [mol]", bounds = (0.0, Inf)]
    p(t),   [guess = 1.0e5, bounds = (0.0, Inf)]
  end

  @equations begin
    # equations for block variables
    D(NH)  ~ a1.xH * a1.G  + a2.xH * a2.G  + a3.xH * a3.G  + b1.xH * b1.G  + b2.xH * b2.G
    D(ND)  ~ a1.xD * a1.G  + a2.xD * a2.G  + a3.xD * a3.G  + b1.xD * b1.G  + b2.xD * b2.G
    D(NT)  ~ a1.xT * a1.G  + a2.xT * a2.G  + a3.xT * a3.G  + b1.xT * b1.G  + b2.xT * b2.G - λT*NT
    D(NHe) ~ a1.xHe * a1.G + a2.xHe * a2.G + a3.xHe * a3.G + b1.xHe * b1.G + b2.xHe * b2.G
    D(NXe) ~ a1.xXe * a1.G + a2.xXe * a2.G + a3.xXe * a3.G + b1.xXe * b1.G + b2.xXe * b2.G
    D(NI)  ~ a1.xI * a1.G  + a2.xI * a2.G  + a3.xI * a3.G  + b1.xI * b1.G  + b2.xI * b2.G
    p ~ (NH/2 + ND/2 + NT/2 + NHe + NXe + NI) * Rg * T / V
    #rho ~ 0.001 * (NH*1 + ND*2 + NT*3 + NHe*4 + NXe*131.293 + NI*40.0)

    # equations for state outlet bs1
    b1.xH  ~ NH  / (NH + ND + NT + NHe + NXe + NI)
    b1.xD  ~ ND  / (NH + ND + NT + NHe + NXe + NI)
    b1.xT  ~ NT  / (NH + ND + NT + NHe + NXe + NI)
    b1.xHe ~ NHe / (NH + ND + NT + NHe + NXe + NI)
    b1.xXe ~ NXe / (NH + ND + NT + NHe + NXe + NI)
    b1.xI  ~ NI  / (NH + ND + NT + NHe + NXe + NI)
    #b1.p   ~ p
    #b1.T   ~ T

    # equation for state outlet bs2
    b2.xH  ~ NH  / (NH + ND + NT + NHe + NXe + NI)
    b2.xD  ~ ND  / (NH + ND + NT + NHe + NXe + NI)
    b2.xT  ~ NT  / (NH + ND + NT + NHe + NXe + NI)
    b2.xHe ~ NHe / (NH + ND + NT + NHe + NXe + NI)
    b2.xXe ~ NXe / (NH + ND + NT + NHe + NXe + NI)
    b2.xI  ~ NI  / (NH + ND + NT + NHe + NXe + NI)
    #b2.p   ~ p
    #b2.T   ~ T

    # Control equations
    ND_ref.u ~ p_op*V/(Rg*T) - (NH/2 + NHe + NXe + NI)
    NT_ref.u ~ p_op*V/(Rg*T) - (NH/2 + NHe + NXe + NI)
    ND_meas.u ~ ND
    NT_meas.u ~ NT
    
    # Temporary
    #a2.xH  ~ 0.0
    #a2.xD  ~ 0.5
    #a2.xT  ~ 0.5
    #a2.xHe ~ 0.0
    #a2.xI  ~ 0.0
    #a2.G   ~ 0.0
    #a3.xH  ~ 0.0
    #a3.xD  ~ 0.5
    #a3.xT  ~ 0.5
    #a3.xHe ~ 0.0
    #a3.xI  ~ 0.0
    #a3.G   ~ 0.0

    b1.dummy ~ 0.0
    b2.dummy ~ 0.0
  end
end



# ============================================================================ #
#tag Vessel1I1O 
# ============================================================================ #
@doc raw"""
    Vessel1I1O(;name)

Ideal gas storage vessel with one inlet port and one outlet port, tracking species
H, D, T, He, Xe, and I. The vessel pressure is computed from the ideal gas law. The outlet
composition is set by the instantaneous mole fractions in the vessel. Control outputs provide
the operating-pressure setpoint and the measured pressure for use with an external PID
controller.

Tritium decay via ``{}^3\mathrm{H} \to {}^3\mathrm{He}`` is included in the `D(NT)` balance.

# States:
- `NH(t)`: [``\mathrm{mol}``] monoatomic H inventory
- `ND(t)`: [``\mathrm{mol}``] monoatomic D inventory
- `NT(t)`: [``\mathrm{mol}``] monoatomic T inventory
- `NHe(t)`: [``\mathrm{mol}``] monoatomic He inventory
- `NXe(t)`: [``\mathrm{mol}``] monoatomic Xe inventory
- `NI(t)`: [``\mathrm{mol}``] monoatomic inert inventory
- `p(t)`: [``\mathrm{Pa}``] vessel pressure

# Connectors:
- `a`: `AFPInlet` [inlet]
- `b`: `AFPOutlet` [outlet]
- `b_p_ref`: `RealOutput` [output] — pressure setpoint
- `b_p_meas`: `RealOutput` [output] — measured pressure

# Parameters:
- `T`: [K] vessel temperature
- `V`: [``\mathrm{m^3}``] vessel volume
- `p_op`: [``\mathrm{Pa}``] operating pressure (absolute)
"""
@mtkmodel Vessel1I1O begin
  @components begin
    b = AFPOutlet() # [outlet]
    a = AFPInlet() # [inlet]
    b_p_ref = RealOutput() # [output] pressure reference
    b_p_meas = RealOutput() # [output] pressure measurement
  end

  @parameters begin
    T = 298.15 # [K] vessel temperature
    V = 1.0 # [m3] vessel volume
    p_op = 0.8e5 # [Pa] operational pressure # WARNING: absolute pressure!
  end

  @variables begin
    #
    NH(t),  [guess = 1.0e-20, description = "moles (monoatomic) [mol]", bounds = (0.0, Inf)]
    ND(t),  [guess = 1.0e-20, description = "moles (monoatomic) [mol]", bounds = (0.0, Inf)]
    NT(t),  [guess = 1.0e-20, description = "moles (monoatomic) [mol]", bounds = (0.0, Inf)]
    NHe(t), [guess = 1.0e-20, description = "moles (monoatomic) [mol]", bounds = (0.0, Inf)]
    NXe(t), [guess = 1.0e-20, description = "moles (monoatomic) [mol]", bounds = (0.0, Inf)]
    NI(t),  [guess = 1.0e-20, description = "moles (monoatomic) [mol]", bounds = (0.0, Inf)]
    p(t),   [guess = 1.0e5, bounds = (0.0, Inf)]
  end

  @equations begin
    # equations for block variables
    D(NH)  ~ a.xH  * a.G + b.xH  * b.G
    D(ND)  ~ a.xD  * a.G + b.xD  * b.G
    D(NT)  ~ a.xT  * a.G + b.xT  * b.G - λT*NT
    D(NHe) ~ a.xHe * a.G + b.xHe * b.G
    D(NXe) ~ a.xHe * a.G + b.xHe * b.G
    D(NI)  ~ a.xI  * a.G + b.xI  * b.G
    p ~ (NH/2 + ND/2 + NT/2 + NHe + NXe + NI) * Rg * T / V
    #rho ~ 0.001 * (nH * 1 + nD * 2 + nT * 3 + nHe * 4 + nI * 40.0)

    # equations for state outlet bs1
    b.xH  ~ NH  / (NH + ND + NT + NHe + NXe + NI)
    b.xD  ~ ND  / (NH + ND + NT + NHe + NXe + NI)
    b.xT  ~ NT  / (NH + ND + NT + NHe + NXe + NI)
    b.xHe ~ NHe / (NH + ND + NT + NHe + NXe + NI)
    b.xXe ~ NXe / (NH + ND + NT + NHe + NXe + NI)
    b.xI  ~ NI  / (NH + ND + NT + NHe + NXe + NI)
    #b.p ~ p
    #b.T ~ T

    # Control equations
    b_p_ref.u ~ p_op
    b_p_meas.u ~ p

    b.dummy ~ 0.0
  end
end



# ============================================================================ #
#tag Vessel2I1O 
# ============================================================================ #
@doc raw"""
    Vessel2I1O(;name)

Ideal gas storage vessel with two inlet ports and one outlet port, tracking species
H, D, T, He, Xe, and I. The vessel pressure is computed from the ideal gas law. The outlet
composition is set by the instantaneous mole fractions in the vessel. Control outputs provide
the operating-pressure setpoint and the measured pressure for use with an external PID
controller.

Tritium decay via ``{}^3\mathrm{H} \to {}^3\mathrm{He}`` is included in the `D(NT)` balance.

# States:
- `NH(t)`: [``\mathrm{mol}``] monoatomic H inventory
- `ND(t)`: [``\mathrm{mol}``] monoatomic D inventory
- `NT(t)`: [``\mathrm{mol}``] monoatomic T inventory
- `NHe(t)`: [``\mathrm{mol}``] monoatomic He inventory
- `NXe(t)`: [``\mathrm{mol}``] monoatomic Xe inventory
- `NI(t)`: [``\mathrm{mol}``] monoatomic inert inventory
- `p(t)`: [``\mathrm{Pa}``] vessel pressure

# Connectors:
- `a1`, `a2`: `AFPInlet` [inlet]
- `b`: `AFPOutlet` [outlet]
- `b_p_ref`: `RealOutput` [output] — pressure setpoint
- `b_p_meas`: `RealOutput` [output] — measured pressure

# Parameters:
- `T`: [K] vessel temperature
- `V`: [``\mathrm{m^3}``] vessel volume
- `p_op`: [``\mathrm{Pa}``] operating pressure (absolute)
"""
@mtkmodel Vessel2I1O begin
  @components begin
    b = AFPOutlet() # [outlet]
    a1 = AFPInlet() # [inlet]
    a2 = AFPInlet() # [inlet]
    b_p_ref = RealOutput() # [output] pressure reference
    b_p_meas = RealOutput() # [output] pressure measurement
  end

  @parameters begin
    T = 298.15 # [K] vessel temperature
    V = 1.0 # [m3] vessel volume
    p_op = 0.8e5 # [Pa] operational pressure # WARNING: absolute pressure!
  end

  @variables begin
    #
    NH(t),  [guess = 1.0e-20, description = "moles (monoatomic) [mol]", bounds = (0.0, Inf)]
    ND(t),  [guess = 1.0e-20, description = "moles (monoatomic) [mol]", bounds = (0.0, Inf)]
    NT(t),  [guess = 1.0e-20, description = "moles (monoatomic) [mol]", bounds = (0.0, Inf)]
    NHe(t), [guess = 1.0e-20, description = "moles (monoatomic) [mol]", bounds = (0.0, Inf)]
    NXe(t), [guess = 1.0e-20, description = "moles (monoatomic) [mol]", bounds = (0.0, Inf)]
    NI(t),  [guess = 1.0e-20, description = "moles (monoatomic) [mol]", bounds = (0.0, Inf)]
    p(t),   [guess = 1.0e5, bounds = (0.0, Inf)]
  end

  @equations begin
    # equations for block variables
    D(NH)  ~ a1.xH  * a1.G + a2.xH  * a2.G + b.xH  * b.G
    D(ND)  ~ a1.xD  * a1.G + a2.xD  * a2.G + b.xD  * b.G
    D(NT)  ~ a1.xT  * a1.G + a2.xT  * a2.G + b.xT  * b.G - λT*NT
    D(NHe) ~ a1.xHe * a1.G + a2.xHe * a2.G + b.xHe * b.G
    D(NXe) ~ a1.xXe * a1.G + a2.xXe * a2.G + b.xXe * b.G
    D(NI)  ~ a1.xI  * a1.G + a2.xI  * a2.G + b.xI  * b.G
    p ~ (NH/2 + ND/2 + NT/2 + NHe + NXe + NI) * Rg * T / V
    #rho ~ 0.001 * (nH * 1 + nD * 2 + nT * 3 + nHe * 4 + nI * 40.0)

    # equations for state outlet bs1
    b.xH  ~ NH  / (NH + ND + NT + NHe + NXe + NI)
    b.xD  ~ ND  / (NH + ND + NT + NHe + NXe + NI)
    b.xT  ~ NT  / (NH + ND + NT + NHe + NXe + NI)
    b.xHe ~ NHe / (NH + ND + NT + NHe + NXe + NI)
    b.xXe ~ NXe / (NH + ND + NT + NHe + NXe + NI)
    b.xI  ~ NI  / (NH + ND + NT + NHe + NXe + NI)
    #b.p ~ p
    #b.T ~ T
    #b.rho ~ 0.001 * (nH * 1 + nD * 2 + nT * 3 + nHe * 4 + nI * 40.0)

    # Control equations
    b_p_ref.u ~ p_op
    b_p_meas.u ~ p

    b.dummy ~ 0.0
  end
end



# ============================================================================ #
#tag fixedBC 
# ============================================================================ #
@doc raw"""
    fixedBC(;name)

Fixed boundary condition block. Prescribes a flow rate and composition on the connected
`AFPInlet` port via freely settable variable states. Intended for use at system boundaries
where an external time-varying signal drives the inlet condition.

# States:
- `G(t)`: [``\mathrm{s^{-1}}``] flow rate (output)
- `xH(t)`: [-] H atomic fraction (output)
- `xD(t)`: [-] D atomic fraction (output)
- `xT(t)`: [-] T atomic fraction (output)
- `xHe(t)`: [-] He atomic fraction (output)
- `xXe(t)`: [-] Xe atomic fraction (output)
- `xI(t)`: [-] inert atomic fraction (output)

# Connectors:
- `a`: `AFPInlet` [inlet]

# Parameters:
- `p`: [``\mathrm{Pa}``] operating pressure
- `T`: [K] temperature
"""
@mtkmodel fixedBC begin
  @components begin
    a = AFPInlet()
  end

  @parameters begin
    p = 0.7e5, [description = "Operating pressure [Pa]", input = true, bounds = (0.0, Inf)]
    T = 300.15, [description = "Temperature [K]", input = true, bounds = (0.0, Inf)]
  end
  @variables begin
    G(t), [guess = 0.0, description = "Flow Rate [mol s-1]"]
    xH(t), [guess = 0.0, description = "Molar fraction [-]", output = true, bounds = (0.0, 1.0)]
    xD(t), [guess = 0.0, description = "Molar fraction [-]", output = true, bounds = (0.0, 1.0)]
    xT(t), [guess = 0.0, description = "Molar fraction [-]", output = true, bounds = (0.0, 1.0)]
    xHe(t), [guess = 0.0, description = "Molar fraction [-]", output = true, bounds = (0.0, 1.0)]
    xXe(t), [guess = 0.0, description = "Molar fraction [-]", output = true, bounds = (0.0, 1.0)]
    xI(t), [guess = 0.0, description = "Molar fraction [-]", output = true, bounds = (0.0, 1.0)]
  end

  @equations begin
    a.G ~ G
    a.xH ~ xH
    a.xD ~ xD
    a.xT ~ xT
    a.xHe ~ xHe
    a.xXe ~ xXe
    a.xI ~ xI
  end
end




# ============================================================================ #
#tag Cap 
# ============================================================================ #
@doc raw"""
    Cap(;name)

Zero-flow terminal cap for an `AFPOutlet` port. Sets the flow rate to zero and fixes the
composition to an inert D/T equimolar mixture. Use this block to seal unused `AFPOutlet`
connections.

# Connectors:
- `b`: `AFPOutlet` [outlet]
"""
@mtkmodel Cap begin
  @components begin
    b = AFPOutlet()
  end

  @parameters begin
  end

  @variables begin
  end

  @equations begin
    b.G   ~ 0.0
    b.xH  ~ 0.0
    b.xD  ~ 0.5
    b.xT  ~ 0.5
    b.xHe ~ 0.0
    b.xXe ~ 0.0
    b.xI  ~ 0.0

    b.dummy ~ 0.0
  end
end




# ============================================================================ #
#tag Cap_withO 
# ============================================================================ #
@doc raw"""
    Cap_withO(;name)

Zero-flow terminal cap for an `AFPOutlet_O` port. Sets the flow rate to zero and fixes the
composition to an inert D/T equimolar mixture with zero oxygen content. Use this block to
seal unused `AFPOutlet_O` connections.

# Connectors:
- `b`: `AFPOutlet_O` [outlet]
"""
@mtkmodel Cap_withO begin
  @components begin
    b = AFPOutlet_O()
  end

  @parameters begin
  end

  @variables begin
  end

  @equations begin
    b.G   ~ 0.0
    b.xH  ~ 0.0
    b.xD  ~ 0.5
    b.xT  ~ 0.5
    b.xHe ~ 0.0
    b.xXe ~ 0.0
    b.xI  ~ 0.0
    b.xO  ~ 0.0

    b.dummy ~ 0.0
  end
end




# ============================================================================ #
#tag CapOCN 
# ============================================================================ #
@doc raw"""
    CapOCN(;name)

Zero-flow terminal cap for an `AFPOutlet_OCN` port. Sets the flow rate to zero and fixes the
composition to an inert D/T equimolar mixture with zero O, C, and N content. Use this block
to seal unused `AFPOutlet_OCN` connections.

# Connectors:
- `b`: `AFPOutlet_OCN` [outlet]
"""
@mtkmodel CapOCN begin
  @components begin
    b = AFPOutlet_OCN()
  end

  @parameters begin
  end

  @variables begin
  end

  @equations begin
    b.G   ~ 0.0
    b.xH  ~ 0.0
    b.xD  ~ 0.5
    b.xT  ~ 0.5
    b.xHe ~ 0.0
    b.xXe ~ 0.0
    b.xI  ~ 0.0
    b.xO  ~ 0.0
    b.xC  ~ 0.0
    b.xN  ~ 0.0

    b.dummy ~ 0.0
  end
end




# ============================================================================ #
#tag Cap_HDTO 
# ============================================================================ #
@doc raw"""
    Cap_HDTO(;name)

Zero-flow terminal cap for an `AFPOutlet_HDTO` port. Sets the flow rate to zero and fixes
the composition to an inert D/T equimolar mixture with zero oxygen content. Use this block
to seal unused `AFPOutlet_HDTO` connections.

# Connectors:
- `b`: `AFPOutlet_HDTO` [outlet]
"""
@mtkmodel Cap_HDTO begin
  @components begin
    b = AFPOutlet_HDTO()
  end

  @parameters begin
  end

  @variables begin
  end

  @equations begin
    b.G   ~ 0.0
    b.xH  ~ 0.0
    b.xD  ~ 0.5
    b.xT  ~ 0.5
    b.xO  ~ 0.0

    b.dummy ~ 0.0
  end
end




# ============================================================================ #
#tag CapOutlet 
# ============================================================================ #
@doc raw"""
    CapOutlet(;name)

Terminal sink block for an `AFPInlet` port. Sets the inlet flow rate to zero and exposes
the local composition variables as observable outputs. Use this block to terminate unused
`AFPInlet` connections while retaining the ability to inspect the composition at that point.

# States:
- `xH(t)`: [-] H atomic fraction (output)
- `xD(t)`: [-] D atomic fraction (output)
- `xT(t)`: [-] T atomic fraction (output)
- `xHe(t)`: [-] He atomic fraction (output)
- `xXe(t)`: [-] Xe atomic fraction (output)
- `xI(t)`: [-] inert atomic fraction (output)

# Connectors:
- `a`: `AFPInlet` [inlet]
"""
@mtkmodel CapOutlet begin
  @components begin
    a = AFPInlet()
  end

  @parameters begin
  end
  
  @variables begin
    xH(t), [guess = 0.0, description = "Molar fraction [-]", output = true, bounds = (0.0, 1.0)]
    xD(t), [guess = 0.0, description = "Molar fraction [-]", output = true, bounds = (0.0, 1.0)]
    xT(t), [guess = 0.0, description = "Molar fraction [-]", output = true, bounds = (0.0, 1.0)]
    xHe(t), [guess = 0.0, description = "Molar fraction [-]", output = true, bounds = (0.0, 1.0)]
    xXe(t), [guess = 0.0, description = "Molar fraction [-]", output = true, bounds = (0.0, 1.0)]
    xI(t), [guess = 0.0, description = "Molar fraction [-]", output = true, bounds = (0.0, 1.0)]
  end

  @equations begin
    a.G   ~ 0.0
    a.xH  ~ xH
    a.xD  ~ xD
    a.xT  ~ xT
    a.xHe ~ xHe
    a.xXe ~ xXe
    a.xI  ~ xI
  end
end



# ============================================================================ #
#tag DTMassFlowController 
# ============================================================================ #
@doc raw"""
    DTMassFlowController(;name)

D/T Mass Flow Controller. Combines separate D and T storage inlet streams into a single
outlet stream whose total flow rate and tritium mole fraction are prescribed by external
control signals. The D flow is ``G_\mathrm{tot}(1 - x_T)`` and the T flow is
``G_\mathrm{tot} \, x_T``. Composition is propagated conservatively through the mix.

# Connectors:
- `aT`: `AFPInlet` [inlet] — tritium stream
- `aD`: `AFPInlet` [inlet] — deuterium stream
- `a_Gtot`: `RealInput` [input] — total flow rate setpoint [``\mathrm{s^{-1}}``]
- `a_xT`: `RealInput` [input] — tritium mole fraction setpoint [-]
- `b`: `AFPOutlet` [outlet] — mixed D/T stream
"""
@mtkmodel DTMassFlowController begin
  @components begin
    aT = AFPInlet()  # [inlet]
    aD = AFPInlet()  # [inlet]
    a_Gtot = RealInput() # [input] output from pressure control PID
    a_xT = RealInput() # [input] output from DT composition control PID
    b = AFPOutlet() # [outlet]
  end

  @equations begin
    aT.G + aD.G + b.G ~ 0.0
    aD.G ~ a_Gtot.u * (1 - a_xT.u)
    aT.G ~ a_Gtot.u * (a_xT.u)
    #
    aT.G * aT.xH  + aD.G * aD.xH  + b.G * b.xH  ~ 0.0
    aT.G * aT.xD  + aD.G * aD.xD  + b.G * b.xD  ~ 0.0
    aT.G * aT.xT  + aD.G * aD.xT  + b.G * b.xT  ~ 0.0
    aT.G * aT.xHe + aD.G * aD.xHe + b.G * b.xHe ~ 0.0
    aT.G * aT.xXe + aD.G * aD.xXe + b.G * b.xXe ~ 0.0
    aT.G * aT.xI  + aD.G * aD.xI  + b.G * b.xI  ~ 0.0

    b.dummy ~ 0.0
  end
end



# ============================================================================ #
#tag Storage 
# ============================================================================ #
@doc raw"""
    Storage(;name)

Ideal storage cylinder with a single outlet port carrying a fixed composition set by
parameters. The total inventory `M` decreases monotonically as the connected system
draws from the cylinder. The outlet flow rate is determined by the connected network.

# States:
- `M(t)`: [``\mathrm{mol}``] total moles remaining in the storage

# Connectors:
- `b`: `AFPOutlet` [outlet]

# Parameters:
- `xH`: [-] H species fraction
- `xD`: [-] D species fraction
- `xT`: [-] T species fraction
- `xHe`: [-] He species fraction
- `xXe`: [-] Xe species fraction
- `xI`: [-] inert species fraction
"""
@mtkmodel Storage begin
  @components begin
    b = AFPOutlet() # [outlet]
  end

  @parameters begin
    xH = 0.0, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xD = 0.0, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xT = 1.0, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xHe = 0.0, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xXe = 0.0, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xI = 0.0, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
  end

  @variables begin
    M(t), [guess = 100.0, description = "Total number of moles in the storage [mol]"]
  end

  @equations begin
    D(M) ~ -b.G
    xH ~ b.xH
    xD ~ b.xD
    xT ~ b.xT
    xHe ~ b.xHe
    xXe ~ b.xXe
    xI ~ b.xI

    b.dummy ~ 0.0
  end
end




# ============================================================================ #
#tag Storage2 
# ============================================================================ #
@doc raw"""
    Storage2(;name)

Ideal storage cylinder with a single outlet port and an externally controlled flow rate.
The composition is fixed by parameters; the flow rate is prescribed via the `a_G` control
input. The total inventory `M` evolves accordingly.

# States:
- `M(t)`: [``\mathrm{mol}``] total moles remaining in the storage

# Connectors:
- `b`: `AFPOutlet` [outlet]
- `a_G`: `RealInput` [input] — flow rate setpoint [``\mathrm{s^{-1}}``]

# Parameters:
- `xH`: [-] H species fraction
- `xD`: [-] D species fraction
- `xT`: [-] T species fraction
- `xHe`: [-] He species fraction
- `xXe`: [-] Xe species fraction
- `xI`: [-] inert species fraction
"""
@mtkmodel Storage2 begin
  @components begin
    b    = AFPOutlet() # [outlet]
    a_G  = RealInput() # [inlet]
  end

  @parameters begin
    xH  = 0.0, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xD  = 0.5, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xT  = 0.5, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xHe = 0.0, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xXe = 0.0, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xI  = 0.0, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
  end

  @variables begin
    M(t), [guess = 100.0, description = "Total number of moles in the storage [mol]"]
  end

  @equations begin
    D(M)   ~ -a_G.u
    b.G    ~ -a_G.u

    xH ~ b.xH
    xD ~ b.xD
    xT ~ b.xT
    xHe ~ b.xHe
    xXe ~ b.xXe
    xI ~ b.xI

    b.dummy ~ 0.0
  end
end






# ============================================================================ #
#tag ControlledDTStorage2 
# ============================================================================ #
@doc raw"""
    ControlledDTStorage2(;name)

Controlled D/T storage cylinder with an externally requested tritium flow rate. The outlet
tritium flow is limited to the available inventory divided by the characteristic time `τ`:
``\dot{N}_T^{\mathrm{out}} = \min(G_{\mathrm{req}},\, N_T / \tau)``. The `error_GT` output
reports the unsatisfied demand; a positive value indicates the storage cannot meet the
requested flow. The composition is fixed by parameters.

# States:
- `M(t)`: [``\mathrm{mol}``] total moles remaining in the storage

# Connectors:
- `b`: `AFPOutlet` [outlet]
- `a_G`: `RealInput` [input] — requested tritium flow rate [``\mathrm{s^{-1}}``]
- `error_GT`: `RealOutput` [output] — unsatisfied flow demand [``\mathrm{s^{-1}}``]

# Parameters:
- `τ`: [s] cylinder characteristic delivery time
- `xH`: [-] H species fraction
- `xD`: [-] D species fraction
- `xT`: [-] T species fraction
- `xHe`: [-] He species fraction
- `xXe`: [-] Xe species fraction
- `xI`: [-] inert species fraction
"""
@mtkmodel ControlledDTStorage2 begin
  @components begin
    b    = AFPOutlet() # [outlet]
    a_G  = RealInput() # [inlet]
    error_GT  = RealOutput() # [outlet]
  end

  @parameters begin
    τ   = 1.0, [description = "cylinder characteristic time [s]", input = true, bounds = (1.0e-4, Inf)]
    xH  = 0.0, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xD  = 0.5, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xT  = 0.5, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xHe = 0.0, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xXe = 0.0, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xI  = 0.0, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
  end

  @variables begin
    M(t), [guess = 100.0, description = "Total number of moles in the storage [mol]"]
  end

  @equations begin
    D(M)     ~ b.G
    b.G*b.xT ~ -min(a_G.u, M*xT/τ)
    error_GT.u    ~ a_G.u + b.G*b.xT # if > 0 then something else has to provide this tritium flow

    xH ~ b.xH
    xD ~ b.xD
    xT ~ b.xT
    xHe ~ b.xHe
    xXe ~ b.xXe
    xI ~ b.xI

    b.dummy ~ 0.0
  end
end




# ============================================================================ #
#tag ControlledDTStorage2_tank 
# ============================================================================ #
@doc raw"""
    ControlledDTStorage2_tank(;name)

Controlled D/T storage tank with two upstream inlet ports (`a1`, `a2`) that replenish the
inventory, one outlet port (`b`), and an externally requested tritium flow rate. The tank
integrates species-wise conservation equations for all six species. The outlet tritium flow is
limited to `min(G_req, NT/τ)`. The `error_GT` output reports the unsatisfied tritium demand.
The outlet composition reflects the instantaneous mole fractions in the tank.

Tritium decay via ``{}^3\mathrm{H} \to {}^3\mathrm{He}`` is included in the `D(NT)` balance.

# States:
- `NH(t)`: [``\mathrm{mol}``] monoatomic H inventory
- `ND(t)`: [``\mathrm{mol}``] monoatomic D inventory
- `NT(t)`: [``\mathrm{mol}``] monoatomic T inventory
- `NHe(t)`: [``\mathrm{mol}``] monoatomic He inventory
- `NXe(t)`: [``\mathrm{mol}``] monoatomic Xe inventory
- `NI(t)`: [``\mathrm{mol}``] monoatomic inert inventory

# Connectors:
- `a1`, `a2`: `AFPOutlet` [inlet] — replenishment streams (note: typed as `AFPOutlet` to match upstream convention)
- `b`: `AFPOutlet` [outlet] — delivery stream
- `a_G`: `RealInput` [input] — requested tritium flow rate [``\mathrm{s^{-1}}``]
- `error_GT`: `RealOutput` [output] — unsatisfied flow demand [``\mathrm{s^{-1}}``]

# Parameters:
- `τ`: [s] cylinder characteristic delivery time
"""
@mtkmodel ControlledDTStorage2_tank begin
  @components begin
    a1    = AFPOutlet() # [inlet]
    a2    = AFPOutlet() # [inlet]
    b    = AFPOutlet() # [outlet]
    a_G  = RealInput()       # [inlet]
    error_GT  = RealOutput() # [outlet]
  end

  @parameters begin
    τ   = 1.0, [description = "cylinder characteristic time [s]", input = true, bounds = (1.0e-4, Inf)]
  end

  @variables begin
    NH(t),  [guess = 1.0e-20, description = "moles (monoatomic) [mol]", bounds = (0.0, Inf)]
    ND(t),  [guess = 1.0e-20, description = "moles (monoatomic) [mol]", bounds = (0.0, Inf)]
    NT(t),  [guess = 1.0e-20, description = "moles (monoatomic) [mol]", bounds = (0.0, Inf)]
    NHe(t), [guess = 1.0e-20, description = "moles (monoatomic) [mol]", bounds = (0.0, Inf)]
    NXe(t), [guess = 1.0e-20, description = "moles (monoatomic) [mol]", bounds = (0.0, Inf)]
    NI(t),  [guess = 1.0e-20, description = "moles (monoatomic) [mol]", bounds = (0.0, Inf)]
  end

  @equations begin
    # equations for block variables
    D(NH)  ~ b.xH  * b.G + a1.xH  * a1.G + a2.xH  * a2.G
    D(ND)  ~ b.xD  * b.G + a1.xD  * a1.G + a2.xD  * a2.G
    D(NT)  ~ b.xT  * b.G + a1.xT  * a1.G + a2.xT  * a2.G - λT*NT
    D(NHe) ~ b.xHe * b.G + a1.xHe * a1.G + a2.xHe * a2.G
    D(NXe) ~ b.xHe * b.G + a1.xHe * a1.G + a2.xHe * a2.G
    D(NI)  ~ b.xI  * b.G + a1.xI  * a1.G + a2.xI  * a2.G

    # equations for state outlet b
    b.xH  ~ NH  / (NH + ND + NT + NHe + NXe + NI)
    b.xD  ~ ND  / (NH + ND + NT + NHe + NXe + NI)
    b.xT  ~ NT  / (NH + ND + NT + NHe + NXe + NI)
    b.xHe ~ NHe / (NH + ND + NT + NHe + NXe + NI)
    b.xXe ~ NXe / (NH + ND + NT + NHe + NXe + NI)
    b.xI  ~ NI  / (NH + ND + NT + NHe + NXe + NI)

    # Eq. for outlet flow
    b.G*b.xT ~ -min(a_G.u, NT/τ)
    error_GT.u    ~ a_G.u + b.G*b.xT # if > 0 then something else has to provide this tritium flow

    b.dummy ~ 0.0
  end
end



# ============================================================================ #
#tag Storage2_withO 
# ============================================================================ #
@doc raw"""
    Storage2_withO(;name)

Ideal storage cylinder with oxygen tracking, a single outlet port (`AFPOutlet_O`), and an
externally controlled flow rate. The composition is fixed by parameters including `xO`. The
flow rate is prescribed via the `a_G` control input and the total inventory `M` evolves
accordingly.

# States:
- `M(t)`: [``\mathrm{mol}``] total moles remaining in the storage

# Connectors:
- `b`: `AFPOutlet_O` [outlet]
- `a_G`: `RealInput` [input] — flow rate setpoint [``\mathrm{s^{-1}}``]

# Parameters:
- `xH`: [-] H species fraction
- `xD`: [-] D species fraction
- `xT`: [-] T species fraction
- `xHe`: [-] He species fraction
- `xXe`: [-] Xe species fraction
- `xI`: [-] inert species fraction
- `xO`: [-] O species fraction
"""
@mtkmodel Storage2_withO begin
  @components begin
    b = AFPOutlet_O() # [outlet]
    a_G = RealInput() # [inlet]
  end

  @parameters begin
    xH  = 0.0, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xD  = 0.5, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xT  = 0.5, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xHe = 0.0, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xXe = 0.0, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xI  = 0.0, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xO  = 0.0, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
  end

  @variables begin
    M(t), [guess = 100.0, description = "Total number of moles in the storage [mol]"]
  end

  @equations begin
    D(M) ~ -a_G.u
    b.G ~ -a_G.u
    xH ~ b.xH
    xD ~ b.xD
    xT ~ b.xT
    xHe ~ b.xHe
    xXe ~ b.xXe
    xI ~ b.xI
    xO ~ b.xO

    b.dummy ~ 0.0
  end
end






# ============================================================================ #
#tag StorageHDTO 
# ============================================================================ #
@doc raw"""
    StorageHDTO(;name)

Ideal storage cylinder for H, D, T, O mixtures, with a single `AFPOutlet_HDTO` outlet and
an externally controlled flow rate. The composition is fixed by parameters; the flow rate is
prescribed via the `a_G` control input. The total inventory `M` evolves accordingly.

# States:
- `M(t)`: [``\mathrm{mol}``] total moles remaining in the storage

# Connectors:
- `b`: `AFPOutlet_HDTO` [outlet]
- `a_G`: `RealInput` [input] — flow rate setpoint [``\mathrm{s^{-1}}``]

# Parameters:
- `xH`: [-] H species fraction
- `xD`: [-] D species fraction
- `xT`: [-] T species fraction
- `xO`: [-] O species fraction
"""
@mtkmodel StorageHDTO begin
  @components begin
    b = AFPOutlet_HDTO() # [outlet]
    a_G = RealInput() # [inlet]
  end

  @parameters begin
    xH  = 0.0, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xD  = 0.5, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xT  = 0.5, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xO  = 0.0, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
  end

  @variables begin
    M(t), [guess = 100.0, description = "Total number of moles in the storage [mol]"]
  end

  @equations begin
    D(M) ~ -a_G.u
    b.G ~ -a_G.u
    xH ~ b.xH
    xD ~ b.xD
    xT ~ b.xT
    xO ~ b.xO

    b.dummy ~ 0.0
  end
end



# ============================================================================ #
#tag Storage3 
# ============================================================================ #
@doc raw"""
    Storage3(;name)

Ideal storage cylinder for H, D, T, He, Xe, I, O, C mixtures, with a single `AFPOutlet_OC`
outlet and an externally controlled flow rate. Intended as the CO source connected to a
`MembraneReactor`. The composition is fixed by parameters; the flow rate is prescribed via
the `a_G` control input. The total inventory `M` evolves accordingly.

# States:
- `M(t)`: [``\mathrm{mol}``] total moles remaining in the storage

# Connectors:
- `b`: `AFPOutlet_OC` [outlet]
- `a_G`: `RealInput` [input] — flow rate setpoint [``\mathrm{s^{-1}}``]

# Parameters:
- `xH`: [-] H species fraction
- `xD`: [-] D species fraction
- `xT`: [-] T species fraction
- `xHe`: [-] He species fraction
- `xXe`: [-] Xe species fraction
- `xI`: [-] inert species fraction
- `xO`: [-] O species fraction
- `xC`: [-] C species fraction
"""
@mtkmodel Storage3 begin
  @components begin
    b = AFPOutlet_OC() # [outlet]
    a_G = RealInput() # [inlet]
  end

  @parameters begin
    xH  = 0.0, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xD  = 0.5, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xT  = 0.5, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xHe = 0.0, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xXe = 0.0, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xI  = 0.0, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xO  = 0.0, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xC  = 0.0, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
  end

  @variables begin
    M(t), [guess = 100.0, description = "Total number of moles in the storage [mol]"]
  end

  @equations begin
    D(M) ~ -a_G.u
    b.G ~ -a_G.u
    xH ~ b.xH
    xD ~ b.xD
    xT ~ b.xT
    xHe ~ b.xHe
    xXe ~ b.xXe
    xI ~ b.xI
    xO ~ b.xO
    xC ~ b.xC

    b.dummy ~ 0.0
  end
end




# ============================================================================ #
#tag StorageOCN 
# ============================================================================ #
@doc raw"""
    StorageOCN(;name)

Ideal storage cylinder for the full H, D, T, He, Xe, I, O, C, N species set, with a single
`AFPOutlet_OCN` outlet and an externally controlled flow rate. The composition is fixed by
parameters; the flow rate is prescribed via the `a_G` control input. The total inventory `M`
evolves accordingly.

# States:
- `M(t)`: [``\mathrm{mol}``] total moles remaining in the storage

# Connectors:
- `b`: `AFPOutlet_OCN` [outlet]
- `a_G`: `RealInput` [input] — flow rate setpoint [``\mathrm{s^{-1}}``]

# Parameters:
- `xH`: [-] H species fraction
- `xD`: [-] D species fraction
- `xT`: [-] T species fraction
- `xHe`: [-] He species fraction
- `xXe`: [-] Xe species fraction
- `xI`: [-] inert species fraction
- `xO`: [-] O species fraction
- `xC`: [-] C species fraction
- `xN`: [-] N species fraction
"""
@mtkmodel StorageOCN begin
  @components begin
    b = AFPOutlet_OCN() # [outlet]
    a_G = RealInput() # [inlet]
  end

  @parameters begin
    xH  = 0.0, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xD  = 0.5, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xT  = 0.5, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xHe = 0.0, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xXe = 0.0, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xI  = 0.0, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xO  = 0.0, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xC  = 0.0, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xN  = 0.0, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
  end

  @variables begin
    M(t), [guess = 100.0, description = "Total number of moles in the storage [mol]"]
  end

  @equations begin
    D(M) ~ -a_G.u
    b.G ~ -a_G.u
    xH ~ b.xH
    xD ~ b.xD
    xT ~ b.xT
    xHe ~ b.xHe
    xXe ~ b.xXe
    xI ~ b.xI
    xO ~ b.xO
    xC ~ b.xC
    xN ~ b.xN

    b.dummy ~ 0.0
  end
end




@doc raw"""
    Storage2temp(;name)

Temporary storage cylinder with a fixed constant outlet flow rate of ``10^{19}\,N_A^{-1}``
[``\mathrm{mol\,s^{-1}}``]. Composition is fixed by parameters. This block is intended for
use during model development and testing where a constant-rate source is required without
an external control input.

# States:
- `M(t)`: [``\mathrm{mol}``] total moles remaining in the storage

# Connectors:
- `b`: `AFPOutlet` [outlet]
- `a_G`: `RealInput` [input] — flow rate input (not used in current equations; reserved for compatibility)

# Parameters:
- `xH`: [-] H species fraction
- `xD`: [-] D species fraction
- `xT`: [-] T species fraction
- `xHe`: [-] He species fraction
- `xXe`: [-] Xe species fraction
- `xI`: [-] inert species fraction
"""
@mtkmodel Storage2temp begin
  @components begin
    b = AFPOutlet() # [outlet]
    a_G = RealInput() # [inlet]
  end

  @parameters begin
    xH  = 0.0, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xD  = 0.5, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xT  = 0.5, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xHe = 0.0, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xXe = 0.0, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xI  = 0.0, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
  end

  @variables begin
    M(t), [guess = 100.0, description = "Total number of moles in the storage [mol]"]
  end

  @equations begin
    D(M) ~ -1.0e19/Na
    b.G ~ -1.0e19/Na
    xH ~ b.xH
    xD ~ b.xD
    xT ~ b.xT
    xHe ~ b.xHe
    xXe ~ b.xXe
    xI ~ b.xI

    b.dummy ~ 0.0
  end
end


# ============================================================================ #
#tag IdealSource
# ============================================================================ #
@doc raw"""
    IdealSource(;name)

Ideal flow source with a fixed composition and an externally prescribed flow rate. The
outlet flow rate equals the signal on `G`; there is no internal inventory tracked. Use
this block to impose a prescribed time-varying flow boundary condition with a known
composition.

# Connectors:
- `b`: `AFPOutlet` [outlet]
- `G`: `RealInput` [input] — flow rate [``\mathrm{s^{-1}}``]

# Parameters:
- `xH`: [-] H species fraction
- `xD`: [-] D species fraction
- `xT`: [-] T species fraction
- `xHe`: [-] He species fraction
- `xXe`: [-] Xe species fraction
- `xI`: [-] inert species fraction
"""
@mtkmodel IdealSource begin
  @components begin
    b = AFPOutlet() # [outlet]
    G = RealInput() # [inlet]
  end

  @parameters begin
    xH = 0.0, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xD = 0.5, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xT = 0.5, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xHe = 0.0, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xXe = 0.0, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
    xI = 0.0, [description = "species fraction [-]", input = true, bounds = (0.0, 1.0)]
  end

  @variables begin
  end

  @equations begin
    b.G ~ -G.u
    xH ~ b.xH
    xD ~ b.xD
    xT ~ b.xT
    xHe ~ b.xHe
    xXe ~ b.xXe
    xI ~ b.xI

    b.dummy ~ 0.0
  end
end

# ============================================================================ #
#tag SimpleControlFlowMeter 
# ============================================================================ #
@doc raw"""
    SimpleControlFlowMeter(;name)

Flow meter combined with a simple proportional flow controller. The block passes composition
through without modification and exposes the measured flow rate and the fixed setpoint as
`RealOutput` signals for use with an external PID controller. The actual outlet flow rate
is set to the measured value of the passing flow (pass-through meter).

# Connectors:
- `a`: `AFPInlet` [inlet]
- `b`: `AFPOutlet` [outlet]
- `a_ctrl`: `RealInput` [input] — control signal (unused in current equations; reserved)
- `G_ctrl`: `RealOutput` [output] — measured flow rate [``\mathrm{s^{-1}}``]
- `ctrl_setpoint`: `RealOutput` [output] — flow rate setpoint [``\mathrm{s^{-1}}``]

# Parameters:
- `setpoint`: [``\mathrm{s^{-1}}``] desired flow rate setpoint
"""
@mtkmodel SimpleControlFlowMeter begin
  @components begin
    a = AFPInlet() # [inlet]
    b = AFPOutlet() # [outlet]
    a_ctrl = RealInput()
    G_ctrl = RealOutput() # [output]
    ctrl_setpoint = RealOutput() # [output]
  end

  @parameters begin
    setpoint = 0.0, [description = "Desired flow rate [mol s-1]", input = true, bounds = (0.0, Inf)]
  end
  @equations begin
    a.xH ~ b.xH
    a.xD ~ b.xD
    a.xT ~ b.xT
    a.xHe ~ b.xHe
    a.xXe ~ b.xXe
    a.xI ~ b.xI
    0.0 ~ a.G + b.G
    #
    a.G ~ G_ctrl.u
    setpoint ~ ctrl_setpoint.u

    b.dummy ~ 0.0
  end
end



# ============================================================================ #
#tag FlowMeter 
# ============================================================================ #
@doc raw"""
    FlowMeter(;name)

In-line flow and composition meter. Passes all species through without modification and
exposes the full flow state (flow rate and all mole fractions) on a `ControlPortOut`
connector for use with downstream control logic.

# States:
- `G(t)`: [``\mathrm{s^{-1}}``] measured flow rate
- `xH(t)`: [-] measured H mole fraction
- `xD(t)`: [-] measured D mole fraction
- `xT(t)`: [-] measured T mole fraction
- `xHe(t)`: [-] measured He mole fraction
- `xXe(t)`: [-] measured Xe mole fraction
- `xI(t)`: [-] measured inert mole fraction

# Connectors:
- `a`: `AFPInlet` [inlet]
- `b`: `AFPOutlet` [outlet]
- `ctrl`: `ControlPortOut` [outlet] — full flow state for control logic
"""
@mtkmodel FlowMeter begin
  @components begin
    a = AFPInlet() # [inlet]
    b = AFPOutlet() # [outlet]
    ctrl = ControlPortOut() # [outlet]
  end

  @variables begin
    G(t), [guess = 0.0, description = "Flow Rate [mol s-1]", bounds = (0.0, Inf)]
    xH(t), [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xD(t), [guess = 0.5, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xT(t), [guess = 0.5, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xHe(t), [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xXe(t), [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xI(t), [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
  end

  @equations begin
    xH ~ a.xH
    xD ~ a.xD
    xT ~ a.xT
    xHe ~ a.xHe
    xXe ~ a.xXe
    xI ~ a.xI
    G ~ a.G
    #
    0.0 ~ xH * G + b.xH * b.G
    0.0 ~ xD * G + b.xD * b.G
    0.0 ~ xT * G + b.xT * b.G
    0.0 ~ xHe * G + b.xHe * b.G
    0.0 ~ xXe * G + b.xXe * b.G
    0.0 ~ xI * G + b.xI * b.G
    0.0 ~ G + b.G
    # Control
    ctrl.G ~ G
    ctrl.xH ~ xH
    ctrl.xD ~ xD
    ctrl.xT ~ xT
    ctrl.xHe ~ xHe
    ctrl.xXe ~ xXe
    ctrl.xI ~ xI

    b.dummy ~ 0.0
  end
end




# ============================================================================ #
#tag FlowMeterOCN 
# ============================================================================ #
@doc raw"""
    FlowMeterOCN(;name)

In-line flow and composition meter for the full OCN species set (H, D, T, He, Xe, I, O, C,
N). Passes all species through without modification. The flow rate is a fixed parameter
rather than a state, which constrains the network flow at this point.

# States:
- `xH(t)`: [-] measured H mole fraction
- `xD(t)`: [-] measured D mole fraction
- `xT(t)`: [-] measured T mole fraction
- `xHe(t)`: [-] measured He mole fraction
- `xXe(t)`: [-] measured Xe mole fraction
- `xI(t)`: [-] measured inert mole fraction
- `xO(t)`: [-] measured O mole fraction
- `xC(t)`: [-] measured C mole fraction
- `xN(t)`: [-] measured N mole fraction

# Connectors:
- `a`: `AFPInlet_OCN` [inlet]
- `b`: `AFPOutlet_OCN` [outlet]

# Parameters:
- `G`: [``\mathrm{s^{-1}}``] flow rate constraint
"""
@mtkmodel FlowMeterOCN begin
  @components begin
    a = AFPInlet_OCN() # [inlet]
    b = AFPOutlet_OCN() # [outlet]
  end

  @parameters begin
    G = 0.0, [description = "flow rate [s-1]", input = true, bounds = (0.0, Inf)]
  end
  
  @variables begin
    xH(t), [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xD(t), [guess = 0.5, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xT(t), [guess = 0.5, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xHe(t), [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xXe(t), [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xI(t), [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xO(t), [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xC(t), [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xN(t), [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
  end

  @equations begin
    xH ~ a.xH
    xD ~ a.xD
    xT ~ a.xT
    xHe ~ a.xHe
    xXe ~ a.xXe
    xI ~ a.xI
    xO ~ a.xO
    xC ~ a.xC
    xN ~ a.xN
    G ~ a.G
    #
    0.0 ~ xH * G + b.xH * b.G
    0.0 ~ xD * G + b.xD * b.G
    0.0 ~ xT * G + b.xT * b.G
    0.0 ~ xHe * G + b.xHe * b.G
    0.0 ~ xXe * G + b.xXe * b.G
    0.0 ~ xI * G + b.xI * b.G
    0.0 ~ xO * G + b.xO * b.G
    0.0 ~ xC * G + b.xC * b.G
    0.0 ~ xN * G + b.xN * b.G
    0.0 ~ a.G + b.G

    b.dummy ~ 0.0
  end
end




# ============================================================================ #
#tag SpeciesFlowMeter 
# ============================================================================ #
@doc raw"""
    SpeciesFlowMeter(;name)

In-line selective species flow meter. Passes all species through without modification and
outputs on `ctrl` the partial flow rate of a user-selected subset of species. The selection
is controlled by the binary flags `jH`, `jD`, `jT`, `jHe`, `jXe`, `jI`; setting a flag
to 1 includes the corresponding species in the output signal.

The control output is:
```math
G_{\mathrm{ctrl}} = G \sum_i j_i \, x_i, \quad i \in \{H, D, T, He, Xe, I\}
```

# States:
- `G(t)`: [``\mathrm{s^{-1}}``] total flow rate
- `xH(t)`, `xD(t)`, `xT(t)`, `xHe(t)`, `xXe(t)`, `xI(t)`: [-] mole fractions

# Connectors:
- `a`: `AFPInlet` [inlet]
- `b`: `AFPOutlet` [outlet]
- `ctrl`: `RealOutput` [outlet] — selected partial flow rate [``\mathrm{s^{-1}}``]

# Parameters:
- `jH`, `jD`, `jT`, `jHe`, `jXe`, `jI`: [-] species selection flags (0 = exclude, 1 = include)
"""
@mtkmodel SpeciesFlowMeter begin
  @components begin
    a = AFPInlet() # [inlet]
    b = AFPOutlet() # [outlet]
    ctrl = RealOutput() # [outlet]
  end

  @variables begin
    G(t), [guess = 0.0, description = "Flow Rate [mol s-1]", bounds = (0.0, Inf)]
    xH(t), [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xD(t), [guess = 0.5, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xT(t), [guess = 0.5, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xHe(t), [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xXe(t), [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xI(t), [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
  end

  @parameters begin
    jH  = 0, [description = "flag for the measurement of H, 0 = flase, 1=true [-]",  input = true, bounds = (0, 1)]
    jD  = 0, [description = "flag for the measurement of D, 0 = flase, 1=true [-]",  input = true, bounds = (0, 1)]
    jT  = 0, [description = "flag for the measurement of T, 0 = flase, 1=true [-]",  input = true, bounds = (0, 1)]
    jHe = 0, [description = "flag for the measurement of He, 0 = flase, 1=true [-]", input = true, bounds = (0, 1)]
    jXe = 0, [description = "flag for the measurement of He, 0 = flase, 1=true [-]", input = true, bounds = (0, 1)]
    jI  = 0, [description = "flag for the measurement of I, 0 = flase, 1=true [-]",  input = true, bounds = (0, 1)]
  end

  @equations begin
    xH ~ a.xH
    xD ~ a.xD
    xT ~ a.xT
    xHe ~ a.xHe
    xXe ~ a.xXe
    xI ~ a.xI
    G ~ a.G
    #
    0.0 ~ xH * G + b.xH * b.G
    0.0 ~ xD * G + b.xD * b.G
    0.0 ~ xT * G + b.xT * b.G
    0.0 ~ xHe * G + b.xHe * b.G
    0.0 ~ xXe * G + b.xXe * b.G
    0.0 ~ xI * G + b.xI * b.G
    0.0 ~ G + b.G
    # Control
    ctrl.u ~ G * (jH*xH + jD*xD + jT*xT + jHe*xHe + jXe*xXe + jI*xI)

    b.dummy ~ 0.0
  end
end


# ============================================================================ #
#tag modPelletInjector PI
# ============================================================================ #
@doc raw"""
    modPelletInjector(;name)

Pellet Injector (PI) model. Separates a mixed D/T/He gas stream into a purified Q stream
delivered to the plasma (`b2`) and a helium-rich purge stream (`b1`). The helium separation
efficiency `eta` governs what fraction of He is diverted to the purge outlet; all Q and
impurity species exit through `b2`. Characteristic times `τ_Q` and `τ_He` set the Q and He
residence times respectively.

Tritium decay via ``{}^3\mathrm{H} \to {}^3\mathrm{He}`` is included in the `D(NT)` balance.

# States:
- `NH(t)`: [``\mathrm{s^{-1}}``] total H inventory
- `ND(t)`: [``\mathrm{s^{-1}}``] total D inventory
- `NT(t)`: [``\mathrm{s^{-1}}``] total T inventory
- `NHe(t)`: [``\mathrm{s^{-1}}``] total He inventory
- `NXe(t)`: [``\mathrm{s^{-1}}``] total Xe inventory
- `NI(t)`: [``\mathrm{s^{-1}}``] total inert inventory

# Connectors:
- `a`: `AFPInlet` [inlet] — mixed fuel + He stream
- `b1`: `AFPOutlet` [outlet] — He-rich purge stream
- `b2`: `AFPOutlet` [outlet] — Q-rich stream to plasma

# Parameters:
- `τ_Q`: [s] Q species characteristic time
- `τ_He`: [s] He characteristic time
- `eta`: [-] He separation efficiency
"""
@mtkmodel modPelletInjector begin
  @components begin
    a  = AFPInlet() # [inlet]
    b1 = AFPOutlet() # [outlet]
    b2 = AFPOutlet() # [outlet]
  end

  @variables begin
    NH(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    ND(t),  [guess = 1.0e20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NT(t),  [guess = 1.0e20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NHe(t), [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NXe(t), [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NI(t),  [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
  end

  @parameters begin
    τ_Q  = 1.0, [description = "Q characteristic time in the system [s]", input = true, bounds = (0.0, Inf)]
    τ_He = 1.0, [description = "He characteristic time in the system [s]", input = true, bounds = (0.0, Inf)]
    eta  = 0.95, [description = "Pellet Injector Helium separation efficiency [-]", input = true, bounds = (0.0, 1.0)]
  end

  @equations begin
    D(NH)  ~ a.xH * a.G + b1.xH * b1.G + b2.xH * b2.G
    D(ND)  ~ a.xD * a.G + b1.xD * b1.G + b2.xD * b2.G
    D(NT)  ~ a.xT * a.G + b1.xT * b1.G + b2.xT * b2.G - λT*NT
    D(NHe) ~ a.xHe * a.G + b1.xHe * b1.G + b2.xHe * b2.G
    D(NXe) ~ a.xXe * a.G + b1.xXe * b1.G + b2.xXe * b2.G
    D(NI)  ~ a.xI * a.G + b1.xI * b1.G + b2.xI * b2.G
    #D(NH) + D(ND) + D(NT) + D(NHe) ~ a.G*(a.xH + a.xD + a.xT + a.xHe) + b1.G*(b1.xH + b1.xD + b1.xT + b1.xHe) + b2.G*(b2.xH + b2.xD + b2.xT + b2.xHe)
    #
    b2.xH * b2.G  ~ -NH/τ_Q
    b2.xD * b2.G  ~ -ND/τ_Q
    b2.xT * b2.G  ~ -NT/τ_Q
    b2.xHe * b2.G ~ -(1 - eta)*NHe/τ_Q
    b2.xXe * b2.G ~ -NXe/τ_Q
    b2.xI * b2.G  ~ -NI/τ_Q
    b2.G ~ -(NH + ND + NT + NXe + NI)/τ_Q  +  NHe*(1-eta)/τ_He
    #
    b1.xH  ~ 0.0
    b1.xD  ~ 0.0
    b1.xT  ~ 0.0
    b1.xHe * b1.G ~ -eta*NHe/τ_He
    b1.xXe ~ 0.0
    b1.xI  ~ 0.0
    b1.G ~ -eta*NHe/τ_He

    b1.dummy ~ 0.0
    b2.dummy ~ 0.0
  end
end



# ============================================================================ #
#tag Plasma 
# ============================================================================ #
@doc raw"""
    Plasma(;name)

Torus model with with an inlet port `a`, one outlet port `b`, and a fluid flow with 4 species.

# States:
- `xH(t)`:  [-] outlet H molar fraction
- `xD(t)`:  [-] outlet D molar fraction
- `xT(t)`:  [-] otulet T molar fraction
- `xHe(t)`: [-] outlet He molar fraction
- `xI(t)`:  [-] outlet He molar fraction
- `G_H(t)`: [``\mathrm{s^{-1}}``] protium production rate
- `G_alpha(t)`: [``\mathrm{s^{-1}}``] ``\alpha`` particle production rate
- `DNH(t)`: [``\mathrm{s^{-1}}``] time derivative of protium population

# Connectors:
- `a`: inlet fluid port
- `b` : outlet fluid port

# Parameters:
- `a`: inlet fluid port
- `b` : outlet fluid port

- `τp`:  [-] particle confinement time
- `τp_α`:  [-] alpha particle confinement time
- `Vₚₗₐₛₘₐ`: [-] plasma volume
- `n`: [-] initial fuel density
- `rR_DD_1`: [``\mathrm{m^6 s^{-1}}``] reaction rate
- `rR_DD_2`: [``\mathrm{m^6 s^{-1}}``] reaction rate
- `rR_DT`: [``\mathrm{m^6 s^{-1}}``] reaction rate
- `rR_DHe3`: [``\mathrm{m^6 s^{-1}}``] reaction rate
- `rR_TT`: [``\mathrm{m^6 s^{-1}}``] reaction rate
- `rR_THe3_1`: [``\mathrm{m^6 s^{-1}}``] reaction rate
- `rR_THe3_2`: [``\mathrm{m^6 s^{-1}}``] reaction rate
- `rR_THe3_3`: [``\mathrm{m^6 s^{-1}}``] reaction rate
- `Wcorrection`: [-] density correction factor


# Equations:

```math
\begin{aligned}
	V_{pl}\frac{\partial n_H}{\partial t} =\ & x_H^{in} \Gamma^{in} +
											x_H^{out} \Gamma^{out} +
											w_c n_D n_{He} R_{DHe^3} +
											w_c n_D n_D R_{DD-1} + \\
											& + w_c n_T n_{He} (R_{THe^3-1} + R_{THe^3-3})\\

	V_{pl}\frac{\partial n_D}{\partial t} =\ & x_D^{in} \Gamma^{in} +
											x_D^{out} \Gamma^{out} +
											w_c n_T n_{He} R_{THe^3-2} -
											w_c n_D n_{He} R_{DHe^3} - \\
											& - w_c n_D n_T R_{DT} -
										    w_c 2 n_D n_D (R_{DD-1} + R_{DD-2})\\

	V_{pl}\frac{\partial n_T}{\partial t} =\ & x_T^{in} \Gamma^{in} +
											x_T^{out} \Gamma^{out} +
											w_c n_D n_D R_{DD-1} -
											w_c n_D n_T R_{DT} - \\
											& - w_c 2 n_T n_T R_{TT} -
											w_c n_T n_{He} (R_{THe^3-1} + R_{THe^3-2} + R_{THe^3-3})\\

	V_{pl}\frac{\partial n_{He}}{\partial t} =\ & x_{He}^{in} \Gamma^{in} +
											   x_{He}^{out} \Gamma^{out} +
											   w_c n_D n_D R_{DD-2} +
											   w_c n_D n_T R_{DT} +
											   w_c n_T n_T R_{TT}\\

  x_i \Gamma^{out} =\ & - n_i V_{pl} / \tau_p,\ \ \ \ i = H,D,T,He\\

  \sum_{i} x_i =\ & 1,\ \ \ \ i = H,D,T,H\\
\end{aligned}
```
"""
@mtkmodel Plasma begin
  @components begin
    a = AFPInlet() # [inlet]
    b = AFPOutlet() # [outlet]
  end


  @parameters begin
    τp = 1.0, [description = "particle confinement time [s]", input = true, bounds = (0.0, 100000.0)]
    τp_α = 1.0, [description = "alpha particle confinement time [s]", input = true, bounds = (0.0, 100000.0)]
    Vₚₗₐₛₘₐ = 2000.0, [description = "plasma volume [m3]", input = true, bounds = (0.0, 5000.0)]
    n = 1.0e20, [description = "D, T fuel density, averaged over a [DT particles/m3]", input = true, bounds = (0.0, 1.0e40)]
    rR_DD_1 = 0.0, [description = "reaction rate, must be multiplied by n^2 [m6 s-1]", input = true, bounds = (0.0, Inf)]
    rR_DD_2 = 0.0, [description = "reaction rate, must be multiplied by n^2 [m6 s-1]", input = true, bounds = (0.0, Inf)]
    rR_DT = 0.0, [description = "reaction rate, must be multiplied by n^2 [m6 s-1]", input = true, bounds = (0.0, Inf)]
    rR_DHe3 = 0.0, [description = "reaction rate, must be multiplied by n^2 [m6 s-1]", input = true, bounds = (0.0, Inf)]
    rR_TT = 0.0, [description = "reaction rate, must be multiplied by n^2 [m6 s-1]", input = true, bounds = (0.0, Inf)]
    rR_THe3_1 = 0.0, [description = "reaction rate, must be multiplied by n^2 [m6 s-1]", input = true, bounds = (0.0, Inf)]
    rR_THe3_2 = 0.0, [description = "reaction rate, must be multiplied by n^2 [m6 s-1]", input = true, bounds = (0.0, Inf)]
    rR_THe3_3 = 0.0, [description = "reaction rate, must be multiplied by n^2 [m6 s-1]", input = true, bounds = (0.0, Inf)]
    wc = 0.0, [description = "correction factor for the averaged density [-]", input = true, bounds = (1.0, 2.0)]
  end


  @variables begin
    #
    nH(t), [guess = 1.0e20, description = "species density [m-3]", bounds = (0.0, Inf)]
    nD(t), [guess = 1.0e20, description = "species density [m-3]", bounds = (0.0, Inf)]
    nT(t), [guess = 1.0e20, description = "species density [m-3]", bounds = (0.0, Inf)]
    nHe(t), [guess = 1.0e20, description = "species density [m-3]", bounds = (0.0, Inf)]
    nXe(t), [guess = 1.0e20, description = "species density [m-3]", bounds = (0.0, Inf)]
    nI(t), [guess = 1.0e20, description = "species density [m-3]", bounds = (0.0, Inf)]
    #nHe45(t),   [guess = 1.0e20, description = "species density [m-3]", bounds = (0.0, Inf)]
    G_H(t), [guess = 0.0, description = "Protium production rate [s-1]", bounds = (0.0, Inf)]
    G_alpha(t), [guess = 0.0, description = "Alpha particles from D-T production rate [s-1]", bounds = (0.0, Inf)]
    DNH(t), [guess = 1.0e20, description = "D(N_H) = d(zN_H)/dt, time derivative of H population [1/s]", bounds = (0.0, Inf)]
  end


  @equations begin
    # Eqs. for ni(t)
    # IMPLICIT EQUATIONS ----------------------------------------------------------------------------------------------------------
    G_H ~ wc*(nD*nD)*rR_DD_1 + wc*(nD*nHe)*rR_DHe3 + wc*(n*n)*rR_DD_2*(rR_THe3_1 + rR_THe3_3)
    G_alpha ~ wc*(nD*nT)*rR_DT
    DNH ~ Vₚₗₐₛₘₐ*D(nH)
    Vₚₗₐₛₘₐ * D(nH) ~ a.xH*a.G + b.xH*b.G + wc*(nD*nD)*rR_DD_1    + wc*(nD*nHe)*rR_DHe3 + wc*(nT*nHe)*(rR_THe3_1 + rR_THe3_3)
    Vₚₗₐₛₘₐ * D(nD) ~ a.xD*a.G + b.xD*b.G + wc*(nT*nHe)*rR_THe3_2 - wc*(nD*nHe)*rR_DHe3 - wc*(nD*nT)*rR_DT   - wc*2*(nD*nD)*(rR_DD_1 + rR_DD_2)
    Vₚₗₐₛₘₐ * D(nT) ~ a.xT*a.G + b.xT*b.G + wc*(nD*nD)*rR_DD_1    - wc*(nD*nT)*rR_DT    - wc*2*(nT*nT)*rR_TT - wc*(nT*nHe)*(rR_THe3_1 + rR_THe3_2 + rR_THe3_3) - λT*nT*Vₚₗₐₛₘₐ
    #Vₚₗₐₛₘₐ * D(nHe) ~ a.xHe*a.G + b.xHe*b.G   + wc*(nD*nD)*rR_DD_2    - wc*(nD*nHe)*rR_DHe3 - wc*(nT*nHe)*(rR_THe3_1 + rR_THe3_2 + rR_THe3_3)
    #Vₚₗₐₛₘₐ * D(nHe45) ~ -nHe45*Vₚₗₐₛₘₐ/τp       + wc*(nD*nT)*rR_DT      + wc*(nD*nHe)*rR_DHe3 + wc*(nT*nT)*rR_TT + wc*(nT*nHe)*(rR_THe3_1 + rR_THe3_2 + rR_THe3_3)
    # Simplification: He is only one, species (consider He3, He4 and He5 as one)
    Vₚₗₐₛₘₐ * D(nHe) ~ a.xHe*a.G + b.xHe*b.G + wc*(nD*nD)*rR_DD_2 + wc*(nD*nT)*rR_DT + wc*(nT*nT)*rR_TT
    Vₚₗₐₛₘₐ * D(nXe) ~ a.xXe*a.G + b.xXe*b.G
    Vₚₗₐₛₘₐ * D(nI)  ~ a.xI*a.G  + b.xI*b.G
    # -----------------------------------------------------------------------------------------------------------------------------
    # EXPLICIT EQUATIONS ----------------------------------------------------------------------------------------------------------
    #G_H ~ wc*(n*n)*rR_DD_1    + wc*(n*n)*rR_DD_2*rR_DHe3      + wc*(n*n)*rR_DD_2*(rR_THe3_1 + rR_THe3_3)
    #G_alpha ~ wc*(n*n)*rR_DT
    #DNH ~ Vₚₗₐₛₘₐ * D(nH)
    #Vₚₗₐₛₘₐ * D(nH)  ~ a.xH*a.G  + b.xH*b.G    + wc*(n*n)*rR_DD_1    + wc*(n*n)*rR_DD_2*rR_DHe3      + wc*(n*n)*rR_DD_2*(rR_THe3_1 + rR_THe3_3)
    #Vₚₗₐₛₘₐ * D(nD)  ~ a.xD*a.G  + b.xD*b.G    - wc*n*n*rR_DT
    #Vₚₗₐₛₘₐ * D(nT)  ~ a.xT*a.G  + b.xT*b.G    - wc*n*n*rR_DT
    #Vₚₗₐₛₘₐ * D(nHe) ~ a.xHe*a.G + b.xHe*b.G   + wc*n*n*rR_DT
    #He45_prodRate  ~ 0.0
    # -----------------------------------------------------------------------------------------------------------------------------
    # Eqs. for b.xi
    b.xH*b.G ~ -nH*Vₚₗₐₛₘₐ/τp
    b.xD*b.G ~ -nD*Vₚₗₐₛₘₐ/τp
    b.xT*b.G ~ -nT*Vₚₗₐₛₘₐ/τp
    #b.xHe*b.G ~ -(nHe + nHe45) * Vₚₗₐₛₘₐ / τp # SIMPLIFICATION: consider He3, He4 and He5 as one
    b.xHe*b.G ~ -nHe*Vₚₗₐₛₘₐ/τp_α # SIMPLIFICATION: consider He3, He4 and He5 as one
    b.xXe*b.G ~ -nXe*Vₚₗₐₛₘₐ/τp
    b.xI*b.G  ~ -nI*Vₚₗₐₛₘₐ/τp
    # Eq. for b.G
    #b.G ~ -(nH + nT + nD + nHe + nHe45) * Vₚₗₐₛₘₐ / τp # SIMPLIFICATION: consider He3, He4 and He5 as one
    b.G ~ -(nH + nT + nD + nXe + nI)*Vₚₗₐₛₘₐ/τp - nHe*Vₚₗₐₛₘₐ/τp_α # SIMPLIFICATION: consider He3, He4 and He5 as one

    b.dummy ~ 0.0
  end
end



# ============================================================================ #
#tag FirstWall 
# ============================================================================ #
@doc raw"""
    FirstWall(;name)

First Wall (FW) outgassing model. Emits a pure protium (H) stream at a rate set by the
`H_in` input signal. The outlet composition is fixed to `xH = 1` with all other species at
zero. This block represents H outgassing from plasma-facing components under ion bombardment.

# Connectors:
- `b`: `AFPOutlet` [outlet] — H outgassing stream
- `H_in`: `RealInput` [input] — H outgassing rate [``\mathrm{s^{-1}}``]

# Parameters:
- `OG_H`: [``\mathrm{s^{-1}}``] nominal H outgassing rate (informational; actual rate is driven by `H_in`)
"""
@mtkmodel FirstWall begin
  @components begin
    b = AFPOutlet() # [outlet]
    H_in = RealInput() # [input]
  end

  @parameters begin
    OG_H = 0.0, [description = "H outgassing from First Wall [s-1]", input = true, bounds = (0.0, Inf)]
  end

  @variables begin end

  @equations begin
    b.xH ~ 1.0
    b.xD ~ 0.0
    b.xT ~ 0.0
    b.xHe ~ 0.0
    b.xXe ~ 0.0
    b.xI ~ 0.0
    b.G ~ -H_in.u

    b.dummy ~ 0.0
  end
end



# ============================================================================ #
#tag TorusMixer 
# ============================================================================ #
@doc raw"""
    TorusMixer(;name)

Mixer volume with two inlet ports 'a1' and 'a2' and one outlet port 'b', a fluid flow and 6 species composition.

# States:
- 'dm(t)': [mol s-1] The inlet molar flow rate
- 'xH(t)':  [-] The inlet molar fraction of H2
- 'xD(t)':  [-] The inlet molar fraction of D2
- 'xT(t)':  [-] The inlet molar fraction of T2
- 'xHe(t)': [-] The inlet molar fraction of He
- 'xXe(t)': [-] The inlet molar fraction of He
- 'xI(t)':  [-] The inlet molar fraction of Ne

# Connectors:
- 'a_FW': inlet fluid port
- 'a_plasma': inlet fluid port
- 'b' : outlet fluid port
"""
@mtkmodel TorusMixer begin
  @components begin
    a_FW = RealInput() # [input]
    a_plasma = AFPInlet() # [inlet]
    a_DC = AFPInlet() # [inlet]
    b = AFPOutlet() # [outlet]
  end

  @variables begin
    G(t), [guess = 0.0, description = "Flow Rate [mol s-1]", bounds = (0.0, Inf)]
    xH(t), [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xD(t), [guess = 0.5, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xT(t), [guess = 0.5, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xHe(t), [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xXe(t), [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xI(t), [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
  end

  @equations begin
    xH * G  ~ a_FW.u + a_plasma.xH  * a_plasma.G + a_DC.xH  * a_DC.G
    xD * G  ~          a_plasma.xD  * a_plasma.G + a_DC.xD  * a_DC.G
    xT * G  ~          a_plasma.xT  * a_plasma.G + a_DC.xT  * a_DC.G
    xHe * G ~          a_plasma.xHe * a_plasma.G + a_DC.xHe * a_DC.G
    xXe * G ~          a_plasma.xXe * a_plasma.G + a_DC.xXe * a_DC.G
    xI * G  ~          a_plasma.xI  * a_plasma.G + a_DC.xI  * a_DC.G
    G ~ a_FW.u + a_plasma.G * (a_plasma.xH + a_plasma.xD + a_plasma.xT + a_plasma.xHe + a_plasma.xXe + a_plasma.xI) + a_DC.G * (a_DC.xH + a_DC.xD + a_DC.xT + a_DC.xHe + a_DC.xXe + a_DC.xI)
    #
    b.xH ~ xH
    b.xD ~ xD
    b.xT ~ xT
    b.xHe ~ xHe
    b.xXe ~ xXe
    b.xI ~ xI
    b.G ~ -G

    b.dummy ~ 0.0
  end
end



# ============================================================================ #
#tag TorusMixer2
# ============================================================================ #
@doc raw"""
    TorusMixer_withFuelEff(;name)

Mixer volume with two inlet ports 'a1' and 'a2' and one outlet port 'b', a fluid flow and 6 species composition.

# States:
- 'dm(t)': [mol s-1] The inlet molar flow rate
- 'xH(t)':  [-] The inlet molar fraction of H2
- 'xD(t)':  [-] The inlet molar fraction of D2
- 'xT(t)':  [-] The inlet molar fraction of T2
- 'xHe(t)': [-] The inlet molar fraction of He
- 'xXe(t)': [-] The inlet molar fraction of He
- 'xI(t)':  [-] The inlet molar fraction of Ne

# Connectors:
- 'a_FW': inlet fluid port
- 'a_plasma': inlet fluid port
- 'b' : outlet fluid port
"""
@mtkmodel TorusMixer_withFuelEff begin
  @components begin
    a_FW = RealInput() # [input]
    a_plasma = AFPInlet() # [inlet]
    a_DC = AFPInlet() # [inlet]
    a_nf = AFPInlet() # [inlet] fuel not delivered to the plasma
    b = AFPOutlet() # [outlet]
  end

  @variables begin
    G(t),   [guess = 0.0, description = "Flow Rate [mol s-1]", bounds = (0.0, Inf)]
    xH(t),  [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xD(t),  [guess = 0.5, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xT(t),  [guess = 0.5, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xHe(t), [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xXe(t), [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xI(t),  [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
  end

  @equations begin
    xH * G  ~ a_FW.u + a_plasma.xH *a_plasma.G  +  a_nf.xH *a_nf.G  +  a_DC.xH *a_DC.G
    xD * G  ~          a_plasma.xD *a_plasma.G  +  a_nf.xD *a_nf.G  +  a_DC.xD *a_DC.G
    xT * G  ~          a_plasma.xT *a_plasma.G  +  a_nf.xT *a_nf.G  +  a_DC.xT *a_DC.G
    xHe * G ~          a_plasma.xHe*a_plasma.G  +  a_nf.xHe*a_nf.G  +  a_DC.xHe*a_DC.G
    xXe * G ~          a_plasma.xXe*a_plasma.G  +  a_nf.xXe*a_nf.G  +  a_DC.xXe*a_DC.G
    xI * G  ~          a_plasma.xI *a_plasma.G  +  a_nf.xI *a_nf.G  +  a_DC.xI *a_DC.G
    G ~ a_FW.u + a_plasma.G*(a_plasma.xH + a_plasma.xD + a_plasma.xT + a_plasma.xHe + a_plasma.xXe + a_plasma.xI) + a_nf.G*(a_nf.xH + a_nf.xD + a_nf.xT + a_nf.xHe + a_nf.xXe + a_nf.xI) + a_DC.G*(a_DC.xH + a_DC.xD + a_DC.xT + a_DC.xHe + a_DC.xXe + a_DC.xI)
    #
    b.xH ~ xH
    b.xD ~ xD
    b.xT ~ xT
    b.xHe ~ xHe
    b.xXe ~ xXe
    b.xI ~ xI
    b.G ~ -G

    b.dummy ~ 0.0
  end
end



# ============================================================================ #
#tag MFP 
# ============================================================================ #
@doc raw"""
    MFP(;name)

Metal Foil Pump (MFP) model. It has one inlet port 'a' and two outlet ports 'b1' and 'b2'. It has a fluid flow and 6 species composition.
A portion eta of the inlet Q (H,D,T) is sent to the outlet 'b1', the rest + (1-eta)Q is sent to the outlet port 'b2'.

# States:
- 'dm(t)':  [mol s-1] The inlet molar flow rate
- 'xH(t)':  [-] The inlet molar fraction of H2
- 'xD(t)':  [-] The inlet molar fraction of D2
- 'xT(t)':  [-] The inlet molar fraction of T2
- 'xHe(t)': [-] The inlet molar fraction of He
- 'xXe(t)': [-] The inlet molar fraction of He
- 'xI(t)':  [-] The inlet molar fraction of I

# Connectors:
- 'a': inlet fluid port
- 'b1' : outlet fluid port
- 'b2' : outlet fluid port
"""
@mtkmodel MFP begin
  @components begin
    a = AFPInlet() # [inlet]
    b1 = AFPOutlet() # [outlet]
    b2 = AFPOutlet() # [outlet]
  end

  @variables begin end

  @parameters begin
    eta = 1.0, [description = "Metal Foil Pump hydrogen (Q) separation efficiency [-]", input = true, bounds = (0.0, 1.0)]
  end

  @equations begin
    0.0 ~ a.xH * a.G + b1.xH * b1.G + b2.xH * b2.G
    0.0 ~ a.xD * a.G + b1.xD * b1.G + b2.xD * b2.G
    0.0 ~ a.xT * a.G + b1.xT * b1.G + b2.xT * b2.G
    0.0 ~ a.xHe * a.G + b1.xHe * b1.G + b2.xHe * b2.G
    0.0 ~ a.xXe * a.G + b1.xXe * b1.G + b2.xXe * b2.G
    0.0 ~ a.xI * a.G + b1.xI * b1.G + b2.xI * b2.G
    #
    0.0 ~ eta * a.xH * a.G + b1.xH * b1.G
    0.0 ~ eta * a.xD * a.G + b1.xD * b1.G
    0.0 ~ eta * a.xT * a.G + b1.xT * b1.G
    0.0 ~ b1.xHe * b1.G
    0.0 ~ b1.xXe * b1.G
    0.0 ~ b1.xI * b1.G
    b1.G ~ -a.G * eta * (a.xH + a.xD + a.xT)
    b2.G ~ -a.G * (1 - eta) * (a.xH + a.xD + a.xT) - a.G * (a.xHe + a.xXe + a.xI)

    b2.dummy ~ 0.0
    b1.dummy ~ 0.0
  end
end



# ============================================================================ #
#tag Splitter 
# ============================================================================ #
@doc raw"""
    Splitter(;name)

Isocompositional flow splitter. Divides the inlet stream into two outlet streams `b1` and
`b2` with identical compositions. A fraction `eta` of the inlet flow exits through `b1` and
the remainder `(1 - eta)` through `b2`.

# States:
- `G(t)`: [``\mathrm{s^{-1}}``] total inlet flow rate
- `xH(t)`, `xD(t)`, `xT(t)`, `xHe(t)`, `xXe(t)`, `xI(t)`: [-] mole fractions

# Connectors:
- `a`: `AFPInlet` [inlet]
- `b1`: `AFPOutlet` [outlet] — fraction `eta` of inlet flow
- `b2`: `AFPOutlet` [outlet] — fraction `(1 - eta)` of inlet flow

# Parameters:
- `eta`: [-] fraction of flow directed to `b1`
"""
@mtkmodel Splitter begin
  @components begin
    a = AFPInlet() # [inlet]
    b1 = AFPOutlet() # [outlet]
    b2 = AFPOutlet() # [outlet]
  end

  @variables begin
    G(t), [guess = 0.0, description = "Flow Rate [mol s-1]", bounds = (0.0, Inf)]
    xH(t), [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xD(t), [guess = 0.5, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xT(t), [guess = 0.5, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xHe(t), [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xXe(t), [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xI(t), [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
  end

  @parameters begin
    eta = 0.8, [description = "Splitter separation fraction [-]. Eta goes to b1, (1-eta) goes to b2.", input = true, bounds = (0.0, 1.0)]
  end

  @equations begin
    xH ~ a.xH
    xD ~ a.xD
    xT ~ a.xT
    xHe ~ a.xHe
    xXe ~ a.xXe
    xI ~ a.xI
    G ~ a.G
    #
    b1.xH ~ xH
    b1.xD ~ xD
    b1.xT ~ xT
    b1.xHe ~ xHe
    b1.xXe ~ xXe
    b1.xI ~ xI
    b1.G ~ -eta * G
    #
    b2.xH ~ xH
    b2.xD ~ xD
    b2.xT ~ xT
    b2.xHe ~ xHe
    b2.xXe ~ xXe
    b2.xI ~ xI
    b2.G ~ -(1 - eta) * G

    b2.dummy ~ 0.0
    b1.dummy ~ 0.0
  end
end






# ============================================================================ #
#tag LimitedSplitter 
# ============================================================================ #
@doc raw"""
    LimitedSplitter(;name)

Isocompositional flow splitter with a maximum flow rate constraint on `b1`. Nominally a
fraction `eta` of the inlet flow exits through `b1`, capped at `Gmax`. Any excess flow
above the cap is redirected to `b2`, so that ``\dot{G}_{b1} = \min(\eta G, G_{\max})`` and
``\dot{G}_{b2} = G - \dot{G}_{b1}``. Both outlets carry the same composition as the inlet.

# States:
- `G(t)`: [``\mathrm{s^{-1}}``] total inlet flow rate
- `xH(t)`, `xD(t)`, `xT(t)`, `xHe(t)`, `xXe(t)`, `xI(t)`: [-] mole fractions

# Connectors:
- `a`: `AFPInlet` [inlet]
- `b1`: `AFPOutlet` [outlet] — limited fraction of inlet flow
- `b2`: `AFPOutlet` [outlet] — remainder of inlet flow

# Parameters:
- `eta`: [-] nominal fraction of flow directed to `b1`
- `Gmax`: [``\mathrm{s^{-1}}``] maximum allowable flow rate through `b1`
"""
@mtkmodel LimitedSplitter begin
  @components begin
    a = AFPInlet() # [inlet]
    b1 = AFPOutlet() # [outlet]
    b2 = AFPOutlet() # [outlet]
  end

  @variables begin
    G(t), [guess = 0.0, description = "Flow Rate [mol s-1]", bounds = (0.0, Inf)]
    xH(t), [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xD(t), [guess = 0.5, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xT(t), [guess = 0.5, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xHe(t), [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xXe(t), [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xI(t), [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
  end

  @parameters begin
    eta = 0.8, [description = "Splitter separation fraction [-]. Eta goes to b1, (1-eta) goes to b2.", input = true, bounds = (0.0, 1.0)]
    Gmax = 0.0, [description = "Maximum flow rate for b1 [s-1].", input = true, bounds = (0.0, 1.0e25)]
  end

  @equations begin
    xH ~ a.xH
    xD ~ a.xD
    xT ~ a.xT
    xHe ~ a.xHe
    xXe ~ a.xXe
    xI ~ a.xI
    G ~ a.G
    #
    b1.xH ~ xH
    b1.xD ~ xD
    b1.xT ~ xT
    b1.xHe ~ xHe
    b1.xXe ~ xXe
    b1.xI ~ xI
    #b1.G ~ -eta * G
    b1.G ~ ifelse(eta*G - Gmax >= 0.0, -Gmax, -eta*G) # b1.G can be at most Gmax
    #
    b2.xH ~ xH
    b2.xD ~ xD
    b2.xT ~ xT
    b2.xHe ~ xHe
    b2.xXe ~ xXe
    b2.xI ~ xI
    #b2.G ~ -(1 - eta) * G
    b2.G ~ -(G + b1.G)

    b1.dummy ~ 0.0
    b2.dummy ~ 0.0
  end
end



# ============================================================================ #
#tag Pump 
# ============================================================================ #
@doc raw"""
    Pump(;name)

Idealised vacuum pump model. Holds an internal inventory and empties it through outlet `b`
with a species-independent characteristic pumping time `τ`. All species (H, D, T, He, Xe, I)
are pumped at the same rate ``N_i / \tau``.

Tritium decay via ``{}^3\mathrm{H} \to {}^3\mathrm{He}`` is included in the `D(NT)` balance.

# States:
- `NH(t)`: [``\mathrm{s^{-1}}``] total H inventory
- `ND(t)`: [``\mathrm{s^{-1}}``] total D inventory
- `NT(t)`: [``\mathrm{s^{-1}}``] total T inventory
- `NHe(t)`: [``\mathrm{s^{-1}}``] total He inventory
- `NXe(t)`: [``\mathrm{s^{-1}}``] total Xe inventory
- `NI(t)`: [``\mathrm{s^{-1}}``] total inert inventory

# Connectors:
- `a`: `AFPInlet` [inlet]
- `b`: `AFPOutlet` [outlet]

# Parameters:
- `τ`: [s] characteristic pumping time
"""
@mtkmodel Pump begin
  @components begin
    a = AFPInlet() # [inlet]
    b = AFPOutlet() # [outlet]
  end


  @parameters begin
    τ = 40.0, [description = "characteristic time [s]", input = true, bounds = (0.0, 100000.0)]
  end


  @variables begin
    #
    NH(t), [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    ND(t), [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NT(t), [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NHe(t), [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NXe(t), [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NI(t), [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
  end


  @equations begin
    D(NH) ~ a.xH * a.G + b.xH * b.G
    D(ND) ~ a.xD * a.G + b.xD * b.G
    D(NT) ~ a.xT * a.G + b.xT * b.G - λT*NT
    D(NHe) ~ a.xHe * a.G + b.xHe * b.G
    D(NXe) ~ a.xXe * a.G + b.xXe * b.G
    D(NI) ~ a.xI * a.G + b.xI * b.G
    #
    b.xH * b.G ~ -NH / τ
    b.xD * b.G ~ -ND / τ
    b.xT * b.G ~ -NT / τ
    b.xHe * b.G ~ -NHe / τ
    b.xXe * b.G ~ -NXe / τ
    b.xI * b.G ~ -NI / τ
    b.G ~ -(NH + ND + NT + NHe + NXe + NI) / τ

    b.dummy ~ 0.0
  end
end




# ============================================================================ #
#tag PumpHDTO
# ============================================================================ #
@doc raw"""
    PumpHDTO(;name)

Idealised vacuum pump model for H, D, T, O mixtures. Holds an internal inventory and empties
it through outlet `b` with a species-independent characteristic pumping time `τ`. All species
are pumped at the rate ``N_i / \tau``.

Tritium decay via ``{}^3\mathrm{H} \to {}^3\mathrm{He}`` is included in the `D(NT)` balance.

# States:
- `NH(t)`: [``\mathrm{s^{-1}}``] total H inventory
- `ND(t)`: [``\mathrm{s^{-1}}``] total D inventory
- `NT(t)`: [``\mathrm{s^{-1}}``] total T inventory
- `NO(t)`: [``\mathrm{s^{-1}}``] total O inventory

# Connectors:
- `a`: `AFPInlet_HDTO` [inlet]
- `b`: `AFPOutlet_HDTO` [outlet]

# Parameters:
- `τ`: [s] characteristic pumping time
"""
@mtkmodel PumpHDTO begin
  @components begin
    a = AFPInlet_HDTO() # [inlet]
    b = AFPOutlet_HDTO() # [outlet]
  end


  @parameters begin
    τ = 40.0, [description = "characteristic time [s]", input = true, bounds = (0.0, 100000.0)]
  end


  @variables begin
    #
    NH(t), [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    ND(t), [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NT(t), [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NO(t), [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
  end


  @equations begin
    D(NH) ~ a.xH * a.G + b.xH * b.G
    D(ND) ~ a.xD * a.G + b.xD * b.G
    D(NT) ~ a.xT * a.G + b.xT * b.G - λT*NT
    D(NO) ~ a.xO * a.G + b.xO * b.G
    #
    b.xH * b.G ~ -NH / τ
    b.xD * b.G ~ -ND / τ
    b.xT * b.G ~ -NT / τ
    b.xO * b.G ~ -NO / τ
    b.G ~ -(NH + ND + NT + NO) / τ

    b.dummy ~ 0.0
  end
end


# ============================================================================ #
#tag WallHPumping 
# ============================================================================ #
@doc raw"""
    WallHPumping(;name)

First-wall H-pumping model. Takes a scalar H source rate from the `a` input signal and
partitions it into a pumped fraction directed to the fuel cycle via outlet `b1` (pure H
stream) and an un-pumped fraction reported as a scalar signal on `b2`. The parameter `eta`
controls the pumping fraction.

# Connectors:
- `a`: `RealInput` [inlet] — total H wall-pumping rate [``\mathrm{s^{-1}}``]
- `b1`: `AFPOutlet` [outlet] — pumped H stream (pure H composition)
- `b2`: `RealOutput` [outlet] — un-pumped H rate [``\mathrm{s^{-1}}``]

# Parameters:
- `eta`: [-] pumped fraction of the total H rate
"""
@mtkmodel WallHPumping begin
  @components begin
    a = RealInput() # [inlet]
    b1 = AFPOutlet() # [outlet]
    b2 = RealOutput() # [outlet]
  end

  @variables begin
  end

  @parameters begin
    eta = 1.0, [description = "Pumping fraction [-].", input = true, bounds = (0.0, 1.0)]
  end

  @equations begin
    b1.xH  ~ 1.0
    b1.xD  ~ 0.0
    b1.xT  ~ 0.0
    b1.xHe ~ 0.0
    b1.xXe ~ 0.0
    b1.xI  ~ 0.0
    b1.G   ~ - a.u
    #
    b2.u ~ (1 - eta) * a.u

    b1.dummy ~ 0.0
  end
end


# ============================================================================ #
#tag FluidOpenBC 
# ============================================================================ #
@doc raw"""
    FluidOpenBC(;name)

Simple BC for a fluid port such that the mass flow 'dm' is non-zero.
It is the opposite from an un-connected port.

# States:
- 'dm(t)': [mol s-1] The molar flow rate passing through the 'a' port
- 'xH(t)': [-] The molar fraction of H2 passing thorugh the 'a' port
- 'xD(t)': [-] The molar fraction of D2 passing thorugh the 'a' port
- 'xT(t)': [-] The molar fraction of T2 passing thorugh the 'a' port
- 'xHe(t)': [-] The molar fraction of He passing thorugh the 'a' port
- 'xXe(t)': [-] The molar fraction of He passing thorugh the 'a' port
- 'xI(t)': [-] The molar fraction of Inert passing thorugh the 'a' port

# Connectors:
- 'a': inlet fluid port
"""
@mtkmodel FluidOpenBC begin
  @components begin
    a = AFPInlet()
  end

  @variables begin
    G(t),   [guess = 0.0, description = "Flow Rate [mol s-1]"]
    xH(t),  [guess = 0.0, description = "Molar fraction [-]", output = true, bounds = (0.0, 1.0)]
    xD(t),  [guess = 0.0, description = "Molar fraction [-]", output = true, bounds = (0.0, 1.0)]
    xT(t),  [guess = 0.0, description = "Molar fraction [-]", output = true, bounds = (0.0, 1.0)]
    xHe(t), [guess = 0.0, description = "Molar fraction [-]", output = true, bounds = (0.0, 1.0)]
    xXe(t), [guess = 0.0, description = "Molar fraction [-]", output = true, bounds = (0.0, 1.0)]
    xI(t),  [guess = 0.0, description = "Molar fraction [-]", output = true, bounds = (0.0, 1.0)]
  end

  @equations begin
    a.G ~ G
    a.xH ~ xH
    a.xD ~ xD
    a.xT ~ xT
    a.xHe ~ xHe
    a.xXe ~ xXe
    a.xI ~ xI
  end
end




# ============================================================================ #
#tag FluidOpenBC_withO 
# ============================================================================ #
@doc raw"""
    FluidOpenBC_withO(;name)

Simple BC for a fluid port such that the mass flow 'dm' is non-zero.
It is the opposite from an un-connected port.
THe ports contain also Oxygen 'O'

# States:
- 'dm(t)': [mol s-1] The molar flow rate passing through the 'a' port
- 'xH(t)': [-] The molar fraction of H2 passing thorugh the 'a' port
- 'xD(t)': [-] The molar fraction of D2 passing thorugh the 'a' port
- 'xT(t)': [-] The molar fraction of T2 passing thorugh the 'a' port
- 'xHe(t)': [-] The molar fraction of He passing thorugh the 'a' port
- 'xXe(t)': [-] The molar fraction of He passing thorugh the 'a' port
- 'xI(t)': [-] The molar fraction of Inert passing thorugh the 'a' port
- 'xO(t)': [-] The molar fraction of O passing thorugh the 'a' port

# Connectors:
- 'a': inlet fluid port
"""
@mtkmodel FluidOpenBC_withO begin
  @components begin
    a = AFPInlet_O()
  end

  @variables begin
    G(t),   [guess = 0.0, description = "Flow Rate [mol s-1]"]
    xH(t),  [guess = 0.0, description = "Molar fraction [-]", output = true, bounds = (0.0, 1.0)]
    xD(t),  [guess = 0.0, description = "Molar fraction [-]", output = true, bounds = (0.0, 1.0)]
    xT(t),  [guess = 0.0, description = "Molar fraction [-]", output = true, bounds = (0.0, 1.0)]
    xHe(t), [guess = 0.0, description = "Molar fraction [-]", output = true, bounds = (0.0, 1.0)]
    xXe(t), [guess = 0.0, description = "Molar fraction [-]", output = true, bounds = (0.0, 1.0)]
    xI(t),  [guess = 0.0, description = "Molar fraction [-]", output = true, bounds = (0.0, 1.0)]
    xO(t),  [guess = 0.0, description = "Molar fraction [-]", output = true, bounds = (0.0, 1.0)]
  end

  @equations begin
    a.G ~ G
    a.xH ~ xH
    a.xD ~ xD
    a.xT ~ xT
    a.xHe ~ xHe
    a.xXe ~ xXe
    a.xI ~ xI
    a.xO ~ xO
  end
end




# ============================================================================ #
#tag FluidOpenBC_OCN 
# ============================================================================ #
@doc raw"""
    FluidOpenBC_OCN(;name)

Simple BC for a fluid port such that the mass flow 'dm' is non-zero.
It is the opposite from an un-connected port.
THe ports contain also Oxygen 'O', 'C', 'N'

# States:
- 'dm(t)': [mol s-1] The molar flow rate passing through the 'a' port
- 'xH(t)': [-] The molar fraction of H2 passing thorugh the 'a' port
- 'xD(t)': [-] The molar fraction of D2 passing thorugh the 'a' port
- 'xT(t)': [-] The molar fraction of T2 passing thorugh the 'a' port
- 'xHe(t)': [-] The molar fraction of He passing thorugh the 'a' port
- 'xXe(t)': [-] The molar fraction of He passing thorugh the 'a' port
- 'xI(t)': [-] The molar fraction of Inert passing thorugh the 'a' port
- 'xO(t)': [-] The molar fraction of O passing thorugh the 'a' port
- 'xC(t)': [-] The molar fraction of C passing thorugh the 'a' port
- 'xN(t)': [-] The molar fraction of N passing thorugh the 'a' port

# Connectors:
- 'a': inlet fluid port
"""
@mtkmodel FluidOpenBC_OCN begin
  @components begin
    a = AFPInlet_OCN()
  end

  @variables begin
    G(t),   [guess = 0.0, description = "Flow Rate [mol s-1]"]
    xH(t),  [guess = 0.0, description = "Molar fraction [-]", output = true, bounds = (0.0, 1.0)]
    xD(t),  [guess = 0.0, description = "Molar fraction [-]", output = true, bounds = (0.0, 1.0)]
    xT(t),  [guess = 0.0, description = "Molar fraction [-]", output = true, bounds = (0.0, 1.0)]
    xHe(t), [guess = 0.0, description = "Molar fraction [-]", output = true, bounds = (0.0, 1.0)]
    xXe(t), [guess = 0.0, description = "Molar fraction [-]", output = true, bounds = (0.0, 1.0)]
    xI(t),  [guess = 0.0, description = "Molar fraction [-]", output = true, bounds = (0.0, 1.0)]
    xO(t),  [guess = 0.0, description = "Molar fraction [-]", output = true, bounds = (0.0, 1.0)]
    xC(t),  [guess = 0.0, description = "Molar fraction [-]", output = true, bounds = (0.0, 1.0)]
    xN(t),  [guess = 0.0, description = "Molar fraction [-]", output = true, bounds = (0.0, 1.0)]
  end

  @equations begin
    a.G ~ G
    a.xH ~ xH
    a.xD ~ xD
    a.xT ~ xT
    a.xHe ~ xHe
    a.xXe ~ xXe
    a.xI ~ xI
    a.xO ~ xO
    a.xC ~ xC
    a.xN ~ xN
  end
end




# ============================================================================ #
#tag FluidOpenBC_HDTO 
# ============================================================================ #
@doc raw"""
    FluidOpenBC_HDTO(;name)

Simple BC for a fluid port such that the mass flow 'dm' is non-zero.
It is the opposite from an un-connected port.

# States:
- 'dm(t)': [mol s-1] The flow rate passing through the 'a' port
- 'xH(t)': [-] The fraction of H passing thorugh the 'a' port
- 'xD(t)': [-] The fraction of D passing thorugh the 'a' port
- 'xT(t)': [-] The fraction of T passing thorugh the 'a' port
- 'xO(t)': [-] The fraction of O passing thorugh the 'a' port

# Connectors:
- 'a': inlet fluid port
"""
@mtkmodel FluidOpenBC_HDTO begin
  @components begin
    a = AFPInlet_HDTO()
  end

  @variables begin
    G(t),   [guess = 0.0, description = "Flow Rate [mol s-1]"]
    xH(t),  [guess = 0.0, description = "Molar fraction [-]", output = true, bounds = (0.0, 1.0)]
    xD(t),  [guess = 0.0, description = "Molar fraction [-]", output = true, bounds = (0.0, 1.0)]
    xT(t),  [guess = 0.0, description = "Molar fraction [-]", output = true, bounds = (0.0, 1.0)]
    xO(t),  [guess = 0.0, description = "Molar fraction [-]", output = true, bounds = (0.0, 1.0)]
  end

  @equations begin
    a.G ~ G
    a.xH ~ xH
    a.xD ~ xD
    a.xT ~ xT
    a.xO ~ xO
  end
end




# ============================================================================ #
#tag MixerTwoWay 
# ============================================================================ #
@doc raw"""
    MixerTwoWay(;name)

Mixer volume with two inlet ports 'a1' and 'a2' and one outlet port 'b', a fluid flow and 6 species composition.

# States:
- 'dm(t)':  [mol s-1] The inlet molar flow rate
- 'xH(t)':  [-] The inlet molar fraction of H2
- 'xD(t)':  [-] The inlet molar fraction of D2
- 'xT(t)':  [-] The inlet molar fraction of T2
- 'xHe(t)': [-] The inlet molar fraction of He
- 'xXe(t)': [-] The inlet molar fraction of He
- 'xI(t)':  [-] The inlet molar fraction of Inert

# Connectors:
- 'a1': inlet fluid port
- 'a2': inlet fluid port
- 'b' : outlet fluid port
"""
@mtkmodel MixerTwoWay begin
  @components begin
    a1 = AFPInlet() # [inlet]
    a2 = AFPInlet() # [inlet]
    b = AFPOutlet() # [outlet]
  end

  @variables begin
    G(t), [guess = 0.0, description = "Flow Rate [mol s-1]", bounds = (0.0, Inf)]
    xH(t), [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xD(t), [guess = 0.5, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xT(t), [guess = 0.5, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xHe(t), [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xXe(t), [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xI(t), [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
  end

  @equations begin
    xH * G  ~ a1.xH * a1.G + a2.xH * a2.G
    xD * G  ~ a1.xD * a1.G + a2.xD * a2.G
    xT * G  ~ a1.xT * a1.G + a2.xT * a2.G
    xHe * G ~ a1.xHe * a1.G + a2.xHe * a2.G
    xXe * G ~ a1.xXe * a1.G + a2.xXe * a2.G
    xI * G  ~ a1.xI * a1.G + a2.xI * a2.G
    G       ~ a1.G * (a1.xH + a1.xD + a1.xT + a1.xHe + a1.xXe + a1.xI) + a2.G * (a2.xH + a2.xD + a2.xT + a2.xHe + a2.xXe + a2.xI)
    #
    b.xH  ~ xH
    b.xD  ~ xD
    b.xT  ~ xT
    b.xHe ~ xHe
    b.xXe ~ xXe
    b.xI  ~ xI
    b.G   ~ -G

    b.dummy ~ 0.0
  end
end




# ============================================================================ #
#tag MixerTwoWay_OCN 
# ============================================================================ #
@doc raw"""
    MixerTwoWayOCN(;name)

Mixer volume with two inlet ports 'a1' and 'a2' and one outlet port 'b', a fluid flow and 6 species composition.

# States:
- 'dm(t)':  [mol s-1] The inlet molar flow rate
- 'xH(t)':  [-] The inlet molar fraction of H2
- 'xD(t)':  [-] The inlet molar fraction of D2
- 'xT(t)':  [-] The inlet molar fraction of T2
- 'xHe(t)': [-] The inlet molar fraction of He
- 'xXe(t)': [-] The inlet molar fraction of He
- 'xI(t)':  [-] The inlet molar fraction of Inert
- 'xO(t)':  [-] The inlet molar fraction of O
- 'xC(t)':  [-] The inlet molar fraction of C
- 'xN(t)':  [-] The inlet molar fraction of N

# Connectors:
- 'a1': inlet fluid port
- 'a2': inlet fluid port
- 'b' : outlet fluid port
"""
@mtkmodel MixerTwoWayOCN begin
  @components begin
    a1 = AFPInlet_OCN() # [inlet]
    a2 = AFPInlet_OCN() # [inlet]
    b = AFPOutlet_OCN() # [outlet]
  end

  @variables begin
    G(t),   [guess = 0.0, description = "Flow Rate [mol s-1]", bounds = (0.0, Inf)]
    xH(t),  [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xD(t),  [guess = 0.5, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xT(t),  [guess = 0.5, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xHe(t), [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xXe(t), [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xI(t),  [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xO(t),  [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xC(t),  [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xN(t),  [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
  end

  @equations begin
    xH *G ~ a1.xH*a1.G  + a2.xH*a2.G 
    xD *G ~ a1.xD*a1.G  + a2.xD*a2.G 
    xT *G ~ a1.xT*a1.G  + a2.xT*a2.G 
    xHe*G ~ a1.xHe*a1.G + a2.xHe*a2.G
    xXe*G ~ a1.xXe*a1.G + a2.xXe*a2.G
    xI *G ~ a1.xI*a1.G  + a2.xI*a2.G 
    xO *G ~ a1.xO*a1.G  + a2.xO*a2.G 
    xC *G ~ a1.xC*a1.G  + a2.xC*a2.G 
    xN *G ~ a1.xN*a1.G  + a2.xN*a2.G 
    G     ~ a1.G*(a1.xH + a1.xD + a1.xT + a1.xHe + a1.xXe + a1.xI + a1.xO + a1.xC + a1.xN) + a2.G*(a2.xH + a2.xD + a2.xT + a2.xHe + a2.xXe + a2.xI + a2.xO + a2.xC + a2.xN)
    #
    b.xH  ~ xH
    b.xD  ~ xD
    b.xT  ~ xT
    b.xHe ~ xHe
    b.xXe ~ xXe
    b.xI  ~ xI
    b.xO  ~ xO
    b.xC  ~ xC
    b.xN  ~ xN
    b.G   ~ -G

    b.dummy ~ 0.0
  end
end




# ============================================================================ #
#tag MixerThreeWay
# ============================================================================ #
@doc raw"""
    MixerThreeWay(;name)

Mixer volume with two inlet ports 'a1' and 'a2' and one outlet port 'b', a fluid flow and 6 species composition.

# States:
- 'dm(t)':  [mol s-1] The inlet molar flow rate
- 'xH(t)':  [-] The inlet molar fraction of H2
- 'xD(t)':  [-] The inlet molar fraction of D2
- 'xT(t)':  [-] The inlet molar fraction of T2
- 'xHe(t)': [-] The inlet molar fraction of He
- 'xXe(t)': [-] The inlet molar fraction of He
- 'xI(t)':  [-] The inlet molar fraction of Inert

# Connectors:
- 'a1': inlet fluid port
- 'a2': inlet fluid port
- 'a3': inlet fluid port
- 'b' : outlet fluid port
"""
@mtkmodel MixerThreeWay begin
  @components begin
    a1 = AFPInlet() # [inlet]
    a2 = AFPInlet() # [inlet]
    a3 = AFPInlet() # [inlet]
    b = AFPOutlet() # [outlet]
  end

  @variables begin
    G(t),   [guess = 0.0, description = "Flow Rate [mol s-1]", bounds = (0.0, Inf)]
    xH(t),  [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xD(t),  [guess = 0.5, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xT(t),  [guess = 0.5, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xHe(t), [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xXe(t), [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xI(t),  [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
  end

  @equations begin
    xH *G ~ a1.xH*a1.G  + a2.xH*a2.G   + a3.xH*a3.G 
    xD *G ~ a1.xD*a1.G  + a2.xD*a2.G   + a3.xD*a3.G 
    xT *G ~ a1.xT*a1.G  + a2.xT*a2.G   + a3.xT*a3.G 
    xHe*G ~ a1.xHe*a1.G + a2.xHe*a2.G  + a3.xHe*a3.G
    xXe*G ~ a1.xXe*a1.G + a2.xXe*a2.G  + a3.xXe*a3.G
    xI *G ~ a1.xI*a1.G  + a2.xI*a2.G   + a3.xI*a3.G 
    G     ~ a1.G*(a1.xH + a1.xD + a1.xT + a1.xHe + a1.xXe + a1.xI) + a2.G*(a2.xH + a2.xD + a2.xT + a2.xHe + a2.xXe + a2.xI) + a3.G*(a3.xH + a3.xD + a3.xT + a3.xHe + a3.xXe + a3.xI)
    #
    b.xH  ~ xH
    b.xD  ~ xD
    b.xT  ~ xT
    b.xHe ~ xHe
    b.xXe ~ xXe
    b.xI  ~ xI
    b.G   ~ -G

    b.dummy ~ 0.0
  end
end




# ============================================================================ #
#tag MixerThreeWayOCN 
# ============================================================================ #
@doc raw"""
    MixerThreeWayOCN(;name)

Mixer volume with two inlet ports 'a1' and 'a2' and one outlet port 'b', a fluid flow and 6 species composition.

# States:
- 'dm(t)':  [mol s-1] The inlet molar flow rate
- 'xH(t)':  [-] The inlet molar fraction of H2
- 'xD(t)':  [-] The inlet molar fraction of D2
- 'xT(t)':  [-] The inlet molar fraction of T2
- 'xHe(t)': [-] The inlet molar fraction of He
- 'xXe(t)': [-] The inlet molar fraction of He
- 'xI(t)':  [-] The inlet molar fraction of Inert
- 'xO(t)':  [-] The inlet molar fraction of O
- 'xC(t)':  [-] The inlet molar fraction of C
- 'xN(t)':  [-] The inlet molar fraction of N

# Connectors:
- 'a1': inlet fluid port
- 'a2': inlet fluid port
- 'a3': inlet fluid port
- 'b' : outlet fluid port
"""
@mtkmodel MixerThreeWayOCN begin
  @components begin
    a1 = AFPInlet_OCN() # [inlet]
    a2 = AFPInlet_OCN() # [inlet]
    a3 = AFPInlet_OCN() # [inlet]
    b = AFPOutlet_OCN() # [outlet]
  end

  @variables begin
    G(t),   [guess = 0.0, description = "Flow Rate [mol s-1]", bounds = (0.0, Inf)]
    xH(t),  [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xD(t),  [guess = 0.5, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xT(t),  [guess = 0.5, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xHe(t), [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xXe(t), [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xI(t),  [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xO(t),  [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xC(t),  [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
    xN(t),  [guess = 0.0, description = "Molar fraction [-]", bounds = (0.0, 1.0)]
  end

  @equations begin
    xH *G ~ a1.xH*a1.G  + a2.xH*a2.G   + a3.xH*a3.G 
    xD *G ~ a1.xD*a1.G  + a2.xD*a2.G   + a3.xD*a3.G 
    xT *G ~ a1.xT*a1.G  + a2.xT*a2.G   + a3.xT*a3.G 
    xHe*G ~ a1.xHe*a1.G + a2.xHe*a2.G  + a3.xHe*a3.G
    xXe*G ~ a1.xXe*a1.G + a2.xXe*a2.G  + a3.xXe*a3.G
    xI *G ~ a1.xI*a1.G  + a2.xI*a2.G   + a3.xI*a3.G 
    xO *G ~ a1.xO*a1.G  + a2.xO*a2.G   + a3.xO*a3.G 
    xC *G ~ a1.xC*a1.G  + a2.xC*a2.G   + a3.xC*a3.G 
    xN *G ~ a1.xN*a1.G  + a2.xN*a2.G   + a3.xN*a3.G 
    G     ~ a1.G*(a1.xH + a1.xD + a1.xT + a1.xHe + a1.xXe + a1.xI + a1.xO + a1.xC + a1.xN) + a2.G*(a2.xH + a2.xD + a2.xT + a2.xHe + a2.xXe + a2.xI + a2.xO + a2.xC + a2.xN) + a3.G*(a3.xH + a3.xD + a3.xT + a3.xHe + a3.xXe + a3.xI + a3.xO + a3.xC + a3.xN)
    #
    b.xH  ~ xH
    b.xD  ~ xD
    b.xT  ~ xT
    b.xHe ~ xHe
    b.xXe ~ xXe
    b.xI  ~ xI
    b.xO  ~ xO
    b.xC  ~ xC
    b.xN  ~ xN
    b.G   ~ -G

    b.dummy ~ 0.0
  end
end




# ============================================================================ #
#tag IRPR 
# ============================================================================ #
@doc raw"""
    IRPR(;name)

Isotope-selective Recirculation and Purification Reactor (IRPR) model. Receives a mixed
hydrogen-isotope stream at inlet `a` and separates it into a Q-enriched product stream
(`b1`) and a waste/recycle stream (`b2`). He and inerts are fully removed to `b2`. The
separation of H, D, T exploits their mass difference via the efficiency factor `eta`:
the purified fractions are ``\eta\, N_H``, ``(\eta/\sqrt{2})\, N_D``, and
``(\eta/2)\, N_T`` respectively.

Tritium decay via ``{}^3\mathrm{H} \to {}^3\mathrm{He}`` is included in the `D(NT)` balance.

# States:
- `NH(t)`: [``\mathrm{s^{-1}}``] total H inventory
- `ND(t)`: [``\mathrm{s^{-1}}``] total D inventory
- `NT(t)`: [``\mathrm{s^{-1}}``] total T inventory
- `NHe(t)`: [``\mathrm{s^{-1}}``] total He inventory
- `NXe(t)`: [``\mathrm{s^{-1}}``] total Xe inventory
- `NI(t)`: [``\mathrm{s^{-1}}``] total inert inventory

# Connectors:
- `a`: `AFPInlet` [inlet]
- `b1`: `AFPOutlet` [outlet] — Q-enriched purified stream
- `b2`: `AFPOutlet` [outlet] — waste/recycle stream (He, Xe, I, and residual Q)

# Parameters:
- `τ`: [s] characteristic residence time
- `eta`: [-] H separation efficiency; D and T efficiencies are scaled by ``1/\sqrt{2}`` and ``1/2`` respectively
"""
@mtkmodel IRPR begin
  @components begin
    a = AFPInlet() # [inlet]
    b1 = AFPOutlet() # [outlet]
    b2 = AFPOutlet() # [outlet]
  end

  @variables begin
    NH(t), [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    ND(t), [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NT(t), [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NHe(t), [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NXe(t), [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
    NI(t), [guess = 1.0e-20, description = "total # of particles [-]", bounds = (0.0, Inf)]
  end

  @parameters begin
    τ = 1.0, [description = "Characteristic time of the system [s]", input = true, bounds = (0.0, 100000.0)]
    eta = 0.95, [description = "Pellet Injector Helium separation efficiency [-]", input = true, bounds = (0.0, 1.0)]
  end

  @equations begin
    D(NH)  ~ a.xH * a.G + b1.xH * b1.G + b2.xH * b2.G
    D(ND)  ~ a.xD * a.G + b1.xD * b1.G + b2.xD * b2.G
    D(NT)  ~ a.xT * a.G + b1.xT * b1.G + b2.xT * b2.G - λT*NT
    D(NHe) ~ a.xHe * a.G + b1.xHe * b1.G + b2.xHe * b2.G
    D(NXe) ~ a.xXe * a.G + b1.xXe * b1.G + b2.xXe * b2.G
    D(NI)  ~ a.xI * a.G + b1.xI * b1.G + b2.xI * b2.G
    #D(NH) + D(ND) + D(NT) + D(NHe) ~ a.G*(a.xH + a.xD + a.xT + a.xHe) + b1.G*(b1.xH + b1.xD + b1.xT + b1.xHe) + b2.G*(b2.xH + b2.xD + b2.xT + b2.xHe)
    #
    b2.xI * b2.G  ~ -NI / τ
    b2.xXe * b2.G ~ -NXe / τ
    b2.xHe * b2.G ~ -NHe / τ
    b2.xD * b2.G  ~ -(1 - eta/1.414) * ND / τ
    b2.xT * b2.G  ~ -(1 - eta/2.0) * NT / τ
    b2.xH * b2.G  ~ -(1 - eta) * NH / τ
    b2.G ~ -(NI + NXe + NHe + ND*(1 - eta/1.414) + NT*(1 - eta/2.0) + NH*(1 - eta)) / τ
    #
    b1.xI ~ 0.0
    b1.xXe ~ 0.0
    b1.xHe ~ 0.0
    b1.xD*b1.G ~ -(eta/1.414) * ND / τ
    b1.xT*b1.G ~ -(eta/2.0) * NT / τ
    b1.xH*b1.G ~ -eta * NH / τ
    b1.G ~ -eta*NH/τ -(eta/1.414)*ND/τ -(eta/2.0)*NT/τ

    b1.dummy ~ 0.0
    b2.dummy ~ 0.0
  end
end

