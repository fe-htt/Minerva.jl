"""
Functions related to fusion reactions and the plasma component.
"""


@doc raw"""
    get_nProfile(n0, np, a, ap, αn, βn)

returns the plasma density profile

```math
n(x, n_p, n_0) = n_p + (n_0 - n_p) \cdot [1 - (\frac{a}{ap})^2 \cdot x^2]^{\alpha n}
```

# Arguments
- `n0`: density at the plasma center [m-3]
- `np`: density at the pedestel [m-3]
- `a`: coefficient [-]
- `ap`: coefficient [-] 
- `αn`: coefficient [-]
- `βn`: coefficient [-]
"""
function get_nProfile(n0, np, a, ap, αn, βn)
  return n(x) = @. np + (n0 - np) * (1.0 - ((a / ap)^βn) * ((x)^βn))^αn
end


@doc raw"""
    get_TProfile(n0, np, a, ap, αn, βn)

returns the plasma temperature profile.

```math
T(x, T_p, T_0) = T_p + (T_0 - T_p) \cdot [1 - (\frac{a}{ap})^{\beta t} \cdot x^{\beta t}]^{\alpha t}
```

# Arguments
- `T0`: temperature at the plasma center [keV]
- `Tp`: temperature at the pedestel [keV]
- `a`: coefficient [-]
- `ap`: coefficient [-]
- `αt`: coefficient [-]
- `βt`: coefficient [-]
"""
function get_TProfile(T0, Tp, a, ap, αt, βt)
  return T(x) = @. Tp + (T0 - Tp) * (1.0 - ((a / ap)^βt) * ((x)^βt))^αt
end


@doc raw"""
    get_nAveraged(n0, np, a, ap, αn, βn)

Returns the line-average plasma density (average over the radius).
Takes as parameters the coefficients of a profile with the shape:
```math
n(x, n_p, n_0) = n_p + (n_0 - n_p) \cdot [1 - (\frac{a}{ap})^2 \cdot x^2]^{\alpha n}
```

# Arguments
- `n0`: density at the plasma center [m-3]
- `np`: density at the pedestel [m-3]
- `a`: coefficient [-]
- `ap`: coefficient [-]
- `αn`: coefficient [-]
- `βn`: coefficient [-]
"""
function get_nAveraged(n0, np, a, ap, αn, βn)
  n = get_nProfile(n0, np, a, ap, αn, βn)
  n_integrated = quadgk(n, 0, 1)[1]
  return n_integrated / 1
end


@doc raw"""
    get_nAveraged_radius(n0, np, a, ap, αn, βn)

Returns the line-average plasma density (average over the radius).
Takes as parameters the coefficients of a profile with the shape:
```math
n(x, n_p, n_0) = n_p + (n_0 - n_p) \cdot [1 - (\frac{a}{ap})^2 \cdot x^2]^{\alpha n}
```

# Arguments
- `n0`: density at the plasma center [m-3]
- `np`: density at the pedestel [m-3]
- `a`: coefficient [-]
- `ap`: coefficient [-]
- `αn`: coefficient [-]
- `βn`: coefficient [-]
"""
function get_nAveraged_radius(n0, np, a, ap, αn, βn) # average over the radius!
  n = get_nProfile(n0, np, a, ap, αn, βn)
  n_integrated = quadgk(n, 0, 1)[1]
  return n_integrated / 1
end


@doc raw"""
    get_nAveraged_disk(n0, np, a, ap, αn, βn)

Returns the area-average plasma density (average over the cross-section).
Takes as parameters the coefficients of a profile with the shape:
```math
n(x, n_p, n_0) = n_p + (n_0 - n_p) \cdot [1 - (\frac{a}{ap})^2 \cdot x^2]^{\alpha n}
```

# Arguments
- `n0`: density at the plasma center [m-3]
- `np`: density at the pedestel [m-3]
- `a`: coefficient [-]
- `ap`: coefficient [-]
- `αn`: coefficient [-]
- `βn`: coefficient [-]
"""
function get_nAveraged_disk(n0, np, a, ap, αn, βn) # average over the disk!
  n = get_nProfile(n0, np, a, ap, αn, βn)
  tmp(x) = @. n(x) * x
  n_integrated = (2/1^2) * quadgk(tmp, 0, 1)[1]
  return n_integrated
end


@doc raw"""
    get_nAveraged_radius_2018(n)

Returns the line-average plasma density (average over the radius).

# Arguments
- `n`: function representing the density profile. n(x), where 0 <= x <= 1 [m-3]
"""
function get_nAveraged_radius_2018(n) # average over the radius!
  n_integrated = quadgk(n, 0, 1, rtol=1e-20, atol=1e-30)[1] / 1
  return n_integrated
end


@doc raw"""
    get_nAveraged_disk_2018(n)

Returns the area-average plasma density (average over the cross-section).

# Arguments
- `n`: function representing the density profile. n(x), where 0 <= x <= 1 [m-3]
"""
function get_nAveraged_disk_2018(n::Function) # average over the disk!
  tmp(x) = @. n(x) * x
  n_averaged = (2/1^2) * quadgk(tmp, 0, 1)[1]
  return n_averaged
end


@doc raw"""
    get_nAveraged2(n0, np, a, ap, αn, βn, xaz)

Returns the line-average plasma density (average over the radius) over the Active Zone.
```math
n(x, n_p, n_0) = n_p + (n_0 - n_p) \cdot [1 - (\frac{a}{ap})^2 \cdot x^2]^{\alpha n}
```

# Arguments
- `n0`: density at the plasma center [m-3]
- `np`: density at the pedestel [m-3]
- `a`: coefficient [-]
- `ap`: coefficient [-]
- `αn`: coefficient [-]
- `βn`: coefficient [-]
- 'xaz': active zone normalized radius [-]
"""
function get_nAveraged2(n0, np, a, ap, αn, βn, xaz)
  n = get_nProfile(n0, np, a, ap, αn, βn)
  n_integrated = quadgk(n, 0, xaz)[1]
  return n_integrated / xaz
end



@doc raw"""
    get_nAvg_powerCorrectionFactor(n0, np, a, ap, αn, βn, xaz)

When computing fusion reactions over pre-integrated reactions rates, a correction factor needs to be used to account for the use of averaged densities.
This functions computes and returns this correction factor for computing reactions rates over the active zone.

It takes as input a density profile with the following shape:
```math
n(x, n_p, n_0) = n_p + (n_0 - n_p) \cdot [1 - (\frac{a}{ap})^2 \cdot x^2]^{\alpha n}
```

# Arguments
- `n0`: density at the plasma center [m-3]
- `np`: density at the pedestel [m-3]
- `a`: coefficient [-]
- `ap`: coefficient [-]
- `αn`: coefficient [-]
- `βn`: coefficient [-]
- 'xaz': active zone normalized radius [-]
"""
function get_nAvg_powerCorrectionFactor(n0, np, a, ap, αn, βn, xaz)
  n = get_nProfile(n0, np, a, ap, αn, βn)
  n_integrated_entire_plasma = quadgk(n, 0, 1)[1] / 1
  n_integrated_active_zone = quadgk(n, 0, xaz)[1] / xaz
  return (n_integrated_active_zone / n_integrated_entire_plasma)^2
end


@doc raw"""
    get_nAvg_radius_powerCorrectionFactor(n0, np, a, ap, αn, βn, xaz)

When computing fusion reactions over pre-integrated reactions rates, a correction factor needs to be used to account for the use of averaged densities.
This functions computes and returns this correction factor for computing reactions rates over the active zone.

It takes as input a density profile with the following shape:
```math
n(x, n_p, n_0) = n_p + (n_0 - n_p) \cdot [1 - (\frac{a}{ap})^2 \cdot x^2]^{\alpha n}
```
Use this function if working with line-averaged densities
# Arguments
- `n0`: density at the plasma center [m-3]
- `np`: density at the pedestel [m-3]
- `a`: coefficient [-]
- `ap`: coefficient [-]
- `αn`: coefficient [-]
- `βn`: coefficient [-]
- 'xaz': active zone normalized radius [-]
"""
function get_nAvg_radius_powerCorrectionFactor(n0, np, a, ap, αn, βn, xaz) # average over the radius, not the disk!
  n = get_nProfile(n0, np, a, ap, αn, βn)
  n_integrated_entire_plasma = quadgk(n, 0, 1)[1] / 1
  n_integrated_active_zone = quadgk(n, 0, xaz)[1] / xaz
  return (n_integrated_active_zone / n_integrated_entire_plasma)^2
end


@doc raw"""
    get_nAvg_disk_powerCorrectionFactor(n0, np, a, ap, αn, βn, xaz)

When computing fusion reactions over pre-integrated reactions rates, a correction factor needs to be used to account for the use of averaged densities.
This functions computes and returns this correction factor for computing reactions rates over the active zone.

It takes as input a density profile with the following shape:
```math
n(x, n_p, n_0) = n_p + (n_0 - n_p) \cdot [1 - (\frac{a}{ap})^2 \cdot x^2]^{\alpha n}
```
Use this function if working with area-averaged densities
# Arguments
- `n0`: density at the plasma center [m-3]
- `np`: density at the pedestel [m-3]
- `a`: coefficient [-]
- `ap`: coefficient [-]
- `αn`: coefficient [-]
- `βn`: coefficient [-]
- 'xaz': active zone normalized radius [-]
"""
function get_nAvg_disk_powerCorrectionFactor(n0, np, a, ap, αn, βn, xaz) # average over the disk!
  n = get_nProfile(n0, np, a, ap, αn, βn)
  tmp(x) = n(x) * x
  n_integrated_entire_plasma = (2/1^2)   * quadgk(tmp, 0, 1)[1]
  n_integrated_active_zone   = (2/xaz^2) * quadgk(tmp, 0, xaz)[1]
  return (n_integrated_active_zone / n_integrated_entire_plasma)^2
end


@doc raw"""
    get_nAvg_radius_powerCorrectionFactor_2018(n, xaz)

When computing fusion reactions over pre-integrated reactions rates, a correction factor needs to be used to account for the use of averaged densities.
This functions computes and returns this correction factor for computing reactions rates over the active zone.

It takes as input a density profile expressed as a function of the normalized plasma radius.
Use this function if working with line-averaged densities

# Arguments
- `n`: density profile as a function of coordinate x, with 0 <= x <= 1 [m-3]
- 'xaz': active zone normalized radius [-]
"""
function get_nAvg_radius_powerCorrectionFactor_2018(n, xaz) # average over the radius, not the disk!
  n_integrated_entire_plasma = quadgk(n, 0, 1)[1] / 1
  n_integrated_active_zone = quadgk(n, 0, xaz)[1] / xaz
  return (n_integrated_active_zone / n_integrated_entire_plasma)^2
end


@doc raw"""
    get_nAvg_disk_powerCorrectionFactor_2018(n, T, xaz)

When computing fusion reactions over pre-integrated reactions rates, a correction factor needs to be used to account for the use of averaged densities.
This functions computes and returns this correction factor for computing reactions rates over the active zone.

It takes as input a density and a temperature profiles expressed as a function of the normalized plasma radius.
Use this function if working with area-averaged densities

# Arguments
- `n`: density profile as a function of coordinate x, with 0 <= x <= 1
- `T`: temperature profile as a function of coordinate x, with 0 <= x <= 1
- 'xaz': active zone normalized radius
"""
function get_nAvg_disk_powerCorrectionFactor_2018(n, T, xaz) # average over the disk!
	u_DT   = Float64[3.35419e-7, 2.47828, 31.8505, 0.292999];
  σv_DT(x) = @. u_DT[1] * T(x)^(-u_DT[2]) * exp(-u_DT[3] * T(x)^(-u_DT[4])) * 1.0e-6           # [m3 s-1]
  intArg_DT(u) = @. σv_DT(u) * u
  intArg_tot(u) = @. n(u) * n(u) * σv_DT(u) * u
  intArg_navg(u) = @. n(u) * u

  # Volume averaged n (over entire plasma volume)
  navg = 2*quadgk(intArg_navg, 0, 1)[1]

  # What we should be doing every time: integrate n*n*σν over the AZ
  nnσν_AZ =quadgk(intArg_tot, 0, xaz, rtol=1e-20, atol=1e-30)[1]

  # What we actually do: <n>_vol * <n>_vol * integral(<σν>)_AZ
  nnσν_computed = navg * navg * quadgk(intArg_DT, 0, xaz, rtol=1e-20, atol=1e-30)[1]

  # Correction factor
  cf = nnσν_AZ / nnσν_computed

  return cf
end


@doc raw"""
    get_reactionRates(xaz::Float64, k::Float64, a::Float64, R::Float64, T::Function)

Integral of reaction rates over the active zone
Returns 8 values related to 8 reaction branches:
- `rR_DD_1`: [``\mathrm{m^6 s^{-1}}``]    D + D --(50%)--> T + p
- `rR_DD_2`: [``\mathrm{m^6 s^{-1}}``]    D + D --(50%)--> He3 + n
- `rR_DT`: [``\mathrm{m^6 s^{-1}}``]      D + T --> He4 + n
- `rR_DHe3`: [``\mathrm{m^6 s^{-1}}``]    D + He3 --> He4 + p
- `rR_TT`: [``\mathrm{m^6 s^{-1}}``]      T + T --> He4 + 2n + 11.3 MeV
- `rR_THe3_1`: [``\mathrm{m^6 s^{-1}}``]  T + He3 --(51%)--> He4 + p + n + 12.1 MeV
- `rR_THe3_2`: [``\mathrm{m^6 s^{-1}}``]  T + He3 --(43%)--> He4 + D
- `rR_THe3_3`: [``\mathrm{m^6 s^{-1}}``]  T + He3 --( 6%)--> He5 + p

# Arguments
- 'xaz': active zone normalized radius [-]
- 'k': plasma shape factor [-]
- 'a': minor radius [m]
- 'R': major radius [m]
- 'T': temperature profile (function of normalized radius) [keV]
"""
function get_reactionRates(xaz::Float64, k::Float64, a::Float64, R::Float64, T::Function)
  # Fit parameters for the reaction rates
	u_DD   = Float64[2.24642e-13, 0.739338, 21.1067, 0.306581];
	u_DT   = Float64[3.35419e-7, 2.47828, 31.8505, 0.292999];
	u_DHe3 = Float64[5.32595e-13 ,0.528691 ,31.4826 ,0.374262];
	u_TT   = Float64[1.01487e-13 ,0.730876 ,19.5003 ,0.283029];
	u_THe3 = Float64[2.20089e-10 ,1.2165 ,42.3469 ,0.30675];

  # Reaction rates for the various fusion reactions
  σv_DD(x) = @. u_DD[1] * T(x)^(-u_DD[2]) * exp(-u_DD[3] * T(x)^(-u_DD[4])) * 1.0e-6           # [m3 s-1]
  σv_DT(x) = @. u_DT[1] * T(x)^(-u_DT[2]) * exp(-u_DT[3] * T(x)^(-u_DT[4])) * 1.0e-6           # [m3 s-1]
  σv_DHe3(x) = @. u_DHe3[1] * T(x)^(-u_DHe3[2]) * exp(-u_DHe3[3] * T(x)^(-u_DHe3[4])) * 1.0e-6 # [m3 s-1]
  σv_TT(x) = @. u_TT[1] * T(x)^(-u_TT[2]) * exp(-u_TT[3] * T(x)^(-u_TT[4])) * 1.0e-6           # [m3 s-1]
  σv_THe3(x) = @. u_THe3[1] * T(x)^(-u_THe3[2]) * exp(-u_THe3[3] * T(x)^(-u_THe3[4])) * 1.0e-6 # [m3 s-1]

  # Integration Arguments
  intArg_DD_1(u) = @. 0.5 * σv_DD(u) * u
  intArg_DD_2(u) = @. 0.5 * σv_DD(u) * u
  intArg_DT(u) = @. σv_DT(u) * u
  intArg_DHe3(u) = @. σv_DHe3(u) * u
  intArg_TT(u) = @. σv_TT(u) * u
  intArg_THe3_1(u) = @. 0.51 * σv_THe3(u) * u
  intArg_THe3_2(u) = @. 0.43 * σv_THe3(u) * u
  intArg_THe3_3(u) = @. 0.06 * σv_THe3(u) * u

  # Integration
  Γ_DD_1   = @. 4 * k * a^2 * R * (pi^2) * quadgk(intArg_DD_1, 0, xaz, rtol=1e-20, atol=1e-30)[1]   # [m6 s-1] to be multiplied for n^2 to obtain the particles per second
  Γ_DD_2   = @. 4 * k * a^2 * R * (pi^2) * quadgk(intArg_DD_2, 0, xaz, rtol=1e-20, atol=1e-30)[1]   # [m6 s-1] to be multiplied for n^2 to obtain the particles per second
  Γ_DT     = @. 4 * k * a^2 * R * (pi^2) * quadgk(intArg_DT, 0, xaz, rtol=1e-20, atol=1e-30)[1]     # [m6 s-1] to be multiplied for n^2 to obtain the particles per second
  Γ_DHe3   = @. 4 * k * a^2 * R * (pi^2) * quadgk(intArg_DHe3, 0, xaz, rtol=1e-20, atol=1e-30)[1]   # [m6 s-1] to be multiplied for n^2 to obtain the particles per second
  Γ_TT     = @. 4 * k * a^2 * R * (pi^2) * quadgk(intArg_TT, 0, xaz, rtol=1e-20, atol=1e-30)[1]     # [m6 s-1] to be multiplied for n^2 to obtain the particles per second
  Γ_THe3_1 = @. 4 * k * a^2 * R * (pi^2) * quadgk(intArg_THe3_1, 0, xaz, rtol=1e-20, atol=1e-30)[1] # [m6 s-1] to be multiplied for n^2 to obtain the particles per second
  Γ_THe3_2 = @. 4 * k * a^2 * R * (pi^2) * quadgk(intArg_THe3_2, 0, xaz, rtol=1e-20, atol=1e-30)[1] # [m6 s-1] to be multiplied for n^2 to obtain the particles per second
  Γ_THe3_3 = @. 4 * k * a^2 * R * (pi^2) * quadgk(intArg_THe3_3, 0, xaz, rtol=1e-20, atol=1e-30)[1] # [m6 s-1] to be multiplied for n^2 to obtain the particles per second

  return Γ_DD_1, Γ_DD_2, Γ_DT, Γ_DHe3, Γ_TT, Γ_THe3_1, Γ_THe3_2, Γ_THe3_3
end


@doc raw"""
    get_DT_reactionRate(xaz::Float64, k::Float64, a::Float64, R::Float64, T::Function, n)

Return the D-T reaction rate [``\mathrm{#reactions s^{-1}}``]

# Arguments
- 'xaz': active zone normalized radius [-]
- 'k': plasma shape factor [-]
- 'a': minor radius [m]
- 'R': major radius [m]
- 'T': temperature profile (function of normalized radius) [keV]
- 'n': density profile (function of normalized radius) [m-3]

"""
function get_DT_reactionRate(xaz::Float64, k::Float64, a::Float64, R::Float64, T::Function, n)
  # Fit parameters for the reaction rates
  u_DT = Float64[3.35419e-7, 2.47828, 31.8505, 0.292999];

  # Reaction rates for the various fusion reactions
  σv_DT(x) = @. u_DT[1] * T(x)^(-u_DT[2]) * exp(-u_DT[3] * T(x)^(-u_DT[4])) * 1.0e-6           # [m3 s-1]

  # Integration Arguments
  intArg_DT(u) = @. n(u) * n(u) * σv_DT(u) * u

  # Integration
  Γ_DT     = @. 4 * k * a^2 * R * (pi^2) * quadgk(intArg_DT, 0, xaz, rtol=1e-20, atol=1e-30)[1]     # [m6 s-1] to be multiplied for n^2 to obtain the particles per second

  return Γ_DT
end



#=
# ============================================================================ #
function get_Q2O_vaporPressure(T::Float64, Q2O::String)
  # Dalla tesi di Jonas, che riporta le seguenti reference:
  # 1) W. A. Van Hook: Vapor pressures of the isotopic waters and ices. The Journal of Physical Chemistry, Volume 72, Issue 4, 1234, (1968). DOI: 10.1021/j100850a028.
  #    https://pubs.acs.org/doi/abs/10.1021/j100850a028
  # 2) D. R. Stull: Vapor Pressure of Pure Substances. Organic and Inorganic Compounds. Industrial & Engineering Chemistry Research, Volume 39, Issue 4, 517, (1947). DOI: 10.1021/ie50448a022
  #    https://pubs.acs.org/doi/abs/10.1021/ie50448a022

  if Q2O=="HDO"
    A = 26398.8
    B = -89.6065
    C = 0.075802
  else
    if Q2O=="HTO"
      A = 37813.2
      B = -136.751
      C = 0.124096
    else
      if Q2O=="D2O"
        A = 49314.9
        B = -164.266
        C = 0.140049components
      else
        if Q2O=="DTO"
          A = 59313.4
          B = -204.941
          C = 0.182686
        else
          if Q2O =="T2O"
            A = 68702.3
            B = -244.687
            C = 0.224388
          end
        end
      end
    endcomponents
  end

  A_H2O = 4.6543
  B_H2O = 1435.264
  C_H2O = -64.848
  p_H2O = 1e5*10^(A_H2O - B_H2O/(C_H2O + T)) # [Pa]

  if Q2O=="H2O"
    return p_H2O
  else
    p_Q2O_H2O = 1/exp(A/T^2 + B/T +C) # [-]
    return p_Q2O_H2O * p_H2O
  end
end




# ============================================================================ #
function get_Q2O_over_H2O_vaporPressure(T::Float64, Q2O::String)
  # Dalla tesi di Jonas, che riporta le seguenti reference:
  # 1) W. A. Van Hook: Vapor pressures of the isotopic waters and ices. The Journal of Physical Chemistry, Volume 72, Issue 4, 1234, (1968). DOI: 10.1021/j100850a028.
  #    https://pubs.acs.org/doi/abs/10.1021/j100850a028
  # 2) D. R. Stull: Vapor Pressure of Pure Substances. Organic and Inorganic Compounds. Industrial & Engineering Chemistry Research, Volume 39, Issue 4, 517, (1947). DOI: 10.1021/ie50448a022
  #    https://pubs.acs.org/doi/abs/10.1021/ie50448a022
  
  if Q2O=="H2O"
    A = 4.6543
    B = 1435.264
    C = -64.848
  else
    if Q2O=="HDO"
      A = 26398.8
      B = -89.6065
      C = 0.075802
    else
      if Q2O=="HTO"
        A = 37813.2
        B = -136.751
        C = 0.124096
      else
        if Q2O=="D2O"
          A = 49314.9
          B = -164.266
          C = 0.140049
        else
          if Q2O=="DTO"
            A = 59313.4
            B = -204.941
            C = 0.182686
          else # T2O
            A = 68702.3
            B = -244.687
            C = 0.224388
          end
        end
      end
    end
  end
  
  # Ratio of vapor pressure p_Q2O/p_H2O
  p_Q2O_H2O = 1/exp(A/T^2 + B/T +C)

  return p_Q2O_H2O

end
=#

function saveSim(all_comb, plot_quantities, sol, sens_name, current_directory)
  for (i,val) in enumerate(all_comb)
    local S = DataFrame()
    # Populate DataFrame with the simulation's results
    for mydata in plot_quantities
      insertcols!(S, Symbol(mydata) => sol[i][mydata])
    end
    # Rename the df column names to simpler names
    for n in names(S)
      rename!(S, n => replace(n, "(t)"=>"", "₊" => "."))
    end
    # Write to .csv file
    CSV.write(current_directory*"/"*sens_name*"/"*join(val,"-")*".csv", S)
  end
end
