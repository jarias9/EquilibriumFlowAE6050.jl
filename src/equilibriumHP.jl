function equilibriumHP(model, h_set, P)
    # Function for evaluating equilibrium properties as
    # a function enthalpy and pressure (h,P)
    # Temperature must be evaluated using a rootfinding approach

    # Generate equation to solve
    function f_enthalpy(model, h_set, T, P)
        h_, ~, ~, ~ = equilibriumTP(model, T, P)
        f = h_ - h_set
        return ustrip(f)
    end

    # Execute rootfinder
    f_sys(T) = f_enthalpy(model, h_set, T, P)
    T_start = 1000u"K"
    T_new = find_zero(f_sys, T_start, Roots.Order0())

    # Evaluate final composition at converged temperature
    h_mix, rho, s_, X = equilibriumTP(model, T_new, P)

    return T_new, rho,  s_, X

end

