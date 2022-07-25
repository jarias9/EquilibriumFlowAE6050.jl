function equilibriumTP(model::Air5s, T, P)
    # Function for evaluating equilibrium properties as
    # a function temperature and pressure (T,P)
    # Function is overloaded depending on reaction model

    # Model: ["N2","O2","NO","N","O"]
    
    # Define Constants
    @unpack h, k, R, c, hc, Na = constants()
    R_x = uconvert.(u"J/kg/K", R ./ model.weights) # [J/kg-K]

    # Calculate equilibrium 
    X = calculateEquilibrium(model::Air5s, T, P)

    # Generate functions for evaluating ∂Q/∂T
    # ["N2","O2","NO","N","O"]
    logQ(T) = [log(ustrip(Q_N2(T)));
                log(ustrip(Q_O2(T)));
                log(ustrip(Q_NO(T)));
                log(ustrip(Q_N(T)));
                log(ustrip(Q_O(T)))] # remove units for using ForwardDiff

    # Generate derivative of log(Q) and
    # calculate internal energy and enthalpy
    dlogQ = ForwardDiff.derivative(logQ, ustrip(T))u"K^-1" 
    e_ =  R_x .* T.^2 .* dlogQ # [J/kg] Sensible internal energy per species
    h_ =  R_x .* T.^2 .* dlogQ .+ R_x.*T # [J/kg] Sensible enthalpy per species
    H_ =  R .* T.^2 .* dlogQ .+ R.*T # [J/mol] Sensible enthalpy per species

    # Enthalpy of Formation
    Δh = uconvert.(u"J/kg", model.deltaH ./ model.weights) # [J/kg]
    ΔH = model.deltaH

    # Calculate mass fractions
    mix_mass = sum(X*u"mol" .* (model.weights*u"g/mol"))
    Y = (X*u"mol" .* (model.weights*u"g/mol")) ./ mix_mass

    # Calculate mixture enthalpy
    h_mix = sum(Y .* h_) + sum(Y .* Δh)
    h0_mix = sum(X .* H_) + sum(X .* ΔH)

    # Calculate entropy
    s_ = equilibriumentropy(model, X, h_mix, T, P)

    # Calculate density
    R_mix = uconvert(u"J/kg/K", R / sum(model.weights .* X))
    rho = uconvert(u"kg/m^3", P/ (R_mix*T))

    return h_mix, rho, s_, X
end

function equilibriumTP(model::Air7s, T, P)
    # Define Constants
    @unpack h, k, R, c, hc, Na = constants()
    R_x = uconvert.(u"J/kg/K", R ./ model.weights) # [J/kg-K]

    # Calculate equilibrium 
    X = calculateEquilibrium(model::Air7s, T, P)

    # Generate functions for evaluating ∂Q/∂T
    # ["N2","O2","NO","N","O", "NO+", "e-"]
    logQ(T) = [log(ustrip(Q_N2(T)));
                log(ustrip(Q_O2(T)));
                log(ustrip(Q_NO(T)));
                log(ustrip(Q_N(T)));
                log(ustrip(Q_O(T)));
                log(ustrip(Q_NO_plus(T)));
                log(ustrip(Q_e_minus(T)))] # remove units for using ForwardDiff

    dlogQ = ForwardDiff.derivative(logQ, ustrip(T))u"K^-1" 
    e_ =  R_x .* T.^2 .* dlogQ # [J/kg] Sensible internal energy per species
    h_ =  R_x .* T.^2 .* dlogQ .+ R_x.*T # [J/kg] Sensible enthalpy per species
    H_ =  R .* T.^2 .* dlogQ .+ R.*T # [J/mol] Sensible enthalpy per species

    # Enthalpy of Formation
    Δh = uconvert.(u"J/kg", model.deltaH ./ model.weights) # [J/kg]
    ΔH = model.deltaH

    # Calculate mass fractions
    mix_mass = sum(X*u"mol" .* (model.weights*u"g/mol"))
    Y = (X*u"mol" .* (model.weights*u"g/mol")) ./ mix_mass

    # Calculate mixture enthalpy
    h_mix = sum(Y .* h_) + sum(Y .* Δh)
    h0_mix = sum(X .* H_) + sum(X .* ΔH)

    # Calculate entropy
    s_ = equilibriumentropy(model, X, h_mix, T, P)

    # Calculate density
    R_mix = uconvert(u"J/kg/K", R / sum(model.weights .* X))
    rho = uconvert(u"kg/m^3", P/ (R_mix*T))

    return h_mix, rho, s_, X
end

function equilibriumTP(model::Air11s, T, P)
    # Define Constants
    @unpack h, k, R, c, hc, Na = constants()
    R_x = uconvert.(u"J/kg/K", R ./ model.weights) # [J/kg-K]

    # Calculate equilibrium 
    X = calculateEquilibrium(model::Air11s, T, P)

    # Generate functions for evaluating ∂Q/∂T
    # ["N2", "O2", "NO", "O", "N", "NO+", "e-", "N2+", "O2+", "N+", "O+"]
    logQ(T) = [log(ustrip(Q_N2(T)));
                log(ustrip(Q_O2(T)));
                log(ustrip(Q_NO(T)));
                log(ustrip(Q_N(T)));
                log(ustrip(Q_O(T)));
                log(ustrip(Q_NO_plus(T)));
                log(ustrip(Q_e_minus(T)));
                log(ustrip(Q_N2_plus(T)));
                log(ustrip(Q_O2_plus(T)));
                log(ustrip(Q_N_plus(T)));
                log(ustrip(Q_O_plus(T)))] # remove units for using ForwardDiff

    dlogQ = ForwardDiff.derivative(logQ, ustrip(T))u"K^-1" 
    e_ =  R_x .* T.^2 .* dlogQ # [J/kg] Sensible internal energy per species
    h_ =  R_x .* T.^2 .* dlogQ .+ R_x.*T # [J/kg] Sensible enthalpy per species
    H_ =  R .* T.^2 .* dlogQ .+ R.*T # [J/mol] Sensible enthalpy per species

    # Enthalpy of Formation
    Δh = uconvert.(u"J/kg", model.deltaH ./ model.weights) # [J/kg]
    ΔH = model.deltaH

    # Calculate mass fractions
    mix_mass = sum(X*u"mol" .* (model.weights*u"g/mol"))
    Y = (X*u"mol" .* (model.weights*u"g/mol")) ./ mix_mass

    # Calculate mixture enthalpy
    h_mix = sum(Y .* h_) + sum(Y .* Δh)
    h0_mix = sum(X .* H_) + sum(X .* ΔH)

    # Calculate entropy
    s_ = equilibriumentropy(model, X, h_mix, T, P)

    # Calculate density
    R_mix = uconvert(u"J/kg/K", R / sum(model.weights .* X))
    rho = uconvert(u"kg/m^3", P/ (R_mix*T))

    return h_mix, rho, s_, X
end

function equilibriumTP(model::Air13s, T, P)
    # Define Constants
    @unpack h, k, R, c, hc, Na = constants()
    R_x = uconvert.(u"J/kg/K", R ./ model.weights) # [J/kg-K]

    # Calculate equilibrium 
    X = calculateEquilibrium(model::Air13s, T, P)

    # Generate functions for evaluating ∂Q/∂T
    # ["N2", "O2", "NO", "O", "N", "NO+", "e-", "N2+", "O2+", "N+", "O+", "Ar", "Ar+"]
    logQ(T) = [log(ustrip(Q_N2(T)));
                log(ustrip(Q_O2(T)));
                log(ustrip(Q_NO(T)));
                log(ustrip(Q_N(T)));
                log(ustrip(Q_O(T)));
                log(ustrip(Q_NO_plus(T)));
                log(ustrip(Q_e_minus(T)));
                log(ustrip(Q_N2_plus(T)));
                log(ustrip(Q_O2_plus(T)));
                log(ustrip(Q_N_plus(T)));
                log(ustrip(Q_O_plus(T)));
                log(ustrip(Q_Ar(T)));
                log(ustrip(Q_Ar_plus(T)))] # remove units for using ForwardDiff

    dlogQ = ForwardDiff.derivative(logQ, ustrip(T))u"K^-1" 
    e_ =  R_x .* T.^2 .* dlogQ # [J/kg] Sensible internal energy per species
    h_ =  R_x .* T.^2 .* dlogQ .+ R_x.*T # [J/kg] Sensible enthalpy per species
    H_ =  R .* T.^2 .* dlogQ .+ R.*T # [J/mol] Sensible enthalpy per species

    # Enthalpy of Formation
    Δh = uconvert.(u"J/kg", model.deltaH ./ model.weights) # [J/kg]
    ΔH = model.deltaH

    # Calculate mass fractions
    mix_mass = sum(X*u"mol" .* (model.weights*u"g/mol"))
    Y = (X*u"mol" .* (model.weights*u"g/mol")) ./ mix_mass

    # Calculate mixture enthalpy
    h_mix = sum(Y .* h_) + sum(Y .* Δh)
    h0_mix = sum(X .* H_) + sum(X .* ΔH)

    # Calculate entropy
    s_ = equilibriumentropy(model, X, h_mix, T, P)

    # Calculate density
    R_mix = uconvert(u"J/kg/K", R / sum(model.weights .* X))
    rho = uconvert(u"kg/m^3", P/ (R_mix*T))

    return h_mix, rho, s_, X
end

function equilibriumTP(model::CO2_6s, T, P)
    # Define Constants
    @unpack h, k, R, c, hc, Na = constants()
    R_x = uconvert.(u"J/kg/K", R ./ model.weights) # [J/kg-K]

    # Calculate equilibrium 
    X = calculateEquilibrium(model::CO2_6s, T, P)

    # Generate functions for evaluating ∂Q/∂T
    # ["CO2", "CO", "N2", "N", "O2", "O"]
    logQ(T) = [log(ustrip(Q_CO2(T)));
                log(ustrip(Q_CO(T)));
                log(ustrip(Q_N2(T)));
                log(ustrip(Q_N(T)));
                log(ustrip(Q_O2(T)));
                log(ustrip(Q_O(T)))] # remove units for using ForwardDiff

    dlogQ = ForwardDiff.derivative(logQ, ustrip(T))u"K^-1" 
    e_ =  R_x .* T.^2 .* dlogQ # [J/kg] Sensible internal energy per species
    h_ =  R_x .* T.^2 .* dlogQ .+ R_x.*T # [J/kg] Sensible enthalpy per species
    H_ =  R .* T.^2 .* dlogQ .+ R.*T # [J/mol] Sensible enthalpy per species

    # Enthalpy of Formation
    Δh = uconvert.(u"J/kg", model.deltaH ./ model.weights) # [J/kg]
    ΔH = model.deltaH

    # Calculate mass fractions
    mix_mass = sum(X*u"mol" .* (model.weights*u"g/mol"))
    Y = (X*u"mol" .* (model.weights*u"g/mol")) ./ mix_mass

    # Calculate mixture enthalpy
    h_mix = sum(Y .* h_) + sum(Y .* Δh)
    h0_mix = sum(X .* H_) + sum(X .* ΔH)

    # Calculate entropy
    s_ = equilibriumentropy(model, X, h_mix, T, P)

    # Calculate density
    R_mix = uconvert(u"J/kg/K", R / sum(model.weights .* X))
    rho = uconvert(u"kg/m^3", P/ (R_mix*T))

    return h_mix, rho, s_, X
end

