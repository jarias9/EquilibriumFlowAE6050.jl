function equilibriumentropy(model::Air5s, X, h, T, P)
    # Function for evaluating equilibrium entropy as
    # a function of composition, enthalpy, 
    # temperature, and pressure (X, h, T, P)
    # Function is overloaded depending on reaction model

    # Model: ["N2","O2","NO","N","O"]

    # Import Constants
    @unpack R, k, Na = constants()

    # In per mole basis
    V = uconvert(u"m^3/mol", Na * k * T / P) # [m^3/mol]
    mw_mix = sum(model.weights .* X)
    R/mw_mix # "J/g-K"

    # Evaluate ln(Q/N)
    logQ(T) = [log(Q_N2(T) * V / (Na * X[1]));
                log(Q_O2(T) * V / (Na * X[2]));
                log(Q_NO(T) * V / (Na * X[3]));
                log(Q_N(T) * V / (Na * X[4]));
                log(Q_O(T) * V / (Na * X[5]))] 

    # Evaluate entropy 
    s_mix = uconvert(u"J/kg/K", R/mw_mix) * sum(X .* logQ(T)) + h/T # "J/kg-K"

    return s_mix
end

function equilibriumentropy(model::Air7s, X, h, T, P)
    # Function for evaluating equilibrium entropy as
    # a function of composition, enthalpy, 
    # temperature, and pressure (X, h, T, P)
    # Function is overloaded depending on reaction model

    # Model: ["N2","O2","NO","N","O", "NO+", "e-"]

    # Import Constants
    @unpack R, k, Na = constants()

    # In per mole basis
    V = uconvert(u"m^3/mol", Na * k * T / P) # [m^3/mol]
    mw_mix = sum(model.weights .* X)
    R/mw_mix # "J/g-K"

    # Evaluate ln(Q/N)
    logQ(T) = [log(Q_N2(T) * V / (Na * X[1]));
                log(Q_O2(T) * V / (Na * X[2]));
                log(Q_NO(T) * V / (Na * X[3]));
                log(Q_N(T) * V / (Na * X[4]));
                log(Q_O(T) * V / (Na * X[5]));
                log(Q_NO_plus(T) * V / (Na * X[6]));
                log(Q_e_minus(T) * V / (Na * X[7]))] 

    # Evaluate entropy 
    s_mix = uconvert(u"J/kg/K", R/mw_mix) * sum(X .* logQ(T)) + h/T # "J/kg-K"

    return s_mix
end

function equilibriumentropy(model::Air11s, X, h, T, P)
    # Function for evaluating equilibrium entropy as
    # a function of composition, enthalpy, 
    # temperature, and pressure (X, h, T, P)
    # Function is overloaded depending on reaction model

    # Model: # ["N2", "O2", "NO", "O", "N", "NO+", "e-", "N2+", "O2+", "N+", "O+"]

    # Import Constants
    @unpack R, k, Na = constants()

    # In per mole basis
    V = uconvert(u"m^3/mol", Na * k * T / P) # [m^3/mol]
    mw_mix = sum(model.weights .* X)
    R/mw_mix # "J/g-K"

    # Evaluate ln(Q/N)
    logQ(T) = [log(Q_N2(T) * V / (Na * X[1]));
                log(Q_O2(T) * V / (Na * X[2]));
                log(Q_NO(T) * V / (Na * X[3]));
                log(Q_N(T) * V / (Na * X[4]));
                log(Q_O(T) * V / (Na * X[5]));
                log(Q_NO_plus(T) * V / (Na * X[6]));
                log(Q_e_minus(T) * V / (Na * X[7]));
                log(Q_N2_plus(T) * V / (Na * X[8]));
                log(Q_O2_plus(T) * V / (Na * X[9]));
                log(Q_N_plus(T) * V / (Na * X[10]));
                log(Q_O_plus(T) * V / (Na * X[11]))] 

    # Evaluate entropy 
    s_mix = uconvert(u"J/kg/K", R/mw_mix) * sum(X .* logQ(T)) + h/T # "J/kg-K"

    return s_mix
end

function equilibriumentropy(model::Air13s, X, h, T, P)
    # Function for evaluating equilibrium entropy as
    # a function of composition, enthalpy, 
    # temperature, and pressure (X, h, T, P)
    # Function is overloaded depending on reaction model

    # Model: # ["N2", "O2", "NO", "O", "N", "NO+", "e-", "N2+", "O2+", "N+", "O+", "Ar", "Ar+"]

    # Import Constants
    @unpack R, k, Na = constants()

    # In per mole basis
    V = uconvert(u"m^3/mol", Na * k * T / P) # [m^3/mol]
    mw_mix = sum(model.weights .* X)
    R/mw_mix # "J/g-K"

    # Evaluate ln(Q/N)
    logQ(T) = [log(Q_N2(T) * V / (Na * X[1]));
                log(Q_O2(T) * V / (Na * X[2]));
                log(Q_NO(T) * V / (Na * X[3]));
                log(Q_N(T) * V / (Na * X[4]));
                log(Q_O(T) * V / (Na * X[5]));
                log(Q_NO_plus(T) * V / (Na * X[6]));
                log(Q_e_minus(T) * V / (Na * X[7]));
                log(Q_N2_plus(T) * V / (Na * X[8]));
                log(Q_O2_plus(T) * V / (Na * X[9]));
                log(Q_N_plus(T) * V / (Na * X[10]));
                log(Q_O_plus(T) * V / (Na * X[11]));
                log(Q_Ar(T) * V / (Na * X[12]));
                log(Q_Ar_plus(T) * V / (Na * X[13]))] 

    # Evaluate entropy 
    s_mix = uconvert(u"J/kg/K", R/mw_mix) * sum(X .* logQ(T)) + h/T # "J/kg-K"

    return s_mix
end

function equilibriumentropy(model::CO2_6s, X, h, T, P)
    # Function for evaluating equilibrium entropy as
    # a function of composition, enthalpy, 
    # temperature, and pressure (X, h, T, P)
    # Function is overloaded depending on reaction model

    # Model: # ["CO2", "CO", "N2", "N", "O2", "O"]

    # Import Constants
    @unpack R, k, Na = constants()

    # In per mole basis
    V = uconvert(u"m^3/mol", Na * k * T / P) # [m^3/mol]
    mw_mix = sum(model.weights .* X)
    R/mw_mix # "J/g-K"

    # Evaluate ln(Q/N)
    logQ(T) = [log(Q_CO2(T) * V / (Na * X[1]));
                log(Q_CO(T) * V / (Na * X[2]));
                log(Q_N2(T) * V / (Na * X[3]));
                log(Q_N(T) * V / (Na * X[4]));
                log(Q_O2(T) * V / (Na * X[5]));
                log(Q_O(T) * V / (Na * X[6]))] 

    # Evaluate entropy 
    s_mix = uconvert(u"J/kg/K", R/mw_mix) * sum(X .* logQ(T)) + h/T # "J/kg-K"

    return s_mix
end

