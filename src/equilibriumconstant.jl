function equilibriumconstant(species::Air5s, T)
    @unpack h, k, R, c, hc, Na = constants()
    # Species
    # ["N2", "O2", "NO", "N", "O"]
    # Eqns
    #
    # N2 = 2N
    # O2 = 2O
    # NO = N + O

    # Enthalpy of Formation
    deltaH_N2 = 0.0u"J/mol" # [J/mol]
    deltaH_N = 470820.0u"J/mol" # [J/mol]
    deltaH_O2 = 0.0u"J/mol" # [J/mol]
    deltaH_O = 246790.0u"J/mol" # [J/mol]
    deltaH_NO = 89775.0u"J/mol"
    deltaH_NO_plus = 983995.0u"J/mol"
    deltaH_e_minus = 0.0u"J/mol"
    Δϵ_0_N2 = (2 * deltaH_N - deltaH_N2)/Na
    Δϵ_0_O2 = (2 * deltaH_O - deltaH_O2)/Na
    Δϵ_0_NO = (deltaH_N + deltaH_O - deltaH_NO)/Na
    Δϵ_0_NO_plus = (deltaH_NO_plus + deltaH_e_minus - deltaH_NO)/Na

    Kp_N2 = uconvert(u"Pa", (k * T) * exp(-Δϵ_0_N2 / (k*T)) * (Q_N(T)^2) / Q_N2(T)) # [Pa]
    Kp_O2 = uconvert(u"Pa", (k * T) * exp(-Δϵ_0_O2 / (k*T)) * (Q_O(T)^2) / Q_O2(T)) # [Pa]
    Kp_NO = uconvert(u"Pa", (k * T) * exp(-Δϵ_0_NO / (k*T)) * (Q_N(T) * Q_O(T)) / Q_NO(T)) # [Pa]
    #Kp_NO_plus = uconvert(u"Pa", (k * T) * exp(-Δϵ_0_NO_plus / (k*T)) * (Q_NO_plus(T) * Q_e_minus(T)) / Q_NO(T)) # [Pa]
    Kp = [Kp_N2; Kp_O2; Kp_NO]

    return Kp
end

function equilibriumconstant(species::Air7s, T)
    @unpack h, k, R, c, hc, Na = constants()
    # Species
    # ["N2", "O2", "NO", "N", "O", "NO+", "e-", "N2+", "O2+", "N+", "O+"]
    # Eqns
    #
    # N2 = 2N
    # O2 = 2O
    # NO = N + O
    # NO = NO+ + e- 

    # Enthalpy of Formation
    deltaH_N2 = 0.0u"J/mol" # [J/mol]
    deltaH_N = 470820.0u"J/mol" # [J/mol]
    deltaH_O2 = 0.0u"J/mol" # [J/mol]
    deltaH_O = 246790.0u"J/mol" # [J/mol]
    deltaH_NO = 89775.0u"J/mol"
    deltaH_NO_plus = 983995.0u"J/mol"
    deltaH_e_minus = 0.0u"J/mol"
    Δϵ_0_N2 = (2 * deltaH_N - deltaH_N2)/Na
    Δϵ_0_O2 = (2 * deltaH_O - deltaH_O2)/Na
    Δϵ_0_NO = (deltaH_N + deltaH_O - deltaH_NO)/Na
    Δϵ_0_NO_plus = (deltaH_NO_plus + deltaH_e_minus - deltaH_NO)/Na

    Kp_N2 = uconvert(u"Pa", (k * T) * exp(-Δϵ_0_N2 / (k*T)) * (Q_N(T)^2) / Q_N2(T)) # [Pa]
    Kp_O2 = uconvert(u"Pa", (k * T) * exp(-Δϵ_0_O2 / (k*T)) * (Q_O(T)^2) / Q_O2(T)) # [Pa]
    Kp_NO = uconvert(u"Pa", (k * T) * exp(-Δϵ_0_NO / (k*T)) * (Q_N(T) * Q_O(T)) / Q_NO(T)) # [Pa]
    Kp_NO_plus = uconvert(u"Pa", (k * T) * exp(-Δϵ_0_NO_plus / (k*T)) * (Q_NO_plus(T) * Q_e_minus(T)) / Q_NO(T)) # [Pa]
    Kp = [Kp_N2; Kp_O2; Kp_NO; Kp_NO_plus]

    return Kp
    
end

function equilibriumconstant(model::Air11s, T)
    @unpack h, k, R, c, hc, Na = constants()
    # Species
    # ["N2", "O2", "NO", "N", "O", "NO+", "e-", "N2+", "O2+", "N+", "O+"]
    # Eqns
    #
    # N2 = 2N
    # O2 = 2O
    # NO = N + O
    # NO = NO+ + e- 

    # Enthalpy of Formation
    deltaH_N2 = model.deltaH[1] # [J/mol]
    deltaH_O2 = model.deltaH[2] # [J/mol]
    deltaH_NO = model.deltaH[3]
    deltaH_N = model.deltaH[4] # [J/mol]
    deltaH_O = model.deltaH[5] # [J/mol]
    deltaH_NO_plus = model.deltaH[6]
    deltaH_e_minus = model.deltaH[7]
    deltaH_N2_plus = model.deltaH[8]
    deltaH_O2_plus = model.deltaH[9]
    deltaH_N_plus = model.deltaH[10]
    deltaH_O_plus = model.deltaH[11]
    Δϵ_0_N2 = (2 * deltaH_N - deltaH_N2)/Na
    Δϵ_0_O2 = (2 * deltaH_O - deltaH_O2)/Na
    Δϵ_0_NO = (deltaH_N + deltaH_O - deltaH_NO)/Na
    Δϵ_0_NO_plus = (deltaH_NO_plus + deltaH_e_minus - deltaH_NO)/Na
    Δϵ_0_N2_plus = (deltaH_N2_plus + deltaH_e_minus - deltaH_N2)/Na
    Δϵ_0_O2_plus = (deltaH_O2_plus + deltaH_e_minus - deltaH_O2)/Na
    Δϵ_0_N_plus = (deltaH_N_plus + deltaH_e_minus - deltaH_N)/Na
    Δϵ_0_O_plus = (deltaH_O_plus + deltaH_e_minus - deltaH_O)/Na

    Kp_N2 = uconvert(u"Pa", (k * T) * exp(-Δϵ_0_N2 / (k*T)) * (Q_N(T)^2) / Q_N2(T)) # [Pa]
    Kp_O2 = uconvert(u"Pa", (k * T) * exp(-Δϵ_0_O2 / (k*T)) * (Q_O(T)^2) / Q_O2(T)) # [Pa]
    Kp_NO = uconvert(u"Pa", (k * T) * exp(-Δϵ_0_NO / (k*T)) * (Q_N(T) * Q_O(T)) / Q_NO(T)) # [Pa]
    Kp_NO_plus = uconvert(u"Pa", (k * T) * exp(-Δϵ_0_NO_plus / (k*T)) * (Q_NO_plus(T) * Q_e_minus(T)) / Q_NO(T)) # [Pa]
    Kp_N2_plus = uconvert(u"Pa", (k * T) * exp(-Δϵ_0_N2_plus / (k*T)) * (Q_N2_plus(T) * Q_e_minus(T)) / Q_N2(T)) # [Pa]
    Kp_O2_plus = uconvert(u"Pa", (k * T) * exp(-Δϵ_0_O2_plus / (k*T)) * (Q_O2_plus(T) * Q_e_minus(T)) / Q_O2(T)) # [Pa]
    Kp_N_plus = uconvert(u"Pa", (k * T) * exp(-Δϵ_0_N_plus / (k*T)) * (Q_N_plus(T) * Q_e_minus(T)) / Q_N(T)) # [Pa]
    Kp_O_plus = uconvert(u"Pa", (k * T) * exp(-Δϵ_0_O_plus / (k*T)) * (Q_O_plus(T) * Q_e_minus(T)) / Q_O(T)) # [Pa]
    Kp = [Kp_N2; Kp_O2; Kp_NO; Kp_NO_plus; Kp_N2_plus; Kp_O2_plus; Kp_N_plus; Kp_O_plus]

    if any(ustrip.(Kp).==0)
        Kp[ustrip.(Kp).==0] .= floatmin()*u"Pa"
    end

    return Kp
end

function equilibriumconstant(model::Air13s, T)
    @unpack h, k, R, c, hc, Na = constants()
    # Species
    # ["N2", "O2", "NO", "N", "O", "NO+", "e-", "N2+", "O2+", "N+", "O+", "Ar", "Ar+"]
    # Eqns
    #
    # N2 = 2N
    # O2 = 2O
    # NO = N + O
    # NO = NO+ + e- 

    # Enthalpy of Formation
    deltaH_N2 = model.deltaH[1] # [J/mol]
    deltaH_O2 = model.deltaH[2] # [J/mol]
    deltaH_NO = model.deltaH[3]
    deltaH_N = model.deltaH[4] # [J/mol]
    deltaH_O = model.deltaH[5] # [J/mol]
    deltaH_NO_plus = model.deltaH[6]
    deltaH_e_minus = model.deltaH[7]
    deltaH_N2_plus = model.deltaH[8]
    deltaH_O2_plus = model.deltaH[9]
    deltaH_N_plus = model.deltaH[10]
    deltaH_O_plus = model.deltaH[11]
    deltaH_Ar = model.deltaH[12]
    deltaH_Ar_plus = model.deltaH[13]
    Δϵ_0_N2 = (2 * deltaH_N - deltaH_N2)/Na
    Δϵ_0_O2 = (2 * deltaH_O - deltaH_O2)/Na
    Δϵ_0_NO = (deltaH_N + deltaH_O - deltaH_NO)/Na
    Δϵ_0_NO_plus = (deltaH_NO_plus + deltaH_e_minus - deltaH_NO)/Na
    Δϵ_0_N2_plus = (deltaH_N2_plus + deltaH_e_minus - deltaH_N2)/Na
    Δϵ_0_O2_plus = (deltaH_O2_plus + deltaH_e_minus - deltaH_O2)/Na
    Δϵ_0_N_plus = (deltaH_N_plus + deltaH_e_minus - deltaH_N)/Na
    Δϵ_0_O_plus = (deltaH_O_plus + deltaH_e_minus - deltaH_O)/Na
    Δϵ_0_Ar_plus = (deltaH_Ar_plus + deltaH_e_minus - deltaH_Ar)/Na

    Kp_N2 = uconvert(u"Pa", (k * T) * exp(-Δϵ_0_N2 / (k*T)) * (Q_N(T)^2) / Q_N2(T)) # [Pa]
    Kp_O2 = uconvert(u"Pa", (k * T) * exp(-Δϵ_0_O2 / (k*T)) * (Q_O(T)^2) / Q_O2(T)) # [Pa]
    Kp_NO = uconvert(u"Pa", (k * T) * exp(-Δϵ_0_NO / (k*T)) * (Q_N(T) * Q_O(T)) / Q_NO(T)) # [Pa]
    Kp_NO_plus = uconvert(u"Pa", (k * T) * exp(-Δϵ_0_NO_plus / (k*T)) * (Q_NO_plus(T) * Q_e_minus(T)) / Q_NO(T)) # [Pa]
    Kp_N2_plus = uconvert(u"Pa", (k * T) * exp(-Δϵ_0_N2_plus / (k*T)) * (Q_N2_plus(T) * Q_e_minus(T)) / Q_N2(T)) # [Pa]
    Kp_O2_plus = uconvert(u"Pa", (k * T) * exp(-Δϵ_0_O2_plus / (k*T)) * (Q_O2_plus(T) * Q_e_minus(T)) / Q_O2(T)) # [Pa]
    Kp_N_plus = uconvert(u"Pa", (k * T) * exp(-Δϵ_0_N_plus / (k*T)) * (Q_N_plus(T) * Q_e_minus(T)) / Q_N(T)) # [Pa]
    Kp_O_plus = uconvert(u"Pa", (k * T) * exp(-Δϵ_0_O_plus / (k*T)) * (Q_O_plus(T) * Q_e_minus(T)) / Q_O(T)) # [Pa]
    Kp_Ar_plus = uconvert(u"Pa", (k * T) * exp(-Δϵ_0_Ar_plus / (k*T)) * (Q_Ar_plus(T) * Q_e_minus(T)) / Q_Ar(T)) # [Pa]
    Kp = [Kp_N2; Kp_O2; Kp_NO; Kp_NO_plus; Kp_N2_plus; Kp_O2_plus; Kp_N_plus; Kp_O_plus; Kp_Ar_plus]

    if any(ustrip.(Kp).==0)
        Kp[ustrip.(Kp).==0] .= floatmin()*u"Pa"
    end

    return Kp
end

function equilibriumconstant(species::CO2_6s, T)
    @unpack h, k, R, c, hc, Na = constants()

    # Enthalpy of Formation
    deltaH_CO = -113805.0u"J/mol" # [J/mol]
    deltaH_CO2 = -393151u"J/mol" # [J/mol]
    deltaH_N2 = 0.0u"J/mol" # [J/mol]
    deltaH_N = 470820.0u"J/mol" # [J/mol]
    deltaH_O2 = 0.0u"J/mol" # [J/mol]
    deltaH_O = 246790.0u"J/mol" # [J/mol]
    Δϵ_0_CO2 = (deltaH_CO + 1/2 * deltaH_O2 - deltaH_CO2)/Na
    Δϵ_0_N2 = (2 * deltaH_N - deltaH_N2)/Na
    Δϵ_0_O2 = (2 * deltaH_O - deltaH_O2)/Na

    Kp_CO2 = uconvert(u"Pa^(1/2)", sqrt(k * T) * exp(-Δϵ_0_CO2 / (k*T)) * (sqrt(Q_O2(T)) * Q_CO(T)) / Q_CO2(T)) # [Pa^(1/2)]
    Kp_N2 = uconvert(u"Pa", (k * T) * exp(-Δϵ_0_N2 / (k*T)) * (Q_N(T)^2) / Q_N2(T)) # [Pa]
    Kp_O2 = uconvert(u"Pa", (k * T) * exp(-Δϵ_0_O2 / (k*T)) * (Q_O(T)^2) / Q_O2(T)) # [Pa]
    Kp = [Kp_CO2; Kp_N2; Kp_O2]

    return Kp
end