function reactionsystem(model::Air5s, T, P)
    # Define each equation 
    # X is a vector containing in each location the following
    # [X_N2, X_O2, X_NO, X_N, X_O]
    # Solving in logspace
    # Y = log10(X)
    # X = 10 .^(Y)

    # Kp = [Kp_N2; Kp_O2; Kp_NO]

    Kp = equilibriumconstant(model::Air5s, T)

    # f(X, Kp, P) = [P * X[4]^2 - X[1] * Kp[1] 
    #     P * X[5]^2 - X[2] * Kp[2] 
    #     P * X[4] * X[5] - X[3] * Kp[3] 
    #     1 - sum(X)
    #     (2 * X[1] + X[3] + X[4]) - (2 * X[2] + X[3] + X[5]) * 79/21]

    f_log(Y, Kp, P) = [log10(ustrip(P)) + 2 * Y[4] - Y[1] - log10(ustrip(Kp[1])) 
        log10(ustrip(P)) + 2*Y[5] - Y[2] - log10(ustrip(Kp[2])) 
        log10(ustrip(P)) + Y[4] + Y[5] - Y[3] - log10(ustrip(Kp[3])) 
        log10(1) - log10(sum(10 .^Y))
        log10(2 * 10^Y[1] + 10^Y[3] + 10^Y[4]) - log10(2 * 10^Y[2] + 10^Y[3] + 10^Y[5]) - log10(79/21)]

    return f_log, Kp
end

function reactionsystem(model::Air7s, T, P)
    # Define each equation 
    # X is a vector containing in each location the following
    # ["N2", "O2", "NO", "N", "O", "NO+", "e-"]

    # Kp = [Kp_N2; Kp_O2; Kp_NO; Kp_NO_plus]

    Kp = equilibriumconstant(model::Air7s, T)

    # f(X, Kp, P) = [P * X[4]^2 - X[1] * Kp[1] 
    #     P * X[5]^2 - X[2] * Kp[2] 
    #     P * X[4] * X[5] - X[3] * Kp[3] 
    #     1 - sum(X)
    #     (2 * X[1] + X[3] + X[4] + X[6]) - (2 * X[2] + X[3] + X[5] + X[6]) * 79/21
    #     P * X[6] * X[7] - X[3] * Kp[4] 
    #     X[7] - (X[6])]

    f_log(Y, Kp, P) = [log10(ustrip(P)) + 2 * Y[4] - Y[1] - log10(ustrip(Kp[1])) 
        log10(ustrip(P)) + 2*Y[5] - Y[2] - log10(ustrip(Kp[2])) 
        log10(ustrip(P)) + Y[4] + Y[5] - Y[3] - log10(ustrip(Kp[3])) 
        log10(1) - log10(sum(10 .^Y))
        log10(2 * 10^Y[1] + 10^Y[3] + 10^Y[4] + 10^Y[6]) - log10(2 * 10^Y[2] + 10^Y[3] + 10^Y[5] + 10^Y[6]) - log10(79/21)
        log10(ustrip(P)) + Y[6] + Y[7] - Y[3] - log10(ustrip(Kp[4])) 
        Y[7] - log10(10^Y[6])]

    return f_log, Kp
    
end

function reactionsystem(model::Air11s, T, P)
    # Define each equation 
    # X is a vector containing in each location the following
    # ["N2", "O2", "NO", "N", "O", "NO+", "e-", "N2+", "O2+", "N+", "O+"]

    # Kp = [Kp_N2; Kp_O2; Kp_NO; Kp_NO_plus; Kp_N2_plus; Kp_O2_plus; Kp_N_plus; Kp_O_plus]

    Kp = equilibriumconstant(model::Air11s, T)

    # f(X, Kp, P) = [P * X[4]^2 - X[1] * Kp[1] 
    #     P * X[5]^2 - X[2] * Kp[2] 
    #     P * X[4] * X[5] - X[3] * Kp[3] 
    #     1 - sum(X)
    #     (2 * X[1] + X[3] + X[4] + X[6] + 2 * X[8] + X[10]) - (2 * X[2] + X[3] + X[5] + X[6] + 2 * X[9] + X[11]) * 79/21
    #     P * X[6] * X[7] - X[3] * Kp[4] 
    #     X[7] - (X[6] + X[8] + X[9] + X[10] + X[11])
    #     P * X[8] * X[7] - X[1] * Kp[5] 
    #     P * X[9] * X[7] - X[2] * Kp[6] 
    #     P * X[10] * X[7] - X[4] * Kp[7] 
    #     P * X[11] * X[7] - X[5] * Kp[8]]

    f_log(Y, Kp, P) = [log10(ustrip(P)) + 2 * Y[4] - Y[1] - log10(ustrip(Kp[1])) 
        log10(ustrip(P)) + 2*Y[5] - Y[2] - log10(ustrip(Kp[2])) 
        log10(ustrip(P)) + Y[4] + Y[5] - Y[3] - log10(ustrip(Kp[3])) 
        log10(1) - log10(sum(10 .^Y))
        log10(2 * 10^Y[1] + 10^Y[3] + 10^Y[4] + 10^Y[6] + 2 * 10^Y[8] + 10^Y[10]) - log10(2 * 10^Y[2] + 10^Y[3] + 10^Y[5] + 10^Y[6] + 2 * 10^Y[9] + 10^Y[11]) - log10(79/21)
        log10(ustrip(P)) + Y[6] + Y[7] - Y[3] - log10(ustrip(Kp[4])) 
        Y[7] - log10(10^Y[6] + 10^Y[8] + 10^Y[9] + 10^Y[10] + 10^Y[11])
        log10(ustrip(P)) + Y[8] + Y[7] - Y[1] - log10(ustrip(Kp[5])) 
        log10(ustrip(P)) + Y[9] + Y[7] - Y[2] - log10(ustrip(Kp[6])) 
        log10(ustrip(P)) + Y[10] + Y[7] - Y[4] - log10(ustrip(Kp[7])) 
        log10(ustrip(P)) + Y[11] + Y[7] - Y[5] - log10(ustrip(Kp[8]))]

    return f_log, Kp
end

function reactionsystem(model::Air13s, T, P)
    # Define each equation 
    # X is a vector containing in each location the following
    # ["N2", "O2", "NO", "N", "O", "NO+", "e-", "N2+", "O2+", "N+", "O+", "Ar", "Ar+"]

    # Kp = [Kp_N2; Kp_O2; Kp_NO; Kp_NO_plus; Kp_N2_plus; Kp_O2_plus; Kp_N_plus; Kp_O_plus; Kp_Ar_plus]

    Kp = equilibriumconstant(model::Air13s, T)

    # f(X, Kp, P) = [P * X[4]^2 - X[1] * Kp[1] 
    #     P * X[5]^2 - X[2] * Kp[2] 
    #     P * X[4] * X[5] - X[3] * Kp[3] 
    #     1 - sum(X)
    #     (2 * X[1] + X[3] + X[4] + X[6] + 2 * X[8] + X[10]) - (2 * X[2] + X[3] + X[5] + X[6] + 2 * X[9] + X[11]) * 79/21
    #     P * X[6] * X[7] - X[3] * Kp[4] 
    #     X[7] - (X[6] + X[8] + X[9] + X[10] + X[11])
    #     P * X[8] * X[7] - X[1] * Kp[5] 
    #     P * X[9] * X[7] - X[2] * Kp[6] 
    #     P * X[10] * X[7] - X[4] * Kp[7] 
    #     P * X[11] * X[7] - X[5] * Kp[8]]

    f_log(Y, Kp, P) = [log10(ustrip(P)) + 2 * Y[4] - Y[1] - log10(ustrip(Kp[1])) 
        log10(ustrip(P)) + 2*Y[5] - Y[2] - log10(ustrip(Kp[2])) 
        log10(ustrip(P)) + Y[4] + Y[5] - Y[3] - log10(ustrip(Kp[3])) 
        log10(1) - log10(sum(10 .^Y))
        log10(2 * 10^Y[1] + 10^Y[3] + 10^Y[4] + 10^Y[6] + 2 * 10^Y[8] + 10^Y[10]) - log10(2 * 10^Y[2] + 10^Y[3] + 10^Y[5] + 10^Y[6] + 2 * 10^Y[9] + 10^Y[11]) - log10(0.78090/0.20950)
        log10(ustrip(P)) + Y[6] + Y[7] - Y[3] - log10(ustrip(Kp[4])) 
        Y[7] - log10(10^Y[6] + 10^Y[8] + 10^Y[9] + 10^Y[10] + 10^Y[11] + 10^Y[13])
        log10(ustrip(P)) + Y[8] + Y[7] - Y[1] - log10(ustrip(Kp[5])) 
        log10(ustrip(P)) + Y[9] + Y[7] - Y[2] - log10(ustrip(Kp[6])) 
        log10(ustrip(P)) + Y[10] + Y[7] - Y[4] - log10(ustrip(Kp[7])) 
        log10(ustrip(P)) + Y[11] + Y[7] - Y[5] - log10(ustrip(Kp[8]))
        log10(ustrip(P)) + Y[13] + Y[7] - Y[12] - log10(ustrip(Kp[9]))
        log10(10^Y[12] + 10^Y[13]) - log10(2 * 10^Y[2] + 10^Y[3] + 10^Y[5] + 10^Y[6] + 2 * 10^Y[9] + 10^Y[11]) - log10(0.00960/0.20950)]

    return f_log, Kp
end



function reactionsystem(model::CO2_6s, T, P)
    # Define each equation 
    # X is a vector containing in each location the following
    # [X_CO2, X_CO, X_N2, X_N, X_O2, X_O]

    # Kp = [Kp_CO2; Kp_N2; Kp_O2]

    Kp = equilibriumconstant(model::CO2_6s, T)

    # f(X, Kp, P) = [sqrt(P * X[5]) * X[2] / Kp[1] - X[1]
    #     P * X[4]^2 / Kp[2] - X[3]
    #     P * X[6]^2 / Kp[3] - X[5]
    #     1 - (X[1] + X[2] + X[3] + X[4] + X[5] + X[6])
    #     (5 * 2 / 95) * (X[1] + X[2] - 95 / 5 * X[3]) - X[4]
    #     (95 / 5) * (2 * X[3] + X[4]) - 2 * X[5] - 2 * X[1] - X[2] - X[6]]

    f_log(Y, Kp, P) = [(1/2)*log10(ustrip(P)) + (1/2) * Y[5] + Y[2] - Y[1] - log10(ustrip(Kp[1])) 
        log10(ustrip(P)) + 2*Y[4] - Y[3] - log10(ustrip(Kp[2])) 
        log10(ustrip(P)) + 2*Y[6] - Y[5] - log10(ustrip(Kp[3])) 
        log10(1) - log10(10^Y[1] + 10^Y[2] + 10^Y[3] + 10^Y[4] + 10^Y[5] + 10^Y[6])
        log10(10^Y[1] + 10^Y[2]) - log10(2*10^Y[3] + 10^Y[4]) - log10(95/10)
        log10(2*10^Y[5] + 10^Y[6] + 2*10^Y[1] + 10^Y[2]) - log10(2*10^Y[3] + 10^Y[4]) - log10(95/5)]

    return f_log, Kp

end

