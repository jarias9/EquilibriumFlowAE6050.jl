function equilibriumflow(problem::NormalShock, model, u1, T1, P1)
    # Generate initial state
    h1, rho1, s1, X1 = equilibriumTP(model, T1, P1)

    # Calculate post shock state
    max_iter = 10
    ϵ = 0.1 # Initial guess ρ₁/ρ₂
    tol = 1e-7
    P2i = deepcopy(P1)
    h2i = deepcopy(h1)
    T2 = deepcopy(T1)
    s2 = deepcopy(s1)
    X2 = deepcopy(X1)
    u2 = deepcopy(u1)
    for i = 1:max_iter
    P2i = uconvert(u"Pa", P1 + rho1 * u1^2 * (1 - ϵ))
    h2i = uconvert(u"J/kg", h1 + u1^2/2 * (1 - ϵ^2))

    # Solve for updated temperature to match new h2i 
    T2, rho2, s2, X2 = equilibriumHP(model, h2i, P2i)
    u2 = rho1 * u1 / rho2 # [m/s]
    ϵ_last = ϵ
    ϵ = rho1 / rho2
    res_ϵ = abs(ϵ - ϵ_last)/abs(ϵ)
    #@printf "ϵ: %f, res_ϵ: %f, i: %f \n" ϵ res_ϵ i
    if res_ϵ < tol
        break
    end
    end

    # Export Final State 
    P2 = P2i
    h2, rho2, s2, X2 = equilibriumTP(model, T2, P2)

    return T2, h2, P2, rho2, u2, s2, X2
end

function equilibriumflow(problem::ObliqueShock, angleType::DeflectionAngle, model, angle, u1, T1, P1)
    # Note: Due to the nature of this solver
    # This function will find the wave angle 
    # corresponding to the weak solution
    # Set starting angles
    theta = angle 
    beta = angle # need starting angle for beta 

    # Generate initial state
    #h1, rho1, s1, X1 = equilibriumTP(model, T1, P1)

    # Since we do not know wave angle, we must iterate
    # until converging on beta value
    u2 = deepcopy(u1)
    u2n = deepcopy(u1)
    tol = 1e-7
    for i = 1:10
    # Calculate normal velocity value
    u1n = u1 * sind(beta)

    # Calculate post shock state
    T2, h2, P2, rho2, u2n, s2, X2 = equilibriumflow(NormalShock(), model, u1n, T1, P1)

    # Calculate β angle
    beta_last = beta
    calcangle(theta, beta, u2n, u1n) = tand(beta - theta)/tand(beta) - u2n/u1n
    angleres(beta) = calcangle(theta, beta, u2n, u1n) # uses beta as initial guess
    beta = find_zero(angleres, theta, Roots.Order0())
    #beta = find_zero(angleres, (0.001, 89.999), Bisection())
    res_θ = abs(beta - beta_last)/abs(beta)
    #@printf "beta: %f, res_ϵ: %f, i: %f \n" beta res_θ i
    if res_θ < tol
        break
    end
    end

    # Calculate Oblique velocity 
    u2 = u2n / sind(beta - theta)

    return T2, h2, P2, rho2, u2, s2, X2, beta
end

function equilibriumflow(problem::ObliqueShock, angleType::ShockAngle, model, angle, u1, T1, P1)
    # Set starting angles
    beta = angle 

    # Generate initial state
    #h1, rho1, s1, X1 = equilibriumTP(model, T1, P1)

    # Calculate normal velocity value
    u1n = u1 * sind(beta)

    # Calculate post shock state
    T2, h2, P2, rho2, u2n, s2, X2 = equilibriumflow(NormalShock(), model, u1n, T1, P1)

    # Calculate θ angle
    calcangle(theta, beta, u2n, u1n) = tand(beta - theta)/tand(beta) - u2n/u1n
    angleres(theta) = calcangle(theta, beta, u2n, u1n)
    #theta = find_zero(angleres, beta, Roots.Order0())
    theta = find_zero(angleres, (-0.00001, 90.00001))

    # Calculate Oblique velocity 
    u2 = u2n / sind(beta - theta)

    return T2, h2, P2, rho2, u2, s2, X2, theta
end