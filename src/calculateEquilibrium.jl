function calculateEquilibrium(model, T, P)
    # Solve species model for equilibrium composition

    # Generate initial composition guess in logspace
    X_start = ones(Float64, length(model.species))/length(model.species)
    Y = log10.(X_start)

    # Generate reaction system of equations to solve 
    # and associated equilibrium constants
    f, Kp = reactionsystem(model, T, P)

    # Execute NL solver using LevenbergMarquardt method
    f_sys(Y) = f(Y, ustrip.(Kp), ustrip(P))
    Y_new = optimize(f_sys, Y, LevenbergMarquardt(LeastSquaresOptim.LSMR()), autodiff = :forward)
    Y_min = Y_new.minimizer

    # Rescale to linspace
    X_min = 10 .^Y_min

    return X_min
end

# function calculateEquilibrium(model, T, P)
#     #X_start = zeros(Float64, length(model.species)) * eps()
#     #X_start[1:2] = [0.78, 0.22]
#     X = X_start

#     Y = log10.(X_start)

#     f, Kp = reactionsystem(model, T, P)

#     f_sys(Y) = f(Y, ustrip.(Kp), ustrip(P))
#     ∂f_sys(Y) = ForwardDiff.jacobian(Y -> f(Y, ustrip.(Kp), ustrip(P)), Y)

#     Y_new = nlsolve(f_sys, ∂f_sys, Y)
#     Y_min = Y_new.zero
#     X_min = 10 .^Y_min
#     #X = optimize(f_optim, X_start, Dogleg(LeastSquaresOptim.QR()))

#     return X_min
# end

# function generateJacobian(model, f, X, T, P, ϵ)
#     m = length(X)
#     du = zeros(Float64, m, m)
#     for j = 1:m
#         eps = zeros(Float64, m)
#         eps[j] = ϵ
#         du[:,j] = (f(model, X + eps, T, P) .- f(model, X, T, P))./ϵ
#     end
#     return du
# end
