# For simplicity we will generate separate functions for calculating each partition function as functions of temperature

# Primitive Partition Functions
Q_rot(T, Theta, sigma) = T / (sigma * Theta)
Q_trans(T, m, k, h) = ((2 * pi * m * k * T)/(h^2))^(3/2)
Q_vib(T, Theta) = 1 / (1 - exp(-Theta/T))
Q_el(T, Theta, g) = g * exp(-Theta / T)

function Q_CO2(T)
    # Constants
    @unpack h, k, R, c, hc, Na = constants()
    # Check units
    if dimension(T) == NoDims
        T = T*u"K"
    end
    # CO2 Parameters
    M_ = 44.01u"g/mol" # [g/mol]
    m_ = uconvert(u"kg",M_ / Na) # CO2 particle mass [kg]
    B_rot = uconvert(u"m^-1",0.389u"cm^-1") # [m^-1]
    omega_vib_ = uconvert.(u"m^-1",[1330 667.3 667.3 2349]u"cm^-1") # [m^-1]
    g_ele_ = [1 1]
    theta_ele_ = [0 66184.17]u"K" # [K]
    theta_rot_ = B_rot * hc / k # [K]
    theta_vib_ = (hc/k) .* omega_vib_ # [K]
    # CO2 Partition Function
    sigma = 2 # CO2 is symmetric
    Q_rot_ = Q_rot(T, theta_rot_, sigma)
    Q_tr_ = uconvert(u"m^-3",Q_trans(T, m_, k, h))
    #Q_tr_ = Q_trans(T, m_, k, h) #[m^-3]
    Q_vib_ = prod(Q_vib.(T, theta_vib_))
    Q_el_ = sum(Q_el.(T, theta_ele_, g_ele_))

    Q = Q_rot_ * Q_tr_ * Q_vib_ * Q_el_ #[m^-3]
    return Q
end

function Q_CO(T)
    # Constants
    @unpack h, k, R, c, hc, Na = constants()
    # Check units
    if dimension(T) == NoDims
        T = T*u"K"
    end
    # CO Parameters
    M_ = 28.010u"g/mol"
    m_ = uconvert(u"kg",M_ / Na) # CO particle mass [kg]
    B_rot = uconvert(u"m^-1",1.945u"cm^-1") # [m^-1]
    omega_vib_ = uconvert.(u"m^-1",2167.6u"cm^-1") # [m^-1]
    g_ele_ = [1 2]
    theta_ele_ = [0 93455.07]u"K" # [K]
    theta_rot_ = B_rot * hc / k # [K]
    theta_vib_ = (hc/k) .* omega_vib_ # [K]
    # CO Partition Function
    sigma = 1 # CO is asymmetric
    Q_rot_ = Q_rot(T, theta_rot_, sigma)
    Q_tr_ = uconvert(u"m^-3",Q_trans(T, m_, k, h)) #[m^-3]
    #Q_tr_ = Q_trans(T, m_, k, h) #[m^-3]
    Q_vib_ = prod(Q_vib.(T, theta_vib_))
    Q_el_ = sum(Q_el.(T, theta_ele_, g_ele_))

    Q = Q_rot_ * Q_tr_ * Q_vib_ * Q_el_ #[m^-3]
    return Q
end

function Q_O2(T)
    # Constants
    @unpack h, k, R, c, hc, Na = constants()
    # Check units
    if dimension(T) == NoDims
        T = T*u"K"
    end
    # O2 Parameters
    M_ = 31.99u"g/mol"
    m_ = uconvert(u"kg",M_ / Na) # O2 particle mass [kg]
    g_ele_ = [3 2]
    theta_ele_ = [0 11400]u"K" # [K]
    theta_rot_ = 2.08u"K" # [K]
    theta_vib_ = 2274.5u"K" # [K]
    # O2 Partition Function
    sigma = 2 # O2 is symmetric
    Q_rot_ = Q_rot(T, theta_rot_, sigma)
    Q_tr_ = uconvert(u"m^-3",Q_trans(T, m_, k, h)) #[m^-3]
    Q_vib_ = prod(Q_vib.(T, theta_vib_))
    Q_el_ = sum(Q_el.(T, theta_ele_, g_ele_))

    Q = Q_rot_ * Q_tr_ * Q_vib_ * Q_el_ #[m^-3]
    return Q
end

function Q_O(T)
    # Constants
    @unpack h, k, R, c, hc, Na = constants()
    # Check units
    if dimension(T) == NoDims
        T = T*u"K"
    end
    # O Parameters
    M_ = 15.999u"g/mol"
    m_ = uconvert(u"kg",M_ / Na) # O particle mass [kg]
    g_ele_ = [5 3 1]
    theta_ele_ = [0.0 227.812 326.718]u"K" # [K]
    # O Partition Function
    Q_rot_ = 1
    Q_tr_ = uconvert(u"m^-3",Q_trans(T, m_, k, h)) #[m^-3]
    Q_vib_ = 1
    Q_el_ = sum(Q_el.(T, theta_ele_, g_ele_))

    Q = Q_rot_ * Q_tr_ * Q_vib_ * Q_el_ #[m^-3]
    return Q
end

function Q_N2(T)
    # Constants
    @unpack h, k, R, c, hc, Na = constants()
    # Check units
    if dimension(T) == NoDims
        T = T*u"K"
    end
    # N2 Parameters
    M_ = 28.013u"g/mol"
    m_ = uconvert(u"kg",M_ / Na) # N2 particle mass [kg]
    g_ele_ = [1 2]
    theta_ele_ = [0 99600]u"K" # [K]
    theta_rot_ = 2.90u"K" # [K]
    theta_vib_ = 3390.0u"K" # [K]
    # N2 Partition Function
    sigma = 2 # N2 is symmetric
    Q_rot_ = Q_rot(T, theta_rot_, sigma)
    Q_tr_ = uconvert(u"m^-3",Q_trans(T, m_, k, h)) #[m^-3]
    Q_vib_ = prod(Q_vib.(T, theta_vib_))
    Q_el_ = sum(Q_el.(T, theta_ele_, g_ele_))

    Q = Q_rot_ * Q_tr_ * Q_vib_ * Q_el_ #[m^-3]
    return Q
end

function Q_N(T)
    # Constants
    @unpack h, k, R, c, hc, Na = constants()
    # Check units
    if dimension(T) == NoDims
        T = T*u"K"
    end
    # N Parameters
    M_ = 14.007u"g/mol"
    m_ = uconvert(u"kg",M_ / Na) # N particle mass [kg]
    g_ele_ = [4 6 4]
    theta_ele_ = [0.0 27691.58 27704.14]u"K" # [K]
    # N Partition Function
    Q_rot_ = 1
    Q_tr_ = uconvert(u"m^-3",Q_trans(T, m_, k, h)) #[m^-3]
    Q_vib_ = 1
    Q_el_ = sum(Q_el.(T, theta_ele_, g_ele_))

    Q = Q_rot_ * Q_tr_ * Q_vib_ * Q_el_ #[m^-3]
    return Q
end

function Q_NO(T)
    # Constants
    @unpack h, k, R, c, hc, Na = constants()
    # Check units
    if dimension(T) == NoDims
        T = T*u"K"
    end
    # NO Parameters
    M_ = 30.0061u"g/mol"
    m_ = uconvert(u"kg",M_ / Na) # NO particle mass [kg]
    B_rot = uconvert(u"m^-1",1.7042u"cm^-1") # [m^-1]
    omega_vib_ = uconvert.(u"m^-1",1903.60u"cm^-1") # [m^-1]
    epsil_ele_ = uconvert.(u"m^-1",[0 , 121.1]u"cm^-1") # [m^-1]
    g_ele_ = [2 2]
    theta_ele_ = (hc/k) .* epsil_ele_ # [K]
    theta_rot_ = B_rot * hc / k # [K]
    theta_vib_ = (hc/k) .* omega_vib_ # [K]
    # NO Partition Function
    sigma = 1 # NO is asymmetric
    Q_rot_ = Q_rot(T, theta_rot_, sigma)
    Q_tr_ = uconvert(u"m^-3",Q_trans(T, m_, k, h)) #[m^-3]
    #Q_tr_ = Q_trans(T, m_, k, h) #[m^-3]
    Q_vib_ = prod(Q_vib.(T, theta_vib_))
    Q_el_ = sum(Q_el.(T, theta_ele_, g_ele_))

    Q = Q_rot_ * Q_tr_ * Q_vib_ * Q_el_ #[m^-3]
    return Q
end

function Q_NO_plus(T)
    # Constants
    @unpack h, k, R, c, hc, Na = constants()
    # Check units
    if dimension(T) == NoDims
        T = T*u"K"
    end
    # NO+ Parameters
    M_ = 30.0056u"g/mol"
    m_ = uconvert(u"kg",M_ / Na) # NO particle mass [kg]
    B_rot = uconvert(u"m^-1",2.002u"cm^-1") # [m^-1]
    omega_vib_ = uconvert.(u"m^-1",2371.1u"cm^-1") # [m^-1]
    epsil_ele_ = uconvert.(u"m^-1",0u"cm^-1") # [m^-1]
    g_ele_ = [1]
    theta_ele_ = (hc/k) .* epsil_ele_ # [K]
    theta_rot_ = B_rot * hc / k # [K]
    theta_vib_ = (hc/k) .* omega_vib_ # [K]
    # NO+ Partition Function
    sigma = 1 # NO+ is asymmetric
    Q_rot_ = Q_rot(T, theta_rot_, sigma)
    Q_tr_ = uconvert(u"m^-3",Q_trans(T, m_, k, h)) #[m^-3]
    #Q_tr_ = Q_trans(T, m_, k, h) #[m^-3]
    Q_vib_ = prod(Q_vib.(T, theta_vib_))
    Q_el_ = sum(Q_el.(T, theta_ele_, g_ele_))

    Q = Q_rot_ * Q_tr_ * Q_vib_ * Q_el_ #[m^-3]
    return Q
end

function Q_e_minus(T)
    # Constants
    @unpack h, k, R, c, hc, Na = constants()
    # Check units
    if dimension(T) == NoDims
        T = T*u"K"
    end
    # e- Parameters
    M_ = 0.0005485799088728282u"g/mol"
    m_ = uconvert(u"kg",M_ / Na) # e- particle mass [kg]
    epsil_ele_ = uconvert.(u"m^-1",0u"cm^-1") # [m^-1]
    g_ele_ = [2]
    theta_ele_ = (hc/k) .* epsil_ele_ # [K]
    # e- Partition Function
    sigma = 1 # e- is asymmetric
    Q_rot_ = 1
    Q_tr_ = uconvert(u"m^-3",Q_trans(T, m_, k, h)) #[m^-3]
    #Q_tr_ = Q_trans(T, m_, k, h) #[m^-3]
    Q_vib_ = 1
    Q_el_ = sum(Q_el.(T, theta_ele_, g_ele_))

    Q = Q_rot_ * Q_tr_ * Q_vib_ * Q_el_ #[m^-3]
    return Q
end

function Q_N2_plus(T)
    # Constants
    @unpack h, k, R, c, hc, Na = constants()
    # Check units
    if dimension(T) == NoDims
        T = T*u"K"
    end
    # N2+ Parameters
    M_ = 28.01285u"g/mol"
    m_ = uconvert(u"kg",M_ / Na) # N2+ particle mass [kg]
    B_rot = uconvert(u"m^-1",1.9319u"cm^-1") # [m^-1]
    omega_vib_ = uconvert.(u"m^-1",[2207.0, 1903.53, 2419.84]u"cm^-1") # [m^-1]
    epsil_ele_ = uconvert.(u"m^-1",[0.0, 9016.4, 25566.0]u"cm^-1") # [m^-1]
    g_ele_ = [2, 4, 2]
    theta_ele_ = (hc/k) .* epsil_ele_ # [K]
    theta_rot_ = B_rot * hc / k # [K]
    theta_vib_ = (hc/k) .* omega_vib_ # [K]
    # N2+ Partition Function
    sigma = 2 # N2+ is symmetric
    Q_rot_ = Q_rot(T, theta_rot_, sigma)
    Q_tr_ = uconvert(u"m^-3",Q_trans(T, m_, k, h)) #[m^-3]
    #Q_tr_ = Q_trans(T, m_, k, h) #[m^-3]
    Q_vib_ = prod(Q_vib.(T, theta_vib_))
    Q_el_ = sum(Q_el.(T, theta_ele_, g_ele_))

    Q = Q_rot_ * Q_tr_ * Q_vib_ * Q_el_ #[m^-3]
    return Q
end

function Q_N_plus(T)
    # Constants
    @unpack h, k, R, c, hc, Na = constants()
    # Check units
    if dimension(T) == NoDims
        T = T*u"K"
    end
    # N+ Parameters
    M_ = 14.00615u"g/mol"
    m_ = uconvert(u"kg",M_ / Na) # N+ particle mass [kg]
    epsil_ele_ = uconvert.(u"m^-1",[0.0, 48.7, 130.8, 15316.2]u"cm^-1") # [m^-1]
    g_ele_ = [1, 3, 5, 5]
    theta_ele_ = (hc/k) .* epsil_ele_ # [K]
    # N+ Partition Function
    Q_rot_ = 1
    Q_tr_ = uconvert(u"m^-3",Q_trans(T, m_, k, h)) #[m^-3]
    #Q_tr_ = Q_trans(T, m_, k, h) #[m^-3]
    Q_vib_ = 1
    Q_el_ = sum(Q_el.(T, theta_ele_, g_ele_))

    Q = Q_rot_ * Q_tr_ * Q_vib_ * Q_el_ #[m^-3]
    return Q
end

function Q_O2_plus(T)
    # Constants
    @unpack h, k, R, c, hc, Na = constants()
    # Check units
    if dimension(T) == NoDims
        T = T*u"K"
    end
    # O2+ Parameters
    M_ = 31.99825u"g/mol"
    m_ = uconvert(u"kg",M_ / Na) # O2+ particle mass [kg]
    B_rot = uconvert(u"m^-1",1.68912u"cm^-1") # [m^-1]
    omega_vib_ = uconvert.(u"m^-1",[1905.13, 1035.534, 898.17, 119.913]u"cm^-1") # [m^-1]
    epsil_ele_ = uconvert.(u"m^-1",[0.0, 32524.0, 40070.0, 49191.0]u"cm^-1") # [m^-1]
    g_ele_ = [4, 8, 4, 4]
    theta_ele_ = (hc/k) .* epsil_ele_ # [K]
    theta_rot_ = B_rot * hc / k # [K]
    theta_vib_ = (hc/k) .* omega_vib_ # [K]
    # O2+ Partition Function
    sigma = 2 # O2+ is symmetric
    Q_rot_ = Q_rot(T, theta_rot_, sigma)
    Q_tr_ = uconvert(u"m^-3",Q_trans(T, m_, k, h)) #[m^-3]
    #Q_tr_ = Q_trans(T, m_, k, h) #[m^-3]
    Q_vib_ = prod(Q_vib.(T, theta_vib_))
    Q_el_ = sum(Q_el.(T, theta_ele_, g_ele_))

    Q = Q_rot_ * Q_tr_ * Q_vib_ * Q_el_ #[m^-3]
    return Q
end

function Q_O_plus(T)
    # Constants
    @unpack h, k, R, c, hc, Na = constants()
    # Check units
    if dimension(T) == NoDims
        T = T*u"K"
    end
    # O+ Parameters
    M_ = 15.99885u"g/mol"
    m_ = uconvert(u"kg",M_ / Na) # O+ particle mass [kg]
    epsil_ele_ = uconvert.(u"m^-1",[0.0, 26810.7, 26830.5]u"cm^-1") # [m^-1]
    g_ele_ = [4, 6, 4]
    theta_ele_ = (hc/k) .* epsil_ele_ # [K]
    # O+ Partition Function
    Q_rot_ = 1
    Q_tr_ = uconvert(u"m^-3",Q_trans(T, m_, k, h)) #[m^-3]
    #Q_tr_ = Q_trans(T, m_, k, h) #[m^-3]
    Q_vib_ = 1
    Q_el_ = sum(Q_el.(T, theta_ele_, g_ele_))

    Q = Q_rot_ * Q_tr_ * Q_vib_ * Q_el_ #[m^-3]
    return Q
end

function Q_Ar(T)
    # Constants
    @unpack h, k, R, c, hc, Na = constants()
    # Check units
    if dimension(T) == NoDims
        T = T*u"K"
    end
    # Ar Parameters
    M_ = 39.948u"g/mol"
    m_ = uconvert(u"kg",M_ / Na) # Ar particle mass [kg]
    epsil_ele_ = uconvert.(u"m^-1",0.0u"cm^-1") # [m^-1]
    g_ele_ = 1
    theta_ele_ = (hc/k) .* epsil_ele_ # [K]
    # Ar Partition Function
    Q_rot_ = 1
    Q_tr_ = uconvert(u"m^-3",Q_trans(T, m_, k, h)) #[m^-3]
    #Q_tr_ = Q_trans(T, m_, k, h) #[m^-3]
    Q_vib_ = 1
    Q_el_ = sum(Q_el.(T, theta_ele_, g_ele_))

    Q = Q_rot_ * Q_tr_ * Q_vib_ * Q_el_ #[m^-3]
    return Q
end

function Q_Ar_plus(T)
    # Constants
    @unpack h, k, R, c, hc, Na = constants()
    # Check units
    if dimension(T) == NoDims
        T = T*u"K"
    end
    # Ar+ Parameters
    M_ = 39.94745u"g/mol"
    m_ = uconvert(u"kg",M_ / Na) # O+ particle mass [kg]
    epsil_ele_ = uconvert.(u"m^-1",[0.0, 1431.41]u"cm^-1") # [m^-1]
    g_ele_ = [4, 2]
    theta_ele_ = (hc/k) .* epsil_ele_ # [K]
    # Ar+ Partition Function
    Q_rot_ = 1
    Q_tr_ = uconvert(u"m^-3",Q_trans(T, m_, k, h)) #[m^-3]
    #Q_tr_ = Q_trans(T, m_, k, h) #[m^-3]
    Q_vib_ = 1
    Q_el_ = sum(Q_el.(T, theta_ele_, g_ele_))

    Q = Q_rot_ * Q_tr_ * Q_vib_ * Q_el_ #[m^-3]
    return Q
end