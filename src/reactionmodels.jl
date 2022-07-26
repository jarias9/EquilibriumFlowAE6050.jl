@with_kw struct Air5s
    species::Vector{String}=["N2","O2","NO","N","O"]
    weights=[28.013, 31.99, 30.0061, 14.0067, 15.9994]u"g/mol" # [g/mol]
    deltaH=[0.0, 0.0, 89775.0, 470820.0, 246790.0]u"J/mol"
    deltaH298=[0.0, 0.0, 90291.0, 472683.0, 249173.0]u"J/mol"
end

@with_kw struct Air7s
    species::Vector{String}=["N2","O2","NO","N","O", "NO+", "e-"]
    weights=[28.013, 31.99, 30.0061, 14.0067, 15.9994, 30.0056, 0.0005485799088728282]u"g/mol" # [g/mol]
    deltaH=[0.0, 0.0, 89775.0, 470820.0, 246790.0, 983995.0, 0.0]u"J/mol"
    deltaH298=[0.0, 0.0, 90291.0, 472683.0, 249173.0, 990185.0, 0.0]u"J/mol"
end

@with_kw struct Air11s
    species::Vector{String}=["N2", "O2", "NO", "N", "O", "NO+", "e-", "N2+", "O2+", "N+", "O+"]
    weights=[28.013, 31.99, 30.0061, 14.0067, 15.9994, 30.0056, 0.0005485799088728282, 28.01285, 31.99825, 14.00615, 15.99885]u"g/mol" # [g/mol]
    deltaH=[0.0, 0.0, 89775.0, 470820.0, 246790.0, 983995.0, 0.0, 1503303.0, 1164700.0, 1873152.0, 1560733.0]u"J/mol"
    deltaH298=[0.0, 0.0, 90291.0, 472683.0, 249173.0, 990185.0, 0.0, 1509499.0, 1170888.0, 1882130.0, 1568786.0]u"J/mol"
end

@with_kw struct Air13s
    species::Vector{String}=["N2", "O2", "NO", "N", "O", "NO+", "e-", "N2+", "O2+", "N+", "O+", "Ar", "Ar+"]
    weights=[28.013, 31.99, 30.0061, 14.0067, 15.9994, 30.0056, 0.0005485799088728282, 28.01285, 31.99825, 14.00615, 15.99885, 39.948, 39.94745]u"g/mol" # [g/mol]
    deltaH=[0.0, 0.0, 89775.0, 470820.0, 246790.0, 983995.0, 0.0, 1503303.0, 1164700.0, 1873152.0, 1560733.0, 0.0, 1520573.0]u"J/mol"
    deltaH298=[0.0, 0.0, 90291.0, 472683.0, 249173.0, 990185.0, 0.0, 1509499.0, 1170888.0, 1882130.0, 1568786.0, 0.0, 1526778.0]u"J/mol"
end

@with_kw struct CO2_6s
    species::Vector{String}=["CO2", "CO", "N2", "N", "O2", "O"]
    weights=[44.01, 28.010, 28.013, 14.007, 31.99, 15.999]u"g/mol" # [g/mol]
    deltaH=[-393151.0, -113805.0, 0.0, 470820.0, 0.0, 246790.0]u"J/mol"
    deltaH298=[-393522.0, -110527.0, 0.0, 472683.0, 0.0, 249173.0]u"J/mol"
end
