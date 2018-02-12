import numpy as np


#Defines the constant mass fractons
X = 0.7
Y_3 = 1e-10
Y = 0.29

Z = 0.01
Z_Li7 = 1e-7
Z_Be7 = 1e-7

#Defines other usefull constants
m_u = 1.66053904e-27
N_A = 6.022e23


#Defines the reaction rates as lambda functions (Trenger N_A)
lambda_pp = lambda T: 1/N_A*4.01*1e-15*T**(-2/3.)*np.exp(-3.380*T**(-1/3.))*\
                    (1+0.123*T**(1/3.) + 1.09*T**(2/3.) + 0.938*T)

lambda_33 = lambda T: 1/N_A*6.04*1e10*T**(-2/3.)*np.exp(-12.276*T**(-1/3.))*\
                    (1+0.034*T**(1/3.) - 0.522*T**(2/3.) - 0.124*T + \
                    0.353*T**(4/3.) + 0.213*T**(-5/3.))

lambda_34 = lambda T: 1/N_A*5.61e6*(T/(1+4.95e-2*T))**(5/6.)*T**(-3/2.)\
                        *np.exp(-12.826*(T/(1+4.95e-2*T))**(-1/3.))

lambda_e7 = lambda T: 1/N_A*1.34e-10*T**(-1/2.)*(1-0.537*T**(1/3.) + 3.86*T**(2/3.)\
                     + 0.0027*T**(-1)*np.exp(2.515e-3*T**(-1)))

lambda_17_mark = lambda T: 1/N_A*(1.096e9*T**(-2/3.)*np.exp(-8.473*T**(-1/3.))\
                        -4.830*(T/(1+0.759*T))**(5/6.)*T**(-3/2.)*\
                        np.exp(-8.472*(T/(1+0.759*T))**(-1/3.))+1.06e10*T**(-3/2)*\
                        np.exp(-30.442*T**(-1)))

lambda_17 = lambda T: 1/N_A*(3.11e5*T**(-2/3.)*np.exp(-10.262*T**(-1/3.))+\
                    2.53e3*T**(-3/2.)*np.exp(-7.306*T**(-1)))





#Defines the energies for the reactions in Joule
Q_pp = (0.15 + 1.02 + 5.49)*1.602e-13

#PPI
Q_33 = 12.86*1.602e-13

#PPII
Q_34 = 1.59*1.602e-13
Q_e7 = 0.05*1.602e-13
Q_17_mark = 17.35*1.602e-13

#PPIII
Q_17 = (0.14)*1.602e-13# + (6.88  + 3.00 + 1.02)*1.602e-13


#Calculates the nuclear energy for a given density rho and temperature T
def calculate_nuclear_energy(rho,T):
    T = T*1e-9

    #Calculates the particle number from the mass fractions
    n_p = rho*X/(1*m_u)
    n_He3 = rho*Y_3/(3*m_u)
    n_He4 = rho*(Y-Y_3)/(4*m_u)

    n_Li7 = rho*Z_Li7/(7*m_u)
    n_Be7 = rho*Z_Be7/(7*m_u)


    n_e = rho/m_u*(X + 2/3.*Y_3 + 1/2*(Y-Y_3))

    #Defines the rates
    r_pp = n_p*n_p/(2*rho)*lambda_pp(T)*1e-6

    r_33 = n_He3*n_He3/(2*rho)*lambda_33(T)*1e-6
    r_34 = n_He3*n_He4/(rho)*lambda_34(T)*1e-6

    r_e7 = n_e*n_Be7/(rho)*lambda_e7(T)*1e-6
    r_17_mark = n_Li7*n_p/(rho)*lambda_17_mark(T)*1e-6

    r_17 = n_Be7*n_p/(rho)*lambda_17(T)*1e-6

    r_33_34 = (r_33 + r_34)
    if r_33_34 > r_pp:
        r_33 = (r_33/r_33_34)*r_pp
        r_34 = (r_34/r_33_34)*r_pp

    r_e7_17 = (r_e7 + r_17)
    if r_e7_17 > r_34:
        r_e7 = (r_e7/r_e7_17)*r_34
        r_17 = (r_17/r_e7_17)*r_34

    if r_17_mark > r_e7:
        r_17_mark = r_e7




    #Calculates the Energy
    e_pp = Q_pp*r_pp
    e_33 = Q_33*r_33
    e_34 = Q_34*r_34

    e_e7 = Q_e7*r_e7
    e_17_mark = Q_17_mark*r_17_mark

    e_17 = Q_17*r_17

    print(e_pp*rho)
    print(e_33*rho)
    print(e_34*rho)
    print(e_e7*rho)
    print(e_17_mark*rho)
    print(e_17*rho)

if __name__ == '__main__':
    rho = 1.62e5
    T = 1e8
    T = 1.57e7
    calculate_nuclear_energy(rho,T)
