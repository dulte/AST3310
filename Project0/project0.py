import numpy as np


#Defines the constant mass fractons
X = 0.7
Y_3 = 1e-10
Y = 0.29

Z = 0.01
Z_Li7 = 1e-13
Z_Be7 = 1e-13

#Defines other usefull constants
m_u = 1.66053904e-27
N_A = 1


#Defines the reaction rates as lambda functions (Trenger N_A)
lambda_pp = lambda T: 1/N_A*4.01*1e-15*T**(-2/3.)*np.exp(-3.380*T**(-1/3.))*\
                    (1+0.123*T**(1/3.) + 1.09*T**(2/3.) + 0.938*T)

lambda_33 = lambda T: 1/N_A*6.04*1e10*T**(-2/3.)*np.exp(-12.276*T**(-1/3.))*\
                    (1+0.034*T**(1/3.) - 0.522*T**(2/3.) - 0.124*T + \
                    0.353*T**(4/3.) + 0.213*T**(-5/3.))

lambda_34 = lambda T: 1/N_A*5.06e6*(T/(1+4.95e-2*T))**(5/6.)*T**(-3/2.)\
                        *np.exp(-12.826*(T/(1+4.95e-2*T))**(-1/3.))

lambda_e7 = lambda T: 1/N_A*1.34e-10*T**(-1/2.)*(1-0.537*T**(1/3.) + 3.86*T**(2/3.))\
                     + 0.0027*T**(-1)*np.exp(2.515e-3*T**(-1))

lambda_17_mark = lambda T: 1/N_A*1.096e9*T**(-2/3.)*np.exp(-8.473*T**(-1/3.))*\
                        -4.830*(T/(1+0.759*T))**(5/6.)*T**(-3/2.)*\
                        np.exp(-8.472*(T/(1+0.759*T))**(-1/3.))+1.06e10*T**(-3/2)*\
                        np.exp(-30.442*T**(-1))

lambda_17 = lambda T: 1/N_A*3.11e5*T**(-2/3.)*np.exp(-10.262*T**(-1/3.))+\
                    2.53e3*T**(-3/2.)*np.exp(-7.306*T**(-1))

lambda_p14 = lambda T: 1/N_A*4.90e7*T**(-2/3.)*np.exp(-15.228*T**(-1/3.)) - 0.092*T**2)*\
                    (1+0.027*T**(1/3.) -0.778*T**(2/3.) -0.149*T + 0.261*T**(4/3.) +\
                    0.127*T**(5/3.)) + 2.37e3*T**(-3/2)*np.exp(3.011*T**(-1)) + \
                    2.19e4*np.exp(-12.53*T**(-1)) 


#Calculates the nuclear energy for a given density rho and temperature T
def calculate_nuclear_energy(rho,T):
    T_9 = T*1e-9

    #Calculates the particle number from the mass fractions

    n_p = rho*X/(1*m_u)
    n_He3 = rho*Y_3/(3*m_u)
    n_He = rho*Y/(4*m_u)

    n_Li7 = rho*Z_Li7/(7*m_u)
    n_Be7 = rho*Z_Be7/(7*m_u)


    #Calculates the
