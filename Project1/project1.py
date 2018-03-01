import numpy as np
import matplotlib.pyplot as plt

from project0 import calculate_nuclear_energy
from scipy.interpolate import interp2d

#Defines the constant mass fractons
X = 0.7
Y_3 = 1e-10
Y = 0.29

Z = 0.01
Z_Li7 = 1e-7
Z_Be7 = 1e-7


class stellar_core:
    def __init__(self,filename="opacity.txt"):

        #Reads and interpolates the opacity file
        self.file = filename
        self.read_opacity()

        #Sets the constants needed for the simulation
        self.G = 6.67408e-11
        self.k = 1.3806505e-23
        self.m_u = 1.660539040e-27
        self.a = 7.5657e-16
        self.mu = 1/2.*1/(X+Y/4.+Z/2.)
        self.sigma = 5.670367e-8


    def read_opacity(self):
            self.opacity = []
            self.log_Rs = []
            self.log_Ts = []
            with open(self.file) as f:
                self.log_Rs = np.array(f.readline().split()[1:]).astype(np.float)
                f.readline() #Skips the empty second row
                for line in f.readlines():
                    words = line.split()
                    self.log_Ts.append(words[0])
                    self.opacity.append(words[1:])

            self.opacity = np.array(self.opacity)
            self.log_Ts = np.array(self.log_Ts).astype(np.float)

            self.opacity_interp2d = interp2d(self.log_Rs,self.log_Ts,self.opacity)



    def get_opacity(self,T,rho):
        #Makes R from rho and T. Takes the log_10 of R and T
        log_R = np.log10(rho*.001/(T*1e-6)**3)
        log_T = np.log10(T)

        #Gets log_10(kappa) in cgs units
        log_kappa = self.opacity_interp2d(log_R,log_T)

        #Returns kappa in SI
        return 10**log_kappa*0.1



    def solve(self,step,L0,R0,M0,rho0,T0,adaptive=True):
        dm = step

        #Initiaze empty lists
        M = []; R = []; L = []; rho = []; T = []; P = []

        #Appends the initial values to the empty lists
        M.append(M0)
        R.append(R0)
        L.append(L0)
        rho.append(rho0)
        T.append(T0)
        P.append(self.P(rho0,T0))


        #Makes arrays for the use in the adaptive steps
        dms = np.zeros(4)
        dV = np.zeros(4)
        V = np.zeros(4)

        #Tolerence for the change dV/V, used for adaptive steps
        p = 1e-4

        while M[-1] >= 0.05*M0:
            try:
                #Integrates using first order Euler
                R.append(R[-1] + dm*self.dr_dm(R[-1],rho[-1]))
                P.append(P[-1] + dm*self.dP_dm(M[-1],R[-2]))
                L.append(L[-1] + dm*self.dL_dm(rho[-1],T[-1]))
                T.append(T[-1] + dm*self.dT_dm(L[-2],R[-2],T[-1],rho[-1]))

                rho.append(self.rho(P[-1],T[-1]))

                M.append(M[-1] + dm)


                #Handles adaptive steps
                if adaptive:

                    dV[0] = R[-2]-R[-1]
                    dV[1] = P[-2]-P[-1]
                    dV[2] = L[-2]-L[-1]
                    dV[3] = T[-2]-T[-1]

                    V[0] = R[-1]
                    V[1] = P[-1]
                    V[2] = L[-1]
                    V[3] = T[-1]

                    dV_V = np.abs(dV)/V

                    if (True in (dV_V>p)):
                        g = np.abs(dV)/np.abs(dm)
                        dm = -np.min(p*V/g)
                        #print("Changing dm to {}".format(dm))

                # if abs(dm) < 1e20:
                #     print("To small dm")
                #     break
                #

                # if R[-1] < 0.05*R0:
                #     break

                #Checks if certain values are below zero to stop non-phyical behavior
                if (R[-1]<=0.0*R0) or (L[-1] <= 0.0*L0):
                    print("At mass {} either temperature, luminocity or density is zeros".format(M[-1]))
                    print("Rho: ",rho[-1])
                    print("T: ",T[-1])
                    print("L: ",L[-1])
                    break

            except KeyboardInterrupt:
                leng = min([len(M),len(R),len(L),len(rho),len(T)])
                plt.plot(np.array(M[:leng])/M0,np.array(R[:leng])/R0)
                plt.show()
                plt.plot(np.array(M[:leng])/M0,np.array(L[:leng])/L0)
                plt.show()
                plt.plot(np.array(M[:leng])/M0,np.array(T[:leng])/T0)
                plt.show()
                plt.plot(np.array(M[:leng])/M0,np.array(rho[:leng])/rho0)
                plt.show()
                exit()
        return np.array(M),np.array(R),np.array(L),np.array(rho),\
                            np.array(T),np.array(P)



    #Functions for the values and derivatives of the needed quantities
    def dr_dm(self,r,rho):
        return 1./(4*np.pi*r**2*rho)

    def dP_dm(self,m,r):
        return -self.G*m/(4*np.pi*r**4)

    def dL_dm(self,rho,T):
        eps = calculate_nuclear_energy(rho,T)
        return eps

    def dT_dm(self,L,r,T,rho):
        kappa = self.get_opacity(T,rho)
        return -3*kappa*L/(256*np.pi**2*self.sigma*r**4*T**3)

    def P(self,rho,T):
        return rho/(self.mu*self.m_u)*self.k*T + self.a/3.*T**4

    def rho(self,P,T):
        return (P - self.a/3.*T**4)*(self.mu*self.m_u)/(self.k*T)

    #Method for optimizing initial values
    def optimize_initial_condition(self,step,L0,R0,M0,rho0,T0):
        # L = L0
        # R = R0
        # M = M0
        # rho = rho0
        # T = T0
        init_values = np.array([R0,L0,rho0,T0])
        prev_init_values = np.copy(init_values)

        change = np.zeros(4)
        prev_change = np.zeros(4)

        last_values = np.zeros(4)
        prev_last_values = np.zeros(4)


        search = True
        number_of_stagnations = 0
        number_of_searches = 0

        eps = 1e-6

        while search:

            R0 = init_values[0]
            L0 = init_values[1]
            rho0 = init_values[2]
            T0 = init_values[3]

            print("Search nr. {} with {}, {}, {}, {} and {}".format(number_of_searches,M0,R0,L0,rho0,T0))
            number_of_searches += 1



            M,R,L,rho,T,P = self.solve(step,L0,R0,M0,rho0,T0,adaptive=False)


            last_values[0] = R[-1]
            last_values[1] = L[-1]
            last_values[2] = rho[-1]
            last_values[3] = T[-1]

            change = (last_values - prev_last_values)/(last_values)
            r = np.random.random()
            if np.sum(prev_last_values)/(np.sum(last_values) + eps) >= r:#np.sum(change) > np.sum(prev_change):
                prev_init_values  = init_values
                prev_last_values = last_values
                prev_change = change

                number_of_stagnations = 0
            else:
                number_of_stagnations += 1

            # if number_of_stagnations >= 10:
            #     print("Minumum reached with: ", init_values)
            #     search = False


            init_values += (np.random.normal(size=4)*0.05)*init_values #0.5*change*init_values



        M0 = init_values[0]
        R0 = init_values[1]
        L0 = init_values[2]
        rho0 = init_values[3]
        T0 = init_values[4]
        M,R,L,rho,T,P = self.solve(step,L0,R0,M0,rho0,T0,adaptive=False)
        plt.plot(M/M0,R/R0)
        plt.show()
        plt.plot(M/M0,L/L0)
        plt.show()
        plt.plot(M/M0,T/T0)
        plt.show()
        plt.plot(M/M0,rho/rho0)
        plt.show()











if __name__ == '__main__':
    sc = stellar_core()
    M_sun = 1.989e30
    L_sun = 3.828e26
    R_sun = 695700e3
    avg_rho = 1.408e3

    M0 = 0.8*M_sun
    L0 = L_sun
    R0 = 0.72*R_sun

    rho0 = 5.1*avg_rho
    T0 = 5.7e6
    M_step = -1e-4*M_sun

    #M0 *=0.5
    L0 *=.9
    R0 *=0.5
    rho0 *=0.5
    T0 *=0.8

    #sc.optimize_initial_condition(M_step,L0,R0,M0,rho0,T0)
    M,R,L,rho,T,P = sc.solve(M_step,L0,R0,M0,rho0,T0,adaptive=False)
    #print(R/R0)
    #plt.semilogy(M/M0,rho/avg_rho)
    print(sc.rho(sc.P(rho0,T0),T0),rho0) #Check for Rho(P(rho0,T0),T0) = rho0
    print(M[-1]/M0)
    print(R[-1]/R0)
    print(L[-1]/L0)

    plt.plot(M/M0,R/R0)
    plt.title("R")
    plt.show()
    plt.plot(M/M0,L/L0)
    plt.title("L")
    plt.show()
    exit()
    plt.plot(M/M0,rho/rho0)
    plt.title("Rho")
    plt.show()
    plt.plot(M/M0,T/T0)
    plt.title("T")
    plt.show()
    plt.plot(M/M0,P/sc.P(rho0,T0))
    plt.title("P")
    plt.show()
