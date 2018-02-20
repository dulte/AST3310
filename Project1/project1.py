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
        self.file = filename
        self.read_opacity()
        self.G = 6.67408e-11
        self.k = 1.3806505e-23
        self.m_u = 1.660539040e-27
        self.a = 7.5657e-16
        self.mu = 1/(X+Y/4.+Z/2.)
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
        log_R = np.log10(rho*.001/(T*1e-6)**3)
        log_T = np.log10(T)

        log_kappa = self.opacity_interp2d(log_R,log_T)

        return 10**log_kappa*0.1




    def solve(self,step,L0,R0,M0,rho0,T0,adaptive=True):
        dm = step

        M = []
        R = []
        L = []
        rho = []
        T = []
        P = []

        M.append(M0)
        R.append(R0)
        L.append(L0)
        rho.append(rho0)
        T.append(T0)
        P.append(self.P(rho0,T0))

        dms = np.zeros(4)
        dV = np.zeros(4)

        V = np.zeros(4)

        p = 1e-4

        while M[-1] >= 0:
            R.append(R[-1] + dm*self.dr_dm(R[-1],rho[-1]))
            P.append(P[-1] + dm*self.dP_dm(M[-1],R[-1]))
            L.append(L[-1] + dm*self.dL_dm(rho[-1],T[-1]))
            T.append(T[-1] + dm*self.dT_dm(L[-1],R[-1],T[-1],rho[-1]))





            dV[0] = R[-2]-R[-1]
            dV[1] = P[-2]-P[-1]
            dV[2] = L[-2]-L[-1]
            dV[3] = T[-2]-T[-1]

            V[0] = R[-1]
            V[1] = P[-1]
            V[2] = L[-1]
            V[3] = T[-1]

            dV_V = np.abs(dV)/V

            if True in (dV_V>p):
                g = np.abs(dV)/np.abs(dm)
                dm = -np.min(p*V/g)
                print("Changing dm to {}".format(dm))


            rho.append(self.rho(P[-1],T[-1]))

            M.append(M[-1] + dm)

            if (rho[-1]<=0) or (T[-1]<=0):
                print("At mass {} either temperature or density is zeros".format(M[-1]))
                break


        return np.array(M),np.array(R),np.array(L),np.array(rho),\
                            np.array(T),np.array(P)



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

    M,R,L,rho,T,P = sc.solve(M_step,L0,R0,M0,rho0,T0)
    #print(R/R0)
    plt.semilogy(M/M0,rho/avg_rho)
    plt.show()
