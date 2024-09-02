import sys

from numpy import log, exp, pi
import numpy as np
import sympy as sy
import matplotlib.pyplot as plt
import scipy.optimize
import functools
import scipy.integrate
from scipy.special import gamma as GammaFunc

def lnOmega_ls(Tstar, *, n, l, s):
    """
    Fokin correlation, with erratum fixed
    """
    A11 = [0,0,0,0,0,-1.45269,29.4682,2.42508,0.0107782,0.587725,-180.714,59.5694,0.0546646,-6.51465,374.457,-137.807,0.485352,24.5523,-336.782,81.4187,-0.385355,-20.6868,132.246,0,0.0847532,5.21812,-18.114,-7.47215]
    A22 = [0,0,0,0,0,1.13086,23.4799,3.10127,0,5.51559,-137.023,18.5848,0.0325909,-29.2925,243.741,0,0.697682,59.0192,-143.67,-123.518,-0.564238,-43.0549,0,137.282,0.126508,10.4273,15.0601,-40.8911]
    if (l,s) == (1,1):
        deltakk = 0
        Deltakk = 0.5
        A = A11
    elif (l,s) == (2,2):
        deltakk = 1
        Deltakk = 0.5
        A = A22
    else:
        raise((l, s))
    sum1 = 0
    for i in range(1, 7):
        a_i = (A[4*i] + A[4*i+1]/n + A[4*i+2]/n**2 + A[4*i+3]/n**2*log(n))
        # print(i, a_i)
        sum1 += a_i*(1/Tstar)**((i-1)*Deltakk)
    return -2/n*log(Tstar) + deltakk*log(1-2.0/(3.0*n)) + sum1

def get_Bstar_scipy(Tstar, *, n, m):
    """ Calculate B^* by numerical integration, as a check """
    C = n/(n-m)*(n/m)**(m/(n-m))
    def integrand(rstar):
        return rstar**2*(1-np.exp(-1/Tstar*C*(rstar**(-n)-rstar**(-m))))
    return 2*np.pi*scipy.integrate.quad(integrand, 0, np.inf)[0]

def get_Bstar_Sadus(Tstar, *, n, m):
    """
    # Sadus, https://doi.org/10.1063/1.5041320, erratum: missing exponent of m
    m is the attractive exponent (the 6 of 12-6 LJ)
    n is the repulsive exponent (the 12 of 12-6 LJ)
    """
    def F(y):
        the_sum = 0
        for i in range(1, 200):
            def my_factorial(k):
                return GammaFunc(k+1)
            the_sum += GammaFunc((i*m-3.0)/n)*y**i/my_factorial(i)
        return y**(3/(n-m))*(GammaFunc((n-3.0)/n) -3/n*the_sum)
    yn = (n/(n-m))**n*((n-m)/m)**m*Tstar**(-(n-m)) # y**n, Eq. 9
    y = yn**(1/n)
    return 2*pi/3*F(y)

def get_B_plus_TdBdT(Tstar, *, n, m):
    """ Get Bstar + Tstar*(dBstar/dTstar) """
    Bstar = get_Bstar_scipy(Tstar, n=n, m=m)
    h = 1e-100
    dBstardTstar = (get_Bstar_Sadus(Tstar+1j*h, n=n, m=m)/h).imag
    return (Bstar + Tstar*dBstardTstar)

def get_eta_plus(Tstar, *, n):
    """ Calculate \eta^+ """
    Omega22 = exp(lnOmega_ls(Tstar, n=n, l=2, s=2))
    etastar = 5.0*Tstar**0.5/(16.0*pi**0.5*Omega22)
    val = etastar/Tstar**(0.5)*get_B_plus_TdBdT(Tstar,n=n,m=6)**(2/3)
    return val

def get_D_plus(Tstar, *, n):
    """ Calculate D^+ """
    Omega11 = exp(lnOmega_ls(Tstar, n=n, l=1, s=1))
    rhostarDstar = 3.0*Tstar**0.5/(8.0*pi**0.5*Omega11)
    val = rhostarDstar/Tstar**(0.5)*get_B_plus_TdBdT(Tstar,n=n,m=6)**(2/3)
    return val    

def get_Boyle(*,n,m):
    """ Calculate the Boyle temperature """
    f = functools.partial(get_Bstar_Sadus, n=n, m=m)
    return scipy.optimize.newton(lambda Tstar: float(f(Tstar)), 1)

def plot_Mie_contributions():
    """ """
    Tstar = np.linspace(0.4, 20, 1000)
    fig, (ax1,ax2) = plt.subplots(1,2,figsize=(6,4))
    for n in [12, 24, 48]:
        Omega22 = np.array([exp(lnOmega_ls(T, n=n, l=2, s=2)) for T in Tstar])
        ax1.plot(Tstar, 1/Omega22,label=str(n)+'-6')
        Bterm = np.array([get_B_plus_TdBdT(T,n=n,m=6)**(2/3) for T in Tstar])
        ax2.plot(Tstar, Bterm, label=str(n)+'-6')
    ax1.set_xlabel(r'$T^*$')
    ax2.set_xlabel(r'$T^*$')
    ax1.set_ylabel(r'$1/\Omega^{(2,2)*}$')
    ax2.set_ylabel(r'$(B^*+T^*dB^*/dT^*)^{2/3}$')
    ax1.legend()
    ax2.legend()
    fig.tight_layout(pad=0.2)
    plt.savefig('Mie_contributions.pdf')
    plt.close()

if __name__ == '__main__':
    print('These values should be identical; just a quick sanity check')
    print('T^* B2^*(Sadus) B2^*(numerical integration)')
    for Tstar in [0.515501912627262,2.45633644062211,10.1106138574776]:
        print(Tstar, get_Bstar_Sadus(Tstar, n=12, m=6), get_Bstar_scipy(Tstar, n=12, m=6))
    # plot_Mie_contributions()