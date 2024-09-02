"""
Stockmayer zero-density properties (Source: 10.1021/acs.jced.9b00455)
-> modified to include only necessary parts and removed ChebTools dependency
"""

# Python standard library
import os
from math import pi
import json

# Scipy stack packages
import pandas
import scipy.optimize
import scipy.interpolate 
import numpy as np
from scipy.special import gamma as GammaFunc

import matplotlib.pyplot as plt

# Pip installable
import numpy.polynomial.chebyshev as cheb

# My files
import Mie_zero_density
# import IPL 

here = os.path.dirname(__file__)

def my_factorial(k):
    return GammaFunc(k+1)

def B2star_Bartke(Tstar):
    """ 
    http://dx.doi.org/10.1103/PhysRevE.75.061503
    """
    def integrand(rstar):
        return rstar**2*(1-np.exp(-4/Tstar*(rstar**(-12)-rstar**(-6))))
    I, uI = scipy.integrate.quad(integrand, 0, np.inf)
    return 2*np.pi*I

def hj_Bartke(Tstar,*,j):
    """ 
    http://dx.doi.org/10.1103/PhysRevE.75.061503
    """
    def integrand(rstar):
        return rstar**(2-6*j)*np.exp(-4/Tstar*(rstar**(-12)-rstar**(-6)))
    I, uI = scipy.integrate.quad(integrand, 0, np.inf)
    return -2*np.pi*I

def B2star_Stockmayer_Bartke(Tstar, mustarsquared, and_increment=False):
    """ 
    http://dx.doi.org/10.1103/PhysRevE.75.061503
    """
    x = mustarsquared/Tstar
    coeffs = {1: 1/3, 2:1/25, 3: 29/11025, 4: 11/99225, 5: 13/4002075, 6: 17/243486243, 7: 523/456536705625, 8: 31/2094271554375, 9: 66197/428670161650355625, 10: 83651/63014513762602276875, 11: 21253/2222311852027773631125, 12:3660541/62502520838281133375390625}
    s = float(Mie_zero_density.get_Bstar_Sadus(Tstar,n=12,m=6))
    for j in range(1,13):
        increment = x**(2*j)*coeffs[j]*hj_Bartke(Tstar, j=j)
        s += increment
    if and_increment:
        return s, increment
    else:
        return s

# Some validation data from Table II-A in Hirschfelder; tabulated value is B^* divided by (2*pi/3)
def verify_Hirschfelder():
    df = pandas.read_csv('Stockmayer_Hirschfelder_TableII_A.csv', comment='#')
    for col in df.keys()[1::]:
        tstar = float(col.split('=')[1])
        mustarsquared = tstar/8**(-1/2)
        newcol = '(Hirsch)'+col
        df[newcol] = df.apply(lambda row: B2star_Stockmayer_Bartke(row['T^*'], mustarsquared=mustarsquared)/(2*np.pi/3), axis=1)
        devs = df[newcol]-df[col]
        err = devs[np.isfinite(df[col])].abs().sum()
        print(col, err)

def get_TdBdTB(Tstar, *, mustarsquared):
    B2star = B2star_Stockmayer_Bartke(Tstar, mustarsquared=mustarsquared)
    dT = 0.00001
    dB2stardTstar = (B2star_Stockmayer_Bartke(Tstar+dT, mustarsquared=mustarsquared)-B2star_Stockmayer_Bartke(Tstar-dT, mustarsquared=mustarsquared))/(2*dT)
    val = B2star + Tstar*dB2stardTstar
    return val

TdBdTB_expansions = {} 
for mustarsquared in np.arange(1,11):
    here = os.path.dirname(__file__)
    fname = here+'/Stockmayer_mustar2_Bterm_'+str(mustarsquared)+'.json'
    if not os.path.exists(fname) or False:
        print('Building (\mu*)^2:', mustarsquared)
        ce = cheb.Chebyshev.fit(lambda lnTstar: np.log(get_TdBdTB(np.exp(lnTstar), mustarsquared=mustarsquared)**(2/3)), [np.log(0.2), np.log(200)], deg=100)

        c = ce.coef().tolist()
        if any(~np.isfinite(c)):
            raise ValueError()
        with open(fname,'w') as fp:
            fp.write(json.dumps(dict(c=c)))

    jj = json.load(open(fname))
    ce = cheb.Chebyshev(jj['c'], domain=[np.log(0.2), np.log(200)])
    TdBdTB_expansions[mustarsquared] = ce

# https://aip.scitation.org/doi/pdf/10.1063/1.1732130
# Monchick and Mason, Table IV, Table V
dfOmega22 = pandas.read_csv(here+'/Stockmayer_MonchickMason_Omega22.csv',sep=r'\s+',engine='python', dtype=float)
dfOmega11 = pandas.read_csv(here+'/Stockmayer_MonchickMason_Omega11.csv',sep=r'\s+',engine='python', dtype=float)

def get_lambdaplus_vec(Tstarvec, *, mustarsquared):
    T_B = get_Boyle(mustarsquared=mustarsquared)
    if abs(mustarsquared-int(mustarsquared)) < 1e-10 and mustarsquared < 11:
        Tstar_dB2dT_plus_B2_23 = np.exp(TdBdTB_expansions[int(mustarsquared)].y(np.log(Tstarvec)))
    else:
        B2star = np.array([B2star_Stockmayer_Bartke(Tstar, mustarsquared=mustarsquared) for Tstar in Tstarvec])
        dT = 0.001
        dB2stardTstar = np.array([(B2star_Stockmayer_Bartke(Tstar+dT, mustarsquared=mustarsquared)-B2star_Stockmayer_Bartke(Tstar-dT, mustarsquared=mustarsquared))/(2*dT)  for Tstar in Tstarvec])
        Tstar_dB2dT_plus_B2_23 = (B2star + Tstarvec*dB2stardTstar)**(2/3)
    deltastar_S = mustarsquared/4 # I think, though this is not at all clear.
    k = 'delta*={0:0.2f}'.format(deltastar_S)
    splOmega22 = scipy.interpolate.splrep(np.log(dfOmega22['T^*']), np.log(dfOmega22[k]))
    Omega22 = np.exp(scipy.interpolate.splev(np.log(Tstarvec), splOmega22))
    eta0_over_sqrtTstar = 5/(16*pi**0.5*Omega22)
    eta0 = eta0_over_sqrtTstar*Tstar_dB2dT_plus_B2_23
    return Tstarvec/T_B, 15/4*eta0

def _Boyle(mustarsquared):
    def f(Tstar):
        val = B2star_Stockmayer_Bartke(Tstar, mustarsquared=mustarsquared)
        return val
    return scipy.optimize.newton(f, 4.0)

Boyles = {}
for mustarsquared in np.arange(1,11):
    Boyles[mustarsquared] = _Boyle(mustarsquared)

def get_Boyle(mustarsquared):
    if mustarsquared in Boyles:
        return Boyles[mustarsquared]
    else:
        return _Boyle(mustarsquared)

def calc_Stockmayer_transport(Tstarvec,mustarsquared):    
    deltastar_S = mustarsquared/4
    k = 'delta*={0:0.2f}'.format(deltastar_S)

    if abs(mustarsquared-int(mustarsquared)) < 1e-10 and mustarsquared < 11:
        Tstar_dB2dT_plus_B2_23 = np.exp(
            TdBdTB_expansions[int(mustarsquared)](np.log(Tstarvec))
        )
    else:
        B2star = np.array([B2star_Stockmayer_Bartke(Tstar, mustarsquared=mustarsquared) for Tstar in Tstarvec])
        dT = 0.001
        dB2stardTstar = np.array([(B2star_Stockmayer_Bartke(Tstar+dT, mustarsquared=mustarsquared)-B2star_Stockmayer_Bartke(Tstar-dT, mustarsquared=mustarsquared))/(2*dT)  for Tstar in Tstarvec])
        Tstar_dB2dT_plus_B2_23 = (B2star + Tstarvec*dB2stardTstar)**(2/3)

    splOmega22 = scipy.interpolate.splrep(np.log(dfOmega22['T^*']), np.log(dfOmega22[k]))
    Omega22 = np.exp(scipy.interpolate.splev(np.log(Tstarvec), splOmega22))
    splOmega11 = scipy.interpolate.splrep(np.log(dfOmega11['T^*']), np.log(dfOmega11[k]))
    Omega11 = np.exp(scipy.interpolate.splev(np.log(Tstarvec), splOmega11))

    eta0 = 5*Tstarvec**0.5/(16*pi**0.5*Omega22)
    eta0_plus = eta0/Tstarvec**0.5*Tstar_dB2dT_plus_B2_23
    rhoD0 = 3*Tstarvec**0.5/(8*pi**0.5*Omega11)
    rhoD0_plus = rhoD0/Tstarvec**0.5*Tstar_dB2dT_plus_B2_23
    lambda0 = 15/4*eta0
    lambda0_plus = 15/4*eta0_plus

    return [eta0_plus, rhoD0_plus, lambda0_plus, eta0, rhoD0, lambda0]

if __name__ == '__main__':
    xT = np.linspace(1., 20, 1000)
    mu2 = 1.0
    out = calc_Stockmayer_transport(xT, mu2)
    plt.figure()
    plt.plot(xT, out[1], label='rhoD0')
    plt.legend()
    plt.show()