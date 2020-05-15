# kinetic L from Clem paper
# https://aip.scitation.org/doi/10.1063/1.4773070
from scipy import interpolate
import scipy.constants as const
from scipy.special import ellipkinc
from scipy.special import ellipk
from scipy.integrate import dblquad
import numpy as np

try:
    pdata = np.loadtxt('pparameter.dat',
                       usecols=list(range(0, 3)), delimiter=',')
except:
    pdata = np.loadtxt('src/pparameter.dat',
                       usecols=list(range(0, 3)), delimiter=',')

zs = pdata[:, 0]
onemps = pdata[:, 1]
ps = pdata[:, 2]
funcp = interpolate.interp1d(zs, ps)
func1mp = interpolate.interp1d(zs, onemps)


def parp(z):
    if z < 0.01:
        return 1-0.67*z
    elif z > 100.:
        return 0.63/np.sqrt(z)
    elif z > 1.:
        return funcp(z)
    elif z <= 1.:
        return 1-func1mp(z)


def Lk(s, w, t, london):
    Lambda = 2*london**2./t
    a = s/2.
    b = s/2.+w
    z = Lambda/a
    p = parp(z)
    k = a/b
    return const.mu_0*Lambda/4/a*gkp(k, p)


def gkp(k, p):
    atanh = np.arctanh
    return ((k+p**2)*atanh(p)-(1+k*p**2)*atanh(k*p))/(p*(1-k**2)*(atanh(p))**2)


def Kz(x, s, w, t, london):
    Lambda = 2*london**2./t
    a = s/2.
    b = s/2.+w
    z = Lambda/a
    p = parp(z)
    k = a/b

    if abs(x) > a and abs(x) < b:
        return 0

    C = b*p/(4*ellipkinc(np.arcsin(p), k**2))
    if abs(x) < a:
        return 2*C/np.sqrt((a**2.-p**2*x**2)*(b**2-p**2*x**2))
    elif abs(x) > b:
        return -2*C/np.sqrt((x**2.-p**2*a**2)*(x**2-p**2*b**2))


def Lg(s, w):
    a = s/2.
    b = s/2.+w
    k = a/b
    kp = np.sqrt(1-k**2)
    return const.mu_0/4.*ellipk(kp**2)/ellipk(k**2)


def Cg(s, w, epsr):
    a = s/2.
    b = s/2.+w
    k = a/b
    kp = np.sqrt(1-k**2)
    return 4*const.epsilon_0*(epsr+1)/2.*ellipk(k**2)/ellipk(kp**2)


'''
def y1(x):
    return  -np.inf

def y2(x):
    return  np.inf

def intLg(x,y,s,w,t,london):
    if x == y:
        result = 0.
    else:
        result = np.log(abs(x-y))*Kz(x,s,w,t,london)*Kz(y,s,w,t,london)
    return result

def Lg(s,w,t,london):
    return dblquad(intLg,-np.inf,np.inf,y1,y2,args=(s,w,t,london))
'''
