import numpy as np
import scipy.constants as const
import scipy.optimize
import copy
try:
    from . import kineticL
except:
    import kineticL
pi = np.pi


def Zl(Z0, Z2, beta, l):
    return Z0*(Z2+1j*Z0*np.tan(beta*l))/(Z0+1j*Z2*np.tan(beta*l))


def Vl(Z0, Z2, beta, l):
    return np.exp(-1j*beta*l)*(1+(Z2-Z0)/(Z2+Z0))/(1+(Zl(Z0, Z2, beta, l)-Z0)/(Zl(Z0, Z2, beta, l)+Z0))


class TLResonator:
    def __init__(self, pars):
        self.setpars(pars)

    def setpars(self, pars):
        self.epseff = pars['epseff']
        self.tand = pars['tand']
        self.L = pars['L']
        self.Lc = pars['Lc']
        self.Cin = pars['Cin']
        self.Z0 = pars['Z0']
        self.Z1 = pars['Z1']

    def w0(self):
        return const.c/np.sqrt(self.epseff)/2./self.L*2.*pi

    def lumpedL(self):
        return self.Z1*pi/2/self.w0()

    def f0(self):
        return self.w0()/2/pi

    def lumpedC(self):
        return 2./pi/self.w0()/self.Z1

    def lumpedR(self):
        k = self.w0()/const.c*np.sqrt(self.epseff)
        alpha = k*self.tand/2.
        return self.Z1*alpha*self.L

    def w1(self):
        C = self.lumpedC()
        L = self.lumpedL()
        Cs = self.Cin
        return np.sqrt((C+Cs)/L/C/Cs)

    def f1(self):
        return self.w1()/2/pi

    def S11(self, f):
        epseff = self.epseff
        tand = self.tand
        L = self.L
        Lc = self.Lc
        Cin = self.Cin
        Z0 = self.Z0
        Z1 = self.Z1

        Zc = 1j*2*pi*f*Lc
        omega = 2*pi*f
        k = omega/const.c*np.sqrt(epseff)
        alpha = k*tand/2
        beta = k-1j*alpha

        Zin = Zl(Z1, Zc, beta, L)
        Zinp = 1/(1j*omega*Cin + 1/Zin)
        result = (Zinp - Z0)/(Zinp + Z0)
        return result

    def absS11(self, f):
        return abs(self.S11(f))


def f0vsLc(Lc, myline):
    fstart = const.c/np.sqrt(myline.epseff)*1./2./myline.L*2./3.
    line = copy.copy(myline)
    line.Lc = Lc
    res = scipy.optimize.fmin(line.absS11, fstart, disp=0)
    return res[0]


def df0vsLc(Lc, myline, sign=1.0, factor=2.):
    f1 = f0vsLc(Lc, myline)
    f2 = f0vsLc(factor*Lc, myline)
    return (f1-f2)*sign


def getlinepars(geo):
    s = geo['s']
    w = geo['w']
    t = geo['t']
    london = geo['london']
    epsr = geo['epsr']

    c = const.c
    LL = kineticL.Lg(s, w)
    LLk = kineticL.Lk(s, w, t, london)
    CC = kineticL.Cg(s, w, epsr)
    Z0 = np.sqrt(LL/CC)
    Z0l = np.sqrt((LL+LLk)/CC)
    vph = 1/np.sqrt(CC*(LL+LLk))
    vph0 = 1/np.sqrt(CC*LL)
    epseff = (c/vph)**2
    epseff0 = (c/vph0)**2

    pars = {}
    pars['epseff0'] = epseff0
    pars['epseff_kinetic'] = epseff
    pars['Z0'] = Z0
    pars['Z0_kinetic'] = Z0l
    pars['Lg'] = LL
    pars['Lk'] = LLk
    pars['Cg'] = CC
    pars['vph'] = vph
    pars['vph0'] = vph0
    return pars
