#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt

Ra_max = 1.5
ABS_EPS = 1e-15

LEN = 500

_Ra = 0
_z_avg = 0
_dz_avg = 0

_dA = 0

_Rq = 0
_dRq = 0

_Rsk = 0

def smooth_abs(x):
    return np.sqrt(x**2 + ABS_EPS)

def smooth_abs_deriv(x):
    return np.divide(x, np.sqrt(np.power(x, 2) + ABS_EPS))

def A(x, z):
    A = 0
    for i in range(1, LEN):
        A += ((z[i]-z[i-1])/(x[i]-x[i-1]))*((x[i]**2)/2-(x[i-1]**2)/2) + z[i-1]*x[i] - z[i]*x[i-1]

    return A

def dA_dz(x, z):
    dA = np.zeros_like(x)
    for i in range(1, LEN):
        dA[i-1] += ((-1)/(x[i]-x[i-1]))*((x[i]**2)/2-(x[i-1]**2)/2) + x[i]
        dA[i-1] += (( 1)/(x[i]-x[i-1]))*((x[i]**2)/2-(x[i-1]**2)/2) - x[i-1]

    return dA

def z_avg(x, z):
    global _z_avg

    _A = A(x, z)
    _z_avg = _A/(x[-1] - x[0])

def dz_avg_dz(x, z):
    global _dz_avg

    _dz_avg = _dA/(x[-1] - x[0])


def Ra(x, z):
    global _Ra

    z_avg(x, z)

    _Ra = np.sum(smooth_abs(z-_z_avg))/LEN
    return (_Ra - Ra_max)**2

def dRa_dz(x, z):

    dz_avg_dz(x, z)

    return 2*(_Ra - Ra_max)*np.multiply(smooth_abs_deriv(z - _z_avg), 1 - _dz_avg)

def Rq(z):
    return np.sqrt(np.sum(np.power(z - _z_avg, 2))/LEN)

def dRq_dz(z):
    return np.multiply((z - _z_avg), 1 - _dz_avg)/(LEN*_Rq)

def Rsk(z):
    return np.sum(np.power(z - _z_avg, 3))/(LEN*_Rq**3)

def dRsk_dz(z):
    return 3*(np.multiply(np.power(z - _z_avg, 2), 1 - _dz_avg)/(LEN*_Rq**3) - (_Rsk/_Rq)*_dRq)
   
delta_max = 1e10
delta_min = 1e-30
delta = 1e-1

z_min = -4
z_max = 4
rng = np.random.Generator(np.random.PCG64())
z = rng.uniform(z_min, z_max, LEN)
x = np.linspace(0, LEN-1, LEN)
zold1 = z
zold2 = z
delta_e = np.ones_like(x)*delta

# This is constant
_dA = dA_dz(x, z)

lRa = Ra(x, z)
lRa_old1 = lRa
lRa_old2 = lRa

_Rq = Rq(z)
_dRq = dRq_dz(z)

_Rsk = Rsk(z)
Rsk_old1 = _Rsk
Rsk_old2 = _Rsk

ETA0 = 1
ETA1 = 5e2
L0 = lRa
L1 = _Rsk**2

plt.ion()

fig, ax = plt.subplots()
ax.set_xlabel("x")
ax.set_ylabel("z")
line, = plt.plot(x, z)

plt.show()

while abs(lRa) > 1e-6 or abs(_Rsk) > 1e-6:
    S = ETA0*dRa_dz(x, z) + ETA1*L1*2*_Rsk*dRsk_dz(z)

    Ra_e = (lRa-lRa_old1)*(lRa_old1-lRa_old2)
    Rsk_e = (_Rsk-Rsk_old1)*(Rsk_old1-Rsk_old2)

    for i in range(LEN):
        d_e = (z[i] - zold1[i])*(zold1[i] - zold2[i])
        if d_e < 0 or Ra_e < 0 or Rsk_e < 0:
            delta_e[i] = max(0.9*delta_e[i], delta_min)
        elif d_e > 0:
            delta_e[i] = min(2*delta_e[i], delta_max)

        z_inf = z[i] - delta_e[i]
        z_sup = z[i] + delta_e[i]
        L = z[i] - 2*delta_max
        U = z[i] + 2*delta_max
        a = np.max([z_inf, z_min, 0.9*L+0.1*z[i]])
        b = np.min([z_sup, z_max, 0.9*U+0.1*z[i]])

        p = ((U-z[i])**2)*(max(S[i], 0.0) + 0.001*abs(S[i]) + 0.5*1e-6*(U-L))
        q = ((z[i]-L)**2)*(max(-S[i], 0.0) + 0.001*abs(S[i]) + 0.5*1e-6*(U-L))
        sqrtp = np.sqrt(p)
        sqrtq = np.sqrt(q)
        zold2[i] = zold1[i]
        zold1[i] = z[i]
        z[i] = max(a, min(b, (L*sqrtp+U*sqrtq)/(sqrtp+sqrtq)))

    line.set_ydata(z)

    fig.canvas.draw()
    fig.canvas.flush_events()

    lRa_old2 = lRa_old1
    lRa_old1 = lRa
    lRa = Ra(x, z)

    _Rq = Rq(z)
    _dRq = dRq_dz(z)

    L1 += _Rsk**2
    Rsk_old1 = Rsk_old2
    Rsk_old2 = _Rsk**2
    _Rsk = Rsk(z)

    print(_z_avg)
    print(_Ra)
    print(_Rsk)
    print("")

input()

