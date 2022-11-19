#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt

Ra_max = 1.5
ABS_EPS = 1e-15

_Ra = 0
_z_avg = 0

def smooth_abs(x):
    return np.sqrt(x**2 + ABS_EPS)

def smooth_abs_deriv(x):
    return np.divide(x, np.sqrt(np.power(x, 2) + ABS_EPS))

def z_avg(z):
    global _z_avg

    _z_avg = np.sum(z)/len(z)

def Ra(z):
    global _Ra

    z_avg(z)

    _Ra = np.sum(smooth_abs(z-_z_avg))/len(z)
    return (_Ra - Ra_max)**2

def dRa_dz(z):
    return 2*(_Ra - Ra_max)*smooth_abs_deriv(z - _z_avg)*(1 - 1/len(z))

LEN = 500
   
delta_max = 1e10
delta_min = 1e-30
delta = 1

z_min = -4
z_max = 4
rng = np.random.Generator(np.random.PCG64())
z = rng.uniform(z_min, z_max, LEN+1)
x = np.linspace(0, LEN, LEN+1)
zold1 = z
zold2 = z
delta_e = np.ones_like(x)*delta

lRa = Ra(z)
lRa_old1 = lRa
lRa_old2 = lRa

ETA0 = 1
ETA1 = 1
L0 = lRa
L1 = 0

plt.ion()

fig, ax = plt.subplots()
ax.set_xlabel("x")
ax.set_ylabel("z")
line, = plt.plot(x, z)

plt.show()

while abs(lRa) > 1e-6:
    S = ETA0*dRa_dz(z)

    Ra_e = (lRa-lRa_old1)*(lRa_old1-lRa_old2)

    for i in range(len(z)):
        d_e = (z[i] - zold1[i])*(zold1[i] - zold2[i])
        if d_e < 0 or Ra_e < 0:
            delta_e[i] = max(0.1*delta_e[i], delta_min)
        elif d_e > 0:
            delta_e[i] = min(1.1*delta_e[i], delta_max)

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

    L0 += lRa
    lRa_old2 = lRa_old1
    lRa_old1 = lRa
    lRa = Ra(z)

    print(_z_avg)
    print(_Ra)
    print("")

print(_z_avg)
print(_Ra)

input()

