### UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO
### POSGRADO EN CIENCIAS DE LA TIERRA
### DOCTORADO EN CIENCIAS DE LA TIERRA
### AUTHOR: MSC. KEVIN AXEL VARGAS-ZAMUDIO

### PROGRAM: mod_GreenFunctions.py

import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.special import hankel2
from scipy.fft import ifft


def G22_SH_FP(ki,rij,mu):
    kr = ki*rij
    h0kr = hankel2(0,kr)
    g22 = 1/(4*1j*mu) * h0kr

    return g22

def dG22_SH(ki,rij,mu):
    kr = ki*rij
    h1kr = hankel2(1,kr)
    dg22 = ki*1j/(4*mu) * h1kr

    return dg22

def dG22_SH_mjk(ki,rij,mu):
    kr = ki*rij
    h1kr = hankel2(1,kr)
    dg22 = ki*1j/(4*mu) * h1kr * mu
    
    return dg22
    
def Gij_PSV(ks,kp,rij,g,C):
    kpr = kp*rij ; ksr = ks*rij
    h0P = hankel2(0,kpr)
    h0S = hankel2(0,ksr)
    h2P = hankel2(2,kpr)
    h2S = hankel2(2,ksr)

    if C[5,5] != 0:
        A = h0P/C[0,0] + h0S/C[5,5]
        B = h2P/C[0,0] - h2S/C[5,5]
    else:
        A = h0P/C[0,0]
        B = h2P/C[0,0]
    
    G = np.zeros((2,2),dtype=complex)
    d = np.eye(2)

    for j in range(2):
        G[0,j] = (1/8j)*(d[0,j]*A-(2*g[0]*g[j]-d[0,j])*B)

    i=1 ; j=1
    G[i,j] = (1/8j) * (d[i,j]*A-(2*g[i]*g[j]-d[i,j])*B)
    G[1,0] = G[0,1]

    return G

def Gij_PSV_2(ks,kp,rij,g,C,alph,bet): #-------------------------- Kausel,2006
    kpr = kp*rij ; ksr = ks*rij
    h0S = hankel2(0,ksr)
    h1S = hankel2(1,ksr)
    h1P = hankel2(1,kpr)
    h2P = hankel2(2,kpr)
    h2S = hankel2(2,ksr)
    
    A = (1j/4) * ((h1S/ksr - (bet/alph)**2 * h1P/kpr) - h0S)
    B = (1j/4) * ((bet/alph)**2 * h2P - h2S)

    G = np.zeros((2,2),dtype=complex)
    d = np.eye(2)
    
    G[0,0] = (1/C[5,5]) * (A + B*g[0]**2)
    G[1,1] = (1/C[5,5]) * (A + B*g[1]**2)
    
    G[0,1] = (1/C[5,5]) * (B*g[0]*g[1])
    G[1,0] = (1/C[5,5]) * (B*g[1]*g[0])

    return G

def Gxx_DP_3(ks,kp,rij,g,C,alph,bet):   # Kausel's formula dgij/dxk pag 38
    kpr = kp*rij 
    ksr = ks*rij
    h0S = hankel2(0,ksr)
    h1S = hankel2(1,ksr)
    h1P = hankel2(1,kpr)
    h2S = hankel2(2,ksr)
    h2P = hankel2(2,kpr)
    
    A = (1j/4)*((h1S/ksr - (bet/alph)**2 * h1P/kpr)-h0S)
    B = (1j/4)*((bet/alph)**2 * h2P - h2S)
    
    dA = (B/rij) + (1j/(4*rij)*ksr*h1S)
    dB = (1j/(4*rij))*((bet/alph)**2 * kpr*h1P - ksr*h1S) - 2*(B/rij)
    
    Gxx = np.zeros((2),dtype=complex)

    Gxx[0] = -(1/C[5,5]) * (g[0]* (1 + g[0]**2)*dA + 2*(B/rij) * (1-g[1]**2))
    Gxx[1] = -(1/C[5,5]) * g[1]*(dA*g[0]**2 + (B/rij) * (1-2*g[0]**2))
    
    return Gxx

def Gxx_DP_2(ks,kp,rij,g,Ci,rho): 
    kpr = kp*rij
    ksr = ks*rij
    h1S = hankel2(1,ksr) ; h1P = hankel2(1,kpr)
    h2S = hankel2(2,ksr) ; h2P = hankel2(2,kpr)
    h3S = hankel2(3,ksr) ; h3P = hankel2(3,kpr)

    
    B = h2P/Ci[0,0] - h2S/Ci[5,5]
    C = h1S/Ci[5,5] - h1P/Ci[0,0]
    D = h3S/Ci[5,5] - h3P/Ci[0,0]
   
    Gxx = np.zeros((2),dtype=complex)
    
    Gxx[0] = -1/8j/rho * g[0] * (C - 2*B/rij + D*(1 - 2*g[0]**2))
    Gxx[1] = 1/4j/rho * g[1] * (B/rij + g[0]**2*D)
    return Gxx

def Gzz_DP_2(ks,kp,rij,g,Ci,rho):
    kpr = kp*rij ; ksr = ks*rij
    h1S = hankel2(1,ksr) ; h1P = hankel2(1,kpr)
    h2S = hankel2(2,ksr) ; h2P = hankel2(2,kpr)
    h3S = hankel2(3,ksr) ; h3P = hankel2(3,kpr)
    B = h2P/Ci[0,0] - h2S/Ci[5,5]
    C = h1S/Ci[5,5] - h1P/Ci[0,0]
    D = h3S/Ci[5,5] - h3P/Ci[0,0]

    Gzz = np.zeros((2),dtype=complex)
    
    Gzz[0] = 1/4j/rho * g[0] * (B/rij + g[1]**2 * D)
    Gzz[1] = -1/8j/rho * g[1] * (C - 2*B/rij + D*(1 - 2*g[1]**2))
    
    return Gzz

def Gxx_DP_pol(ks,kp,rij,g,C,alph,bet,gamma):
    kpr = kp*rij 
    ksr = ks*rij
    h0S = hankel2(0,ksr)
    h1S = hankel2(1,ksr)
    h1P = hankel2(1,kpr)
    h2S = hankel2(2,ksr)
    h2P = hankel2(2,kpr)
    gam = 2*gamma
    A = (1j/4)*((h1S/ksr - (bet**2/alph**2)*h1P/kpr)-h0S)
    B = (1j/4)*((bet**2/alph**2)*h2P - h2S)
    dA = (B/rij) + (1j/(4*rij)*ksr*h1S)
    dB = (1/(4*rij))*((bet**2/alph**2)*kpr*h1P - ksr*h1S) - 2*(B/rij)
    Gxx = np.zeros((2),dtype=complex)
    Gxx[0] = -(1/2/C[5,5]) * (dA + dB + (B/rij) + np.cos(gam)*(dA + dB - (B/rij)))
    Gxx[1] = -(1/2/C[5,5]) * -(np.sin(gam)*(dB + (B/rij)))
    
    return Gxx

def Gxx_DP(ks,kp,rij,g,C,alph,bet): #   KAUSEL
    kpr = kp*rij 
    ksr = ks*rij
    h0S = hankel2(0,ksr)
    h1S = hankel2(1,ksr)
    h1P = hankel2(1,kpr)
    h2S = hankel2(2,ksr)
    h2P = hankel2(2,kpr)

    A = (1j/4)*((h1S/ksr - (bet**2/alph**2)*h1P/kpr)-h0S)
    B = (1j/4)*((bet**2/alph**2)*h2P - h2S)

    dA = (B/rij) + (1j/(4*rij)*ksr*h1S)
    dB = (1j/(4*rij))*((bet/alph)**2 * kpr*h1P - ksr*h1S) - 2*(B/rij)

    Gxx = np.zeros((2),dtype=complex)

    Gxx[0] = -(1/C[5,5]) * (g[0]*(dA + dB*g[0]**2 + 2*(B/rij)*g[1]**2))
    Gxx[1] = -(1/C[5,5]) * g[1]*((dB - 2*(B/rij))*g[0]**2 + B/rij)

    return Gxx

def Gzz_DP(ks,kp,rij,g,C,alph,bet):
    kpr = kp*rij ; ksr = ks*rij
    h0S = hankel2(0,ksr)
    h1S = hankel2(1,ksr)
    h1P = hankel2(1,kpr)
    h2S = hankel2(2,ksr)
    h2P = hankel2(2,kpr)

    A = (1j/4)*((h1S/ksr - (bet**2/alph**2)*h1P/kpr)-h0S)
    B = (1j/4)*((bet**2/alph**2)*h2P - h2S)

    dA = (B/rij) + (1j/(4*rij)*ksr*h1S)
    dB = (1j/(4*rij))*((bet/alph)**2 *kpr*h1P - ksr*h1S) - 2*(B/rij)

    Gzz = np.zeros((2),dtype=complex)

    Gzz[0] = -(1/C[5,5]) * g[0]*((dB - 2*(B/rij))*g[1]**2 + B/rij)
    Gzz[1] = -(1/C[5,5]) * (g[1]*(dA + dB*g[1]**2 + 2*(B/rij)*g[0]**2))

    return Gzz


def Gzx_DP(ks,kp,rij,g,C,alph,bet):
    kpr = kp*rij ; ksr = ks*rij
    h0S = hankel2(0,ksr)
    h1S = hankel2(1,ksr)
    h1P = hankel2(1,kpr)
    h2S = hankel2(2,ksr)
    h2P = hankel2(2,kpr)

    A = (1j/4)*((h1S/ksr - (bet**2/alph**2)*h1P/kpr)-h0S)
    B = (1j/4)*((bet**2/alph**2)*h2P - h2S)

    dA = (B/rij) + (1j/(4*rij)*ksr*h1S)
    dB = (1j/(4*rij))*((bet**2/alph**2)*kpr*h1P - ksr*h1S) - 2*(B/rij)

    Gzx = np.zeros((2),dtype=complex)

    Gzx[0] = -(1/C[5,5]) * g[1]*((dB - 2*(B/rij))*g[0]**2 + (B/rij))
    Gzx[1] = -(1/C[5,5]) * g[0]*(dA + (dB - 2*(B/rij))*g[1]**2)

    return Gzx

def Gxz_DP(ks,kp,rij,g,C,alph,bet):
    kpr = kp*rij ; ksr = ks*rij
    h0S = hankel2(0,ksr)
    h1S = hankel2(1,ksr)
    h1P = hankel2(1,kpr)
    h2S = hankel2(2,ksr)
    h2P = hankel2(2,kpr)

    A = (1j/4)*((h1S/ksr - (bet**2/alph**2)*h1P/kpr)-h0S)
    B = (1j/4)*((bet**2/alph**2)*h2P - h2S)

    dA = (B/rij) + (1j/(4*rij)*ksr*h1S)
    dB = (1j/(4*rij))*((bet**2/alph**2)*kpr*h1P - ksr*h1S) - 2*(B/rij)

    Gxz = np.zeros((2),dtype=complex)

    Gxz[0] = -(1/C[5,5]) * g[1] * (dA + (dB - 2*(B/rij))*g[0]**2)
    Gxz[1] = -(1/C[5,5]) * g[0] * ((dB - 2*(B/rij))*g[1]**2 + (B/rij))

    return Gxz


def dgxxk_dxk_ux(ks,kp,rij,g,C,alph,bet,mij):
    
    kpr = kp*rij ; ksr = ks*rij
    h0S = hankel2(0,ksr)
    h1S = hankel2(1,ksr)
    h1P = hankel2(1,kpr)
    h2S = hankel2(2,ksr)
    h2P = hankel2(2,kpr)
    
    mu = C[5,5]

    A = (1j/4)*((h1S/ksr - (bet**2/alph**2)*h1P/kpr)-h0S)
    B = (1j/4)*((bet**2/alph**2)*h2P - h2S)

    dA = (B/rij) + (1j/(4*rij)*ksr*h1S)
    dB = (1j/(4*rij))*((bet**2/alph**2)*kpr*h1P - ksr*h1S) - 2*(B/rij)
    

    dgxx_dx = -g[0]/mu * (dA + dB*g[0]**2 + (2*B/rij) * g[1]**2) 
    dgxx_dz = -g[1]/mu * (dA + (dB - (2*B/rij)*g[0]**2))         
    
    dgxz_dx = -g[1]/mu * ((dB - (2*B/rij))*g[0]**2 + (B/rij)) 
    dgxz_dz = -g[0]/mu * ((dB - (2*B/rij))*g[1]**2 + (B/rij)) 
    
    dgxxkux = np.zeros((4),dtype=complex)
    
    dgxxkux[0] = dgxx_dx
    dgxxkux[1] = dgxx_dz
    dgxxkux[2] = dgxz_dx
    dgxxkux[3] = dgxz_dz

    return dgxxkux

def dgzxk_dxk_uz(ks,kp,rij,g,C,alph,bet,mij):
    
    kpr = kp*rij ; ksr = ks*rij
    h0S = hankel2(0,ksr)
    h1S = hankel2(1,ksr)
    h1P = hankel2(1,kpr)
    h2S = hankel2(2,ksr)
    h2P = hankel2(2,kpr)
    
    mu = C[5,5]

    A = (1j/4)*((h1S/ksr - (bet**2/alph**2)*h1P/kpr)-h0S)
    B = (1j/4)*((bet**2/alph**2)*h2P - h2S)

    dA = (B/rij) + (1j/(4*rij)*ksr*h1S)
    dB = (1j/(4*rij))*((bet**2/alph**2)*kpr*h1P - ksr*h1S) - 2*(B/rij)    

    dgzx_dx = -g[1]/mu * ((dB - (2*B/rij))*g[0]**2 + (B/rij)) #* mij[0,0]
    dgzx_dz = -g[0]/mu * ((dB - (2*B/rij))*g[1]**2 + (B/rij)) #* mij[1,0]
    
    dgzz_dx = -g[0]/mu * (dA + (dB - (2*B/rij)*g[1]**2))         #* mij[0,1]
    dgzz_dz = -g[1]/mu * (dA + dB*g[1]**2 + (2*B/rij) * g[0]**2) #* mij[1,1]

    dgxxkuz = np.zeros((4),dtype=complex)
    
    dgxxkuz[0] = dgzx_dx
    dgxxkuz[1] = dgzx_dz
    dgxxkuz[2] = dgzz_dx
    dgxxkuz[3] = dgzz_dz

    
    return dgxxkuz

def Gij_k_PSV(ks,kp,rij,g,C,alph,bet,mij):
    
    kpr = kp*rij ; ksr = ks*rij
    h0S = hankel2(0,ksr)
    h1S = hankel2(1,ksr)
    h1P = hankel2(1,kpr)
    h2S = hankel2(2,ksr)
    h2P = hankel2(2,kpr)
    
    mu = C[5,5]

    A = (1j/4)*((h1S/ksr - (bet**2/alph**2)*h1P/kpr)-h0S)
    B = (1j/4)*((bet**2/alph**2)*h2P - h2S)

    dA = (B/rij) + (1j/(4*rij)*ksr*h1S)
    dB = (1j/(4*rij))*((bet**2/alph**2)*kpr*h1P - ksr*h1S) - 2*(B/rij)

    Gijk = np.zeros((2,2,2),dtype=complex)    
    
    Gijk[0,0,0] = -g[0]/mu * (dA + dB*g[0]**2 + (2*B/rij) * g[1]**2)            #dGxx_dx 
    Gijk[0,1,0] = -g[1]/mu * ((dB - (2*B/rij))*g[0]**2 + (B/rij))               #dGxz_dx
    Gijk[1,0,0] = -g[1]/mu * ((dB - (2*B/rij))*g[0]**2 + (B/rij))               #dGzx_dx
    Gijk[1,1,0] = -g[0]/mu * (dA + (dB - (2*B/rij)*g[1]**2))                    #dGzz_dx
    
    Gijk[0,0,1] = -g[1]/mu * (dA + (dB - (2*B/rij)*g[0]**2))                    #dgxx_dz
    Gijk[0,1,1] = -g[0]/mu * ((dB - (2*B/rij))*g[1]**2 + (B/rij))               #dgxz_dz
    Gijk[1,0,1] = -g[0]/mu * ((dB - (2*B/rij))*g[1]**2 + (B/rij))               #dgzx_dz
    Gijk[1,1,1] = -g[1]/mu * (dA + dB*g[1]**2 + (2*B/rij) * g[0]**2)            #dgzz_dz
    
    return Gijk

def My_TM(ks,kp,rij,g,C,alph,bet):
    
    kpr = kp*rij ; ksr = ks*rij
    h0S = hankel2(0,ksr)
    h1S = hankel2(1,ksr)
    h1P = hankel2(1,kpr)
    h2S = hankel2(2,ksr)
    h2P = hankel2(2,kpr)
    
    A = (1j/4)*((h1S/ksr - (bet**2/alph**2)*h1P/kpr)-h0S)
    B = (1j/4)*((bet**2/alph**2)*h2P - h2S)

    dA = (B/rij) + (1j/(4*rij)*ksr*h1S)
    dB = (1j/(4*rij))*((bet**2/alph**2)*kpr*h1P - ksr*h1S) - 2*(B/rij)
    
    My = np.zeros((2),dtype=complex)
    
    My[0] = -(1/C[5,5]) * g[1] * (dA + 2*dB*g[0]**2 - 4*B/rij*g[0]**2 + B/rij)
    My[1] = -(1/C[5,5]) * g[0] * (dA + 2*dB*g[1]**2 - 4*B/rij*g[1]**2 + B/rij)
    
    return My

def correct_spectre_ricker(nfrec,df,Rtp,delay):
    dt = (1/(df*2*nfrec)*(nfrec/(nfrec-1)))
    tp = 0.5*Rtp
    n = np.arange(1,2*nfrec-1,1)

    rick = np.zeros((len(n)))
    crick = np.zeros((len(n)),dtype=complex)
    cr = np.zeros(len(n))

    rick = ricker(n,dt,tp)
    crick[:] = rick[:] + 1j

    cr = 2*(2*nfrec-2) * ifft(0.5*rick)

    spec = np.zeros((nfrec))
    spec = cr[0:nfrec]
    spec[0] = 0
    
    w = np.zeros((nfrec))
    w[:] = [2*np.pi*(df*((i))) for i in range(nfrec)]

    spec *= np.exp(1j*w*delay)

    return np.real(spec)


def ricker(n,dt,tp):
    n0 = len(n)
    a = np.zeros((n0))
    ric = np.zeros((n0))

    for i in range(n0):
        a[i] = (np.pi * (dt*(n[i]-1)/tp))**2
        ric[i] = (0.5-a[i])*np.exp(-a[i])
        a[i] = (np.pi*(dt*(n0+n[i]-1))/tp)**2
        ric[i] += (0.5-a[i])*np.exp(-a[i])
        a[i] = (np.pi*(dt*(-n0+n[i]-1))/tp)**2
        ric[i] += (0.5-a[i])*np.exp(-a[i])

    return ric[:]