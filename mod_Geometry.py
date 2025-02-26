# -*- coding: utf-8 -*-
"""
UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO
POSGRADO EN CIENCIAS DE LA TIERRA
DOCTORADO EN CIENCIAS DE LA TIERRA
Created on Mon May  9 13:36:58 2022

@author: M.Sc. Kevin Axel Vargas-Zamudio
PROGRAM: mod_Geometry.py

"""
import numpy as np
import matplotlib.pyplot as plt

def source_geometry(ninc,geom,h0,k0):
    if geom == 'Circ':
        # R = []
        # for i in range(ngeom):
        #     R.append(int(input(f'Type {i+1} array Radius: \n')))
        
        R = 31.819 #float(input('Type the circumference radius: \n'))
            
        # Circumference center coordinates

        h = np.mean(h0)
        k = np.mean(k0)
        g2rad = np.pi / 180.0
        
        theta = np.linspace(0,360,ninc)
        theta *= g2rad

        
        x = np.zeros((ninc))
        z = np.zeros((ninc))
        
        x[:] = [R*np.cos(theta[i]) + h for i in range(ninc)]
        z[:] = [R*np.sin(theta[i]) + k for i in range(ninc)]
        
        # arch = [f'farfield_'+str(ninc)+'.dat',
        # 'middlefield_'+str(ninc)+'.dat',
        # 'nearfield_'+str(ninc)+'.dat']
        
        arch = f'receivers_{ninc}_{geom}_Rad_{R:.2f}.dat'
        
        # for i in range(3):
        #     ofile = open(arch[i],'w')
        
        #     for j in range(ninc):
        #         #ofile.write(str(x[j,i])+'\t'+str(y[j,i]) +'\n')
        #         coord = '{0:3.3f} \t {1:3.3f} \n'.format(x[j,i],y[j,i])
        #         ofile.write(coord)
        #     ofile.close()
        
        # opc = input('Add source position randomness? y/n \n')
        # if opc == 'y' or opc == 'Y':
        #     np.random.seed(7)
        #     x += (max(x)-min(x))/(R)* np.random.randn(ninc)
            # z += (max(z)-min(z))/(R)* np.random.randn(ninc)

        
        # ofile = open(arch,'w')
        # for j in range(ninc):
        #     coord = '{0:3.3f} \t {1:3.3f} \n'.format(x[j],z[j])
        #     ofile.write(coord)
            
        # ofile.close()
        
        return x,z,R
        
    elif geom == 'Rect':
        
        xv = [-50,50,50,-50]
        zv = [10,10,110,110]
        
        nv = len(xv)
        npl = int(ninc/nv)
        nseg = npl + 1
        
        xseg = np.zeros((nv,npl))
        zseg = np.zeros((nv,npl))
        d = np.zeros((nv))
        
        #d[0] = np.sqrt((xv[1] - xv[0])**2 + ((zv[0] - zv[0])**2))
        
        for i in range(nv):
            if i != nv-1:
                d[i] = np.sqrt((xv[i+1] - xv[i])**2 + ((zv[i+1] - zv[i])**2))
                if (zv[i+1] - zv[i])**2 == 0:
                    
                    segx = d[i]/nseg
                    segz = 0
                    auxx = 0
                    auxz = 0
                    for j in range(npl):
                        auxx += segx
                        if i == 0:
                            xseg[i,j] = xv[i] + auxx
                        elif i == 2:
                            xseg[i,j] = xv[i] - auxx
                        zseg[i,j] = zv[i]
                        
                elif (xv[i+1] - xv[i])**2 == 0:
                    segz = d[i]/nseg
                    segx = 0
                    auxx = 0
                    auxz = 0
                    for j in range(npl):
                        auxz += segz
                        zseg[i,j] = zv[i] + auxz
                        xseg[i,j] = xv[i]
        
            else:
                d[i] = np.sqrt((xv[0] - xv[i])**2 + ((zv[0] - zv[i])**2))
                if (zv[0] - zv[i])**2 == 0:
                    segx = d[i]/nseg
                    segz = 0
                    auxx = 0
                    auxz = 0
                    for j in range(npl):
                        auxx += segx
                        xseg[i,j] = xv[i] + auxx
                        zseg[i,j] = zv[i]
                        
                elif (xv[0] - xv[i])**2 == 0:
                    segz = d[i]/nseg
                    segx = 0
                    auxx = 0
                    auxz = 0
                    for j in range(npl):
                        auxz += segz
                        zseg[i,j] = zv[i] - auxz
                        xseg[i,j] = xv[i]
        
        ntot = nv + npl*nv
        xtot = np.zeros((ntot+1))
        ztot = np.zeros((ntot+1))
        
        # Filling Vertices
        for i in range(0,ntot-1,nseg):
            #print(i)
            xtot[i] = xv[int(i/nseg)]
            ztot[i] = zv[int(i/nseg)]
        
        # Filling x.z segment coordinates
        for i in range(nv):
            for j in range(npl):
                xtot[i*nseg + (j+1)] = xseg[i,j]
                ztot[i*nseg + (j+1)] = zseg[i,j] 
        
        # Completing the polygon
        
        xtot[ntot] = xv[0]
        ztot[ntot] = zv[0]
        
        # arch = f'receivers_{ninc}_{geom}.dat'
        # ofile = open(arch,'w')
        # for j in range(ntot):
        #     coord = '{0:3.3f} \t {1:3.3f} \n'.format(xtot[j],ztot[j])
        #     ofile.write(coord)
             
        # ofile.close()
        
                
        return xtot,ztot,ntot+1

def receiver_geometry(xs,zs,nrec,geom):
    
    ''' Surface arrays for X coordinate distribution when  Z = 0'''
    if geom == 'Surface_Array':
        #xini = float(input('Enter initial X: \n'))
        #xfin = float(input('Enter final X: \n'))
        xini = -30.0
        xfin = 30.0
        dx = (xfin - xini)/nrec
        
        xrec = np.zeros((nrec))
        xrec = np.arange(xini,xfin,nrec)
        xrec = np.linspace(xini,xfin,nrec)
        # for i in range(nrec):
        #     xrec[i] = xini + dx*i
        #print(xrec)
        return xrec
    
    elif geom == 'Circ':
        
        R = int(input('Type the circumference radius: \n'))
            
        # Circumference center coordinates
            
        h = xs
        k = -zs
        g2rad = np.pi / 180.0
        
        theta = np.linspace(0,360,nrec)
        theta *= g2rad
            
        x = np.zeros((nrec))
        z = np.zeros((nrec))

        
        x[:] = [R*np.cos(theta[i]) - h for i in range(nrec)]
        z[:] = [R*np.sin(theta[i]) - k for i in range(nrec)]
              
        arch = f'receivers_{nrec}_{geom}_Rad_{R:.2f}.dat'
       
        # ofile = open(arch,'w')
        # for j in range(ninc):
        #     coord = '{0:3.3f} \t {1:3.3f} \n'.format(x[j],z[j])
        #     ofile.write(coord)
            
        # ofile.close()
        
        return x,z,R     


    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        