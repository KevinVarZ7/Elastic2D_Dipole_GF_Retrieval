### UNIVERSIDAD NACIONAL AUTÓNOMA DE MÉXICO
### POSGRADO EN CIENCIAS DE LA TIERRA
### DOCTORADO EN CIENCIAS DE LA TIERRA
### AUTHOR: MSC. KEVIN AXEL VARGAS-ZAMUDIO

### PROGRAM: mod_Visualization.py
### 18/04/2022

import numpy as np
# Basic Visualization tools
import matplotlib.pyplot as plt
import seaborn as sns
import pylab as pl
# 3D visualization
import matplotlib.image as mpimg
from   mpl_toolkits.axes_grid1 import make_axes_locatable
from   mpl_toolkits.mplot3d import Axes3D
# SVG format image tools
import matplotlib as mpl
from svgpathtools import svg2paths
from svgpath2mpl import parse_path

def model_geometry(fig,xr,zr,xs,zs,model,source,wave,gam,ninc):
    '''Function: model_geometry
    input: fig: Figure object with geological model
        xr,zr,xs,zs : Float array with receiver and source XZ coordinates
    output: PNG image displaying model geometry '''

    plt.scatter(xr[:],zr[:],color='red',marker='v',s=100,label='Stations')
        
    if source == 'FP' and wave == 'PSV':
        u = np.zeros((ninc))
        v = np.zeros((ninc))
        
        for i in range(ninc):
            
            u[i] = np.sin(gam[i])
            v[i] = np.cos(gam[i])
            print(gam,u,v)
            plt.quiver(xs[:],zs[:],u[i],v[i],color='green',scale=10,label = 'Sources')

    elif source == 'FP' and wave == 'SH':
        shfp_marker = SHFP_marker()
        plt.scatter(xs[:],zs[:],color='green',marker=shfp_marker,s=150,label = 'Sources')

    elif source == 'DP':
        plt.scatter(xs[:],zs[:],color='green',marker='*',s=200,label = 'Source')

    # elif source == 'DP' and wave == 'PSV':
    #     if dipole = 'Gxx':
    #         plt.scatter(xs[:],zs[:],color='green',marker='*',s=100,label='Sources')
    
    # if geo_mod == 'layers':
    #     n = len(h)

    plt.title('Model Geometry: %s_%s ' %(wave,model),fontsize=14)
    plt.xlabel('x[km]')
    plt.ylabel('z[km]')
    plt.gca().invert_yaxis()
    plt.legend()
    plt.grid(alpha=0.5,color='gray')
    fig.savefig('geometry_model_%s_%s.png' %(wave,model),dpi=500,bbox_inches='tight',
        pad_inches = 0.05)
    plt.show()

def model_geometry_mult_source_in_point(xr,zr,xs,zs,model,source,wave,gam,ninc,nfpp):
    '''Function: model_geometry_mult_source_in_point
    input: xr,zr,xs,zs : Float array with receiver and source XZ coordinates
    output: PNG image displaying model geometry for multiple sources in an array point '''

    fig = plt.figure(figsize=(10,10))
    plt.scatter(xr[:],zr[:],color='red',marker='v',s=100,label='Stations')

    if source == 'FP' and wave == 'PSV':
        u = np.zeros((ninc,nfpp))
        v = np.zeros((ninc,nfpp))
        
        for n in range(nfpp):
            u[:,n] = np.sin(gam[n])
            v[:,n] = np.cos(gam[n])
            plt.quiver(xs[:],zs[:],u[:,n],v[:,n],color='green',scale=15,label = 'Sources')
        
    if source == 'FP' and wave == 'SH':
        shfp_marker = SHFP_marker()
        plt.scatter(xs[:],zs[:],color='green',marker=shfp_marker,s=150,label = 'Sources')        
            
    elif source == 'DP' and wave == 'PSV':
        for n in range(nfpp):
            plt.scatter(xs[:],zs[:],color='green',marker='*',s=200,label = 'Source')
            
    elif source == 'DP' and wave == 'SH':
        for n in range(nfpp):
            print('into plot DP mult sh')
            plt.scatter(xs[:],zs[:],color='green',marker='*',s=200,label = 'Source')

    plt.title('Model Geometry: %s_%s ' %(wave,model),fontsize=14)
    plt.xlabel('x[km]')
    plt.ylabel('z[km]')
    plt.gca().invert_yaxis()
    #plt.legend()
    plt.grid(alpha=0.5,color='gray')
    fig.savefig('geometry_model_%s_%s.png' %(wave,model),dpi=500,bbox_inches='tight',
        pad_inches = 0.05)
    plt.show()

def geological_model(mod_geom,properties):
    
    color_layer = ['yellow','gray','black']
    
    fig = plt.figure(figsize=(10,10))
    
    if mod_geom == 'layers':
        ne = len(properties['h'])
        xini = properties['xi']
        xfin = properties['xf']

        x = np.arange(xini,xfin+1)
        y = np.zeros((len(x),ne))
        y0 = np.zeros((len(x)))
        
        for i in range(ne):
            h = properties['h'][i]
            y[:,i] = h
        
            plt.plot(x,y,color = 'brown',linewidth=2.5,linestyle='--')
            if i==0:
                plt.fill_between(x,y0,y[:,i],color=color_layer[i],alpha=0.5,hatch='.')
            else:
                plt.fill_between(x,y[:,i-1],y[:,i],color=color_layer[i],alpha=0.5)
                
    
    return fig
            

    
# Custom marker
def SHFP_marker(*args,**kwargs):
    path = '../AuxImgs/CustomMarkers/'
    #path = '../AuxImgs/CustomMarkers/'
    shfp_path, attributes = svg2paths(path+'SHFP_OutPlane3.svg')
    part1 = shfp_path[0].d()
    part2 = shfp_path[1].d()
    total = part1 + part2
    shfp_marker = parse_path(total)
    shfp_marker = shfp_marker.transformed(mpl.transforms.Affine2D().rotate_deg(180))
    shfp_marker.vertices -= shfp_marker.vertices.mean(axis=0)

    return shfp_marker

def array_geometry(xstat,zstat,x,z,R):
    ninc = len(x)

    fig = plt.figure(figsize=(10,10))
    plt.scatter(x[:],z[:],color='red',marker='*', s=100,label = 'Sources')
    plt.scatter(xstat,zstat,color='blue',marker=7,s=150,label='Stations')

    for j in range(len(xstat)):
        for i in range(ninc):
            xval = [x[i],xstat[j]]
            zval = [z[i],zstat[j]]
            plt.plot(xval,zval,color='black',linestyle='dotted',linewidth=0.75)

    # plt.title(f'Crosscorrelation Array, 2 Receivers, {ninc} Sources',font='Candara',
    #     fontsize=24,fontweight='bold')
    
    plt.title(f'Source Array',font='Candara',
        fontsize=24,fontweight='bold')
    
    plt.xlabel('X km',font='Candara',fontsize=16,weight='bold')
    plt.ylabel('Z km',font='Candara',fontsize=16,weight='bold')
    plt.yticks(fontsize=16)
    plt.xticks(fontsize=16)
    plt.gca().invert_yaxis()
    plt.legend()
    plt.grid(alpha=0.0)
    fig.savefig('Crosscorrelation_Array_%s_Radius_%s.png' %(ninc,R),dpi=300,bbox_inches='tight',pad_inches = 0.05)
    plt.show()
    
    
###############################################################################
def vis_transf_function(est,comp,frec,us,vs,ws):
    
    if comp == 'RTV':
        tit1 = 'Radial'
        tit2 = 'Transversal'
        tit3 = 'Vertical'
    if comp == 'NEV':
        tit1 = 'North - South'
        tit2 = 'East - West'
        tit3 = 'Vertical' 
    
    mpl.rcParams['figure.figsize'] = (15,10)
    fig,ax = plt.subplots(3,1)
    fig.suptitle(f'Transference Functions {comp} components, Soft Layer STATION: {est}',
                      font= 'Candara',fontsize=16,weight = 'bold')
    ax[0].plot(frec,abs(us),color = 'blue', linewidth = 1.2)
    ax[0].set_title(f'{tit1} Component',font='Candara',fontsize=14,weight = 'bold')
    #ax[0].set_xlabel('t[s]',font='Candara',fontsize=12)
    ax[0].set_ylabel('Amplitude',font='Candara',fontsize=12,weight = 'bold')
    ax[0].grid(alpha=0.5)
    ax[1].plot(frec,abs(vs),color = 'green', linewidth = 1.2)
    ax[1].set_title(f'{tit2} Component',font='Candara',fontsize=14,weight = 'bold')
    #ax[1].set_xlabel('t[s]',font='Candara',fontsize=12)
    ax[1].grid(alpha=0.5)
    ax[1].set_ylabel('Amplitude',font='Candara',fontsize=12,weight = 'bold')
    ax[2].plot(frec,abs(ws),color = 'red', linewidth = 1.2)
    ax[2].set_title(f'{tit3} Component',font='Candara',fontsize=14,weight = 'bold')
    ax[2].set_xlabel('Frequency [Hz]',font='Candara',fontsize=12,weight = 'bold')
    ax[2].set_ylabel('Amplitude',font='Candara',fontsize=12,weight = 'bold')
    ax[2].grid(alpha=0.5)


###############################################################################

def vis_synthseis(t,u,v,w,nrec,N,model,gama,fi,pol,comp):
    
    if comp == 'RTV':
        tit1 = 'Radial'
        tit2 = 'Transversal'
        tit3 = 'Vertical'
        
    if comp == 'NEV':
        tit1 = 'North - South'
        tit2 = 'East - West'
        tit3 = 'Vertical' 
        
    scaleu = 2
    scalev = 2
    scalew = 1
    offset = 1
    nrec2 = nrec/2
    
    SisU,ax1 = plt.subplots(figsize=(15,7))
    for i in range(nrec):
        iin = N*i
        iif = N*(i+1)
        ax1.plot(t,u[iin:iif]*scaleu+offset*(i-nrec2),color='blue',lw=1.5, ls='-', marker='s', markersize=0.01)
        ax1.set_xlabel('Tiempo [s]' ,fontsize=12,font = 'Candara')
        ax1.set_ylabel('Estaciones [km]' ,fontsize = 12,font = 'Candara')
        ax1.set_title('Synth Seismograms %s component %s, 'r'$\gamma=%3i$,'r'$\phi=%3i$ 'r'$Pol=%3i$' \
                      %(tit1,model,gama,fi,pol),font = 'Candara',fontsize = 14,weight='bold')
        #ax1.legend(loc=1)
        ax1.grid(color='black', alpha=0.5, linestyle='dashed', linewidth=0.3)
                
        plt.show()
    
        SisU.savefig(f'SynSeis_UxT_{model}_{comp}.png',bbox_inches='tight',pad_inches = 0.05)
    
    
    SisV,ax1 = plt.subplots(figsize=(15,7))
    for i in range(nrec):
        iin = N*i
        iif = N*(i+1)
        ax1.plot(t,v[iin:iif]*scalev+offset*(i-nrec2),color='green',lw=1.5, ls='-', marker='s', markersize=0.01)
        ax1.set_xlabel('Tiempo [s]' ,fontsize=12,font = 'Candara')
        ax1.set_ylabel('Estaciones [km]' ,fontsize = 12,font = 'Candara')
        ax1.set_title('Synth Seismograms %s component %s, 'r'$\gamma=%3i$,'r'$\phi=%3i$ 'r'$Pol=%3i$' \
                      %(tit2,model,gama,fi,pol),font = 'Candara',fontsize = 14,weight='bold')
        #ax1.legend(loc=1)
        ax1.grid(color='black', alpha=0.5, linestyle='dashed', linewidth=0.3)
                
        plt.show()
    
        SisV.savefig(f'SynSeis_VxT_{model}_{comp}.png',bbox_inches='tight',pad_inches = 0.05)
    
    SisW,ax1 = plt.subplots(figsize=(15,7))
    for i in range(nrec):
        iin = N*i
        iif = N*(i+1)
        ax1.plot(t,w[iin:iif]*scalew+offset*(i-nrec2),color='red',lw=1.5, ls='-', marker='s', markersize=0.01)
        ax1.set_xlabel('Tiempo [s]' ,fontsize=12,font = 'Candara')
        ax1.set_ylabel('Estaciones [km]' ,fontsize = 12,font = 'Candara')
        ax1.set_title('Synth Seismograms %s component %s, 'r'$\gamma=%3i$,'r'$\phi=%3i$ 'r'$Pol=%3i$' \
                      %(tit3,model,gama,fi,pol),font = 'Candara',fontsize = 14,weight='bold')
        #ax1.legend(loc=1)
        ax1.grid(color='black', alpha=0.5, linestyle='dashed', linewidth=0.3)
                
        plt.show()
    
        SisW.savefig(f'SynSeis_WxT_{model}_{comp}.png',bbox_inches='tight',pad_inches = 0.05)
