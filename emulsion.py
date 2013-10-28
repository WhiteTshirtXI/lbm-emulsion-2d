import os
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import findlist
import circ
import vtk
import time
os.system("cls")
#Nombre archivo para guardar las imagenes
nombre = 'Cortante_lamda3.5_vel0.5_sep28'
#Numero de iteraciones
maxT   = 5000
#Frecuencia de imagen y guardar resultados
itPlot = 50
#Constantes generales
ly = 100
lx = 100
#Amplitud de la fuerza de interaccion
G = -1.5
#Velocidad de las paredes (flujo cortante)
UW = 0
#Densidad del fluido 1
rho1_0 = 1
#Densidad del fluido 2
rho2_0 = 1
#Numero de gotas
numdrops = 1
#Coordenadas de las gotas
coordinates = np.array([[-20,0,15],[20, 14, 15]]) #Coordenadas
#Relaxation parameter for fluid 1
omega1 = 1
#Relaxation parameter for fluid 2
omega2 = 1.3
#Numero de Reynolds
Re = 0.16

##Parametros por defecto
#Tao fluido 1
tao1 = 1/omega1
#Tao fluido 1
tao2 = 1/omega2
#Viscosidad de lattice boltzmann fluido 1
vlb1 = (tao1-0.5)*(1./3)
#Viscosidad de lattice boltzmann fluido 2
vlb2 = (tao2-0.5)*(1./3)
#Delta de distancia
dx = 1./lx
#Delta de tiempo
dt = Re*(dx**2)*vlb1
#Gravedad
grav = 9.81*(dt**2)/dx

##Comienzo de programa

obst_x =np.zeros([numdrops,1])
obst_y =np.zeros([numdrops,1])
obst_r =np.zeros([numdrops,1])
 
for i in xrange(numdrops):
    #Posicion de las gotas
    obst_x[i,0] = coordinates[i,0]   
    obst_y[i,0] = coordinates[i,1]
    #Radio
    obst_r[i,0] = coordinates[i,2]   

obst = np.zeros([lx+1,ly+1])
# Location of top/bottom boundary
obst[:,0] = 1
### Location of top/bottom boundary
obst[:,ly] = 1
bbRegion = findlist.find(obst,1) # Boolean mask for bounce-back cells
bbRegiona = list(bbRegion[:,0])
bbRegionb = list(bbRegion[:,1])
del(bbRegion)

Gomega1 = G/omega1
Gomega2 = G/omega2
counter = 0
itsaved = 0
paso = 1
 
# D2Q9 LATTICE CONSTANTS
tNS   = np.mat([4./9,1./9,1./9,1./9,1./9,1./36,1./36,1./36,1./36])   #weigh factors wi
cxNS  = np.mat([0,1,0,-1,0,1,-1,-1,1])   #velocity ex
cyNS  = np.mat([0,0,1,0,-1,1,1,-1,-1])   #velocity ey
oppNS = np.mat([0,3,4,1,2,7,8,5,6])   #opposing node
 
[x,y] = np.mgrid[-lx/2:lx/2+1,-ly/2:ly/2+1]
 
jx2 = np.zeros([lx+1,ly+1])
jy2 = np.zeros([lx+1,ly+1])
 
# INITIAL CONDITION: Velocity
jx1 = np.zeros([lx+1,ly+1])
jy1 = np.zeros([lx+1,ly+1])
# Density distribution
fIn= np.zeros([9,lx+1,ly+1])
gIn= np.zeros([9,lx+1,ly+1])

##Caso en que no se tienen 2 gotas sino una distribucion aleatoria
#drho = 0.001
#import random
#rand = np.zeros([lx+1,ly+1])
#for i in xrange(lx):
#    for j in xrange(ly):
#        rand[i,j] = random.random()
#        
#delta_rho = -drho*(1-2.0*rand)

### INITIAL CONDITION FOR BOTH DISTRIBUTION FUNCTIONS: (T=0) ==> TIn(i) = t(i)
##for i in xrange(9):
##    fIn[i,:,:] = tNS[0,i]*(1.0 + delta_rho)
##    gIn[i,:,:] = tNS[0,i]*(1.0 - delta_rho)

# Calculates density distribution
for i in xrange(9):
    cuNS1 = 3*(cxNS[0,i]*jx1+cyNS[0,i]*jy1)
    cuNS2 = 3*(cxNS[0,i]*jx2+cyNS[0,i]*jy2)
    for j in xrange(lx+1):
        for k in xrange(ly+1):
            t = 0
            for z in xrange(numdrops):
                if (abs((x[j,k]-obst_x[z]))**2+abs((y[j,k]-obst_y[z]))**2)<((obst_r[z])**2) or (abs((x[j,k]-obst_x[z]))**2+abs((y[j,k]-obst_y[z]))**2)==((obst_r[z])**2):
                    fIn[i,j,k] = 0
                    gIn[i,j,k] = rho2_0*tNS[0,i]*(1 + cuNS2[j,k]+(1./2)*(cuNS2[j,k]*cuNS2[j,k])-(3./2)*(jx2[j,k]**2+jy2[j,k]**2))
                    t = t+1
                elif (abs((x[j,k]-obst_x[z]))**2+abs((y[j,k]-obst_y[z]))**2) > ((obst_r[z])**2) and (t<1):
                    fIn[i,j,k] = rho1_0*tNS[0,i]*( 1 + cuNS1[j,k] + 1./2*(cuNS1[j,k]*cuNS1[j,k])-(3./2)*(jx1[j,k]**2+jy1[j,k]**2))
                    gIn[i,j,k] = 0
                    t = t+1

t = 0
## MAIN LOOP (TIME CYCLES)
for cycle in xrange(paso,maxT+1):
        tic = time.clock()
    ##  MACROSCOPIC VARIABLES
        counter = counter+1
        rho1 = sum(fIn).reshape(lx+1,ly+1)
        rho2 = sum(gIn).reshape(lx+1,ly+1)
    ##  e*f
        temp1 = (((fIn).reshape(9,(lx+1)*(ly+1))))
        temp1 = np.mat(temp1)
        tempf = cxNS*temp1
        jx1 = tempf.reshape(lx+1,ly+1) 
        del(tempf)
        tempf = cyNS*temp1
        jy1 = tempf.reshape(lx+1,ly+1)
        del(tempf)
        del(temp1)
        temp2 = (((gIn).reshape(9,(ly+1)*(lx+1))))
        temp2=np.mat(temp2)
        tempf=cxNS*temp2
        jx2 = tempf.reshape(lx+1,ly+1)
        del(tempf)
        tempf=cyNS*temp2
        jy2 = tempf.reshape(lx+1,ly+1) + (dt/2)*((rho2-rho1))*grav #Este ultimo termino es de fuerza necesario para la inclusion de gravedad (ver Fuerzas Lattice)
        del(tempf)
        del(temp2)

        rhoTot_OMEGA = rho1*omega1 + rho2*omega2 
        uTotX = (jx1*omega1+jx2*omega2) / rhoTot_OMEGA # Equation 95 Sukop
        uTotY = (jy1*omega1+jy2*omega2) / rhoTot_OMEGA
      
        rhoContrib1x = np.zeros([rho1.shape[0],rho1.shape[1]])
        rhoContrib2x = np.zeros([rho2.shape[0],rho2.shape[1]])
        
        rhoContrib1y = np.zeros([rho1.shape[0],rho1.shape[1]])
        rhoContrib2y = np.zeros([rho2.shape[0],rho2.shape[1]])
   
    ##  Calculates the force therm (sum) explained on chapter 7 of Sukop Eq 98
        for i in xrange(9):
            rhoContrib1x = rhoContrib1x + circ.circshift(rho1*tNS[0,i],cxNS[0,i],cyNS[0,i])*cxNS[0,i]
            rhoContrib1y = rhoContrib1y + circ.circshift(rho1*tNS[0,i],cxNS[0,i],cyNS[0,i])*cyNS[0,i]
            
            rhoContrib2x = rhoContrib2x + circ.circshift(rho2*tNS[0,i],cxNS[0,i],cyNS[0,i])*cxNS[0,i]
            rhoContrib2y = rhoContrib2y + circ.circshift(rho2*tNS[0,i],cxNS[0,i],cyNS[0,i])*cyNS[0,i]

##    Total velocity with force term
        uTotX1 = uTotX - Gomega1*rhoContrib2x #POTENTIAL CONTRIBUTION OF FLUID 2 ON 1  Equation 100 Sukop
        uTotY1 = uTotY - Gomega1*rhoContrib2y 

        uTotX2 = uTotX - Gomega2*rhoContrib1x #POTENTIAL CONTRIBUTION OF FLUID 2 ON 1
        uTotY2 = uTotY - Gomega2*rhoContrib1y      

    #  COLLISION STEP FLUID 1 AND 2
        fEq = np.zeros([9,lx+1,ly+1])
        gEq = np.zeros([9,lx+1,ly+1])
        fOut = np.zeros([9,lx+1,ly+1])
        gOut = np.zeros([9,lx+1,ly+1])
        fg2 = np.zeros([9,lx+1,ly+1])
        for i in xrange(9):
           cuNS1 = 3*(cxNS[0,i]*uTotX1+cyNS[0,i]*uTotY1)
           cuNS2 = 3*(cxNS[0,i]*uTotX2+cyNS[0,i]*uTotY2)
           fEq[i,:,:] = (rho1*tNS[0,i]) + np.multiply((rho1*tNS[0,i]),cuNS1) + np.multiply((rho1*tNS[0,i]),(0.5*(np.multiply(cuNS1,cuNS1))-(1.5*(np.multiply(uTotX1,uTotX1) + np.multiply(uTotY1,uTotY1)))))
           gEq[i,:,:] = (rho2*tNS[0,i]) + np.multiply((rho2*tNS[0,i]),cuNS2) + np.multiply((rho2*tNS[0,i]),(0.5*(np.multiply(cuNS2,cuNS2))-(1.5*(np.multiply(uTotX2,uTotX2) + np.multiply(uTotY2,uTotY2)))))               
           fOut[i,:,:] = fIn[i,:,:] - omega1*(fIn[i,:,:]-fEq[i,:,:])
           fOut[i,:,:] = fOut[i,:,:] 
           fg2[i,:,:] = np.multiply((1-(omega2/2))*tNS[0,i]*((3*(cyNS[0,i]-uTotY2))+(9*(cuNS2*cyNS[0,i]))),(rho2-rho1)*grav)#Gravedad Fluido 2
           gOut[i,:,:] = gIn[i,:,:] - omega2*(gIn[i,:,:]-gEq[i,:,:])
           gOut[i,:,:] = gOut[i,:,:] + dt*fg2[i,:,:]

        if(len(bbRegiona)>0):
    ##      OBSTACLE (BOUNCE-BACK) Including Moving Walls for Shear Flow
            for i in xrange(9):
                cuNS1 = 3*(cxNS[0,i]*uTotX1+cyNS[0,i]*uTotY1)
                cuNS2 = 3*(cxNS[0,i]*uTotX2+cyNS[0,i]*uTotY2)
                fOut[i,bbRegiona[0:2*(lx+1):2],bbRegionb[0:2*(lx+1):2]] = fIn[oppNS[0,i],bbRegiona[0:2*(lx+1):2],bbRegionb[0:2*(lx+1):2]] + (2*rho1[bbRegiona[0:2*(lx+1):2],bbRegionb[0:2*(lx+1):2]]*tNS[0,i]*(UW)*cxNS[0,i])
                fOut[i,bbRegiona[1:2*(lx+1):2],bbRegionb[1:2*(lx+1):2]] = fIn[oppNS[0,i],bbRegiona[1:2*(lx+1):2],bbRegionb[1:2*(lx+1):2]] - (2*rho1[bbRegiona[1:2*(lx+1):2],bbRegionb[1:2*(lx+1):2]]*tNS[0,i]*(UW)*cxNS[0,i])
                gOut[i,bbRegiona[0:2*(lx+1):2],bbRegionb[0:2*(lx+1):2]] = gIn[oppNS[0,i],bbRegiona[0:2*(lx+1):2],bbRegionb[0:2*(lx+1):2]] + (2*rho2[bbRegiona[0:2*(lx+1):2],bbRegionb[0:2*(lx+1):2]]*tNS[0,i]*(UW)*cxNS[0,i])
                gOut[i,bbRegiona[1:2*(lx+1):2],bbRegionb[1:2*(lx+1):2]] = gIn[oppNS[0,i],bbRegiona[1:2*(lx+1):2],bbRegionb[1:2*(lx+1):2]] - (2*rho2[bbRegiona[1:2*(lx+1):2],bbRegionb[1:2*(lx+1):2]]*tNS[0,i]*(UW)*cxNS[0,i])

            
    ## STREAMING STEP FLUID 1 AND 2
        for i in xrange(9):
            fIn[i,:,:] = circ.circshift(fOut[i,:,:],cxNS[0,i],cyNS[0,i])
            gIn[i,:,:] = circ.circshift(gOut[i,:,:],cxNS[0,i],cyNS[0,i])

        t=t+dt

        print np.min(rho1)
        print np.max(rho1)
        
    ## VISUALIZATION
        if(cycle%itPlot==0):
            print 'iteracion:', cycle
            toc = time.clock()
            print toc-tic
            #Calculates Pressure where P=Cs^2(rho1+rho2)+1/Cs^2*G*rho1*rho2  and Cs^2=1/3
            Pin = ((1/3)*np.max(np.max(rho2))+ np.min(np.min(rho2)) + 3*G*np.max(np.max(rho2))*np.min(np.min(rho2)))
            Pout =((1/3)*np.max(np.max(rho1))+ np.min(np.min(rho1)) + 3*G*np.max(np.max(rho1))*np.min(np.min(rho1)))
            deltaP = Pout-Pin
            print deltaP
            curv = 1/coordinates[1,2]
            u = np.zeros([uTotX1.shape[0],uTotX1.shape[1]])
            u = sp.sqrt(np.multiply(uTotX1,uTotX1)+np.multiply(uTotY1,uTotY1)).reshape(lx+1,ly+1)
            #u[bbRegiona,bbRegionb] = sp.nan;
            name = str(cycle)
            name = nombre + name
            name2 = 'vel' + name
            plt.figure()
            plt.imshow(rho1.T, cmap=None, norm=None, aspect=None, interpolation=None,alpha=1.0, vmin=None, vmax=None, origin=None, extent=None)
            plt.colorbar()
            plt.savefig(name+'tiempo'+str(t)+'.png', dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format='png',transparent=False)
            plt.figure()
            plt.imshow(u.T, cmap=None, norm=None, aspect=None, interpolation=None,alpha=1.0, vmin=None, vmax=None, origin=None, extent=None)
            plt.colorbar()
            plt.savefig(name2+'.png', dpi=None, facecolor='w', edgecolor='w',orientation='portrait', papertype=None, format='png',transparent=False)
##            u = np.zeros([uTotX1.shape[0]*uTotX1.shape[1],3])
##            u[:,0]=uTotX1.reshape((lx+1)*(ly+1))
##            u[:,1]=uTotY1.reshape((lx+1)*(ly+1))
##            u[:,2]=0
##            pts = np.zeros([(lx+1)*(ly+1),3])
##            pts[:,0]=x.reshape((lx+1)*(lx+1),1).T
##            pts[:,1]=y.reshape((lx+1)*(lx+1),1).T
##            Den = pyvtk.StructuredGrid(dimensions=(lx+1,ly+1,1),points=pts)
##            Den.CellData = pyvtk.CellData(pyvtk.Scalars(rho1.T.ravel(),name='Densidad',lookup_table='default'))
##            #d = pyvtk.VtkData(Den,'Datos Densidad',Den.CellData)
##            Den.CellData = (pyvtk.CellData(pyvtk.Vectors(u,name='Velocidad')))
##            d = pyvtk.VtkData(Den,'Datos Velocidad',Den.CellData)
##            d.tofile(name+'.vtk', format='ascii')
