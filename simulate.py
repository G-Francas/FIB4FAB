#Allows notebook to access other files within drive
from google.colab import drive
drive.mount('/content/drive') #This will require you enter passcode to access your drive
import sys

#Path to import extra modules from may need changing to match the path within your drive
path='/content/drive/My Drive/FIB 4 FAB/RCWA code/rcwa/examples'

sys.path.insert(0,f'{path}')

#Additional modules that RCWA code requires
import context
import matplotlib.pyplot as plt
import unittest
from shorthand import *
from matrices import *
from fresnel import *
from matrixParser import *
from source import Source
from layer import Layer
from solver import *
from crystal import Crystal
import cmath
import json
#Modules I've created
from slopecsv import Slopecsv
from checkcsv import *
from profileMkr import makeProf
from plot import plotter
from squarecsv import *
from Blazed import Blazedcsv
from polarization import *
def sim(gType,depthF,slope,profile,wavelength,period,har,res,material,theta,phi,pTEM, loop):
    tests=np.size(locals()[loop])
    loops=locals()[loop]
    
    depthFs=depthF
    slopes=slope
    wavelengths=wavelength
    periods=period
    hars=har
    ress=res
    
    depthF=depthFs[0]
    slope=slopes[0]
    wavelength=wavelengths[0]
    period=periods[0]
    har=hars[0]
    res=ress[0]
    
    

    print(f"\nSimulating a {gType} diffraction grating made of {material}....")
    ##refractive index vals at 400nm
    Si_n = 5.5674
    #Si_n=6.4732
    #Si_n=5.4754
    #Si_n=4.9760
    #Si_n=4.6784
    Al_n = 0.48787
    Ag_n = 0.05
    Au_n=1.4684
    Ti_n=2.0913
    
    #extinction coefficient vals at 400nm
    Si_k = 0.38612
    #Si_k=1.0686
    #Si_k=3.0024
    #Si_k=4.1990
    #Si_k=0.14851
    Al_k= 4.8355
    Ag_k = 2.1035
    Au_k=1.9530
    Ti_k=2.9556
    
    er_si=(Si_n)**2-(Si_k)**2+2*Si_n*Si_k*1j
    if material == 'Silicon':
        er=er_si
    elif material == 'Aluminium':
        er=(Al_n)**2-(Al_k)**2+2*Al_n*Al_k*1j
    elif material == 'Silver':
        er=(Ag_n)**2-(Ag_k)**2+2*Ag_n*Ag_k*1j
    elif material == 'Gold':
        er=(Au_n)**2-(Au_k)**2+2*Au_n*Au_k*1j
    elif material == 'Titanium':
        er=(Ti_n)**2-(Ti_k)**2+2*Ti_n*Ti_k*1j
 
    print(f'Permittivity = {er}')
    
    #Creating basic layers
    reflectionLayer = Layer(er=1.0006, ur=1, L=100)
    baseLayer = Layer(er=er_si,ur=1, L=1)
    source = Source(wavelength=wavelength, theta=theta, phi=phi,pTEM=pTEM, layer=reflectionLayer)

    
    if gType == 'Checkerboard':
        print("Creating layers for checkerboard simulation...\n")
        nLayers=Check(res,er, material, gType)
    if gType == 'Rectangular' and slope == 0:
        print("Creating layers for rectangular simulation with vertical sidewalls...\n")
        nLayers = Slopecsv(res,slope,depthF, material, gType, er)
    if gType == 'Square':
       print("Creating layers for square simulation...\n")
       nLayers=Square(res,er, material, gType) 
    if gType == 'CoatedChecker':
       print("Creating layers for coated checkerboard simulation...\n")
       nLayers=CheckCoat(res,er,material,gType, er_si)
    if gType == 'CheckerError':
        err=input("What percentage error would you like on the checkerboard? Please enter as a decimal, eg 0.3 for 30% error.")
        print(f"Creating layers for checkerboard simulation with error {err}...\n")
        nLayers=CheckErr(res,er, material, gType, err)
    if gType=='CheckDiag':
        print(f"Creating layers for checkerboard simulation on diagonal...\n")
        nLayers=CheckDiag(res,er,material,gType)
    if gType=='CheckDiagErr':
        err=input("What percentage error would you like on the checkerboard? Please enter as a decimal, eg 0.3 for 30% reduction in square size.")
        print(f"Creating layers for checkerboard simulation with error {err}...\n")
        nLayers=CheckDiagErr(res,er, material, gType, err)
    if gType=='Circ':
        print(f"Creating layers for circular grating simulation...\n")
        nLayers=Circ(res,er,material,gType)
        
    #setting up arrays to hold results of zeroth and first diffraction orders from simulation
    firstxs=np.zeros(tests)
    firstys=np.zeros(tests)
    zeros=np.zeros(tests)
    diags=np.zeros(tests)
    
    
    
    count=0
    while count < tests:
        depthF=depthFs.take([count],mode='clip')[0]
        
        slope=slopes.take([count],mode='clip')[0]
        wavelength=wavelengths.take([count],mode='clip')[0]
        period=periods.take([count],mode='clip')[0]
        har=hars.take([count],mode='clip')[0]
        res=ress.take([count],mode='clip')[0]
        depth=depthF*wavelength
        #generate csv files
        if slope != 0:
            if gType=="Rectangular":
                print(f"Creating layers for rectangular grating with sloped sidewalls, slope={slope} degrees...\n")
                nLayers = Slopecsv(res,slope,depthF, material, gType, er)
            elif gType == "Blazed":
                print(f"Creating layers for blazed grating of blaze angle {slope} degrees...\n")
                nLayers = Blazedcsv(res,slope,depthF, material, gType, er)
        #set period in each direction
        t1, t2 = complexArray([period, 0, 0]), complexArray([0, period, 0])
        n=0
        layers=np.empty(nLayers, Layer)
        print(f"Performing simulation trial {count+1} of {tests}, with {loop}={loops[count]}...")

        print(f"This simulation requires {nLayers} layers... ")
        #slope=90-slope

        while n < nLayers:
            
            print(f"creating layer {n+1} of {nLayers}...")
        
            if slope != 0:
              eCellData = np.transpose(np.genfromtxt(f'{path}/csvs/{gType}_{material}_slope{slope}_{depthF}_layer{count}_{res}.csv', delimiter=',', dtype=complex))
            else:
                eCellData = np.transpose(np.genfromtxt(f'{path}/csvs/{gType}_{material}_slope0_layer{n}_{res}.csv', delimiter=',', dtype=complex))
                print(f'file used:{gType}_{material}_slope0_layer{n}_{res}.csv ')
            uCellData = np.transpose(np.genfromtxt(f'{path}/csvs/U{res}.csv', delimiter=','))
            
            
            #setting depth of grating
            if gType == "Blazed":
                depth=period*np.tan(slope)
                
            crystalThickness = depth/nLayers
            
            if gType=="CoatedChecker":
                crystalThickness = depth

            numberHarmonics = (har,har)
            #numberHarmonics=(har,har)
            
            #creating crystal (grating) and grating layer
            deviceCrystal = Crystal(eCellData, uCellData, t1, t2)
            layers[n] = Layer(crystal=deviceCrystal, L=crystalThickness, numberHarmonics=numberHarmonics)
            n=n+1
        

        print(f"There are {nLayers} Layers of thickness {crystalThickness}")
        layerStack=LayerStack(reflectionLayer,*layers,baseLayer)

        print("Running simulation...\n")
        #perform simulation
        solver = Solver(layerStack, source, numberHarmonics)
        solver.Solve()
        sResults=solver.results
        rx=np.reshape(solver.rx,(har,har))
        ry=np.reshape(solver.ry,(har,har))
        rz=np.reshape(solver.rz,(har,har))
        # np.savetxt(f"{path}/results/RX_{gType}_{loop}_tests({tests})_{wavelength}_{period}.csv",rx , delimiter=",")
        # np.savetxt(f"{path}/results/RY_{gType}_{loop}_tests({tests})_{wavelength}_{period}.csv",ry , delimiter=",")
        # np.savetxt(f"{path}/results/RZ_{gType}_{loop}_tests({tests})_{wavelength}_{period}.csv",rz , delimiter=",")
        RX=np.ones([3,3], dtype='complex')
        RY=np.ones([3,3], dtype='complex')
        RZ=np.ones([3,3], dtype='complex')

        i,j=int(har/2)-1,int(har/2)-1
        a,b=0,0
        while i>=int(har/2)-1 and i <=int(har/2)+1:
            j=int(har/2)-1
            b=0
            while j>=int(har/2)-1 and j <=int(har/2)+1:
                RX[a,b]=complex(rx[i,j])
                RY[a,b]=complex(ry[i,j])
                RZ[a,b]=complex(rz[i,j])
                j=j+1
                b=b+1
            i=i+1
            a=a+1
        
        
        
        stokes=transform(RX,RY,RZ,wavelength,period,'x',False)
        print("STOKES:\n")
        print(stokes)
        
        
        # Get the diffraction efficiencies R and T and overall reflection and transmission coefficients R and T
        (R, RTot) = (solver.R, solver.RTot)
        print(f"Total Reflection {RTot}")
        
        
        # if profile:
        #     makeProf(slope,depthF,nLayers,res, material)
        firstx=R[int(har/2)][int(har/2 - 1)]
        firsty=R[int(har/2-1)][int(har/2)]
        zeroth=R[int(har/2)][int(har/2)]
        diag=R[int(har/2+count)][int(har/2+count)]
        diags[count]=diag
        firstxs[count]=firstx
        firstys[count]=firsty
        zeros[count]=zeroth
        print("Plotting...")
        plotter(firstxs,firstys,zeros,diags,R,loop, gType, tests,count, profile, slope, depthF, nLayers, res, material, wavelength, period, har, loops[count])
        count=count+1
        
    plt.figure()
    plt.plot(loops,zeros,label="0th order efficiency")   #simulated 0th orders at varying depths
    plt.plot(loops,firstxs,label="1st order efficiency - x direction")  #simulated 1st orders at varying depths
    #plt.scatter(loops,firstys,label="1st order efficiency - y direction")
    plt.plot(loops,diags,label="diagonal efficiency modes")
    plt.ylabel("Intensity")
    plt.xlabel(f"{loop}")
    plt.legend(loc='upper left', bbox_to_anchor=(1.05, 1))
    plt.title(f"Efficiency of {gType} {material} grating at different {loop} ")
    tuple=(loops, zeros, firstxs, firstys, diags)
    results=np.vstack(tuple)
    np.savetxt(f"{path}/results/{gType}_{loop}_tests({tests})_{wavelength}_{period}.csv",results , delimiter=",")
    return 0
