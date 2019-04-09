#!/usr/bin/env python
import sys
from math import *
import random
import numpy as np
from scipy.linalg import expm
from astropy.io import fits
import os, sys, shutil
from astropy.convolution import AiryDisk2DKernel
import params

pi=np.pi
def myso_logo(wh):
    if (wh == 'logo'):        
        print ' ================================================================================= '
        print '|   ===============   Make   Your   Own Synthetic   ObservaTIonS  =============   |'
        print '| =============================================================================== |'
        print '||                                                                               ||'
        print '|| |%|       |%| |%|   |%| |%||%||%| |%||%||%| |%||%||%| |%||%||%| |%| |%||%||%| ||'
        print '|| |%||%| |%||%| |%|   |%| |%|   |%| |%|       |%|   |%|    |%|    |%| |%|       ||'
        print '|| |%|  |%|  |%|    |%|    |%|   |%| |%||%||%| |%|   |%|    |%|    |%| |%||%||%| ||'
        print '|| |%|  |%|  |%|    |%|    |%|   |%|       |%| |%|   |%|    |%|    |%|       |%| ||'
        print '|| |%|  |%|  |%|    |%|    |%||%||%| |%||%||%| |%||%||%|    |%|    |%| |%||%||%  ||'
        print '||                                                                               ||'
        print '| =============================================================================== |'
        print '||        MMMMMMMMMMMMMMMMMWN00KX0OOO0KNMMMMMMMMMMMMMMMMMMMM      ||'
        print '||        MMMMMMMMMMMMMMMMW0dloddolcccld0NMMMMMMMMMMMMMMMMMM      ||'
        print '||        MMMMMMMMMMMMMMMNklcloooooollccckWMMMMMMMMMMMMMMMMM      ||'
        print '||        MMMMMMMMMMMMMMMKo::clooooollccclOWMMMMMMMMMMMMMMMM      ||'
        print '||        MMMMMMMMMMMMMMWOl:ccclllllllccccoKMMMMMMMMMMMMMMMM      ||'
        print '||        MMMMMMMMMMMMMMNxccccccclllllcccclkNMMMMMWWNXNNWMMM      ||'
        print '||        MMMMMMMWWWWWWWNxcccc::::ccccc::ccxNMWNKOkxdoodxONM      ||'
        print '||        MMMWX0OkkkxxxkOxcc::::;;::cc:::clONX0kxdxxxdoolcxN      ||'
        print '||        MMXkllloooooooollloc::;;;:::::ldkOOkxkkkkkxxxddodK      ||'
        print '||        MNd::cclloodddoolcoxd:,;:c:;:dOxddddddxdddddddoooO      ||'
        print '||        Xxlccclllloooooll::dOkodkkxoxko:clooolccccccccc:cO      ||'
        print '||        kllllccccc:cccclllldkkxdllooxxollcc::::cccccccclkN      ||'
        print '||        klllcccc:::::::::cokxdc,...llool;,;;::::::::cldKWM      ||'
        print '||        Xklcccc:::;;;;;;;;ldool;...llolddol:,,,;;::lxKWMMM      ||'
        print '||        MWX0xlc::;;;,,,;cdkxdxdolcloxkxlllc;,,:oxk0XWMMMMM      ||'
        print '||        MMMMWNK0Okxol;;::;;:lxolxkOxdkkoc;,,,,c0WMMMMMMMMM      ||'
        print '||        MMMMMMMMMMWNxc:::;:cllccx00klldddolc:::l0WMMMMMMMM      ||'
        print '||        MMMMMMMMMMWOlccccccllcccccllodddddolcccco0WMMMMMMM      ||'
        print '||        MMMMMMMMMM0occ::clccclcccclodddddddoollccoKMMMMMMM      ||'
        print '||        MMMMMMMMMWkccccccccccclolcdkxdxxxdddooolllxNMMMMMM      ||'
        print '||        MMMMMMMMMNxcc::ccccclloxO0XWX0Oxxddddoooolo0WMMMMM      ||'
        print '||        MMMMMMMMMNx::::ccclllokXWMMMMMWXKOxddooolllkWMMMMM      ||'
        print '||        MMMMMMMMMMXxl::cllclxKWMMMMMMMMMMWNK0OxdodONMMMMMM      ||'
        print '||        MMMMMMMMMMMWNOlccokKWMMMMMMMMMMMMMMMMMWNNWMMMMMMMM      ||'
        print '||        MMMMMMMMMMMMMW0xkKWMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM      ||'
        print ' =======================     FORGET ME NOT     ===================  '

    if (wh == 'star'): 
        print ' =============================================================== '
        print '|               READING YOUR STELLAR SOURCES FILE               |'
        print ' =============================================================== '

    if (wh == 'cloud'):
        print ' =============================================================== '
        print '|                   READING YOUR CLOUD  FILE                    |'
        print ' =============================================================== '

    if (wh  == 'filter'):
        print ' =============================================================== '
        print '|                         READING FILTER                        |'
        print ' =============================================================== '
 
    if (wh  == 'iso'):
        print ' =============================================================== '
        print '|                  READING EVOLUTIONARY MODELS                  |'
        print ' =============================================================== '

    if (wh  == 'spec'):
        print ' =============================================================== '
        print '|                  SPECTRA CHANNEL is CREATED                   |'
        print ' =============================================================== '
 
    if (wh  == 'tlusty'):
        print ' =============================================================== '
        print '|                USING TLUSTY for MASSIVE STARs                 |'
        print ' =============================================================== '

    if (wh  == 'ao'):
        print ' =============================================================== '
        print '|              ADOPTIVE OPTICS is being APPLIED                 |'
        print ' =============================================================== '
 
    if (wh  == 'outim'):
        print ' =============================================================== '
        print '|            OUTPUT SYNTHETIC IMAGE is WRITTEN in :             |'
        print ' =============================================================== '
 
    if (wh  == 'outspec'):
        print ' =============================================================== '
        print '|          OUTPUT SYNTHETIC SPECTRA CUBE is WRITTEN in :        |'
        print ' =============================================================== '

def FINDCLOSE(par1,pararr):
    nclose=np.argsort(abs(np.add(-par1,pararr)))
    return nclose[0]

def rot_euler(v, xyz):
    ''' Rotate vector v (or array of vectors) by the euler angles xyz '''
    # https://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector
    for theta, axis in zip(xyz, np.eye(3)):
        v = np.dot(np.array(v), expm(np.cross(np.eye(3), axis*-theta)))
    return v

def Kernel_fn(dr, h):
    pi=np.pi
    if (dr/h < 1.0):
        w = (1/(pi*h**3.)) * (1-1.5*(dr/h)**2. + 0.75*(dr/h)**3.)
    elif (dr/h < 2.0):
        w = (1/(pi*h**3.)) * (0.25 * (2.0-(dr/h))**3.)
    else: w = 0.0

    return w
  
def COLUMN_DENSITY(xstar,ystar,zstar,xarr,yarr,zarr,marr,harr):
    nc=len(marr)
    den2dtot=0.0
    for ii in range(nc):
        if ((sqrt((xstar-xarr[ii])**2+(ystar-yarr[ii])**2) < 2.*harr[ii]) and (zstar < zarr[ii])):
            dr=sqrt((xstar-xarr[ii])**2+(ystar-yarr[ii])**2)
            lc=2*sqrt(4.*harr[ii]*harr[ii]-dr*dr)
            myres=20
            deltalc=lc/myres
            denarr=np.zeros(myres)
            lcarr=np.zeros(myres)
            for jj in range(myres): denarr[jj]= marr[ii]*Kernel_fn(sqrt(dr*dr+(jj*lc/myres-lc/2.)**2.),harr[ii])*deltalc
            den2dtot = den2dtot+ np.sum(denarr)
    return den2dtot

def extinctions(lam,Av,Rv):
    Alam=0.0  
    xarr=1./lam
    if ((xarr >= 0.3) and (xarr < 1.1)):
      ax1=0.574*(xarr**1.61)
      bx1=-0.527*(xarr**1.61)
      Alam=Av*(ax1+bx1/Rv)
    elif ((xarr >= 1.1) and (xarr < 3.3)):
      y = xarr-1.82
      ax2 = 1.0+0.17699*y-0.50447*(y**2.0)-0.02427*(y**3.0)+0.72085*(y**4.0) +0.01979*(y**5.0)-0.77530*(y**6.0)+0.32999*(y**7.0)
      bx2 = 1.41338*y+2.28305*(y**2.0)+1.07233*(y**3.0)-5.38434*(y**4.0)-0.62251*(y**5.0)+5.30260*(y**6.0) -2.09002*(y**7.0)
      Alam=Av*(ax2+bx2/Rv)
    elif ((xarr >= 3.3) and (xarr < 5.9)):
      ax3 = 1.752 -0.316*xarr-0.104/((xarr-4.67)**2.0 + 0.341)
      bx3 = (-3.090)+1.825*xarr+1.206/((xarr-4.62)**2.0+0.263)
      Alam=Av*(ax3+bx3/Rv)
    elif ((xarr >= 5.9) and (xarr < 8.0)):
      fax4 = -0.04473*(xarr - 5.9)**2.0 - 0.009779*(xarr -5.9)**3.0
      fbx4 = 0.2130*(xarr-5.9)**2.0 + 0.1207*(xarr-5.9)**3.0
      ax4 = 1.752 - 0.316*xarr - (0.104/((xarr - 4.67)**2.0 + 0.341)) +fax4
      bx4 = -3.090 + 1.825*xarr + (1.206/((xarr - 4.62)**2.0 + 0.263)) +fbx4
      Alam=Av*(ax4+bx4/Rv)
    else:  print 'Lambda,',lam,', is out of range for Fitzpatrick (1999) model!!!'
    return Alam

def BCcal(wavelength,flux,lambdaF,weight,AVstar,Rv,Teff,par2vega,EXTmodel,DrainearrLam,DrainearrK): # sed, filter
#  BC= Mbol_sun -2.5 log10 [4!PI 10pc^2 Sigma Teff^4 / Lsun] + 2.5 .... Girardi + 2002
    Mbolsun = 4.77
    bolconstant = -2.5*log10(4*pi*(3.0857**2.0)*5.67051/3.826)
 
    n_wave=len(wavelength)
# I put the zero weight at the edges just to solve the problem of edges in linear interpolation
    weight[0]=0.0
    weight[len(weight)-1]=0.0
    wavelength_weight=np.interp(wavelength,lambdaF,weight)
    par1=0.0
    par2=0.0
    alamarr=np.zeros(len(wavelength)) #fltarr((size(wavelength))[-1])

    if (EXTmodel == 'Dmodel'):
     kappased=np.interp(wavelength,DrainearrLam,DrainearrK)
     alamarr=AVstar*kappased

  
    for i in range(1, len(wavelength)):
      ltemp=wavelength[i]*1.0e-4 
      if (EXTmodel == 'Fmodel' and wavelength[i] >= 1250. and wavelength[i] <= 33333.):
          alamarr[i]=extinctions(ltemp,AVstar,Rv)
      par1 += wavelength[i]*flux[i]*(10.0**(-0.4*alamarr[i]))*wavelength_weight[i]*(wavelength[i]-wavelength[i-1])
#     par2 += wavelength[i]*wavelength_weight[i]*(wavelength[i]-wavelength[i-1])  ;!!!! par2 will be calculating from Vega flux
    BCfilter = Mbolsun + bolconstant - 10.*log10(Teff)+2.5*log10(par1/par2vega)
    return BCfilter

def BCcals(wavelength,flux,lambdaF,weight,AVstar,Rv,Teff,EXTmodel,DrainearrLam,DrainearrK): # sed, filter
#  BC= Mbol_sun -2.5 log10 [4!PI 10pc^2 Sigma Teff^4 / Lsun] + 2.5 .... Girardi + 2002
    Mbolsun = 4.77
    bolconstant = -2.5*log10(4*pi*(3.0857**2.0)*5.67051/3.826)
 
    n_wave=len(wavelength)
# I put the zero weight at the edges just to solve the problem of edges in linear interpolation
    weight[0]=0.0
    weight[len(weight)-1]=0.0
    wavelength_weight=np.interp(wavelength,lambdaF,weight)
    par1=0.0
    par2=0.0
    alamarr=np.zeros(len(wavelength)) #fltarr((size(wavelength))[-1])

    if (EXTmodel == 'Dmodel'):
     kappased=np.interp(wavelength,DrainearrLam,DrainearrK)
     alamarr=AVstar*kappased

  
    for i in range(1, len(wavelength)):
      ltemp=wavelength[i]*1.0e-4 
      if (EXTmodel == 'Fmodel' and wavelength[i] >= 1250. and wavelength[i] <= 33333.):
          alamarr[i]=extinctions(ltemp,AVstar,Rv)
      par1 += wavelength[i]*flux[i]*(10.0**(-0.4*alamarr[i]))*wavelength_weight[i]*(wavelength[i]-wavelength[i-1])
      par2 += wavelength[i]*wavelength_weight[i]*(wavelength[i]-wavelength[i-1])  #!!!! par2 will be calculating from Vega flux
    BCfilter = Mbolsun + bolconstant - 10.*log10(Teff)+2.5*log10(par1/par2)
    return BCfilter

def makeGaussian(size, fwhm, center):
    """ Make a square gaussian kernel.
    size is the length-array of a side of the square
    fwhm is full-width-half-maximum, which
    can be thought of as an effective radius.
    """
    x = np.arange(0, size[0], 1, float)
    y = np.arange(0, size[1], 1, float)
    y = y[:,np.newaxis]
    if center is None:
        x0 = y0 = size // 2
    else:
        x0 = center[0]
        y0 = center[1]
    return np.exp(-4*np.log(2) * ((x-x0)**2 + (y-y0)**2) / fwhm**2)

project_name=params.project_name
EXTmodel=params.EXTmodel
filterfile=params.filterfile
filestar=params.filestar
filecloud=params.filecloud
Columndensities=params.Columndensities
OBtreatment=params.OBtreatment
alphai=params.alphai
bettai=params.bettai
gammai=params.gammai
res=params.res
fovx=params.fovx
fovy=params.fovy
fwhm=params.fwhm
SNR=params.SNR
Rv=params.Rv
distance=params.distance

spectroscopy = params.spectroscopy     # Spectroscopy output, choose 'yes' or 'no'
lminspec     = params.lminspec    # Minimum wavelength [A] should be set within your filter transparency
lmaxspec     = params.lmaxspec   # Maximum wavelength [A] 
Rspec        = params.Rspec #700      # Spectral resolution (please check your SED library, should not be larger than the resolution of the SEDs)

velocitydis  = params.velocitydis

PSFtype=params.PSFtype
Adaptiveoptics=params.Adaptiveoptics
seeing=params.seeing
SR=params.SR
vegafile='vegaf'
evolutionary='Evolutionary/Z0p015.dat'
foldersed='SEDs/'
Mbolsun = 4.77

outputim     = project_name+'/'+project_name+'_image.fits' 
outputimnoise=project_name+'/'+project_name+'_imageNoise.fits'
outputspecFL = project_name+'/'+project_name+'_cube_spectra.fits'
outputspecL  = project_name+'/'+project_name+'_Lambda.txt'
outputstarinfo  = project_name+'/'+project_name+'_star_info.txt'


if not os.path.exists(project_name):
    os.makedirs(project_name)
else:
    shutil.rmtree(project_name)           
    os.makedirs(project_name)

myso_logo('logo')

lunstarinfo=open(outputstarinfo,"w")

lunstarinfo.write("#     mass[mo] ,   logage[yr]  ,     Z      , log[Teff[k]] ,    logg    ,   logL/Lo    ,  Av_star   ,      mag     ,    Xpix     ,   Ypix       ,       assignedSED \n")

lunstarinfo.write("#       1              2             3            4              5              6           7               8            9            10                   11 \n")


lambdaF,weight=np.loadtxt(filterfile,unpack=True)

DrainearrLam, DrainearrK=[0.0],[0.0]
if (EXTmodel == 'Dmodel'):
    if (Rv == 3.1): 
        Drainemodel='source/Draine3p1.txt'
        DraineKappaV=8.551E+03
    elif (Rv == 4.0): 
        Drainemodel='source/Draine4.txt'
        DraineKappaV=8.492E+03
    elif (Rv == 5.5): 
        Drainemodel='source/Draine5p5.txt'
        DraineKappaV=7.313E+03
    else: print 'For Dmodel, R_V should be 3.1 or 4.0 or 5.5. If you need other Rv values please choose Fmodel'
    DrainearrLamu,drainealbedo,drainecos,draineC,DrainearrKu,drainecos2=np.loadtxt(Drainemodel,usecols=(0,1,2,3,4,5),unpack=True)
    DrainearrLamu=DrainearrLamu*1.0E+4
    DrainearrKu=DrainearrKu/DraineKappaV
    DrainearrLam, DrainearrK = zip(*sorted(zip(DrainearrLamu,DrainearrKu)))

xpix=round(fovx/res)
ypix=round(fovy/res)

print 'FoV: ',fovx,'" x ',fovy,'" =',int(xpix),'[pix] x',int(ypix),'[pix]'

sceneim=np.full((int(ypix),int(xpix)),0.0)

wavelengthvega,fluxvega=np.loadtxt(foldersed+vegafile,unpack=True)
wavelength_weightvega=np.interp(wavelengthvega,lambdaF,weight)
par2vega=0.0
for i in range(1, (len(wavelengthvega))):
  par2vega += fluxvega[i]*wavelengthvega[i]*wavelength_weightvega[i]*(wavelengthvega[i]-wavelengthvega[i-1])


Lspecarr=[]
if (spectroscopy == 'yes'):
    myso_logo('spec')
    Lspecarr.append(lminspec)
    while (Lspecarr[-1] <= lmaxspec): Lspecarr.append(Lspecarr[-1]+Lspecarr[-1]/Rspec)
    nspec=len(Lspecarr)
    sceneL=np.zeros(nspec)
    sceneimFL=np.full((nspec,int(ypix),int(xpix)),0.0)
    lun=open(outputspecL,"w")
    lun.write('# Lambda[A] \n')
    for ll in range(nspec): lun.write("%f\n" %(Lspecarr[ll]))  
    lun.close()
    
xstar,ystar,zstar,vxstar,vystar,vzstar,massstar,logagestar,kzstar,rhostar=np.loadtxt(filestar,unpack=True)
nstar=len(xstar)
positionvector=zip(xstar,ystar,zstar)
velocityvector=zip(vxstar,vystar,vzstar)

zeinab=rot_euler(positionvector,np.multiply([alphai,bettai,gammai],pi/180.))
xstar=zeinab[0:nstar,0]
ystar=zeinab[0:nstar,1]
zstar=zeinab[0:nstar,2]

zeinabv=rot_euler(velocityvector,np.multiply([alphai,bettai,gammai],pi/180.))
vxstar=zeinabv[0:nstar,0]
vystar=zeinabv[0:nstar,1]
vzstar=zeinabv[0:nstar,2]

distancestar=np.add(distance,-zstar)
pc2pixstar=206264.806247/distancestar/res
Teffstar=np.zeros(nstar)
loggstar=np.zeros(nstar)
loglstar=np.zeros(nstar)
sedstar=np.zeros(nstar,dtype=np.uint64)
fluxstar=np.zeros(nstar)
newx=xstar*pc2pixstar #convert x[pc] into pixel position
newy=ystar*pc2pixstar
newz=zstar*pc2pixstar
columncloud=np.zeros(nstar) #column density of the cloud in fron of each star
print 'READ STARs: ',nstar

if (Columndensities == 'sph'):
    myso_logo('cloud')
    xcloud,ycloud,zcloud,vxcloud,vycloud,vzcloud,masspar,hpar=np.loadtxt(filecloud,unpack=True)
    ncloud=len(xcloud)
    
    positioncvector=zip(xcloud,ycloud,zcloud)
    zeinabc=rot_euler(positioncvector,np.multiply([alphai,bettai,gammai],pi/180.))
    xcloud=zeinabc[0:nstar,0]
    ycloud=zeinabc[0:nstar,1]
    zcloud=zeinabc[0:nstar,2]

    distancecloud=np.add(distance,-zcloud)
    pc2pixcloud=206264.806247/distancecloud/res
    newxcloud=xcloud*pc2pixcloud #convert x[pc] into pixel position
    newycloud=ycloud*pc2pixcloud
    newzcloud=zcloud*pc2pixcloud
    newhcloud=hpar*pc2pixcloud

#  reading the isochrones
myso_logo('iso')
ziso,logageiso,miniiso,mactiso,logliso,logteff,loggiso=np.loadtxt(evolutionary,usecols=(0,1,2,3,4,5,6),unpack=True)
teffiso=10.**logteff
niso=len(miniiso)


#  ;reading list of SEDs:
teffsed,loggsed,metallicity,lh,vtur,sedname=np.loadtxt(foldersed+'kseds.dat',dtype={'names': ('col1', 'col2', 'col3','col4', 'col5', 'col6'), 'formats':(np.float,np.float,np.float,np.float,np.float,'|S60')},usecols=(0,1,2,3,4,5),unpack=True)
nseds=len(teffsed)

if (OBtreatment == 'yes'):
    teffsedOB,loggsedOB,metallicityOB,sednameOB=np.loadtxt(foldersed+'Tseds.dat',dtype={'names': ('col1', 'col2', 'col3','col4'), 'formats':(np.float,np.float,np.float,'|S60')},usecols=(0,1,2,3),unpack=True)
    nsedsOB=len(teffsedOB)
    myso_logo('tlusty')


#wavelength,flux=np.loadtxt(foldersed+sedname[10],comments=['fn:', '#'],unpack=True)

nstars=0
for ii in range(nstar):
    if ((abs(newx[ii]) < (xpix/2)-1) and (abs(newy[ii]) < (ypix/2)-1)):
        nstars += 1
        if (Columndensities == 'sph'):
            columncloud[ii]=COLUMN_DENSITY(newx[ii],newy[ii],newz[ii],newxcloud,newycloud,newzcloud,masspar,newhcloud)
            cloudden=columncloud[ii]*((1.989/1.673534)/(((3.0857/pc2pixstar[ii])**2.))) #*1.0e21 ;convert column density unit [Msun/pix^2] --> [10^21 * Mhydrogen/cm^2]
        else: cloudden=rhostar[ii]*((1.989/1.673534)/((3.0857**2.)))#*1.0e21 ;convert column density unit [Msun/pc^2] --> [10^21 * Mhydrogen/cm^2]
        
#        cloudden=rhostar[ii]*((1.989/1.673534)/((3.0857**2))) #*1.0e21 #convert column density unit [Msun/pc^2] --> [10^21 * Mhydrogen/cm^2]
        AVstar=cloudden/2.21 #e21 ;Guver&Ozel2009: The relation between Optical Extinction and Hydrogen column density
        nage=FINDCLOSE(logagestar[ii],logageiso)
        selectedage=logageiso[nage]
        nmetalicity=FINDCLOSE(kzstar[ii],ziso)
        selectedz=ziso[nmetalicity]
        marrtemp=np.full(niso,999.99)
        for kk in range(niso):  
            if ((ziso[kk] == selectedz) and (logageiso[kk] == selectedage)):  marrtemp[kk]=mactiso[kk]
        ns=FINDCLOSE(massstar[ii],marrtemp)
        Teffstar[ii]=teffiso[ns]
        loggstar[ii]=loggiso[ns]
        loglstar[ii]=logliso[ns]

        deltaT=abs(Teffstar[ii]-teffsed)

        deltagarr=np.full(nseds,99.)


        for jj in range(nseds): 
            if (deltaT[jj] == min(deltaT)):  deltagarr[jj]=abs(loggstar[ii]-loggsed[jj])
#        sedstar[ii]=where(deltagarr eq min(deltagarr))
        sedstar[ii]=FINDCLOSE(min(deltagarr),deltagarr)
        readsed=sedname[sedstar[ii]]
#        readcol,foldersed+sedname[sedstar[ii]],wavelength,flux,/silent
        wavelength,flux=np.loadtxt(foldersed+sedname[sedstar[ii]],comments=['fn:', '#'],unpack=True)

        if ((OBtreatment == 'yes') and (Teffstar[ii] >= 15000.)):
            deltaT=abs(Teffstar[ii]-teffsedOB)
            deltagarr=np.full(nsedsOB,99.)
            for jj in range(nsedsOB):  
                if (deltaT[jj] == min(deltaT)):  deltagarr[jj]=abs(loggstar[ii]-loggsedOB[jj])
#            sedstar[ii]=where(deltagarr eq min(deltagarr))
            sedstar[ii]=FINDCLOSE(min(deltagarr),deltagarr)

            readsed=sednameOB[sedstar[ii]]
#            readcol,foldersed+sednameOB[sedstar[ii]],wavelength,flux,/silent
            wavelength,flux=np.loadtxt(foldersed+sednameOB[sedstar[ii]],comments=['fn:', '#'],unpack=True)

        bc1=BCcal(wavelength,flux,lambdaF,weight,AVstar,Rv,Teffstar[ii],par2vega,EXTmodel,DrainearrLam,DrainearrK)
        mag=Mbolsun-2.5*loglstar[ii]-(bc1)+5.0*log10(distancestar[ii]/10.0)
        fluxstar[ii]=10.**(mag/(-2.5))

        lunstarinfo.write("%13.4f %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f %13.4f %60s \n" %(massstar[ii], logagestar[ii],kzstar[ii],log10(Teffstar[ii]),loggstar[ii],loglstar[ii],AVstar,mag,newx[ii]+xpix/2.,newy[ii]+ypix/2.,readsed))

        if (PSFtype == 'gaussian'):
            airy1=makeGaussian((xpix,ypix),fwhm/res,center=(newx[ii]+xpix/2.,newy[ii]+ypix/2.))
        
        airy1=airy1/np.sum(airy1) #to be sure about normaized total flux across FOV 

        if (Adaptiveoptics == 'yes'):
            halo=makeGaussian((xpix,ypix),seeing/res,center=(newx[ii]+xpix/2.,newy[ii]+ypix/2.))
            halo=halo/np.sum(halo)
            sceneim += fluxstar[ii]*(SR*airy1+(1.0-SR)*halo)
        if (Adaptiveoptics == 'no'): sceneim += fluxstar[ii]*airy1

        if (spectroscopy == 'yes'):
            for ll in range(2,nspec-3): 
#                linterp,lambda,weight,wavelength,wavelength_weight
                wavelength_weight=np.interp(wavelength,lambdaF,weight)
                bc3=BCcals(wavelength,flux*wavelength_weight,[Lspecarr[ll-1],Lspecarr[ll],Lspecarr[ll+1],Lspecarr[ll+2]],[0.0,1.0,1.0,0.0],AVstar,Rv,Teffstar[ii],EXTmodel,DrainearrLam,DrainearrK)
                mag3=Mbolsun-2.5*loglstar[ii]-bc3+5.0*log10(distancestar[ii]/10.0)
                fluxstar3=float(10.**(mag3/(-2.5)))
                if (PSFtype == 'gaussian'):
                    airy3=makeGaussian((xpix,ypix),fwhm/res,center=(newx[ii]+xpix/2.,newy[ii]+ypix/2.))

                airy3=airy3/(np.sum(airy3))

                if (velocitydis == 'yes'):
                    shiftv=vzstar[ii]*Lspecarr[ll]/3.0E+5 #shift is in A
                    lambdashift=Lspecarr[ll]-shiftv        #if vz>0 ==> source comes toward observer ==> lambdashift<lamda0 (blue-shift)
                    llchannel=FINDCLOSE(lambdashift,Lspecarr)
                    if (Adaptiveoptics == 'yes'): sceneimFL[llchannel,:,:] += fluxstar3*(SR*airy3+ (1.0-SR)*halo)
                    if (Adaptiveoptics == 'no' ): sceneimFL[llchannel,:,:] += fluxstar3*airy3

                elif (velocitydis == 'no'):
                    if (Adaptiveoptics == 'yes'): sceneimFL[ll,:,:] += fluxstar3*(SR*airy3+ (1.0-SR)*halo)
                    if (Adaptiveoptics == 'no' ): sceneimFL[ll,:,:] += fluxstar3*airy3




print 'Number of stars in the FoV: ',nstars


hdu = fits.PrimaryHDU(sceneim)
hdu.writeto(outputim)

faintestflux=min(i for i in fluxstar if i > 0)
noise=faintestflux*4.*log(2.0)/(pi*(fwhm/res)*(fwhm/res))/SNR
noise2add=noise*np.random.rand(int(ypix),int(xpix))

hdu = fits.PrimaryHDU(sceneim+noise2add)
hdu.writeto(outputimnoise)

if (spectroscopy == 'yes'):
    hdu = fits.PrimaryHDU(sceneimFL)
    hdu.writeto(outputspecFL)

lunstarinfo.write('#faintestflux: %f \n' %(faintestflux))
lunstarinfo.write('#noise:  %f \n' %(noise))
lunstarinfo.write('#SNR:  %f \n' %(SNR))
lunstarinfo.write('#Filter:  %s \n' %(filterfile))
lunstarinfo.write('#distance:  %f \n' %(distance))
lunstarinfo.write('#res:  %f \n' %(res))
lunstarinfo.write('#FoV:  %d %d \n' %(fovx,fovy))
lunstarinfo.write('#fwhm:  %f \n' %(fwhm))
lunstarinfo.write('#OBtreatment:  %s \n' %(OBtreatment))
lunstarinfo.write('#Adaptiveoptics:  %s \n' %(Adaptiveoptics))
lunstarinfo.write('#SR:  %f \n' %(SR))
lunstarinfo.write('#Seeing:  %f \n' %(seeing))
lunstarinfo.write('#Spectroscopic Resolution:  %f \n' %(Rspec))
lunstarinfo.write('#Lambda_min:  %f \n' %(lminspec))
lunstarinfo.write('#Lambda_max:  %f \n' %(lmaxspec))
lunstarinfo.write('#filestar:  %f \n' %(filestar))
lunstarinfo.write('#Columndensities:  %f \n' %(Columndensities))
lunstarinfo.write('#filecloud:  %f \n' %(filecloud))
lunstarinfo.write('#EXTmodel:  %f \n' %(EXTmodel))
lunstarinfo.write('#Rv:  %f \n' %(Rv))
lunstarinfo.close()

myso_logo('outim')
print outputim
print '   '
myso_logo('outspec')
print outputspecFL
print outputspecL
print  '   '
