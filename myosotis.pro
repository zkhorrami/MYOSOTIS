; Make Your Own Synthetic ObservaTIonS (MYOSOTIS)
; !!! IDL version 8.5.1 !!!
; This code creates synthetic observations both imaging and spectroscopic data.
; The stellar and interstellar medium information (position, mass, velosity, metallicity,age, cloud density)
; should provided by the user.
; The user can choose different filters from the provided list of filters, depend of the observational instrument.
; Then the code will create the image and/or spectra within any FoV.
; this code is able to simulate ground-based telescopes data in addition to the space telescopes.
; The user can change the seeing and Strehl-Ratio values.
; The user can change the resolution of the image since most of the instruments do not reach their optimum resolution (~Lambda/Diameter).
; The user can choose her/his own pixel-scale for the given detector.
; The extinction can be applied on the output data, knowing the column density of the gas infront of each source.
;
;Inputs:
;   Project_name: any name chosen by the user
;   filestar: name of the file contains stellar sources information: X[pc],Y[pc],Z[pc],Vx[km/s],Vy[km/s],Vz[km/s],Mass[Mo],logage[yr],Z,cloud_density[Mo/pc^2]
;   filecloud: name of the file contains the cloud's information: X[pc],Y[pc],Z[pc],Vx[km/s],Vy[km/s],Vz[km/s],Mass[Mo],smoothing-length[pc]
;   filter: name of the instrument's filter which contains: Lambda[A],filter_trasparency
;   res: Pixel sampling of the instrument [arcsec/pix]
;   fovx,fovy: Field of view in x and y [arcsec]
;   fwhm: angular resolution of the observational instrument [arcsec]. The optimum resolution for a single dish telescope is about Lambda/diametr
;   PSFtype:type of PSF for distribution of the stellar fluxes on the detector: 'gassian' or 'airy'
;   distance: distance of the center of the simulation [X=0,Y=0,Z=0] from the observer [pc]
;   EXTmodel = 'Dmodel'    ; Extinction model, choose either 'Dmodel' OR 'Fmodel' (see the manual for more information on these models)
;   Rv: Av/E(B-V) constant value (e.g. 3.1 for Galaxy). Rv should be 3.1 or 4.0 or 5.5 if EXTmodel='Dmodel', otherwise can be set to any value.
;   
;   spectroscopy: 'yes' providing a cube data in wavelength range of lminspec[A]-lmaxspec[A] with the spectral resolution of Rspec
;                 'no' does not provide any spectroscopic data
;   velocitydis : 'yes' or 'no' : adding Doppler shift on the spectra according to the velocity of the stars
;   OBtreatment: 'yes' uses TLUSTY SEDs for O and B type stars
;   Adaptiveoptics: 'yes' the flux of the star will distribute partially in the airy pattern and the halo
;   SR: Strehl ratio, value between 0.0-1.0. if SR=1.0 the images are as perfect as space telescope's images
;   seeing: atmospheric resolution (fwhm of the halo) [arcsec]
;   Euler angles for rotation [degree]: 
;       alphai : rotation around x [degree]
;       bettai : rotation around y [degree]
;       gammai : rotation around z [degree]
;       For example [0,0,0] --> X-Y
;                   [90,0,0] --> X-Z
;                   [0,90,0] --> Y-Z
;   SNR : Signal to Noise ratio for the faintest star in the FOV. This will add sky noise to the image. MYOSOTIS will create two images, with and without noise.

; Output:
;   *_image.fits : 2D fits image
;   *_imageNoise.fits : Same as original image (*_image.fits) but with sky noise
;   *_star_info.txt : contains the information of the stellar sources in the FoV (11 columns) and also initial condition which user had used to create the image.
;                     These 11 columns are:
;                     mass[mo], logage[yr], Z, log[Teff[k]], logg, logL/Lo, Av_star, mag, Xpix, Ypix, assignedSED
;   if spectroscopy='yes' the two other outputs:
;       *_cube_spectra.fits : 3D cube, X-Y is the position of stellar sources, z is flux in different wavelengths
;       *_Lambda.txt : the wavelengths [A] which is used for the 3rd dimension of the spectral_cube
;
; Calling sequence:
;   IDL> .r myosotis
;   IDL> myosotis
;   
;  written by Z. Khorrami (KhorramiZ@cardiff.ac.uk)
;  updated Jan. 2019

pro MYOSOTIS

  project_name = 'test2'    ; Name of the project
  filestar = 'Examples/Teststar1.txt' ; Star's information file (should be 10 columns)
 ; filestar = 'mystars.txt' ; Star's information file (should be 10 columns)
  filecloud= 'Examples/NoCloud';'inputmyso/homo-M2000.00-R2.00000' ;'NoCloud'      ; Cloud's information file (should be 8 columns). Will be needed if Columndensities = 'sph'
  Columndensities= 'sph'   ; 'sph': reads the cloud information file, 'user': reads it from the filestar, last column in the unit of [Msun/pc2]

  filter = 'Filters/hst/wfpc2/HST-WFPC2.f814w.dat' ;choose the path of your favourite filter!
  res  = 0.05          ; pixel resolution [arcsec/pix]
  fovx = 35             ; Field-of-view X  [arcsec]
  fovy = 35             ; Field-of-view  Y [arcsec]
  fwhm = 0.11           ; FWHM of the PSF, angular resolution of the instrument
  PSFtype= 'gaussian'   ;'airy' Type of PSF for distribution of the stellar fluxes on the detector: 'gassian' or 'airy'
  
  distance = 50000.      ; Distance of the observer from the reference in filestar and filecloud [pc]
  EXTmodel = 'Dmodel'    ; Extinction model, choose either 'Dmodel' OR 'Fmodel' (see the manual for more information on these models)
  Rv       = 3.1         ; Rv should be 3.1 or 4.0 or 5.5 if EXTmodel='Dmodel', otherwise can be set to any value.
  
  metallicityZ=0.5         ; Z=1 for Galactic or Z=0.5 for LMC  will affect the evolutionary and atmosphere models 
  OBtreatment    = 'yes'  ; Uses TLUSTY SEDs for stars with Teff > 15000 K
  Adaptiveoptics = 'no'   ; Using of adaptive optics
  SR = 1.0                ; Strehl Ratio, should be between 0.0 - 1.0
  seeing = 0.8            ; Seeing [arcsec]

  spectroscopy = 'no'     ; Spectroscopy output, choose 'yes' or 'no'
  lminspec     = 4000.    ; Minimum wavelength [A] should be set within your filter transparency
  lmaxspec     = 6000.    ; Maximum wavelength [A] 
  Rspec        = 600      ; Spectral resolution (please check your SED library, should not be larger than the resolution of the SEDs)

  velocitydis  = 'no'     ; 'yes' or 'no' : Adding Doppler shift on the spectra according to the velocity of the stars 

  ;Euler angles for rotation [degree]:
  ;if all of them are zero, then outout image is X-Y
  ;[0,0,0] --> X-Y
  ;[90,0,0] --> X-Z
  ;[0,90,0] --> Y-Z
  alphai = 90   ; rotation around x [degree]
  bettai = 30   ; rotation around y [degree]
  gammai = 0   ; rotation around z [degree]

  SNR = 0.0;0.1    ; Signal to Noise ratio for the faintest star in the FOV
  noise2add=6.15883e-12;1.55350e-22;6.15883e-12 ; noise in the unit of flux/pix2=erg/cm2/s/A/pix2
  t0=SYSTIME(/SECONDS)
  FILE_MKDIR, project_name
  !PATH=!PATH+':'+Expand_Path('+./')
  myso_logo,'logo'
  mainscene,project_name,filestar,filecloud,filter,distance,Rv,res,fovx,fovy,fwhm,spectroscopy,$
    lminspec,lmaxspec,Rspec,OBtreatment,Adaptiveoptics,SR,seeing,alphai,bettai,gammai,velocitydis,SNR,EXTmodel,Columndensities,PSFtype,noise2add,metallicityZ
  t1=SYSTIME(/SECONDS)
  print,'Simulation time using 1 core: ',(t1-t0)/60.,'[min]'

end