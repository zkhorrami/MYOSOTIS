; This is a main part of the MYOSOTIS code which calculate stellar fluxes, extinction and 
; creates images and spectra.
; It reads the Evolutionary models and SEDs
; Sets up right extinction model according to the Rv value
; Creates output files and arrays and if spectroscopy is ON then sets up the spectroscopic channel
; Reads the "filestar" which should be 10 columns
; If the columndensities='sph' (means gas cloud is provided by the standard sph codes) then 
;   it reads the 'filecloud' which should have 8 columns
; Reads Vega SED (flux observed from Earth) to be used to calibrate bolometric corrections.
; For each star it finds the right stellar parametrs from the evolutionary models and selects a proper SED
; from the SED library, then calculates the stellar flux within a given filter using its extinction along the sight of observer
; 
;
;
;  written by Z. Khorrami (KhorramiZ@cardiff.ac.uk)
;  updated Jan. 2018

pro mainscene,project_name,filestar,filecloud,filter,distance,Rv,res,fovx,fovy,fwhm,spectroscopy,$
  lminspec,lmaxspec,Rspec,OBtreatment,Adaptiveoptics,SR,seeing,alphai,bettai,gammai,velocitydis,SNR,EXTmodel,Columndensities,PSFtype,noise2add,metallicityZ
 
if (metallicityZ eq 1.0) then begin
  evolutionary = 'Evolutionary/Z0p015.dat'
  foldersed = 'SEDs/Z1/'
endif else if (metallicityZ eq 0.5) then begin
  evolutionary = 'Evolutionary/Z0p008.dat'
  foldersed = 'SEDs/Z0p5/'
endif else begin
  print, '!!! metallicityZ should be 1.0 (for Solar) or 0.5 (for LMC)'
  stop
endelse

if (EXTmodel eq 'Dmodel') then begin

  CASE Rv OF
    3.1: begin
            Drainemodel='source/Draine3p1.txt'
            DraineKappaV=8.551E+03
         end
    4.0: begin
           Drainemodel='source/Draine4.txt'
           DraineKappaV=8.492E+03
         end
    5.5: begin
          Drainemodel='source/Draine5p5.txt'
          DraineKappaV=7.313E+03
         end
    ELSE: begin
            PRINT, 'For Dmodel, R_V should be 3.1 or 4.0 or 5.5. If you need other Rv values please choose Fmodel'
            stop
          end
  ENDCASE
  readcol,Drainemodel,DrainearrLamu,drainealbedo,drainecos,draineC,DrainearrKu,drainecos2 
  DrainearrLamu=DrainearrLamu*1.0E+4  ;um to A
  DrainearrKu=DrainearrKu/DraineKappaV ;normalizing kappa according to Draine(2003) definition
  DrainearrLam=DrainearrLamu[sort(DrainearrLamu)] ;Sorting lambda
  DrainearrK=DrainearrKu[sort(DrainearrLamu)]
endif


  
  outputim     = project_name+'/'+project_name+'_image.fits'
  outputpsf     = project_name+'/PSF_'+project_name+'_image.fits'
  outputimnoise=project_name+'/'+project_name+'_imageNoise.fits'
  outputspecFL = project_name+'/'+project_name+'_cube_spectra.fits'
  outputspecL  = project_name+'/'+project_name+'_Lambda.txt'
  outputstarinfo  = project_name+'/'+project_name+'_star_info.txt'
  openw,lunstarinfo,outputstarinfo,/get_lun
  printf,lunstarinfo,'#     mass[mo] ,   logage[yr]  ,     Z      , log[Teff[k]] ,    logg    ,   logL/Lo    ,  Av_star   ,      mag     ,    Xpix     ,   Ypix       ,       assignedSED'
  printf,lunstarinfo,'#       1              2             3            4              5              6           7               8            9            10                   11'
  xpix=round(fovx/res)
  ypix=round(fovy/res)
 ;Fov=xpix*xres
  print,'FoV: ',fovx,'" x ',fovy,'" =',xpix,'[pix] x',ypix,'[pix]'
  sceneim=fltarr(xpix,ypix)
;  pc2pix=206264.806247/distance/res
  
  if (spectroscopy eq 'yes') then begin
    myso_logo,'spec'
    Lspecarr=lminspec
    repeat Lspecarr=[Lspecarr,Lspecarr[-1]+Lspecarr[-1]/Rspec] until (Lspecarr[-1] ge lmaxspec)
    nspec=(size(Lspecarr))[1]
    sceneL=fltarr(nspec)
    sceneimFL=fltarr(xpix,ypix,nspec)
    openw,lun,outputspecL,/get_lun
    printf,lun,'# Lambda[A]'
    for ll=0,(size(lspecarr))[1]-1 do printf,lun,lspecarr[ll]    
    close,lun
    free_lun,lun
  endif
 

  ;reading the stars:
  myso_logo,'star'
  readcol,filestar,xstar,ystar,zstar,vxstar,vystar,vzstar,massstar,logagestar,kzstar,rhostar, /silent
;  logagestar=alog10(agestar)
  nstar=(size(xstar))[1]
  Rott=euler_to_rotmatrix([alphai,bettai,gammai],/DEGREE)
  zeinab=[[xstar],[ystar],[zstar]]#Rott
  xstar=zeinab[*,0]
  ystar=zeinab[*,1]
  zstar=zeinab[*,2]

  zeinabv=[[vxstar],[vystar],[vzstar]]#Rott
  vxstar=zeinabv[*,0]
  vystar=zeinabv[*,1]
  vzstar=zeinabv[*,2]
    
  distancestar=distance-zstar
  pc2pixstar=206264.806247/distancestar/res
  Teffstar=fltarr(nstar)
  loggstar=fltarr(nstar)
  loglstar=fltarr(nstar)
  sedstar=intarr(nstar)
  fluxstar=fltarr(nstar)
  newx=xstar*pc2pixstar ;convert x[pc] into pixel position
  newy=ystar*pc2pixstar
  newz=zstar*pc2pixstar
  columncloud=fltarr(nstar) ;column density of the cloud in fron of each star
  print,'READ STARs: ',nstar

 if (Columndensities eq 'sph') then begin
    myso_logo,'cloud'
    ;reading the cloud:
    readcol,filecloud,xcloud,ycloud,zcloud,vxcloud,vycloud,vzcloud,masspar,hpar,/silent
    ncloud=(size(xcloud))[1]
  
   Rott=euler_to_rotmatrix([alphai,bettai,gammai],/DEGREE)
   zeinab=[[xcloud],[ycloud],[zcloud]]#Rott
   xcloud=zeinab[*,0]
   ycloud=zeinab[*,1]
   zcloud=zeinab[*,2]
  
   distancecloud=distance-zcloud
   pc2pixcloud=206264.806247/distancecloud/res
   newxcloud=xcloud*pc2pixcloud ;convert x[pc] into pixel position
   newycloud=ycloud*pc2pixcloud
   newzcloud=zcloud*pc2pixcloud
   newhcloud=hpar*pc2pixcloud
  ;print,'READ CLOUD: ',ncloud
  endif
  ;
  ;reading the filters:
  myso_logo,'filter'
readcol, filter,lambda,weight,/silent
  ;print,'READ FILTER: ',(size(lambda))[1]

  ;reading the isochrones
  myso_logo,'iso'
  readcol,evolutionary,ziso,logageiso,miniiso,mactiso,logliso,logteff,loggiso,mboliso,UXiso,BXiso,Biso,Viso,Riso,Iiso,Jiso,Hiso,Kiso,Liso,Lpiso,Miso,/silent
  teffiso=10.^logteff
;  logageiso=alog10(ageiso) ;
  niso=(size(miniiso))[1]
  ;print,'READ Isochrones: ',niso

  ;reading list of SEDs:
  readcol,foldersed+'kseds.dat',teffsed,loggsed,metallicity,lh,vtur,sedname,FORMAT='F,F,F,F,F,A',/silent
  nseds=(size(teffsed))[1]
;  print,'READ SEDs: ',nseds

  if (OBtreatment eq 'yes') then begin
    readcol,foldersed+'Tseds.dat',teffsedOB,loggsedOB,metallicityOB,sednameOB,FORMAT='F,F,F,A',/silent
    nsedsOB=(size(teffsedOB))[1]
;    print,'READ TLUSTY SEDs: ',nsedsOB
     myso_logo,'tlusty'
  endif

  if (Adaptiveoptics eq 'yes') then begin
    myso_logo,'ao'
    print,'SR: ',sr,'  Seeing: ',seeing
  endif

;Vega mag system:
readcol,'SEDs/vegaf',wavelengthvega,fluxvega,/silent
linterp,lambda,weight,wavelengthvega,wavelength_weightvega
par2vega=0.0
for i=1, (size(wavelengthvega))[-1]-1 do begin
  par2vega += fluxvega[i]*wavelengthvega[i]*wavelength_weightvega[i]*(wavelengthvega[i]-wavelengthvega[i-1])
endfor
;par2vega=par2vega*6.247e-17 ;; if using Vega surface flux (F) then this will be the conversion: f=F*(R/d)^2


  nstars=0
  for ii=0,nstar-1 do begin
    if ((abs(newx[ii]) lt (xpix/2)-1) and (abs(newy[ii]) lt (ypix/2)-1)) then begin
      nstars += 1
      
      if (Columndensities eq 'sph') then begin
        columncloud[ii]=COLUMN_DENSITY(newx[ii],newy[ii],newz[ii],newxcloud,newycloud,newzcloud,masspar,newhcloud)
        cloudden=columncloud[ii]*((1.989/1.673534)/(((3.0857/pc2pixstar[ii])^2)));*1.0e21 ;convert column density unit [Msun/pix^2] --> [10^21 * Mhydrogen/cm^2]
      endif else cloudden=rhostar[ii]*((1.989/1.673534)/((3.0857^2)));*1.0e21 ;convert column density unit [Msun/pc^2] --> [10^21 * Mhydrogen/cm^2]

      AVstar=cloudden/2.21 ;e21 ;Guver&Ozel2009: The relation between Optical Extinction and Hydrogen column density
      ;      print,'Column density of the cloud infront of the star: ',cloudden, 'x 1.0e21 NH/cm^2  ===> AV star:',AVstar

      nage=FINDCLOSE(logagestar[ii],logageiso)
      selectedage=logageiso[nage]
      nmetalicity=FINDCLOSE(kzstar[ii],ziso)
      selectedz=ziso[nmetalicity]
      marrtemp=fltarr(niso)
      marrtemp[*]=999.99
      for kk=0,niso-1 do if ((ziso[kk] eq selectedz) and (logageiso[kk] eq selectedage)) then marrtemp[kk]=mactiso[kk]
      ns=FINDCLOSE(massstar[ii],marrtemp)
      Teffstar[ii]=teffiso[ns]
      loggstar[ii]=loggiso[ns]
      loglstar[ii]=logliso[ns]

      deltaT=abs(Teffstar[ii]-teffsed)

      deltagarr=fltarr(nseds)
      deltagarr[*]=99.
      
      for jj=0,nseds-1 do if ( deltaT[jj] eq min(deltaT) ) then deltagarr[jj]=abs(loggstar[ii]-loggsed[jj])
      sedstar[ii]=where(deltagarr eq min(deltagarr))      
      ;print,'sed: ',teffsed[jj],loggsed[jj]
      
      readsed=sedname[sedstar[ii]]
      readcol,foldersed+sedname[sedstar[ii]],wavelength,flux,/silent

      if ((OBtreatment eq 'yes') and (Teffstar[ii] gt 15000.)) then begin
        deltaT=abs(Teffstar[ii]-teffsedOB)

        deltagarr=fltarr(nsedsOB)
        deltagarr[*]=99.
        
        for jj=0,nsedsOB-1 do if (deltaT[jj] eq min(deltaT)) then deltagarr[jj]=abs(loggstar[ii]-loggsedOB[jj])
        sedstar[ii]=where(deltagarr eq min(deltagarr))
        
        readsed=sednameOB[sedstar[ii]]
        readcol,foldersed+sednameOB[sedstar[ii]],wavelength,flux,/silent
        ;        print,Teffstar[ii],loggstar[ii],foldersed+sednameOB[sedstar[ii]]
      endif

      bc1=BCcal(wavelength,flux,lambda,weight,AVstar,Rv,Teffstar[ii],par2vega,EXTmodel,DrainearrLam,DrainearrK)
      mag=4.77-2.5*loglstar[ii]-(bc1)+5.0*alog10(distancestar[ii]/10.0)
      ;print,bc1,mag
      fluxstar[ii]=10.^(mag/(-2.5))
      
      ;User may need this information:
      printf,lunstarinfo,massstar[ii], logagestar[ii],kzstar[ii],alog10(teffstar[ii]),loggstar[ii],loglstar[ii],AVstar,mag,newx[ii]+xpix/2.,newy[ii]+ypix/2.,readsed,$
        form='(f13.4,1x,f13.4,1x,f13.4,1x,f13.4,1x,f13.4,1x,f13.4,1x,f13.4,1x,f13.4,1x,f13.4,1x,f13.4,1x,a60)'

      if (PSFtype eq 'gaussian') then $
      airy1=psf_gaussian(Npixel=[xpix,ypix],fwhm=[fwhm/res,fwhm/res],CENTROID =[newx[ii]+xpix/2.,newy[ii]+ypix/2.],/NORMALIZE,/Double) $
      else if (PSFtype eq 'airy') then airy1=airy_pattern(xpix,ypix,newx[ii]+xpix/2.,newy[ii]+ypix/2.,2./(fwhm/res))
      
      airy1=airy1/(total(airy1)) ;to be sure about normaized total flux

      if (Adaptiveoptics eq 'yes') then begin
        halo=psf_gaussian(Npixel=[xpix,ypix],fwhm=[seeing/res,seeing/res],CENTROID =[newx[ii]+xpix/2.,newy[ii]+ypix/2.],/NORMALIZE,/Double)
        halo = halo/total(halo)
        sceneim += fluxstar[ii]*(SR*airy1 + (1.0-SR)*halo)
      endif

      if (Adaptiveoptics eq 'no') then sceneim += fluxstar[ii]*airy1

      if (spectroscopy eq 'yes') then begin
        for ll=2,nspec-3 do begin
          ;          sceneL[ll]=Lspecarr[ll]
          
          linterp,lambda,weight,wavelength,wavelength_weight ;for spectroscopy you need 2 interpolation:1)filter 2)spectra-windows
          bc3=BCcals(wavelength,flux*wavelength_weight,[Lspecarr[ll-1],Lspecarr[ll],Lspecarr[ll+1],Lspecarr[ll+2]],[0.0,1.0,1.0,0.0],AVstar,Rv,Teffstar[ii],EXTmodel,DrainearrLam,DrainearrK)
          mag3=4.77-2.5*loglstar[ii]-bc3+5.0*alog10(distancestar[ii]/10.0)
          fluxstar3=double(10.^(mag3/(-2.5)))
          if (PSFtype eq 'gaussian') then $
            airy3=psf_gaussian(Npixel=[xpix,ypix],fwhm=[fwhm/res,fwhm/res],CENTROID =[newx[ii]+xpix/2.,newy[ii]+ypix/2.],/NORMALIZE,/Double) $
          else if (PSFtype eq 'airy') then airy3=airy_pattern(xpix,ypix,newx[ii]+xpix/2.,newy[ii]+ypix/2.,2./(fwhm/res))


          airy3=airy3/(total(airy3)) ;to be sure about normaized total flux

          if (velocitydis eq 'yes') then begin
            shiftv=vzstar[ii]*Lspecarr[ll]/3.0e+5 ;shift is in A
            lambdashift=Lspecarr[ll]-shiftv  ;if vz>0 ==> source comes toward observer ==> lambdashift<lamda0 (blue-shift)
            llchannel=(sort(abs(lambdashift-Lspecarr)))[0]
            if (Adaptiveoptics eq 'yes') then sceneimFL[*,*,llchannel] += fluxstar3*(SR*airy3+ (1.0-SR)*halo)
            if (Adaptiveoptics eq 'no' ) then sceneimFL[*,*,llchannel] += fluxstar3*airy3
          endif
          
          if (velocitydis eq 'no') then begin
            if (Adaptiveoptics eq 'yes') then sceneimFL[*,*,ll] += fluxstar3*(SR*airy3+ (1.0-SR)*halo)
            if (Adaptiveoptics eq 'no' ) then sceneimFL[*,*,ll] += fluxstar3*airy3
          endif
           
        endfor
      endif

    endif
;    tvscl,-asinh(sceneim^0.1)
  endfor
  print,'Number of stars in the FoV: ',nstars
  writefits,outputim,sceneim
  faintestflux=min(fluxstar(where(fluxstar gt 0.0)))
  if (SNR ne 0.0) then noise=faintestflux*4.*alog(2.0)/(!PI*(fwhm/res)*(fwhm/res))/SNR else noise=noise2add
  ;print,noise,mean(sceneim)
  writefits,outputimnoise,sceneim+noise*randomu(seed,xpix,ypix)
  ;writefits,outputspecL,sceneL
  if (spectroscopy eq 'yes') then  writefits,outputspecFL,sceneimFL  

  printf,lunstarinfo,'#faintestflux: ',faintestflux
  printf,lunstarinfo,'#noise: ',noise
  printf,lunstarinfo,'#SNR: ',SNR
  printf,lunstarinfo,'#Filter: ',filter
  printf,lunstarinfo,'#distance: ',distance
  printf,lunstarinfo,'#res: ',res
  printf,lunstarinfo,'#FoV: ',fovx,fovy
  printf,lunstarinfo,'#fwhm: ',fwhm
  printf,lunstarinfo,'#OBtreatment: ',OBtreatment
  printf,lunstarinfo,'#Adaptiveoptics: ',Adaptiveoptics
  printf,lunstarinfo,'#SR: ',SR
  printf,lunstarinfo,'#Seeing: ',seeing
  printf,lunstarinfo,'#Spectroscopic Resolution: ',Rspec
  printf,lunstarinfo,'#Lambda_min: ',lminspec 
  printf,lunstarinfo,'#Lambda_max: ',lmaxspec
  printf,lunstarinfo,'#filestar: ',filestar
  printf,lunstarinfo,'#Columndensities: ',Columndensities
  printf,lunstarinfo,'#filecloud: ',filecloud
  printf,lunstarinfo,'#EXTmodel: ',EXTmodel
  printf,lunstarinfo,'#Rv: ',Rv

  close,lunstarinfo
  free_lun,lunstarinfo
  myso_logo,'outim'
  print, outputim
  print, '   '
  myso_logo,'outspec'
  print, outputspecFL
  print, outputspecL
  print, '   '
  if (PSFtype eq 'gaussian') then airypsf=psf_gaussian(Npixel=[16*round(fwhm/res),16*round(fwhm/res)],fwhm=[fwhm/res,fwhm/res],CENTROID =[(16*round(fwhm/res))/2.,(16*round(fwhm/res))/2.],/NORMALIZE,/Double)
  if (PSFtype eq 'airy') then airypsf=airy_pattern(16*round(fwhm/res),16*round(fwhm/res),(16*round(fwhm/res))/2.,(16*round(fwhm/res))/2.,2./(fwhm/res))
  writefits,outputpsf,airypsf
end














