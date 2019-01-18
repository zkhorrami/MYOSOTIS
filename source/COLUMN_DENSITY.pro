; Calculate column density in front of a source
; Column density is calculated knowing the cloud's particle position, mass and smoothing lengths
; This code can be easily used for the output of the SPH simulations where we have all the information of a cloud.
;
;Inputs:
;   Source position: Xstar, ystar, Zstar
;   Cloud position: Xarr, Yarr, Zarr
;   Cloud's particle masses and smoothing lengths: marr, harr
;   source position and cloud position array and smoothing length array should have a same unit
;
; Output:
;   the total column density infrint of the point source
;   The unit of the output is the [input mass unit/input length unit^2]
;   For example if mass_array is in Solar_mass and length in pc
;   then the column density unit is [solar_mass/pc^2]
;   
;
; Calling sequence:
;   columndensity = COLUMN_DENSITY(xstar,ystar,zstar,xarr,yarr,zarr,marr,harr)
;
;The output can be used to detemin the Av according to the column density infront of a given source
;
;  written by Z. Khorrami (KhorramiZ@cardiff.ac.uk)
;  updated Feb. 2017
  
FUNCTION COLUMN_DENSITY,xstar,ystar,zstar,xarr,yarr,zarr,marr,harr
  On_error,2
  compile_opt idl2

  nc=(size(marr))[1]
  den2dtot=0.0
  for ii=0,nc-1 do begin
    if ( (sqrt((xstar-xarr[ii])^2+(ystar-yarr[ii])^2) lt 2.*harr[ii]) and (zstar lt zarr[ii]) ) then begin
      dr=sqrt((xstar-xarr[ii])^2+(ystar-yarr[ii])^2)
      lc=2*sqrt(4.*harr[ii]*harr[ii]-dr*dr)
      myres=20
      deltalc=lc/myres
      denarr=fltarr(myres)
      lcarr=fltarr(myres)
      for jj=0,myres-1 do denarr[jj]= marr[ii]*kernel_fn(sqrt(dr*dr+(jj*lc/myres-lc/2.)^2.),harr[ii])*deltalc
      den2dtot = den2dtot+ total(denarr)
    endif
  endfor
  RETURN, den2dtot
END
