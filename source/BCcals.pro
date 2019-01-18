; Calculate Bolometric Correction for a given filter considering the optical extinction
; 
;Inputs:
;   wavelength : input wavelength array [A]
;                to convert frequency [1/s] into wavelength [A] use: wavelength=2.997925e18/nu
;   flux: The stellar flux (SED) in different wavelengths [erg/s/cm2/A]. Sould has a same size as the wavelength.
;            to covert flux unit from [erg/s/cm2/Hz] into [erg/s/cm2/A] use flux = flux*2.997925e18/wavelength/wavelength*4*!pi
;   lambda: wavelength array of the filter [A]
;   weight: weight values of the filter (filter's transparency). Should has a same size as lambda
;   AV: optical extinction
;   Rv: R(V)=AV/E(B-V) both Av and Rv depend on the amount of gas/dust infront of the stellar source. If there is no extinction put Av=0
;            
; Output:
;   The BC of the star in a given filter 
;   For more details see Girardi + 2002
;   
; Element of the closest value in pararr to par1
; 
; Calling sequence:
;   BClambda = BCcal(lambda_SED, F_lambda_SED, lambda_filter, Filter_transparency, Av_value, Rv_value)
;   
;   
;  written by Z. Khorrami (KhorramiZ@cardiff.ac.uk)
;  updated Feb. 2017

FUNCTION BCcals,wavelength,flux,lambda,weight,AVstar,Rv,Teff,EXTmodel,DrainearrLam,DrainearrK;,par2vega; sed, filter
  ;Some needed constants:
  ;BC= Mbol_sun -2.5 log10 [4!PI 10pc^2 Sigma Teff^4 / Lsun] + 2.5 .... Girardi + 2002
  mbolsun = 4.77
;  bolconstant = -2.5*alog10(4*!pi*567.0373/(3.24078^2.*3.844))
  bolconstant = -2.5*alog10(4*!pi*(3.0857^2.0)*5.67051/3.826)
  ;Teff=5750.
 
  n_wave=(size(wavelength))[-1]
  ;reading the filters
  ;      readcol, filter,lambda,weight
  ;I put the zero weight at the edges just to solve the problem of edges in linear interpolation
  weight[0]=0.0
  weight[-1]=0.0
  linterp,lambda,weight,wavelength,wavelength_weight
  par1=0.0
  par2=0.0
  alamarr=fltarr((size(wavelength))[-1])
  
  if (EXTmodel eq 'Dmodel') then begin
    linterp,DrainearrLam,DrainearrK,wavelength,kappased
    alamarr=AVstar*kappased
  endif
  
  for i=1, (size(wavelength))[-1]-1 do begin
    ltemp=wavelength[i]*1.0e-4 ; don't foorget to convert A to um for extinction calculation!
   if (EXTmodel eq 'Fmodel' and wavelength[i] ge 1250. and wavelength[i] le 33333.) then alamarr[i]=extinctions(ltemp,Avstar,Rv);,'Fmodel',DrainearrLam,DrainearrK,DraineKappaV)

    par1 += wavelength[i]*flux[i]*(10.0^(-0.4*alamarr[i]))*wavelength_weight[i]*(wavelength[i]-wavelength[i-1])
    par2 += wavelength[i]*wavelength_weight[i]*(wavelength[i]-wavelength[i-1])  ;!!!! par2 will be calculating from Vega flux
   endfor
  BCfilter = Mbolsun + bolconstant - 10.*alog10(teff)+2.5*alog10(par1/par2);+21.10
  RETURN, BCfilter
END