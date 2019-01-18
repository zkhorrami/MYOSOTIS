;This function calculates the extinction in a given wavelength knowing optical extinction value.
;this function uses The average extinction curve in the optical-through-IR range is
;reproduced with a cubic spline and a set of anchor points from Fitzpatrick (1999)
;Works only in the wavelength range of 0.125 - 3.333 um
;
;Draine model covers [0.0001 - 10000] um and the references are here:
;Draine, B.T. 2003a, "Interstellar Dust Grains", Annu. Rev. Astr. Astrophys.,41, 241-289
;Draine, B.T. 2003b, "Scattering by Interstellar Dust Grains. I. Optical and Ultraviolet", Astrophys. J., 598, 1017-1025
;Draine, B.T. 2003c, "Scattering by Interstellar Dust Grains. II. X-Rays", Astrophys. J., 598, 1026-1037
;Li, A., & Draine, B.T. 2001, "Infrared Emission from Interstellar Dust. II. The Diffuse Interstellar Medium", Astrophys. J., 554, 778-802
;Weingartner, J.C., & Draine, B.T. 2001, Dust Grain Size Distributions and Extinction in the Milky Way, LMC, and SMC", Astrophys. J., 548,
;296-309 (WD01).
;
;INPUTs: Lambda [um]
;        Av
;        Rv
;OUTOUT: A_{lambda} array
;
;Calling sequence: extinctions(lamarr,Av,Rv)
;
;  written by Z. Khorrami (KhorramiZ@cardiff.ac.uk)
;  updated Feb. 2017

FUNCTION extinctions,lam,Av,Rv;,EXTmodel,DrainearrLam,DrainearrK,DraineKappaV
alam=0.0  
;if (EXTmodel eq 'Fmodel') then begin
  xarr=1./lam
  
  if ((xarr ge 0.3) and (xarr lt 1.1)) then begin
    ax1=0.574*(xarr^1.61)
    bx1=-0.527*(xarr^1.61)
    Alam=Av*(ax1+bx1/Rv)
  endif else begin
    if ((xarr ge 1.1) and (xarr lt 3.3)) then begin
      y = xarr-1.82
      ax2 = 1.0+0.17699*y-0.50447*(y^2.0)-0.02427*(y^3.0)+0.72085*(y^4.0) +0.01979*(y^5.0)-0.77530*(y^6.0)+0.32999*(y^7.0)
      bx2 = 1.41338*y+2.28305*(y^2.0)+1.07233*(y^3.0)-5.38434*(y^4.0)-0.62251*(y^5.0)+5.30260*(y^6.0) -2.09002*(y^7.0)
      Alam=Av*(ax2+bx2/Rv)
      ;   print,xarr[ii],alam[ii],y,ax2,bx2
    endif else begin
      if ((xarr ge 3.3) and (xarr lt 5.9)) then begin
        ax3 = 1.752 -0.316*xarr-0.104/((xarr-4.67)^2.0 + 0.341)
        bx3 = (-3.090)+1.825*xarr+1.206/((xarr-4.62)^2.0+0.263)
        Alam=Av*(ax3+bx3/Rv)
        ;   print,xarr[ii],alam[ii],ax3,bx3
      endif else begin
        if ((xarr ge 5.9) and (xarr lt 8.0)) then begin
          fax4 = -0.04473*(xarr - 5.9)^2.0 - 0.009779*(xarr -5.9)^3.0
          fbx4 = 0.2130*(xarr-5.9)^2.0 + 0.1207*(xarr-5.9)^3.0
          ax4 = 1.752 - 0.316*xarr - (0.104/((xarr - 4.67)^2.0 + 0.341)) +fax4
          bx4 = -3.090 + 1.825*xarr + (1.206/((xarr - 4.62)^2.0 + 0.263)) +fbx4
          Alam=Av*(ax4+bx4/Rv)
          ;   print,xarr[ii],alam[ii],ax4,bx4
        endif else  print,'Lambda,',lam,', is out of range for Fitzpatrick (1999) model!!!'
      endelse
    endelse
  endelse
;  print,alam
  return,alam
;endif 

;if (EXTmodel eq 'Dmodel') then begin
;  ncloselam=FINDCLOSE(lam,Drainearrlam)
;  alam=Av*DrainearrK[ncloselam]/DraineKappaV
;  return,alam
;endif

END