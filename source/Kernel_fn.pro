; This subroutine is an implementation of the kernel function(W)
; It takes distance from the center of some particle(dr) and the particle smoothing length(h) and it returns the kernel function (w).
; Uses 3D M4 algorithem 
; Inputs: 
;     dr: distance from the center of the cloud particle
;     h: smoothing length of the cloud particle
; Output: 
;     w: Kernel function value
;
; Calling sequence: 
;     w = Kernel_fn(dr,h)         
;
;  written by Z. Khorrami (KhorramiZ@cardiff.ac.uk)
;  updated Feb. 2017

FUNCTION Kernel_fn,dr, h
  pi=3.14159265358

  if (dr/h lt 1.0) then begin
    w = (1/(pi*h^3)) * (1-1.5*(dr/h)^2 + 0.75*(dr/h)^3)
  endif else if (dr/h lt 2.0) then begin
    w = (1/(pi*h^3)) * (0.25 * (2.0-(dr/h))^3)
  endif else w = 0.0

  return,w

END