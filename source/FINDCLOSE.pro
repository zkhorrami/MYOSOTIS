; Finding closest value (its element) to the input
;Inputs:
; par1 : float value
; pararr: an array
; Output:
; Element of the closest value in pararr to par1
; 
; Calling sequence:
;   n = findclose(par1,pararr)
;   
; Example:
;   $ array = [1, 2, 2.5, 3.5, 1.4]
;   $ n = findclose(1.5, array)
;   $ print, n
;   $ 4
;   $ print, array[n]
;   $ 1.4
;   
;  written by Z. Khorrami (KhorramiZ@cardiff.ac.uk)
;  updated Feb. 2017

FUNCTION FINDCLOSE,par1, pararr
  nclose= (sort(abs(par1-pararr)))[0]
  RETURN, nclose
END