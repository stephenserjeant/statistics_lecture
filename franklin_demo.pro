function expand_image, image

s = size(image)
bigimage = dblarr(3*s[1],3*s[2])
for i=0,2 do begin 
   for j=0,2 do begin 
      bigimage[i*s[1]:(i+1)*s[1]-1,j*s[2]:(j+1)*s[2]-1] = image[*,*]
   endfor 
endfor 

return, bigimage
end 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function dequarter, image

s = size(image)
result = 0*image 

result[0:s[1]/2-1,s[2]/2:s[2]-1] = image[s[1]/2:s[1]-1,0:s[2]/2-1]
result[s[1]/2:s[1]-1,s[2]/2:s[2]-1] = image[0:s[1]/2-1,0:s[2]/2-1]
result[0:s[1]/2-1,0:s[2]/2-1] = image[s[1]/2:s[1]-1,s[2]/2:s[2]-1]
result[s[1]/2:s[1]-1,0:s[2]/2-1] = image[0:s[1]/2-1,s[2]/2:s[2]-1]

return, result
end 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function convolve_periodic, image, kernel, _extra=e
bigimage = expand_image(image)
smooth_image = convolve(bigimage, kernel, _extra=e)
s = size(image)
result = smooth_image[s[1]:2*s[1]-1,s[2]:2*s[2]-1]
return, result 
end 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro franklin_demo 

;;; read in RF image 
read_jpeg, '~/idl/rosalind_franklin.jpg',image0
s = size(image0)
image = dblarr(s[1],s[2]+1)
image[*,1:s[2]] = image0
image[*,0] = image[*,0]

;;; smooth very slightly
kernel = psf_gaussian(fwhm=2.0,npix=[5,5])
image = convolve(image,kernel)

;;; make window
s = size(image)
window,0,xsize=s[1],ysize=s[2]

;;; plot image
tv, bytscl(image)
saveimage,'franklin1.jpg',/jpeg
print,'Press any key...'
jnk=get_kbrd(1)

;;; smooth image with a real-space summation 
fwhm=7.0
kernel = psf_gaussian(fwhm=fwhm,npix=[19,19])
smooth_image = convolve_periodic(image,kernel,/no_ft) 

;;; plot image
tv, bytscl(smooth_image)
saveimage,'franklin2.jpg',/jpeg
print,'Press any key...'
jnk=get_kbrd(1)

;;; restore using Fourier deconvolution 
kernel2 = psf_gaussian(fwhm=fwhm,npix=[s[1],s[2]])
input_image_fft = fft(smooth_image)
input_psf_fft = fft(kernel2)
ok = where(abs(input_psf_fft) gt 1d-7)
restored_image_fft = 0*input_image_fft
restored_image_fft[ok] = input_image_fft[ok] / input_psf_fft[ok]
restored_image = float(fft(restored_image_fft,/inverse))

;;; plot image
tv, bytscl(dequarter(restored_image))
saveimage,'franklin3.jpg',/jpeg

;;; closing 
stop 

end 


