# RadialProfiles

For both radio and X-ray images:

-in ds9 create FirstPanda.reg (you can also modify the example provided here)

-run make_CircularPanda_Spitted.py

Input: FWHM of beam/psf, annulus width in units of the FWHM, maximum radial distance

################

For radio images 

-run Radio_sb.py

Input: image in fits format, rms in Jy/beam

################

For X-ray images

-run Xray_sb.py

Input: instrument ("XMM", "CHANDRA"), count/exposure/background maps in fits format 



