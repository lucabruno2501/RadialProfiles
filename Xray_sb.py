import sys
import os
import numpy as np
import math
from astropy.io import fits
import pyregion
import matplotlib.pyplot as plt
import glob


"""
USAGE

- first run make_CircularPanda_Splitted.py or make_boxreg.py or make_EllipticalPanda_Splitted.py
- INPUT: countmap, exposure map

"""

reg_type = "panda"
side = "minor"
instrument = "XMM"
countmap = "mosaic_a2142.fits.mask.fits"
exposuremap = "mosaic_a2142_expo.fits.mask.fits"
blank_files = glob.glob("mosaic_a2142_bkg.fits.mask.fits")




regions = []
region_list = open('listREG.txt')

for line in region_list:
    regions.append(line)
region_list.close()
    


# read data from count map
def read_pixel_value(reg):
    hdul = fits.open(countmap)
    header = hdul[0].header
    data = hdul[0].data

    x_size = header['NAXIS1']
    y_size = header['NAXIS2']
    
    r = pyregion.open('REGIONS/'+reg).as_imagecoord(header)
    
    
    mask = r.get_mask(shape=(y_size, x_size))
    data_inside_reg = data[mask]
    
    return data_inside_reg


# read data from exposure map
def read_exposure(reg):
    hdul = fits.open(exposuremap)
    header = hdul[0].header
    data = hdul[0].data
    
    x_size = header['NAXIS1']
    y_size = header['NAXIS2']

    r = pyregion.open('REGIONS/'+reg).as_imagecoord(header)


    mask = r.get_mask(shape=(y_size, x_size))
    exposure_inside_reg = data[mask]

    return exposure_inside_reg
  

# read data from combined blanksky files
def read_background(reg):
    Nback = 0.
    for files in blank_files:
    
        hdul = fits.open(files)
        header = hdul[0].header
        data = hdul[0].data

        x_size = header['NAXIS1']
        y_size = header['NAXIS2']

        r = pyregion.open('REGIONS/'+reg).as_imagecoord(header)

        mask = r.get_mask(shape=(y_size, x_size))
        backgound_inside_reg = data[mask]
        
        Nback = Nback + np.nansum(backgound_inside_reg)
        

    
    return Nback



hdul = fits.open(countmap)
header = hdul[0].header
pixscale = abs(header['CDELT1']*3600.) #arcsec

if reg_type=='epanda':
    first_reg = pyregion.open('FirstPanda.reg')
    err_radius = (first_reg[0].coord_list[7] - first_reg[0].coord_list[5])*3600./2.
if reg_type=='panda':
    first_reg = pyregion.open('FirstPanda.reg')
    err_radius = (first_reg[0].coord_list[6] - first_reg[0].coord_list[5])*3600./2.
if reg_type=='box' and side=='minor':
    first_reg = pyregion.open('FirstBox.reg')
    err_radius = first_reg[0].coord_list[3] *3600./2.
if reg_type=='box' and side=='major':
    first_reg = pyregion.open('FirstBox.reg')
    err_radius = first_reg[0].coord_list[2] *3600./2.

if os.path.exists('SBX_output.txt'):
    os.remove('SBX_output.txt')
with open('SBX_output.txt', 'w') as f:
    f.write('#reg_number, meanR, errR, SB, errSB \n ')
    f.write('# ---, arcsec, arcsec,cnt/s cm2 arcsec2, cnt/s cm2 arcsec2 \n ')
f.close()

i=0
for line in regions:
    reg_name = line.strip(' \n')
    reg_name = str(reg_name)

    data_inside_reg = read_pixel_value(reg_name)
    
#    Ncount = np.sum(data_inside_reg)
    
    
#    Npix = len(data_inside_reg)
    Npix_all = len(data_inside_reg)
    check_nan = np.isnan(data_inside_reg)
    Npix_nan = 0
    for j in check_nan:
        if j==True:
            Npix_nan = Npix_nan + 1
    Npix = Npix_all - Npix_nan
    
    Ncount = np.nansum(data_inside_reg)
    
    area = Npix*pixscale**2
    
    exposure_inside_reg = read_exposure(reg_name)
        
    exposure = np.nansum(exposure_inside_reg)/(len(exposure_inside_reg)-Npix_nan)
    
    Nback = read_background(reg_name)
    
    flux = (Ncount-Nback)/exposure
    errflux = np.sqrt(Ncount+Nback)/exposure
    
    
    region = pyregion.open('REGIONS/'+reg_name)
    
    if reg_type=='panda':
        SB = flux/area
        errSB = errflux/area
        mean_radius = (region[0].coord_list[5] + region[0].coord_list[6])*3600./2.
    if reg_type=='epanda':
        SB = flux/area
        errSB = errflux/area
        mean_radius = (region[0].coord_list[7] + region[0].coord_list[5])*3600./2.
    if reg_type=='box' and side=='minor':
        SB = flux/area
        errSB = errflux/area
        mean_radius = (region[0].coord_list[3] + i * region[0].coord_list[3])*3600
    if reg_type=='box' and side=='major':
        SB = flux/area
        errSB = errflux/area
        mean_radius = (region[0].coord_list[2] + i * region[0].coord_list[3])*3600
    
    with open('SBX_output.txt', 'a') as f:
        f.write(str(i) +str(',') + str(mean_radius) +str(',') + str(err_radius) +str(',') + str(SB) +str(',') + str(errSB) + '\n')
    i = i+1
f.close()


fileX = np.loadtxt('SBX_output.txt', delimiter=',', skiprows=2, usecols=(1,2,3,4))
mean_radius = fileX[:,0]
err_radius = fileX[:,1]
SBX = fileX[:,2]
errSBX = fileX[:,3]


plt.errorbar(x=mean_radius, y=SBX, xerr=err_radius, yerr=errSBX, fmt='ro', ecolor='r', capsize=3)

if instrument=='CHANDRA':
    plt.ylabel('SB$_{X}$ ($cnts \; s^{-1} \; cm^{-2} \; arcsec^{-2}$)')
if instrument=='XMM':
    plt.ylabel('SB$_{X}$ ($cnts \; s^{-1} \; arcsec^{-2}$)')
plt.xlabel('r ($arcsec$) ')
plt.title('X-RAY surface brightness profile')


plt.xscale('log')
plt.yscale('log')

plt.savefig('X-ray_SB_profile.pdf',dpi=300)

plt.show()


#sys.exit()



