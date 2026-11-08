import sys
import os
import numpy as np
import math
from astropy.io import fits
import pyregion
import matplotlib.pyplot as plt
import glob
from astropy.wcs import WCS
from astropy.wcs import WCS as pywcs


"""
USAGE

- first run make_CircularPanda_Splitted.py
- INPUT: radio map, rms

"""

radio_map = "map143_SUB_134x134_4axis.fits.mask2.fits"
sigma = 1.00e-3 #Jy/beam

#flatten
def flatten(filename, channel=0, freqaxis=0):
    f = fits.open(filename)
    naxis=f[0].header['NAXIS']
    if naxis<2:
        raise RadioError('Cannot make map from this')
    if naxis==2:
        pass
        #return f[0].header,f[0].data

    w = pywcs(f[0].header)
    wn = pywcs(naxis=2)

    wn.wcs.crpix[0]=w.wcs.crpix[0]
    wn.wcs.crpix[1]=w.wcs.crpix[1]
    wn.wcs.cdelt=w.wcs.cdelt[0:2]
    wn.wcs.crval=w.wcs.crval[0:2]
    wn.wcs.ctype[0]=w.wcs.ctype[0]
    wn.wcs.ctype[1]=w.wcs.ctype[1]

    header = wn.to_header()
    header['NAXIS']=2
    header['NAXIS1']=f[0].header['NAXIS1']
    header['NAXIS2']=f[0].header['NAXIS2']
    copy=('EQUINOX','EPOCH')
    for k in copy:
        r=f[0].header.get(k)
        if r:
            header[k]=r

    dataslice=[]
    for i in range(naxis,0,-1):
        if i<=2:
            dataslice.append(np.s_[:],)
        elif i==freqaxis:
            dataslice.append(channel)
        else:
            dataslice.append(0)
# add freq
    header["FREQ"] = f[0].header['CRVAL3']

    # add beam if present
    try:
        header["BMAJ"]=f[0].header['BMAJ']
        header["BMIN"]=f[0].header['BMIN']
        header["BPA"]=f[0].header['BPA']
    except:
        pass
    return header, f[0].data[tuple(dataslice)]



# read data from radio map
def read_pixel_value(reg):
    hdul = fits.open(radio_map+'temp')
    header = hdul[0].header
    data = hdul[0].data

    x_size = header['NAXIS1']
    y_size = header['NAXIS2']
    
    r = pyregion.open('REGIONS/'+reg).as_imagecoord(header)
    
    
    mask = r.get_mask(shape=(y_size, x_size))
    data_inside_reg = data[mask]
    
    return data_inside_reg



head, data = flatten(radio_map)
pixscale = abs(head['CDELT1']*3600.) #arcsec
bmaj = head['BMAJ']*3600
bmin = head['BMIN']*3600
Abeam_pix = (np.pi * bmaj * bmin) / (4 * np.log(2) * pixscale**2)
Abeam_arcsec = (np.pi * bmaj * bmin) / (4 * np.log(2))
first_reg = pyregion.open('FirstPanda.reg')
err_radius = (first_reg[0].coord_list[6] - first_reg[0].coord_list[5])*3600./2.


if os.path.exists(radio_map+'temp'):
    os.remove(radio_map+'temp')
fits.writeto(radio_map+'temp', data, head)


regions = []
region_list = open('listREG.txt')

for line in region_list:
    regions.append(line)
region_list.close()


if os.path.exists('SBR_output.txt'):
    os.remove('SBR_output.txt')

i = 0
with open('SBR_output.txt', 'a') as f:
    f.write('#  reg_number, meanR, SB, errSB \n ')
    f.write('# ---, arcsec, Jy/arcsec2, Jy/arcsec2 \n ')
    
    for line in regions:
        reg_name = line.strip(' \n')
        reg_name = str(reg_name)
        data_inside_reg = read_pixel_value(reg_name)
        Npix_all = len(data_inside_reg)
        check_nan = np.isnan(data_inside_reg)
        Npix_nan = 0
        for j in check_nan:
            if j==True:
                Npix_nan = Npix_nan + 1
        Npix = Npix_all - Npix_nan
        Nbeam = Npix/Abeam_pix
        calib_error = 0.1
        area = Npix*pixscale**2
        region = pyregion.open('REGIONS/'+reg_name)
        
        flux = np.nansum(data_inside_reg)/Abeam_pix
        flux_err = np.sqrt( (sigma*np.sqrt(Nbeam))**2 + (calib_error*flux)**2 )
        
        SB = flux/area
        errSB = flux_err/area
        mean_radius = (region[0].coord_list[5] + region[0].coord_list[6])*3600./2.
        
        
        f.write(str(i) +str(',') + str(mean_radius) +str(',') + str(SB) +str(',') + str(errSB) + '\n')
        
        i = i+1
        
f.close()
os.remove(radio_map+'temp')





#sys.exit()

fileR = np.loadtxt('SBR_output.txt', delimiter=',', skiprows=2, usecols=(1,2,3))
mean_radius = fileR[:,0]
SBR = fileR[:,1]
errSBR = fileR[:,2]


plt.errorbar(x=mean_radius, y=SBR*1.e6, xerr=err_radius, yerr=errSBR*1.e6, fmt='ro', ecolor='r', capsize=3)



#r_lorentz1D = np.arange(201, 804, 10)
#I_lorentz1D = 3.e-6 * np.exp(-aaa/170 )
#plt.plot(r_lorentz1D, I_lorentz1D*1.e6, 'k--')


plt.ylabel('SB$_{R}$ ($\mu Jy \; arcsec^{-2}$)')
plt.xlabel('r ($arcsec$) ')

plt.title('RADIO azimuthally-averaged surface brightness profile')

plt.xscale('log')
plt.yscale('log')

plt.savefig('Radio_SB_profile.pdf',dpi=300)

plt.show()


#sys.exit()




