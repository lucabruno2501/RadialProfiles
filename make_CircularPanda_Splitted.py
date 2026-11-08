#!/usr/bin/python

import numpy as np
import os, sys
import pyregion
from shutil import copyfile, move, rmtree
from shutil import move


#OUTPUT a number of panda circular regions to be used for SB analysis

"""
REQUIRED INPUT and PRELIMINARIES:
FWHM of beam (sampling with FWHM/2), max radius (in arcsec, multiple of FWHM/2)

in ds9 create a single bin sector of width=FWHM/2 and save it as 'FirstPanda.reg'

OUTPUT: list of regions + REGIONS/*.reg
"""


FWHM=134
width=FWHM/2.
r_max=804./3600. #arcsec




input_reg = 'FirstPanda.reg'
reg = pyregion.open(input_reg)
outputreg = 'CircularPanda_width'+str(round(width,1))+'arcsec.reg'

if os.path.exists(outputreg):
    os.remove(outputreg)
copyfile(input_reg, outputreg)


i_max= np.round(r_max/reg[0].coord_list[6],0)




if os.path.exists('listREG.txt'):
    os.remove('listREG.txt')
if os.path.exists('REGIONS'):
    rmtree('REGIONS')
os.mkdir('REGIONS')


i=0
while i<i_max:
    if os.path.exists('CircularPanda_width'+str(round(width,1))+'arcsec_'+str(i)+'.reg'):
        os.remove('CircularPanda_width'+str(round(width,1))+'arcsec_'+str(i)+'.reg')
    with open('CircularPanda_width'+str(round(width,1))+'arcsec_'+str(i)+'.reg', 'w') as f:
        f.write('# Region file format: DS9 version 4.1 \n')
        f.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1 \n')
        f.write('fk5 \n')
        f.write('panda(' +str(reg[0].coord_list[0]) +str(',') + str(reg[0].coord_list[1]) +str(',') +str(reg[0].coord_list[2]) +str(',') +str(reg[0].coord_list[3]) +str(',') +str(reg[0].coord_list[4]) +str(',') +str(reg[0].coord_list[5]+(i*width)) +str('",') +str(reg[0].coord_list[6]*3600+(i*width)) +str('",') +'1)')
    with open('listREG.txt', 'a') as g:
        g.write('CircularPanda_width'+str(round(width,1))+'arcsec_'+str(i)+'.reg \n')
    move('CircularPanda_width'+str(round(width,1))+'arcsec_'+str(i)+'.reg', './REGIONS')
    i=i+1
f.close()
g.close()

i=1
with open(outputreg, 'a') as h:
    while i<i_max:
        h.write('panda(' +str(reg[0].coord_list[0]) +str(',') + str(reg[0].coord_list[1]) +str(',') +str(reg[0].coord_list[2]) +str(',') +str(reg[0].coord_list[3]) +str(',') +str(reg[0].coord_list[4]) +str(',') +str(reg[0].coord_list[5]+(i*width)) +str('",') +str(reg[0].coord_list[6]*3600+(i*width)) +str('",') +'1)\n')
        i=i+1
h.close()
