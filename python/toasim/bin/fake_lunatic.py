#!/usr/bin/env python


import sys
from astropy import coordinates as coord
from astropy import units as u
from astropy import time as atime
from astropy import constants as const
import numpy as np

import toasim

nreal = 1
header = toasim.header()

header.parfile_name=sys.argv[1]
header.timfile_name=sys.argv[2]


with open(sys.argv[1]) as par:
    header.orig_parfile=par.read()
with open(sys.argv[2]) as tim:
    header.idealised_toas=tim.read()




with open(header.parfile_name) as f:
    for line in f:
        e=line.split()
        if e[0]=="RAJ":
            raj=e[1]
        elif e[0]=="DECJ":
            decj=e[1]

toas=[]
with open(header.timfile_name) as f:
    for line in f:
        if line.startswith("C") or line.startswith("#"):
            continue
        e=line.split()
        if len(e) > 4:
            toas.append(float(e[2]))

psrpos=coord.SkyCoord(raj,decj,unit=(u.hourangle,u.degree))

t = atime.Time(toas,format='mjd')
moon=coord.get_moon(t)

oproduct = psrpos.cartesian.dot(moon.cartesian)

mass_ratio = 7.3476e22 / 5.972e24

path_delay = (oproduct/const.c).to(u.s).value*mass_ratio

header.ntoa=len(path_delay)
header.nrealisations=nreal
real = toasim.correction(header,path_delay,0,0,0,"")
with open(header.timfile_name+".noMoon","wb") as f:
    header.write(f)
    real.write(f)



lin_time = np.linspace(np.amin(toas),np.amax(toas),int(np.amax(toas)-np.min(toas))+1)
t = atime.Time(lin_time,format='mjd')
moon=coord.get_moon(t)

oproduct = psrpos.cartesian.dot(moon.cartesian)

mass_ratio = 7.3476e22 / 5.972e24

path_delay = (oproduct/const.c).to(u.s).value*mass_ratio


with open(header.timfile_name[:-4]+".asc","w") as f:
    for i in range(len(lin_time)):
        f.write("{:9.1f} {:g}\n".format(lin_time[i],path_delay[i]))


