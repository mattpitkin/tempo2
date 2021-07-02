#!/usr/bin/python
import toasim
import numpy as np
from scipy.interpolate import interp1d
import sys


nreal = 1
header = toasim.header()

header.parfile_name=sys.argv[1]
header.timfile_name=sys.argv[2]

par=open(sys.argv[1])
tim=open(sys.argv[2])

header.orig_parfile=par.read()
header.idealised_toas=tim.read()


x,y=np.loadtxt(sys.argv[3],unpack=True,usecols=(0,1))

Y = interp1d(x,y,kind='linear')

errs=list()

offsets=list()
ntoa=0
for line in header.idealised_toas.split("\n"):
    if line.startswith(" "):
        elems=line.strip().split()
        toa=float(elems[2])
        off=0
        print(toa)

        off = Y(toa)


        offsets.append(off)
        ntoa+=1



real = toasim.correction(header,offsets,0,0,0,"")
header.ntoa=ntoa
header.nrealisations=nreal
print("\nWriting....")

header.orig_parfile=""
header.idealised_toas=""
file=open(header.timfile_name+".addSig","w")
header.write(file)
real.write(file)
file.close()

    
