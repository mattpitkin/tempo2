#!/usr/bin/python
import toasim
import sys,random,math


nreal=1
header = toasim.header()

header.parfile_name=sys.argv[1]
header.timfile_name=sys.argv[2]

par=open(sys.argv[1])
tim=open(sys.argv[2])

header.orig_parfile=par.read()
header.idealised_toas=tim.read()

reals=list()

range=1e-4

i=0
while i < len(sys.argv):
    if sys.argv[i] == "--nreal":
        i+=1
        nreal=int(sys.argv[i])
        continue
    if sys.argv[i] == "--seed":
        i+=1
        seed = int(sys.argv[i])
        random.seed(seed)
        continue
    if sys.argv[i] == "--range":
        i+=1
        range = float(sys.argv[i])
        continue

    i+=1

offs=dict()


r=0
while r < nreal:
    offsets=list()
    ntoa=0
    for line in header.idealised_toas.split("\n"):
        if line.startswith(" "):
            elems=line.strip().split()
            freq=float(elems[1])
            toa=elems[2]
            band=3e8/(freq*1e6)*1e2 # cm waveband
            if band > 25 and band < 55:
                band=50
            elif band < 25 and band > 15:
                band=20
            elif band < 15 and band > 8:
                band=10
            else:
                band = int(band) - int(band)%10
            band=str(band)
            if not band in offs:
                offs[band] = random.uniform(-range/2.0,range/2.0)
                print "%s cm => %g"%(band,offs[band])
            
            offsets.append(offs[band])
            ntoa+=1
    r+=1
    print "\b\b\b\b\b\b\b\b",
    print "%d"%r,
    reals.append(toasim.correction(header,offsets,0,0,0,""))
header.ntoa=ntoa
header.nrealisations=nreal
print "\nWriting...."


header.orig_parfile=""
header.idealised_toas=""
file=open(header.timfile_name+".addFreqVar","w")
header.write(file)
for real in reals:
    real.write(file)
file.close()

    
