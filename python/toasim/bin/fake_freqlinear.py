#!/usr/bin/python
import toasim
import sys,random,math


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
nreal = 1
while r < nreal:
    offsets=list()
    ntoa=0
    mult=random.uniform(-range/2.0,range/2.0)
    for line in header.idealised_toas.split("\n"):
        if line.startswith(" "):
            elems=line.strip().split()
            freq=float(elems[1])
            toa=elems[2]
            v=range/freq
            offsets.append(v)
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
file=open(header.timfile_name+".addFreqLinear","w")
header.write(file)
for real in reals:
    real.write(file)
file.close()

    
