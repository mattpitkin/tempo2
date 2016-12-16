#!/usr/bin/python
import toasim
import sys,random


nreal = 1
header = toasim.header()

header.parfile_name=sys.argv[1]
header.timfile_name=sys.argv[2]

par=open(sys.argv[1])
tim=open(sys.argv[2])

header.orig_parfile=par.read()
header.idealised_toas=tim.read()

reals=list()

nbad=10
sig=4.0

i=0
while i < len(sys.argv):
    if sys.argv[i] == "--nreal":
        i+=1
        nreal=int(sys.argv[i])
        continue
    if sys.argv[i] == "--nbad":
        i+=1
        nbad=int(sys.argv[i])
        continue

    if sys.argv[i] == "--sig":
        i+=1
        sig=float(sys.argv[i])
        continue


    i+=1




errs=list()

for line in header.idealised_toas.split("\n"):
    if line.startswith(" "):
        elems=line.strip().split()
        error=float(elems[3])*1e-6
        errs.append(error)


errs.sort()

errlim=errs[nbad]

print errlim

fff=open("bad.tim","w")

r=0
while r < nreal:
    offsets=list()
    ntoa=0
    for line in header.idealised_toas.split("\n"):
        if line.startswith(" "):
            elems=line.strip().split()
            error=float(elems[3])*1e-6
            toa=elems[2]
            off=0
            if error <= errlim:
                off+=error*sig*random.uniform(-1,1)
                print toa,error,off
                fff.write("%s -bad bad\n"%line)

            offsets.append(off)
            
            ntoa+=1
    r+=1
    print "\b\b\b\b\b\b\b\b",
    print "%d"%r,
    reals.append(toasim.correction(header,offsets,0,0,0,""))
header.ntoa=ntoa
header.nrealisations=nreal
print "\nWriting...."

fff.close()

header.orig_parfile=""
header.idealised_toas=""
file=open(header.timfile_name+".addBad","w")
header.write(file)
for real in reals:
    real.write(file)
file.close()

    
