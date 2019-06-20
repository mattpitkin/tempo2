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

jumpflag='-f'
large=0
range=1.0

i=0
while i < len(sys.argv):
    if sys.argv[i] == "--flag":
        i+=1
        jumpflag = (sys.argv[i])
    elif sys.argv[i] == "--seed":
        i+=1
        seed = int(sys.argv[i])
        random.seed(seed)
    elif sys.argv[i] == "--large":
        large=1
    elif sys.argv[i] == "--range":
        i+=1
        range=float(sys.argv[i])
    elif sys.argv[i] == "--nreal":
        i+=1
        nreal=int(sys.argv[i])
    i+=1




reals=list()

toas=list()
flags=dict()
for line in header.idealised_toas.split("\n"):
    if line.startswith(" "):
        elems=line.strip().split()
        error=float(elems[3])*1e-6
        i=0
        while i < len(elems):
            if elems[i] == jumpflag:
                flagval=elems[i+1]
                break
            i+=1
        if not flagval in flags:
            flags[flagval] = dict(sum=0,n=0)
        flags[flagval]['sum']+=error
        flags[flagval]['n']+=1
        toas.append(flagval)



r=0
while r < nreal:
    for flagval in flags:
        flag=flags[flagval]
        mean=flag['sum']/float(flag['n'])
        if large:
            flag['jump'] = random.uniform(-range/2.0,range/2.0)
        else:
            flag['jump'] = random.gauss(0,1.0*mean)


    offsets=list()
    ntoa=0
    for flag in toas:
        offsets.append(flags[flag]['jump'])
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
file=open(header.timfile_name+".addFlagJumps","w")
header.write(file)
for real in reals:
    real.write(file)
file.close()

    
