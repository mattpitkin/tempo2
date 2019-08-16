#!/usr/bin/python
from numpy.fft import fft,ifft
import random,math,sys
import toasim


sindex=-11.0/3.0
sf_cal_scale=200
sf_cal_value=math.pow(10,-7.6)
nreal=1
niter=512
isc=50
osc=5000

kDM = 4.148808e-3

header = toasim.header()

header.parfile_name=sys.argv[1]
header.timfile_name=sys.argv[2]

par=open(sys.argv[1])
tim=open(sys.argv[2])

i=0
while i < len(sys.argv):
    if sys.argv[i] == "--nreal":
        i+=1
        nreal = int(sys.argv[i])
        continue
    if sys.argv[i] == "--os":
        i+=1
        osc = float(sys.argv[i])
        continue
    if sys.argv[i] == "--is":
        i+=1
        isc = float(sys.argv[i])
        continue
    if sys.argv[i] == "--itr":
        i+=1
        niter = int(sys.argv[i])
        continue
    if sys.argv[i] == "--sf_cal_scale":
        i+=1
        sf_cal_scale = float(sys.argv[i])
        continue
    if sys.argv[i] == "--sf_cal_value":
        i+=1
        sf_cal_value = float(sys.argv[i])
        continue
    if sys.argv[i] == "--seed":
        i+=1
        seed = int(sys.argv[i])
        random.seed(seed)
        continue


    i+=1

header.orig_parfile=par.read()
header.idealised_toas=tim.read()
header.nrealisations=nreal

toadata=header.idealised_toas.split("\n")

toas=list()
for line in toadata:
    if line.startswith(" "):
        elems=line.strip().split()
        mjd=float(elems[2])
        freq=float(elems[1])
        toas.append([mjd,freq])
toadata=""
mjd0 = toas[0][0]
ndays=toas[-1][0] - toas[0][0]
ntoas=len(toas)
header.ntoa=ntoas


nsamp = int(2.0*ndays/isc)

xscale = float(ndays)/float(nsamp)

print "nsamp=%d xscale=%f"%(nsamp,xscale)

v=osc/xscale
minf=1.0/float(v*2)
v=isc/xscale
maxf=v
sf_cal_scale = int(sf_cal_scale/xscale + 0.5)

cut_f=minf*2

print sf_cal_scale

ofile=open(header.timfile_name+".addDMvar","w")

header.write(ofile)

print mjd0,ndays,xscale

dmfile=open("dmvar.asc","w")
dmrawf=open("dmraw.asc","w")
dmra2f=open("dmra2.asc","w")
specf=open("spec.asc","w")
sf_file=open("sf.asc","w")

r=0
while r < nreal:

    
    data=list()
    spec=list()
    specx=list()
    ddl=list()
    i=0
    while i < nsamp:
        data.append(0)
        i+=1
    i=0
    while i < ntoas:
        ddl.append(0)
        i+=1

    f=minf
    while f < maxf:
        spec.append(0)
        specx.append(f/xscale)
        f+=minf


    for n in range(0,niter):
        freq=(fmax-fmin)*(float(n)/float(niter))+fmin 
        bin=int(freq/minf)
        if bin > len(spec)-1:
            bin=len(spec)-1
        if bin < 0:
            bin=0
        spec[bin] += 1
        a = math.pow(freq,sindex)
        if freq < cut_f:
            a *= freq/cut_f
        phase = random.random()
        i=0
        while i < nsamp:
            data[i] += a*math.sin(2*math.pi*(i*freq + phase))
            i+=1
        i=0
        while i < ntoas:
            mjd=toas[i][0]
            t=(mjd-mjd0)/xscale
            ddl[i] += a*math.sin(2*math.pi*(t*freq + phase))
            i+=1

        if n%100 == 0:
            print "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",
            print "%07d/%07d"%(n,niter),
            sys.stdout.flush()

    print "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",
    print "%07d/%07d\n"%(n,niter),
    sys.stdout.flush()


    if r==0:
        i=0
        while i < nsamp:
            dmrawf.write("%f %g\n"%(xscale*i,data[i]))
            i+=1

        i=0
        f=minf
        while f < maxf:
            p = 1.0/specx[i]
            specf.write("%f %g %f\n"%(specx[i],spec[i],p))
            f+=minf
            i+=1

    #compute structure function
    i=0
    sf=0
    n=0
    while i < nsamp-sf_cal_scale:
        v=abs(data[sf_cal_scale+i]-data[i])
        sf+=v*v
        n+=1
        i+=1
    sf/=n
    scale=math.sqrt(sf_cal_value/sf)
    print sf,scale

    i=0
    sum=0
    while i < nsamp:
        data[i]*=scale
        sum+=data[i]
        i+=1
    mean=sum/nsamp

    i=0
    while i < nsamp:
        data[i]-=mean
        if r ==0:
            dmra2f.write("%f %g\n"%(mjd0+xscale*i,data[i]))
        i+=1

    i=0
    while i < ntoas:
        ddl[i]*=scale
        ddl[i]-=mean
        i+=1


    i=0
    sf=0
    while i < nsamp-sf_cal_scale:
        v=abs(data[i]-data[i+sf_cal_scale])
        sf+=v*v
        i+=1



    sf_test_scale=2
    while sf_test_scale < nsamp/1.5:
        i=0
        sf2=0
        n=0
        while i < int(nsamp-sf_test_scale):
            v=abs(data[int(i+sf_test_scale)]-data[i])
            sf2+=v*v
            i+=1
            n+=1
        sf2/=n
        sf_file.write("%f %g\n"%(sf_test_scale*xscale,sf2))
        sf_test_scale = sf_test_scale*1.2




    offsets=list()
    i=0
    t=0
    for toa in toas:
        day=toa[0]-mjd0
        freq=toa[1]/1000.0 # GHz
        while day > i*xscale and i < nsamp-1:
            i+=1
        dmoff=ddl[t]
        if r == 0:
            dmfile.write("%f %g\n"%(day+mjd0,dmoff))
        delay=kDM * dmoff * pow(freq,-2)
        offsets.append(delay)
        t+=1
    corr = toasim.correction(header,offsets,0,0,0,"")
    corr.write(ofile)

    r+=1
ofile.close()
dmfile.close()
dmrawf.close()
