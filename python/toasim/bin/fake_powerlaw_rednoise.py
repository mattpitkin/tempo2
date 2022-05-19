#!/usr/bin/env python
import toasim
import numpy as np
from scipy.interpolate import interp1d
import sys
import argparse


parser = argparse.ArgumentParser("Fake pure power-law red noise")
parser.add_argument("-P",type=float, help="Power at fref")
parser.add_argument("-a",type=float, help="index")
parser.add_argument("--fref",type=float, help="fref")
parser.add_argument("--nreal",default=1,type=int)
parser.add_argument("--tn-like",action='store_true')
parser.add_argument("--plot",action='store_true')
parser.add_argument("--seed",type=int)
parser.add_argument("parfile")
parser.add_argument("timfile")


args=parser.parse_args()
print(args)


if args.plot:
    from matplotlib import pyplot as plt

if args.seed is None:
    seed = np.random.randint(low=0,high=np.iinfo(np.uint32).max,dtype=np.uint32)
else:
    seed = args.seed

nreal = args.nreal
header = toasim.header()


header.parfile_name=args.parfile
header.timfile_name=args.timfile

with open(args.parfile) as par, open(args.timfile) as tim:
    header.orig_parfile=par.read()
    header.idealised_toas=tim.read()

with open(header.timfile_name+".addRedNoise","wb") as outfile:
    toas=[]
    for line in header.idealised_toas.split("\n"):
        if line.startswith(" "):
            elems=line.strip().split()
            toa=float(elems[2])
            toas.append(toa)

    ntoa=len(toas)
    toas=np.array(toas)
    header.ntoa=ntoa
    header.nrealisations=nreal
    header.invocation=" ".join(sys.argv)

    header.seed = seed

    print("\nWriting....")
    header.write(outfile)


    itoas = np.argsort(toas)

    tspan=(np.amax(toas)-np.amin(toas))/365.25

    fmax = 365.25
    fmin = 1.0/200.0
    nharm = int(fmax/fmin)+1
    print("nharm={}".format(nharm))

    if args.tn_like:
        Pref = ((10**args.P)**2) / 12.0/np.pi/np.pi
    else:
        Pref=args.P

    freq = np.arange(1,nharm)*fmin
    P = Pref*(freq/args.fref)**-args.a

    rand = np.random.RandomState(seed)
    for ireal in range(nreal):
        print("ireal={}/{}".format(ireal,nreal))
        phases = rand.uniform(0,2*np.pi,size=len(freq))
        white = rand.normal(0,1,size=len(freq))
        offsets = np.zeros_like(toas)
        for i,f in enumerate(freq):
            if i%100 == 0:
                print("{:>10d}\r".format(i),end="")
            # The factor of 2 here is because we are doing A.sin(wt+phi)
            # rather than A.sin(wt) + B.cos(wt)
            A = np.sqrt(2*P[i]*fmin * (86400.0 * 365.25) ** 2)
            omega = 2*np.pi*f/365.25
            dat = A*white[i]*np.sin(omega*toas+phases[i])
            pp = np.poly1d(np.polyfit(toas,dat,2))
            dat -= pp(toas)
            offsets += dat
        print("")
        if args.plot:
            plt.plot(toas,offsets)
            plt.show()
        real = toasim.correction(header,offsets,0,0,0,"")
        real.write(outfile)
