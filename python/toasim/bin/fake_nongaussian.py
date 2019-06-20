#!/usr/bin/env python2
import numpy as np
from sys import argv,exit
from numpy import random
import toasim
import ast
def strbool(v):
    return v.lower() in ("yes", "true", "t", "1")

nreal = 1
default_distribution={'name':"normal", 'args':{}, 'scale':1}

outdata=list()

distributions=[]
default_normal=True
override_default_normal=False
global_rescale = False

header = toasim.header()

header.parfile_name=argv[1]
header.timfile_name=argv[2]

par=open(argv[1])
tim=open(argv[2])

header.orig_parfile=par.read()
header.idealised_toas=tim.read()

i=0
while i < len(argv):
    if argv[i] == "--nreal":
        i+=1
        nreal=int(argv[i])

    if argv[i] == "--default-normal":
        override_default_normal=True

    if argv[i] == "--rescale":
        global_rescale=True

    if argv[i] == "-d":
        default_normal=False
        i+=1
        name=argv[i]
        i+=1
        args=ast.literal_eval(argv[i])
        i+=1
        scale=float(argv[i])

        distributions.append({'name':name, 'args':args, 'scale':scale})
        continue
    i+=1


if default_normal or override_default_normal:
    distributions.append(default_distribution)


errs=list()
mjds=list()

for line in header.idealised_toas.split("\n"):
    if line.startswith(" "):
        elems=line.strip().split()
        error=float(elems[3])*1e-6
        errs.append(error)
        mjds.append(float(elems[2]))


errs = np.array(errs);
mjds = np.array(mjds);

real=0
while real < nreal:
    corn = np.zeros(len(errs))

    for distribution in distributions:
        r = getattr(random, distribution['name'])
        distribution_args = distribution['args']
        distribution_args['size'] = len(errs)

        vals = r(**distribution_args)
        if distribution['scale']:
            vals *= errs*distribution['scale']
        print distribution['name'],distribution['args']
        print "mean: %.3g"%np.mean(vals)
        print "sig : %.3g"%np.sqrt(np.var(vals))
        corn+=vals

    print ""
    print "==========="
    print "mean: %.3g"%np.mean(corn)
    print "sig : {0:.3g}   ({1:.3f} us)".format(np.sqrt(np.var(corn)),np.sqrt(np.var(corn))*1e6)


    if global_rescale:

        m = np.mean(corn)
        corn -=m
        s = np.mean(np.power(corn / errs,2))
        corn /= np.sqrt(s)
        corn+=m

        print ""
        print "Scale by ",(1.0/s)
        print "==========="
        print "mean: %.3g"%np.mean(corn)
        print "sig : {0:.3g}   ({1:.3f} us)".format(np.sqrt(np.var(corn)),np.sqrt(np.var(corn))*1e6)


    outdata.append(toasim.correction(header,corn,0,0,0,""))


    real+=1


header.ntoa=len(errs)
header.nrealisations=nreal
print "\nWriting...."

with open(header.timfile_name+".addNonGauss","w") as f:
    header.write(f)
    for r in outdata:
        r.write(f)

 
