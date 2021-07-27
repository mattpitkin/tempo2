#!/usr/bin/env python
import toasim
import numpy as np
from scipy.interpolate import interp1d
import sys
import argparse


parser = argparse.ArgumentParser("Fake nudot variations")
parser.add_argument("--alt-nudot-factor", default=0.1, type=float)
parser.add_argument("--alt-time-N", nargs=2, type=float, help="Timescale of 'alt' mode given by a normal distribition given by mu/sigma in days")
parser.add_argument("--std-time-N", nargs=2, type=float, help="Timescale of 'std' mode given by a normal distribition given by mu/sigma in days")
parser.add_argument("--plot",action='store_true')
parser.add_argument("--subtract-quad","-s",action='store_true')
parser.add_argument("--nreal",default=1,type=int)
parser.add_argument("parfile")
parser.add_argument("timfile")



def integrate_phase(nudot,t0,t1, nu0):
    t0*=86400.0
    t1*=86400.0
    return 0.5*nudot*(t1**2-t0**2) + nu0*(t1-t0)

def integrate_nu(nudot,t0,t1):
    t0*=86400.0
    t1*=86400.0
    return nudot*(t1-t0)


args=parser.parse_args()
print(args)


if args.plot:
    from matplotlib import pyplot as plt

nreal = args.nreal
header = toasim.header()


header.parfile_name=args.parfile
header.timfile_name=args.timfile

with open(args.parfile) as par, open(args.timfile) as tim:
    header.orig_parfile=par.read()
    header.idealised_toas=tim.read()



with open(header.timfile_name+".addNudot","wb") as outfile:

    f1=None
    f0=None
    with open(args.parfile) as ff:
        for line in ff:
            e=line.split()
            if len(e) > 1:
                if e[0]=="F1":
                    f1=float(e[1])
                if e[0]=="F0":
                    f0=float(e[1])

    if f1 is None:
        print("No F1 found in par file")
        print(header.orig_parfile)
        sys.exit(1)
    if f0 is None:
        print("No F0 found in par file")
        print(header.orig_parfile)
        sys.exit(1)



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
    print("\nWriting....")
    header.write(outfile)


    itoas = np.argsort(toas)



    for ireal in range(nreal):
        print("ireal={}/{}".format(ireal,nreal))
        t = toas[itoas[0]] ## The time we have accumulated phase until
        t0=t
        accumulated_phase=0 ## the accumulated phase at t
        accumulated_nu=0 ## the accumulated change in nu at t.

        STD=(0,args.std_time_N[0],args.std_time_N[1])
        ALT=(args.alt_nudot_factor*f1, args.alt_time_N[0],args.alt_time_N[1])

        state=STD if np.random.uniform() < 0.5 else ALT
        next_state= ALT if state==STD else STD


        cur_nudot,mu,sigma = state
        next_switch = toas[itoas[0]] + np.random.normal(mu,sigma) * np.random.uniform() ## start a random way through the first interval

        phases=np.zeros_like(toas)
        state_lag = t
        other_state_lag = t

        for i in itoas:
            while toas[i] > next_switch:
                ## The next ToA occurs after a switch, so integrate phase up to the end of this switch.
                #accumulated_phase += integrate_phase(nudot=cur_nudot, t0=t-state_lag, t1=next_switch-state_lag,nu0=accumulated_nu)
                #accumulated_nu    += integrate_nu(nudot=cur_nudot, t0=t-state_lag, t1=next_switch-state_lag)

                accumulated_phase += integrate_phase(nudot=cur_nudot, t0=0, t1=next_switch-t,nu0=accumulated_nu)
                accumulated_nu    += integrate_nu(nudot=cur_nudot, t0=0, t1=next_switch-t)
                other_state_lag += next_switch - t
                t = next_switch
                state,next_state = next_state,state # swap state and next state
                cur_nudot,mu,sigma = state
                other_state_lag,state_lag = state_lag,other_state_lag
                next_switch = np.random.normal(mu,sigma) + t
            # Now integrate to the ToA
            #accumulated_phase += integrate_phase(nudot=cur_nudot, t0=t-state_lag, t1=toas[i]-state_lag,nu0=accumulated_nu)
            #accumulated_nu    += integrate_nu(nudot=cur_nudot, t0=t-state_lag, t1=toas[i]-state_lag)

            accumulated_phase += integrate_phase(nudot=cur_nudot, t0=0, t1=toas[i]-t,nu0=accumulated_nu)
            accumulated_nu    += integrate_nu(nudot=cur_nudot, t0=0, t1=toas[i]-t)
            other_state_lag   += toas[i]-t
            t = toas[i]
            phases[i] = accumulated_phase
        if args.subtract_quad:
            ## fit and remove quadratic
            pp = np.poly1d(np.polyfit(toas,phases,2))
            phases -= pp(toas)




        if args.plot:
            plt.subplot(311)
            plt.plot(toas[itoas],phases[itoas],ls=':',marker='x')
            plt.subplot(312)
            d1 = np.diff(phases[itoas]) / np.diff(toas[itoas])
            plt.plot(toas[itoas][:-1],d1,ls=':',marker='x')
            plt.subplot(313)
            d2 = np.diff(d1) / np.diff(toas[itoas])[:-1]
            plt.plot(toas[itoas][:-2],d2,ls=':',marker='x')
            plt.show()

        offsets=phases/f0
        

        real = toasim.correction(header,offsets,0,0,0,"")
        real.write(outfile)
