#!/usr/bin/env python


import sys
from astropy import coordinates as coord
from astropy import units as u
from astropy import time as atime
from astropy import constants as const
import numpy as np

from scipy import optimize

import toasim

#from matplotlib import pyplot as plt

import argparse
parser=argparse.ArgumentParser(description="Inject a pulsar signal into a filterbank file")
parser.add_argument('par', type=str, help="Input par file")
parser.add_argument('tim', type=str, help="Input tim file")

parser.add_argument('-m','--mass',type=float, help='planet mass (earth mass)',default=1.0)
parser.add_argument('-M','--psrmass',type=float, help='pulsar mass (solar mass)',default=1.4)
parser.add_argument('-p','--period',type=float, help='orbital period (day)',default=100.0)
parser.add_argument('-i','--inclination', type=float, help='orbital inclination (deg)',default=90.0)
parser.add_argument('-T','--t0', type=float, help='time of periastron (mjd)',default=55000.0)
parser.add_argument('-O','--omega', type=float, help='angle of periastron (deg)',default=0.0)
parser.add_argument('-e','--ecc', type=float, help='eccintricity',default=0.0)

args=parser.parse_args()

nreal = 1
header = toasim.header()

header.parfile_name=args.par
header.timfile_name=args.tim

with open(args.par) as par:
    header.orig_parfile=par.read()
with open(args.tim) as tim:
    header.idealised_toas=tim.read()

toas=[]
with open(header.timfile_name) as f:
    for line in f:
        if line.startswith("C") or line.startswith("#"):
            continue
        e=line.split()
        if len(e) > 4:
            toas.append(float(e[2]))


Omega_b = 2.0*np.pi/(u.day*args.period)
inc=args.inclination*u.degree

M1 = args.psrmass * u.M_sun
M2 = args.mass * u.M_earth
Mtot = M1+M2
Mr = M2**3 / Mtot**2
t0=args.t0*u.day
e=args.ecc
a1 = np.power(Mr*const.G/Omega_b**2,1.0/3.0).to(u.m)

asini = a1 * np.sin(inc)

om=coord.Angle(args.omega*u.deg)


def get_roemer(t):

    def ecc_anom(E,e,M):
        return (E-e*np.sin(E))-M

    mean_anom = coord.Angle((Omega_b * (t*u.day - t0)).decompose().value*u.rad)
    mean_anom.wrap_at(2*np.pi*u.rad,inplace=True)
    mean_anom = mean_anom.rad

    E=np.zeros_like(mean_anom)


    for i in range(len(E)):
        v1=ecc_anom(0,e,mean_anom[i])
        E1=0
        for trial in np.linspace(0,2*np.pi,8)[1:]:
            E2=trial
            v2 = ecc_anom(trial,e,mean_anom[i])
            if v1*v2 < 0:
                break
            v1=v2
            E1=E2

        sol = optimize.root_scalar(ecc_anom,args=(e,mean_anom[i]),bracket=[E1,E2])
        E[i] = sol.root


    roemer = (asini*(np.cos(E)-e)*np.sin(om) + asini*np.sin(E)*np.sqrt(1.0-e**2)*np.cos(om))/const.c
    return roemer


roemer = get_roemer(toas)

header.ntoa=len(roemer)
header.nrealisations=nreal
real = toasim.correction(header,roemer.to(u.s).value,0,0,0,"")
with open(header.timfile_name+".addPlanet","wb") as f:
    header.write(f)
    real.write(f)


t=np.linspace(np.amin(toas),np.amax(toas),int(np.amax(toas)-np.min(toas))+1)
roemer2 = get_roemer(t)
with open(header.timfile_name[:-4]+".asc","w") as f:
    for i in range(len(t)):
        f.write("{:9.1f} {:g}\n".format(t[i],roemer2[i].value))



with open(header.parfile_name+".planet","w") as f:
    f.write(header.orig_parfile)
    f.write("#Binary parameters injected by simulated planet\n")
    f.write("BINARY BT\n")
    f.write("PB {}\n".format(args.period))
    f.write("A1 {}\n".format((asini/const.c).to(u.s).value))
    f.write("T0 {}\n".format(args.t0))
    f.write("OM {}\n".format(om.deg))
    f.write("ECC {}\n".format(args.ecc))

with open(header.parfile_name+".orbit","w") as f:
    f.write(header.orig_parfile)
    f.write("#Binary parameters injected by simulated planet\n")
    f.write("BINARY BT\n")
    f.write("PB {}\n".format(args.period))
    f.write("A1 0\n")
    f.write("T0 {}\n".format(args.t0))
    f.write("OM {}\n".format(om.deg))
    f.write("ECC {}\n".format(args.ecc))

