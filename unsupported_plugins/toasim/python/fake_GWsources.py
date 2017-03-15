#!/usr/bin/python
import toasim
import argparse
from astropy import coordinates as coord
from astropy import units as u
from astropy import constants as const
import numpy as np

def fakeGWsources(args):
    # Set up the output file.
    header = toasim.header()
    header.parfile_name=args.par
    header.timfile_name=args.tim

    # Read the .par and .tim file
    toas,par = header.readParTim()

    header.nrealisations=1
    header.ntoa = len(toas)


    # Just a debug print the content
    print toas[0]
    print par['coord']

    # Where is the source?
    source_pos = coord.SkyCoord(args.source_ra,args.source_dec,unit=(u.hourangle,u.degree))
    print source_pos

    # Extract the GW signal parameters from the program arguments
    gwfreq = float(args.frequency)*u.hertz
    gwamp = float(args.amp) # in what units?
    epoch = np.longdouble(args.epoch)*u.day

    offsets = list() # This will be the perturbations from the GW signal

    polterm = 1 # The polariation term, how to set?
    posterm = 1 # The position term, how to set?

    # Loop through all the ToAs and compute the GW signal at that time.
    for toa in toas:
        time = (toa.toa - epoch).to(u.second)
        perturbation = gwamp * np.cos((time * gwfreq).decompose()*2*np.pi*u.radian) * polterm * posterm
        offsets.append(perturbation) # CHECK units
        print perturbation # Print this out to check

    # Output the simulation
    realisation = toasim.correction(header, offsets,0,0,0,"")
    with open(header.timfile_name+".addGWSingle","w") as f:
        header.write(f)
        realisation.write(f)



if __name__=="__main__":

    # Read arguments from the command line...
    parser = argparse.ArgumentParser(description="Simulate a single GW source.")
    parser.add_argument('--nreal','-N',default=1)
    parser.add_argument('par',help="the .par file")
    parser.add_argument('tim',help="the .tim file")
    parser.add_argument('--amp','-A',help="GW amplitude",required=True) # Units???
    parser.add_argument('--epoch','-E',help="GW epoch of zero phase",required=True) # how to define ???
    parser.add_argument('--frequency','-f',help="GW source frequency",required=True)
    parser.add_argument('--source-ra','-r',help="GW source RA (direction GW proagating FROM)",required=True)
    parser.add_argument('--source-dec','-d',help="GW source DEC",required=True) 
    parser.add_argument('--source-pol','-p',help="GW source polarisation angle",required=True) # how to define
    args=parser.parse_args()

    # Run the GW source simulator
    fakeGWsources(args)

