#!/usr/bin/python
import sys
import matplotlib.pylab as plt
from math import sqrt,pi,degrees
import random
import ephem



def printhelp():
    print "TOA_characterise [options] tempo2.tim"
    print ""
    print "Characterises the given TOAs (tempo2 format tim file)."
    print "Can simulate future observations based on the existing data"
    print ""
    print "--simulate [ndays]    simulate ndays future TOAs."
    print "--shuffle             shuffle observations inside observing blocks."
    print "--seed [seed]         set the random seed for simulating"
    print "--select \"[flags]\"  only passes through given '-f' flags"
    print "--psr [psrname]       give the pulsar name for HA calculations"
    print "--plot                show some plots of the loaded TOAs"
    print ""
    print "e.g. typical usage:"
    print "TOA_characterise.py 0437-4715.tim --simulate 365 --shuffle --select '20CM_DFB3 20CM_APSR' --psr J0437-4715"
    print ""
    sys.exit(0)


def mjd2djd(mjd):
    jd=mjd+2400000.5
    return jd-2415020

def getHA(target,tel,mjd):
    tel.date=mjd2djd(mjd)
    ha_rad = float(tel.sidereal_time()) - float(target.ra)
    while ha_rad > pi:
        ha_rad -= 2*pi
    while ha_rad < -pi:
        ha_rad += 2*pi
    return ha_rad / pi * 12 # hours


class TOA:
    def __init__(self,t,e,f,F,o,xtraf):
        self.t=t
        self.e=e
        self.f=f
        self.F=F
        self.o=o.lower()
        self.siml=list()
        self.ha=0
        self.xtraf=xtraf

    def printxtraf(self):
        ret=""
        for k,v in self.xtraf:
            ret+=" %s %s"%(k,v)
        return ret


observatories = dict()

parkes = ephem.Observer()
parkes.long = '148.263333'
parkes.lat = '-32.999962'
parkes.elevation = 414.8
parkes.horizon = '30.25'


gbt = ephem.Observer()
gbt.long = '-79.839495'
gbt.lat = '38.432783'
gbt.elevation = 826.9
gbt.horizon = '10.0'

ao = ephem.Observer()
ao.long = '-66.752311'
ao.lat = '18.343853'
ao.elevation = 504.2
ao.horizon = '72.0'

observatories['7'] = parkes
observatories['parkes'] = parkes
observatories['arecibo'] = ao
observatories['gbt'] = gbt

if len(sys.argv) < 2:
    printhelp()

target = None
psrn=""
simulate=0
shuffle = 0
plot=0
file=open(sys.argv[1])
select=list()
random.seed()
good_has=list()
i=2
while i < len(sys.argv):
    if sys.argv[i]=="--select":
        i+=1
        select.extend(sys.argv[i].split())
    if sys.argv[i]=="--shuffle":
        shuffle=1
    if sys.argv[i]=="--simulate":
        i+=1
        simulate=int(sys.argv[i])
    if sys.argv[i]=="--plot":
        i+=1
        plot=sys.argv[i]
    if sys.argv[i]=="--seed":
        i+=1
        random.seed(int(sys.argv[i]))
    if sys.argv[i]=='--psr':
        i+=1
        psrn=sys.argv[i]
    if sys.argv[i]=='--help' or sys.argv[i]=='-h':
        printhelp()
    i+=1

if len(psrn) > 6:
    s=psrn.strip(" J")
    ra=s[0:2]+":"+s[2:4]
    dec=s[4:7]+":"+s[7:9]
    target=ephem.Equatorial(ra,dec,epoch='2000')
    print ra,target.ra,float(target.ra)
    print dec,target.dec,float(target.dec)

flags=dict()

total=0
sessions=list()
cur_ses=list()
sessions.append(cur_ses)
sess_seps=list()
prev_toa=None
total+=1
alltoas=list()
for line in file:
    if line[0] == "C":
        continue
    elems = line.strip().split()
    if elems[0] == "FORMAT":
        continue
    i = 5
    flag="??"
    xtraf=list()
    while i < len(elems):
        if elems[i] == "-f":
            flag=elems[i+1]
            i+=1
        if elems[i][0] == '-':
            xtraf.append((elems[i],elems[i+1]))
            i+=1
        i+=1
    if len(select) > 0 and not flag in select:
        continue
    if not flag in flags:
        flags[flag] = list()
    if len(elems) < 4:
        print >> sys.stderr, "ERROR: ",line

    #WARNING: This truncates the TOA! No use for timing
    #         but ok for simulations
    approx_time=float(elems[2])

    error = float(elems[3])
    freq = float(elems[1])
    obs = elems[4]
    toa = TOA(approx_time,error,freq,flag,obs,xtraf)
    if target!=None:
        toa.ha=getHA(target,observatories[toa.o],approx_time)
        good_has.append(toa.ha)
    if prev_toa!=None and abs(approx_time-prev_toa.t) < 0.005:
        # simlt
        toa.siml.append(prev_toa)
        prev_toa.siml.append(toa)
    if prev_toa!=None and abs(approx_time-prev_toa.t) > 5:
        if approx_time-prev_toa.t < 0:
            print "ERR",total
            print approx_time,error,freq,flag,obs,xtraf
        sess_seps.append(approx_time - cur_ses[0].t)
        cur_ses=list()
        sessions.append(cur_ses)
    prev_toa=toa
    flags[flag].append(toa)
    cur_ses.append(toa)
    total+=1
    alltoas.append(toa)

file.close()

print "There are %d TOAs in %d sessions"%(total,len(sessions))

keys=flags.keys()

fig = plt.figure()
if plot == "sepn" or plot == "errs":
    nx=int(sqrt(len(keys)))
    ny=int(len(keys)/nx + 0.99)
    p=0
    for flag in keys:
        if len(flags[flag]) < 5:
            continue
        print flag, len(flags[flag])
        toas=flags[flag]
        sepns=list()
        errs=list()
        freq=list()
        prev_t = -1
        for toa in toas:
            f=toa.f
            e=toa.e
            t=toa.t
            if prev_t == -1:
                prev_t=t
                continue
            t_sepn = t-prev_t
            errs.append(e)
            freq.append(f)
            sepns.append(t_sepn)
            prev_t=t
        p+=1
        if p==1:
            ax=fig.add_subplot(nx,ny,p)
            ax1=ax
        else:
            ax=fig.add_subplot(nx,ny,p)#, sharex=ax1)
        if plot=="sepn":
            print len(sepns)
            ax.hist(sepns,50)
        if plot=="errs":
            ax.hist(errs,50)
        ax.set_title(flag)
    print plot

if plot == "sessep":
    ax= fig.add_subplot(111)
    ax.hist(sess_seps,20)
if(plot!=0):
    plt.show()


print simulate
if simulate > 0:
    simtoas=list()
    start = sessions[-1][0].t
    print start
    last = start
    while last-start < simulate:
    #    next=random.choice(sessions)
    #    next_st=last+random.choice(sess_seps)
        r=random.randint(0,len(sess_seps)-1)
        next=sessions[r]
        next_st=last + sess_seps[r]
#        print "ST %f"%next_st
        range = next[-1].t - next[0].t
        for toa in next:
            toa.siml_shuf=None
        for toa in next:
            shuf=0
            if toa.siml_shuf != None:
                shuf=toa.siml_shuf
                # we have already done this for the siml obs
            else:
                if (shuffle):
                    shuf=random.uniform(next[0].t,next[-1].t)
                    shuf-=toa.t
                if target != None:
                #   check HA
                    goodHA=0
                    count = 0
                    while not goodHA:
                        count+=1
                        t=toa.t-next[0].t + next_st + shuf
                        ha=getHA(target,observatories[toa.o],t)
                        if shuffle and count < 5:
                            reroll=1
                            random.shuffle(good_has)
                            for test in good_has:
                                v=abs(test-ha)
                                if v > 12:
                                    v = 24 - v
                                if v < 0.2:
                                    goodHA=1
                                    reroll=0
                                    break
                                if v < 2:
                                    v= test-ha
                                    if v > 12:
                                        v = 24 - v
                                    if v < -12:
                                        v = -24 - v
                                    shuf += v/24.0
                                    reroll=0
                                    break
                            if reroll:
                                shuf=random.uniform(next[0].t,next[-1].t)
                                shuf-=toa.t
                        else:
                            goodHA = (abs(ha-toa.ha) < 0.2)
                            shuf+=0.1/24.0
            for tsim in toa.siml:
                tsim.siml_shuf = shuf
            t=toa.t-next[0].t + next_st + shuf
#            print "TT %f"%t
            newtoa=TOA(t,toa.e,toa.f,toa.F,toa.o,toa.xtraf)
            if target != None:
                newtoa.ha = getHA(target,observatories[toa.o],t)
            simtoas.append(newtoa)
        last=next_st
        print "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",
        print "%d"%(last-start),
    print "Done"
    i=0
    file = open("sim.tim","w")
    file.write("FORMAT 1\n")
    for toa in alltoas:
        file.write(" REAL_%d %f %f %f %s -f %s -sim R %s\n"%(i,toa.f,toa.t,toa.e,toa.o,toa.F,toa.printxtraf()))
        i+=1
    i=0
    for toa in simtoas:
        file.write(" SIML_%d %f %f %f %s -f %s -sim F %s\n"%(i,toa.f,toa.t,toa.e,toa.o,toa.F,toa.printxtraf()))
        i+=1
    file.close()
