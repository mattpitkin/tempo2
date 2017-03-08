#!/usr/bin/python
import struct
from array import array
import numpy as np
import astropy.units as u
import astropy.coordinates as coord


class ToA:
    def __init__(self,toa,freq,error,site,flags=dict(),name="X"):
        e=toa.split(".")
        self.toa = (np.longdouble(e[0])+ np.longdouble("0."+e[1]))*u.day
        self.freq = freq.to(u.megahertz)
        self.site = site
        self.flags=flags
        self.error = error.to(u.microsecond)
        self.name=name

    @staticmethod
    def ld_to_string(toa):
        intpart = int(np.floor(toa))
        floatpart = toa - np.longdouble(intpart)
        str1="%d"%intpart
        str2="%.20f"%(floatpart)
        return "%s%s"%(str1,str2[1:])

    def __repr__(self):
        return "%s %f %f"%(ToA.ld_to_string(self.toa.value),self.error.value,self.freq.value)

    def __str__(self):
        return "%s %f %f"%(ToA.ld_to_string(self.toa.value),self.error.value,self.freq.value)

class correction:
    def __init__(self,header,offsets,a0,a1,a2,params):
        self.offsets=offsets
        self.a0=a0
        self.a1=a1
        self.a2=a2
        self.params=params
        self.header=header
    def write(self,stream):
        s=stream
        s.write("CORR")
        s.write(struct.pack("=ddd",self.a0,self.a1,self.a2))
        a=array('d',self.offsets)
        a.tofile(s)
        s.write(struct.pack("=%ds"%self.header.rparam_len,self.params))

class header:
    VERSION=2
    def __init__(self):
        self.version=header.VERSION
        self.writer=""
        self.timfile_name=""
        self.parfile_name=""
        self.invocation=""
        self.short_desc=""
        self.description=""
        self.idealised_toas=""
        self.orig_parfile=""
        self.gparam_desc=""
        self.gparam_vals=""
        self.rparam_desc=""
        self.rparam_len=0
        self.seed=-1
        self.ntoa=0
        self.nrealisations=0
        self.d_start=0
        self.d_offset=0

    def readParTim(self):
        with open(self.parfile_name) as f:
            self.orig_parfile = f.read()
        with open(self.timfile_name) as f:
            self.idealised_toas = f.read()


        toas=list()

        first=True
        t2format=False
        for line in self.idealised_toas.split("\n"):
            if line.startswith("C") or line.startswith("#"):
                continue
            line = line.strip()
            if first and line.startswith("FORMAT 1"):
                first=False
                t2format=True
            first=False
            elems=line.split()
            if len(elems) < 5:
                continue
            if t2format:
                flags=dict()
                i=5
                while (i+1) < len(elems):
                    flags[elems[i]] = elems[i+1]
                    i+=2
                toas.append(ToA(toa=elems[2],  \
                                freq=float(elems[1])*u.megahertz,   \
                                error=float(elems[3])*u.microsecond,\
                                site=elems[4],                      \
                                flags=flags,name=elems[0]))



        par = dict()
        for line in self.orig_parfile.split("\n"):
            line = line.strip()
            if line.startswith("#"):
                continue
            elems = line.split()
            if len(elems) < 2:
                continue
            key=elems[0]
            val=elems[1]

            par[key] = val

        if "RAJ" in par and "DECJ" in par:
            par['coord'] = coord.SkyCoord(ra=par['RAJ'],dec=par['DECJ'], unit=(u.hourangle,u.degree))
        return toas,par

    def write(self,stream):
        self.version=header.VERSION
        s=stream
        s.write("TOASIM")
        write_uint32(s,"VERS",self.version)
        writestr(s,"WRTR","toasim.py")
        writestr(s,"T_NM",self.timfile_name)
        writestr(s,"P_NM",self.parfile_name)
        writestr(s,"INVK",self.invocation)
        writestr(s,"SHRT",self.short_desc)
        writestr(s,"DESC",self.description)
        writestr(s,"TOAS",self.idealised_toas)
        writestr(s,"OPAR",self.orig_parfile)
        writestr(s,"GP_D",self.gparam_desc)
        writestr(s,"GP_V",self.gparam_vals)
        writestr(s,"RP_D",self.rparam_desc)
        write_uint32(s,"RP_L",self.rparam_len)
        write_int64(s,"SEED",self.seed)
        write_uint32(s,"NTOA",self.ntoa)
        write_uint32(s,"NREA",self.nrealisations)
        self.d_start=s.tell()+4*6
        self.d_offset= 4 + self.ntoa*8 + 3*8 + self.rparam_len
        write_uint32(s,"DSTT",self.d_start)
        write_uint32(s,"DOFF",self.d_offset)


def write_int64(stream,str,int):
    stream.write(str)
    stream.write(struct.pack("=Iq",8,int))
def write_int32(stream,str,int):
    stream.write(str)
    stream.write(struct.pack("=Ii",4,int))
def write_uint32(stream,str,int):
    stream.write(str)
    stream.write(struct.pack("=II",4,int))
def write_double(stream,str,double):
    stream.write(str)
    stream.write(struct.pack("=Id",8,double))
def writestr(stream,str,val):
    stream.write(str)
    stream.write(struct.pack("=I",len(val)))
    stream.write(val)
