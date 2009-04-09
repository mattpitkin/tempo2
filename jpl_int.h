//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russel Edwards

/*
*    This file is part of TEMPO2. 
* 
*    TEMPO2 is free software: you can redistribute it and/or modify 
*    it under the terms of the GNU General Public License as published by 
*    the Free Software Foundation, either version 3 of the License, or 
*    (at your option) any later version. 
*    TEMPO2 is distributed in the hope that it will be useful, 
*    but WITHOUT ANY WARRANTY; without even the implied warranty of 
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
*    GNU General Public License for more details. 
*    You should have received a copy of the GNU General Public License 
*    along with TEMPO2.  If not, see <http://www.gnu.org/licenses/>. 
*/

/*
*    If you use TEMPO2 then please acknowledge it by citing 
*    Hobbs, Edwards & Manchester (2006) MNRAS, Vol 369, Issue 2, 
*    pp. 655-672 (bibtex: 2006MNRAS.369..655H)
*    or Edwards, Hobbs & Manchester (2006) MNRAS, VOl 372, Issue 4,
*    pp. 1549-1574 (bibtex: 2006MNRAS.372.1549E) when discussing the
*    timing model.
*/

/* Define a new "long" which is 4 bytes in size -- Addition by George Hobbs */
typedef unsigned int JPLlong;

/* Right now,  DEs 403 and 405 have the maximum kernel size,  of 2036.    */
/* This value may need to be updated the next time JPL releases a new DE: */

#define MAX_KERNEL_SIZE 2036

/***** THERE IS NO NEED TO MODIFY THE REST OF THIS SOURCE (I hope) *********/


            /* A JPL binary ephemeris header contains five doubles and */
            /* (up to) 41 long integers,  so:                          */


#define JPL_HEADER_SIZE (5 * sizeof( double) + 41 * sizeof( JPLlong))



#pragma pack(1)

struct jpl_eph_data {
   double ephem_start, ephem_end, ephem_step;
   JPLlong ncon;
   double au;
   double emrat;
   JPLlong ipt[13][3];
   JPLlong ephemeris_version;
   JPLlong kernel_size, recsize, ncoeff;
   JPLlong swap_bytes;
   JPLlong curr_cache_loc;
   double pvsun[6];
   double *cache;
   void *iinfo;
   FILE *ifile;
   };

struct interpolation_info
   {
   double pc[18],vc[18], twot;
   int np, nv;
   };

#pragma pack()
