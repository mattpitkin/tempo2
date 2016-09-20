#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
//  Copyright (C) 2006,2007,2008,2009, George Hobbs, Russell Edwards

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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tempo2.h"
#include "read_fortran2.h"

double date2mjd(int idat);

void makechars(char* raw, char cbuf[25][9]);
void swap4(char*);
void swap8(char*);
void swap8(double* in){ swap8((char*)in);}

void readJBO_bat(char *fname,pulsar *psr,int p)
{
    int itype;
// UNUSED VARIABLE //     int endit=0;
    int nobs=0;
    char raw_record[25*8];
    char cbuf[25][9];
// UNUSED VARIABLE //     int *ints;
    double *doubles;
    //ints = (int*)raw_record;
    doubles = (double*)raw_record;

    const char *CVS_verNum = "$Id$";

    if (displayCVSversion == 1) CVSdisplayVersion("readJBO_bat.C","readJBO_bat()",CVS_verNum);

    printf("Attempting to read %s\n",fname);
    FILE *batfile = fopen(fname,"r");

    int irecord=0;
    int dmyi=0;
    int nsite=0;
    char sites[256][25];

    // skip 2 bytes

    while(!feof(batfile)) { 
        fread(&dmyi,1,4,batfile);
        fread(&itype,1,4,batfile);
        swap4((char*)(&itype));
        logmsg("Got record type = %d (%d)",itype,++irecord);
        fread(raw_record,25,8,batfile);
        fread(&dmyi,1,4,batfile);
        switch(itype){
            case 3:
                makechars(raw_record,cbuf);
                logmsg("Read BAT file created by %s, ver %s, user %s",cbuf[1],cbuf[2], cbuf[3]);
                break;
            case 0:
                makechars(raw_record,cbuf);
                logmsg("Source name: %s%s%s",cbuf[0],cbuf[1],cbuf[2]);
                break;
            case 1:
                makechars(raw_record,cbuf);
                swap8(doubles+3);
                swap8(doubles+4);
                swap8(doubles+5);
                logmsg("Observatory name: %s%s%s %lf %lf %lf",cbuf[0],cbuf[1],cbuf[2], doubles[3],doubles[4],doubles[5]);
                strcpy(sites[nsite],cbuf[0]);
                strcpy(sites[nsite]+8,cbuf[1]);
                strcpy(sites[nsite]+16,cbuf[2]);

                for (unsigned i=0; i < strlen(sites[nsite]); i++){
                    if(sites[nsite][i] == ' ')sites[nsite][i]='\0';
                }
                nsite++;
                break;
            case 2:
                makechars(raw_record,cbuf);
                for (int i=0; i < 25; i++){
                    if(i==13){
                        logmsg("%d %d %s",i,i+1,cbuf[i]);
                    } else {
                        swap8(doubles+i);
                        logmsg("%d %d %lf",i,i+1,doubles[i]);
                    }
                }
                char* site = sites[(int)doubles[4]-1];
                logmsg("site = %s",site);

                double diff = (double)(date2mjd((int)doubles[23])*86400.0+doubles[24]) -
                    (double)(date2mjd((int)doubles[0])*86400.0+doubles[1]);

                logmsg("sat-bat = %lg",diff);

                longdouble bat = (longdouble)(date2mjd((int)doubles[0])+doubles[1]/86400.0);

                double corr = /*doubles[14] + */doubles[15] + doubles[16] + doubles[17] + doubles[21];
                logmsg("corr =  = %lf",corr);
                psr->obsn[nobs].sat = bat;//bat + diff/86400.0;
                psr->obsn[nobs].bat = bat;
                logmsg("SAT = %.12lf",(double)psr->obsn[nobs].sat);
                logmsg("BAT = %.12lf",(double)psr->obsn[nobs].bat);

                //psr->obsn[nobs].sat = psr->obsn[nobs].bat;
                psr->obsn[nobs].toaErr = doubles[2]; 
                psr->obsn[nobs].origErr = doubles[2]; 
                psr->obsn[nobs].efac = 1;
                psr->obsn[nobs].equad = 0;
                psr->obsn[nobs].freq = doubles[3];
                psr->obsn[nobs].deleted = 0;
                psr->obsn[nobs].clockCorr = 0; 
                psr->obsn[nobs].delayCorr = 0;
                strcpy(psr->obsn[nobs].telID,site);
                strcpy(psr->obsn[nobs].fname,"UNKNOWN"); 
                ++nobs;

                break;
        }

    }
    /*
       for (i=0;i<25;i++)
       { buf.dbl[i] = read_double2(); }
       read_record_int2(); read_record_int2(); 
       if (itype==3)
       {
       if (swap==1)
       {
       }
       printf("Barycentric arrival times written by %8s ver %8s\n",buf.ch[1],buf.ch[2]); 
       }
       else if (itype==0) // New observation 
       {
       printf("Source: %s%s%s\n",buf.ch[1],buf.ch[2],buf.ch[3]);
       }
       else if (itype==1) // Observatory site 
       {
       printf("Observatory: %s%s%s\n",buf.ch[1],buf.ch[2],buf.ch[3]);
       }
       else if (itype==-1)
       endit=1;
       else if (itype==2) //BAT record 
       { 
       psr->obsn[nobs].bat = (longdouble)(date2mjd((int)buf.dbl[0])+buf.dbl[1]/86400.0);
       psr->obsn[nobs].sat = psr->obsn[nobs].bat;
       psr->obsn[nobs].toaErr = buf.dbl[2]; 
       psr->obsn[nobs].freq = buf.dbl[3];
       psr->obsn[nobs].deleted = 0;
       psr->obsn[nobs].clockCorr = 0; 
       psr->obsn[nobs].delayCorr = 0;
       strcpy(psr->obsn[nobs].telID,"BAT");  
       strcpy(psr->obsn[nobs].fname,"UNKNOWN"); 
       nobs++;
       }
       else
       {
       printf("Unknown itype = %d\n",itype);
       exit(1);
       }
       */
    psr->nobs=nobs;
    fclose(batfile);
}

double date2mjd(int idat)
{
    int monthd[13]={0,31,28,31,30,31,30,31,31,30,31,30,31};
    int iyr,imn,idy,i;
    double nd;

    iyr = (int)(idat / 10000) + 1900;
    imn = (int)fortran_mod(idat / 100,100);
    idy = (int)fortran_mod(idat,100);
    nd = (int)(iyr * 365.25);
    if (fortran_mod(iyr,4) == 0) 
    {
        monthd[2] = 29;
        nd = nd - 1;
    }
    else
        monthd[2] = 28;

    for (i=1;i<imn;i++)
        nd = nd + monthd[i];
    return (double)((nd + idy) - 678956);
}

void swap4(char* raw){
    char a = raw[0];
    char b = raw[1];
    char c = raw[2];
    char d = raw[3];
    raw[0]=d;
    raw[1]=c;
    raw[2]=b;
    raw[3]=a;
}
void swap8(char* raw){
    swap4(raw);
    swap4(raw+4);
    int* i = (int*)raw;
    int a = i[0];
    int b = i[1];
    i[1] = a;
    i[0] = b;
}

void makechars(char* raw, char cbuf[25][9]){
    int i;
    for(i=0; i < 25 ; i++){
        memcpy(cbuf[i],raw+i*8,8);
        cbuf[i][8]='\0';
    }
}
