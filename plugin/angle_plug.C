//  Copyright (C) 2006,2007,2008,2009,2010,2011 George Hobbs, Russell Edwards

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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "tempo2.h"
#include "T2toolkit.h"

using namespace std;

double psrangle(double centre_long,double centre_lat,double psr_long,double psr_lat);

void help() /* Display help */
{
}


extern "C" int graphicalInterface(int argc,char *argv[],pulsar *psr,int *npsr) 
{
    char parFile[MAX_PSR][MAX_FILELEN];
    char timFile[MAX_PSR][MAX_FILELEN];
    int i;

    *npsr = 0;  

    printf("Graphical Interface: angle\n");
    printf("Author:              M. Keith\n");
    printf("Version:             1.0\n");
    printf(" --- type 'h' for help information\n");

    /* Obtain all parameters from the command line */
    for (i=2;i<argc;i++)
    {
        if (strcmp(argv[i],"-f")==0)
        {
            strcpy(parFile[*npsr],argv[++i]); 
            strcpy(timFile[*npsr],argv[++i]);
            (*npsr)++;
        }

    }

    readParfile(psr,parFile,timFile,*npsr); /* Load the parameters       */
    FILE *o = fopen("angles","w");
    for (int p1=0; p1 < *npsr ; ++p1){
        for(int p2=0; p2 < p1 ; ++p2){
            double angle = (double)psrangle(psr[p1].param[param_raj].val[0],
                    psr[p1].param[param_decj].val[0],
                    psr[p2].param[param_raj].val[0],
                    psr[p2].param[param_decj].val[0]);
            fprintf(o,"%s\t%s\t%f\n",psr[p1].name,psr[p2].name,angle);
        }
    }

    return 0;
}

double psrangle(double centre_long,double centre_lat,double psr_long,double psr_lat)
{
    double dlon,dlat,a,c;
    double deg2rad = M_PI/180.0;

    /* Apply the Haversine formula */
    dlon = (psr_long - centre_long);
    dlat = (psr_lat  - centre_lat);
    a = pow(sin(dlat/2.0),2) + cos(centre_lat) * 
        cos(psr_lat)*pow(sin(dlon/2.0),2);
    if (a==1)
        c = M_PI/deg2rad;
    else
        c = 2.0 * atan2(sqrt(a),sqrt(1.0-a))/deg2rad;  
    return c;
}



