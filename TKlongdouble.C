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

/* This file contains useful functions that are commonly used throughout the TEMPO2 */
/* software */

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "TKlongdouble.h"

/* print out a longdouble to a std::string */
std::string print_longdouble(const longdouble &ld)
{
    char buf[1024];
#ifdef LONGDOUBLE_IS_FLOAT128
    quadmath_snprintf(buf,1024, "%Qg", ld);
#endif

#ifdef LONGDOUBLE_IS_IEEE754
    ld_sprintf(buf, "%Lg", ld);
#endif

#ifdef LONGDOUBLE_IS_DDREAL
    ld.write(buf);
#endif
    return std::string(buf);

}

longdouble parse_longdouble(const char *str)
{
    longdouble ld;
#ifdef LONGDOUBLE_IS_FLOAT128
    ld = strtoflt128(str,NULL);
#endif

#ifdef LONGDOUBLE_IS_IEEE754
    sscanf(str, "%Lf", &ld);
#endif

#ifdef LONGDOUBLE_IS_DDREAL
    ld = str;
#endif
    return ld;
}
/* Long double support routines */
#ifdef LONGDOUBLE_IS_DDREAL
dd_real pow(const dd_real &a, const dd_real &b)
{ return exp(b*log(a)); }
// operator float(const dd_real &a) 
// {return (float)(double)a;}
#endif


#ifdef LONGDOUBLE_IS_FLOAT128
#include <stdarg.h>
#define BUFSIZE 4096
int ld_vsprintf(char *buf, const char *__format, va_list args){
    char fmt[BUFSIZE];
    char qfmt[BUFSIZE];
    va_list oargs;
    const char *c;
    const char *e;
    char* o=buf;
    c=__format;
    do {
        if (*c!='%'){
            *o=*c;
            ++o;
            continue;
        }

        e=c+1;
        do {
            if (*e=='%'){
                *o='%';
                ++o;
                ++c;
                break;
            }
            if (*e=='L' && (*(e+1)=='f' || *(e+1)=='g' || *(e+1) == 'e')){
                longdouble ld = va_arg(args,longdouble);
                size_t n=e-c;
                ++e;
                memcpy(qfmt,c,n);
                qfmt[n]='Q';
                qfmt[n+1]=*e;
                qfmt[n+2]='\0';
                quadmath_snprintf(o,BUFSIZE-(o-buf),qfmt, ld);
                while(*(++o)!='\0'){
                    continue;
                }
                c+=n+1;
                break;
            }
            if ( *e=='u' || *e=='d' || *e=='x' || *e=='X' || *e=='o' || *e=='i' || *e=='c'){ 
                size_t n=e-c;
                memcpy(qfmt,c,n+1);
                qfmt[n+1]='\0';

                if(*(e-1)=='l' && *(e-2)=='l'){
                    long long int str = va_arg(args,long long int);
                    sprintf(o,qfmt,str);
                } else {
                    int str = va_arg(args,int);
                    sprintf(o,qfmt,str);
                }
                --c;
                while(*(++o)!='\0'){
                    continue;
                }
                c+=n+1;

                break;

            }
            if ( *e=='f' || *e=='g' || *e=='e') {
                size_t n=e-c;
                memcpy(qfmt,c,n+1);
                qfmt[n+1]='\0';
                double str = va_arg(args,double);
                sprintf(o,qfmt,str);
                --c;
                while(*(++o)!='\0'){
                    continue;
                }
                c+=n+1;

                break;
            }
            if ( *e=='s' ) {
                size_t n=e-c;
                memcpy(qfmt,c,n+1);
                qfmt[n+1]='\0';
                const char* str = va_arg(args,const char*);
                sprintf(o,qfmt,str);
                --c;
                while(*(++o)!='\0'){
                    continue;
                }
                c+=n+1;

                break;
            }
        } while(*(++e) != '\0');
    } while(*(++c) != '\0');

    *o='\0';
    int ret_status = 0;

    va_end(args);
    va_end(oargs);
    return ret_status;
}

int ld_fprintf(FILE* __stream, const char *__format, ...) {
    char buf[BUFSIZE];
    va_list args;
    va_start(args,__format);
    int ret = ld_vsprintf(buf,__format,args);
    va_end(args);
    fprintf(__stream,"%s",buf);
    return ret;
}


int ld_printf(const char *__format, ...) {
    char buf[BUFSIZE];
    va_list args;
    va_start(args,__format);
    int ret = ld_vsprintf(buf,__format,args);
    va_end(args);
    printf("%s",buf);
    return ret;
}

int ld_sprintf(char *buf, const char *__format, ...){
    va_list args;
    va_start(args,__format);
    int ret = ld_vsprintf(buf,__format,args);
    va_end(args);
    return ret;
}


#undef BUFSIZE

#endif

