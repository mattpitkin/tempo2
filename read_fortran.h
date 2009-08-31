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

/* Routines to read a FORTAN binary file */

#include <stdio.h>
#include <string.h>

FILE *c_fileptr;
int swapByte;

int open_file(char *fname)
{
  int i;
  
  /* NOTE MUST PUT BACK TO 0 FOR SOLARIS!!!!! */
  union Convert
  {
    char asChar[sizeof(int)];
    int asInt;
  } intcnv;
  intcnv.asInt = 1;
  if (intcnv.asChar[0]==1)
    swapByte = 1;
  else
    swapByte = 0;
  /* Look for blank character in filename */
  for (i=0;fname[i];i++)
    {
      if (fname[i]==' ') 
      {
        fname[i]='\0';
        break;
      }
    }
  if (!(c_fileptr = fopen(fname,"rb")))
    {
      printf("C PROGRAM: Unable to open filename >>%s<<\n",fname);
      return 1;
    }
  return 0;
}
void close_file()
{
  fclose(c_fileptr);
}

void read_character(int len,char *str)
{
  int i;
  for (i=0;i<len;i++)
    str[i] = fgetc(c_fileptr);
}
char read_char()
{
  return fgetc(c_fileptr);
}

int read_int()
{
  int i;
  
  union Convert
  {
    char asChar[sizeof(int)];
    int asInt;
  } intcnv;

  if (swapByte==1)
    {
      for (i=sizeof(int)-1;i>=0;i--)
	intcnv.asChar[i]=fgetc(c_fileptr);
    }
  else
    {
      for (i=0;i<(int)sizeof(int);i++)
	intcnv.asChar[i]=fgetc(c_fileptr);
    }
  return intcnv.asInt;
}

float read_float()
{
  int i;

  union Convert
  {
    char asChar[sizeof(float)];
    float asFloat;
  } intcnv;
  if (swapByte==1)
    {
      for (i=sizeof(float)-1;i>=0;i--)
	intcnv.asChar[i]=fgetc(c_fileptr);
    }
  else
    {
      for (i=0;i<(int)sizeof(float);i++)
	intcnv.asChar[i]=fgetc(c_fileptr);
    }
  return intcnv.asFloat;
}

double read_double()
{
  int i;

  union Convert
  {
    char asChar[sizeof(double)];
    double asDouble;
  } dblcnv;

  if (swapByte==1)
    {
      for (i=sizeof(double)-1;i>=0;i--)
	dblcnv.asChar[i]=fgetc(c_fileptr);
    }
  else
    {
      for (i=0;i<(int)sizeof(double);i++)
	dblcnv.asChar[i]=fgetc(c_fileptr);
    }
  return dblcnv.asDouble;
}
int read_record_int()
{
  int ival;
  fread(&ival,sizeof(int),1,c_fileptr);
  return ival;
}
