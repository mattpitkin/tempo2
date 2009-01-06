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
