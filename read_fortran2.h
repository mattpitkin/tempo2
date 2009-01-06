/* Routines to read a FORTAN binary file */

#include <stdio.h>
#include <string.h>

FILE *c_fileptr2;
int swapByte2;

void open_file2(char *fname,int *swap)
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
    swapByte2 = 1;
  else
    swapByte2 = 0;
  *swap = swapByte2;

  printf("Byte swapping = %d\n",swapByte2);
  /* Look for blank character in filename */
  for (i=0;fname[i];i++)
    {
      if (fname[i]==' ') 
      {
        fname[i]='\0';
        break;
      }
    }
  if (!(c_fileptr2 = fopen(fname,"rb")))
    {
      printf("C PROGRAM: Unable to open filename >>%s<<\n",fname);
      exit(1);
    }
}
void close_file2()
{
  fclose(c_fileptr2);
}

void read_character2(int len,char *str)
{
  int i;
  if (swapByte2==1)
    {
      for (i=len-1;i>=0;i--)
	str[i] = fgetc(c_fileptr2);
    }
  else
    {
      for (i=0;i<len;i++)
	str[i] = fgetc(c_fileptr2);
    }
}

int read_int2()
{
  int i;
  
  union Convert
  {
    char asChar[sizeof(int)];
    int asInt;
  } intcnv;

  if (swapByte2==1)
    {
      for (i=sizeof(int)-1;i>=0;i--)
	intcnv.asChar[i]=fgetc(c_fileptr2);
    }
  else
    {
      for (i=0;i<(int)sizeof(int);i++)
	{
	  intcnv.asChar[i]=fgetc(c_fileptr2);

	}
    }
  return intcnv.asInt;
}

float read_float2()
{
  int i;

  union Convert
  {
    char asChar[sizeof(float)];
    float asFloat;
  } intcnv;
  if (swapByte2==1)
    {
      for (i=sizeof(float)-1;i>=0;i--)
	intcnv.asChar[i]=fgetc(c_fileptr2);
    }
  else
    {
      for (i=0;i<(int)sizeof(float);i++)
	{
	  intcnv.asChar[i]=fgetc(c_fileptr2);
	}
    }
  return intcnv.asFloat;
}

double read_double2()
{
  int i;
  static int count=0;

  union Convert
  {
    char asChar[sizeof(double)];
    double asDouble;
  } dblcnv;

  if (swapByte2==1)
    {
      for (i=sizeof(double)-1;i>=0;i--)
	{
	  dblcnv.asChar[i]=fgetc(c_fileptr2);      
	  /*	  printf("a) Have read %c (%d)\n",dblcnv.asChar[i],dblcnv.asChar[i]); */
	}
    }
  else 
    {
      for (i=0;i<(int)sizeof(double);i++)
	{
	  dblcnv.asChar[i]=fgetc(c_fileptr2);
	  /*	  printf("b) Have read %c (%d)\n",dblcnv.asChar[i],dblcnv.asChar[i]);  */
	}
    }
  /*    printf("Ret %d\n",count); */
  count++;
  return dblcnv.asDouble;
}
int read_record_int2()
{
  int ival;
  fread(&ival,sizeof(int),1,c_fileptr2);
  return ival;
}
