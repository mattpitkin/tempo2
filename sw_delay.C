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

/* Solar wind model developed by You Xiaopeng, Bill Coles and George Hobbs */
/* This model is written up in You et al. (2007)                           */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "tempo2.h"

const double a1=1.155e11;
const double a2=32.3e11;
const double a3=3254e11;
const double a4=4.1e11;
const double a5=23.534e11;
const double a6=1.5e14;
const double a7=2.99e14;

using namespace std;

// float solar(pulsar *psr, int npsr, int p);

void mcl2(float eclon,float eclat,int iyr,int iday,float secs,float vel,float helat[],float crlon[],float rots[],float *helate,float *crlne,float *rote,float elong,float beta,float dlon,float delng,float cl,float zlonan,float e[],float ble,float delcrle,int lp);
float elsun2(int iyr,int iday,float secs,float gst,float sra,float sdec);
void calcRotN(float crlne, float rote, int *irot1, int *irot2, float *bcrlon);
double fast(double Rmin,double theta2,double theta1);
double slow(double Rmin,double theta2,double theta1);
void JulToGreg(long jd, long *year, long *month, long *day);

int calmjd(int y,int m,int d);
int calwd(int mjd);
int calwy(int mjd);
int calwn(int mjd);
int ocalmjd(int wy,int wn,int wd);
void mjd2date(int mjd,int *iyr,int *yy, int *mm, int *dd, int *iday);
void convertEcliptic(double raj,double decj,float *elong,float *elat);
 
float solarWindModel(pulsar psr,int iobs)
{
  int mjd;
  int yy,mm,dd,iyr,iday;
  float pi=3.1415927,theta,dp,inc,cl0,e[2],h[3],delcrle;
  int last,finished;
  float a[3],b[3],c[3],x[3][3],helat[36],crlon[36],rots[36],pe1[36],pe2[36],pe3[36];
  float eclon,eclat,vel=400.0;
  float secs;
  int lp;
  float elong,beta,dlon,delng,cl,zlonan,ble;
  float crlne,rote,helate;
  int i,j;
  float shiftangle=0.0;
  int istart;
  int irot1,irot2;
  float bcrlon;
  int introts[36];
  float crlon2[36]; /*using compare and give the minimum */
  int ijump; /*find where the longitude has about 360 jump*/
  ijump=1000;
  int npt=35; /*number of point with rots!=0*/
  float elsun,gst,sra,sdec;
  x[0][2]=0.0;
  x[1][2]=0.1262;
  x[2][2]=0.9920;
  if(psr.eclCoord==1)
    {
      eclat=psr.param[param_decj].val[0];
      eclon=psr.param[param_raj].val[0];
    }
  else if(psr.eclCoord==0)
    {
      convertEcliptic(psr.param[param_raj].val[0],psr.param[param_decj].val[0],&eclon,&eclat);
    }

  //  printf("%d, %30.20Lf\n",iobs, psr.obsn[iobs].sat);
  mjd=(int)psr.obsn[iobs].sat;
  mjd2date(mjd, &iyr, &yy, &mm, &dd, &iday);
  secs=(psr.obsn[iobs].sat-mjd)*86400;
  //_new3.tim      printf("%d %d %d %d %d\n",yy,mm,dd,iyr,iday);
  mcl2(eclon,eclat,iyr,iday,secs,vel,helat,crlon,rots,&helate,&crlne,&rote,elong,beta,dlon,delng,cl,zlonan,e,ble,delcrle,lp);    
  elsun=elsun2(iyr,iday,secs,gst,sra,sdec);
  //  printf("elsun=%f\n",elsun);
  calcRotN(crlne, rote, &irot1, &irot2, &bcrlon);
  for(i=0;i<36;i++)
    {
      introts[i]=(int)rots[i];
    }
  for(i=0;i<36;i++)crlon2[i]=crlon[i];
  for(i=1;i<35;i++)
    {
      if((rots[i]!=rots[i+1] || crlon[i]!=crlon[i+1] || helat[i]!=helat[i+1]) && introts[i]>1000 && introts[i+1]>1000)
	{
	  if(crlon2[i+1]-crlon2[i]>4.0)
	    {
	      for(j=0;j<=i;j++)
		{
		  crlon2[j]=crlon2[j]+2*pi;  
		}
	    }
	  else if(crlon2[i]-crlon2[i+1]>4.0)
	    {
	      crlon2[i+1]=crlon2[i+1]+2*pi;
	    }
	}
      else
	{
	  npt=i-1;
	  break;
	}
    }
  float midp,minp=1000,minp2=1000;
  float minrote=rots[1];
  int minnumber;
  for(i=1;i<=npt;i++)
    {
      if(rots[i]!=rots[i+1])
	{
	  if(crlon2[i]<minp)
	    {
	      minp=crlon2[i];
	      minrote=rots[i];
	    }
	  if(crlon[i]<minp2)minp2=crlon[i];
	}
    }
//  if(crlon2[16]<crlon2[18]) irot1=(int)minrote;
//  else irot1=(int)minrote-1;
  if(eclon>pi)
    {
      if(elsun>eclon-pi && elsun < eclon)
	irot1=(int)minrote-1;
      else
	irot1=(int)minrote;
    }
  else
    {
      if(elsun>eclon && elsun < eclon+pi)
	irot1=(int)minrote-1;
      else
	irot1=(int)minrote;
    }
  int rotn;
  shiftangle=minp-0.5;
  if (shiftangle < 0.0) 
    {
      shiftangle=shiftangle+2*pi;
      irot1=irot1+1;
    }
  irot2=irot1+1;
  FILE *fin,*fin1,*fin2;
  char str1[1000],str[100];
  float data[100];
  float lon[1000],lat[1000],m,d,x1,x2,nullx[1000],nully[1000];
  int n=0,type;
  float psrx[1000],psry[1000],xpos,px[2],py[2];
  shiftangle=shiftangle*180.0/pi;
  midp=midp*180.0/pi; 
  int npp=0;
  float *cont,contArray[5];
  float vals[200][30];
  int Ni = 49;
  int Nj = 29,pos;
  float tr[6];
  char name1[100],name2[100];
  char uplabel[100];
  int CTn;
  CTn=360/5-(int)(shiftangle/5);
  sprintf(name1,"%s/solarWindModel/CT%d.dat",getenv(TEMPO2_ENVIRON),irot1);
  sprintf(name2,"%s/solarWindModel/CT%d.dat",getenv(TEMPO2_ENVIRON),irot2);
  if(shiftangle/360+0.5>1)
    sprintf(uplabel,"PSR %s,  Rotation:%7.1f",psr.name,irot1+shiftangle/360-0.5);
  else
    sprintf(uplabel,"PSR %s,  Rotation:%7.1f",psr.name,irot1+shiftangle/360+0.5);
  
  fin1 = fopen(name1,"r");
  fin2 = fopen(name2,"r");
  for(j=1;j<=CTn;j++)
    {
      if (fscanf(fin1,"%s",&str1)==1)
	{
	  for (i=0;i<29;i++)
	    {
	      fscanf(fin1,"%f",&data[i]);
	      vals[n][i] = data[i];
	    }
	}
    }
  while (!feof(fin1))
    {
      if (fscanf(fin1,"%s",&str1)==1)
	{
	  for (i=0;i<29;i++)
	    {
	      fscanf(fin1,"%f",&data[i]);
	      vals[n][i] = data[i];
	    }
	  sscanf(str1+strlen(str1)-3,"%f",&lon[n]);
	  lon[n] = lon[n]-shiftangle;
	  if (lon[n]<0) lon[n]+=360.0;
	  if (lon[n]>360) lon[n]-=360.0;
	  /* Find neutral line */
	  if (data[0]>0) type=1;
	  else type=-1;
	  for (i=0;i<29;i++)
	    {
	      if (type==1 && data[i]<=0) break;
	      if (type==-1 && data[i]>0) break;
	    }
	  x1 = 70.0-5.0*(i-1);
	  x2 = 70.0-5.0*(i);
	  m = (data[i-1]-data[i])/(x1-x2);
	  d = data[i]-m*x2;
	  lat[n] = -d/m;
	  n++;
	}
    }
  fclose(fin1);
  if (fscanf(fin2,"%s",&str1)==1)
    {
      for (i=0;i<29;i++)
	{
	  fscanf(fin2,"%f",&data[i]);
	  vals[n][i] = data[i];
	}
    }
  while(n<72)
    {
      if (fscanf(fin2,"%s",&str1)==1)
	{
	  for (i=0;i<29;i++)
	    {
	      fscanf(fin2,"%f",&data[i]);
	      vals[n][i] = data[i];
	    }
	  sscanf(str1+strlen(str1)-3,"%f",&lon[n]);
	  lon[n] = lon[n]-shiftangle;
	  if (lon[n]<0) lon[n]+=360.0;
	  if (lon[n]>360) lon[n]-=360.0;
	  /* Find neutral line */
	  if (data[0]>0) type=1;
	  else type=-1;
	  for (i=0;i<29;i++)
	    {
	      if (type==1 && data[i]<=0) break;
	      if (type==-1 && data[i]>0) break;
	    }
	  x1 = 70.0-5.0*(i-1);
	  x2 = 70.0-5.0*(i);
	  m = (data[i-1]-data[i])/(x1-x2);
	  d = data[i]-m*x2;
	  lat[n] = -d/m;
	  n++;
	}
      
    }
  fclose(fin2);
  Ni = n;
  cont = (float *)malloc(Ni*Nj*sizeof(float));
  for (j=0;j<Nj;j++)
    {
      for (i=0;i<Ni;i++)
	{
	  pos = j*Ni+(int)((lon[i])/360.0*72.0);
	  cont[pos]=vals[i][j];
	}
    }
  
  for (i=0;i<360;i+=30)
    {
      xpos = i-shiftangle;
      if (xpos<0) xpos+=360.0;
      if (xpos>360) xpos-=360.0;
      sprintf(str,"%0d",i);
    }
  px[0]=0; px[1]=360;
  for (i=-90;i<90;i+=30)
    {
      py[0]=i;py[1]=i;
    }
  py[0]=-90; py[1]=90;
  for (i=0;i<360;i+=30)
    {
      px[0] = i-shiftangle;
      if (px[0]<0) px[0]+=360.0;      
      if (px[0]>360) px[0]-=360.0;      
      px[1]=px[0];
    }
  /* Now plot line-of-sight to pulsar */
  for(npp=0;npp<36;npp++)
    {
      psrx[npp]=crlon[npp];
      psrx[npp]*=180.0/pi;
      psrx[npp] = psrx[npp]-shiftangle;
      if (psrx[npp]<0) psrx[npp]+=360.0;
      if (psrx[npp]>360) psrx[npp]-=360.0;	  
      psry[npp]=helat[npp];
      psry[npp]*=180.0/pi;
    }
  int nsm=5;
  int nli=71*nsm+1;
  float lonli[nli],latli[nli];
  for(i=0;i<71;i++)
    {
      for(j=0;j<nsm;j++)
	{
	  lonli[i*nsm+j]=lon[i]-5.0/nsm*j;
	  if(lonli[i*nsm+j]<0) lonli[i*nsm+j]=lonli[i*nsm+j]+360;
	  latli[i*nsm+j]=lat[i]+(lat[i+1]-lat[i])/nsm*j;
	}
    }
  lonli[nli-1]=lon[71];
  latli[nli-1]=lat[71];
  float lonsm[356+nsm],latsm[356+nsm];
  float lon20p[nli],lat20p[nli];
  float lon20n[nli],lat20n[nli];
  float lonc,latc,slope;
  for(i=1;i<nli-1;i++)
    {
      lonc=(lonli[i]+lonli[i+1])/2;
      latc=(latli[i]+latli[i+1])/2;
      slope=(latli[i+1]-latli[i])/(lonli[i+1]-lonli[i]);
      slope=-1.0/slope;
      if (slope>=0.0)
	{
	  lon20p[i-1]=lonc+sqrt(400.0/(1.0+slope*slope));
	  lat20p[i-1]=latc+sqrt(400.0/(1.0+1.0/slope/slope));
	  lon20n[i-1]=lonc-sqrt(400.0/(1.0+slope*slope));
	  lat20n[i-1]=latc-sqrt(400.0/(1.0+1.0/slope/slope));
	}	   
      else	   
	{	   
	  lon20p[i-1]=lonc-sqrt(400.0/(1.0+slope*slope));
	  lat20p[i-1]=latc+sqrt(400.0/(1.0+1.0/slope/slope));
	  lon20n[i-1]=lonc+sqrt(400.0/(1.0+slope*slope));
	  lat20n[i-1]=latc-sqrt(400.0/(1.0+1.0/slope/slope));
	}
    }
  int l,k;
  float x20p[nli],y20p[nli];
  float x20n[nli],y20n[nli];
  float maxy20p,miny20n;
  float xbg,xed;
  if(lon20n[0]<lon20p[0])
    {xbg=lon20n[0];}
  else
    {xbg=lon20p[0];}
  if(lon20n[nli-3]>lon20p[nli-3])
    {xed=lon20n[nli-3];}
  else
    {xed=lon20p[nli-3];}
  
  for(l=0;l<nli;l++)
    {
      x20n[l]=(xed-xbg)/nli*l+xbg;
      x20p[l]=x20n[l];
      y20p[l]=-10000.0;
      y20n[l]=10000.0;
      for(k=0;k<nli-3;k++)
	{
	  if(x20p[l]<=lon20p[k] && x20p[l]>lon20p[k+1])
	    {
	      maxy20p=lat20p[k]+(lat20p[k+1]-lat20p[k])/(lon20p[k+1]-lon20p[k])*(x20p[l]-lon20p[k]);
	      if(maxy20p>y20p[l])y20p[l]=maxy20p;
	    }
	  if(x20n[l]<=lon20n[k] && x20n[l]>lon20n[k+1])
	    {
	      miny20n=lat20n[k]+(lat20n[k+1]-lat20n[k])/(lon20n[k+1]-lon20n[k])*(x20n[l]-lon20n[k]);
	      if(miny20n<y20n[l])y20n[l]=miny20n;
	    }
	}
    }
  

  tr[0] = -5.0;
  tr[1] = 5.0;
  tr[2] = 0.0;
  tr[3] = 75.0;
  tr[4] = 0.0;
  tr[5] = -5.0;
  contArray[0] = 0.0;
  contArray[1] = 0.0;
  
  /*calculate integration*/
  float thetakey[100];
  for(i=0;i<100;i++) thetakey[i]=0.0;
  int key1,key2,startkey;
  k=0;
  float tht[180];
  for(i=0;i<180;i++) tht[i]=0.0;
  float dm=0.0;
  float exdm=0.0;
  float dm2=0.0;
  for(i=0;i<npt*5;i++) /*degree of tracking line*/
    {
      tht[i]=(85-i*1.0)*pi/180; 
    }
  float psryint[180],psrxint[180];
  for(j=0;j<npt;j++)
    {
      for(i=0;i<5;i++) /*inteperation of tracking line*/
	{
	  psrxint[j*5+i]=psrx[j]+(psrx[j+1]-psrx[j])/5*i; 
	  psryint[j*5+i]=psry[j]+(psry[j+1]-psry[j])/5*i; 
	}  
    }

  float thtinte[2000];
  for(i=0;i<180;i++) thtinte[i]=0.0;
  thetakey[0]=tht[5];
  for(j=1;j<nli-1;j++)
    {
      if(psrx[1]<= (x20p[j-1]+x20p[j])/2 && psrx[1]> (x20p[j]+x20p[j+1])/2)
	{
	  if(psry[1]>=y20n[j] && psry[1]<=y20p[j])
	    startkey=0;
	  else 
	    startkey=1;
	}
    }
  for(i=6;i<npt*5;i++) /*start from 80 degree psrx[1] or psrxint[6]*/
    {
      for(j=1;j<nli-1;j++)
	{
	  if(psrxint[i]<= (x20p[j-1]+x20p[j])/2 && psrxint[i]> (x20p[j]+x20p[j+1])/2)
	    {
	      if(psryint[i]<y20n[j])key1=-1;
	      if(psryint[i]>y20p[j])key1=1;
	      if(psryint[i]>=y20n[j] && psryint[i]<=y20p[j])key1=0;
	    }
	  if(psrxint[i+1]<= (x20p[j-1]+x20p[j])/2 && psrxint[i+1]> (x20p[j]+x20p[j+1])/2)
	    {
	      if(psryint[i+1]<y20n[j])key2=-1;
	      if(psryint[i+1]>y20p[j])key2=1;
	      if(psryint[i+1]>=y20n[j] && psryint[i+1]<=y20p[j])key2=0;
	    }
	}
      if(key1!=key2)
	{
	  k++;
	  thetakey[k]=(tht[i-1]+tht[i])/2;
	}
    }
  float thetaearth;
  double ne_sw,ctheta,freqf,r,rsa[3],posp[3],pospos,delt;
  freqf = psr.obsn[iobs].freqSSB;
  delt = (psr.obsn[iobs].sat-psr.param[param_posepoch].val[0] + 
	  (getCorrectionTT(psr.obsn+iobs)+psr.obsn[iobs].correctionTT_TB)/SECDAY)/36525.0;
  
  for (j=0;j<3;j++)
    {
      rsa[j] = -psr.obsn[iobs].sun_ssb[j] + psr.obsn[iobs].earth_ssb[j] + 
	psr.obsn[iobs].observatory_earth[j];
      posp[j] = psr.posPulsar[j]+delt*psr.velPulsar[j];      
      pospos = sqrt(posp[0]*posp[0] + posp[1]*posp[1] + posp[2]*posp[2]);	    
    }
  for (j=0;j<3;j++)
    posp[j] /= pospos;
  
  r = sqrt(dotproduct(rsa,rsa));
  ctheta = dotproduct(posp,rsa)/r;
  //      printf("%d %d %d %f\n",yy,mm,dd,acos(ctheta)*180/pi);
  thetaearth=pi/2-acos(ctheta);
  //      if(thetaearth>thetakey[k])(k=k-1);
  thetakey[k+1]=thetaearth;
  float Rmin=(AU_DIST/SOLAR_RADIUS)*sin(thetaearth+pi/2);
  int N;
  //      for(j=0;j<=k+1;j++)printf("thetakey=%f\n",thetakey[j]*180/pi);
  float thout[10]; /* 89 degree to 81 dgree */
  float exdmfast=0.0,exdmslow=0.0;
  for(i=0;i<10;i++)
    {
      thout[i]=(float)(89-i)*pi/180;
    }
  for(i=0;i<9;i++)
    {
      exdmslow=exdmslow+slow(Rmin,thout[i],thout[i+1]);
      exdmfast=exdmfast+fast(Rmin,thout[i],thout[i+1]);
    }
  exdmslow=exdmslow+(a5/pow(Rmin,1.7)*pow(cos(89.0*pi/180),0.7)
		     +a6/pow(Rmin,5.0)*pow(cos(89.0*pi/180),4.0)
		     +a7/pow(Rmin,15.0)*pow(cos(89.0*pi/180),14.0))
    *pi/180/2;
  exdmfast=exdmfast+(a2/pow(Rmin,3.39)*pow(cos(89.0*pi/180),2.39)
		     +a3/pow(Rmin,15.25)*pow(cos(89.0*pi/180),14.25))
    *pi/180/2;
  if(thetaearth*180/pi<80)
    {
      if(k==0)
	{
	  N=(int)((thetakey[k]-thetakey[k+1])*180/pi);
	  for(j=0;j<=N;j++)
	    {
	  thtinte[j]=thetakey[k]-(thetakey[k]-thetakey[k+1])/N*j;
	    }
	  if(startkey==0)
	    {
	      dm2=dm2+a4/Rmin*(thetakey[k]-thetakey[k+1]);
	      for(j=0;j<N;j++)
		{
		  exdm=exdm+slow(Rmin,thtinte[j],thtinte[j+1]);
		}
	      exdm=exdm+exdmslow;
	    }
	  else
	    {
	      dm2=dm2+a1/Rmin*(thetakey[k]-thetakey[k+1]);
	      for(j=0;j<N;j++)
		{
		  exdm=exdm+fast(Rmin,thtinte[j],thtinte[j+1]);
		}
	      exdm=exdm+exdmfast;
	    }
	}
      else
	{
	  if(startkey==0)
	    {
	      for(i=0;i<=k;i=i+2)
		{
		  dm2=dm2+a4/Rmin*(thetakey[i]-thetakey[i+1]);
		  N=(int)((thetakey[i]-thetakey[i+1])*180/pi);
		  for(j=0;j<=N;j++)
		    {
		      thtinte[j]=thetakey[i]-(thetakey[i]-thetakey[i+1])/N*j;
		    }
		  for(j=0;j<N;j++)
		    {
		      exdm=exdm+slow(Rmin,thtinte[j],thtinte[j+1]);
		    }
		}
	      for(l=1;l<=k;l=l+2)
		{
		  dm2=dm2+a1/Rmin*(thetakey[l]-thetakey[l+1]);
		  N=(int)((thetakey[l]-thetakey[l+1])*180/pi);
		  for(j=0;j<=N;j++)
		    {
		      thtinte[j]=thetakey[l]-(thetakey[l]-thetakey[l+1])/N*j;
		    }
		  for(j=0;j<N;j++)
		    {
		      exdm=exdm+fast(Rmin,thtinte[j],thtinte[j+1]);
		    }
		}
	      exdm=exdm+exdmslow;
	    }
	  else
	    {
	      for(i=0;i<=k;i=i+2)
		{
		  dm2=dm2+a1/Rmin*(thetakey[i]-thetakey[i+1]);
		  N=(int)((thetakey[i]-thetakey[i+1])*180/pi);
		  for(j=0;j<=N;j++)
		    {
		      thtinte[j]=thetakey[i]-(thetakey[i]-thetakey[i+1])/N*j;
		    }
		  for(j=0;j<N;j++)
		    {
		      exdm=exdm+fast(Rmin,thtinte[j],thtinte[j+1]);
		    }
		}
	      for(l=1;l<=k;l=l+2)
		{
		  dm2=dm2+a4/Rmin*(thetakey[l]-thetakey[l+1]);
		  N=(int)((thetakey[l]-thetakey[l+1])*180/pi);
		  for(j=0;j<=N;j++)
		    {
		      thtinte[j]=thetakey[l]-(thetakey[l]-thetakey[l+1])/N*j;
		    }
		  for(j=0;j<N;j++)
		    {
		      exdm=exdm+slow(Rmin,thtinte[j],thtinte[j+1]);
		    }
		}
	      exdm=exdm+exdmfast; 
	    }
	}
    }
  else
    {
      for(j=1;j<nli-1;j++)
	{
	  if(psrx[npt]<= (x20p[j-1]+x20p[j])/2 && psrx[npt]> (x20p[j]+x20p[j+1])/2)
	    {
	      if(psry[npt]>=y20n[j] && psry[npt]<=y20p[j])
		startkey=0;
	      else 
		startkey=1;
	    }
	}
      N=(89-(int)(thetaearth*180/pi))+1;
      thtinte[0]=thetaearth;
      for(j=1;j<N;j++)
	{
	  thtinte[j]=(float)((int)(thetaearth*180/pi)+j)*pi/180;
	}
      if(startkey==0)
	{
	  dm2=dm2+a4/Rmin*(pi/2-thetaearth);
	  for(j=0;j<N-1;j++)
	    {
	      exdm=exdm+slow(Rmin,thtinte[j+1],thtinte[j]);
	    }
	  exdm=exdm+(a5/pow(Rmin,1.7)*pow(cos(89.0*pi/180),0.7)
		     +a6/pow(Rmin,5.0)*pow(cos(89.0*pi/180),4.0)
		     +a7/pow(Rmin,15.0)*pow(cos(89.0*pi/180),14.0))
		*pi/180/2;
	}
      else
	{
	  dm2=dm2+a1/Rmin*(pi/2-thetaearth);
	  for(j=0;j<N-1;j++)
	    {
	      exdm=exdm+fast(Rmin,thtinte[j+1],thtinte[j]);
	    }
	  exdm=exdm+(a2/pow(Rmin,3.39)*pow(cos(89.0*pi/180),2.39)
		     +a3/pow(Rmin,15.25)*pow(cos(89.0*pi/180),14.25))
	    *pi/180/2;
	}
    }
  dm2=dm2*1.0e-6*(SOLAR_RADIUS/PCM); 
  exdm=exdm*1.0e-6*(SOLAR_RADIUS/PCM); 
  dm=dm2+exdm;
  //      if(thetaearth*180/pi<-30)
  //	{
  //  printf("dm=%f\n",dm);
//  if(dm<0.01)
//    {
//      if(thetaearth*180/pi<-60)
//	{
//       printf("%d %d %d ",yy,mm,dd);
      //  if(iyr%4==1) printf("%d %f ",(iyr-104)*366+iday,dm);
//	  else printf("%d %f ",(iyr-104)*365+iday,dm);
//        printf(" %f\n",4.0*(pi/2-thetaearth)/206265/sin(pi/2-thetaearth));
//	}
//    }

//  if(acos(ctheta)*180/3.14159265<120)
//    dm =1.0e6*AU_DIST*AU_DIST/SPEED_LIGHT/DM_CONST_SI*psr.ne_sw*
//		  acos(ctheta)/r/sqrt(1.0-ctheta*ctheta)*DM_CONST*1e-12; 

  return(dm);
      //	}
    
}

void mcl2(float eclon,float eclat,int iyr,int iday,float secs,float vel,float helat[],float crlon[],float rots[],float *helate,float *crlne,float *rote,float elong,float beta,float dlon,float delng,float cl,float zlonan,float e[],float ble,float delcrle,int lp)
{
  float pi,theta,dp,inc,cl0,h[3];
  int last,finished;
  float gs,ra,de;
  float a[3],b[3],c[3],x[3][3],pe1[36],pe2[36],pe3[36];
  int i,j;
  for(i=0;i++;i<3)
    for(j=0;j++;j<3)
      x[i][j]=0.0;
  x[0][2]=0.0;
  x[1][2]=0.1262;
  x[2][2]=0.9920;
  finished=0;
  last=0;
  float	rad=180.0/3.141593;
  float helatep,crlnep,rotep;
  pi = 3.1415927;
  zlonan=1.28573 + 7.89327 * (iyr + iday/365.+ 50. )/ 32400.;
  float zkl=cos(eclat);
  float ya=elsun2(iyr,iday,secs,gs,ra,de);
  float delta_PA=atan(-.12722*cos(ya-zlonan));
  //  printf("ya,zlonana,delta_PA %f %f %f\n",ya*rad,zlonan*rad,delta_PA*rad);
  cl=cos(ya-eclon)*zkl;
  dp = sqrt(1.-cl*cl);
  elong=atan2(dp,cl);
  lp=0;
  cl0=cl;
  inc = 5./rad;
  theta = pi/2 - inc;
  cl = cl0 + dp*tan(theta);
  while(finished==0)
    {


      float clan=cos(zlonan);
      float slan=sin(zlonan);
      x[0][0]=clan;
      x[0][1]=slan;
      x[1][0]=-0.9920049497*slan;
      x[2][0]=0.1261989691*slan;
      x[1][1]=0.9920049497*clan;
      x[2][1]=-0.1261989691*clan;
      float sl=sqrt((1.-cl0*cl0)+((cl0-cl)*(cl0-cl)));
      float snx=-cos(ya);
      float sny=-sin(ya);
      a[0]=cos(eclon)*zkl;
      a[1]=sin(eclon)*zkl;
      a[2]=sin(eclat);
      b[0]=a[0]*cl+snx;
      b[1]=a[1]*cl+sny;
      b[2]=a[2]*cl;
      pe1[lp]=b[0];
      pe2[lp]=b[1];
      pe3[lp]=b[2];
      for(i=0;i<3;i++)
	{
	  c[i]=0.;
	  for(j=0;j<3;j++)
	    {
	      c[i]=c[i]+x[i][j]*a[j];
	    }
	}
      beta=atan2(c[2],sqrt(c[0]*c[0]+c[1]*c[1]));
      for(i=0;i<3;i++)
	{
	  a[i]=0.;
	  for(j=0;j<3;j++)
	    {
	      a[i]=a[i]+x[i][j]*b[j];
	    }
	}
      helat[lp]=atan2(a[2],sqrt(a[0]*a[0]+a[1]*a[1]));
      b[0]=x[0][0]*snx+x[0][1]*sny;
      b[1]=x[1][0]*snx+x[1][1]*sny;
      b[2]=x[2][0]*snx+x[2][1]*sny;
      float b1b2=sqrt(b[0]*b[0]+b[1]*b[1]);
      e[0]=snx;
      e[1]=sny;
      h[0]=b[0];
      h[1]=b[1];
      h[2]=b[2];
      helatep=atan2(b[2],b1b2);
      float cep=((-b[0]*c[0]-b[1]*c[1])/b1b2)/sqrt(c[0]*c[0]+c[1]*c[1]);
      if (cep > 1.) cep=1.;
      delng=atan2(sqrt(1.-cep*cep),cep);
      float blp=atan2(a[1],a[0]);
      if(blp < 0.)blp=blp+6.283185;
      ble=atan2(b[1],b[0]);
      if(ble<0.)ble=ble+6.283185;
      float zlagtm=1.496e8*sl/3600./vel;
      int ileap=(iyr-69)/4;
      float zhr= 24.*((iyr-69)*365.+ileap+iday-341.0) + secs/3600. -14.77;
      float zle1=228.42-zhr/1.81835;   
      float zle2= 228.42 - zhr/1.692002 + ble*rad;
      crlnep=360.+fmod(zle2,360.);
      //      printf("ileap,zhr,zle1,zle2,crlne %d %f %f %f %f\n",ileap,zhr,zle1,zle2,crlne);
      if(crlnep>360.) crlnep = fmod(crlnep,360.);
      delcrle = zlagtm/1.692002;
      float rote1=1556.-zle1/360.;
      float rote2=1556.-zle2/360.;
      rotep = (int)(rote1) + (rote2 - (int)(rote2));
      //      printf("rote=%f\n",rotep);
      if (rote1-rotep>(4./360.)) rotep = rotep + 1.;
      if (rote1-rotep<(-4./360.)) rotep = rotep - 1.;
      dlon=blp-ble;
      float delcrl = zlagtm/1.692002 + dlon*rad;
      float zlpp= zle2 + delcrl;
      crlon[lp]= 360. + fmod(zlpp,360.);
      if(crlon[lp] > 360.)crlon[lp] = fmod(crlon[lp],360.);
      if(crlon[lp] > crlnep)
	{
	  rots[lp]=(int)(rotep)+(1.-crlon[lp]/360.);
	}
      else
	{
	  rots[lp]=(int)(rotep)-(crlon[lp]/360.);
	  if (crlon[lp]/360. < 0.01) rots[lp]=(int)(rotep)-0.01;
	}
      crlon[lp]=crlon[lp]/rad;
      crlnep=crlnep/rad;
      //      printf("n,theta,cl %d %f %f %f %f %f\n",lp,theta*rad,cl,helat[lp]*rad,crlon[lp]*rad,rots[lp]);     
      theta = theta - inc;
      cl = cl0 + dp*tan(theta);
      if (last==1) finished=1;
      if (cl<0.)
	{
	  cl = 0.;
	  last=1;
	}
      lp=lp+1;
    }
  for(i=lp+1;i<36;i++)
    {
      helat[i] = helat[lp];
      crlon[i] = crlon[lp];
      rots[i] = rots[lp];
    }
  *helate=helatep;
  *crlne=crlnep;
  *rote=rotep;
}

float elsun2(int iyr,int iday,float secs,float gst,float sra,float sdec)
{
  double dj,fday;
  float rad=57.29578;
  float elsun3;
  if(iyr<1 || iyr>199) return elsun3;
  fday= secs/(double)86400.0;
  int idd= 365*iyr + (iyr-1)/4 + iday;
  dj= idd + fday - (double)0.5;
  float t= dj/(double)36525.;
  float vl= fmod((double)279.696678+(double)0.9856473354*dj,(double)360.);
  gst= fmod((double)279.690983+(double)0.9856473354*dj+360.*fday+180.,(double)360.);
  float g= fmod((double)358.475845+(double)0.985600267*dj,(double)360.)/rad;
  float elsun= vl+(1.91946-0.004789*t)*sin(g)+0.020094*sin(2.*g);
  float obliq= (23.45229- 0.0130125*t)/rad;
  float slp= (elsun-0.005686)/rad;
  float sind= sin(obliq)*sin(slp);
  float cosd= sqrt(1.-sind*sind);
  sdec= rad*atan(sind/cosd);
  float cot= cos(obliq)/sin(obliq);
  sra= 180.-rad*atan2(sind/cosd*cot,-cos(slp)/cosd);
  elsun3= elsun/rad;
  sdec= sdec/rad;
  sra= sra/rad;
  return elsun3;
}

void calcRotN(float crlne, float rote, int *irot1, int *irot2, float *bcrlon)
{
  (*bcrlon) = (double)((int)((.5+crlne*180./3.1415927)) + 180.);

  if (crlne<=3.1415927) 
    {
      (*irot1) = (int)rote;
      (*irot2) = (int)rote+1;
    }
  else
    {
      (*irot1) = (int)rote-1;
      (*irot2) = (int)rote;
      if ((*bcrlon) > 360.) 
	(*bcrlon) = (*bcrlon) - 360.;
    }
}

double fast(double Rmin,double theta2,double theta1)
{
  double exfast;
  exfast=(a2/pow(Rmin,3.39)*(pow(cos(theta2),2.39)+pow(cos(theta1),2.39))
	  +a3/pow(Rmin,15.25)*(pow(cos(theta2),14.25)+pow(cos(theta1),14.25)))
          *fabs(theta2-theta1)/2;
  //  exfast=exfast;
  return exfast;
}

double slow(double Rmin,double theta2,double theta1)
{
  double exslow;
  exslow=(a5/pow(Rmin,1.7)*(pow(cos(theta2),0.7)+pow(cos(theta1),0.7))
	  +a6/pow(Rmin,5.0)*(pow(cos(theta2),4.0)+pow(cos(theta1),4.0))
	  +a7/pow(Rmin,15.0)*(pow(cos(theta2),14.0)+pow(cos(theta1),14.0)))
	  *fabs(theta2-theta1)/2;
  //  exslow=exslow*1.0e-6*(r/pc);
  return exslow;
}

void JulToGreg(long jd, long *year, long *month, long *day) /* JD to Date*/
{
     long l,n,i,j;

     l = jd + 68569L;
     n = (4L * l) / 146097L;
     l = l - (146097L * n + 3L) / 4L;
     i = (4000L * (l + 1L)) / 1461001L;
     l = l - (1461L * i) / 4L + 31L;
     j = (80L * l) / 2447L;
     *day = l - (2447L * j) / 80L;
     l = j / 11L;
     *month = j + 2L - (12L * l);
     *year = 100L * (n - 49L) + i + l;
     //     printf("%d %d %d %d\n",l,n,i,j);
     return; 
}

int calmjd(int y,int m,int d) /*Date to MJD y=year-1900*/
{ 
  int l; 
  l=0; 
  if(m==1||m==2) l=1; 
  return(14956+d+(int)(((float)y-(float)l)*365.25)+(int)(((float)m+1+l*12)*30.6001)); 
} 
int calwd(int mjd) /*MJD to day in a week*/
{ 
  return((int)((mjd+2)%7)+1); 
} 
int calwy(int mjd) /*MJD to day the number of week since 1900*/
{ 
  int w; 
  w=(int)(((((float)mjd)/7)-2144.64)); 
  return((int)(((float)w*28/1461)-0.0079)); 
} 
int calwn(int mjd) /* MJD to the number of the week due to ISO 8601*/
{ 
  int w,wy; 
  w=(int)(((((float)mjd)/7)-2144.64)); 
  wy=(int)(((float)w*28/1461)-0.0079); 
  return(w-(int)(((float)wy*1461/28)+0.41)); 
}

int ocalmjd(int wy,int wn,int wd) 
{ 
  return(15012+wd+7*(int)((float)wn+(((float)wy*1461/28)+0.41))); 
}

void mjd2date(int mjd,int *iyr,int *yy, int *mm, int *dd, int *iday) /*MJD to date, to iyr=year-1900, to iday (number of the day in a year) */
{
  int y1,m1,k; 
  y1=(int)(((float)mjd-15078.2)/365.25); 
  m1=(int)(((float)mjd-14956.1-(int)((float)y1*365.25))/30.6001); 
  *dd=mjd-14956-(int)((float)y1*365.25)-(int)((float)m1*30.6001);
  k=0;
  if( m1==14 || m1==15 ) k=1; 
  *mm=m1-1-k*12;
  int mk=m1-1-k*12;
  *iyr=y1+k;
  *yy=*iyr+1900;
  if((*yy%4==0 && *yy%100!=0) || (*yy%400==0))
    {
      if(*mm==1)*iday=    (*dd);
      else if(*mm==2)*iday=31 +(*dd);
      else if(*mm==3)*iday=60 +(*dd);
      else if(*mm==4)*iday=91 +(*dd);
      else if(*mm==5)*iday=121+(*dd);
      else if(*mm==6)*iday=152+(*dd);
      else if(*mm==7)*iday=182+(*dd);
      else if(*mm==8)*iday=213+(*dd);
      else if(*mm==9)*iday=244+(*dd);
      else if(*mm==10)*iday=274+(*dd);
      else if(*mm==11)*iday=305+(*dd);
      else if(*mm==12)*iday=335+(*dd);
      else printf("date error\n");
    }
  else
    {
      if(*mm==1)*iday=    (*dd);
      else if(*mm==2)*iday=31 +(*dd);
      else if(*mm==3)*iday=59 +(*dd);
      else if(*mm==4)*iday=90 +(*dd);
      else if(*mm==5)*iday=120+(*dd);
      else if(*mm==6)*iday=151+(*dd);
      else if(*mm==7)*iday=181+(*dd);
      else if(*mm==8)*iday=212+(*dd);
      else if(*mm==9)*iday=243+(*dd);
      else if(*mm==10)*iday=273+(*dd);
      else if(*mm==11)*iday=304+(*dd);
      else if(*mm==12)*iday=334+(*dd);
      else printf("date error\n");
    }
}


void convertEcliptic(double raj,double decj,float *elong,float *elat)
{
  double sinb,beta,x,y,lambdap,lambda;
  double deg2rad = M_PI/180.0;
  double epsilon = 23.439292*deg2rad;
  /*  double epsilon = 23.441884*deg2rad;*/

  sinb = sin(decj)*cos(epsilon)-cos(decj)*sin(epsilon)*sin(raj);
  beta = asin(sinb);
  y = sin(raj)*cos(epsilon)+tan(decj)*sin(epsilon);
  x = cos(raj);
  
  lambdap = atan2(y,x);
  if (lambdap<0) lambda=lambdap+2*M_PI;
  else lambda = lambdap;

  *elong = lambda;
  *elat  = beta;
}
