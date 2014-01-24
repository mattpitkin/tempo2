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


/*
 * Set of routines necessary for simulating gravitational wave sources on pulsar
 * timing data
 *
 * Most of these routines are described in Hobbs, Jenet, Verbiest, Yardley, Manchester,
 * Lommen, Coles, Edwards & Shettigara (2009) MNRAS: "TEMPO2, a new pulsar timing 
 * package. III: Gravitational wave simulation
 *
 * These routines are used in the following plugins:
 *
 * GWbkgrd: simulation of a GW background
 *
 * This code has been developed by G. Hobbs and D. Yardley with help from F. Jenet and K. J. Lee.
 */

#include "T2toolkit.h"
#include "GWsim.h"

/* Have a structure to define a gravitational wave source */
typedef struct gwgeneralSrc
{
  long double theta_g;     /* Angle of source from "z"-direction (e.g. 90-declination)             */
  long double phi_g;       /* Azimuthal angle of source (e.g. right ascension)                     */
  long double omega_g;     /* Frequency of gravitational wave source (Hz) in observer's rest frame */   
  long double phi_polar_g; /* Polarization angle of the gravitational wave source w.r.t. the unit vector along te direction of increasing elevation*/
  long double phase_g;
  long double aplus_g;
  long double aplus_im_g; /* imagainary part of plus polarization */
  long double across_g;
  long double across_im_g; /*imag part of cross polarization */
  long double ast_g; /* real part of scalar transverse polarization */
  long double ast_im_g; /*imag part of scalar transverse polarization */
  long double asl_g; /* real part of scalar longitudinal polarization */
  long double asl_im_g; /*imag part of scalar longitudinal polarization */
  long double avx_g; /* real part of X vector longitudinal polarization */
  long double avx_im_g; /*imag part of X vector longitudinal polarization */
  long double avy_g; /* real part of Y vector longitudinal polarization */
  long double avy_im_g; /*imag part of Y vector longitudinal polarization */
  
  /*Important if the source is a supermassive blackhole binary. These are parameters of the binary. */
  
  long double phi_bin;    /* orientation of line of nodes on sky, recall line of nodes is intersection of plane of sky with plane of binary orbit. phi is
							 measured relative to the unit vector in the direction of increasing(?) right ascension */
  long double theta_bin;  /* orbital phase, assume time stationary */
  //long double chirp_mass; /* chirp mass of binary = (m_1+m_2)*[(m_1*m_2)/(m_1+m_2)^2]^(3/5) */
  long double inc_bin;    /* orbital inclination angle w.r.t. plane of the sky*/
  long double dist_bin;   /* proper distance to the binary system */

  /* Derived parameters                                                             */
  long double h[3][3];     /* The gravitational wave strain                         */
  long double h_im[3][3];  /* same as h but with im amplitudes. Ideally would have h = h + i*h_im */
  long double kg[3];       /* The gravitational wave unit vector                    */
}gwgeneralSrc;



/* Set up GW:                                                              
 * Sets up the vector pointing at the GW source                            
 * Sets up the strain matrix
 */

void setupgeneralGW(gwgeneralSrc *gw)
{
  long double deg2rad = M_PI/180.0;
  long double convert[3][3],trans[3][3];
  long double out[3][3];
  long double out_im[3][3];
  int i,j;
  long double st,ct,cp,sp,cpp,spp;

  st = sin(gw->theta_g);
  ct = cos(gw->theta_g);
  sp = sin(gw->phi_g);
  cp = cos(gw->phi_g);
  cpp = cos(gw->phi_polar_g);
  spp = sin(gw->phi_polar_g);

  gw->kg[0] = st*cp;
  gw->kg[1] = st*sp;
  gw->kg[2] = ct;

  /* Set up the GW strain: do not assume general relativity */
  gw->h[0][0] = gw->aplus_g+gw->ast_g;
  gw->h[0][1] = gw->across_g;
  gw->h[1][0] = gw->across_g;
  gw->h[1][1] = -gw->aplus_g+gw->ast_g+gw->asl_g;
  gw->h[2][0] = gw->avx_g;
  gw->h[2][1] = gw->avy_g;
  gw->h[2][2] = gw->asl_g;
  gw->h[0][2] = gw->avx_g;
  gw->h[1][2] = gw->avy_g;
  
  gw->h_im[0][0] = gw->aplus_im_g+gw->ast_im_g; // ALL SHOULD HAVE sqrt(-1) *  out the front
  gw->h_im[0][1] = gw->across_im_g; // only used for calcuateResidualGW and setupGW
  gw->h_im[1][0] = gw->across_im_g;
  gw->h_im[1][1] = -gw->aplus_im_g+gw->ast_im_g+gw->asl_im_g;
  gw->h_im[2][0] = gw->avx_im_g;
  gw->h_im[2][1] = gw->avy_im_g;
  gw->h_im[2][2] = gw->asl_im_g;
  gw->h_im[0][2] = gw->avx_im_g;
  gw->h_im[1][2] = gw->avy_im_g;
 
  /*printf("%10.8Le %10.8Le %10.8Le\n",gw->aplus_g,gw->ast_g,gw->asl_g);
  for (int jj=0;jj<3;jj++) {
                for (int kk=0;kk<3;kk++) {
                        printf("%10.8Le ",gw->h[jj][kk]);
                }
        }
        printf("\n");*/

  /* Now carry out a coordinate conversion (to convert to same coordinates 
     as the pulsar and GW source coordinates?)                                */

  /* These are from Euler's angles with pitch-roll-yaw convention */
  /* Also different definition of theta in MathWorld (mw) => sin(theta_mw) = -cos(theta_g),
     cos(theta_mw) = sin(theta_g) */
  /* Should check as there seems to be a minus sign difference in the top line
     with http://mathworld.wolfram.com/EulerAngles.html. Also the top line
     of this matrix is the bottom line of the wolfram matrix.  Probably doesn't
     matter, but should check */

  convert[0][0] = cp*ct*cpp-sp*spp; 
  convert[0][1] = sp*ct*cpp+cp*spp;
  convert[0][2] = -st*cpp;
       
  convert[1][0]= -cp*ct*spp-sp*cpp;
  convert[1][1]= -sp*ct*spp+cp*cpp;
  convert[1][2]= st*spp;
  
  convert[2][0]= cp*st;
  convert[2][1]= sp*st;
  convert[2][2]= ct;

  /* Calculate h.convert */
  matrixMult(gw->h,convert,out);
  matrixMult(gw->h_im,convert,out_im);

  /* Calculate transpose of 'convert' */
  for (i=0;i<3;i++)
    {
      for (j=0;j<3;j++)
	{
	  //	  printf("%g (%g, %g) ",(double)convert[j][i],(double)gw->h[j][i],(double)out[j][i]);
	  trans[i][j] = convert[j][i];
	}
      //      printf("\n");
    }
  /* Pre-multiply "out" by transpose */
  /* why multiply here -- we have already rotated it???? */
  matrixMult(trans,out,gw->h);
  matrixMult(trans,out_im,gw->h_im);
  //  matrixMult(trans,convert,out);
  /*  for (i=0;i<3;i++)
    {
      for (j=0;j<3;j++)
	{
	  printf("%g ",(double)gw->h[j][i]);
	}
      printf("\n");
      }*/

}

typedef struct gwgenSpec {
	double tensor_amp;
	double st_amp;
	double sl_amp;
	double vl_amp;
	double tensor_alpha;
	double st_alpha;
	double sl_alpha;
	double vl_alpha;
}gwgenSpec;

/* Produce a background of gravitational waves */
void GWgeneralbackground(gwgeneralSrc *gw,int *numberGW,long *idum,long double flo,long double fhi,
		  gwgenSpec gwAmps, int loglin)
{
  int j,k,l;
  k=0;
  for (j=0;j<4;j++) {
  	for (l=0;l<numberGW[j];l++)
    	{
      	gw[k].theta_g     = acos((TKranDev(idum)-0.5)*2);  
      	gw[k].phi_g       = TKranDev(idum)*2*M_PI;   
      	gw[k].phi_polar_g = 0.0; 
      	gw[k].phase_g     = TKranDev(idum)*2*M_PI;   
      	if (loglin==1)  /* Use equal sampling in log */
		{
	  	gw[k].omega_g  = 2*M_PI*exp(log(flo)+TKranDev(idum)*(log(fhi/flo))); 
		if (j == 0) {
	  		gw[k].aplus_g  = gwAmps.tensor_amp*pow(gw[k].omega_g/(2.0*M_PI), gwAmps.tensor_alpha)/sqrt((double)(numberGW[j])/log(fhi/flo))*TKgaussDev(idum);
	  		gw[k].across_g = gwAmps.tensor_amp*pow(gw[k].omega_g/(2.0*M_PI), gwAmps.tensor_alpha)/sqrt((double)(numberGW[j])/log(fhi/flo))*TKgaussDev(idum);
	  		gw[k].ast_g=gw[k].asl_g=gw[k].avx_g=gw[k].avy_g=0.;
	  		gw[k].aplus_im_g=gw[k].across_im_g=gw[k].ast_im_g=gw[k].asl_im_g=gw[k].avx_im_g=gw[k].avy_im_g=0.;
		} else if (j == 1) {
			gw[k].ast_g=gwAmps.st_amp*pow(gw[k].omega_g/(2.0*M_PI), gwAmps.st_alpha)/sqrt((double)(numberGW[j])/log(fhi/flo))*TKgaussDev(idum);
	  		gw[k].aplus_g=gw[k].across_g=gw[k].asl_g=gw[k].avx_g=gw[k].avy_g=0.;
	  		gw[k].aplus_im_g=gw[k].across_im_g=gw[k].ast_im_g=gw[k].asl_im_g=gw[k].avx_im_g=gw[k].avy_im_g=0.;
		} else if (j == 2) {
			gw[k].asl_g=gwAmps.sl_amp*pow(gw[k].omega_g/(2.0*M_PI), gwAmps.sl_alpha)/sqrt((double)(numberGW[j])/log(fhi/flo))*TKgaussDev(idum);
	  		gw[k].aplus_g=gw[k].across_g=gw[k].ast_g=gw[k].avx_g=gw[k].avy_g=0.;
	  		gw[k].aplus_im_g=gw[k].across_im_g=gw[k].ast_im_g=gw[k].asl_im_g=gw[k].avx_im_g=gw[k].avy_im_g=0.;
		} else if (j == 3) {
			gw[k].avx_g=gwAmps.vl_amp*pow(gw[k].omega_g/(2.0*M_PI), gwAmps.vl_alpha)/sqrt((double)(numberGW[j])/log(fhi/flo))*TKgaussDev(idum);
			gw[k].avy_g=gwAmps.vl_amp*pow(gw[k].omega_g/(2.0*M_PI), gwAmps.vl_alpha)/sqrt((double)(numberGW[j])/log(fhi/flo))*TKgaussDev(idum);
	  		gw[k].ast_g=gw[k].asl_g=gw[k].aplus_g=gw[k].across_g=0.;
	  		gw[k].aplus_im_g=gw[k].across_im_g=gw[k].ast_im_g=gw[k].asl_im_g=gw[k].avx_im_g=gw[k].avy_im_g=0.;
		}
	} 
      else
	{
	  gw[k].omega_g = 2*M_PI*(TKranDev(idum)*(fhi-flo)+flo); 
	  	if (j == 0) {
	  		gw[k].aplus_g  = gwAmps.tensor_amp*pow(gw[k].omega_g/(2.0*M_PI), gwAmps.tensor_alpha)*sqrt(1.0/(double)numberGW[j])*(sqrt((fhi-flo)*(2.0*M_PI)/gw[k].omega_g))*TKgaussDev(idum);
	  		gw[k].across_g  = gwAmps.tensor_amp*pow(gw[k].omega_g/(2.0*M_PI), gwAmps.tensor_alpha)*sqrt(1.0/(double)numberGW[j])*(sqrt((fhi-flo)*(2.0*M_PI)/gw[k].omega_g))*TKgaussDev(idum);
	  		gw[k].ast_g=gw[k].asl_g=gw[k].avx_g=gw[k].avy_g=0.;
	  		gw[k].aplus_im_g=gw[k].across_im_g=gw[k].ast_im_g=gw[k].asl_im_g=gw[k].avx_im_g=gw[k].avy_im_g=0.;
		} else if (j == 1) {
			gw[k].ast_g=gwAmps.st_amp*pow(gw[k].omega_g/(2.0*M_PI), gwAmps.st_alpha)*sqrt(1.0/(double)numberGW[j])*(sqrt((fhi-flo)*(2.0*M_PI)/gw[k].omega_g))*TKgaussDev(idum);
	  		gw[k].aplus_g=gw[k].across_g=gw[k].asl_g=gw[k].avx_g=gw[k].avy_g=0.;
	  		gw[k].aplus_im_g=gw[k].across_im_g=gw[k].ast_im_g=gw[k].asl_im_g=gw[k].avx_im_g=gw[k].avy_im_g=0.;
		} else if (j == 2) {
			gw[k].asl_g=gwAmps.sl_amp*pow(gw[k].omega_g/(2.0*M_PI), gwAmps.sl_alpha)*sqrt(1.0/(double)numberGW[j])*(sqrt((fhi-flo)*(2.0*M_PI)/gw[k].omega_g))*TKgaussDev(idum);
	  		gw[k].aplus_g=gw[k].across_g=gw[k].ast_g=gw[k].avx_g=gw[k].avy_g=0.;
	  		gw[k].aplus_im_g=gw[k].across_im_g=gw[k].ast_im_g=gw[k].asl_im_g=gw[k].avx_im_g=gw[k].avy_im_g=0.;
		} else if (j == 3) {
			gw[k].avx_g=gwAmps.vl_amp*pow(gw[k].omega_g/(2.0*M_PI), gwAmps.vl_alpha)*sqrt(1.0/(double)numberGW[j])*(sqrt((fhi-flo)*(2.0*M_PI)/gw[k].omega_g))*TKgaussDev(idum);
			gw[k].avy_g=gwAmps.vl_amp*pow(gw[k].omega_g/(2.0*M_PI), gwAmps.vl_alpha)*sqrt(1.0/(double)numberGW[j])*(sqrt((fhi-flo)*(2.0*M_PI)/gw[k].omega_g))*TKgaussDev(idum);
	  		gw[k].ast_g=gw[k].asl_g=gw[k].aplus_g=gw[k].across_g==0.;
	  		gw[k].aplus_im_g=gw[k].across_im_g=gw[k].ast_im_g=gw[k].asl_im_g=gw[k].avx_im_g=gw[k].avy_im_g=0.;
		}
	}
      //      printf("omega = %g\n",(double)(1.0/(gw[k].omega_g/2.0/M_PI)/86400.0));
	k++;
    }
  }
}

/* Calculate the pulsar timing residual induced by the gravitational wave */
long double calculateResidualgeneralGW(long double *kp,gwgeneralSrc *gw,long double time,long double dist)
{
  long double cosMu;
  long double psrVal1,psrVal2,psrVal1_im;
  long double earthVal1;
  long double tempVal[3];
  long double tempVal_im[3];
  long double geo;
  long double residual;
  long double deg2rad = M_PI/180.0L;
  long double VC = 299792456.200L;
  int i,k;

  for (i=0;i<3;i++)
    {
      tempVal[i] = 0.0;
      tempVal_im[i] = 0.0;
      for (k=0;k<3;k++)
	{
	  tempVal[i]    += gw->h[i][k]*kp[k];
	  tempVal_im[i] += gw->h_im[i][k]*kp[k];
	  //printf("%10.8Le ",gw->h[i][k]);
	}
            //printf("tempVal[%d] = %g ",i,(double)tempVal[i]);
    }
   //printf("\n");
   // printf("kp = %g %g %g\n",(double)kp[0],(double)kp[1],(double)kp[2]);

		
  cosMu = dotProduct(kp,gw->kg);   /* Angle between pulsar and GW source */
  //    printf(" cosmu = %g\n",(double)cosMu);
  if ((1+cosMu)!=0) 
  {
    psrVal1    = 0.5L/(1.0L+cosMu)*dotProduct(kp,tempVal); 
    psrVal1_im = 0.5L/(1.0L+cosMu)*dotProduct(kp,tempVal_im);
    //    printf("psrval1 = %g %g %g\n",(double)psrVal1,(double)(1+cosMu),(double)(1-cosMu));
  }
  else psrVal1 = 0.0L;
  
  geo      = psrVal1;
  //  printf("psrVal = %g %g\n",(double)psrVal1,(double)psrVal1_im);
  earthVal1 = psrVal1*sinl(gw->phase_g+time*gw->omega_g)+ 
    psrVal1_im*cosl(gw->phase_g+time*gw->omega_g); // sin -> cos
  if (dist==0) /* No pulsar term */
    psrVal2=0.0L;
  else
    psrVal2 = psrVal1*sinl(gw->phase_g-(1+cosMu)*dist/VC*gw->omega_g+time*gw->omega_g) + 
      psrVal1_im*cosl(gw->phase_g-(1+cosMu)*dist/VC*gw->omega_g+time*gw->omega_g); 

  //printf("%10.8Le %10.8Le %10.8Le\n",tempVal[0],tempVal[1],tempVal[2]);
  //printf("%10.8Le %10.8Le %10.8Le %10.8Le %10.8Le\n",psrVal1, psrVal1_im,earthVal1,psrVal2,gw->omega_g);
  residual = (earthVal1-psrVal2)/gw->omega_g; 
  //    printf("GWsim %g %g %g %g %g %g %g %g %g %g %g %g %g\n",(double)residual,(double)earthVal1,(double)psrVal2,(double)psrVal1,(double)cosMu,(double)psrVal1_im,(double)kp[0],(double)kp[1],(double)kp[2],(double)cosMu,(double)gw->h_im[0][0],(double)gw->h_im[0][1],(double)gw->h_im[1][0]);
  return residual;
}

int GWgeneralbackground_read(gwgeneralSrc *gw, FILE *file, int ireal){
	char key[13];
	int nreal;
	int ngw,id,igw,i;
	const unsigned int gwsize = 17*sizeof(long double);

	for (i=0;i<=ireal;i++){
		fread(&id,sizeof(int),1,file);
		printf("%d\n",id);
		if(id==ireal)break;
		fread(&ngw,sizeof(int),1,file);
		fseek(file,ngw*gwsize,SEEK_CUR);
		if(feof(file)){
			fprintf(stderr,"Could not read enough file to find realiation required\n");
			return -1;
		}
	}

	fread(&ngw,sizeof(int),1,file);

	printf("Reading %d GW sources from real %d (%d)\n",ngw,id,ireal);
	
	for (igw=0;igw<ngw;igw++){
		fread(&(gw[igw]),gwsize,1,file);
	}

	return ngw;
}



void GWgeneralbackground_write(gwgeneralSrc *gw, FILE *file,int ngw, int ireal){
	int igw;
	const unsigned int gwsize = 17*sizeof(long double);
	fwrite(&ireal,sizeof(int),1,file);
	fwrite(&ngw,sizeof(int),1,file);
	for (igw=0;igw<ngw;igw++){
		fwrite(&(gw[igw]),gwsize,1,file);
	}
	fflush(file);
}
