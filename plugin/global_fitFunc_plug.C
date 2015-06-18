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


/**************
 *
 * A NOTICE TO THE PUBLIC
 *
 * This plugin is now depricated and does not need to be complied ever.
 *
 * Global fitting now works directly within tempo2, just use the -global
 * option as normal without needing to specify a fitfunc.
 *
 * If you have checked out from CVS and are having build problems because
 * of this file - try redoing the ./bootstrap command.
 *
 * Michael Keith 2013. mkeith@pulsarastronomy.net
 */


#include <stdio.h>
#include <stdlib.h>
#include <tempo2.h>
#include <math.h>
#include "TKfit.h"
#include "constraints.h"

void globalFITfuncs(double x,double afunc[],int ma,pulsar *psr,int ipos);
int gnpsr;

extern "C" int pluginFitFunc(pulsar *psr,int npsr,int writeModel) 
{
  int i,j,k,p;
  double tol = 1.0e-40;  /* Tolerence for singular value decomposition routine */
  int npol=1;
  int ip[MAX_OBSN];
  double *val,*error;
  double *x,*y,*sig,**covar;
  long double toffset;
  double chisq;
  int offset;
  int count=0;
  int weightfit=0;
 
  gnpsr = npsr;

  x = (double *)malloc(MAX_OBSN*sizeof(double));
  y = (double *)malloc(MAX_OBSN*sizeof(double));
  sig = (double *)malloc(MAX_OBSN*sizeof(double));
  printf("Update. About to undertake a global fit, number of pulsars = %d\n",npsr);

  // Form pre-fit residuals
  for (p=0;p<npsr;p++)
    {
      //      printf("Fitmode = %d\n",psr[p].fitMode);
      if (psr[p].fitMode==1) weightfit=1;
      if (weightfit==1 && psr[p].fitMode==0)  
	printf("WARNING: A weighted fit is being carried out, but PSR %s does not have MODE 1 in the parameter file\n",psr[p].name);
      for (i=0;i<psr[p].nobs;i++)
	{
	  if (psr[p].obsn[i].deleted!=0)
	    {
	      printf("ERROR: Please remove all deleted files in your TOA files before doing a global fit.  There is a problem with PSR %s\n",psr[p].name);
	      exit(1);	       
	    }
	  psr[p].obsn[i].prefitResidual = psr[p].obsn[i].residual;
	  x[count] = (double)(psr[p].obsn[i].bbat-psr[p].param[param_pepoch].val[0]);
	  y[count] = (double)(psr[p].obsn[i].prefitResidual);
	  if (psr[p].fitMode==0) sig[count] = 1.0;
	  else sig[count] = psr[p].obsn[i].toaErr*1.0e-6;
	  ip[count]=count;
	  //	  printf("Have %g %g %g %d\n",x[count],y[count],sig[count],ip[count]);
	  count++;
	  if (count >= MAX_OBSN)
	    {
	      printf("global_fitFunc: Must increase max obsn, must equal total number of observations for all pulsars\n");
	      exit(1);
	    }
	}
     // add constraints as extra pseudo observations
     // These point to non-existant observations after the nobs array
     // These are later caught by getParamDeriv.
     for (i=0; i < psr[p].nconstraints; i++){
       ip[count] = count; //psr[p].nobs+i;
	x[count]=0;
	y[count]=0;
	sig[count]=1e-12;
	count++;
      }

    }
  for (p=0;p<npsr;p++)
    psr[p].nFit=count;
  printf("Total number of points = %d\n",count);
  if (weightfit==1)printf("Doing a weighted fit (rescale = %d)\n",psr[0].rescaleErrChisq);
  else printf("Doing an unweighted fit\n");
  // Determine number of fit parameters
  npol=0;
  // Add global parameters
  for (i=0;i<MAX_PARAMS;i++)
    {
      for (k=0;k<psr[0].param[i].aSize;k++)
	{
	  if (psr[0].param[i].fitFlag[k]==2) {
	    {
	      printf("Found a global parameter: %s %d %d\n",psr[0].param[i].label[k],npol,psr[0].ifuncN);
	      npol++;
	    }
	  }
	}
    }
  if (psr[0].param[param_wave_om].fitFlag[0]==2)	  
    {
      npol+=psr[0].nWhite*2-1;
      printf("Adding to npol %d (%d)\n",psr[0].nWhite*2-1,npol); 
      for (p=1;p<npsr;p++)
	psr[p].param[param_waveepoch].val[0] = psr[0].param[param_waveepoch].val[0];
    }
  printf("Checking here %d\n",psr[0].param[param_quad_om].fitFlag[0]);
  printf("1] Adding %d global parameters\n",npol);
  if (psr[0].param[param_quad_om].fitFlag[0]==2)
      npol+=psr[0].nQuad*4-1;
  if (psr[0].param[param_ifunc].fitFlag[0]==2)	  
    npol+=psr[0].ifuncN-1;
  if (psr[0].param[param_clk_offs].fitFlag[0]==2)	  
    npol+=psr[0].clkOffsN-2;
  printf("2] Adding %d global parameters\n",npol);
  if (psr[0].param[param_tel_dx].fitFlag[0]==2 && psr[0].param[param_tel_dx].val[0] < 2)	  
    npol+=(psr[0].nTelDX-1);
  else if (psr[0].param[param_tel_dx].fitFlag[0]==2 && psr[0].param[param_tel_dx].val[0] == 2)	  
    npol+=(psr[0].nTelDX-2);
  printf("3] Adding %d global parameters\n",npol);
  if (psr[0].param[param_tel_dz].fitFlag[0]==2 && psr[0].param[param_tel_dz].val[0] < 2)	  
    npol+=(psr[0].nTelDZ-1);
  else if (psr[0].param[param_tel_dz].fitFlag[0]==2 && psr[0].param[param_tel_dz].val[0] == 2)	  
    npol+=(psr[0].nTelDZ-2);
  printf("4] Adding %d global parameters\n",npol);
  if (psr[0].param[param_tel_dy].fitFlag[0]==2 && psr[0].param[param_tel_dy].val[0] < 2)	  
    npol+=(psr[0].nTelDY-1);
  else if (psr[0].param[param_tel_dy].fitFlag[0]==2 && psr[0].param[param_tel_dy].val[0] == 2)	  
    npol+=(psr[0].nTelDY-2);
  printf("5] Adding %d global parameters\n",npol);
  if (psr[0].param[param_quad_ifunc_p].fitFlag[0]==2)	  
    npol+=psr[0].quad_ifuncN_p-1;
  if (psr[0].param[param_quad_ifunc_c].fitFlag[0]==2)	  
    npol+=psr[0].quad_ifuncN_c-1;
  printf("6] Adding %d global parameters\n",npol);
  if (psr[0].param[param_gwsingle].fitFlag[0]==2)
    npol+=(4-1); 

  printf("7] Adding %d global parameters\n",npol);
  for (p=0;p<npsr;p++)
    {
      npol++; // For the offset
      //      printf("Adding fitting func for the pulsar offset %d (%d)\n",p,npol);
      for (i=0;i<MAX_PARAMS;i++)
	{
	  for (k=0;k<psr[p].param[i].aSize;k++)
	    {
	      if (psr[p].param[i].fitFlag[k]==1) {
		if (i!=param_start && i!=param_finish)
		  {
		    npol++;
		    //		    printf("Adding fitting func %d/%d for pulsar %d (%d)\n",i,k,p,npol);
		  }
	      }
	    }
	}
      // Should do the same for jumps
      /* Add extra parameters for jumps */
      for (i=1;i<=psr[p].nJumps;i++)
	{
	  if (psr[p].fitJump[i]==1)
	    {
	      npol++;
	      printf("Adding jump fit for pulsar %d (%d)\n",p,npol);
	    }
	}
      /* Add extra parameters for sinusoidal whitening */
      if (psr[p].param[param_wave_om].fitFlag[0]==1)
	npol+=psr[p].nWhite*2-1;
      if (psr[p].param[param_quad_om].fitFlag[0]==1)
	npol+=psr[p].nQuad*4-1;
      if (psr[p].param[param_ifunc].fitFlag[0]==1)
	npol+=psr[p].ifuncN-1;
      if (psr[p].param[param_clk_offs].fitFlag[0]==1)
	npol+=psr[p].clkOffsN-2;
      if (psr[p].param[param_tel_dx].fitFlag[0]==1 && psr[p].param[param_tel_dx].val[0] < 2)
	npol+=(psr[p].nTelDX-1);
      else if (psr[p].param[param_tel_dx].fitFlag[0]==1 && psr[p].param[param_tel_dx].val[0] == 2)
	npol+=(psr[p].nTelDX-2);

      if (psr[p].param[param_tel_dy].fitFlag[0]==1 && psr[p].param[param_tel_dy].val[0] < 2)
	npol+=(psr[p].nTelDY-1);
      else if (psr[p].param[param_tel_dy].fitFlag[0]==1 && psr[p].param[param_tel_dy].val[0] == 2)
	npol+=(psr[p].nTelDY-2);

      if (psr[p].param[param_tel_dz].fitFlag[0]==1 && psr[p].param[param_tel_dz].val[0] < 2)
	npol+=(psr[p].nTelDZ-1);
      else if (psr[p].param[param_tel_dz].fitFlag[0]==1 && psr[p].param[param_tel_dz].val[0] == 2)
	npol+=(psr[p].nTelDZ-2);

      if (psr[p].param[param_quad_ifunc_p].fitFlag[0]==1)
	npol+=psr[p].quad_ifuncN_p-1;
      if (psr[p].param[param_quad_ifunc_c].fitFlag[0]==1)
	npol+=psr[p].quad_ifuncN_c-1;
      if (psr[p].param[param_gwsingle].fitFlag[0]==1)
	npol+=(4-1);
    }
  if (!(val   = (double *)malloc(npol*sizeof(double))))
    {
      printf("Unable to allocate enough memory\n");
      exit(1);
    }
  if (!(error = (double *)malloc(npol*sizeof(double))))
    {
      printf("Unable to allocate enough memory\n");
      exit(1);
    }
  
  if (!(covar = (double **)malloc(npol*sizeof(double *))))
    {
      printf("Unable to allocate enough memory\n");
      exit(1);      
    }
  for (i=0;i<npol;i++)
    {
      if (!(covar[i] = (double *)malloc(npol*sizeof(double))))
	{
	  printf("Unable to allocate enough memory\n");
	  exit(1);      
	}
    }
  printf("global_fitFunc: Number of fit parameters = %d\n",npol);
  for (i=0;i<npol;i++)
    val[i]=0.0;
  printf("Doing the fit: npol = %d\n",npol);
  TKleastSquares_svd_psr(x,y,sig,count,val,error,npol,covar,&chisq,globalFITfuncs,weightfit,psr,tol,ip);
  for (p=0;p<npsr;p++)
    {
      psr[p].fitChisq= chisq;
      psr[p].fitNfree = count-npol;
    }
  printf("Done the fit\n");
  //  for (i=0;i<npol;i++)
  //    printf("covar diag = %g\n",covar[i][i]);
    printf("Here 2, chisq = %g, red chisq = %g\n",chisq,chisq/(double)(count-npol));
    //  if (npol > MAX_PARAMS)
    //    {
    //      printf("ERROR: nterms=%d > MAX_PARAMS=%d\n",npol,MAX_PARAMS);
    //      exit(1);       
    //    }
    for (i=0;i<npol;i++)
      {
	for (j=0;j<npol;j++)
	  psr[0].covar[i][j]=covar[i][j];
      }
    printf("Now updating the parameters\n");
  //    for (i=0;i<npol;i++)
  //      printf("val %d = %g %g\n",i,val[i],error[i]);
  // now update the parameters
  offset=0;
  // update global parameters
  for (i=0;i<MAX_PARAMS;i++)
    {
      for (k=0;k<psr[0].param[i].aSize;k++)
	{
	  if (psr[0].param[i].fitFlag[k] == 2)
	    {
	      if (i==param_wave_om)
		{		  
		  for (j=0;j<psr[0].nWhite;j++)
		    {
		      for (p=0;p<npsr;p++)
			{
			  psr[p].wave_cos[j] -= val[offset];
			  psr[p].wave_cos_err[j] = error[offset];
			}
		      offset++;
		      for (p=0;p<npsr;p++)
			{
			  psr[p].wave_sine[j] -= val[offset];
			  psr[p].wave_sine_err[j] = error[offset];
			}
		      offset++;		      
		    }
		}
	      else if (i==param_quad_om)
		{
		  for (j=0;j<psr[0].nQuad;j++)
		    {
		      for (p=0;p<npsr;p++)
			{
			  psr[p].quad_aplus_r[j]    -= val[offset];		      
			  psr[p].quad_aplus_i[j]    -= val[offset+1];		      
			  psr[p].quad_across_r[j]   -= val[offset+2];		      
			  psr[p].quad_across_i[j]   -= val[offset+3];		      
			  psr[p].quad_aplus_r_e[j]   = error[offset];		      
			  psr[p].quad_aplus_i_e[j]   = error[offset+1];		      
			  psr[p].quad_across_r_e[j]  = error[offset+2];		      
			  psr[p].quad_across_i_e[j]  = error[offset+3];		      
			}	      
		      offset+=4;
		    }
		}
	      else if (i==param_gwsingle)
		{
		  for (p=0;p<npsr;p++)
		    {
		      psr[p].gwsrc_aplus_r -= val[offset];		      
		      psr[p].gwsrc_across_r -= val[offset+1];		      
		      psr[p].gwsrc_aplus_r_e = error[offset];		      
		      psr[p].gwsrc_across_r_e = error[offset+1];		      
		      psr[p].gwsrc_aplus_i -= val[offset+2];		      
		      psr[p].gwsrc_across_i -= val[offset+3];		      
		      psr[p].gwsrc_aplus_i_e = error[offset+2];		      
		      psr[p].gwsrc_across_i_e = error[offset+3];		      
		    }
		  offset+=4;
		}
	      else if (i==param_ifunc)
		{
		  printf("Updating %d point\n",psr[0].ifuncN);
		  for (j=0;j<psr[0].ifuncN;j++)
		    {
		      printf("Updating %g\n",val[offset]);
		      for (p=0;p<npsr;p++)
			{
			  psr[p].ifuncV[j]-=val[offset];
			  psr[p].ifuncE[j]=error[offset];
			}
		      offset++;
		    }
		}
	      else if (i==param_clk_offs)
		{
		  for (j=0;j<psr[0].clkOffsN-1;j++)
		    {
		      for (p=0;p<npsr;p++)
			{
			  psr[p].clk_offsV[j]-=val[offset];
			  psr[p].clk_offsE[j]=error[offset];
			}
		      offset++;
		    }
		}
	      else if (i==param_tel_dx)
		{
		  if (psr[0].param[param_tel_dx].val[0]==-1)
		    {
		      printf("UPDATING IN HERE %g %g\n",(double)psr[0].telDX_v[0],(double)val[offset]);
		      for (p=0;p<npsr;p++)
			{
			  psr[p].telDX_v[0]-=val[offset];
			  psr[p].telDX_e[0]=error[offset];
			}
		      offset++;
		    }
		  else if (psr[0].param[param_tel_dx].val[0] < 2)
		    {
		      for (j=0;j<psr[0].nTelDX;j++)
			{
			  printf("Setting: %d %g\n",j,val[offset]);
			  for (p=0;p<npsr;p++)
			    {
			      psr[p].telDX_v[j]-=val[offset];
			      psr[p].telDX_e[j]=error[offset];
			    }
			  offset++;
			}
		    }
		  else
		    {
		      for (j=0;j<psr[0].nTelDX-1;j++)
			{
			  printf("Setting: %d %g\n",j,val[offset]);
			  for (p=0;p<npsr;p++)
			    {
			      psr[p].telDX_v[j]-=val[offset];
			      psr[p].telDX_e[j]=error[offset];
			    }
			  offset++;
			}
		    }
		}
	      else if (i==param_tel_dy)
		{
		  if (psr[0].param[param_tel_dy].val[0]==-1)
		    {
		      for (p=0;p<npsr;p++)
			{
			  psr[p].telDY_v[0]-=val[offset];
			  psr[p].telDY_e[0]=error[offset];
			}
		      offset++;
		    }
		  else if (psr[0].param[param_tel_dy].val[0] < 2)
		    {
		      for (j=0;j<psr[0].nTelDY;j++)
			{
			  for (p=0;p<npsr;p++)
			    {
			      psr[p].telDY_v[j]-=val[offset];
			      psr[p].telDY_e[j]=error[offset];
			    }
			  offset++;
			}
		    }
		  else
		    {
		      for (j=0;j<psr[0].nTelDY-1;j++)
			{
			  for (p=0;p<npsr;p++)
			    {
			      psr[p].telDY_v[j]-=val[offset];
			      psr[p].telDY_e[j]=error[offset];
			    }
			  offset++;
			}
		    }
		}
	      else if (i==param_tel_dz)
		{
		  if (psr[0].param[param_tel_dz].val[0]==-1)
		    {
		      for (p=0;p<npsr;p++)
			{
			  psr[p].telDZ_v[0]-=val[offset];
			  psr[p].telDZ_e[0]=error[offset];
			}
		      offset++;
		    }
		  else if (psr[0].param[param_tel_dz].val[0] < 2)
		    {
		      for (j=0;j<psr[0].nTelDZ;j++)
			{
			  for (p=0;p<npsr;p++)
			    {
			      psr[p].telDZ_v[j]-=val[offset];
			      psr[p].telDZ_e[j]=error[offset];
			    }
			  offset++;
			}
		    }
		  else
		    {
		      for (j=0;j<psr[0].nTelDZ-1;j++)
			{
			  for (p=0;p<npsr;p++)
			    {
			      psr[p].telDZ_v[j]-=val[offset];
			      psr[p].telDZ_e[j]=error[offset];
			    }
			  offset++;
			}
		    }
		}
	      else if (i==param_tel_x0)
		{
		  for (p=0;p<npsr;p++)
		    {
		      psr[p].param[param_tel_x0].val[0]-=val[offset];
		      psr[p].param[param_tel_x0].err[0]=error[offset];
		    }
		  offset++;

		}
	      else if (i==param_tel_y0)
		{
		  for (p=0;p<npsr;p++)
		    {
		      psr[p].param[param_tel_y0].val[0]-=val[offset];
		      psr[p].param[param_tel_y0].err[0]=error[offset];
		    }
		  offset++;

		}
	      else if (i==param_tel_z0)
		{
		  for (p=0;p<npsr;p++)
		    {
		      psr[p].param[param_tel_z0].val[0]-=val[offset];
		      psr[p].param[param_tel_z0].err[0]=error[offset];
		    }
		  offset++;

		}
	      else if (i==param_tel_vx)
		{
		  for (p=0;p<npsr;p++)
		    {
		      psr[p].param[param_tel_vx].val[0]-=val[offset];
		      psr[p].param[param_tel_vx].err[0]=error[offset];
		    }
		  offset++;

		}
	      else if (i==param_tel_vy)
		{
		  for (p=0;p<npsr;p++)
		    {
		      psr[p].param[param_tel_vy].val[0]-=val[offset];
		      psr[p].param[param_tel_vy].err[0]=error[offset];
		    }
		  offset++;

		}
	      else if (i==param_tel_vz)
		{
		  for (p=0;p<npsr;p++)
		    {
		      psr[p].param[param_tel_vz].val[0]-=val[offset];
		      psr[p].param[param_tel_vz].err[0]=error[offset];
		    }
		  offset++;

		}
	      else if (i==param_quad_ifunc_p)
		{
		  printf("Updating %d point\n",psr[0].quad_ifuncN_p);
		  for (j=0;j<psr[0].quad_ifuncN_p;j++)
		    {
		      printf("Updating %g\n",val[offset]);
		      for (p=0;p<npsr;p++)
			{
			  psr[p].quad_ifuncV_p[j]-=val[offset];
			  psr[p].quad_ifuncE_p[j]=error[offset];
			}
		      offset++;
		    }
		}
	      else if (i==param_quad_ifunc_c)
		{
		  printf("Updating %d point\n",psr[0].quad_ifuncN_c);
		  for (j=0;j<psr[0].quad_ifuncN_c;j++)
		    {
		      printf("Updating %g\n",val[offset]);
		      for (p=0;p<npsr;p++)
			{
			  psr[p].quad_ifuncV_c[j]-=val[offset];
			  psr[p].quad_ifuncE_c[j]=error[offset];
			}
		      offset++;
		    }
		}
	      else
		{
		  for (p=0;p<npsr;p++)
		    {
		      printf("Fitting (1): %d %d %d %d %g\n",i,k,p,offset,val[offset]);
		      if (i==param_telx || i==param_tely || i==param_telz || i==param_gwm_amp)
			psr[p].param[i].val[k] -= val[offset];
		      else
			psr[p].param[i].val[k] += val[offset];
		      psr[p].param[i].err[k] = error[offset];
		    }
		  offset++;
		}
	    }
	}
    }
  for (p=0;p<npsr;p++)
    {
      updateParameters(psr,p,val+offset,error+offset);
      for (i=0;i<MAX_PARAMS;i++)
	{
	  for (k=0;k<psr[p].param[i].aSize;k++)
	    {
	      if (psr[p].param[i].fitFlag[k]==1) {
		if (i!=param_start && i!=param_finish)
		  offset++;
	      }
	    }
	}
      /* Add extra parameters for jumps */
      for (i=1;i<=psr[p].nJumps;i++)
	{
	  if (psr[p].fitJump[i]==1)
	    offset++;
	}
      /* Add extra parameters for sinusoidal whitening */
      if (psr[p].param[param_wave_om].fitFlag[0]==1)
	offset+=psr[p].nWhite*2-1;
      if (psr[p].param[param_quad_om].fitFlag[0]==1)
	offset+=psr[p].nQuad*4-1;
      if (psr[p].param[param_ifunc].fitFlag[0]==1)
	offset+=psr[p].ifuncN-1;
      if (psr[p].param[param_clk_offs].fitFlag[0]==1)
	offset+=psr[p].clkOffsN-2;
      if (psr[p].param[param_tel_dx].fitFlag[0]==1 &&
	  psr[p].param[param_tel_dx].val[0] < 2)
	offset+=(psr[p].nTelDX-1);
      else if (psr[p].param[param_tel_dx].fitFlag[0]==1 &&
	       psr[p].param[param_tel_dx].val[0] == 2)
	offset+=(psr[p].nTelDX-2);

      if (psr[p].param[param_tel_dy].fitFlag[0]==1 &&
	  psr[p].param[param_tel_dy].val[0] < 2)
	offset+=(psr[p].nTelDY-1);
      else if (psr[p].param[param_tel_dy].fitFlag[0]==1 &&
	       psr[p].param[param_tel_dy].val[0] == 2)
	offset+=(psr[p].nTelDY-2);

      if (psr[p].param[param_tel_dz].fitFlag[0]==1 &&
	  psr[p].param[param_tel_dz].val[0] < 2)
	offset+=(psr[p].nTelDZ-1);
      else if (psr[p].param[param_tel_dz].fitFlag[0]==1 &&
	       psr[p].param[param_tel_dz].val[0] == 2)
	offset+=(psr[p].nTelDZ-2);

      if (psr[p].param[param_quad_ifunc_p].fitFlag[0]==1)
	offset+=psr[p].quad_ifuncN_p-1;
      if (psr[p].param[param_quad_ifunc_c].fitFlag[0]==1)
	offset+=psr[p].quad_ifuncN_c-1;
      if (psr[p].param[param_gwsingle].fitFlag[0]==1)
	offset+=4;
      offset++; // For arbitrary phase
    }
  //  printf("At this point\n");
  if (psr[0].param[param_wave_om].fitFlag[0]==2)
    {
      double sx[MAX_OBSN],sy[MAX_OBSN],sy2[MAX_OBSN],sye[MAX_OBSN];
      double xval;
      double earliestTime,latestTime;
      int npt;

      FILE *fout,*fout2;
      printf("Outputting wave.dat file\n");

      fout = fopen("wave.dat","w");
      fout2 = fopen("ept2tt.clk","w");
      fprintf(fout2,"# TT(TAI) TT(EPT)\n# EPT derived using the pte plugin package\n#\n");

      earliestTime = psr[0].obsn[i].sat;
      latestTime = psr[0].obsn[psr[0].nobs-1].sat;
      for (p=1;p<npsr;p++)
	{
	  if ((double)psr[p].obsn[i].sat < earliestTime) earliestTime = (double)psr[p].obsn[i].sat-1;
	  if ((double)psr[p].obsn[psr[p].nobs-1].sat > latestTime) latestTime = (double)psr[p].obsn[psr[p].nobs-1].sat+1;
	}
      npt = (int)((latestTime-earliestTime)/14.0);

      for (i=0;i<npt;i++)
	{
	  sx[i] = earliestTime+(i*14.0)-psr[0].param[param_waveepoch].val[0];
	  sy[i] = 0.0;
	  xval = sx[i];
	  sye[i] = 0.0;
	  for (j=0;j<psr[0].nWhite;j++)
	    {
	      sy[i]+=(psr[0].wave_cos[j]*cos((j+1)*psr[0].param[param_wave_om].val[0]*xval)+psr[0].wave_sine[j]*sin((j+1)*psr[0].param[param_wave_om].val[0]*xval));
	      sye[i] += pow(psr[0].wave_cos_err[j]*cos((j+1)*psr[0].param[param_wave_om].val[0]*xval),2);
	      sye[i] += pow(psr[0].wave_sine_err[j]*sin((j+1)*psr[0].param[param_wave_om].val[0]*xval),2);
	    }
	  sy2[i] = sy[i];
	  sye[i] = sqrt(sye[i]);
	}
      TKremovePoly_d(sx,sy2,npt,3);
      for (i=0;i<npt;i++)
	{
	  fprintf(fout,"%.5f %g %g %g\n",(double)earliestTime+14*i,(double)sy[i],(double)sy2[i],(double)sye[i]);
	  fprintf(fout2,"%.5f %g\n",(double)earliestTime+14*i,(double)sy2[i]);
	}
      fclose(fout);
      fclose(fout2);
    }
  
  
  printf("Freeing the memory\n");
  free(x);
  free(y);
  free(sig);
  free(val);
  free(error);
  for (i=0;i<npol;i++)
    free(covar[i]);
  free(covar);
}


void globalFITfuncs(double x,double afunc[],int ma,pulsar *psr,int counter)
{
  int i;
  int n=0;
  int p,pp;
  int new_ma;
  int ipos;
  int j,k;
  int tot=0;
  int nglobal=0;

  // MUST DO SOMETHING WITH CONSTRAINTS
  //  printf("Here at the top\n");
  //  printf("Here with joe %g %d %d %s\n",x,ma,counter,psr->name);
  for (i=0;i<ma;i++) afunc[i]=0.0;

  for (p=0;p<gnpsr;p++)
    {
      if (counter < tot+psr[p].nobs+psr[p].nconstraints) break;
      tot+=psr[p].nobs+psr[p].nconstraints;
    }
  //  printf("pos 1\n");
  ipos = counter-tot;
  //  printf("Here with joe2 %s %d %d %d %d\n",psr[p].name,ipos,counter,tot,psr[p].nobs);
  new_ma=1;
  // Add global parameters
  for (i=0;i<MAX_PARAMS;i++)
    {
      for (k=0;k<psr[0].param[i].aSize;k++)
	{
	  if (psr[0].param[i].fitFlag[k] == 2)
	    nglobal++;
	}
    }
  //  printf("pos 2\n");
  /* Add extra parameters for sinusoidal whitening */
  if (psr[p].param[param_wave_om].fitFlag[0]==2)
    nglobal+=psr[p].nWhite*2-1;
  if (psr[p].param[param_quad_om].fitFlag[0]==2)
    nglobal+=psr[p].nQuad*4-1;
  if (psr[p].param[param_ifunc].fitFlag[0]==2)
    nglobal+=psr[p].ifuncN-1;
  if (psr[p].param[param_clk_offs].fitFlag[0]==2)
    nglobal+=psr[p].clkOffsN-2;
  if (psr[p].param[param_tel_dx].fitFlag[0]==2 &&
      psr[p].param[param_tel_dx].val[0] < 2)
    nglobal+=(psr[p].nTelDX-1);
  else if (psr[p].param[param_tel_dx].fitFlag[0]==2 &&
      psr[p].param[param_tel_dx].val[0] == 2)
    nglobal+=(psr[p].nTelDX-2);

  if (psr[p].param[param_tel_dy].fitFlag[0]==2 &&
      psr[p].param[param_tel_dy].val[0] < 2)
    nglobal+=(psr[p].nTelDY-1);
  else if (psr[p].param[param_tel_dy].fitFlag[0]==2 &&
      psr[p].param[param_tel_dy].val[0] == 2)
    nglobal+=(psr[p].nTelDY-2);

  if (psr[p].param[param_tel_dz].fitFlag[0]==2 &&
      psr[p].param[param_tel_dz].val[0] < 2)
    nglobal+=(psr[p].nTelDZ-1);
  else if (psr[p].param[param_tel_dz].fitFlag[0]==2 &&
      psr[p].param[param_tel_dz].val[0] == 2)
    nglobal+=(psr[p].nTelDZ-2);

  if (psr[p].param[param_quad_ifunc_p].fitFlag[0]==2)
    nglobal+=psr[p].quad_ifuncN_p-1;
  if (psr[p].param[param_quad_ifunc_c].fitFlag[0]==2)
    nglobal+=psr[p].quad_ifuncN_c-1;
  if (psr[p].param[param_gwsingle].fitFlag[0]==2)
    nglobal+=(4-1);
  //  printf("pos 3\n");
  new_ma+=nglobal;
  for (i=0;i<MAX_PARAMS;i++)
    {
      for (k=0;k<psr[p].param[i].aSize;k++)
	{
	  if (psr[p].param[i].fitFlag[k]==1) {
	    if (i!=param_start && i!=param_finish)
	      new_ma++;
	  }
	}
    }
  //  printf("pos 4\n");
  /* Add extra parameters for jumps */
  for (i=1;i<=psr[p].nJumps;i++)
    {
      if (psr[p].fitJump[i]==1)
	new_ma++;
    }
  //  printf("pos 5\n");
  /* Add extra parameters for sinusoidal whitening */
  if (psr[p].param[param_wave_om].fitFlag[0]==1)
    new_ma+=psr[p].nWhite*2-1;
  if (psr[p].param[param_quad_om].fitFlag[0]==1)
    new_ma+=psr[p].nQuad*4-1;
  if (psr[p].param[param_ifunc].fitFlag[0]==1)
    new_ma+=psr[p].ifuncN-1;
  if (psr[p].param[param_clk_offs].fitFlag[0]==1)
    new_ma+=psr[p].clkOffsN-2;
  if (psr[p].param[param_tel_dx].fitFlag[0]==1 && 
      psr[p].param[param_tel_dx].val[0] < 2)
    new_ma+=(psr[p].nTelDX-1);
  else if (psr[p].param[param_tel_dx].fitFlag[0]==1 && 
      psr[p].param[param_tel_dx].val[0] == 2)
    new_ma+=(psr[p].nTelDX-2);

  if (psr[p].param[param_tel_dy].fitFlag[0]==1 && 
      psr[p].param[param_tel_dy].val[0] < 2)
    new_ma+=(psr[p].nTelDY-1);
  else if (psr[p].param[param_tel_dy].fitFlag[0]==1 && 
      psr[p].param[param_tel_dy].val[0] == 2)
    new_ma+=(psr[p].nTelDY-2);

  if (psr[p].param[param_tel_dz].fitFlag[0]==1 && 
      psr[p].param[param_tel_dz].val[0] < 2)
    new_ma+=(psr[p].nTelDZ-1);
  else if (psr[p].param[param_tel_dz].fitFlag[0]==1 && 
      psr[p].param[param_tel_dz].val[0] == 2)
    new_ma+=(psr[p].nTelDZ-2);

  if (psr[p].param[param_quad_ifunc_p].fitFlag[0]==1)
    new_ma+=psr[p].quad_ifuncN_p-1;
  if (psr[p].param[param_quad_ifunc_c].fitFlag[0]==1)
    new_ma+=psr[p].quad_ifuncN_c-1;
  if (psr[p].param[param_gwsingle].fitFlag[0]==1)
    new_ma+=4;
  //  printf("pos 6\n");
  // Now calculate position in afunc array
  n=0;
  // Add global parameters
  n+=nglobal;
  for (pp=0;pp<p;pp++)
    {
      n++; // For the offset
      for (i=0;i<MAX_PARAMS;i++)
	{
	  for (k=0;k<psr[pp].param[i].aSize;k++)
	    {
	      if (psr[pp].param[i].fitFlag[k]==1) {
		if (i!=param_start && i!=param_finish)
		  n++;
	      }
	    }
	}
      // Should do the same for jumps
      /* Add extra parameters for jumps */
      for (i=1;i<=psr[pp].nJumps;i++)
	{
	  if (psr[pp].fitJump[i]==1)
	    n++;
	}
      //  printf("pos 7\n");
      /* Add extra parameters for sinusoidal whitening */
      if (psr[pp].param[param_wave_om].fitFlag[0]==1)
	n+=psr[pp].nWhite*2-1;
      if (psr[pp].param[param_quad_om].fitFlag[0]==1)
	n+=psr[pp].nQuad*4-1;
      if (psr[pp].param[param_ifunc].fitFlag[0]==1)
	n+=psr[pp].ifuncN-1;
      if (psr[pp].param[param_clk_offs].fitFlag[0]==1)
	n+=psr[pp].clkOffsN-2;
      if (psr[p].param[param_tel_dx].fitFlag[0]==1 && 
	  psr[p].param[param_tel_dx].val[0] < 2)
	n+=(psr[p].nTelDX-1);
      else if (psr[p].param[param_tel_dx].fitFlag[0]==1 && 
	       psr[p].param[param_tel_dx].val[0] == 2)
	n+=(psr[p].nTelDX-2);

      if (psr[p].param[param_tel_dy].fitFlag[0]==1 && 
	  psr[p].param[param_tel_dy].val[0] < 2)
	n+=(psr[p].nTelDY-1);
      else if (psr[p].param[param_tel_dy].fitFlag[0]==1 && 
	       psr[p].param[param_tel_dy].val[0] == 2)
	n+=(psr[p].nTelDY-2);

      if (psr[p].param[param_tel_dz].fitFlag[0]==1 && 
	  psr[p].param[param_tel_dz].val[0] < 2)
	n+=(psr[p].nTelDZ-1);
      else if (psr[p].param[param_tel_dz].fitFlag[0]==1 && 
	       psr[p].param[param_tel_dz].val[0] == 2)
	n+=(psr[p].nTelDZ-2);

      if (psr[pp].param[param_quad_ifunc_p].fitFlag[0]==1)
	n+=psr[pp].quad_ifuncN_p-1;
      if (psr[pp].param[param_quad_ifunc_c].fitFlag[0]==1)
	n+=psr[pp].quad_ifuncN_c-1;
      if (psr[pp].param[param_gwsingle].fitFlag[0]==1)
	n+=4-1;
      //  printf("pos 8\n");
    }

  // Global fit
  int c=0;
  //  printf("pos 9\n");
  for (i=0;i<MAX_PARAMS;i++)
    {
      //      printf("Checking param %d %s\n",i,psr[p].param[i].label[0]);
      for (k=0;k<psr[p].param[i].aSize;k++)
	{
	  if (psr[p].param[i].fitFlag[k] == 2)  
	    {
	      //	      afunc[c] = dotproduct(psr[p].posPulsar,psr[p].obsn[ipos].planet_ssb[4]);
	      if (i==param_wave_om)
		{
		  for (j=0;j<psr[p].nWhite*2;j++)
		    {
		      // Note that a global fit to white data must be related to the same time
		      afunc[c] = getParamDeriv(&psr[p],ipos,x+(double)psr[p].param[param_pepoch].val[0] - (double)psr[0].param[param_waveepoch].val[0],i,j);
		      c++;		      
		    }
		}
	      else if (i==param_quad_om) 
		{
		  for (j=0;j<psr[p].nQuad*4;j++)
		    {
		      //		      printf("Calling fit\n");
		      afunc[c] = getParamDeriv(&psr[p],ipos,x+(double)psr[p].param[param_pepoch].val[0],i,j);
		      c++;		      
		    }
		}
	      else if (i==param_gwsingle)
		{
		  afunc[c] = getParamDeriv(&psr[p],ipos,x,i,0);
		  c++;
		  afunc[c] = getParamDeriv(&psr[p],ipos,x,i,1);
		  c++;
		  afunc[c] = getParamDeriv(&psr[p],ipos,x,i,2);
		  c++;
		  afunc[c] = getParamDeriv(&psr[p],ipos,x,i,3);
		  c++;
		}
	      else if (i==param_ifunc)
		{
		  for (j=0;j<psr[p].ifuncN;j++)
		    {
		      afunc[c] = getParamDeriv(&psr[p],ipos,x,i,j);
		      c++;
		    }
		}
	      else if (i==param_clk_offs)
		{
		  for (j=0;j<psr[p].clkOffsN-1;j++)
		    {
		      afunc[c] = getParamDeriv(&psr[p],ipos,x,i,j);
		      c++;
		    }
		}
	      else if (i==param_tel_dx)
		{
		  if (psr[p].param[param_tel_dx].val[0] < 2)
		    {
		      for (j=0;j<psr[p].nTelDX;j++)
			{
			  afunc[c] = getParamDeriv(&psr[p],ipos,x,i,j);
			  c++;
			}
		    }
		  else if (psr[p].param[param_tel_dx].val[0] == 2)
		    {
		      for (j=0;j<psr[p].nTelDX-1;j++)
			{
			  afunc[c] = getParamDeriv(&psr[p],ipos,x,i,j);
			  c++;
			}
		    }
		}
	      else if (i==param_tel_dy)
		{
		  if (psr[p].param[param_tel_dy].val[0] < 2)
		    {
		      for (j=0;j<psr[p].nTelDY;j++)
			{
			  afunc[c] = getParamDeriv(&psr[p],ipos,x,i,j);
			  c++;
			}
		    }
		  else if (psr[p].param[param_tel_dy].val[0] == 2)
		    {
		      for (j=0;j<psr[p].nTelDY-1;j++)
			{
			  afunc[c] = getParamDeriv(&psr[p],ipos,x,i,j);
			  c++;
			}
		    }
		}
	      else if (i==param_tel_dz)
		{
		  if (psr[p].param[param_tel_dz].val[0] < 2)
		    {
		      for (j=0;j<psr[p].nTelDZ;j++)
			{
			  afunc[c] = getParamDeriv(&psr[p],ipos,x,i,j);
			  c++;
			}
		    }
		  else if (psr[p].param[param_tel_dz].val[0] == 2)
		    {
		      for (j=0;j<psr[p].nTelDZ-1;j++)
			{
			  afunc[c] = getParamDeriv(&psr[p],ipos,x,i,j);
			  c++;
			}
		    }
		}
	      else if (i==param_quad_ifunc_p)
		{
		  for (j=0;j<psr[p].quad_ifuncN_p;j++)
		    {
		      afunc[c] = getParamDeriv(&psr[p],ipos,x,i,j);
		      c++;
		    }
		}
	      else if (i==param_quad_ifunc_c)
		{
		  for (j=0;j<psr[p].quad_ifuncN_c;j++)
		    {
		      afunc[c] = getParamDeriv(&psr[p],ipos,x,i,j);
		      c++;
		    }
		}
	      else
		{
		  afunc[c] = getParamDeriv(&psr[p],ipos,x,i,k);
		  c++;
		}
	    }
	} 
    }
  //  printf("pos 10\n");
  //  printf("Global fit n = %d, nglobal = %d\n",n,nglobal);
  //  printf("Global fit = %g\n",afunc[0]);
  //  printf("Now calling fitfuncs\n");
  FITfuncs(x,afunc+n,new_ma-nglobal,&psr[p],ipos);
  //  printf("Finished\n");
  //  printf("-----------------------------\n");
  //  for (i=0;i<ma;i++)
  //    printf("have %g\n",afunc[i]);
  //  n+=new_ma;
  //  exit(1);
}
char * plugVersionCheck = TEMPO2_h_VER;
