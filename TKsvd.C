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

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "tempo2.h"
#include "TKsvd.h"
#include "T2toolkit.h"


/* Calculates SVD by following technique given in wikipedia */
void TKsingularValueDecomposition_lsq(longdouble **designMatrix,int n,int nf,longdouble **v,longdouble *w,longdouble **u)
{
   longdouble an;
   int i,j,k,its,l,nm,jj;
   int max_its = 40,pos1;
   longdouble c,s,f,g,h,y,z,x;
   longdouble rv1[nf];
   /* For A = U.W.V^T - obtain U, W and V */  

   /* Step 1: Reduce the matrix to a bidiagonal matrix */
   TKbidiagonal(designMatrix,&an,n,nf,v,w,u,rv1);
   /* Step 2: Compute the SVD of the bidiagonal matrix */

   // Diagonalisation of the bidiagonal form: Loop over singular values
   // Code based on the numerical recipes routine
   for (k=nf-1;k>=0;k--)
   {
	  for (its=1;its<=max_its;its++)
	  {
		 pos1=0;
		 for (l=k;l>=0;l--)
		 {
			nm=l-1;
			if ((fabs(rv1[l])+an)==an) {pos1=2; break;}
			if ((fabs(w[nm])+an)==an)  {pos1=1; break;}
		 }
		 if (pos1!=2)
		 {
			c=0.0;
			s=1.0;
			for (i=l;i<=k;i++) // Check <= sign
			{
			   f=s*rv1[i];		  
			   rv1[i]=c*rv1[i];
			   if ((fabs(f)+an)==an) break;
			   g=w[i];
			   h=TKpythag(f,g);
			   w[i]=h;
			   h=1.0/h;
			   c= (g*h);
			   s= -(f*h);
			   for (j=0;j<n;j++)
			   {
				  y = designMatrix[j][nm];
				  z = designMatrix[j][i];
				  designMatrix[j][nm] = (y*c)+(z*s);
				  designMatrix[j][i] = -(y*s)+(z*c);
			   }
			}
		 }	      
		 z = w[k];
		 if (l==k)
		 {
			if (z < 0.0)
			{
			   w[k] = -z;
			   for (j=0;j<nf;j++)
				  v[j][k] = -v[j][k];
			}
			break;
		 }
		 if (its==30)
		 {
			printf("1. No convergence in singular value decomposition after 30 iterations\n");
			exit(1);
		 }
		 x = w[l];
		 nm = k-1;
		 y = w[nm];
		 g =rv1[nm];
		 h= rv1[k];
		 f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
		 g = TKpythag(f,1.0);
		 f = ((x-z)*(x+z)+h*((y/(f+TKsign_d(g,f)))-h))/x;
		 c=1.0;
		 s=1.0;
		 for (j=l;j<=nm;j++)
		 {
			i = j+1;
			g = rv1[i];
			y = w[i];
			h = s*g;
			g = c*g;
			z = TKpythag(f,h);
			rv1[j] = z;
			c = f/z;
			s = h/z;
			f = (x*c)+(g*s);
			g = -(x*s)+(g*c);
			h = y*s;
			y = y*c;
			for (jj=0;jj<nf;jj++)
			{
			   x = v[jj][j];
			   z = v[jj][i];
			   v[jj][j] = (x*c)+(z*s);
			   v[jj][i] = -(x*s)+(z*c);		  
			}
			z = TKpythag(f,h);
			w[j] = z;
			if (z != 0.0)
			{
			   z = 1.0/z;
			   c=f*z;
			   s=h*z;
			}
			f = (c*g)+(s*y);
			x = -(s*g)+(c*y);
			for (jj=0;jj<n;jj++)
			{
			   y = designMatrix[jj][j];
			   z = designMatrix[jj][i];
			   designMatrix[jj][j] = (y*c)+(z*s);
			   designMatrix[jj][i] = -(y*s)+(z*c);		  
			}
		 }
		 rv1[l] = 0.0;
		 rv1[k] = f;
		 w[k]=x;
	  }
   }

}


/* Use Householder reflections to do reduce the matrix to bidiagonal form       */
/* This is converted from the Fortran EISPACK code which is very similar        */
/* to the more recent numerical recipes code - originally this came from        */
/* the algol procedure svd, num. math. 14, 403-420 (1970) by Golub and Reinsch  */
/* The Fortran code was developed by Burton S. Garbow from the Argonne National */
/* laboratory (1983)                                                            */

void TKbidiagonal(longdouble **a,longdouble *an,int ndata,int nfit,longdouble **v,longdouble *w,longdouble **u,longdouble *rv1)
{
   int i,j,k,l;
   longdouble g=0.0;
   longdouble scale=0.0,s,f,h;

   *an=0.0;
   for (i=0;i<nfit;i++)
   {
	  l=i+1;
	  rv1[i] = scale*g;
	  g=0.0; s=0.0; scale = 0.0;
	  if (i <= ndata)
	  {
		 for (k=i;k<ndata;k++)
			scale+=fabs(a[k][i]);
		 if (scale != 0.0)
		 {
			for (k=i;k<ndata;k++)
			{
			   a[k][i]/=scale;
			   s+=a[k][i]*a[k][i];
			}
			f=a[i][i];
			g=-TKsign_d(sqrt(s),f);
			h = f*g-s;
			a[i][i]=f-g;
			for (j=l;j<nfit;j++)
			{
			   s=0.0;
			   for (k=i;k<ndata;k++)
				  s+=a[k][i]*a[k][j];
			   f = s/h;
			   for (k=i;k<ndata;k++)
				  a[k][j]+=f*a[k][i];
			}
			for (k=i;k<ndata;k++)
			   a[k][i]*=scale;	      

		 }
	  }

	  w[i]=scale*g;
	  g=0.0;
	  s=0.0;
	  scale=0.0;
	  if (i<=ndata && i != nfit)
	  {
		 for (k=l;k<nfit;k++)
			scale+=fabs(a[i][k]);
		 if (scale != 0.0)
		 {
			for (k=l;k<nfit;k++)
			{
			   a[i][k] = a[i][k]/scale;
			   s+=a[i][k]*a[i][k];
			}
			f=a[i][l];
			g=-TKsign_d(sqrt(s),f);
			h=f*g-s;
			a[i][l] = f-g;
			for (k=l;k<nfit;k++)
			   rv1[k]=a[i][k]/h;
			for (j=l;j<ndata;j++)
			{
			   s=0.0;
			   for (k=l;k<nfit;k++)
				  s+=a[j][k]*a[i][k];
			   for (k=l;k<nfit;k++)
				  a[j][k]+=s*rv1[k];
			}
			for (k=l;k<nfit;k++)
			   a[i][k]=scale*a[i][k];
		 }
	  }
	  *an=TKretMax_d(*an,(fabs(w[i])+fabs(rv1[i])));
   }

   // Accumulation of right-hand transformations
   for (i=nfit-1;i>=0;i--)
   {
	  if (i < nfit-1)
	  {
		 if (g != 0.0)
		 {
			for (j=l;j<nfit;j++) v[j][i]=(a[i][j]/a[i][l])/g;
			for (j=l;j<nfit;j++)
			{
			   s=0.0;
			   for (k=l;k<nfit;k++) s+=a[i][k]*v[k][j];
			   for (k=l;k<nfit;k++) v[k][j]+=s*v[k][i];
			}			 
		 }
		 for (j=l;j<nfit;j++)
			v[i][j] = v[j][i] = 0.0;
	  }
	  v[i][i]=1.0;
	  g = rv1[i];
	  l=i;
   }

   // Accumulation of left-hand transformations
   for (i=TKretMin_i(ndata,nfit)-1;i>=0;i--)
   {
	  l=i+1;
	  g=w[i];
	  for (j=l;j<nfit;j++) a[i][j]=0.0;
	  if (g != 0.0)
	  {
		 g=1.0/g;
		 for (j=l;j<nfit;j++)
		 {
			s=0.0;
			for (k=l;k<ndata;k++) s+=a[k][i]*a[k][j];
			f = (s/a[i][i])*g;
			for (k=i;k<ndata;k++) a[k][j]=a[k][j]+f*a[k][i];
		 }
		 for (j=i;j<ndata;j++) a[j][i]*=g;
	  }
	  else
	  {
		 for (j=i;j<ndata;j++) a[j][i]=0.0;
	  }
	  a[i][i]++;
   }

}

/* Solves A.X = B for vector X using equation 2.6.7 in numerical recipes */
/* equation: x = V . [diag(1/w_j)] . (U^T.b)                             */ 
/*                                                                       */
/* Returns 'x'                                                           */
void TKbacksubstitution_svd(longdouble **V, longdouble *w,longdouble **U,longdouble *b,longdouble *x,int n,int nf)
{
   int i,j;
   longdouble uTb[nf];

   /* Calculate [diag(1/w_j)] . U^T.b */
   for (i=0;i<nf;i++)
   {
	  uTb[i]=0.0;
	  if (w[i]!=0.0)
	  {
		 for (j=0;j<n;j++)
			uTb[i]+=U[j][i]*b[j];
		 uTb[i]/=w[i];
	  }
   }
   /* Now multiply by V as in equation 2.6.7 */
   for (i=0;i<nf;i++)
   {
	  x[i]=0.0;
	  for (j=0;j<nf;j++)
		 x[i]+=V[i][j]*uTb[j];
   }
}

/* Computes (a^2 + b^2)^1/2 */
longdouble TKpythag(longdouble a,longdouble b)
{
   longdouble ret=0.0;
   longdouble absa,absb;

   absa = fabs(a);
   absb = fabs(b);
   if (absa > absb)
	  ret = absa*sqrt(1.0+pow(absb/absa,2));
   else
   {
	  if (absb==0) ret = 0.0;
	  else ret = absb*sqrt(1.0+pow(absa/absb,2));
   }
   return ret;
}

