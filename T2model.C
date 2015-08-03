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

/* Generalised timing model for tempo2                 
 *
 * Based on the DD model, but includes
 *
 *  conversion to ELL1 model if EPS1 and EPS2 are set
 *  BT model                                          (set allTerms = 0)
 *  jumps from BTJ model
 *  use of SHAPMAX (i.e. DDS model)
 *  extra terms implemented in DDK model
 *  flag to convert to DDGR model
 *  multiple binary systems (e.g. planetary systems)
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tempo2.h"


static double calcDH( double ae, double h3, double h4, int nharm, int sel);
long double getParameter(pulsar *psr,int p,int k);

void calcGR(double mtot,double m2,double x,double ecc,double an,double afac,
            double f0, double *dr,double *dth,double *er,double *eth,
            double *xk,double *si,double *gamma, double *pbdot,double *a0,
            double *b0);
void getKeplerian(pulsar *psr,int com,double *pb,long double *t0,double *ecc,
                  double *omz,double *x,double *eps1,double *eps2,
                  long double *t0asc,double *shapmax,double *kom,double *kin);
void addKeplerianJumps(pulsar *psr,int ipos,double *torb,double *x,double *ecc,
                       double *omz,double *pb);
void getPostKeplerian(pulsar *psr,int com,double an,double *si,double *m2,
                      double *mtot,double *omdot, double *gamma,double *xdot,
                      double *xpbdot, double *pbdot, double *edot,double *pmra,
                      double *pmdec,double *dpara, double *dr,double *dth,
                      double *a0,double *b0,double *xomdot,double *afac,
                      double *eps1dot,double *eps2dot,double *daop);
void updateParameters(double edot,double xdot,double eps1dot,double eps2dot,
                      long double tt0,double *ecc,double *x,double *eps1,
                      double *eps2);
void deriveKeplerian(double pb,double kom,double *an,double *sin_omega,
                     double *cos_omega);
void derivePostKeplerian(double mtot,double m2,double dr,double dth,
                         double ecc,double *m1,double *er,double *eth);
void KopeikinTerms(pulsar *psr,int ipos,double ki,double pmra,double sin_omega,
                   double pmdec, double cos_omega,long double tt0,double dpara, 
                   double daop, double si,double *x, long double *DK011, 
                   long double *DK012, long double *DK021,long double *DK022,
                   long double *DK031, long double *DK032, long double *DK041, 
                   long double *DK042, long double *DK013, long double *DK014, 
                   long double *DK023, long double *DK024, long double *DK033, 
                   long double *DK034, long double *DK043, long double *DK044);
void computeU(double phase,double ecc,double *u);

double T2model(pulsar *psr,int p,int ipos,int param,int arr)
{
  double an;
  double pb,omdot;
  double rad2deg = 180.0/M_PI;
  double SUNMASS = 4.925490947e-6;
  long double tt0,t0,ct,t0asc;
  double m2,x,ecc,er,xdot,edot,dr,dth,eth;
  double pbdot,xpbdot,phase,u,gamma;
  double orbits;
  int    norbits;
  double cu,onemecu=0,cae,sae,ae,omega,omz,sw,cw,alpha,beta,bg,dre,drep,drepp,
    anhat,su=0;
  double sqr1me2,cume,brace,si,dlogbr,ds,da,a0,b0,d2bar,torb;
  double csigma,ce,cx,comega,cgamma,cdth,cm2,csi, ckom, ckin;
  double eps1,eps2,eps1dot,eps2dot;
  double ceps1,ceps2;
  double shapmax,cshapmax,sdds;
  int    com,com1,com2;
  int    allTerms=1;            /* = 0 to emulate BT model */
  double dpara;
  double pmra,pmdec;
  double sin_omega,cos_omega,ki;
  double mtot,m1,xk,xomdot,afac,kom;
  longdouble DK011,DK012, DK021, DK022, DK031,DK032, DK041,DK042,C,S;
  longdouble DK013, DK014, DK023, DK024, DK033, DK034, DK043, DK044;
  longdouble DAOP=0.0L, DSR=0.0L;
  long double DOP=0.0L; // Orbital parallax time delay
  double daop;// DAOP is the time delay due to annual orbital
              // parallax. daop is the aop distance.

  long double h3,h4,stig;
  long double lgf,TrueAnom;
  long double lsc, fs;
  int nharm=4;
  int mode = -1; // See ELL1Hmodel.C

  torb = 0.0;
  const char *CVS_verNum = "$Revision: 1.18 $";


  if (displayCVSversion == 1) 
    CVSdisplayVersion("T2model.C","T2model()",CVS_verNum);

  if (param==-1) 
    {
      com1 = 0;
      com2 = psr[p].nCompanion;
    }
  else
    {
      com1 = arr;
      com2 = arr+1;
    }

  //    printf("Number of companions = %d %d\n",com1,com2);
  
  for (com = com1; com < com2;com++)
    {      
      /* Obtain Keplerian parameters */   
      getKeplerian(&psr[p],com,&pb,&t0,&ecc,&omz,&x,&eps1,&eps2,&t0asc,
                   &shapmax,&kom,&ki);
      /* Now add in the jumps */
      addKeplerianJumps(&psr[p],ipos,&torb,&x,&ecc,&omz,&pb);
      /* Parameters derived from the Keplerian parameters */
      deriveKeplerian(pb,kom,&an,&sin_omega,&cos_omega);
      /* Obtain post-Keplerian parameters */
      getPostKeplerian(&psr[p],com,an,&si,&m2,&mtot,&omdot,&gamma,&xdot,&xpbdot,
                       &pbdot,&edot,&pmra,&pmdec,&dpara,&dr,&dth,&a0,&b0,
                       &xomdot,&afac,&eps1dot,&eps2dot,&daop);

      /* If the beta-prime parameters are set then omdot, gamma, si, dr, er, 
         dth, a0 and b0 can be calculated from the beta-prime values
       * - see papers of Taylor, Wolszczan, Damour & Weisberg 1992 - Nature
       *                 Damour & Taylor 1992 - Phys Rev D
       */
      if (psr[p].param[param_bp].paramSet[0]==1 && 
          psr[p].param[param_bpp].paramSet[0]==1)
        {
          /*	  useBeta(psr[p]); */
          printf("Beta model not implemented yet\n");
        }

      /* If general relativity is assummed to be correct and 
       * the total system mass has been determined then the values 
       * of dr,dth,er, eth, sini, gamma and pbdot can be calculated
       */
      if (psr[p].param[param_mtot].paramSet[com]==1)
        {
          calcGR(mtot,m2,x,ecc,an,afac,(double)psr[p].param[param_f].val[0],
                 &dr,&dth,&er,&eth,&xk,&si,&gamma,&pbdot,&a0,&b0);
          omdot   = xomdot + xk; 
        }
      /* Derive parameters from the post-Keplerian parameters */
      derivePostKeplerian(mtot,m2,dr,dth,ecc,&m1,&er,&eth);
      /* Obtain delta T */
      ct  = psr[p].obsn[ipos].bbat;      
      if (psr[p].param[param_t0].paramSet[com]==1)        
        tt0 = (ct-t0)*SECDAY;
      else if( psr[p].param[param_tasc].paramSet[com]==1) 
        tt0 = (ct-t0asc)*SECDAY;
      else {
        printf("ERROR [T2] T0 or TASC needs to be set in the parameter file\n");
        exit(1);
      }
      //      logdbg("going to update Parameters");
      /* Update parameters with their time derivatives */
      updateParameters(edot,xdot,eps1dot,eps2dot,tt0,&ecc,&x,&eps1,&eps2);
      //      logdbg("updated parameters");

      /* Do some checks */
      if (ecc < 0.0 || ecc > 1.0)
        {
          displayMsg(1,"BIN2","Eccentricity of orbit < 0, setting to 0","",
                     psr[p].noWarnings);
          psr[p].param[param_ecc].val[com]=0.0;
          ecc = 0.0;
        }

      /* Obtain number of orbits in tt0 */
      orbits  = tt0/pb - 0.5*(pbdot+xpbdot)*pow(tt0/pb,2);
      norbits = (int)orbits;
      if (orbits<0.0) norbits--;
      
      /* Obtain phase of orbit */
      phase=2.0*M_PI*(orbits-norbits);
      //      printf("Orbit phase = %.15Lf %.15g\n",psr[p].obsn[ipos].bbat,phase);
      if (psr[p].param[param_ecc].paramSet[com]==1)
        {
	  //          logdbg("going to compute U");
          /* Compute eccentric anomaly u by iterating Kepler's equation */	 
          computeU(phase,ecc,&u);
	  //          logdbg("computed U");
          /*  DD equations 17b and 17c */
          su=sin(u); cu=cos(u);
          onemecu=1.0-ecc*cu;
          cae=(cu-ecc)/onemecu;                         /* Equation 17b */
          sae=sqrt(1.0-pow(ecc,2))*su/onemecu;          /* Equation 17c */
          ae=atan2(sae,cae);
          if(ae<0.0) ae=ae+2.0*M_PI;
          ae=2.0*M_PI*orbits + ae - phase;
          omega=omz/rad2deg + omdot*ae;
          sw=sin(omega);
          cw=cos(omega);
	  //          logdbg("In the middle of DD");
          /* DD equations 26, 27, 57: */
          sqr1me2=sqrt(1-pow(ecc,2));
          cume=cu-ecc;
	  //          logdbg("going to Kopeikin");
          /* Update parameters due to proper motion - Kopeikin 1996 */
          /* And annual-orbital and orbital parallax - Kopeikin 1995 */
          if (psr[p].param[param_kin].paramSet[com]==1 && 
              psr[p].param[param_kom].paramSet[com]==1 &&
              (psr[p].param[param_pmra].paramSet[com]==1 || 
               psr[p].param[param_pmdec].paramSet[com]==1))
            {
	      //              logdbg("going to do KopeikinTerms");
              KopeikinTerms(&psr[p],ipos,ki,pmra,sin_omega,pmdec,cos_omega,tt0,
                            dpara,daop,si,&x,&DK011,&DK012,&DK021,&DK022,
                            &DK031,&DK032,&DK041,&DK042,&DK013,&DK014,&DK023,
                            &DK024,&DK033,&DK034,&DK043,&DK044);
	      //              logdbg("did KopeikinTerms");
              C = (longdouble)(cw*(cu-er)-sqrt(1.0-pow(eth,2.0))*sw*su);
              S = (longdouble)(sw*(cu-er)+cw*sqrt(1.0-pow(eth,2.0))*su);
              DAOP = (DK011+DK012)*C-(DK021+DK022)*S;
              DSR = (DK031+DK032)*C+(DK041+DK042)*S;
              DOP = dpara/AULTSC/2.0*pow(x,2.0)*
                ( pow( sin( ki ), -2.0 ) - 0.5 + 0.5 * pow( ecc, 2.0 ) *
                  ( 1 + pow( sw, 2.0 ) - 3.0 / pow( sin( ki ), 2.0 ) ) -
                  2.0 * ecc * ( pow( sin( ki ), -2.0 ) - pow( sw, 2.0 ) ) * 
                  cume - sqr1me2 * 2 * sw * cw * su * cume + 0.5 * 
                  ( cos( 2.0 * omega ) + pow( ecc, 2.0 ) * 
                    ( pow( sin( ki ), -2.0 ) + pow( cu, 2.0 ) ) ) *
                  cos( 2.0 * u ) );

	      //	      printf("T2model: %g DAOP = %g DSR = %g DOP = %g\n",(double)psr[p].obsn[ipos].bbat,(double)DAOP,(double)DSR,(double)DOP);

	      //                logdbg("DAOP is %g and DSR is %g\n", (double)DAOP, (double)DSR);
	      //                logdbg("DAOP is %g, DK011 and DK021 are %f and %f\n",
	      //                       (double)DAOP,(double)DK011,(double)DK021);
            }
	  
	  
          if (psr[p].param[param_shapmax].paramSet[com]==1)  /* DDS model */
            {
              sdds  = 1.0 - exp(-1.0*shapmax);
              brace = onemecu-sdds*(sw*cume+sqr1me2*cw*su);
            }
          else
            brace=onemecu-si*(sw*cume+sqr1me2*cw*su);

          da=a0*(sin(omega+ae) + ecc*sw) + b0*(cos(omega+ae) + 
                                               ecc*cw); /* Equation 27 */
	  
          /* DD equations 46 to 51 */	  
          alpha=x*sw;                                   /* Equation 46  */
          beta=x*sqrt(1-pow(eth,2))*cw;                 /* Equation 47  */
          bg=beta+gamma;
          dre=alpha*(cu-er) + bg*su;                    /* Equation 48  */
          drep=-alpha*su + bg*cu;                       /* Equation 49  */
          drepp=-alpha*cu - bg*su;                      /* Equation 50  */
          anhat=an/onemecu;                             /* Equation 51  */

	  dlogbr=log(brace);
	  ds=-2*m2*dlogbr;        /* Equation 26 */
	  
        }
      else if (psr[p].param[param_eps1].paramSet[com]==1)  /* ELL1 model */
        {
          dre  = x*(sin(phase)-0.5*(eps1*cos(2.0*phase)-eps2*sin(2.0*phase)));
          drep = x*cos(phase);
          drepp=-x*sin(phase);
	  //          logdbg("going to Kopeikin");
          /* Update parameters due to proper motion - Kopeikin 1996 */
          if (psr[p].param[param_kin].paramSet[com]==1 && 
              psr[p].param[param_kom].paramSet[com]==1 &&
              (psr[p].param[param_pmra].paramSet[com]==1 || 
               psr[p].param[param_pmdec].paramSet[com]==1))
            {
              S = (sin(phase)-0.5*(eps1*cos(2.0*phase)-eps2*sin(2.0*phase)));
              C = cos(phase)+0.5*(eps2*cos(2.0*phase)+eps1*sin(2.0*phase));
              KopeikinTerms(&psr[p],ipos,ki,pmra,sin_omega,pmdec,cos_omega,tt0,
                            dpara,daop,si,&x,&DK011,&DK012,&DK021,&DK022,
                            &DK031,&DK032,&DK041,&DK042,&DK013,&DK014,&DK023,
                            &DK024,&DK033,&DK034,&DK043,&DK044);
              DAOP = (DK011+DK012)*C-(DK021+DK022)*S;
              DSR = (DK031+DK032)*C+(DK041+DK042)*S;
            }
	  
          brace=1-si*sin(phase);
          da=a0*sin(phase)+b0*cos(phase);  
	  
          anhat = an; ecc = 0.0;
	  
	  
	  /* Shapiro delay */
	  if ( psr[p].param[param_h3].paramSet[0] * psr[p].param[param_stig].paramSet[0] == 1 
	       || psr[p].param[param_h3].paramSet[0] * psr[p].param[param_h4].paramSet[0] == 1){
	    
	    //	    printf("Using the Friere & Wex formalism for the Shapiro delay\n");
	    // Based on ELL1Hmodel.C
	    
	    //h3 = psr[p].param[param_h3].val[0];
	    h3 = getParameterValue( &psr[p], param_h3, 0 );
	    
	    // Determine fw10 mode
	    if( psr[p].param[param_h4].paramSet[0] == 1 ){
	      //	      h4 = psr[p].param[param_h4].val[0];
	      h4 = getParameterValue( &psr[p], param_h4, 0 );
	      // mode 2 or 3 take preference over mode 1 as they are more stable
	      if( psr[p].param[param_nharm].paramSet[0] == 1 ){
		nharm = (int)psr[p].param[param_nharm].val[0];
		//nharm = (int)getParameterValue( &psr[p], param_nharm, 0 );
		if( nharm > 4 )
		  mode = 3;
		else
		  mode = 2;
	      }
	      if( psr[p].param[param_stig].paramSet[0] == 1 ){
		// Conflict. Unsure whether to select mode 1 or modes 2/3, so will default
		// to the most stable one.
		printf( "WARNING! You specified both H4 and STIG.\n" );
		printf( "We will ignore STIG and perform the approx. H4 fit instead.\n" );
		printf( "If you want to perform the exact fit for H3 and STIG, then " );
		printf( "please remove H4 from your parameter file.\n");
	      }
	      // Have H3, H4, but no NHARM
	      mode = 2;
	    }else{ 
	      // Have H3, but no H4
	      if( psr[p].param[param_stig].paramSet[0] == 1 ){
		//		stig = psr[p].param[param_stig].val[0];
		stig = getParameterValue( &psr[p], param_stig, 0 );
		mode = 1;
	      }else{
		mode = 0;
		h4 = 0;
		nharm = 3;
	      }
	    }// fw10 mode determined.
	    
	    // Define sin(i) and m2 for calculation of the orbital phases etc.
	    if( mode == 1 ){
	      // fw10, Eq. 22:
	      si = 2.0 * stig / ( 1.0 + pow( stig, 2.0 ) );
	      // fw10, Eq. 20:
	      m2 = h3 / pow( stig, 3.0 ); // Shapiro r, not just M2.
	      
	      if( si > 1.0 ){
		displayMsg(1,"BIN1",
			   "SIN I > 1.0, setting to 1: should probably use DDS model",
			   "",psr[p].noWarnings);
		si = 1.0;
		psr[p].param[param_sini].val[0] = 1.0L;
	      }
	    }else if( mode == 2 || mode == 3 ){
	      // fw10, Eq. 25:
	      si = 2.0 * h3 * h4 / ( h3 * h3 + h4 * h4 );
	      // fw10, Eq. 26:
	      m2 = pow( h3, 4.0 ) / pow( h4, 3.0 );
	      if( si > 1.0 ){
		displayMsg(1,"BIN1",
			   "SIN I > 1.0, setting to 1: should probably use DDS model",
			   "",psr[p].noWarnings);
		si = 1.0;
		psr[p].param[param_sini].val[0] = 1.0L;
	      }
	    }else if( mode == 0 ){
	      // Cannot determine m2 and/or sini. Will have to determine the
	      // Shapiro delay based on h3 alone.
	    }else{
	      printf( "This should not be possible. Go Away.\n" );
	      printf( "And tell someone about it: joris.verbiest@gmail.com, e.g.\n" );
	    }
	    brace=1-si*sin(phase);
	    dlogbr=log(brace);
	    
	    ecc = sqrt( eps1 * eps1 + eps2 * eps2 );
	    TrueAnom = phase;
	    //TrueAnom = 2.0 * atan2( sqrt( 1.0 + ecc ) * sin( phase / 2.0 ), 
	    //                       sqrt( 1.0 - ecc ) * cos( phase / 2.0 ) );
	    omega = atan2( eps1, eps2 );
	    //lgf = log( 1.0 + stig * stig - 2.0 * stig * sin( TrueAnom + omega ) );

	    fs = 1.0 + stig * stig - 2.0 * stig * sin( TrueAnom );
	    lgf = log( fs );
	    lsc = lgf + 2.0 * stig * sin( TrueAnom ) - stig * stig * cos( 2.0 * TrueAnom );
	    
	    if( mode == 0 ){
	      // mode 0: only h3 is known. 
	      ds = -4.0 / 3.0 * h3 * sin( 3.0 * TrueAnom );
	    }else if( mode == 1 ){
	      ds = -2.0 *m2* lsc;
	      //ds = -2.0 * m2 * dlogbr;
	    }else{ // modes 2 and 3
	      ds = calcDH( TrueAnom, h3, h4, nharm, 0 );
	    }
	    
	  }
	  else{
	    dlogbr=log(brace);
	    ds=-2*m2*dlogbr;        /* Equation 26 */
	  }
	  
        }
      else
        {
          printf("Require eccentricity set or EPS1/EPS2 parameters for companion %d\n",com+1);
          exit(1);
        }
      // printf("T2: %g %g %g %g %g %g %g\n",brace,phase,a0,b0,dlogbr,ds,m2);

      /* Now compute d2bar, the orbital time correction in DD equation 42. */
      /* Equation 52 */
      if (onemecu != 0.0)
        {
          d2bar=dre*(1-anhat*drep+allTerms*pow(anhat,2)*
                     (pow(drep,2) + 0.5*dre*drepp - 
                      0.5*ecc*su*dre*drep/onemecu)) + allTerms*(ds+da+DAOP+DSR
                                                                + DOP);
        }
      else
        {
          d2bar=dre*(1-anhat*drep+allTerms*pow(anhat,2)*
                     (pow(drep,2) + 0.5*dre*drepp))
            + allTerms*(ds+da+DAOP+DSR+DOP);
        }    
      //      printf("T2a: %g %g %g %g %g drepp=%g ecc=%g su =%g ome =%g ds = %g da = %g %g %g\n",(double)d2bar,(double)dre,(double)anhat,(double)drep,(double)allTerms,(double)drepp,(double)ecc,(double)su,(double)onemecu,(double)ds,(double)da,(double)DAOP,(double)DSR);
      torb-=d2bar;                                  /* Equation 42  */

      if (param==-1 && com == psr[p].nCompanion-1) return torb;
      else if (param!=-1 && com==arr)
        {
          // Now we need the partial derivatives. Use DD equations 62a - 62k.
          if (psr[p].param[param_ecc].paramSet[com]==1)
            {
              csigma=x*(-sw*su+sqr1me2*cw*cu)/onemecu;      /* Equation 62a */
              ce=su*csigma-x*sw-ecc*x*cw*su/sqr1me2;        /* Equation 62c */
              cx=sw*cume+sqr1me2*cw*su;                     /* Equation 62d */
              comega=x*(cw*cume-sqr1me2*sw*su);             /* Equation 62e */
              cgamma=su;                                    /* Equation 62g */
              cdth=-ecc*ecc*x*cw*su/sqr1me2;                /* Equation 62i */
              cm2=-2*dlogbr;                                /* Equation 62j */
              if (psr[p].param[param_shapmax].paramSet[com]==1)
                cshapmax = 2*m2*(sw*cume+sqr1me2*cw*su)/brace * (1.0-sdds);
              else if (psr[p].param[param_sini].nLinkTo>0){
                csi= 2*m2*(sw*cume+sqr1me2*cw*su)*cos(ki)/brace;
              }
              else
                csi=2*m2*(sw*cume+sqr1me2*cw*su)/brace;     /* Equation 62k */
            }
          else if (psr[p].param[param_eps1].paramSet[com]==1) /* ELL1 model */
            {
              csigma = x*cos(phase);
              cx     = sin(phase);
              ceps1  = -0.5*x*cos(2*phase);
              ceps2  =  0.5*x*sin(2*phase);
              cm2    = -2*dlogbr;
              csi    = 2*m2*sin(phase)/brace;
            }

          if (param==param_pb)	         return -csigma*an*SECDAY*tt0/pb; 
          else if (param==param_a1)      return cx;
          else if (param==param_ecc)     return ce;
          else if (param==param_edot)     return ce*tt0;
          else if (param==param_om)      return comega;
          else if (param==param_omdot)   
            return ae*comega/(an*360.0/(2.0*M_PI)*365.25*SECDAY);
          else if (param==param_t0)      return -csigma*an*SECDAY;
          else if (param==param_pbdot){
            if(psr[p].param[param_pbdot].nLinkFrom>0){
              return 0.5*tt0*(-csigma*an*SECDAY*tt0/(pb*SECDAY));
              /*- SPEED_LIGHT/(getParameterValue(&psr[p],param_pb,0)*SECDAY*
                (pow(getParameterValue(&psr[p],param_pmra,0)*MASYR2RADS,2)+
                pow(getParameterValue(&psr[p],param_pmdec,0)*MASYR2RADS,2))*
                getParameterValue(&psr[p],param_daop,0))*
                (C*(-DK011-DK012)+S*(DK021+DK022));*/
            }
            else  return 0.5*tt0*(-csigma*an*SECDAY*tt0/(pb*SECDAY));
          }
          else if (param==param_sini)    return csi;
          else if (param==param_gamma)   return cgamma;
          else if (param==param_m2)      return cm2*SUNMASS;
          else if (param==param_a1dot)   return cx*tt0;
          else if (param==param_eps1)    return ceps1;
          else if (param==param_eps1dot) return ceps1*tt0;
          else if (param==param_eps2dot) return ceps2*tt0;
          else if (param==param_eps2)    return ceps2;
          else if (param==param_tasc)    return -csigma*an*SECDAY;
          else if (param==param_shapmax) return cshapmax;
	  else if (param==param_stig){
	    return( -2.0 * m2 / stig * ( 1.0 - 3.0 * lgf - ( 1.0 - stig * stig ) / fs ) 
		    + 2.0 * m2 * ( 4.0 * sin( TrueAnom ) - stig * cos( 2.0 * TrueAnom ) ) );	    
	  }
	  else if (param==param_h3) {
	    if( mode == 0 || mode == 2)
	      return( -4.0 / 3.0 * sin( 3.0 * TrueAnom ) );
	    else if( mode == 1 ){
	      return( -2.0 * lsc / pow( stig, 3.0 ) );
	      //return( 2.0 / pow( stig, 3.0 ) * 
	      //        ( lgf + 2.0 * stig * sin( TrueAnom ) - stig * stig * cos( TrueAnom ) ) );
	    }else if( mode == 3 )
	      return( calcDH( TrueAnom, h3, h4, nharm, 3 ) );
	    else{
	      printf( "ERROR in ELLH model in T2. This really shouldn't happen.\n" );
	    }
	  }else if( param == param_h4 ){
	    if( mode == 2 )
	      return( cos( 4.0 * TrueAnom ) );
	    else
	      return( calcDH( TrueAnom, h3, h4, nharm, 4 ) );
	  }	  
          else if (param==param_kom){
            ckom = C* (DK033+DK034+DK013+DK014) + S*(DK043+DK044+DK023+DK024);
            return ckom;
          }
          else if (param==param_kin){
            ckin = C/sin(ki)*(DK043+DK044+DK023+DK024)-
              S/sin(ki)*(DK013+DK014+DK033+DK034);
            if( psr[p].param[param_ecc].paramSet[com]== 1 )
                ckin += dpara/AULTSC/2.0*pow( x, 2.0 ) * cos( ki ) *
                pow( sin( ki ), -3.0 ) * 
                ( pow( ecc, 2.0 ) * ( 3.0 - cos( 2.0 * u ) ) + 
                  4.0 * ecc * cume - 2.0 );
            if(psr[p].param[param_kin].nLinkFrom>0)
              ckin += csi; //  * cos( ki );
            return ckin;
          }
          /* Update the binary parameter jumps */
          if (psr[p].param[param_bpjep].paramSet[arr]==1 && 
              psr[p].obsn[ipos].bbat > psr[p].param[param_bpjep].val[arr])
            {
              if (param==param_bpjph)
                return 1.0/psr[p].param[param_f].val[0];
              else if (param==param_bpja1) return cx;
              else if (param==param_bpjec) return ce;
              else if (param==param_bpjom) return comega;
              else if (param==param_bpjpb)
                return -csigma*an*SECDAY*tt0/(pb*SECDAY); 
            }
          else
            return 0.0;
        }
    }
  return 0.0;
}


void updateT2(pulsar *psr,double val,double err,int pos,int arr){
  if (pos==param_pb || pos==param_a1 || pos==param_ecc || pos==param_t0 || 
      pos==param_sini || pos==param_m2 || pos == param_gamma || 
      pos==param_eps1 || pos==param_eps2 || pos==param_tasc ||
      pos == param_bpjph || pos==param_bpja1 || pos==param_bpjec || 
      pos==param_bpjom || pos == param_bpjpb || pos==param_shapmax || pos==param_h3 || pos==param_h4 || pos==param_stig || pos==param_edot){
    psr->param[pos].val[arr] += val;
    psr->param[pos].err[arr]  = err;
  }
  else if (pos==param_om || pos==param_kom || pos==param_kin)
    {
      psr->param[pos].val[arr] += val*180.0/M_PI;
      psr->param[pos].err[arr]  = err*180.0/M_PI;
    }
  else if (pos==param_pbdot)
    {
      psr->param[pos].val[arr] += val;
      psr->param[pos].err[arr]  = err;
    }
  else if (pos==param_a1dot || pos == param_eps1dot || pos==param_eps2dot)
    {
      psr->param[pos].val[arr] += val;
      psr->param[pos].err[arr]  = err;
    }
  else if (pos==param_omdot)
    {
      psr->param[pos].val[arr] += val; /* *(SECDAY*365.25)*180.0/M_PI; */
      psr->param[pos].err[arr]  = err; /* *(SECDAY*365.25)*180.0/M_PI; */
    }
}

long double getParameter(pulsar *psr,int p,int k)
{  
  if (k > psr->param[p].aSize) return 0.0;
  //  if (psr->param[p].paramSet[k]==1) return(double)psr->param[p].val[k];

 
  if (psr->param[p].paramSet[k]==1){
    if(psr->param[p].nLinkTo>0)
      return getParameterValue(psr,p,k);

    else 
      return psr->param[p].val[k];
  }
 
  return 0.0;
}


/* Given system masses of m,m2 and Keplerian parameters x,ecc and an, this 
 * routine calculates values of dr,dth,er,eth, si, gamma and pbdot under GR */

void calcGR(double mtot,double m2,double x,double ecc,double an,double afac,
            double f0, double *dr,double *dth,double *er,double *eth,
            double *xk,double *si,double *gamma, double *pbdot,double *a0,
            double *b0){
  double ARRTOL = 1.0e-10;
  double m1,arr0,arrold;
  double arr,ar;
  double a0aligned;

  m1 = mtot-m2;

  if (mtot<0)
    {
      printf("ERROR: problem in DDGR model (mtot < 0)\n");
      exit(1);
    }
  arr0 = pow(mtot/(an*an),1.0/3.0);
  arr  = arr0;
  do {
    arrold = arr;
    arr = arr0*pow(1.0+(m1*m2/pow(mtot,2) - 9.0)*0.5*mtot/arr,2.0/3.0);
  } while (fabs((arr-arrold)/arr) > ARRTOL);


  arr = arr0*pow(1.0+(m1*m2/pow(mtot,2) - 9.0)*0.5*mtot/arr,2.0/3.0);
  ar  = arr*m2/mtot;


  *si=x/ar;
  *xk=3.0*mtot/(arr*(1.0-ecc*ecc));
  *gamma = ecc*m2*(m1+2*m2)/(an*arr*mtot);
  *pbdot = -(96.0*2.0*M_PI/5.0)*pow(an,5.0/3.0)*pow(1.0-pow(ecc,2),-3.5)
    * (1+(73.0/24)*pow(ecc,2) + (37.0/96)*pow(ecc,4))
    * m1*m2*pow(mtot,-1.0/3.0);

  *dr  = (3.0*pow(m1,2)+6.0*m1*m2 + 2.0*pow(m2,2))/(arr*mtot);
  *er  = ecc*(1.0+(*dr)); 
  // Shouldn't this be "1.0+(*dr)"??? instead of "1.0*(*dr)".

  *dth = (3.5*m1*m1 + 6*m1*m2 + 2*m2*m2)/(arr*mtot);
  *eth = ecc*(1.0+(*dth));

  a0aligned = an*ar/(2.0*M_PI*f0*(*si)*(sqrt(1-ecc*ecc)));
  *a0 = afac*a0aligned;
  *b0 = 0;
}

/*
 *
 * pb  = orbital period in seconds 
 * t0  = epoch of periastron (MJD)
 * ecc = eccentricity
 * omz = omega (deg)
 * x   = a1 (lt-s)
 *
 * eps1= ELL1 binary model parameter 1 
 * eps2= ELL1 binary model parameter 2 
 * tasc= Time of ascending node
 */

void getKeplerian(pulsar *psr,int com,double *pb,long double *t0,double *ecc,
                  double *omz,double *x,double *eps1,double *eps2,
                  long double *t0asc,double *shapmax,double *kom,double *kin)
{
  *pb  = getParameter(psr,param_pb,com)*SECDAY;
  *t0  = getParameter(psr,param_t0,com);
  *ecc = getParameter(psr,param_ecc,com); 
  *omz = getParameter(psr,param_om,com);
  *x   = getParameter(psr,param_a1,com);
  
  *eps1= getParameter(psr,param_eps1,com);
  *eps2= getParameter(psr,param_eps2,com);
  *t0asc = getParameter(psr,param_tasc,com);
  
  *shapmax = getParameter(psr,param_shapmax,com);

  *kom = getParameter(psr,param_kom,com)*M_PI/180.0;
  *kin = getParameter(psr,param_kin,com)*M_PI/180.0;
}

/*
 * Following the BTJ model, this function adds jumps to the Keplerian parameters
 * at a specified epoch
 *
 */
void addKeplerianJumps(pulsar *psr,int ipos,double *torb,double *x,double *ecc,
                       double *omz,double *pb)
{
  int i;

  for (i=0;i<psr->param[param_bpjep].aSize;i++)
    {
      if (psr->param[param_bpjep].paramSet[i]==1 && 
          psr->obsn[ipos].bbat > psr->param[param_bpjep].val[i])
        {
          *torb = *torb - (double)(psr->param[param_bpjph].val[i]
                                   / psr->param[param_f].val[0]);  
          *x    = *x    + (double)psr->param[param_bpja1].val[i];
          *ecc  = *ecc  + (double)psr->param[param_bpjec].val[i];
          *omz  = *omz  + (double)psr->param[param_bpjom].val[i];    
          *pb   = *pb   + (double)psr->param[param_bpjpb].val[i]*SECDAY; 
        }
    } 
}

/* Post-Keplerian parameters 
 *
 * si    = sine of inclination angle
 * m2    = companion mass (kg)
 * omdot = rate of periastron advance (orbits/s)
 * gamma = gamma term
 * xdot  = rate of change of projected semi-major axis of orbit
 * pbdot = rate of change of orbital period
 * edot  = rate of change of eccentricity
 * xpbdot= rate of change of orbital period minus GR prediction
 *
 * eps1dot
 * eps2dot
 */

void getPostKeplerian(pulsar *psr,int com,double an,double *si,double *m2,
                      double *mtot,double *omdot, double *gamma,double *xdot,
                      double *xpbdot,double *pbdot, double *edot,double *pmra,
                      double *pmdec,double *dpara, double *dr,double *dth,
                      double *a0,double *b0,double *xomdot,double *afac,
                      double *eps1dot,double *eps2dot, double *daop){
  double SUNMASS = 4.925490947e-6;
  double rad2deg = 180.0/M_PI;
  //double pxConv = 1.74532925199432958E-2/3600.0e3;//converts mas to rad
  double pxConv = M_PI/180.0/3600*1e-3; // converts mas to rad
  double daopConv = 3.08568025e16;//pc in m

  //  logdbg("Going to get parameters");
  *si      = getParameter(psr,param_sini,com);
  if (*si > 1.0)
    {
      displayMsg(1,"BIN1","SIN I > 1.0, setting to 1: should probably use DDS model","",psr[0].noWarnings);
      *si = 1.0;
      psr[0].param[param_sini].val[0] = 1.0;
    }
  if (*si < -1.0)
    {
      displayMsg(1,"BIN1","SIN I < -1.0, setting to -1: should probably use DDS model","",psr[0].noWarnings);
      *si = -1.0;
      psr[0].param[param_sini].val[0] = -1.0;
    }
  *m2      = getParameter(psr,param_m2,com)*SUNMASS;
  *mtot    = getParameter(psr,param_mtot,com)*SUNMASS;
  *omdot   = getParameter(psr,param_omdot,com)/(rad2deg*365.25*SECDAY*an);
  *gamma   = getParameter(psr,param_gamma,com);
  *xdot    = getParameter(psr,param_a1dot,com);
  *xpbdot  = getParameter(psr,param_xpbdot,com);
  *pbdot   = getParameter(psr,param_pbdot,com);
  *edot    = getParameter(psr,param_edot,com);
  *pmra    = getParameter(psr,param_pmra,com)
    * M_PI/(180.0*3600.0e3)/(365.25*86400.0);
  *pmdec   = getParameter(psr,param_pmdec,com)
    * M_PI/(180.0*3600.0e3)/(365.25*86400.0);
  *dpara   = getParameter(psr,param_px,com)*pxConv;
  *dr      = getParameter(psr,param_dr,com);
  *dth     = getParameter(psr,param_dth,com);
  *a0      = getParameter(psr,param_a0,com);
  *b0      = getParameter(psr,param_b0,com);
  *xomdot  = getParameter(psr,param_xomdot,com)/(an*rad2deg*365.25*86400.0);
  *afac    = getParameter(psr,param_afac,com);
  *eps1dot = getParameter(psr,param_eps1dot,com);
  *eps2dot = getParameter(psr,param_eps2dot,com);
  *daop    = getParameter(psr,param_daop,com)*1e-3/pxConv;
}


void updateParameters(double edot,double xdot,double eps1dot,double eps2dot,
                      long double tt0,double *ecc,double *x,double *eps1,
                      double *eps2){
  (*ecc)  += edot*tt0;
  (*x)    += xdot*tt0;
  (*eps1) += eps1dot*tt0;
  (*eps2) += eps2dot*tt0;
}

void deriveKeplerian(double pb,double kom,double *an,double *sin_omega,
                     double *cos_omega){
  *an    = 2.0*M_PI/pb;
  *sin_omega = sin(kom);
  *cos_omega = cos(kom);
}

void derivePostKeplerian(double mtot,double m2,double dr,double dth,
                         double ecc,double *m1,double *er,double *eth)
{
  *m1  = mtot - m2;  /* Pulsar mass */
  *er  = ecc*(1.0+dr);
  *eth = ecc*(1.0+dth);
}

void KopeikinTerms(pulsar *psr,int ipos,double ki,double pmra,double sin_omega,
                   double pmdec, double cos_omega,long double tt0,double dpara, 
                   double daop, double si,double *x, long double *DK011, 
                   long double *DK012, long double *DK021,long double *DK022, 
                   long double *DK031,long double *DK032, long double *DK041, 
                   long double *DK042,long double *DK013, long double *DK014, 
                   long double *DK023, long double *DK024, long double *DK033, 
                   long double *DK034, long double *DK043, long double *DK044){
  double ki_dot,sini,cosi,tani;
  double sin_delta,cos_delta,sin_alpha,cos_alpha;
  double delta_i0,delta_j0,xpr,ypr;
  si = sin(ki);
  /* Equation 10 in Kopeikin 1996 */

  //  ki    += ki_dot*tt0;
  sini = sin(ki);
  cosi = cos(ki);
  tani = sini/cosi;
  ki_dot = -pmra * sin_omega + pmdec*cos_omega;  
  /* Equation 8 in Kopeikin 1996 */
  //  (*x) += ((*x)*ki_dot/tani)*tt0;
  /* Equation 9 in Kopeikin 1996 */
  //(*omz) += (pmra*cos_omega+pmdec*sin_omega)/sini*tt0;

  /* Now modify x and omega due to the annual-orbital parallax term 
   * as described in Kopeikin 1995 
   *
   * Require knowledge of the barycentric earth position vector - earth_ssb
   */

  /* Obtain vector pointing at the pulsar */
  sin_delta = psr->obsn[ipos].psrPos[2];
  cos_delta = cos(asin(sin_delta));
  sin_alpha = psr->obsn[ipos].psrPos[1]/cos_delta;
  cos_alpha = psr->obsn[ipos].psrPos[0]/cos_delta;
  
  /* Equation 15 in Kopeikin 1995 */
  delta_i0 = -psr->obsn[ipos].earth_ssb[0]/AULTSC*sin_alpha+
    psr->obsn[ipos].earth_ssb[1]/AULTSC*cos_alpha;
  /* Equation 16 in Kopeikin 1995 */
  delta_j0 = -psr->obsn[ipos].earth_ssb[0]/AULTSC*sin_delta*cos_alpha-
    psr->obsn[ipos].earth_ssb[1]/AULTSC*sin_delta*sin_alpha+
    psr->obsn[ipos].earth_ssb[2]/AULTSC*cos_delta;
  
  xpr = delta_i0*sin_omega - delta_j0*cos_omega;
  ypr = delta_i0*cos_omega + delta_j0*sin_omega;
  
  /* Equations 18 and 19 in Kopeikin 1995 */

  if(psr->param[param_daop].paramSet[0]==1){
    //    logdbg("Using daop for par file for Kopeikin delays!");
      
    *DK011 = (longdouble)(-(*x)/daop/si*delta_i0*sin_omega);
    *DK012 = (longdouble)(-(*x)/daop/si*delta_j0*cos_omega);
    *DK013 = (longdouble)(-(*x)/daop/si*delta_i0*cos_omega);
    *DK014 = (longdouble)((*x)/daop/si*delta_j0*sin_omega);
      
    *DK021 = (longdouble)((*x)/daop/tani*delta_i0*cos_omega);
    *DK022 = (longdouble)(-(*x)/daop/tani*delta_j0*sin_omega);
    *DK023 = (longdouble)((*x)/daop/tani*delta_i0*sin_omega);
    *DK024 = (longdouble)((*x)/daop/tani*delta_j0*cos_omega);
  }
  else{
    *DK011 = (longdouble)(-(*x)*dpara/si*delta_i0*sin_omega);
    *DK012 = (longdouble)(-(*x)*dpara/si*delta_j0*cos_omega);
    *DK013 = (longdouble)(-(*x)*dpara/si*delta_i0*cos_omega);
    *DK014 = (longdouble)((*x)*dpara/si*delta_j0*sin_omega);

    *DK021 = (longdouble)((*x)*dpara/tani*delta_i0*cos_omega);
    *DK022 = (longdouble)(-(*x)*dpara/tani*delta_j0*sin_omega);
    *DK023 = (longdouble)((*x)*dpara/tani*delta_i0*sin_omega);
    *DK024 = (longdouble)((*x)*dpara/tani*delta_j0*cos_omega);
  }  
  *DK031 = (longdouble)((*x)*tt0/si*pmra*sin_omega);
  *DK032 = (longdouble)((*x)*tt0/si*pmdec*cos_omega);
  *DK033 = (longdouble)((*x)*tt0/si*pmra*cos_omega);
  *DK034 = (longdouble)(-(*x)*tt0/si*pmdec*sin_omega);

  *DK041 = (longdouble)((*x)*tt0/tani*pmra*cos_omega);
  *DK042 = (longdouble)(-(*x)*tt0/tani*pmdec*sin_omega);
  *DK043 = (longdouble)(-(*x)*tt0/tani*pmra*sin_omega);
  *DK044 = (longdouble)(-(*x)*tt0/tani*pmdec*cos_omega);
    
    
  //  logdbg("DK011 %g, DK021 %g, x %f, dpara %g, ypr %f and si %f and tani %f, xpr %f",
  //                         (double)(*DK011),(double)(*DK021),(double)(*x),
  //                         (double)(dpara),(double)(ypr),(double)(si),
  //                         (double)(tani),(double)(xpr));
}  


/*  Compute eccentric anomaly u by iterating Kepler's equation if
    eccentricity is set. The equation is solved using a Newton-Raphson
    technique and the S9 * starting value in Odell & Gooding 1986
    CeMec 38 307 */
void computeU(double phase,double ecc,double *u)
{
  double du;

  /*	  u=phase+ecc*sin(phase)*(1.0+ecc*cos(phase));*/
  *u = phase+ecc*sin(phase)/sqrt(1.0-2*ecc*cos(phase)+ecc*ecc);
  do {
    du=(phase-(*u-ecc*sin(*u)))/(1.0-ecc*cos(*u));
    (*u)+=du;
  } while (fabs(du)>1.0e-14);
}

/* Based on the DDT model and the equations in Damour & Taylor (1992) */
/* DDT model includes a loop around d2bar --- MUST INCLUDE */

/* void useBeta(pulsar psr,int iteration)
   {
   int i;
   //  double TSUN=4.925490947e-6; // should be defined in tempo2.h

   am1z=am;
   am2z=am2;
   deltam1=1.0e-8;
   deltam2=1.0e-8;

   am1=am1z;
   am2=am2z;
   if(iteration==1) am1=am1z+deltam1;
   if(iteration==2) am2=am2z+deltam2;
  
   // Compute quantities depending on M, m2, bp, and bpp. 
   am=am1+am2;
   m=am*TSUN;
   m1=am1*TSUN;
   m2=am2*TSUN;
   aa=2.1569176;
   bb=1.0261529;
   cp=0.21;
   c1=cp*am1;
   c2=cp*am2;
   a1a2=0.5*bp*bb*(c1*c1 + c2*c2);
   bk1=-c2-bb*c1*c1 + (aa-3*bb)*c2*c2 - 2*(aa-bb)*c1*c1 * c2 +
   (2*aa*aa -7*aa*bb + 5*bb*bb)*c1*c1 * c2*c2;
   bk2=-3*c2*c2*c2 + 2*c1*c1 * c2*c2 + c2*c2*c2*c2 + 0.5*c1*c1*c1*c1 * c2 + 
   aa*c1*c1*c1*c1 * c2*c2;
   a1b2a1=bp*bk1 + pow(bp*bb,2) * bk2 + 0.5*bpp*bb*c2*c1;
   bk3=-c1-bb*c2*c2 + (aa-3*bb)*c1*c1 - 2*(aa-bb)*c2*c2 * c1 +
   (2*aa*aa -7*aa*bb + 5*bb*bb)*c2*c2 * c1*c1;
   bk4=-3*c1*c1*c1 + 2*c2*c2 * c1*c1 + c1*c1*c1*c1 + 0.5*c2*c2*c2*c2 * c1 + 
   aa*c2*c2*c2*c2 * c1*c1;
   a2b1a2=bp*bk3 + pow(bp*bb,2) * bk4 + 0.5*bpp*bb*c1*c1;
   a0a2=0.5*bp*bb*c2*c2;

   brk1=(1.d0-a1a2/3.d0)*pow(1.d0+a1a2,-1.0/3.0) - 
   (1.0/(6.0*am)) *
   (am1*a1b2a1+am2*a2b1a2)*pow(1.d0+a1a2,-4.0/3.0);
   k=(3.d0/(1.d0-ecc*ecc)) * pow(an*TSUN*am,2.0/3.0) * brk1;
   brk2=1.d0 + (am2/am)*(1.0+a1a2)+a0a2;
   gamma=(ecc/an)*pow(an*TSUN*am,2.0/3.0) * (am2/am) *
   pow(1.0+a1a2,-1.0/3.0) * brk2;
   si=x*pow(an,2.d0/3.d0) * pow(TSUN*am,-1.0/3.0) * (am/am2) *
   pow(1.d0+a1a2,-1.0/3.0);

   arr=pow(m/an*an,1.0/3.0);
   ar=arr*m2/m;

   //  DD equations 36, 37: 
   dr=(3*m1*m1 + 6*m1*m2 + 2*m2*m2)/(arr*m);
   er=ecc*(1+dr);
   dth=(3.5d0*m1*m1 + 6*m1*m2 + 2*m2*m2)/(arr*m);
   eth=ecc*(1+dth);
   }
*/

static double calcDH( double ae, double h3, double h4, int nharm, int sel ){
  // The final input argument "sel", selects the output for this function.
  // sel = 3 outputs the derivative for H3;
  // sel = 4 outputs the derivative for H4;
  // sel = 0 outputs the shapiro delay.
  double s = h4 / h3;
  double FirstH3 = -4.0 / 3.0 * sin( 3.0 * ae );
  double SecondH3 = 0.0;
  double FirstH4 = cos( 4.0 * ae );
  double SecondH4 = 0.0;
  double sd3 = -4.0 / 3.0 * h3 * sin( 3.0 * ae );
  double sd4 = h4 * cos( 4.0 * ae );
  double sd5 = 0.0;

  double fs = s * sin( 5.0 * ae ) / 5.0;
  double dfsds = sin( 5.0 * ae ) / 5.0;
  double count;

  if( nharm > 5)
    for( int EvenCt = 6; EvenCt <= nharm; EvenCt += 2 ){
      count = (double)EvenCt;
      // Add even harmonics
      fs += pow( -1.0, count / 2.0 )/
        count * pow( s, count - 4.0 ) * 
        cos( count * ae );
      dfsds += pow( -1.0, count / 2.0 ) / 
        count * pow( s, count - 5.0 ) *
        cos( count * ae ) * ( count - 4.0 );
    }
  
  if( nharm > 6 )
    for( int OddCt = 7; OddCt <= nharm; OddCt += 2 ){
      count = (double)OddCt;
      // Add odd harmonics
      fs += pow( -1.0, ( count - 1.0 ) / 2.0 ) / 
        count * pow( s, count - 4.0 ) * 
        sin( count * ae );                                    
      dfsds += pow( -1.0, ( count - 1.0 ) / 2.0 ) / 
        count * pow( s, count - 5.0 ) * 
        sin( count * ae ) * ( count - 4.0 );
    }

  if( nharm > 4 ){
    SecondH3 = -4.0 * dfsds * s * s;
    SecondH4 = 4.0 * ( fs + dfsds * s );
    sd5 = 4.0 * h4 * fs;
  }

  if( sel == 3 )
    return( FirstH3 + SecondH3 );
  else if( sel == 4 )
    return( FirstH4 + SecondH4 );
  else if( sel == 0 )
    return( sd3 + sd4 + sd5 );
  else{
    printf( "ERROR in ELL1Hmodel! This shouldn't be happening.\n" );
  }
}
