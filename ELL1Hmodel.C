#include <stdio.h>
#include <math.h>
#include "tempo2.h"

/* ------------------------------------------------------------------------- */
/*  Timing model for small-eccentricity binary pulsars, e<<1 (Wex 1998)      */
/*  Expanded with harmonic Shapiro delay fitting (Freire & Wex, 2010) 
    by JPWV in Parkes observatory, May 2011.                                 */
/*                                                                           */
/*  Instead of e and omega the Laplace parameters                            */
/*     epsilon1 = e*sin(omega)	                                             */
/*     epsilon2 = e*cos(omega)                                               */
/* are used as new parameters. T0 is related to the ascending node (not to   */
/* periastron as in BT, DD, ...)                                             */
/*                                                                           */
/*  Time derivatives:                                                        */
/*     nell1=0 -> fit for eps1dot,eps2dot                                    */
/*     nell1=1 -> fit for omdot,edot                                         */
/*                                                                           */
/*  Computes pulsar orbit time, torb, at time of observation t=ct(n)-pepoch. */
/*  Pulsar proper time is then TP=T+TORB.                                    */
/*  Units are such that c=G=1. Thus masses have units of seconds, with       */
/*  one solar mass = 4.925490947 usec.                                       */
/*                                                                           */
/*  Also computes the binary orbit-related values of fctn: partial           */
/*  derivatives of each arrival time residual with respect to the model      */
/*  parameters.                                                              */
/*                                                                           */
/*  Based on bnryell1.f                                                      */
/* ------------------------------------------------------------------------- */

static double calcDH( double ae, double h3, double h4, int nharm, int sel);

double ELL1Hmodel(pulsar *psr,int p,int ipos,int param){
  double an; // orbital angular velocity
  double x0,m2,tt0,orbits,phase,e1,e2,dre,drep,drepp,brace,dlogbr,ds,da,pb;
  double eps1,eps2,eps1dot,eps2dot,si,a0,b0,d2bar;
  double torb,Csigma,Cx,Ceps1,Ceps2,Cm2,Csi,ct,t0asc,pbdot,xpbdot,x,xdot,am2;
  int norbits;
  double SUNMASS = 4.925490947e-6;
  double ecc; // eccentricity
  double TrueAnom; // true anomaly
  double omega; // omega
  // fw10 parameters:
  double h3, h4, stig;
  double lgf; // intermediate derivative calculation
  int nharm=4;
  int mode = -1; // How to implement the fw10 model. There are four modes:
  // mode 0: Only h3 is fitted. This is an incomplete parameterisation of the 
  //         Shapiro delay.
  // mode 1: h3 and stigma are fitted. This is a complete and exact parameterisation 
  //         of the Shapiro delay and should only be used if the SD is well 
  //         constrained.
  // mode 2: h3 and h4 fitted. This is an approximate parameterisation of the SD
  //         and can be translated in m2 and sini.
  // mode 3: h3, h4, nharm (>=5). This determines h3 and h4 and in doing so, 
  //         takes nharm harmonics into account. In the lower limit (nharm = 4), 
  //         this mode is identical to mode 2. In the higher limit (nharm->infty),
  //         this mode approaches mode 1. NHARM is a constant, integer input. 
  //         Be aware that nharm>=5 can only be used in case h3 is clearly 
  //         inconsistent with zero.
  
  a0 = 0.0; // WHAT SHOULD THESE BE? -- They should be zero because
  b0 = 0.0; // the Aberration delay is not yet implemented. (See below;
            // JPWV; 15 Jan 2010.)
  
  pb    = psr[p].param[param_pb].val[0]*SECDAY;
  
  if( psr[p].param[param_pbdot].paramSet[0] == 1 ) 
    pbdot = psr[p].param[param_pbdot].val[0];
  else 
    pbdot=0.0;
  
  an    = 2.0*M_PI/pb;

  // Sanity check
  if( psr[p].param[param_h3].paramSet[0] != 1 ){
    printf( "ERROR! Cannot use the ELL1H model without H3.\n" );
  }else
    h3 = psr[p].param[param_h3].val[0];
    //h3 = getParameterValue( &psr[p], param_h3, 0 );

  // Determine fw10 mode
  if( psr[p].param[param_h4].paramSet[0] == 1 ){
    h4 = psr[p].param[param_h4].val[0];
    //    h4 = getParameterValue( &psr[p], param_h4, 0 );
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
      stig = psr[p].param[param_stig].val[0];
      //stig = getParameterValue( &psr[p], param_stig, 0 );
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

  
  x0    = psr[p].param[param_a1].val[0];
  if( psr[p].param[param_a1dot].paramSet[0] == 1 ) 
    xdot  = psr[p].param[param_a1dot].val[0];
  else 
    xdot = 0.0;

  t0asc = psr[p].param[param_tasc].val[0];
  
  xpbdot = 0.0;
  eps1  = psr[p].param[param_eps1].val[0];
  eps2  = psr[p].param[param_eps2].val[0];
  if( psr[p].param[param_eps1dot].paramSet[0] == 1 ) 
    eps1dot = psr[p].param[param_eps1dot].val[0];
  else 
    eps1dot=0;

  if( psr[p].param[param_eps2dot].paramSet[0] == 1 ) 
    eps2dot = psr[p].param[param_eps2dot].val[0];
  else 
    eps2dot=0;

  ct = psr[p].obsn[ipos].bbat;      
  tt0 = (ct-t0asc)*SECDAY;
  orbits = tt0/pb-0.5*(pbdot+xpbdot)*pow(tt0/pb,2);
  norbits = (int)orbits;
  if( orbits < 0.0 )
    norbits = norbits-1;

  phase = 2.0*M_PI*(orbits-norbits);
  
  x = x0+xdot*tt0;
  
  e1 = eps1+eps1dot*tt0;
  e2 = eps2+eps2dot*tt0;
  dre  = x*(sin(phase)-0.5*(e1*cos(2.0*phase)-e2*sin(2.0*phase)));
  drep = x*cos(phase);
  drepp=-x*sin(phase);

  brace = 1.0 - si * sin( phase );
  dlogbr = log( brace );
  
  ecc = sqrt( e1 * e1 + e2 * e2 );
  TrueAnom = phase;
  //TrueAnom = 2.0 * atan2( sqrt( 1.0 + ecc ) * sin( phase / 2.0 ), 
  //                       sqrt( 1.0 - ecc ) * cos( phase / 2.0 ) );
  omega = atan2( e1, e2 );
  //lgf = log( 1.0 + stig * stig - 2.0 * stig * sin( TrueAnom + omega ) );
  double lsc, fs;
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

  /* ================
     ABERRATION DELAY
     ================ */
  /* NOTE: a0 and b0 are always zero -- they are not set in the
     original TEMPO!!!!! */
  da=a0*sin(phase)+b0*cos(phase);  
  
  /*  Now compute d2bar (cf. DD 52) */
  
  // I suspect this following equation gives the expansion of the
  // R\"omer delay, then the Shapiro delay (ds) and finally the
  // aberration delay (da). The Aberration dealy is, however,
  // unmeasurably small in all presently known systems - this is why
  // da == 0, implied by a0 == 0 and b0 == 0. (JPWV; 15 Jan 2010)
  d2bar=dre*(1-an*drep+pow(an*drep,2)+0.5*pow(an,2)*dre*drepp)+ds+da;
  torb=-d2bar;
  if( param == -1 )
    return torb;

  /* Now we need the partial derivatives. */
  Csigma   = x*cos(phase);
  Cx       = sin(phase);
  Ceps1    = -0.5*x*cos(2*phase);
  Ceps2    =  0.5*x*sin(2*phase);
  Cm2      = -2*dlogbr;
  Csi      = 2*m2*sin(phase)/brace; 
  if( param == param_pb ){
    return -Csigma*an*SECDAY*tt0/(pb*SECDAY); /* Pb    */

  }else if( param == param_a1 ){
    return Cx;

  }else if( param == param_eps1 ){
    return Ceps1;

  }else if( param == param_tasc ){
    return -Csigma*an*SECDAY;

  }else if( param == param_eps2 ){
    return Ceps2;

  }else if( param == param_eps1dot ){
    return Ceps1*tt0;

  }else if( param == param_eps2dot ){
    return Ceps2*tt0;

  }else if( param == param_pbdot ){
    return 0.5*tt0*(-Csigma*an*SECDAY*tt0/(pb*SECDAY));

  }else if( param == param_a1dot ){
    return Cx*tt0;  

  }else if( param == param_sini ){
    return Csi;

  }else if( param == param_m2 ){
    return Cm2*SUNMASS;

  }else if( param == param_stig ){
    return( -2.0 * m2 / stig * ( 1.0 - 3.0 * lgf - ( 1.0 - stig * stig ) / fs )
            + 2.0 * m2 * ( 4.0 * sin( TrueAnom ) - stig * cos( 2.0 * TrueAnom ) ) );
    //return( -2.0 * m2 / stig * 
    //        ( 1.0 + 3.0 * lgf - ( 1.0 - stig * stig ) / 
    //          ( 1.0 + stig * stig - 2.0 * stig * sin( TrueAnom + omega ) ) ) );
    //       //+ 2.0 * m2 * ( 4.0 * sin( TrueAnom ) - stig * cos( 2.0 * TrueAnom ) ) );
  }else if( param == param_h3 ){
    if( mode == 0 || mode == 2)
      return( -4.0 / 3.0 * sin( 3.0 * TrueAnom ) );
    else if( mode == 1 ){
      return( -2.0 * lsc / pow( stig, 3.0 ) );
      //return( 2.0 / pow( stig, 3.0 ) * 
      //        ( lgf + 2.0 * stig * sin( TrueAnom ) - stig * stig * cos( TrueAnom ) ) );
    }else if( mode == 3 )
      return( calcDH( TrueAnom, h3, h4, nharm, 3 ) );
    else{
      printf( "ERROR in ELL1H. This really shouldn't happen.\n" );
    }
  }else if( param == param_h4 ){
    if( mode == 2 )
      return( cos( 4.0 * TrueAnom ) );
    else
      return( calcDH( TrueAnom, h3, h4, nharm, 4 ) );
  }

  return 0.0;
}

void updateELL1H(pulsar *psr,double val,double err,int pos){
  if( pos == param_pb ){
    psr->param[param_pb].val[0] += val/SECDAY;
    psr->param[param_pb].err[0]  = err/SECDAY;
  }else if( pos == param_a1 || pos == param_eps1 || pos == param_eps2 ||
            pos == param_tasc || pos == param_sini || pos == param_m2 || 
            pos == param_eps1dot || pos == param_eps2dot || 
            pos == param_h3 || pos == param_h4 || pos == param_stig ){
    psr->param[pos].val[0] += val;
    psr->param[pos].err[0]  = err;
  }else if( pos == param_pbdot ){
    psr->param[pos].val[0] += val;
    psr->param[pos].err[0]  = err;
  }else if( pos == param_a1dot ){
    // JPWV 15 March 2010
    psr->param[pos].val[0] += val;
    psr->param[pos].err[0] += err;
  }else if( pos == param_omdot ){
    // JPWV 15 March 2010
    psr->param[pos].val[0] += val;//*(SECDAY*365.25)*180.0/M_PI;
    psr->param[pos].err[0]  = err;//*(SECDAY*365.25)*180.0/M_PI;
  }
}

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
