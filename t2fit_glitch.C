#include <tempo2.h>
#include <math.h>
#include <assert.h>
#include "enum_str.h"






double t2FitFunc_stdGlitch(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k){

    if (label==param_glph) {
        if (psr[ipsr].obsn[ipos].bbat >= psr[ipsr].param[param_glep].val[k])
            return 1.0/psr[ipsr].param[param_f].val[0];
        else
            return  0.0;
    }
    else if (label==param_glf0d || label==param_glf0d2 || label==param_glf0d3)
    {
        longdouble dt1,expf,tp,tgl;
        param_label our_gltd;
        if (label==param_glf0d) our_gltd = param_gltd;
        if (label==param_glf0d2) our_gltd = param_gltd2;
        if (label==param_glf0d3) our_gltd = param_gltd3;

        tp = (psr[ipsr].obsn[ipos].bbat-psr[ipsr].param[param_pepoch].val[0])*86400.0;
        tgl = (psr[ipsr].param[param_glep].val[k] - psr[ipsr].param[param_pepoch].val[0])*86400.0;

        dt1 = tp-tgl;

        if (psr[ipsr].param[our_gltd].val[k]!=0.0)
            expf = exp(-dt1/86400.0/psr[ipsr].param[our_gltd].val[k]);
        else
            expf = 1.0;

        if (psr[ipsr].obsn[ipos].bbat >= psr[ipsr].param[param_glep].val[k])
        {
            return  psr[ipsr].param[our_gltd].val[k]*SECDAY*(1.0-expf)/psr[ipsr].param[param_f].val[0]; ///psr[ipsr].param[param_f].val[0];
            //	  printf("Glitch diff = %d %.10f %.10f %.10f\n",k+1,afunc,(double)tp,(double)tgl,(double)psr[ipsr].param[param_gltd].val[0]);

        }
        else
            return  0.0;
    }
    else if (label==param_gltd || label == param_gltd2 || label == param_gltd3 )
    {
        longdouble dt1,expf,tp,tgl;
        param_label our_glf0d;
        if (label==param_gltd) our_glf0d=param_glf0d;
        if (label==param_gltd2) our_glf0d=param_glf0d2;
        if (label==param_gltd3) our_glf0d=param_glf0d3;
        param_label our_gltd = label;

        tp = (psr[ipsr].obsn[ipos].bbat-psr[ipsr].param[param_pepoch].val[0])*86400.0L;
        tgl = (psr[ipsr].param[param_glep].val[k] - psr[ipsr].param[param_pepoch].val[0])*86400.0L;

        dt1 = tp-tgl;

        if (psr[ipsr].param[our_gltd].val[k]!=0.0)
            expf = exp(-dt1/86400.0L/psr[ipsr].param[our_gltd].val[k]);
        else
            expf = 1.0;

        if (psr[ipsr].obsn[ipos].bbat >= psr[ipsr].param[param_glep].val[k])
            return  psr[ipsr].param[our_glf0d].val[k]*
                (1.0-(1.0+dt1/SECDAY/(psr[ipsr].param[our_gltd].val[k]))*expf)/psr[ipsr].param[param_f].val[0]*SECDAY;
        else
            return  0.0;
    }
    else if (label==param_glf0)
    {
        if (psr[ipsr].obsn[ipos].bbat >= psr[ipsr].param[param_glep].val[k]){
            return  (psr[ipsr].obsn[ipos].bbat-psr[ipsr].param[param_glep].val[k])*86400.0/psr[ipsr].param[param_f].val[0];
        }
        else
            return  0.0;
    }
    else if (label==param_glf1)
    {
        if (psr[ipsr].obsn[ipos].bbat >= psr[ipsr].param[param_glep].val[k])
            return  0.5*pow((psr[ipsr].obsn[ipos].bbat-psr[ipsr].param[param_glep].val[k])*86400.0,2)/psr[ipsr].param[param_f].val[0];
        else
            return  0.0;
    }
    else if (label==param_glf2)
    {
        if (psr[ipsr].obsn[ipos].bbat >= psr[ipsr].param[param_glep].val[k])
        {
            return  (double)(1.0L/6.0L*powl((psr[ipsr].obsn[ipos].bbat-psr[ipsr].param[param_glep].val[k])*86400.0L/1.0e9,3)/psr[ipsr].param[param_f].val[0]);
        }
        else
            return  0.0;
    }

    logerr("Unknown glitch parameter: label=%d str=%s k=%d",(int)label,label_str[label],k);
    assert(false);
}

void t2UpdateFunc_stdGlitch(pulsar *psr, int ipsr ,param_label label,int k, double val, double error){
    double F=1;
    if(label==param_glf2)F=1e-27;
    psr[ipsr].param[label].val[k] -= val*F;
    psr[ipsr].param[label].err[k]  = error*F;
}


// function for fitting exponential dips (profile events and scattering events)

double t2FitFunc_expdip(pulsar *psr, int ipsr, double x, int ipos, param_label label, int k)
{
  
  
  long double val;

  long double dt;
  long double dm;
  long double tau;
  long double freq;
  long double gamma;
  

  
  // reference to 1.4 GHz  to agree with Enterprise
  freq= psr[ipsr].obsn[ipos].freqSSB/1.4e9;
  
  dt=(psr[ipsr].obsn[ipos].bbat - psr[ipsr].param[param_expep].val[k]);
  dm=psr[ipsr].param[param_expph].val[k];

  tau=psr[ipsr].param[param_exptau].val[k];

  if (psr[ipsr].param[param_expindex].paramSet[k] ==1)
    {
      gamma=psr[ipsr].param[param_expindex].val[k];
    }
  else{
    gamma=-2;
    
  }
  

  if (label ==param_expph)
    {
      if (dt> 0)
	{	   
	  
	  return val =  powl(freq,gamma)*exp(-dt/tau);
	}
      else
	{
	  return 0.0;
	}
    }
  if( label ==param_exptau)
    {
      if (dt > 0)
	{
	  return val = dm*powl(freq,gamma)*expl(-dt/tau)*dt/tau/tau;
	}
      else
	{
	  return 0.0;
	}
      
    }
  if (label == param_expindex)
    {
      if (dt  > 0)
	{
	  return val = dm*log(freq)*powl(freq,gamma)*expl(-dt/tau);
	}
      else{
	return 0.0;
      }
    }
  if (label == param_expep)
    {
      if (dt  > 0)
	{
	  return val = dm*powl(freq,gamma)*expl(-dt/tau)/tau;
	}
      else
	{
	  return 0.0;
	}
    }
      
    
      

	
	    

  



  return 0;
}




double t2FitFunc_gausdip(pulsar *psr, int ipsr, double x, int ipos, param_label label, int k)
{
  
  
  long double val;

  long double dt;
  long double amp;
  long double sig;
  long double freq;
  long double gamma;
  

  
  // reference to 1.4 GHz  to agree with Enterprise
  freq= psr[ipsr].obsn[ipos].freqSSB/1.4e9;
  
  dt=(psr[ipsr].obsn[ipos].bbat - psr[ipsr].param[param_gausep].val[k]);
  amp=psr[ipsr].param[param_gausamp].val[k];
  sig=psr[ipsr].param[param_gaussig].val[k];


 // if index is not set assume the Gaussian is achromatic
  if (psr[ipsr].param[param_gausindex].paramSet[k] ==1)
    {
      gamma=psr[ipsr].param[param_gausindex].val[k];
    }
  else{
    gamma=0;
    
  }
  
// model  is
// amp*powl(freq, gamma)*exp(-dt*dt/2./sig/sig);

    val= amp*powl(freq, gamma)*exp(-dt*dt/2./sig/sig);     

  if (label ==param_gausamp)
    {
        return val/amp;
    }
  if( label ==param_gausep)
    {
        return val*(-dt/sig/sig); 
    }
    if (label==param_gaussig)
    {
        return val*(dt*dt/6./sig/sig/sig);    
    }
    if (label==param_gausindex)
    {
        return val*log(freq);
    }
  return 0;
}

