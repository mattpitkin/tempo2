#include <assert.h>
#include "t2fit_gw.h"


double t2FitFunc_gwm_amp(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k){
    assert(label==param_gwm_amp);
    longdouble dt;
    double res;

    if (psr[ipsr].param[param_gwm_amp].paramSet[1]==1){
        dt = (psr[ipsr].obsn[ipos].bbat - psr[ipsr].gwm_epoch)*longdouble(86400.0);
        if (dt > 0)
        {
            if (k==0)
                res = psr[ipsr].quad_ifunc_geom_p*dt;
            else if (k==1)
                res = psr[ipsr].quad_ifunc_geom_c*dt;
        }
        else res = 0;
    }
    else
    {
        double n1,n2,n3;
        double cosTheta;
        double lambda_p,beta_p,lambda,beta;
        double g1,g2,g3;


        if (psr[ipsr].param[param_raj].paramSet[1] == 1)
            lambda_p = (double)psr[ipsr].param[param_raj].val[1];
        else
            lambda_p = (double)psr[ipsr].param[param_raj].val[0];

        if (psr[ipsr].param[param_decj].paramSet[1] == 1)
            beta_p   = (double)psr[ipsr].param[param_decj].val[1];
        else
            beta_p   = (double)psr[ipsr].param[param_decj].val[0];

        lambda   = psr[ipsr].gwm_raj;
        beta     = psr[ipsr].gwm_decj;

        // GW vector
        g1 = -cosl(lambda)*cosl(beta);
        g2 = -sinl(lambda)*cosl(beta);
        g3 = -sinl(beta);

        // Pulsar vector
        n1 = cosl(lambda_p)*cosl(beta_p);
        n2 = sinl(lambda_p)*cosl(beta_p);
        n3 = sinl(beta_p);
        cosTheta = -(cosl(beta)*cosl(beta_p)*cosl(lambda-lambda_p)+
                sinl(beta)*sinl(beta_p));

        /* Only has effect after the glitch epoch */
        if (psr[ipsr].obsn[ipos].sat >= psr[ipsr].gwm_epoch)
        {
            longdouble dt,scale;
            double cos2Phi;
            double cosPhi;
            double l1,l2,l3,m1,m2,m3;
            double d1,d2,d3,md;
            double a1,a2,a3,ma;

            //   if  (g3 != 0) 
            //	   {beta_m = atan2(-cos(beta)*cos(lambda-psr[ipsr].gwm_phi),sin(beta));}
            //  else  
            //      {beta_m = atan2(sinl(psr[ipsr].gwm_phi),cosl(psr[ipsr].gwm_phi));
            //       psr[ipsr].gwm_phi = lambda + 1.5708;}
            //  m1 = cosl(psr[ipsr].gwm_phi)*cosl(beta_m);
            //  m2 = sinl(psr[ipsr].gwm_phi)*cosl(beta_m);
            //  m3 = sinl(beta_m);

            if (beta == 0.0 )
            {
                d1 = 0.0;
                d2 = 0.0;
                d3 = 1.0;
            }

            if ( beta > 0)
            {
                d1 = g1*cosl(0.5*M_PI - beta);
                d2 = g2*cosl(0.5*M_PI - beta);
                d3 = 1.0 + g3*cos(0.5*M_PI - beta);
                md = sqrt(d1*d1 + d2*d2 + d3*d3);
                d1 = d1/md;
                d2 = d2/md;
                d3 = d3/md;
                /*covert d to unit vector */
            } 
            else if (beta < 0) 
            {
                d1 = g1*cosl(-0.5*M_PI - beta);
                d2 = g2*cosl(-0.5*M_PI - beta);
                d3 = -1.0 + g3*cos(-0.5*M_PI - beta);
                md = sqrt(d1*d1 + d2*d2 + d3*d3);
                d1 = d1/md;
                d2 = d2/md;
                d3 = d3/md;
            } 

            //if (g2*d3-d2*g3 != 0)
            // {
            //  a1 = 1.0; 
            //  a2 = (d1*g3-g1*d3)/(g2*d3-d2*g3);
            //  a3 = (g2*d1-g1*d2)/(g3*d2-g2*d3); 
            // }
            //else if (g1*d3-d1*g3 != 0)
            // {
            //  a1 = (g3*d2-d3*g2)/(g1*d3-g3*d1); 
            //  a2 = 1.0;
            //  a3 = (g1*d2-d1*g2)/(g3*d1-d1*d3);
            // }
            //else if (d2*g1-g2*d1 != 0)			
            // {
            //  a1 = (g2*d3-d2*g3)/(d2*g1-g2*d1); 
            //  a2 = (g1*d3-d1*g3)/(d1*g2-g1*d2);
            //  a3 =1.0; 
            // }
            a1 =  (d2*g3-d3*g2);
            a2 =  (d3*g1-d1*g3);
            a3 =  (d1*g2-d2*g1);

            /* conver it to unit vector */
            ma = sqrt(a1*a1 +a2*a2 + a3*a3);
            a1 = a1/ma;
            a2 = a2/ma;
            a3 = a3/ma;

            /* polarisation vector of GW source */
            m1 = d1*cosl(psr[ipsr].gwm_phi)	+ a1*sinl(psr[ipsr].gwm_phi);   
            m2 = d2*cosl(psr[ipsr].gwm_phi)	+ a2*sinl(psr[ipsr].gwm_phi);
            m3 = d3*cosl(psr[ipsr].gwm_phi)	+ a3*sinl(psr[ipsr].gwm_phi);

            if  (cosTheta != 1.0 && cosTheta != -1.0)
            {g1 = g1*cosTheta; 
                g2 = g2*cosTheta;
                g3 = g3*cosTheta;
                l1 = n1 - g1;
                l2 = n2 - g2;
                l3 = n3 - g3;
                cosPhi = (l1*m1 + l2*m2 + l3*m3)/sqrt(l1*l1 + l2*l2 + l3*l3);
                //		 if  (cosPhi >= 1.0/sqrt(2.0))
                cos2Phi = 2*cosPhi*cosPhi - 1.0;
                //		 else
                //		     cos2Phi = 2*sqrt(1.0 - cosPhi*cosPhi)*sqrt(1.0 - cosPhi*cosPhi) - 1.0;
            }
            else 
            {cos2Phi = 0;}

            dt = (psr[ipsr].obsn[ipos].sat - psr[ipsr].gwm_epoch)*86400.0;
            scale = -0.5*cos2Phi*(1-cosTheta);
            //	    scale=1;
            res = scale*dt;
        }
        else
            res = 0;
    }
    return res;
}




double t2FitFunc_gwb_amp(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k){
    double res=0;
    longdouble dt;
    longdouble prefac;
    if (psr[ipsr].param[param_gwb_amp].paramSet[1]==1)
    {
        dt = (psr[ipsr].obsn[ipos].bbat - psr[ipsr].gwb_epoch)/psr[ipsr].gwb_width;
        prefac = dt*exp( (double) -dt*dt/2.);

        if (k==0)
        {
            res = psr[ipsr].gwb_geom_p*prefac;
        }
        else if (k==1)
        {
            res = psr[ipsr].gwb_geom_c*prefac;
        }
        else
        {
            res = 0;
        }
    }
    return res;

}


double t2FitFunc_gwcs_amp(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k){
    double res=0;
    longdouble dt;
    longdouble prefac;
    double extra;
    longdouble width,width_day;
    
    if (psr[ipsr].param[param_gwcs_amp].paramSet[0]==1)
    {
        dt = (psr[ipsr].obsn[ipos].bbat - psr[ipsr].gwcs_epoch)*86400.0;
	width = psr[ipsr].gwcs_width*86400.0;
	width_day = psr[ipsr].gwcs_width;
	
	
	if (psr[ipsr].obsn[ipos].sat < psr[ipsr].gwcs_epoch-width_day/2.0)
	  res=0;
	else if (psr[ipsr].obsn[ipos].sat < psr[ipsr].gwcs_epoch)
	  {
	    extra =   (3.0/4.0*(pow(0.5*width,4.0/3.0)-pow(fabs(dt),4.0/3.0))-
		       pow(0.5*width,1.0/3.0)*(dt+0.5*width));
	    
	  }
	else if (psr[ipsr].obsn[ipos].sat < psr[ipsr].gwcs_epoch+width_day/2.0)
	  {
	    extra =   (3.0/4.0*(pow(0.5*width,4.0/3.0)+pow(fabs(dt),4.0/3.0))-
		       pow(0.5*width,1.0/3.0)*(dt+0.5*width));
	    
	  }
	else
	  {
	    extra=-0.25*(pow(0.5,1.0/3.0)*pow(width,4.0/3.0));
	  }
    }

    if (k==0)
      res = psr[ipsr].gwcs_geom_p*extra;
    else if (k==1)
      res = psr[ipsr].gwcs_geom_c*extra;

    return res;

}




double t2FitFunc_quad_om(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k) {
    assert(label==param_quad_om);
    double afunc = 0;
    double n1,n2,n3;
    double e11p,e21p,e31p,e12p,e22p,e32p,e13p,e23p,e33p;
    double e11c,e21c,e31c,e12c,e22c,e32c,e13c,e23c,e33c;
    double cosTheta,omega_g;
    longdouble resp,resc;
    double lambda_p,beta_p,lambda,beta;
    longdouble time;
    time    = (psr[ipsr].obsn[ipos].bbat - psr[ipsr].quadEpoch)*longdouble(86400.0);
    if (psr[ipsr].param[param_raj].paramSet[1] == 1)
        lambda_p = (double)psr[ipsr].param[param_raj].val[1];
    else
        lambda_p = (double)psr[ipsr].param[param_raj].val[0];

    if (psr[ipsr].param[param_raj].paramSet[1] == 1)
        beta_p   = (double)psr[ipsr].param[param_decj].val[1];
    else
        beta_p   = (double)psr[ipsr].param[param_decj].val[0];

    lambda   = psr[ipsr].quadRA;
    beta     = psr[ipsr].quadDEC;
    // Pulsar vector
    n1 = cosl(lambda_p)*cosl(beta_p);
    n2 = sinl(lambda_p)*cosl(beta_p);
    n3 = sinl(beta_p);
    cosTheta = cosl(beta)*cosl(beta_p)*cosl(lambda-lambda_p)+
        sinl(beta)*sinl(beta_p);

    e11p = pow(sinl(lambda),2)-pow(cosl(lambda),2)*pow(sinl(beta),2);
    e21p = -sinl(lambda)*cosl(lambda)*(pow(sinl(beta),2)+1);
    e31p = cosl(lambda)*sinl(beta)*cosl(beta);

    e12p = -sinl(lambda)*cosl(lambda)*(pow(sinl(beta),2)+1);
    e22p = pow(cosl(lambda),2)-pow(sinl(lambda),2)*pow(sinl(beta),2);
    e32p = sinl(lambda)*sinl(beta)*cosl(beta);

    e13p = cosl(lambda)*sinl(beta)*cosl(beta);
    e23p = sinl(lambda)*sinl(beta)*cosl(beta);
    e33p = -powl(cosl(beta),2);

    resp = (n1*(n1*e11p+n2*e12p+n3*e13p)+
            n2*(n1*e21p+n2*e22p+n3*e23p)+
            n3*(n1*e31p+n2*e32p+n3*e33p));

    e11c = sin(2*lambda)*sin(beta);
    e21c = -cos(2*lambda)*sin(beta);
    e31c = -sin(lambda)*cos(beta);

    e12c = -cos(2*lambda)*sin(beta);
    e22c = -sin(2*lambda)*sin(beta);
    e32c = cos(lambda)*cos(beta);

    e13c = -sin(lambda)*cos(beta);
    e23c = cos(lambda)*cos(beta);
    e33c  = 0;

    resc = (n1*(n1*e11c+n2*e12c+n3*e13c)+
            n2*(n1*e21c+n2*e22c+n3*e23c)+
            n3*(n1*e31c+n2*e32c+n3*e33c));

    omega_g = (double)psr[ipsr].param[param_quad_om].val[0]*((int)(k/4.0)+1);

    //            printf("In fitting with k = %d %d\n",k,k%4);
    if ((longdouble(1.0)-cosTheta)==0)
        afunc=0;
    else
    {
        if (k%4==0)      afunc = resp*sin(omega_g*time)/(longdouble(2.0)*omega_g*(longdouble(1.0)-cosTheta)); // aplus_re
        else if (k%4==1) afunc = resp*cos(omega_g*time)/(longdouble(2.0)*omega_g*(longdouble(1.0)-cosTheta)); // aplus_im
        else if (k%4==2) afunc = resc*sin(omega_g*time)/(longdouble(2.0)*omega_g*(longdouble(1.0)-cosTheta)); // across_re
        else if (k%4==3) afunc = resc*cos(omega_g*time)/(longdouble(2.0)*omega_g*(longdouble(1.0)-cosTheta)); // across_im

        if (psr[ipsr].gwsrc_psrdist > 0) // Add in the pulsar term
        {
            if (k%4==0) afunc      -=  resp*sinl(omega_g*time-(1-cosTheta)*psr[ipsr].gwsrc_psrdist/SPEED_LIGHT*omega_g)/(longdouble(2.0)*omega_g*(longdouble(1.0)-cosTheta)); // aplus_im
            else if (k%4==1) afunc -=  resp*cosl(omega_g*time-(1-cosTheta)*psr[ipsr].gwsrc_psrdist/SPEED_LIGHT*omega_g)/(longdouble(2.0)*omega_g*(longdouble(1.0)-cosTheta)); // aplus_im
            else if (k%4==2) afunc -=  resc*sinl(omega_g*time-(1-cosTheta)*psr[ipsr].gwsrc_psrdist/SPEED_LIGHT*omega_g)/(longdouble(2.0)*omega_g*(longdouble(1.0)-cosTheta)); // aplus_im
            else if (k%4==3) afunc -=  resc*cosl(omega_g*time-(1-cosTheta)*psr[ipsr].gwsrc_psrdist/SPEED_LIGHT*omega_g)/(longdouble(2.0)*omega_g*(longdouble(1.0)-cosTheta)); // aplus_im
        }


    }
    return afunc;
}


void t2UpdateFunc_quad_om(pulsar *psr, int ipsr ,param_label label,int k, double val, double err) {

    int quad_j = k / 4;
    int quad_k = k % 4;
    switch(quad_k){
        case 0:
            psr[ipsr].quad_aplus_r[quad_j] -= val;
            psr[ipsr].quad_aplus_r_e[quad_j] = err;
            break;
        case 1:
            psr[ipsr].quad_aplus_i[quad_j] -= val;
            psr[ipsr].quad_aplus_i_e[quad_j] = err;
            break;
        case 2:
            psr[ipsr].quad_across_r[quad_j] -= val;
            psr[ipsr].quad_across_r_e[quad_j]  = err;
            break;
        case 3:
            psr[ipsr].quad_across_i[quad_j] -= val;
            psr[ipsr].quad_across_i_e[quad_j] = err;
            break;
    }
}


double t2FitFunc_gwsingle(pulsar *psr, int ipsr ,double x ,int ipos ,param_label label,int k) {
    assert(label==param_gwsingle);
    double afunc = 0;
    double n1,n2,n3;
    double e11p,e21p,e31p,e12p,e22p,e32p,e13p,e23p,e33p;
    double e11c,e21c,e31c,e12c,e22c,e32c,e13c,e23c,e33c;
    double cosTheta,omega_g;
    longdouble resp,resc;
    double lambda_p,beta_p,lambda,beta;
    longdouble time;

    time    = (psr[ipsr].obsn[ipos].bbat - psr[ipsr].gwsrc_epoch)*longdouble(86400.0);
    lambda_p = (double)psr[ipsr].param[param_raj].val[0];
    beta_p   = (double)psr[ipsr].param[param_decj].val[0];
    lambda   = psr[ipsr].gwsrc_ra;
    beta     = psr[ipsr].gwsrc_dec;
    // Pulsar vector
    n1 = cosl(lambda_p)*cosl(beta_p);
    n2 = sinl(lambda_p)*cosl(beta_p);
    n3 = sinl(beta_p);
    cosTheta = cosl(beta)*cosl(beta_p)*cosl(lambda-lambda_p)+
        sinl(beta)*sinl(beta_p);

    e11p = pow(sinl(lambda),2)-pow(cosl(lambda),2)*pow(sinl(beta),2);
    e21p = -sinl(lambda)*cosl(lambda)*(pow(sinl(beta),2)+1);
    e31p = cosl(lambda)*sinl(beta)*cosl(beta);

    e12p = -sinl(lambda)*cosl(lambda)*(pow(sinl(beta),2)+1);
    e22p = pow(cosl(lambda),2)-pow(sinl(lambda),2)*pow(sinl(beta),2);
    e32p = sinl(lambda)*sinl(beta)*cosl(beta);

    e13p = cosl(lambda)*sinl(beta)*cosl(beta);
    e23p = sinl(lambda)*sinl(beta)*cosl(beta);
    e33p = -powl(cosl(beta),2);

    omega_g = (double)psr[ipsr].param[param_gwsingle].val[0];

    resp = (n1*(n1*e11p+n2*e12p+n3*e13p)+
            n2*(n1*e21p+n2*e22p+n3*e23p)+
            n3*(n1*e31p+n2*e32p+n3*e33p));

    e11c = sin(2*lambda)*sin(beta);
    e21c = -cos(2*lambda)*sin(beta);
    e31c = -sin(lambda)*cos(beta);

    e12c = -cos(2*lambda)*sin(beta);
    e22c = -sin(2*lambda)*sin(beta);
    e32c = cos(lambda)*cos(beta);

    e13c = -sin(lambda)*cos(beta);
    e23c = cos(lambda)*cos(beta);
    e33c  = 0;

    resc = (n1*(n1*e11c+n2*e12c+n3*e13c)+
            n2*(n1*e21c+n2*e22c+n3*e23c)+
            n3*(n1*e31c+n2*e32c+n3*e33c));

    if ((longdouble(1.0)-cosTheta)==0)
        afunc=0;
    else
    {
        if (k==0)      afunc = resp*sin(omega_g*time)/(longdouble(2.0)*omega_g*(longdouble(1.0)-cosTheta)); // aplus_re
        else if (k==1) afunc = resc*sin(omega_g*time)/(longdouble(2.0)*omega_g*(longdouble(1.0)-cosTheta)); // across_re
        else if (k==2) afunc = resp*(cos(omega_g*time))/(longdouble(2.0)*omega_g*(longdouble(1.0)-cosTheta)); // aplus_im
        else if (k==3) afunc = resc*(cos(omega_g*time))/(longdouble(2.0)*omega_g*(longdouble(1.0)-cosTheta)); // across_im

    }
    return afunc;
}

void t2UpdateFunc_gwsingle(pulsar *psr, int ipsr ,param_label label,int k, double val, double err) {
    switch(k){
        case 0:
            psr[ipsr].gwsrc_aplus_r -= val;
            break;
        case 1:
            psr[ipsr].gwsrc_across_r -= val;
            break;
        case 2:
            psr[ipsr].gwsrc_aplus_i -= val;
            break;
        case 3:
            psr[ipsr].gwsrc_across_i -= val;
            break;
    }
}
