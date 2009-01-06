#include "tempo2.h"
//#include "tempo2Util.h"
//#include "tempo2pred.h"
//#include "tempo2pred_int.h"

char TEMPO2_ENVIRON[MAX_STRLEN]="TEMPO2";
char TEMPO2_ERROR[MAX_STRLEN]="";

int MAX_PSR   = MAX_PSR_VAL;    /* Maximum number of pulsars to fit simultaneously  */
int MAX_OBSN  = MAX_OBSN_VAL;
double ECLIPTIC_OBLIQUITY = ECLIPTIC_OBLIQUITY_VAL;
int debugFlag = 0;
int veryFast = 0;


void extra_delays(pulsar *psr,int npsr)
{  
  calculate_bclt(psr,npsr);/* 3. Calculate bclt  */
  /*  shapiro_delay(psr,npsr); */ /* 1. Calculate the Shapiro delay */
  /* dm_delays(psr,npsr); */    /* 2. Extra dispersion measure delays */  
}

void clock_corrections(pulsar *psr,int npsr)
{  
  if (debugFlag==1) printf("Calling toa2utc\n");
  toa2utc(psr,npsr);        /* 1. UTC(Observatory) -> UTC(NIST) */
  if (debugFlag==1) printf("Calling tai2ut1\n");
//   utc2tai(psr,npsr);     /* 2. UTC(NIST) -> TAI              */
  tai2ut1(psr,npsr);        /* 3. TAI -> UT1                    */
//   tai2tt(psr,npsr);      /* 4. TAI -> TT                     */
  if (debugFlag==1) printf("Calling tt2tb\n");
  tt2tb(psr,npsr);          /* 5. Rough estimate of TT-TB (+-2.2 microsec) */
  if (debugFlag==1) printf("Done clock corrections\n");
}

void ephemeris_routines(pulsar *psr,int npsr)
{ 
  vectorPulsar(psr,npsr);   /* 1. Form a vector pointing at the pulsar */
  readEphemeris(psr,npsr,0);/* 2. Read the ephemeris */
  get_obsCoord(psr,npsr);   /* 3. Get Coordinate of observatory relative to Earth's centre */
  tt2tb(psr,npsr);          /* Observatory/time-dependent part of TT-TB */
  readEphemeris(psr,npsr,0);  /* Re-evaluate ephemeris with correct TB */ 
}

void formBatsAll(pulsar *psr,int npsr)
{
  if (debugFlag==1) printf("Calling clock corrections\n");
  clock_corrections(psr,npsr);          /* Clock corrections  ... */  
  if (debugFlag==1) printf("Reading ephemeris routines\n");
  ephemeris_routines(psr,npsr);         /* Ephemeris routines ... */
  if (debugFlag==1) printf("Reading extra delays\n");
  extra_delays(psr,npsr);               /* Other time delays  ... */
  formBats(psr,npsr);                   /* Form Barycentric arrival times */
  secularMotion(psr,npsr); 
}
