// redwards, chebyshev-related code NOT for external linkage in libtempo2pred

#include "tempo2pred.h"
#include "tempo2pred_int.h"
#include "tempo2.h"
#include <math.h>
#include <assert.h>
 
// struct to hold info needed to evaluate timing model during chebyshev
// construction
typedef struct
{
  ChebyModel *model;
  const pulsar *psr;
  bool compute_dispersion_constant;
} ChebyModelInfo;

void
chebyModelFunc(long double *x, long double *y, 
	       int nx, int ny, long double *z, void *info_in)
{
  ChebyModelInfo *info = (ChebyModelInfo *)info_in;
  pulsar *psr = const_cast<pulsar*>( info->psr );
  int ix, iy, iobs=1;

  if (psr->param[param_tzrmjd].paramSet[0] == 0) {
    fprintf (stderr, "chebyModelFunc: TZRMJD not set in pulsar parameters\n");
    exit(-1);
  }
  if (psr->param[param_tzrfrq].paramSet[0] == 0) {
    fprintf (stderr, "chebyModelFunc: TZRFRQ not set in pulsar parameters\n");
    exit(-1);
  }

  // set up a bunch of observations in the psr struct
  psr->nobs = nx*ny+1;
  psr->param[param_track].paramSet[0] = 0;
  for (iy=0; iy < ny; iy++)
    for (ix=0; ix < nx; ix++)
    {
      psr->obsn[iobs].sat = info->model->mjd_start + 
	(x[ix]+1.0L)*0.5L*(info->model->mjd_end-info->model->mjd_start);
      psr->obsn[iobs].freq =  info->model->freq_start + 
       (y[iy]+1.0L)*0.5L*(info->model->freq_end-info->model->freq_start);
      psr->obsn[iobs].deleted = 0;
      psr->obsn[iobs].nFlags = 0;
      psr->obsn[iobs].delayCorr = 1;
      psr->obsn[iobs].clockCorr = 1;
      psr->obsn[iobs].phaseOffset = 0.0;
      strcpy(psr->obsn[iobs].telID, info->model->sitename);
      if (strcmp(psr->obsn[iobs].telID,"@")==0 || strcasecmp(psr->obsn[iobs].telID,"bat")==0)
	{
	  psr->obsn[iobs].clockCorr=0;  /* therefore don't do clock corrections */
	  psr->obsn[iobs].delayCorr=0;
	}

      iobs++;
    }
  // stick in a fake obs to get the reference phase
  psr->obsn[0].sat = psr->param[param_tzrmjd].val[0];
  psr->obsn[0].freq = psr->param[param_tzrfrq].val[0];

  strcpy(psr->obsn[0].telID,  info->model->sitename);//XXX!!!
  psr->obsn[0].deleted = 0;
  psr->obsn[0].nFlags = 0;
  psr->obsn[0].delayCorr = 1;
  psr->obsn[0].clockCorr = 1;
  psr->obsn[0].phaseOffset = 0;

  if (strcmp(psr->obsn[0].telID,"@")==0 || strcasecmp(psr->obsn[0].telID,"bat")==0)
    {
      psr->obsn[0].clockCorr=0;  /* therefore don't do clock corrections */
      psr->obsn[0].delayCorr=0;
    }

  // Compute the phase for each observation
  formBatsAll(psr, 1);
  formResiduals(psr, 1, 0.0);
  // Get a nominal apparent dispersion constant
  if (info->compute_dispersion_constant)
  {
    int iref0 = nx/2, iref1 = nx/2 + (ny-1)*nx;
    info->model->dispersion_constant = 
      (psr->obsn[iref0].phase - psr->obsn[iref1].phase)
      / (1.0/(psr->obsn[iref0].freq*psr->obsn[iref0].freq)
	 - 1.0/(psr->obsn[iref1].freq*psr->obsn[iref1].freq));
  }
  //  info->model->dispersion_constant = 0.0;

  // fill in the "z" array with the phases
  long double rphase = psr->obsn[0].phase;
  for (iobs=1; iobs < psr->nobs; iobs++)
  {
    z[iobs-1] = psr->obsn[iobs].phase-rphase
      - info->model->dispersion_constant 
      / (psr->obsn[iobs].freq*psr->obsn[iobs].freq);
  }
//   printf("%Lg %Lg %Lg\n", z[0]+rphase, rphase, z[0]);
//   exit(1);
    // this is what we should use, only it gets set to something wierd
    //    z[iobs] = psr->obsn[iobs].phase;
}


void
ChebyModel_Construct(ChebyModel *cm, const pulsar *psr)
{
  // construct the phase chebyshev polynomial
  ChebyModelInfo info;
  strcpy(cm->psrname, psr->name);
  info.model = cm;
  info.psr = psr;
  info.compute_dispersion_constant = true;

  Cheby2D_Construct(&cm->cheby, chebyModelFunc, (void*)&info);
  // compute the derivatives
  Cheby2D_Construct_x_Derivative(&cm->frequency_cheby, &cm->cheby);
}


void
ChebyModel_Test(ChebyModel *cm, const pulsar *psr, int nmjd, int nfreq,
		long double *residualRMS, long double *residualMAV)
{
  ChebyModelInfo info;
  info.model = cm;
  info.psr = psr;
  info.compute_dispersion_constant = false; // use the one already there!
  Cheby2D_Test(&cm->cheby, nmjd, nfreq, chebyModelFunc, (void *)&info,
	       residualRMS, residualMAV);
}


void ChebyModelSet_Construct(ChebyModelSet *cms, const pulsar *psr,
			     const char *sitename,
			     long double mjd_start, long double mjd_end,
			     long double segment_length, long double overlap,
			     long double freq_start,long double freq_end,
			     int ntimecoeff, int nfreqcoeff)
{
  cms->nsegments = (int)ceil((mjd_end-mjd_start)/(segment_length-overlap));

  if (debugFlag)
    fprintf (stderr, "ChebyModelSet_Construct MJD start=%f end=%f \n"
                     "\t length=%f overlap=%f nsegments=%d\n", 
                     (double)mjd_start, (double)mjd_end,
                     (double)segment_length, (double)overlap, cms->nsegments);

  assert (cms->nsegments > 0);

  if (psr->param[param_tzrmjd].paramSet[0]==0)
  {
    if (debugFlag)
      printf("WARNING: tzrmjd not set.  Setting to %g\n",(double)mjd_start);
    psr->param[param_tzrmjd].paramSet[0]=1;
    psr->param[param_tzrmjd].val[0] = mjd_start;
  }
  if (psr->param[param_tzrfrq].paramSet[0]==0)
  {
    if (debugFlag)
      printf("WARNING: tzrfrq not set.  Setting to %g\n",(double)freq_start);
    psr->param[param_tzrfrq].paramSet[0]=1;
    psr->param[param_tzrfrq].val[0] = freq_start;
  }
  if (strlen(psr->tzrsite)<1)
  {
    if (debugFlag)
      printf("WARNING: tzrsite not set.  Setting to %s\n",sitename);
    strcpy(const_cast<char*>(psr->tzrsite),sitename);
  }
  
  cms->segments = (ChebyModel *)malloc(cms->nsegments*sizeof(ChebyModel));

  int iseg;

  for (iseg=0; iseg < cms->nsegments ; iseg++)
  {
    strcpy(cms->segments[iseg].sitename, sitename);
    cms->segments[iseg].mjd_start = mjd_start+iseg*(segment_length-overlap);
    cms->segments[iseg].mjd_end = cms->segments[iseg].mjd_start
      +segment_length;
    if (cms->segments[iseg].mjd_end > mjd_end)
      cms->segments[iseg].mjd_end = mjd_end;
    cms->segments[iseg].freq_start = freq_start;
    cms->segments[iseg].freq_end = freq_end;
    ChebyModel_Init(&cms->segments[iseg], ntimecoeff, nfreqcoeff);
    ChebyModel_Construct(&cms->segments[iseg], psr);
  }
}

void ChebyModelSet_Test(ChebyModelSet *cms, const pulsar *psr,
			int nmjd, int nfreq,
			long double *residualRMS, long double *residualMAV)
{
  *residualRMS = *residualMAV = 0.0;
  long double rms, mav;
  int iseg;
  for (iseg=0; iseg < cms->nsegments ; iseg++)
  {
    ChebyModel_Test(&cms->segments[iseg], psr, nmjd/cms->nsegments,
		    nfreq/cms->nsegments, &rms, &mav);
    *residualRMS += rms*rms;
    if (mav > *residualMAV)
      *residualMAV = mav;
  }

  *residualRMS = sqrt(*residualRMS/cms->nsegments);
}

