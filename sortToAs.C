#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <tempo2.h>

int compareObs(const void *o1, const void *o2);
/**
 * Sort ToAs for one pulsar.
 */
void sortToAs(pulsar* psr){
   if (!psr->sorted){
	  logmsg("Sorting ToAs for psr %s",psr->name);
	  logtchk("call qsort()");
	  qsort(psr->obsn,psr->nobs,sizeof(observation),compareObs);
	  psr->sorted=1;
	  logtchk("exit qsort()");
   }

}

int compareObs(const void *o1, const void *o2)
{
   const observation *elem1 = (const observation*) o1;    
   const observation *elem2 = (const observation*) o2;

   if ( elem1->sat < elem2->sat)
	  return -1;
   else if (elem1->sat > elem2->sat)
	  return 1;
   else
	  return 0;
}
