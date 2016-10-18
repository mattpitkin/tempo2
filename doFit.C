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


#include "tempo2.h"
#include "t2fit.h"
#include <string.h>
#include <dlfcn.h>
#include <stdlib.h>


/**
 * Master fitting routine with or without cholesky, global or not.
 */
void doFitAll(pulsar *psr,int npsr, const char *covarFuncFile) {
    t2Fit(psr,npsr,covarFuncFile);
}

double getParamDeriv(pulsar *psr,int ipos,double x,int i,int k){
    return t2Fit_getParamDeriv(psr,ipos,x,i,k);
}



void callFitFuncPlugin(pulsar *psr,int npsr, const char *covarFuncFile){
    if (strcmp(psr[0].fitFunc,"default")!=0)
    {
        char *(*entry)(pulsar *,int,const char*);
        void * module;
        char str[100];

        logmsg("Calling fitting plugin: %s",psr[0].fitFunc);

        logmsg("%d",tempo2_plug_path_len);
        for (int iplug=0; iplug < tempo2_plug_path_len; iplug++) {
            sprintf(str,"%s/%s_fitFunc_%s_plug.t2",tempo2_plug_path[iplug],
                    psr[0].fitFunc,tempo2MachineType);
            logmsg("Looking for '%s'",str);
            module = dlopen(str, RTLD_NOW); 
            if (module==NULL) {
                logerr("dlerror() = %s",dlerror());
            } else break;
        }

        module = dlopen(str, RTLD_NOW); 
        if(!module)  {
            logerr("dlopen() unable to open plugin: %s.",str);
            exit(1);
        }
        /*
         * Check that the plugin is compiled against the same version of tempo2.h
         */
        char ** pv  = (char**)dlsym(module, "plugVersionCheck");
        if(pv!=NULL){
            // there is a version check for this plugin
            if(strcmp(TEMPO2_h_VER,*pv)){
                logerr("Plugin version mismatch");
                logerr(" '%s' != '%s'",TEMPO2_h_VER,*pv);
                logerr(" Please recompile plugin against same tempo2 version!");
                dlclose(module);
                exit(1);
            }
        }


        entry = (char*(*)(pulsar *,int,const char*))dlsym(module, "pluginFitFunc");
        if( entry == NULL ) {
            dlclose(module);
            logerr("dlerror() failed while retrieving address.");
            exit(1);
        }
        entry(psr,npsr,covarFuncFile);
        logmsg("FitFunc plugin Returning\n");

        return;
    }
}
