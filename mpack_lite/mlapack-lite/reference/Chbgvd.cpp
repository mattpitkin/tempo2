/*
 * Copyright (c) 2008-2010
 *      Nakata, Maho
 *      All rights reserved.
 *
 *  $Id: Chbgvd.cpp,v 1.3 2010/08/07 04:48:32 nakatamaho Exp $ 
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 */
/*
Copyright (c) 1992-2007 The University of Tennessee.  All rights reserved.

$COPYRIGHT$

Additional copyrights may follow

$HEADER$

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

- Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer. 
  
- Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer listed
  in this license in the documentation and/or other materials
  provided with the distribution.
  
- Neither the name of the copyright holders nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.
  
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT  
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT  
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
*/

#include <mblas.h>
#include <mlapack.h>

void Chbgvd(const char *jobz, const char *uplo, INTEGER n, INTEGER ka,
	    INTEGER kb, COMPLEX * AB, INTEGER ldab, COMPLEX * BB,
	    INTEGER ldbb, REAL * w, COMPLEX * z, INTEGER ldz, COMPLEX * work, INTEGER lwork, REAL * rwork, INTEGER lrwork, INTEGER * iwork, INTEGER liwork, INTEGER * info)
{
    INTEGER inde;
    char vect;
    INTEGER llwk2;
    INTEGER iinfo;
    INTEGER lwmin;
    INTEGER upper;
    INTEGER llrwk;
    INTEGER wantz;
    INTEGER indwk2;
    INTEGER lrwmin;
    INTEGER lquery;
    INTEGER liwmin;
    INTEGER indwrk;
    REAL Zero = 0.0, One = 1.0;

//Test the input parameters.
    wantz = Mlsame(jobz, "V");
    upper = Mlsame(uplo, "U");
    lquery = lwork == -1 || lrwork == -1 || liwork == -1;

    *info = 0;
    if (n <= 1) {
	lwmin = 1;
	lrwmin = 1;
	liwmin = 1;
    } else if (wantz) {
	lwmin = n * n * 2;
	lrwmin = n * 5 + 1 + n * n * 2;
	liwmin = n * 5 + 3;
    } else {
	lwmin = n;
	lrwmin = n;
	liwmin = 1;
    }
    if (!(wantz || Mlsame(jobz, "N"))) {
	*info = -1;
    } else if (!(upper || Mlsame(uplo, "L"))) {
	*info = -2;
    } else if (n < 0) {
	*info = -3;
    } else if (ka < 0) {
	*info = -4;
    } else if (kb < 0 || kb > ka) {
	*info = -5;
    } else if (ldab < ka + 1) {
	*info = -7;
    } else if (ldbb < kb + 1) {
	*info = -9;
    } else if (ldz < 1 || (wantz && ldz < n)) {
	*info = -12;
    }
    if (*info == 0) {
	work[1] = lwmin;
	rwork[1] = lrwmin;
	iwork[1] = liwmin;
	if (lwork < lwmin && !lquery) {
	    *info = -14;
	} else if (lrwork < (lrwmin && !lquery)) {
	    *info = -16;
	} else if (liwork < (liwmin && !lquery)) {
	    *info = -18;
	}
    }
    if (*info != 0) {
	Mxerbla("Chbgvd", -(*info));
	return;
    } else if (lquery) {
	return;
    }
//Quick return if possible
    if (n == 0) {
	return;
    }
//Form a split Cholesky factorization of B.
    Cpbstf(uplo, n, kb, &BB[0], ldbb, info);
    if (*info != 0) {
	*info = n + *info;
	return;
    }
//Transform problem to standard eigenvalue problem.
    inde = 1;
    indwrk = inde + n;
    indwk2 = n * n + 1;
    llwk2 = lwork - indwk2 + 2;
    llrwk = lrwork - indwrk + 2;
    Chbgst(jobz, uplo, n, ka, kb, &AB[0], ldab, &BB[0], ldbb, &z[0], ldz, &work[0], &rwork[indwrk], &iinfo);
//Reduce Hermitian band matrix to tridiagonal form.
    if (wantz) {
	vect = 'U';
    } else {
	vect = 'N';
    }
    Chbtrd(&vect, uplo, n, ka, &AB[0], ldab, &w[1], &rwork[inde], &z[0], ldz, &work[0], &iinfo);
//For eigenvalues only, call DSTERF.  For eigenvectors, call ZSTEDC.
    if (!wantz) {
	Rsterf(n, &w[1], &rwork[inde], info);
    } else {
	Cstedc("I", n, &w[1], &rwork[inde], &work[0], n, &work[indwk2], llwk2, &rwork[indwrk], llrwk, &iwork[1], liwork, info);
	Cgemm("N", "N", n, n, n, (COMPLEX) One, &z[0], ldz, &work[0], n, (COMPLEX) Zero, &work[indwk2], n);
	Clacpy("A", n, n, &work[indwk2], n, &z[0], ldz);
    }
    work[1] = lwmin;
    rwork[1] = lrwmin;
    iwork[1] = liwmin;
    return;
}
