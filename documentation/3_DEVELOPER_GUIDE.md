Developer Guide                                                 {#devguide}
===============
Tempo2 Developer Guide                                                 {#dgtop}
======================
[TOC]

About this guide        {#about}
----------------

This guide has been developed to encourage development of tempo2, and to improve the consistency between developers. The majority of this guide has been written by \ref MJK, although all are welcome to contribute.


General code guidelines       {#generalguide}
-----------------------

Tempo2 is, for historical reasons, mostly written in C but compiled using a C++ compiler. However, be aware that a few parts of tempo2 use C++ clases or other C++ extensions. There is no particular C or C++ version in use, but for now assume that we are using C++98 with GNU extensions (i.e. --std=gnu++98)

@todo determine if we should migrate to C++ 11. It has lots of good features, but we need to check that all compilers support it.


### Core tempo2 code
As a general rule, we try to minimise the libraries needed to build the core of tempo2 (not plugins). This means you can't link against libfftw, libpgplot, etc. from the core code. Some linear algebra features from BLAS/LAPACK are made avaliable to the code code via the T2toolkit, and fallback routines have been generated to ensure that the code still works without BLAS/LAPACK. These routines are being expanded all the time.

### plugins
For plugins, the rules are much less strict. Currently we compile plugins with links to cfitsio, fftw and pgplot as part of the main plugin distribution.

### libt2toolkit
\ref MJK is attempting to introduce a little more rigour in the coding standards for the code that makes up libt2toolkit, but in general this is treated exactly the same as code temo2.


Development workflow  {#workflow}
--------------------

### Recommended workflow
The recommended workflow is as follows.

Step 1: create a new branch:

    git checkout -b myfeature

Step 2: Make and commit your changes to that branch

    git commit -a

Step 3: Build, test, run your code.

    make
    make check

Step 4: If the new features seem good, promote them to the "master" branch.

    # if the first time
    git push --set-upstream origin docs
    # otherwise
    git push origin

and go to https://bitbucket.org/mkeith/tempo2/pull-requests/new to make a new pull request. The code will be reviewed by the core developers to check that the changes do not break any important features. If the modification is accepted (almost always) then it will be merged.

### Alternative workflow
If you can't be bothered with branches, you can simply work directly on the "dev" branch:

    git checkout dev

And commit as you want.

    git commit -a && git push origin

The dev branch will be merged into master, after code review, as and when required.
The drawbacks of this method are that you have to deal with conflicts yourself.



Coding style                                             {#style}
------------
Tempo2 does not have a strict coding style. However, it is recommended to adopt the following practice, as illustrated by the snippet below:

~~~~{.cpp}
// copyright statement up here.
#ifdef HAVE_CONFIG_H
#include <config.h> // make sure to include config.h
#endif

#include <cstdint>   // standard libries are included first
#include <fftw.h>    // then external libraries
#include "TKlog.h"  // then internal libraries


// functions are prefarably camelCase with small first letter.
// strings should be declared as const char* (or std::string) as they are immutable.
void myFunction(int anInt, const char *str, double **matrix) {
    // indent is 4 spaces.

    // use stdint types where possible to avoid confusion on 32-bit vs 64-bit machines.
    // use unsigned types whre sutable
    // use const when a variable will not change
    const uint64_t myconst = 1024;


    // keywords have a space before parenthesis (e.g. if, for, while).
    if (anInt < 10) { // always use braces, even if one line!
        // use TKlog for logging debug messages and warnings.

        // debug for statements that are to be printed when debug flag is set
        logdbg("anInt = %d",anInt);

        // warnings when problem might be an issue but can continue
        logwarn("anInt should be less than 10"); // adds a message to the warning stack

        // messages always appear
        logmsg("Print to terminal")

        // errors for when the operation is likely to fail.
        logerr("aborting because anInt was too large (%d)",anInt);

        // prefer to return on error rather than exit
        return;
    }

    // best to declare variables in for loops, but give them a proper name (not i, j, k) if possible.
    for (size_t iVal = 0; ival < myconst; ival++) {
        // ...
    }

}
~~~~

\note Core tempo2 code should be copyright George Hobbs and Russell Edwards until we decide to change this.

Headers should declare the functions and have documentation! Please avoid globals as much as possible, but sometimes they are required. Use any doxygen markup required to document the interface, ESPECIALLY if it is to be called from outside tempo2.

\code{.cpp}
// use defines to prevent double declaration
#ifndef myHeader_h 
#define myHeader_h
\endcode
\verbatim
/*!
 *  @brief A brief description of the function
 *  @param anInt[in]   description of this parameter
 *  @param str[in]     description of this parameter
 *  @param matrix[out] descrition, note if it is an "output" parameter!
 *
 *  More description if required
 */
\endverbatim
\code{.cpp}
void myFunction(int anInt, const char* str, double** matrix);
#endif
\endcode










