User Guide                                      {#userguide}
==========

Tempo2 User Manual                              {#ug-top}
==================
[TOC]

About tempo2                                    {#ug-about}
------------

Tempo2 is a pulsar timing package, based on the old fortran tempo code to address some shortcomings in that code for high precision pulsar timing. Over the years tempo2 has been expanded my many \ref devs "developers", and has grown to become the premier package for all kinds of pulsar timing experiments.

For more details on pulsar timing in general, you may wish to read the Tempo2 paper series:

+ I. An overview <http://adsabs.harvard.edu/abs/2006MNRAS.369..655H>
+ II. The timing model and precision estimates <http://adsabs.harvard.edu/abs/2006MNRAS.372.1549E>
+ III. Gravitational wave simulation <http://adsabs.harvard.edu/abs/2009MNRAS.394.1945H>

There is also a lot of useful information on the tempo2 wiki, <http://www.atnf.csiro.au/research/pulsar/tempo2/index.php?n=Main.HomePage>. Some of the details are outdated as of 2015, but the general principles are sound. **The wiki is the best place for tutorials and basic introduction to tempo2.**



Terminology and basic usage                     {#ug-basic}
---------------------------
This documentation will focus on providing some basic overview of the many functions of tempo2 and is mostly intended as a reference for those who have mastered the basics. However, for completeness, here we will cover the most basic functions of tempo2.

Tempo2 brings together time-of-arival measurements (ToAs), stored in a `.tim` file, and a pulsar ephemeris stored in a `.par` file to produce the difference between the pulsar ephemeris model and the actual arrival times. This step is generally known as "forming residuals", and depends on having accurate models of the Earth ephemeris and of the clocks used to measure the ToAs. If all is well, these differences will be consistent with the uncertanty in the measurements. This is not generally the case, therefore the second part of tempo2 is a fitting routine that attempts to update the model parameters to get the best-fit model.

The basic usage of tempo2 is to feed in a `.par` and a `.tim` file, form residuals, do the fit and write out the best-fit parameters.

    tempo2 -f example1.par example1.tim -newpar

This will write out `new.par` file with the updated parameters, as well as printing to the console the pre and post-fit parameters. Note any warnings that are printed. One of the side-effects of tempo2 is that it sometimes prints a lot of warnings, some you can ignore and some that you can't, so you have to read them!


If you compiled the `pgplot` plugins, you can run the graphical interface `plk`

    tempo2 -gr plk -f example1.par example1.tim


Running plugins
---------------
There are many, many plugins. Some plugins are better supported than others. To get a list of the plugins you have installed try `tempo2 -h`. The majority of plugins are "graphical" plugins, even if they do not use graphics. This is to do with the way that the plugin is called, rather than anything to do with it being graphical. Graphical plugins are run with the `-gr` option. A few other types of plugins exist:

+ `-gr <plugin_name>` for so-called "graphical" plugins. This is most plugins.
+ `-output <plugin_name>` for "output" plugins, like `general` and `general2`
+ `-fitfunc <plugin_name>` for alternative fit routines. This is likely to be removed in a future release.
+ `-select <plugin_name>` for ToA filtering plugins.


You may have to review the source code if you can't find documentation for the plugin you desire. See the @ref plugin for more details on the avaliable plugins.

