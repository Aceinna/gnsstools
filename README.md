# gnsstools
GNSS Data Analysis Utility Tools

GNSS pseudorange multipath skyplot

This tool is based on the original program developed by Stephen Hilla of National Geodetic Survey. Due to the upgrade of TEQC and GMT, the anticipated plots could not be generated with the original program. Specifically, the current version of TEQC produces COMPACT3 plot files instead of COMPACT plots; and some GMT options have been deprecated. The new development adds functions to adapt to these changes. Moreover, functionality of analyzing observation gaps in the RINEX file has been implemented.

In order to use the program, TEQC, GMT, and Ghostview have to be installed. After generating plot files with TEQC and changing some parmenters in cf2sky.inp file, the program can be run to get the skyplot in .ps format, which can be viewed and converted to other formats with Ghostview. 
