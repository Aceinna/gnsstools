/* @(#)consts.h	1.1 09/16/98 */

/* ======== CONSTANTS FOR PROGRAM TDSOLVE ========== */
/* S. Hilla, 22 Nov 1996 */

#if !defined(_consts_h)
#define _consts_h

#define  MAX_PRN_ID     36
#define  c_in_vac       299792458.0
#define  freq_L1        1575420000.0
#define  freq_L2        1227600000.0
#if !defined(PI)
#  define  PI             3.1415926535898
#endif
#define  roundoff(x)  ( ((x)<0.0) ? ((long)((x)-0.5)) : ((long)((x)+0.5)) )

#define AE84            6378137.0
#define FLAT84          0.00335281066474
#define FLATINV84       298.257223563

#endif /* _consts_h */
