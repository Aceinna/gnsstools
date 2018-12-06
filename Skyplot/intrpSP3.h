// intrpSP3.h

// Class definitions for SP3 file objects
// and prototypes for precise orbit interpolation methods

#if !defined( __SP3FILE__ )

#define __SP3FILE__

#include <fstream>
#include <string>
#include "datetime.h"
using namespace std;
using namespace NGSdatetime;


//========================================== constants


const unsigned short  MAXSVSEPH = 36;
const unsigned short  MAXPARAM  =  8;



//========================================== SP3 File Class

class SP3File
{
    public:
      SP3File();
      SP3File(string pathFilename, ios_base::openmode );
      ~SP3File();

      void setPathFilenameMode(string pathFilename,
                               ios_base::openmode );

      void initHeaderInfo();
      int readHeader();
      int getSVPosVel(DateTime tuser, unsigned short PRNid, double rvec[]);

    private:
      string               pathFilename;
      fstream              fileStream;
      ios_base::openmode   fileMode;

      char                 formatVersion;
      char                 modeFlag;
      DateTime             SP3StartTime;
      DateTime             SP3EndTime;
      unsigned long        numberSP3epochs;
      string               dataUsed;
      string               coordFrame;
      string               orbitType;
      string               sourceAgency;
      unsigned long        gpsWeek;
      double               secsOfWeek;
      double               SP3interval;
      long                 SP3mjd;
      double               SP3fmjd;
      unsigned short       numberSP3svs;
      unsigned short       sp3PRNs[MAXSVSEPH + 1];
      unsigned short       svAccur[MAXSVSEPH + 1];

      unsigned long        lastEpochRead;
      unsigned long        currEpoch;
      unsigned short       numberSVparams;
      double               intrpCoeff[5 + 1][2 + 1][MAXPARAM + 1][MAXSVSEPH + 1];
      double               inputValues[10 + 1][MAXPARAM + 1][MAXSVSEPH + 1];
      unsigned short       numberGoodPRNs;
      unsigned short       numberGoodACCURs;

};

#endif
