// intrpSP3.cpp

// Methods for reading and interpolating SP3 precise ephemerides
#include "stdafx.h"
#include <sstream>
#include <iomanip>
#include "datetime.h"
#include "intrpSP3.h"
using namespace std;
using namespace NGSdatetime;

const double EVCF[5][5] = {
{ 1.7873015873015873e+00,
-0.9359567901234568e+00,
 0.1582175925925926e+00,
-0.9755291005291005e-02,
 0.1929012345679012e-03},
{-0.4960317460317460e+00,
 0.6057098765432098e+00,
-0.1171296296296296e+00,
 0.7605820105820106e-02,
-0.1543209876543210e-03},
{ 0.1206349206349206e+00,
-0.1632716049382716e+00,
 0.4606481481481481e-01,
-0.3505291005291005e-02,
 0.7716049382716048e-04},
{-0.1984126984126984e-01,
 0.2779982363315696e-01,
-0.8796296296296296e-02,
 0.8597883597883598e-03,
-0.2204585537918871e-04},
{ 0.1587301587301587e-02,
-0.2259700176366843e-02,
 0.7523148148148148e-03,
-0.8267195767195767e-04,
 0.2755731922398589e-05} };  // the Everett coefficients


const double MFACT[9] = { 1.0, 3.0, 5.0, 7.0, 9.0,
     2.0, 4.0, 6.0, 18.0 }; // multiplicative factors



//================== SP3File class =============================

SP3File::SP3File() 
{
    pathFilename = "nofilename.out";
    fileMode = ios_base::out;
}

SP3File::SP3File(string inputFilePath, ios_base::openmode mode)
{
    string record;
    pathFilename = inputFilePath;
    fileMode = mode;
    fileStream.open( pathFilename.c_str(), fileMode );
    if( ! fileStream )
    {
	cerr << "Error: unable to open SP3 file in constructor: "
	     << pathFilename << " using mode: " << fileMode << endl;
//        SP3FileException  excep( "Error opening file " );
//	excep.ErrorMessage.append( pathFilename );
//	throw excep;
    }
}

// Destructor
SP3File::~SP3File()
{
    pathFilename = "nofilename.out";
    fileMode = ios_base::out;
    fileStream.close();
}


void SP3File::setPathFilenameMode(string inputFilePath,
                                    ios_base::openmode mode)
{
    pathFilename = inputFilePath;
    fileMode = mode;
    fileStream.open( pathFilename.c_str(), fileMode );
    if( ! fileStream )
    {
	cerr << "Error: unable to open SP3 file in setPathFilenameMode: "
	     << pathFilename << " using mode: " << fileMode << endl;
//        SP3FileException  excep( "Error opening file " );
//	excep.ErrorMessage.append( pathFilename );
//	throw excep;
    }
}



void SP3File::initHeaderInfo()
{
   int i,j,k,l;

   formatVersion = ' ';
   modeFlag = ' ';
   SP3StartTime.SetYMDHMS(1999,1,1,0,0,0.0);
   SP3EndTime.SetYMDHMS(1999,1,1,0,0,0.0);
   numberSP3epochs = 0;
   dataUsed = "";
   coordFrame = "";
   orbitType = "";
   sourceAgency = "";
   gpsWeek = 9999;
   secsOfWeek = 0.0;
   SP3interval = 0.0;
   SP3mjd = 999999;
   SP3fmjd = 0.0;
   numberSP3svs = 0;
   numberGoodPRNs = 0;
   numberGoodACCURs = 0;
   for( i = 0; i < MAXSVSEPH + 1; i++)
   {
     sp3PRNs[i] = 0;
     svAccur[i] = 0;
   }

   lastEpochRead = 0;
   currEpoch = 0;
   numberSVparams = 0;

   for(i = 0; i < 5 + 1;  i++ )
     for(j = 0; j < 2 + 1; j++ )
       for(k = 0; k < MAXPARAM + 1; k++ )
         for(l = 0; l < MAXSVSEPH + 1; l++ )
            intrpCoeff[i][j][k][l] = 0.0;
   for(i = 0; i < 10 + 1; i++ )
     for(j = 0; j < MAXPARAM + 1; j++ )
       for(k = 0; k < MAXSVSEPH + 1; k++ )
            inputValues[i][j][k] = 0.0;

};


int SP3File::readHeader()
{
   string      recType;
   string      inputRec;
   string      temp;
   char        nextFirstChar;
   char        inputRecC[ 256 ];
   long        year, month, day, hour, minute;
//   long        mjd;
   double      sec;
//   double      fmjd;
//   int         numsvs = 0;
   int         indexPRN = 0;
   int         indexACC = 0;
   int         lineLength = 0;
   int         i, prn, accur;
//   int         reply;

   prn = 0;
   initHeaderInfo();

   // position the stream at the beginning of the file
   fileStream.seekg(0);

   while( fileStream.getline(inputRecC, 255, '\n') )
   {
      inputRec = inputRecC;
      recType = inputRec.substr( 0, 2 );
      if( recType[0] == '#' && recType[1] != '#' )
      {
	  formatVersion = recType[1];
//cout << formatVersion << endl;  cin >> reply;
	  if( formatVersion != 'a' && formatVersion != 'A' )
	  {
	      cerr << "This program reads only SP3 files with format-mode: aP "
	      << endl << "This file has format version: " << formatVersion << endl;
	      return( -1 );
          }
	  modeFlag = inputRec[2];
//cout << modeFlag << endl;  cin >> reply;
	  if( modeFlag != 'P' && modeFlag != 'p' )
	  {
	      cerr << "This program reads only SP3 files with format-mode: aP "
	      << endl << "This file has mode flag: " << modeFlag << endl;
	      return( -1 );
          }
          temp = inputRec.substr( 3, 4 );
            year = atol( temp.c_str() );
          temp = inputRec.substr( 8, 2 );
            month = atol( temp.c_str() );
          temp = inputRec.substr(11, 2 );
            day = atol( temp.c_str() );
          temp = inputRec.substr(14, 2 );
            hour = atol( temp.c_str() );

          temp = inputRec.substr(17, 2 );
            minute = atol( temp.c_str() );
          temp = inputRec.substr(20, 11 );
            sec = atof( temp.c_str() );
	  SP3StartTime.SetYMDHMS(year, month, day, hour, minute, sec);
	  SP3EndTime.SetYMDHMS(year, month, day, hour, minute, sec); // this will get added to later

          temp = inputRec.substr(32, 7 );
            numberSP3epochs = (unsigned long) atol( temp.c_str() );
          dataUsed = inputRec.substr(40, 5 );
          coordFrame = inputRec.substr(46, 5 );
          orbitType = inputRec.substr(52, 3 );
          sourceAgency = inputRec.substr(56, 4 );

      }
      else if( recType[0] == '#' && recType[1] == '#' )
      {
          temp = inputRec.substr(3, 4 );
            gpsWeek = (unsigned long) atol( temp.c_str() );
          temp = inputRec.substr(8, 15 );
            secsOfWeek = atof( temp.c_str() );
          temp = inputRec.substr(24, 14 );
            SP3interval = atof( temp.c_str() );
          temp = inputRec.substr(39, 5 );
            SP3mjd = atol( temp.c_str() );
          temp = inputRec.substr(45, 15 );
            SP3fmjd = atof( temp.c_str() );

      }
      else if( recType[0] == '+' && recType[1] == ' ' )
      {
	  if( numberSP3svs == 0 )
	  {
            temp = inputRec.substr(4, 2 );
              numberSP3svs = (unsigned short) atoi( temp.c_str() );

            lineLength = inputRec.size();
	    for( i = 9; i < lineLength; i=i+3 )
	    {
              temp = inputRec.substr( i, 3 );
	      prn = atoi( temp.c_str() );
	      if( prn > 0 )
	      {
	        sp3PRNs[indexPRN + 1] = (unsigned short) prn;
	        indexPRN++;
	      }
            }

          }
	  else
	  {
            lineLength = inputRec.size();
	    for( i = 9; i < lineLength; i=i+3 )
	    {
              temp = inputRec.substr( i, 3 );
	      prn = atoi( temp.c_str() );
	      if( prn > 0 )
	      {
	        sp3PRNs[indexPRN + 1] = (unsigned short) prn;
	        indexPRN++;
	      }
            }

          }

      }
      else if( recType[0] == '+' && recType[1] == '+' )
      {
         lineLength = inputRec.size();
	 for( i = 9; i < lineLength; i=i+3 )
	 {
           temp = inputRec.substr( i, 3 );
	   accur = atoi( temp.c_str() );
	   if( accur > 0 )
	   {
	     svAccur[indexACC + 1] = (unsigned short) accur;
	     if( accur <= 0 )
	     {
	       cerr << "WARNING ! accuracy code ZERO for PRN: " <<
	       sp3PRNs[indexACC + 1] << endl;
               sp3PRNs[indexACC + 1] = (unsigned short)(-1.0 * sp3PRNs[indexACC + 1]);
             }
	     indexACC++;
	   }
         }

      }
      nextFirstChar = (char) fileStream.peek();
      if( nextFirstChar == '*' ) break;  // exit while loop

   } // end of while loop to read all SP3 header records


   numberSVparams = 4;
   numberGoodPRNs = (unsigned short) indexPRN;
   numberGoodACCURs = (unsigned short) indexPRN;
   SP3EndTime = SP3EndTime +
     (double)(numberSP3epochs - 1)*SP3interval/86400.0;

   if( numberGoodPRNs <= 0  ||  numberSP3svs != numberGoodPRNs )
   {
      cerr << "Error reading the PRNs from the header for SP3 file: "
      << endl << pathFilename << endl
      << "Number of SVs expected: " << setw(6) << numberSP3svs << endl
      << "Number of PRNs read in: " << setw(6) << numberGoodPRNs << endl;
      return(1);
   }
   else
     return(0);


} // end of method SP3File::readHeader()

//------------------------------------------------------------------------


int SP3File::getSVPosVel(DateTime tuser, unsigned short PRNid, double rvec[])
{
   double    p, q, p2, q2, sum1, sum2, trun;
   int       i, j, k, l, m, n, jsv;
   long      jmax, jmid, jmin, mbeg, mend, muse;
   bool      foundSV = false;
   char      inputRecC[ 256 ];
   string    inputRec;
   string    temp;

   int  iret;
//   fstream   dout;
//   dout.open("hilla.out",ios_base::out);
//   if( !dout ) cerr << "Cannot open hilla.out file !!! " << endl;


   for( i = 0; i < (4 + numberSVparams); i++)
     rvec[i] = 0.0;

//   dout << "numberSP3svs = " << numberSP3svs << endl;
//   dout << "numberSVparams = " << numberSVparams << endl;
//   dout << "numberSP3epochs = " << numberSP3epochs << endl;
//   dout << "lastEpochRead = " << lastEpochRead << endl;

   for( i = 1; i <= numberSP3svs; i++)
   {
//     dout << "i,prn: " << i << " " << sp3PRNs[i] << endl;
     if( PRNid == sp3PRNs[i] )
     {
       jsv = i;
       foundSV = true;
       break;  // exit loop over SP3svs once a match is found
     }
   }


   if( !foundSV )
   {
     //cerr << "Cannot find PRN " << PRNid << " in SP3 file "
     //<< pathFilename << endl << "Terminate Program." << endl;
     return( -2 );
   }


   trun = (tuser - SP3StartTime)*86400.0;  // diff between DateTimes = days
   trun = trun/SP3interval + 1.0;
//   dout << "trun as epoch span = " << trun << endl;
   jmin = (long)trun - 4;
   jmax = (long)trun + 5;
//   dout << "jmin, jmax: " << jmin << " " << jmax << endl;

   if( jmin < 1 )
   {
       jmin = 1;
       jmax = jmin + 9;
       if( jmax > (long) numberSP3epochs)
       {
         return( -3 );  // not enough data at beginning of SP3 file
       }
   }

   if( jmax > (long) numberSP3epochs )
   {
       jmax = numberSP3epochs;
       jmin = jmax - 9;
       if( jmin < 1 )
       {
	   return( -4 );  // not enough data at end of SP3 file
       }
   }


// move in file
   if( jmax != (long) lastEpochRead   &&
       fabs( (double)lastEpochRead - (trun + 5.0) ) > (1.0/SP3interval) )
   {
	 // move in file backward
	 if( jmax < (long) lastEpochRead )
	 {
	     iret = readHeader();
	     if( iret != 0 ) return( iret );
         }

	 // move in file forward
    if( lastEpochRead > 0 )
    {
         for( k = 1; k <= numberSP3svs; k++ )
	 {
	   for( j = 1; j <= numberSVparams; j++ )
	   {
	     for( i = 1; i <= (int) lastEpochRead - jmin + 1; i++ )
	     {
 		 m = i + 10 - (lastEpochRead - jmin + 1);
//                 dout << "file forward ====  m = " << m << endl;
		 inputValues[i][j][k] = inputValues[m][j][k];
//                 dout << " i,j,k,yye: " << i << " " << j << " " << k <<
//                  " " << setw(30) << setprecision(15) << inputValues[i][j][k];
             }
	   }
	 }
     } // do this section only if lastEpochRead != 0


         // skip over lines to get to the correct epoch in the SP3 file
         for( i = lastEpochRead + 1; i <= jmin - 1; i++ )
	 {
	   for( j = 1; j <= numberSP3svs + 1; j++ )
	   {
             //if( fileStream.getline(inputRecC, 255, '\n') == NULL )
		   if (!fileStream.getline(inputRecC,255,'\n'))
	     {
		 cerr << "Error skipping lines in the SP3 file." << endl;
		 return( -1 );
             }
//             dout << " aline = " << inputRecC << endl;
	   }
	 }

	 // read the data into the input values
	 if( (long) (lastEpochRead + 1) > jmin)
	   m = lastEpochRead + 1;
	 else
	   m = jmin;

         for( i = m; i <= jmax; i++ )
	 {
	   l = i - jmin + 1;
           if( !fileStream.getline(inputRecC, 255, '\n') )//== NULL )
	   {
	      cerr << "Error skipping time tag line in the SP3 file." << endl;
	      return( -1 );
           }
	   for( k = 1; k <= numberSP3svs; k++ )
	   {
              if( !fileStream.getline(inputRecC, 255, '\n'))// == NULL )
              {
               cerr << "Error reading input values in the SP3 file." << endl;
	       return( -1 );
              }
	      inputRec = inputRecC;
	      for( j = 1; j <= numberSVparams; j++ )
	      {
		  temp = inputRec.substr( (4 + ((j - 1)*14)), 14);
		  inputValues[l][j][k] = atof( temp.c_str() );
                  if( inputValues[l][j][k] > 999999.9 )
                     inputValues[l][j][k] = 0.0;   // do for bad clocks
// dout << "read new yye (i,j,k,yye): " << setw(8) << (i - jmin + 1) << " "
// << j << " " << k << " " << setw(30) << setprecision(15) <<
// inputValues[l][j][k] << endl;
	      }

	   }
	 }
	 lastEpochRead = jmax;
//         dout << "new jlaste: " << setw(8) << lastEpochRead << endl;


	 // interpolate terms for position and velocity
	 for( i = 1; i <= numberSP3svs; i++ )
	 {
	   n = 4;
	   if( numberSVparams < 4 ) n = numberSVparams;
	   for( j = 1; j <= n; j++ )
	   {
	     for( k = 1; k <= 2; k++ )
	     {
                jmid = 4 + k;
		for( l = 1; l <= 5; l++ )
		{
                    intrpCoeff[l][k][j][i] = 0.0;
		    for( m = 5; m >= 2; m-- )
		    {
		       intrpCoeff[l][k][j][i] = intrpCoeff[l][k][j][i] +
		       EVCF[m - 1][l - 1] *
		       (inputValues[jmid + m - 1][j][i]
		            + inputValues[jmid - m + 1][j][i]);
// dout << " 1st iloop l,k,j,i,m: " << l << " " << k << " " << j << " " << i
//      << " " << m << endl;
                    }
		    intrpCoeff[l][k][j][i] = intrpCoeff[l][k][j][i] +
		    EVCF[1 - 1][l - 1] * inputValues[jmid][j][i];

		    intrpCoeff[l][k][j + numberSVparams][i] =
		      intrpCoeff[l][k][j][i] * MFACT[l - 1]/SP3interval;
                }
             }

           }
	   for( j = 5; j <= numberSVparams; j++ )
	   {
	     for( k = 1; k <= 2; k++ )
	     {
                jmid = 4 + k;
		for( l = 1; l <= 5; l++ )
		{
                    intrpCoeff[l][k][j][i] = 0.0;
		    for( m = 5; m >= 2; m-- )
		    {
		       intrpCoeff[l][k][j][i] = intrpCoeff[l][k][j][i] +
		       EVCF[m - 1][l - 1] *
		       (inputValues[jmid + m - 1][j][i]
		            + inputValues[jmid - m + 1][j][i]);
                    }
		    intrpCoeff[l][k][j][i] = intrpCoeff[l][k][j][i] +
		    EVCF[1 - 1][l - 1] * inputValues[jmid][j][i];
                }
             }

           }


       } // end of loop over satellites


     } // end of if test for jmax != lastEpochRead && ...

     // initialize interpolation parameters
     p = trun - lastEpochRead + 5;
     q = 1.0 - p;
     p2 = p * p;
     q2 = q * q;

//     dout << "p  = " << p  << endl;
//     dout << "q  = " << q  << endl;
//     dout << "p2 = " << p2 << endl;
//     dout << "q2 = " << q2 << endl;

     mbeg = 1;
     mend = numberSVparams + 4;   // do all: x,y,z,clk, xdot,ydot,zdot,clk-dot
//     dout << "mbeg,mend: " << mbeg << " " << mend << endl;

     // interpolate
     for( m = mbeg; m <= mend; m++ )
     {
        if( m <= 4 )
	  muse = m;
	else if( m > 4  &&  m <= numberSVparams )
	  muse = m + 4;

        if( m > numberSVparams  &&  m <= (numberSVparams + 4) )
	  muse = m - numberSVparams + 4;


        sum1 = 0.0;
	sum2 = 0.0;
	for( j = 5; j >= 2; j-- )
	{
//           dout << "3rd iloop j,m,jsv,muse: " << j << " " << m << " " <<
//           " " << jsv << " " << muse << endl;
	   sum1 = p2*(sum1 + intrpCoeff[j][2][m][jsv]);
	   sum2 = q2*(sum2 + intrpCoeff[j][1][m][jsv]);
   	}

	if( m <= numberSVparams )
	  rvec[muse - 1] = p * (intrpCoeff[1][2][m][jsv] + sum1 ) +
	    q * (intrpCoeff[1][1][m][jsv] + sum2 );
	else
	  rvec[muse - 1] = (intrpCoeff[1][2][m][jsv] + sum1 ) -
	    (intrpCoeff[1][1][m][jsv] + sum2 );

   }

//   dout << setw(30) << setprecision(15) << endl << endl;
//   dout << "inside rvec[0] = " << rvec[0] << endl;
//   dout << "inside rvec[1] = " << rvec[1] << endl;
//   dout << "inside rvec[2] = " << rvec[2] << endl;
//   dout << "inside rvec[3] = " << rvec[3] << endl;
//   dout << "inside rvec[4] = " << rvec[4] << endl;
//   dout << "inside rvec[5] = " << rvec[5] << endl;
//   dout << "inside rvec[6] = " << rvec[6] << endl;
//   dout << "inside rvec[7] = " << rvec[7] << endl;

//   dout.close();
   return(0);

} // end of method SP3File::getSVPosVel()

