/***************************************************************
* 
*  FILE   :  parse.c
*  DATE   :  March 31, 1993
*  AUTHOR :  Mark Heng
*  
*  DESCRIPTION :
*    This is a source code listing of a program designed to serve as
*    an example of how to decode and parse :
*      1) the Motorola PVT6 ephemeris message (in either ASCII or binary
*         formats) into the ephemeris components described in the 
*         ICD-GPS-200 for  use in post-processing applications
*      2) the Motorola PVT6 position/status/data message (in binary
*         format only) into its components as described in the Motorola
*         GPS Technical Reference Manual
*
*    This program is intended for use as an application note and is 
*    provided free of charge to customers in order to aid in the use of
*    the Motorola PVT6 Core GPS receiver.  The program coded herein has
*    been compiled, tested, and debugged by the author and to the best of
*    author's knowledge the program performs its function correctly.
*    Motorola, however, does not assume responsibility for any errors or
*    oversights of this program in use by customers in any end use application.
*    
*  NOTES :
*    Other than the procedure main() all procedures in this file
*    are listed in alphabetical order.
*
***************************************************************/

#include <stdio.h>
#include <math.h>

/**/
/********* type definitions *********/

typedef char  ONEBYTE;  /* data type allowing for the use of  8 bits. */
typedef short TWOBYTE;  /* data type allowing for the use of 16 bits. */
typedef long  FOURBYTE; /* data type allowing for the use of 32 bits. */

typedef unsigned char  UNSIGNED_ONEBYTE;          /* unsigned  8 bits */
typedef unsigned short UNSIGNED_TWOBYTE;          /* unsigned 16 bits */
typedef unsigned long  UNSIGNED_FOURBYTE;         /* unsigned 32 bits */

/********* function prototypes *********/

void      main();
void      Convert_Ephemeris();
void      Ephemeris_Decode_Example(UNSIGNED_ONEBYTE eg_type);
TWOBYTE   Hex(ONEBYTE ASCIIvalue);
void      Pos_Status_Data_Decode_Example() ;

/********* constant definitions *********/

#define EG_ASCII  0
#define EG_BINARY 1

#define PI          3.1415926535898      /* GPS value of Pi */
#define WE          7.2921151467E-5      /* WGS84 earths rotation rate (rads/sec) */

#define SQRMU       1.99649818432174E7   /* Square root of MU.               */
#define FCONST     -4.442807633E-10      /* FCONST is used in relativistic clock
                                            correction.  Defined on pg. 73 of
                                            ICD-GPS-200.                     */

#define TWOM5       0.03125000000000     /* 2^-5  */
#define TWOM11      (TWOM5  /  64.0) 
#define TWOM19      (TWOM11 / 256.0) 
#define TWOM20      (TWOM19 /   2.0)
#define TWOM21      (TWOM20 /   2.0)
#define TWOM23      (TWOM21 /   4.0)
#define TWOM24      (TWOM23 /   2.0)
#define TWOM27      (TWOM24 /   8.0)
#define TWOM29      (TWOM23 /  64.0)
#define TWOM30      (TWOM29 /   2.0)
#define TWOM31      (TWOM30 /   2.0)
#define TWOM33      (TWOM31 /   4.0)
#define TWOM38      (TWOM33 /  32.0)
#define TWOM43      (TWOM38 /  32.0)
#define TWOM50      (TWOM43 / 128.0)
#define TWOM55      (TWOM50 /  32.0)
#define TWOM3125    (0.25 * TWOM31)      

#define TWOP11      2048.0               /* 2^11  */
#define TWOP12      (TWOP11 * 2.0)
#define TWOP14      (TWOP12 * 4.0)
#define TWOP16      (TWOP14 * 4.0)

#define TWOP12_INTEGER   4096            /* 2^12  */
#define TWOP16_INTEGER  65536            /* 2^16  */

#define PITWOM31    (PI * TWOM31)        
#define PITWOM43    (PI * TWOM43)

/********* struct and global data definitions *********/

UNSIGNED_ONEBYTE hex_ephemeris[24][3] ; 
    /* [word][byte] : 24 words (3 bytes per word) corresponding to words 3-10
       of subframes 1-3 of the satellite nav message. */

struct T_FLOAT_EPHEMERIS {
   UNSIGNED_ONEBYTE iode;

   /* Field name abbreviations were extracted from ICD-GPS-200. */
   double toe;                    /* Time of ephemeris.                           */
   double m0;                     /* Mean anomaly at refernce time.               */
   double delta_n;                /* Mean motion difference from computed value.  */
   double e;                      /* Eccentricity.                                */
   double sqrt_a;                 /* Square root of semi-major axis.              */
   double omega_0;                /* Longitude of ascending node of orbit plane
                                     at weekly epoch.                             */
   double i0;                     /* Inclination angle at reference time.         */
   double w;                      /* Argument of perigee.                         */
   double omega_dot;              /* Rate of right ascension.                     */
   double i_dot;                  /* Rate of inclination angle.                   */
   double cuc;                    /* Amplitude of the cosine harmonic correction
                                     term to the argument of latitude.            */
   double cus;                    /* Amplitude of the sine harmonic correction
                                     term to the argument of latitude.            */
   double crc;                    /* Amplitude of the cosine harmonic correction
                                     term to the orbit radius.                    */
   double crs;                    /* Amplitude of the sine harmonic correction
                                     term to the orbit radius.                    */
   double cic;                    /* Amplitude of the cosine harmonic correction
                                     term to the angle of inclination.            */
   double cis;                    /* Amplitude of the sine harmonic correction
                                     term to the angle of inclination.            */
   double tgd;                    /* Estimated group delay differential.          */
   double toc;                    /* Clock data reference time.                   */
   double af2;                    /* Polynomial coefficient (SV clock correction) */
   double af1;                    /* Polynomial coefficient (SV clock correction) */
   double af0;                    /* Polynomial coefficient (SV clock correction) */
   double n;                      /* Corrected mean motion.                       */
   double e_sqrt;                 /* Square root of eccentricity^2 subtracted
                                     from 1.                                      */
   double wk1;                    /* Omega_dot - WE.                              */
   double wk0;                    /* Omega_0 - WE * toe.                          */
   double sin_w;                  /* Sine w.                                      */
   double cos_w;                  /* Cosine w.                                    */
   double a;                      /* Semi-major axis.                             */
   double sqrta_e_fconst;         /* Square root of semi-major axis *
                                       eccentricity * FCONST.                     */
} float_ephemeris ;

struct T_RECEIVER_CHANNELS {
   UNSIGNED_ONEBYTE	svid;
   UNSIGNED_ONEBYTE	mode;
   UNSIGNED_ONEBYTE	strength;
   UNSIGNED_ONEBYTE	flags;
} ;

struct T_GEODETIC  {
   TWOBYTE          degrees ;
   UNSIGNED_TWOBYTE minutes ;
   double           seconds ;
} ;
  
#define NUM_CHANNELS 6 

struct T_POS_CHAN_STATUS {
   UNSIGNED_ONEBYTE	month;
   UNSIGNED_ONEBYTE	day;
   UNSIGNED_TWOBYTE	year;
   UNSIGNED_ONEBYTE	hours;
   UNSIGNED_ONEBYTE	minutes;
   double   	        seconds;
   struct T_GEODETIC   	latitude;
   struct T_GEODETIC    longitude;
   double 	        datum_height;  /* meters */
   double	        msl_height;    /* meters */
   double		velocity;      /* m/sec */
   double		heading;       /* degrees */
   double       	current_dop;
   UNSIGNED_ONEBYTE	dop_type;
   UNSIGNED_ONEBYTE	visible_sats;
   UNSIGNED_ONEBYTE	sats_tracked;
   struct T_RECEIVER_CHANNELS channel[NUM_CHANNELS] ;
   UNSIGNED_ONEBYTE	rcvr_status;
} pos_chan ;

FILE *EPHfile ;
FILE *POSfile ;

/**/
/*******************************************************************
*
* FUNCTION NAME:   main
* 
* DESCRIPTION:
*   This procedure calls two routines :
*     1) a routine that parses and decodes the ephemeris record (this
*        routine is called twice, once for each data format), and
*     2) a routine that parses and decodes the position/status/data
*        record.
*
********************************************************************/

void main()
{
   Ephemeris_Decode_Example(EG_ASCII) ;
   Ephemeris_Decode_Example(EG_BINARY) ;
   Pos_Status_Data_Decode_Example() ;
}

/**/
/*******************************************************************
*
* FUNCTION NAME:       Convert_Ephemeris
* 
* DESCRIPTION:
*   This procedure converts the raw binary ephemeris data into a
*   representation usable for post-processing applications.  The
*   equations coded in this procedure are taken from the ICD-GPS-200.
*   (the GPS NAVSTAR ICD).
*
********************************************************************/

void Convert_Ephemeris()
{
/******* this MACRO sign extends a negative valued hex byte into its
   TWOBYTE (integer) equivalent *******/
#define MACRO_SET_UPPER(word_to_set) \
   if ( word_to_set & 0x0080 )  word_to_set |= 0xFF00 ;
   
/******* this MACRO takes two hex bytes and merges them into a single
   TWOBYTE (integer) value *******/
#define MACRO_MAKE_EPH_WORD(e_wd, w1, b1, w2, b2) \
   e_wd = ( ( ( UNSIGNED_TWOBYTE ) hex_ephemeris[w1][b1] << 8 ) | \
              ( UNSIGNED_TWOBYTE ) hex_ephemeris[w2][b2] );
   
   UNSIGNED_TWOBYTE word1, word2;  /* Temporary storage for words built. */
   struct T_FLOAT_EPHEMERIS *eph ;

   eph                      = &float_ephemeris ;
   
   /* Time of ephemeris */
   MACRO_MAKE_EPH_WORD(word1,  15, 0, 15, 1 );
   eph->toe = (double) ( ((UNSIGNED_FOURBYTE) word1) * 16 ) ;
   
   /* Mean anomaly at reference time ( m0 ). */
   MACRO_MAKE_EPH_WORD(word1,  9, 2, 10, 0 );
   MACRO_MAKE_EPH_WORD(word2, 10, 1, 10, 2 );
   eph->m0 = ( ( ( FOURBYTE ) word1 << 16 ) | word2 ) * PITWOM31;

   /* Mean motion difference from computed value ( delta_n ). */
   MACRO_MAKE_EPH_WORD(word1,  9, 0, 9, 1 );
   eph->delta_n = ( TWOBYTE ) word1 * PITWOM43;
   
   /* Eccentricity ( e ). */
   MACRO_MAKE_EPH_WORD(word1,  11, 2, 12, 0 );
   MACRO_MAKE_EPH_WORD(word2,  12, 1, 12, 2 );
   eph->e = ( ( ( UNSIGNED_FOURBYTE ) word1 << 16 ) | word2 ) * TWOM33;
   
   /* Square root of the semi-major axis ( sqrt_a ). */
   MACRO_MAKE_EPH_WORD(word1,  13, 2, 14, 0 );
   MACRO_MAKE_EPH_WORD(word2,  14, 1, 14, 2 );
   eph->sqrt_a = ( ( ( UNSIGNED_FOURBYTE ) word1 << 16 ) | word2 ) * TWOM19;
   
   /* Longitude of ascending mode of orbit plane at weekly epoch ( omega_0 ). */
   MACRO_MAKE_EPH_WORD(word1,  16, 2, 17, 0 );
   MACRO_MAKE_EPH_WORD(word2,  17, 1, 17, 2 );
   eph->omega_0 = ( ( ( FOURBYTE ) word1 << 16 ) | word2 ) * PITWOM31;
   
   /* Inclination angle at reference time ( i0 ). */
   MACRO_MAKE_EPH_WORD(word1,  18, 2, 19, 0 );
   MACRO_MAKE_EPH_WORD(word2,  19, 1, 19, 2 );
   eph->i0 = ( ( ( FOURBYTE ) word1 << 16 ) | word2 ) * PITWOM31;
   
   /* Argument of perigee ( w ). */
   MACRO_MAKE_EPH_WORD(word1,  20, 2, 21, 0 );
   MACRO_MAKE_EPH_WORD(word2,  21, 1, 21, 2 );
   eph->w = ( ( ( FOURBYTE ) word1 << 16 ) | word2 ) * PITWOM31;
   
   /* Rate of right ascension ( omega_dot ). */
   word1 = ( UNSIGNED_TWOBYTE ) hex_ephemeris[22][0];
   MACRO_SET_UPPER(word1) ;   /* Set upper bits if negative. */
   MACRO_MAKE_EPH_WORD(word2,  22, 1, 22, 2 );
   eph->omega_dot = ( ( ( FOURBYTE ) word1 << 16 ) | word2 ) * PITWOM43;
   
   /* Rate of inclination angle ( i_dot ). */
   MACRO_MAKE_EPH_WORD(word1,  23, 1, 23, 2 ) ;
   word1 = word1 >> 2 ;
   if ( word1 & 0x2000 )     /* Set upper bits if negative. */
      word1 |= 0xC000;
   eph->i_dot = ( TWOBYTE ) word1 * PITWOM43;
   
   /* Amplitude of the cosine harmonic correction term to the argument of latitude (cuc). */
   MACRO_MAKE_EPH_WORD(word1,  11, 0, 11, 1 );
   eph->cuc = ( TWOBYTE ) word1 * TWOM29;
   
   /* Amplitude of the sine harmonic correction term to the argument of latitude (cus). */
   MACRO_MAKE_EPH_WORD(word1,  13, 0, 13, 1 );
   eph->cus = ( TWOBYTE ) word1 * TWOM29;
   
   /* Amplitude of the cosine harmonic correction term to the orbit radius (crc). */
   MACRO_MAKE_EPH_WORD(word1,  20, 0, 20, 1 );
   eph->crc = ( TWOBYTE ) word1 * TWOM5;
   
   /* Amplitude of the sine harmonic correction term to the orbit radius (crs). */
   MACRO_MAKE_EPH_WORD(word1,  8, 1, 8, 2 );
   eph->crs = ( TWOBYTE ) word1 * TWOM5;
   
   /* Amplitude of the cosine harmonic correction term to the angle of inclination (cic). */
   MACRO_MAKE_EPH_WORD(word1,  16, 0, 16, 1 );
   eph->cic = ( TWOBYTE ) word1 * TWOM29;
   
   /* Amplitude of the sine harmonic correction term to the angle of inclination (cis). */
   MACRO_MAKE_EPH_WORD(word1,  18, 0, 18, 1 );
   eph->cis = ( TWOBYTE ) word1 * TWOM29;
   
   /* Estimated group delay differential ( tgd ). */
   word1 = ( UNSIGNED_TWOBYTE ) hex_ephemeris[4][2];
   MACRO_SET_UPPER(word1) ;   /* Set upper bits if negative. */
   eph->tgd = ( TWOBYTE ) word1 * TWOM31;
   
   /* Clock data reference time ( toc ). */
   MACRO_MAKE_EPH_WORD(word1,  5, 1, 5, 2 );
   eph->toc = word1 * 16.0;
   
   /* Polynomial coefficient ( af2 ). */
   word1 = ( UNSIGNED_TWOBYTE ) hex_ephemeris[6][0];
   MACRO_SET_UPPER(word1) ;   /* Set upper bits if negative. */
   eph->af2 = ( TWOBYTE ) word1 * TWOM55;
   
   /* Polynomial coefficient ( af1 ). */
   MACRO_MAKE_EPH_WORD(word1,  6, 1, 6, 2 );
   eph->af1 = ( TWOBYTE ) word1 * TWOM43;
   
   /* Polynomial coefficient ( af0 ). */
   word1 = ( UNSIGNED_TWOBYTE ) hex_ephemeris[7][0];
   MACRO_SET_UPPER(word1) ;   /* Set upper bits if negative. */
   MACRO_MAKE_EPH_WORD(word2,  7, 1, 7, 2 ) ;
   word2 = word2 & 0xFFFC ; /* clear 2 bits */
   eph->af0 = ( ( ( FOURBYTE ) word1 << 16 ) | word2 ) * TWOM3125;
   
   /* Semi-major axis ( a ). */
   eph->a = eph->sqrt_a * eph->sqrt_a;
   
   /* Corrected mean motion ( n ). */
   eph->n = SQRMU / ( eph->a * eph->sqrt_a ) + eph->delta_n;
   
   /* Square root of eccentricity^2 subtracted from 1 ( esqrt ). */
   eph->e_sqrt = sqrt ( 1.0 - eph->e * eph->e );
   
   /* Omega_dot - WE ( wk1 ). */
   eph->wk1 = eph->omega_dot - WE;
   
   /* Omega0 - WE * toe ( wk0 ). */
   eph->wk0 = eph->omega_0 - WE * eph->toe;
   
   /* Sine w ( sin_w ). */
   eph->sin_w = sin ( eph->w );
   
   /* Cosine w ( cos_w ). */
   eph->cos_w = cos ( eph->w );
   
   /* Square root of semi-major axis * eccentricity * FCONST. */
   eph->sqrta_e_fconst = eph->sqrt_a * eph->e * FCONST;
   
   /* Issue of data, ephemeris. */
   eph->iode = hex_ephemeris[23][0];
}

/**/
/*******************************************************************
*
* FUNCTION NAME:       Ephemeris_Decode_Example
* 
* DESCRIPTION:
*   This procedure parses and decodes an ephemeris record received from
*   the Motorola PVT6.
*
*   The input to this procedure is a file (either ASCII or binary format)
*   containing the EPHEMERIS DATA OUTPUT MESSAGE (as described in the
*   Motorola GPS Technical Reference Manual.  Both files can be generated
*   by the Motorola PVT6 Controller program.
*
*   The output of this file is the ephemeris message data converted to
*   floating point format (as described in the ICD-GPS-200) useful for
*   post processing applications.
*
*   The flow of this procedure is as follows :
*     IF example type is ASCII :
*       - open ASCII file containing 150 characters representing the message
*         header (@@Bf), the svid descriptor, and the 72 bytes of ephemeris data
*       - read the ASCII ephemeris data only into storage (ascii_buffer[])
*       - convert the ASCII data to hex data and store to hex_ephemeris[] array
*       - close the ASCII file
*       - convert hex_ephemeris[] array according to the ICD-GPS-200
*     ELSE IF example type is BINARY :
*       - open BINARY file containing 80 hex bytes
*       - read the BINARY data into hex_ephemeris[] array
*       - close the BINARY file
*       - convert hex_ephemeris[] array according to the ICD-GPS-200
*
*********************************************************************/

void Ephemeris_Decode_Example(UNSIGNED_ONEBYTE eg_type)
{
#define EPH_NUM_CHARS 144   /* (72 bytes * 2 nibbles/byte) */
#define EPH_NUM_WORDS  24   /* at 3 bytes/word */
   
   /* local variables */
      ONEBYTE tempchar ;
      UNSIGNED_ONEBYTE char_ctr, byte_ctr, word_ctr, i, base ;
      ONEBYTE ascii_buffer[500] ;
      ONEBYTE *char_ptr ;
      
   /* open ephemeris file */
   
      if ( eg_type == EG_ASCII ) EPHfile = fopen("eph.ascii","r") ;
      else                       EPHfile = fopen("eph.bin","r") ;
      
   /* parse data into temporary storage data buffer (skip 1st 6 chars since
      they are not ephemeris data information) */
   
      if ( eg_type == EG_ASCII )  {
         for ( char_ctr = 0 ; char_ctr < 6 ; char_ctr++ )
            tempchar = fgetc(EPHfile) ;
         for ( char_ctr = 0 ; char_ctr < EPH_NUM_CHARS ; char_ctr++ ) 
            ascii_buffer[char_ctr] = fgetc(EPHfile) ;
	 
         /* pack data into hex storage structure */
            for ( i = 0 ; i < EPH_NUM_WORDS ; i++ ) {
               base = i * 6 ;
               hex_ephemeris[i][0] = 
	          (ONEBYTE) (((Hex(ascii_buffer[base])) << 4) | Hex(ascii_buffer[base+1]));
               hex_ephemeris[i][1] = 
   	          (ONEBYTE) (((Hex(ascii_buffer[base+2])) << 4) | Hex(ascii_buffer[base+3]));
               hex_ephemeris[i][2] = 
	          (ONEBYTE) (((Hex(ascii_buffer[base+4])) << 4) | Hex(ascii_buffer[base+5]));
            }
      }
      else {
         for ( byte_ctr = 0 ; byte_ctr < 5 ; byte_ctr++ )
            tempchar = fgetc(EPHfile) ;
         for ( i = 0 ; i < EPH_NUM_WORDS ; i++ ) {
            hex_ephemeris[i][0] = fgetc(EPHfile) ;
            hex_ephemeris[i][1] = fgetc(EPHfile) ;
            hex_ephemeris[i][2] = fgetc(EPHfile) ;
         }
      }
      
   /* close ephemeris file */
      
      fclose(EPHfile) ;
      
   /* decode from temporary storage into floating point array */
      
      Convert_Ephemeris() ;
}

/**/
/***************************************************************
*
* FUNCTION NAME:	Hex
*
* DESCRIPTION:
*   This procedure converts an ASCII value to its hex equivalent.
*
***************************************************************/

TWOBYTE Hex(ONEBYTE ASCIIvalue)
{
   /* 0 is 0x30 ->   9 is 0x39 */
   /* A is 0x41 ->   F is 0x46 */
   /* a is 0x61 ->   f is 0x66 */
   
  if      (((TWOBYTE) ASCIIvalue) < 0x40)   return(((TWOBYTE) ASCIIvalue) - 0x30);
  else if (((TWOBYTE) ASCIIvalue) < 0x60)   return(((TWOBYTE) ASCIIvalue) - 0x37);
  else                                      return(((TWOBYTE) ASCIIvalue) - 0x57);
}   

/**/
/*******************************************************************
*
* FUNCTION NAME:       Pos_Status_Data_Decode_Example()
* 
* DESCRIPTION:
*   This procedure parses and decodes a position/status/data record
*   received from the Motorola PVT6.
*
*   The input to this procedure is a binary format file containing the
*   POSITION/STATUS/DATA OUTPUT MESSAGE as described in the Motorola
*   GPS Technical Reference Manual.  This file can be generated by the
*   Motorola PVT6 Controller program.
*
*   The output of this file is a data structure containing each of the 
*   message components in their most common scale factors (eg, velocity
*   in m/sec, heading in degrees).
*
*   The flow of this procedure is as follows :
*    - open BINARY file containing 68 bytes
*    - read the BINARY data from the file and convert directly into
*      the final data structure 
*    - close the BINARY file
*
***************************************************************/

void Pos_Status_Data_Decode_Example()
{
#define MSECS_TO_DEGREES ( ( 1.0 / 1000.0 ) / 3600.0 )
   
  UNSIGNED_ONEBYTE  i ;
  UNSIGNED_ONEBYTE  tempchar ;
  UNSIGNED_FOURBYTE tempu4byte ;
  FOURBYTE temps4byte ;
  UNSIGNED_FOURBYTE nsecs ;
  double degrees, minutes, seconds ;
  
  /* open the binary file */
  
     POSfile = fopen("pos.bin","r") ;
  
  /* skip first 4 bytes (@@Ba) */
  
     for ( i = 0 ; i < 4 ; i++ )
        tempchar = fgetc(POSfile) ;
  
  /* read and scale the rest of the data */
     
     pos_chan.month = fgetc(POSfile) ;
     pos_chan.day   = fgetc(POSfile) ;
     
     tempchar = fgetc(POSfile) ;
     pos_chan.year = ( tempchar << 8 ) + fgetc(POSfile) ;
     
     pos_chan.hours    = fgetc(POSfile) ;
     pos_chan.minutes  = fgetc(POSfile) ;
     
     tempchar   = fgetc(POSfile) ;  /* integer seconds */
     tempu4byte = fgetc(POSfile) ;
     tempu4byte = ( tempu4byte << 8 ) + fgetc(POSfile) ;
     tempu4byte = ( tempu4byte << 8 ) + fgetc(POSfile) ;
     tempu4byte = ( tempu4byte << 8 ) + fgetc(POSfile) ;
     pos_chan.seconds = (double) tempchar + ( ( (double) tempu4byte ) / 1.0E+9 );
     
     temps4byte = fgetc(POSfile) ;
     temps4byte = ( temps4byte << 8 ) + fgetc(POSfile) ;
     temps4byte = ( temps4byte << 8 ) + fgetc(POSfile) ;
     temps4byte = ( temps4byte << 8 ) + fgetc(POSfile) ;
     degrees = (double) temps4byte * MSECS_TO_DEGREES ;
     
     pos_chan.latitude.degrees = (TWOBYTE) degrees ;
     if ( degrees < 0 )
	degrees = fabs ( degrees ) ;
     minutes =  ( degrees - (TWOBYTE) degrees ) * 60.0 ; 
     pos_chan.latitude.minutes = (TWOBYTE) ( minutes ) ;
     pos_chan.latitude.seconds = ( minutes - (TWOBYTE) minutes ) * 60.0 ;
	
     temps4byte = fgetc(POSfile) ;
     temps4byte = ( temps4byte << 8 ) + fgetc(POSfile) ;
     temps4byte = ( temps4byte << 8 ) + fgetc(POSfile) ;
     temps4byte = ( temps4byte << 8 ) + fgetc(POSfile) ;
     degrees = (double) temps4byte * MSECS_TO_DEGREES ;
     
     pos_chan.longitude.degrees = (TWOBYTE) degrees ;
     if ( degrees < 0 )
	degrees = fabs ( degrees ) ;
     minutes = ( degrees - (TWOBYTE) degrees ) * 60.0 ; 
     pos_chan.longitude.minutes = (TWOBYTE) ( minutes ) ;
     pos_chan.longitude.seconds = ( minutes - (TWOBYTE) minutes ) * 60.0 ;
	
     temps4byte = fgetc(POSfile) ;
     temps4byte = ( temps4byte << 8 ) + fgetc(POSfile) ;
     temps4byte = ( temps4byte << 8 ) + fgetc(POSfile) ;
     temps4byte = ( temps4byte << 8 ) + fgetc(POSfile) ;
     pos_chan.datum_height = (double) temps4byte / 100.0 ;
     
     temps4byte = fgetc(POSfile) ;
     temps4byte = ( temps4byte << 8 ) + fgetc(POSfile) ;
     temps4byte = ( temps4byte << 8 ) + fgetc(POSfile) ;
     temps4byte = ( temps4byte << 8 ) + fgetc(POSfile) ;
     pos_chan.msl_height = (double) temps4byte / 100.0 ;
     
     tempchar = fgetc(POSfile) ;
     pos_chan.velocity = (double) ( ( tempchar << 8 ) + fgetc(POSfile) ) / 100.0 ;
     
     tempchar = fgetc(POSfile) ;
     pos_chan.heading = (double) ( ( tempchar << 8 ) + fgetc(POSfile) ) / 10.0 ;
     
     tempchar = fgetc(POSfile) ;
     pos_chan.current_dop = (double) ( ( tempchar << 8 ) + fgetc(POSfile) ) / 10.0 ;
     
     pos_chan.dop_type     = fgetc(POSfile) ;
     pos_chan.visible_sats = fgetc(POSfile) ;
     pos_chan.sats_tracked = fgetc(POSfile) ;

     for (i = 0; i < NUM_CHANNELS; i++)   {
	pos_chan.channel[i].svid     = fgetc(POSfile);
	pos_chan.channel[i].mode     = fgetc(POSfile);
	pos_chan.channel[i].strength = fgetc(POSfile);
	pos_chan.channel[i].flags    = fgetc(POSfile);
     }
     pos_chan.rcvr_status = fgetc(POSfile);

  /* close the binary file */
  
     fclose(POSfile) ;
}
