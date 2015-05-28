/*
 * This function will be called as follows:
 *
 * makeMultipleNewrayInC   freq    num_rays   center_lat ...
 *                       delt_low  delt_high  final_TG 
 *			lat_spread
 *
 * Rays are launched at 1000 km altitude and have a spread of 
 * lat_spread degrees in latitude, symmetrical about center_lat.  
 * NB - there are actually num_rays+1 rays, so that eg. if we 
 * have 10 rays between lat=20 and =30, so that they don't go 
 * 20-29, but more like what we'd expect, 20-30.
 *
 * This file differs from makeMultipleNewrayInC in that it models 
 * the smooth magnetosphere from the Carpenter & Anderson '92 model
 * with the parameters:
 * R = 90  (13 month average sunspot number)
 * Kp_max = 4 (maximumn Kp value for previous 24hours)
 * d = 0 (day number in the year)
 * t = 2 (time of day)
 * 
 */ 





// Includes:
// ---------

#include <stdio.h>
#include <stdlib.h>
#include <time.h>





// Defines:
// --------

#define		RANGE		10






int main(int argc, char *argv[] )
{
  FILE *filePtr; 
  char *outputFileName;
  int i, num_rays;
  long freq;
  double center_lat, delt_low, delt_high, final_TG;
  double DELT_offset;

  double INTERA, NUMRES, NSUPPR, SPELAT;
  double DISTRE, LATITUD;
  double NUM, KSKIP, MODE, KOUNT, KDUCTS, KTAPE, REFALT;
  double DSRRNG, DSRLAT, DSRDENS;
  double EGFEQ, THERM, HM, ABSB, RELB;
  double RBASE, ANE0, ALPHA0_2, ALPHA0_3, ALPHA0_4;
  double RZERO, SCBOT, RSTOP, RDIV, HMIN;
  double LK, EXPK, DDK, RCONSN, SCR;
  double L0, DEF, DD, RDUCLN, HDUCLN, RDUCUN, HDUCUN, RDUCLS;
  double HDUCLS, RDUCUS, HDUCUS, SIDEDU;
  double PSTALT, PALT1, PALT2, PLATIT, PSTLAT, PLAT1, PLAT2, PALTIT;
  double start_lat, lat_incr, delta_spread, delta_incr;
  double FKC, X0, LATITU, DELT, PSI, TINITI_P, TINITI_G, TGFINA;
  double KPLOT, XSCALE, YSCALE, FLBOT, FLTOP;
  double LP_1, LP_2, LP_3, LP_4, LP_5;
  double LAT_SPREAD;


  
  if(argc != 8) {
    printf("Wrong number of arguments!\n");
    exit(0);
  }

  freq		= atol(argv[1]);
  num_rays	= atoi(argv[2]);
  center_lat	= atof(argv[3]);
  delt_low	= atof(argv[4]);
  delt_high	= atof(argv[5]);
  final_TG	= atof(argv[6]);
  LAT_SPREAD	= atof(argv[7]);



  // Open a file for writing
  outputFileName = "newray.in";
  if((filePtr = fopen(outputFileName,"w"))==NULL) {
    printf("Problem opening %s \n",outputFileName);
    exit(0);
  } else {
    printf("Opened %s successfully for writing\n",outputFileName);
  }





  // Write all the cards:
  // --------------------



  // -------------------   CARD 1:  ---------------------------
  // Interactive: if n>0, enter NUMRES and NSUPPR interactively and
  // override the existing values.  If n<=0 use existing values.
  INTERA	= 0;	

  // Maximum number of raytracing continuations at the resonant
  // cone. If zero, then no continuation.
  NUMRES	= 10;	

  // values 1/0. If =1, then satellite trajectory in the range
  // lambda_min-3 <= lambda <= lambda_max+3 is printed.  
  NSUPPR	= 0;	

  // Special latitude: if SPELAT ~= 0, pogram tries to compute 
  // and print the output at the latitude lambda = 10 + eps
  SPELAT	= 0;	

  fprintf(filePtr,"%g %g %g %g\n", INTERA, NUMRES, NSUPPR, SPELAT);






  // ---------------- CARD 2:  Satellite info ------------------
  // Geoccentric distance of Satellite
  DISTRE = 22400;

  // Latitude of satellite (Don't confuse with LATITU of ray)
  LATITUD = 0;

  fprintf(filePtr,"%g %g\n", DISTRE, LATITUD);






  // ------------------  CARD 3: Separator------------------------
  // This just seems to be a kind of separator card
  fprintf(filePtr, "%g %g \n", -1.0, 0);






  // ------------------- CARD 4: Model input -------------------
  // Number of components in the plasma - number from 2-4, we 
  // have H, He, O, and el.'s
  NUM	= 4;

  // Is 0/1. If 0, integration size automatically adjusted to 
  // maintain accuracy specified by ABSB and RELB.  If 1, 
  // stepsize held constant
  KSKIP	= 0;

  // Is -1/1. If -1 proton whistler mode, if 1, el whistler mode
  MODE	= 1;

  // Spacing of points on line printer and (if KTAPE~=0) tape. 
  // If 0 - all points outputted on printer, otherwise, only 
  // every nth point is outputted. Critical points always printed 
  // i.e <0 critical points only, =0 all points, >0 critical 
  // and spaced ponits.
  KOUNT	= 2;

  // Indicates number of ducts in density model.  We usually have 
  // 1 or 2, because plasmapause is a duct and we had another one 
  // to make the steepness of electron dropoff accurate.
  KDUCTS = 4;

  // File number of results recorded on tape.  If =0, only
  // lineprinter, if =1, the first file is the program itself.
  KTAPE	= 1;

  // No description: maybe reference altitude of some sort?
  REFALT = 200;

  // Density reference point - range in L shell 
  DSRRNG	= 2;

  // Density reference point - latitude in degrees
  DSRLAT	= 0;

  // Density reference point - el/cc @ reference point
  DSRDENS	= 2500;

  fprintf(filePtr,"%g %g %g %g %g %g %g %g %g %g\n", 
	  NUM, KSKIP, MODE, KOUNT, KDUCTS, KTAPE, REFALT, 
	  DSRRNG, DSRLAT, DSRDENS );







  // ---------------------  CARD 5: Model Input ----------------
  // electron gyrofrequency in kHz, at surface of the earth, at
  // equator, for dipole magnetic field model
  EGFEQ	= 880;

  // Temperature used in the diffusive equilibrium density model
  THERM	= 1600;

  // Initial integration stepsize H = 50*HM/sqrt(f), except below
  // RDIV, where H is 1/16 of this value
  HM	= 10;

  // Absolute error in integration.  When any of the 5 variables
  // exceeds 14.2*ABSB and RELB, stepsize is halved.
  ABSB	= 0.001;

  // Relative error in integration - absolute error divided by the
  // value. When relative error exceeds 14.2*RELB and absolute 
  // error, stepsize is halved.
  RELB	= 1e-5;

  fprintf(filePtr, "%g %g %g %g %g\n", EGFEQ, THERM, HM, ABSB, RELB);







  // ----------- CARD 6: Diffusive Equilirium input  -----------
  // Geocentric base (in km) of DE density model
  RBASE	= 7370;

  // Electron density (cm^-3) at RBASE.  Final density here could
  // actually be different due to other things such as ducts, etc.
  ANE0	= 1;

  // relative concentration of H+ at RBASE, i.e. [H+]/[H+ + He+ + O+]
  ALPHA0_2= 0.08;

  // Relative concnetration of He+ at RBASE
  ALPHA0_3= 0.02;
  
  // Relative concentration of O+ at RBASE
  ALPHA0_4= 0.9;

  fprintf(filePtr,"%g %g %g %g %g\n", 
	  RBASE, ANE0, ALPHA0_2, ALPHA0_3, ALPHA0_4);







  // ---------------- CARD 7: Model Input ------------------------
  // Geocentric distance of the bottom of the ionosphere, where
  // density =0 (km)
  RZERO	= 6460;

  // Scale height at the bottom side of the lower ionosphere (km)
  SCBOT	= 140;
 
  // Geocentric distance below which raytracing stops (km)
  RSTOP	= 6470;		// 100km

  // Geocentric distance below which stepsize reset to 1/16 of 
  // normal starting value  (km)
  RDIV	= 6873;

  // Minimum allowed value of integration stepsize
  HMIN	= 0.001;

  fprintf( filePtr, "%g %g %g %g %g\n", 
	  RZERO, SCBOT, RSTOP, RDIV, HMIN);






  // ---------------  CARD 8: Plasmapause input ---------------
  // L-value of inner edge of the plasmapause, center of the 
  // knee is at LK + DDK
  LK	= 5.55;

  // Exponential component of density decrease beyond the knee,
  // R^{-EXPK}.  This is in addition to decrease due to DE model.
  EXPK 	= 0.13;

  // Halfwidth in L of the knee
  DDK	= 0.07;

  // Geocentric distance (km) where density outside the knee is 
  // equal to the density inside (???).  Modeil is not good 
  // below RCONSS outside of the knee.
  RCONSN = 1e-8;

  // Scale height of the radial density decrease above RCONSN 
  // above the knee.
  SCR	= 50;

  fprintf(filePtr, "%g %g %g %g %g\n", LK, EXPK, DDK, RCONSN, SCR);





  // ---------- CARD 9: Duct input (optional, repeatable) ----------
  // repeat for as many ducts as we want to include

  // L-value of the center of the duct
  L0		= 1.1;

  // Enhancement factor of the duct
  DEF		= 30.0;

  // Duct half-width in L.  If -ve, then sinusoidal perturbation, 
  // and DD is the period in L of the sinusoidal perturbation
  DD		= 4.0;

  // Geocentric radius (km) of the power end of the duct.  Below 
  // this height, the duct merges into the background plasma. If 
  // extends to earth, put RDUCLN = 6370 - In the NORTH
  RDUCLN	= 10370;

  // Radial scale height of duct, with which the LOWER end merges 
  // into the background plasma
  HDUCLN	= 20000;

  // Geocentric radius of the upper end of the duct.  If it extends 
  // to equator, set = 6370*L0 - In the NORTH
  RDUCUN	= 70081;

  // Radial scale height with which upper end of duct merges into 
  // the background plasma
  HDUCUN	= 20000;

  // Lower radius of duct in the SOUTH
  RDUCLS	= 10370;

  // Lower scale height of duct in the South
  HDUCLS	= 20000;

  // Upper radius of duct in the South
  RDUCUS	= 70081;

  // Upper scale height of duct in the South
  HDUCUS	= 20000;

  // Allows one-sided ducts.  If >0, all lower L values below L0 
  // stay at the enhancement factor DEF, like a secondary 
  // plasmapause. If =0, then is a normal 2 sides duct, if <0, 
  // all higher L values above L0 stay at the enhancement factor 
  // DEF, like a recovery region.
  SIDEDU	= 1;

  // DUCT 1
  fprintf(filePtr, "%g %g %g %g %g %g %g %g %g %g %g %g \n",
	  L0,  DEF,  DD,  RDUCLN,  HDUCLN,  
	  RDUCUN,  HDUCUN,  RDUCLS,  HDUCLS, 
	  RDUCUS,  HDUCUS,  SIDEDU );

  //DUCT 2
  fprintf(filePtr, "%g %g %g %g %g %g %g %g %g %g %g %g \n",
	  L0,  DEF,  DD,  RDUCLN,  HDUCLN,  
	  RDUCUN,  HDUCUN,  RDUCLS,  HDUCLS, 
	  RDUCUS,  HDUCUS,  SIDEDU );

  // DUCT 3
  fprintf(filePtr, "%g %g %g %g %g %g %g %g %g %g %g %g \n",
	  L0,  1.8,  1.8,  RDUCLN,  HDUCLN,  
	  RDUCUN,  HDUCUN,  RDUCLS,  HDUCLS, 
	  RDUCUS,  HDUCUS,  SIDEDU );





  // ------------------  CARD 10: Profile Input --------------------
  // Altitude step in Earth radii between points on radial profile, 
  // =0 for no profile
  PSTALT	= 0;

  // Geocentric distance (Earth radii) of lower limit of radial 
  // profile 
  PALT1	= 1;

  // Geocentric distance of upper limit of radial profile
  PALT2	= 7;

  // Constant latitude (in degrees) of radial profile
  PLATIT	= 0;

  // Latitude step (in degrees) between points on teh latitude
  // profile. =0 for no profile
  PSTLAT	= 1;

  // Lower latitude limit of profile (in degrees)
  PLAT1	= -60;
  
  // Upper latitude limit of profile (in degrees)
  PLAT2	= 60;

  // Constant altitude at which lat profile is to be computed
  PALTIT	= 1000;
  
  fprintf( filePtr, "%g %g %g %g %g %g %g %g\n", 
	   PSTALT, PALT1, PALT2, PLATIT, PSTLAT, PLAT1, PLAT2, PALTIT);
  



  // --------------  CARD 11: Ray input - repeatable --------------
  srand( (int) time(NULL) );
  for( i=0 ; i<=num_rays ; i++ ) {

    start_lat  = center_lat - LAT_SPREAD/2 ;
    lat_incr   = LAT_SPREAD/num_rays;
    
    delta_spread = delt_high - delt_low ;
    delta_incr   = delta_spread / num_rays ;
  
    // FKC is the frequency of the individual ray in kHz.
    FKC	= freq/1000.0;

    // Initial geocentric distance in km - launch altitude=1000km
    X0	= 7370 ;

    // Latitude of ray launch in degrees
    LATITU = start_lat + i*lat_incr ;

    // Wavenormal angle wrt local zenith
    DELT_offset = ( rand()/((double)RAND_MAX)- 0.5 )*RANGE;
    DELT	= delt_low + i*delta_incr;

    // Wavenormal angle wrt local B-field.  This value ignored 
    // unless DELT=360
    PSI	= 360;
    
    // Initial phase time in sec
    TINITI_P= 0;

    // Initial group delay time in sec
    TINITI_G= 0;

    // Final group time delay
    TGFINA = final_TG;

    fprintf(filePtr, "%g %g %g %.4g %g %g %g %g\n",  
	    FKC,X0,LATITU,DELT,PSI,TINITI_P,TINITI_G,TGFINA);
    
  }	
  
  
  
  
  
  // ---------------  CARD 12: Terminating card  -----------------
  fprintf(filePtr, "%g %d %d %d %d %d %d %d \n", 
	  -1.0, 0, 0, 0, 0, 0, 0, 0);
  
  
  
  
  
  // ----------------- CARD 13: Title of plot    ------------------
  fprintf(filePtr, "%s \n", "My Rays" );

  





  // -------- CARD 14: Plotting function specifications  --------
  // I'm not too sure what all of these do, but as far as I know 
  // they do not change.  Set to the same value as Dave's 
  // automatic scripts.
  KPLOT		= 2;
  XSCALE	= 1;
  YSCALE	= 1;
  FLBOT		= 1;
  FLTOP		= 8;
  LP_1		= 5;
  LP_2		= 4;
  LP_3		= 3;
  LP_4		= 2;
  LP_5		= 6;

  fprintf(filePtr, "%g %g %g %g %g %g %g %g %g %g\n", 
	  KPLOT, XSCALE, YSCALE, FLBOT, FLTOP, 
	  LP_1, LP_2, LP_3, LP_4, LP_5);







  // --------------- CARD 15: Terminating card -----------------
  fprintf(filePtr, "%s \n", "END");
 




  // ---------------------- Close newray.in -----------------
  fclose(filePtr);




  return 0;
}





