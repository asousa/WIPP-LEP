/*---------------------------------------------------------------
* 8/10/06:
  Already made some changes.  Unfortunately, just started documenting
  them now.

  1) Instead of using large number of frequencies (necessary for broad-band
  whistlers), reading in 2 frequency components from command line.  The lower is
  the transmitter fkc, and the upper is a small bandwidth added to the signal.

  Other than this modification, the code is basically unchanged.
------------------------------------------------------------------*/



#include <stdio.h>
#include <math.h>
#include "time.h"
#include <unistd.h>
#include <stdlib.h>




// DEFINITIONS
// -----------


#define		RES_DT		0.02
#define		RES_FINT	5.0
#define		C		2.997956376932163e+08
#define		R_E		6370000.0
#define		H_IONO		1E5
#define		PI		3.141592653589793115997963468544185
#define		D2R		PI/180.0
#define		M_EL		9.1E-31
#define		E_EL		5.105396765648739E5

#define		DE_EXP		0.003
#define		E_EXP_BOT	1.477
#define		E_EXP_TOP	7.477
#define		NUM_E		2200
#define		SQUARE		1


/*
 * FUNCTION: main
 * --------------
 * This is an upgrade of calcFluxE.  It now reads AE8 
 * (L-dependent) flux values converts into differential 
 * flux, and assumes "sin" shaped loss-cone.
 *
 */

int main(int argc, char *argv[])
{
  printf("Floats are %s bytes on this machine\n",sizeof(float));

  return 0;
}