#include <stdio.h>
#include <math.h>
#include "time.h"
#include <unistd.h>
#include <stdlib.h> // Added, APS, 5/2015 




// Global variables:
// -----------------
int		GEOM;




// Definitions:
//-------------

#define		PI		3.14159265358979311599796346854418516159
#define		D2R		PI/180.0
#define		R2D		180.0/PI

#define		NUM_STEPS	15000
#define		T_STEP		0.0004

#define		Q_EL		1.602E-19
#define		M_EL		9.1E-31
#define		MU0		PI*4E-7		
#define		EPS0		8.854E-12
#define		C		2.997956376932163e+08
#define		Z0		377.0
#define		R_E		6370000.0
#define		H_MAGNETO	1E6
#define		H_IONO		1E5

#define		A		5E3	
#define		B		1E5
#define		H_E		5000.0
#define		I0		-10530.0

#define		P_DIST		0.0
#define		Q_DIST		2.0
#define		AN_CM_DIST	2E5
#define		V0_DIST		1.0
#define		M_RES		0








// Data structures:
//-----------------

// Store each ray in rayT
typedef struct {
  float	tg[NUM_STEPS];		// group time
  float	distre[NUM_STEPS];     	// geocentric dist in Re
  float	lat[NUM_STEPS];		// latitude
  float	delta[NUM_STEPS];      	// angle wrt vert
  float	t[NUM_STEPS];		// phase time	
  float	l_sh[NUM_STEPS];       	// local L-shell
  float	psi[NUM_STEPS];		// angle wrt B-field
  float	psiray[NUM_STEPS];     	// ray angle wrt B
  float	psires[NUM_STEPS];     	// res. cone angle wrt B
  float	mu[NUM_STEPS];		// refr. index
  float	dens[NUM_STEPS];       	// electron density /cm^3
  float	anH[NUM_STEPS];		// concentrations of: H+
  float	anHe[NUM_STEPS];       	//      ,,            He+
  float	anO[NUM_STEPS];		//      ,,            O+
  float	fH[NUM_STEPS];		// gyrofreq. in kHz
  float	stixP[NUM_STEPS];		
  float	stixR[NUM_STEPS];
  float	stixL[NUM_STEPS];	
  int	numSteps;		// num of steps in this ray
  long	f;			// ray frequency
  float pwr[NUM_STEPS];	       	// pwr/m^2 @ every point
} rayT;











// Prototypes:
//------------

FILE *readFile(char *fileName, FILE *filePtr, rayT *ray, long f);
void initRay(rayT *ray );
void initStep(rayT *ray, int i);
FILE *writeFile(long lower_freq, char *pref, char *suff);
void damping(rayT *ray, double AN, double v0_sq);
FILE *writeDamping(rayT *ray, FILE *outPtr);










/*
 * FUNCTION: main
 * --------------
 * This function reads in one rayfile of the format newray.dat and 
 * outputs a data file with the ray's Landau damping.
 *
 * Command-line arguments:
 * freq: frequency of ray in Hz
 * rayfile: eg. 'newray1.dat'
 * GEOM: include geometric factor, 0=NO, 1=YES
 *
 */

int main(int argc, char *argv[])
{
  FILE	*outPtr=NULL, *ptr=NULL;
  char sysCmd[128], *file, *freqPref="d", *freqSuff=""; 
  double AN;
  rayT ray;
  int i;
  long freq;






  // --------------  get arguments from command line ---------------
  if(argc != 4) 
    { printf("Wrong number of input arguments\n");  exit(0); }
  freq = atol(argv[1]);
  file = argv[2];
  GEOM = atoi(argv[3]);


  printf("\nCalculating damping for f: %d\n",freq);
  // ----------------------------------------------------------------




  // ----------------  Do some initialization  ----------------------
  AN	= AN_CM_DIST * pow( 10.0 , (12.0-(4.0*Q_DIST)) );
  outPtr = (FILE *)writeFile(freq, freqPref, freqSuff);
  // ----------------------------------------------------------------




  // ---------- Step through all rays in one rayfile ----------------
 
  while( 1 ) {
    ptr = (FILE *)readFile( file, ptr, &ray, freq);
    if(ptr==NULL) break;
    damping(&ray, AN, V0_DIST );
    outPtr = writeDamping(&ray, outPtr);    
  }
  // ------------------------------------------------------------------




  // ------------------ Close all files and return ---------------------
  if(ptr != NULL) fclose(ptr);
  fclose(outPtr);
  return 0;
  // --------------------------------------------------------------------

}









/*
 * FUNCTION: readFile
 * ------------------
 * This will take in the name of a file and return the next ray as 
 * a rayT* until the end of the file, when the file will be closed.
 * 
 */

FILE  *readFile( char *fileName, FILE *filePtr, rayT *ray, long f)
{
  int i=0;
  

  //  Open the file for reading  //
  if(filePtr == NULL) {
    if( (filePtr = fopen( fileName ,"r")) == NULL ) {
      printf("Hey buddy, I can't open the file!\n");
      exit(1);
    } else {
      printf("\nFile %s\nOpened successfully! \n", fileName);
    }
    // Place the pointer after the header for reading in data
    fseek(filePtr, 146*sizeof(float), SEEK_SET);

  }


  // Read in each ray //
  initRay( ray );
  while( 1 ) { 
    if(fscanf(filePtr, 
	      "%e%e%e%e%e%e%e%e%e%e%e%e%e%e%e%e%e%e", 
	      &(ray->tg[i]), &(ray->distre[i]), &(ray->lat[i]), 
	      &(ray->delta[i]), &(ray->t[i]), &(ray->l_sh[i]),
	      &(ray->psi[i]), &(ray->psiray[i]), &(ray->psires[i]), 
	      &(ray->mu[i]), &(ray->dens[i]), &(ray->anH[i]), 
	      &(ray->anHe[i]), &(ray->anO[i]), &(ray->fH[i]), 
	      &(ray->stixP[i]), &(ray->stixR[i]), 
	      &(ray->stixL[i]) ) == EOF ) {
      ray->numSteps=0;
      return NULL;
    }
    i++;
    
    if( ray->tg[(i-1)]>=99998.0 ) {
      initStep( ray, (i-1) );
      break; 
    }

    if(i==NUM_STEPS) {
      i--;	// back up a step, record next step over this one until
		// eventually it is rewritten by last line & initialised
      printf("\nnumber of NUM_STEPS reached\a\n");
    }

  }
  ray->numSteps = i-1;

  // Set the ray frequency 
  ray->f = f;

  printf("\nNow reading, f: %d, lat: %g\n",ray->f, ray->lat[0]);

  return filePtr; 
}

  




/* 
 * FUNCTION initRay
 * ----------------
 * Takes in a ray structure and initialises all the variables
 *
 */

void initRay(rayT *ray )
{ 
  int i;

  for( i=0 ; i<NUM_STEPS ; i++)
    initStep(ray,i);
      
  ray->numSteps = 0;
  ray->f = 0;
}


void initStep(rayT *ray, int i)
{
   ray->tg[i]=0.0;	     	// group time
   ray->distre[i]=0.0;     	// geocentric dist in Re
   ray->lat[i]=0.0;		// latitude
   ray->delta[i]=0.0;      	// angle wrt vert
   ray->t[i]=0.0;      		// phase time	
   ray->l_sh[i]=0.0;       	// local L-shell
   ray->psi[i]=0.0;		// angle wrt B-field
   ray->psiray[i]=0.0;     	// ray angle wrt B
   ray->psires[i]=0.0;     	// res. cone angle wrt B
   ray->mu[i]=0.0;		// refr. index
   ray->dens[i]=0.0;       	// electron density /cm^3
   ray->anH[i]=0.0;		// concentrations of: H+
   ray->anHe[i]=0.0;       	//      ,,            He+
   ray->anO[i]=0.0;		//      ,,            O+
   ray->fH[i]=0.0;		// gyrofreq. in kHz
   ray->stixP[i]=0.0;		
   ray->stixR[i]=0.0;
   ray->stixL[i]=0.0;
   ray->pwr[i]=0.0;
}









/*
 * FUNCTION: writeFile
 * -------------------
 * This function just gets an appropriate name for the output file 
 * and opens it, checking that there is no problem.  Returns a 
 * pointer to the open file.
 *
 */

FILE *writeFile(long freq, char *pref, char *suff)
{
  char outFileName[64];
  FILE *outPtr;

  // make some kind of filename convention, return pointer
  sprintf( outFileName, "%s%d%s.dat", pref, freq, suff);  

  if((outPtr=fopen(outFileName,"w"))== NULL) {
    printf("Problem opening output file!\n");
    exit(1);
  } else {
    printf("File %s \nopened successfully for writing!\n", outFileName);
  }
  return outPtr;
}









/*
 * FUNCTION: damping
 * -----------------
 * This function calculates the damping of the wave at every step 
 * along its path, using Brinca's 1972 formula.
 *
 *
 *
 */

void damping(rayT *ray, double AN, double v0_sq)
{
  double p,r,l,s,d,a,b,n,theta,w,dens,Omega,v_perp, v_step;
  double sin_th, cos_th, sin_th_sq, cos_th_sq,n_sq,v_perp_sq;
  double k,k_para,k_perp,wp,const1,Vm,c1,bessel_arg, bessel_term;
  double distrib_term, velocity_term,Vm_sq,v_tot_sq,c2,term1;
  double twoPPlusOne,g2,c3,term2,term3, chi1, ds, k_im, damp=0.0;
  int i; 
  double n_steps, geom_fact, DELTA=1.0;
  time_t start,end;

  start = time(NULL);
  twoPPlusOne = 2*P_DIST + 1 ;
  n_steps = 500.0;
  //v_step = C/n_steps; <- change back to this!
  v_step = (3.0E8-10.0)/n_steps;

  w = (ray->f)*2*PI;	// already in Hz
  ray->pwr[0] = 1;

  for(i=1 ; i<ray->numSteps ; i++) {
    
    // calc geometric factor  
    if(GEOM==0) {
      geom_fact = 1.0;
    } else {
    geom_fact = ray->distre[0] * cos(ray->lat[0]*D2R) / 
      ( (ray->distre[i]) * cos(ray->lat[i]*D2R) ) ;
    }

    if(DELTA==0 /*|| (ray->pwr[i-1] < 1e-7*MAX_POWER)*/  ) {
      if(DELTA==0) 
	ray->pwr[i] = geom_fact*ray->pwr[0]; 

      // otherwise leave the power value as 0

    } else  {
      
      if( fabs(ray->psi[i])<=90 ) theta = fabs(ray->psi[i])*D2R;
      
      if(fabs(ray->psi[i]) > 90.0 && fabs(ray->psi[i]) <= 180.0 ) {
	theta = (180.0 - fabs(ray->psi[i]) )*D2R ; }
      
      if(fabs(ray->psi[i]) > 180.0 && fabs(ray->psi[i]) <= 270.0 ) {
	theta = (fabs(ray->psi[i]) - 180.0)*D2R ; }
      
      
      
    p = (double)ray->stixP[i];
    r = (double)ray->stixR[i];
    l = (double)ray->stixL[i];    
    s = (r+l)/2 ;
    d = (r-l)/2;
  
    sin_th = sin(theta);
    cos_th = cos(theta);
    sin_th_sq = pow(sin_th,2);
    cos_th_sq = pow(cos_th,2);
    
    a = s*sin_th_sq + p*cos_th_sq;
    b = r*l*sin_th_sq + p*s*(1+cos_th_sq);

    dens = ray->dens[i] * 1e6;	//convert from el/cm^3 -> el/m^3
    DELTA = 1/dens;
    Omega = ray->fH[i]*2000.0*PI;	//convert from f[kHz]->rad/s 
    n = ray->mu[i];
    n_sq = pow(n,2);
    
    k = w/C*n;
    k_para = k*cos_th ;
    k_perp = k*sin_th ;
    
    wp = Q_EL*sqrt( dens / (EPS0*M_EL) );  
    const1  = DELTA / ( 4*n*(2*a*n_sq - b) );
    Vm = (w + M_RES*Omega) / k_para;  
    Vm_sq = pow(Vm,2);
    c1 = (4* AN * pow(PI,2) * pow(wp,2) ) / ( w * k_para ) ;
    
    chi1 = 0.0;

    for(v_perp=10; v_perp<C ; v_perp+=v_step ) {
     
      v_perp_sq = pow(v_perp,2);
      bessel_arg = k_perp*v_perp / Omega ;
      bessel_term =	(r-n_sq) * jn((M_RES-1),bessel_arg) + 
			(l-n_sq) * jn((M_RES+1),bessel_arg) ;
      v_tot_sq = v_perp_sq + Vm_sq + v0_sq ;
      distrib_term = P_DIST*v_tot_sq  -  v_perp_sq*Q_DIST; 
      velocity_term = pow( v_tot_sq , (Q_DIST+1) );


      //TERM1:
      c2 = ( n_sq*sin_th_sq - p ) / ( 2*(s-n_sq) ) ;
      term1 =	c2 * c1 * 
	pow( bessel_term , 2 ) *
	pow( v_perp , twoPPlusOne ) / (w*velocity_term) *
	( w*distrib_term - k_para*Vm*P_DIST*v_tot_sq ) ;


      //TERM2:
      g2 = (1+ M_RES*Omega/w)*(Q_DIST*Vm*v_perp) + 
	(M_RES*Omega*Vm/(w*v_perp)) * distrib_term ;
      term2 = 2*c1*((s-n_sq*cos_th_sq)*(s-n_sq) - pow(d,2) ) *
	pow( jn(M_RES,bessel_arg) , 2 ) * Vm * pow(v_perp_sq,P_DIST)/
	velocity_term * g2 ;


      //TERM3:
      c3 = 2*n_sq*sin_th*cos_th ;
      term3 = c3 * c1 * pow(v_perp,twoPPlusOne)*jn(M_RES,bessel_arg) /
	velocity_term * bessel_term * g2 ;


      chi1 += (term1 + term2 + term3) ;
    }

    ds = sqrt( pow( (ray->distre[i]*cos((ray->lat[i])*D2R) - 
	       ray->distre[i-1]* cos((ray->lat[i-1])*D2R)), 2 ) +
	       pow( (ray->distre[i]*sin((ray->lat[i])*D2R) - 
	       ray->distre[i-1]*sin((ray->lat[i-1])*D2R)), 2))*R_E;

    k_im = chi1*const1*v_step*w/C ;
    damp += k_im*ds;

    // use 2*damp because it's POWER damping, not Bw or E
    ray->pwr[i] = geom_fact * exp( 2.0 * damp ) * ray->pwr[0];
    // printf("%e\n",exp(2*damp) );
    
    }  // else ... (if DELTA !=0) 
  } // for i<NUM_STEPS

  end = time(NULL);
  printf("\ndamping of f: %d, lat: %g took %g sec\n", 
	 ray->f, ray->lat[0], difftime(end,start));  
}






/*
 * FUNCTION: writeDamping
 * ----------------------
 * This function simply writes out the relative power damping of the current
 * ray as one long column, separated by the number 99999 to be consistent 
 * file separators in newray.dat
 *
 */
FILE *writeDamping(rayT *ray, FILE *outPtr)
{
  int i;

  for(i=0; i<ray->numSteps; i++) {
    fprintf(outPtr, "%g\n", ray->pwr[i]);
  }
  fprintf(outPtr, "99999\n");
  return outPtr;
}
