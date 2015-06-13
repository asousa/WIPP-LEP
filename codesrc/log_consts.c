/* log_consts.c
 *
 *  Prints a matlab script that loads all constants used in the simulation
 *  (Useful for plotting data, making sure we're using correct values)
 */


#include <stdio.h>
#include <math.h>
#include "time.h"
#include <unistd.h>
#include <stdlib.h>
#include "consts.h"




int main(int argc, char *argv[])
{
  int numFreqs, i;
  FILE *f = fopen("sim_consts.m", "w");
  if (f == NULL)
  {
      printf("Error opening file!\n");
      exit(1);
  }


fprintf(f, "%% -------- Log of simulation constants ---------\n");
time_t rawtime;
struct tm * timeinfo;

time ( &rawtime );
timeinfo = localtime ( &rawtime );
fprintf (f,"%% Ran at: %s", asctime (timeinfo) );


fprintf(f,"RAYTRACE_TIME=%g;\n",RAYTRACE_TIME);
fprintf(f,"NUM_STEPS=%d;\n",NUM_STEPS);
fprintf(f,"T_STEP=%g;\n",T_STEP);


fprintf(f,"Q_EL=%g;\n",Q_EL);
fprintf(f,"M_EL=%g;\n",M_EL);
fprintf(f,"E_EL=%g;\n",E_EL);
fprintf(f,"MU0=%g;\n",MU0);
fprintf(f,"EPS0=%g;\n",EPS0);
fprintf(f,"C=%g;\n",C);
fprintf(f,"Z0=%g;\n",Z0);
fprintf(f,"R_E=%g;\n",R_E);
fprintf(f,"H_MAGNETO=%g;\n",H_MAGNETO);
fprintf(f,"H_IONO=%g;\n",H_IONO);

fprintf(f,"A=%g;\n",A);
fprintf(f,"B=%g;\n",B);
fprintf(f,"H_E=%g;\n",H_E);
fprintf(f,"I=%g;\n",I0);

fprintf(f,"P_DIST=%g;\n",P_DIST);
fprintf(f,"Q_DIST=%g;\n",Q_DIST);
fprintf(f,"AN_CM_DIST=%g;\n",AN_CM_DIST);
fprintf(f,"V0_DIST=%g;\n",V0_DIST);
fprintf(f,"M_RES=%g;\n",M_RES);

fprintf(f,"RES_DT=%g;\n",RES_DT);
fprintf(f,"RES_FINT=%g;\n",RES_FINT);
fprintf(f,"M_EL=%g;\n",M_EL);
fprintf(f,"E_EL=%e;\n",E_EL);

// Energy scaling (change these inputs):
//fprintf(f,"   E_EXP_BOT 
fprintf(f,"E_MIN=%e;\n",E_MIN);
fprintf(f,"E_MAX=%e;\n",E_MAX);
fprintf(f,"NUM_E=%d;\n",NUM_E);
fprintf(f,"SQUARE=%d;\n",SQUARE);

// Calculated params
fprintf(f,"E_EXP_BOT=%g;\n",E_EXP_BOT);
fprintf(f,"E_EXP_TOP=%g;\n",E_EXP_TOP);
fprintf(f,"DE_EXP=%g;\n",DE_EXP);

fprintf(f,"EALimS=%g;\n",EALimS);
fprintf(f,"EALimN=%g;\n",EALimN);
fprintf(f,"EAIncr=%g;\n",EAIncr);
fprintf(f,"dL0=%g;\n",dL0);
fprintf(f,"DF=%g;\n",DF);
fprintf(f,"DT=%g;\n",DT);
fprintf(f,"NUMLATS=%f;\n",NUMLATS);

fprintf(f,"EA_SPLIT=%d;\n",EA_SPLIT);
fprintf(f,"MULT=%f;\n",MULT);

fprintf(f,"sim_freqs=[ ");//%d,\n",freqs[0]);


numFreqs = (sizeof freqs)/sizeof(long) - 1;

for(i=1; i<=numFreqs; i++) {
  fprintf(f," %d,",freqs[i]);
}
fprintf(f," ];\n");//,freqs[numFreqs]);


// fprintf(f, "Integer: %d, float: %f\n", i, py);
fclose(f); 
}


