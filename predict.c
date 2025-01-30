#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"


/*! \file predict.c
 *  \brief drift particles by a small time interval
 *
 *  This function contains code to implement a drift operation on all the
 *  particles, which represents one part of the leapfrog integration scheme.
 */


/*! This function drifts all particles from the current time to the future:
 *  time0 - > time1
 *
 *  If there is no explicit tree construction in the following timestep, the
 *  tree nodes are also drifted and updated accordingly. Note: For periodic
 *  boundary conditions, the mapping of coordinates onto the interval
 *  [0,All.BoxSize] is only done before the domain decomposition, or for
 *  outputs to snapshot files.  This simplifies dynamic tree updates, and
 *  allows the domain decomposition to be carried out only every once in a
 *  while.
 */
void move_particles(int time0, int time1)
{
  int i, j;
  double dt_drift, dt_gravkick, dt_hydrokick, dt_entr;
  double t0, t1;


  t0 = second();

  if(All.ComovingIntegrationOn)
    {
      dt_drift = get_drift_factor(time0, time1);
      dt_gravkick = get_gravkick_factor(time0, time1);
      dt_hydrokick = get_hydrokick_factor(time0, time1);
    }
  else
    {
      dt_drift = dt_gravkick = dt_hydrokick = (time1 - time0) * All.Timebase_interval;
    }

  // We use the force tree to search for neighborhood and do scattering
  if(dt_drift>0.0 && !(TreeReconstructFlag)){
  //if(dt_drift>0.0){ // results similar, but disturbance exists occasionally...
     // In the beginning of scattering, will clear all nscatt
     sidm(dt_drift);
  }

  //if(ThisTask==0) printf("\n XXXXXX dt_drift XXXXXXXXXX %g\n",dt_drift);
  for(i = 0; i < NumPart; i++)
    {
      
      //if(abs(P[i].Pos[0])>3000||abs(P[i].Pos[1])>3000||abs(P[i].Pos[2])>3000){
      //   double ri = sqrt(P[i].Pos[0]*P[i].Pos[0]+P[i].Pos[1]*P[i].Pos[1]+P[i].Pos[2]*P[i].Pos[2]);
      //   if(ri>4500){
      //      P[i].Mass=0;
      //      P[i].Type=3;
      //      P[i].Pos[0]*=4500./ri;
      //      P[i].Pos[1]*=4500./ri;
      //      P[i].Pos[2]*=4500./ri;
      //      //P[i].Pos[0]=0.;
      //      //P[i].Pos[1]=0.;
      //      //P[i].Pos[2]=0.;
      //      //P[i].Vel[0]*=0.00001;
      //      //P[i].Vel[1]*=0.00001;
      //      //P[i].Vel[2]*=0.00001;
      //      P[i].Vel[0]=0.;
      //      P[i].Vel[1]=0.;
      //      P[i].Vel[2]=0.;
      //   }
      //}
      
      
      for(j = 0; j < 3; j++)
	P[i].Pos[j] += P[i].Vel[j] * dt_drift;

/*
      double rin=sqrt(pow(P[i].Pos[0],2)+pow(P[i].Pos[1],2)+pow(P[i].Pos[2],2));
      if(rin<0.8){
         FILE * fp;
         fp = fopen("./output/innerJ.txt","a");
         double Ji = pow(P[i].Pos[1]*P[i].Vel[2]-P[i].Pos[2]*P[i].Vel[1],2)+pow(P[i].Pos[0]*P[i].Vel[2]-P[i].Pos[2]*P[i].Vel[0],2)+pow(P[i].Pos[0]*P[i].Vel[1]-P[i].Pos[1]*P[i].Vel[0],2);
         Ji=sqrt(Ji);
         fprintf (fp, "%f %f\n",rin,Ji);
         fclose (fp);
      }
*/
      // print out the orbit of a particle
      //if(P[i].ID==1190995456){
      //if(P[i].ID==1201233920){
      //if(P[i].ID==1234){
      //     FILE * fp;
      //     fp = fopen("./distfunc/part1234.txt","a");
      //     fprintf (fp, "%f %f %f %f %f %f\n",P[i].Pos[0],P[i].Pos[1],P[i].Pos[2],P[i].Vel[0],P[i].Vel[1],P[i].Vel[2]);
      //     fclose (fp);
      //}
      //printf("@@@@ particle ID:%llu\n",P[i].ID);

      if(P[i].Type == 0)
	{
#ifdef PMGRID
	  for(j = 0; j < 3; j++)
	    SphP[i].VelPred[j] +=
	      (P[i].GravAccel[j] + P[i].GravPM[j]) * dt_gravkick + SphP[i].HydroAccel[j] * dt_hydrokick;
#else
	  for(j = 0; j < 3; j++)
	    SphP[i].VelPred[j] += P[i].GravAccel[j] * dt_gravkick + SphP[i].HydroAccel[j] * dt_hydrokick;
#endif
	  SphP[i].Density *= exp(-SphP[i].DivVel * dt_drift);
	  SphP[i].Hsml *= exp(0.333333333333 * SphP[i].DivVel * dt_drift);

	  if(SphP[i].Hsml < All.MinGasHsml)
	    SphP[i].Hsml = All.MinGasHsml;

	  dt_entr = (time1 - (P[i].Ti_begstep + P[i].Ti_endstep) / 2) * All.Timebase_interval;

	  SphP[i].Pressure = (SphP[i].Entropy + SphP[i].DtEntropy * dt_entr) * pow(SphP[i].Density, GAMMA);
	}
    }

  /* if domain-decomp and tree are not going to be reconstructed, update dynamically.  */
  if(All.NumForcesSinceLastDomainDecomp < All.TotNumPart * All.TreeDomainUpdateFrequency)
    {
      for(i = 0; i < Numnodestree; i++)
	for(j = 0; j < 3; j++)
	  Nodes[All.MaxPart + i].u.d.s[j] += Extnodes[All.MaxPart + i].vs[j] * dt_drift;

      force_update_len();

      force_update_pseudoparticles();
    }

  t1 = second();

  All.CPU_Predict += timediff(t0, t1);
}



/*! This function makes sure that all particle coordinates (Pos) are
 *  periodically mapped onto the interval [0, BoxSize].  After this function
 *  has been called, a new domain decomposition should be done, which will
 *  also force a new tree construction.
 */
#ifdef PERIODIC
void do_box_wrapping(void)
{
  int i, j;
  double boxsize[3];

  for(j = 0; j < 3; j++)
    boxsize[j] = All.BoxSize;

#ifdef LONG_X
  boxsize[0] *= LONG_X;
#endif
#ifdef LONG_Y
  boxsize[1] *= LONG_Y;
#endif
#ifdef LONG_Z
  boxsize[2] *= LONG_Z;
#endif

  for(i = 0; i < NumPart; i++)
    for(j = 0; j < 3; j++)
      {
	while(P[i].Pos[j] < 0)
	  P[i].Pos[j] += boxsize[j];

	while(P[i].Pos[j] >= boxsize[j])
	  P[i].Pos[j] -= boxsize[j];
      }
}
#endif
