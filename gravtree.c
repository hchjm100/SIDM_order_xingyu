#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"

// Pos: comoving coordinates in units h^-1 kpc
// physical pos: Pos*a/h = Pos*All.Time/All.HubbleParam
double MOverR3(FLOAT r) // comoving r in units kpc/h
{
   double  z = 1.0/All.Time -1.0;
   double  rhos= 0.00008018508808657744 + 0.00001173220974381105*z + 0.00005788062708937496*pow(z,2) - 0.00001091975902862681*pow(z,3) + 4.322228318940011e-6*pow(z,4);
   rhos = rhos/pow(All.HubbleParam,2);
   double  rs = 101.07117877864096 - 8.071187120090881*z - 40.28285233649378*pow(z,2) + 22.552770822134857*pow(z,3) - 4.007626606926756*pow(z,4);

   rs = rs*All.HubbleParam; 
   double  rphys = All.Time*r;
   //----------------------
   double  rhoh=0.104684*pow((1. + z),3.63)*(10.0076 - 2.78014*z + 3.65695*pow(z,2) - 4.04992*pow(z,3) + 1.02619*pow(z,4));
   rhoh = rhoh/pow(All.HubbleParam,2);
   double  rh=1.14986/pow((1. + z),1.21);
   rh = rh*All.HubbleParam;
   //----------------------
   double h = All.ForceSoftening[1]; //SofteningHalo*2.8
   double h_inv,h3_inv,u;
   h_inv = 1.0 / h;
   h3_inv = h_inv * h_inv * h_inv;
   double valHalo=(4*3.14159265358979*rhos*pow(rs,3)*(-(rphys/(rphys + rs)) - log(rs/(rphys + rs)) ));
   double valHern=2*3.14159265358979*rhoh*pow(rphys,2)*pow(rh,3)/(pow(rphys+rh,2));

   //printf("XXXXX: r %g, rho %g, rh %g, M %g \n",rphys/0.7,rhoh*0.7*0.7,rh/0.7,valHern/0.7);

   if(r>h){
      valHalo=valHalo/pow(r,3);
      valHern=valHern/pow(r,3);
   }
   else{
          u = r * h_inv;
          if(u < 0.5){
            valHalo = valHalo * h3_inv * (10.666666666667 + u * u * (32.0 * u - 38.4));
            valHern = valHern * h3_inv * (10.666666666667 + u * u * (32.0 * u - 38.4));
          }
          else{
            valHalo =
              valHalo * h3_inv * (21.333333333333 - 48.0 * u +
                               38.4 * u * u - 10.666666666667 * u * u * u - 0.066666666667 / (u * u * u));
            valHern =
              valHern * h3_inv * (21.333333333333 - 48.0 * u +
                               38.4 * u * u - 10.666666666667 * u * u * u - 0.066666666667 / (u * u * u));
          }
   }
   return valHalo+valHern; 

}


double FR_disk_ZXY(FLOAT M_d, FLOAT a_d, FLOAT b_d, FLOAT R, FLOAT z){

  return -All.G * M_d * R / pow(R*R + pow( a_d + sqrt(b_d*b_d + z*z), 2), 1.5);

}


double Fz_disk_ZXY(FLOAT M_d, FLOAT a_d, FLOAT b_d, FLOAT R, FLOAT z){

  return -All.G * M_d * z * ( a_d + sqrt(b_d*b_d + z*z) ) / sqrt(b_d*b_d + z*z) / pow(R*R + pow( a_d + sqrt(b_d*b_d + z*z), 2), 1.5);

}


double Fr_Hern_ZXY(FLOAT M_h, FLOAT r_h, FLOAT r){

  return -All.G * M_h / pow( r_h + r , 2); 

}

// 2023.04.20 ZXY: Actually all the forces do not diverge as r -> 0, therefore no need of softening ......
//                 I just take the limit of r -> 0 ......

double Fr_NFW_ZXY(FLOAT rhos, FLOAT rs, FLOAT r){

  if(r < 3*All.SofteningHalo){
    return -All.G * 2*3.141592653 * rs * rhos;
  }
  else{
    return -All.G * 4*3.141592653 * pow(rs,3) * rhos * ( -r/(r+rs) + log(1+r/rs) ) / pow(r,2);
  }

}



/*! \file gravtree.c 
 *  \brief main driver routines for gravitational (short-range) force computation
 *
 *  This file contains the code for the gravitational force computation by
 *  means of the tree algorithm. To this end, a tree force is computed for
 *  all active local particles, and particles are exported to other
 *  processors if needed, where they can receive additional force
 *  contributions. If the TreePM algorithm is enabled, the force computed
 *  will only be the short-range part.
 */

/*! This function computes the gravitational forces for all active
 *  particles.  If needed, a new tree is constructed, otherwise the
 *  dynamically updated tree is used.  Particles are only exported to other
 *  processors when really needed, thereby allowing a good use of the
 *  communication buffer.
 */
void gravity_tree(void)
{
  long long ntot;
  int numnodes, nexportsum = 0;
  int i, j, iter = 0;
  int *numnodeslist, maxnumnodes, nexport, *numlist, *nrecv, *ndonelist;
  double tstart, tend, timetree = 0, timecommsumm = 0, timeimbalance = 0, sumimbalance;
  double ewaldcount;
  double costtotal, ewaldtot, *costtreelist, *ewaldlist;
  double maxt, sumt, *timetreelist, *timecommlist;
  double fac, plb, plb_max, sumcomm;

#ifndef NOGRAVITY
  int *noffset, *nbuffer, *nsend, *nsend_local;
  long long ntotleft;
  int ndone, maxfill, ngrp;
  int k, place;
  int level, sendTask, recvTask;
  double ax, ay, az;
  MPI_Status status;
#endif

  /* set new softening lengths */
  if(All.ComovingIntegrationOn)
    set_softenings();


  /* contruct tree if needed */
  tstart = second();
  if(TreeReconstructFlag)
    {
      if(ThisTask == 0)
	printf("Tree construction.\n");

      force_treebuild(NumPart);

      TreeReconstructFlag = 0;

      if(ThisTask == 0)
	printf("Tree construction done.\n");
    }
  tend = second();
  All.CPU_TreeConstruction += timediff(tstart, tend);

  costtotal = ewaldcount = 0;

  /* Note: 'NumForceUpdate' has already been determined in find_next_sync_point_and_drift() */
  numlist = malloc(NTask * sizeof(int) * NTask);
  MPI_Allgather(&NumForceUpdate, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);
  for(i = 0, ntot = 0; i < NTask; i++)
    ntot += numlist[i];
  free(numlist);


#ifndef NOGRAVITY
  if(ThisTask == 0)
    printf("Begin tree force.\n");


#ifdef SELECTIVE_NO_GRAVITY
  for(i = 0; i < NumPart; i++)
    if(((1 << P[i].Type) & (SELECTIVE_NO_GRAVITY)))
      P[i].Ti_endstep = -P[i].Ti_endstep - 1;
#endif


  noffset = malloc(sizeof(int) * NTask);	/* offsets of bunches in common list */
  nbuffer = malloc(sizeof(int) * NTask);
  nsend_local = malloc(sizeof(int) * NTask);
  nsend = malloc(sizeof(int) * NTask * NTask);
  ndonelist = malloc(sizeof(int) * NTask);

  i = 0;			/* beginn with this index */
  ntotleft = ntot;		/* particles left for all tasks together */

  while(ntotleft > 0)
    {
      iter++;

      for(j = 0; j < NTask; j++)
	nsend_local[j] = 0;

      /* do local particles and prepare export list */
      tstart = second();
      for(nexport = 0, ndone = 0; i < NumPart && nexport < All.BunchSizeForce - NTask; i++)
	if(P[i].Ti_endstep == All.Ti_Current)
	  {
	    ndone++;

	    for(j = 0; j < NTask; j++)
	      Exportflag[j] = 0;
#ifndef PMGRID
	    costtotal += force_treeevaluate(i, 0, &ewaldcount);
#else
	    costtotal += force_treeevaluate_shortrange(i, 0);
#endif
	    for(j = 0; j < NTask; j++)
	      {
		if(Exportflag[j])
		  {
		    for(k = 0; k < 3; k++)
		      GravDataGet[nexport].u.Pos[k] = P[i].Pos[k];
#ifdef UNEQUALSOFTENINGS
		    GravDataGet[nexport].Type = P[i].Type;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
		    if(P[i].Type == 0)
		      GravDataGet[nexport].Soft = SphP[i].Hsml;
#endif
#endif
		    GravDataGet[nexport].w.OldAcc = P[i].OldAcc;
		    GravDataIndexTable[nexport].Task = j;
		    GravDataIndexTable[nexport].Index = i;
		    GravDataIndexTable[nexport].SortIndex = nexport;
		    nexport++;
		    nexportsum++;
		    nsend_local[j]++;
		  }
	      }
	  }
      tend = second();
      timetree += timediff(tstart, tend);

      qsort(GravDataIndexTable, nexport, sizeof(struct gravdata_index), grav_tree_compare_key);

      for(j = 0; j < nexport; j++)
	GravDataIn[j] = GravDataGet[GravDataIndexTable[j].SortIndex];

      for(j = 1, noffset[0] = 0; j < NTask; j++)
	noffset[j] = noffset[j - 1] + nsend_local[j - 1];

      tstart = second();

      MPI_Allgather(nsend_local, NTask, MPI_INT, nsend, NTask, MPI_INT, MPI_COMM_WORLD);

      tend = second();
      timeimbalance += timediff(tstart, tend);

      /* now do the particles that need to be exported */

      for(level = 1; level < (1 << PTask); level++)
	{
	  tstart = second();
	  for(j = 0; j < NTask; j++)
	    nbuffer[j] = 0;
	  for(ngrp = level; ngrp < (1 << PTask); ngrp++)
	    {
	      maxfill = 0;
	      for(j = 0; j < NTask; j++)
		{
		  if((j ^ ngrp) < NTask)
		    if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
		      maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
		}
	      if(maxfill >= All.BunchSizeForce)
		break;

	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;

	      if(recvTask < NTask)
		{
		  if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
		      /* get the particles */
		      MPI_Sendrecv(&GravDataIn[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(struct gravdata_in), MPI_BYTE,
				   recvTask, TAG_GRAV_A,
				   &GravDataGet[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct gravdata_in), MPI_BYTE,
				   recvTask, TAG_GRAV_A, MPI_COMM_WORLD, &status);
		    }
		}

	      for(j = 0; j < NTask; j++)
		if((j ^ ngrp) < NTask)
		  nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	    }
	  tend = second();
	  timecommsumm += timediff(tstart, tend);


	  tstart = second();
	  for(j = 0; j < nbuffer[ThisTask]; j++)
	    {
#ifndef PMGRID
	      costtotal += force_treeevaluate(j, 1, &ewaldcount);
#else
	      costtotal += force_treeevaluate_shortrange(j, 1);
              //if(j==0 && ThisTask==0) printf("Grav nbuffer[ThisTask=%d]=%d,level=%d\n",ThisTask,nbuffer[ThisTask],level);

#endif
	    }
	  tend = second();
	  timetree += timediff(tstart, tend);

	  tstart = second();
	  MPI_Barrier(MPI_COMM_WORLD);
	  tend = second();
	  timeimbalance += timediff(tstart, tend);

	  /* get the result */
	  tstart = second();
	  for(j = 0; j < NTask; j++)
	    nbuffer[j] = 0;
	  for(ngrp = level; ngrp < (1 << PTask); ngrp++)
	    {
	      maxfill = 0;
	      for(j = 0; j < NTask; j++)
		{
		  if((j ^ ngrp) < NTask)
		    if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
		      maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
		}
	      if(maxfill >= All.BunchSizeForce)
		break;

	      sendTask = ThisTask;
	      recvTask = ThisTask ^ ngrp;
	      if(recvTask < NTask)
		{
		  if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
		      /* send the results */
		      MPI_Sendrecv(&GravDataResult[nbuffer[ThisTask]],
				   nsend[recvTask * NTask + ThisTask] * sizeof(struct gravdata_in),
				   MPI_BYTE, recvTask, TAG_GRAV_B,
				   &GravDataOut[noffset[recvTask]],
				   nsend_local[recvTask] * sizeof(struct gravdata_in),
				   MPI_BYTE, recvTask, TAG_GRAV_B, MPI_COMM_WORLD, &status);

		      /* add the result to the particles */
		      for(j = 0; j < nsend_local[recvTask]; j++)
			{
			  place = GravDataIndexTable[noffset[recvTask] + j].Index;

			  for(k = 0; k < 3; k++)
			    P[place].GravAccel[k] += GravDataOut[j + noffset[recvTask]].u.Acc[k];

			  P[place].GravCost += GravDataOut[j + noffset[recvTask]].w.Ninteractions;
			}
		    }
		}

	      for(j = 0; j < NTask; j++)
		if((j ^ ngrp) < NTask)
		  nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	    }
	  tend = second();
	  timecommsumm += timediff(tstart, tend);

	  level = ngrp - 1;
	}

      MPI_Allgather(&ndone, 1, MPI_INT, ndonelist, 1, MPI_INT, MPI_COMM_WORLD);
      for(j = 0; j < NTask; j++)
	ntotleft -= ndonelist[j];
    }

  free(ndonelist);
  free(nsend);
  free(nsend_local);
  free(nbuffer);
  free(noffset);

  /* now add things for comoving integration */

#ifndef PERIODIC
#ifndef PMGRID
  if(All.ComovingIntegrationOn)
    {
      fac = 0.5 * All.Hubble * All.Hubble * All.Omega0 / All.G;
      //printf("XXXX %g %g\n",All.Hubble,All.Omega0);

      for(i = 0; i < NumPart; i++)
	if(P[i].Ti_endstep == All.Ti_Current)
	  for(j = 0; j < 3; j++)
	    P[i].GravAccel[j] += fac * P[i].Pos[j];
    }
#endif
#endif

  for(i = 0; i < NumPart; i++)
    if(P[i].Ti_endstep == All.Ti_Current)
      {
#ifdef PMGRID
	ax = P[i].GravAccel[0] + P[i].GravPM[0] / All.G;
	ay = P[i].GravAccel[1] + P[i].GravPM[1] / All.G;
	az = P[i].GravAccel[2] + P[i].GravPM[2] / All.G;
#else
	ax = P[i].GravAccel[0];
	ay = P[i].GravAccel[1];
	az = P[i].GravAccel[2];
#endif
	P[i].OldAcc = sqrt(ax * ax + ay * ay + az * az);
      }


  if(All.TypeOfOpeningCriterion == 1)
    All.ErrTolTheta = 0;	/* This will switch to the relative opening criterion for the following force computations */

  /*  muliply by G */
  for(i = 0; i < NumPart; i++)
    if(P[i].Ti_endstep == All.Ti_Current)
      for(j = 0; j < 3; j++)
	P[i].GravAccel[j] *= All.G;


  /* Finally, the following factor allows a computation of a cosmological simulation 
     with vacuum energy in physical coordinates */
#ifndef PERIODIC
#ifndef PMGRID
  if(All.ComovingIntegrationOn == 0)
    {
      fac = All.OmegaLambda * All.Hubble * All.Hubble;

      for(i = 0; i < NumPart; i++)
	if(P[i].Ti_endstep == All.Ti_Current)
	  for(j = 0; j < 3; j++)
	    P[i].GravAccel[j] += fac * P[i].Pos[j];
    }
#endif
#endif

#ifdef SELECTIVE_NO_GRAVITY
  for(i = 0; i < NumPart; i++)
    if(P[i].Ti_endstep < 0)
      P[i].Ti_endstep = -P[i].Ti_endstep - 1;
#endif

  if(ThisTask == 0)
    printf("tree is done.\n");

#else /* gravity is switched off */

  for(i = 0; i < NumPart; i++)
    if(P[i].Ti_endstep == All.Ti_Current)
      for(j = 0; j < 3; j++)
	P[i].GravAccel[j] = 0;

#endif




  /* Now the force computation is finished */

  /*  gather some diagnostic information */

  timetreelist = malloc(sizeof(double) * NTask);
  timecommlist = malloc(sizeof(double) * NTask);
  costtreelist = malloc(sizeof(double) * NTask);
  numnodeslist = malloc(sizeof(int) * NTask);
  ewaldlist = malloc(sizeof(double) * NTask);
  nrecv = malloc(sizeof(int) * NTask);

  numnodes = Numnodestree;

  MPI_Gather(&costtotal, 1, MPI_DOUBLE, costtreelist, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&numnodes, 1, MPI_INT, numnodeslist, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(&timetree, 1, MPI_DOUBLE, timetreelist, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&timecommsumm, 1, MPI_DOUBLE, timecommlist, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Gather(&NumPart, 1, MPI_INT, nrecv, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(&ewaldcount, 1, MPI_DOUBLE, ewaldlist, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Reduce(&nexportsum, &nexport, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&timeimbalance, &sumimbalance, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      All.TotNumOfForces += ntot;

      fprintf(FdTimings, "Step= %d  t= %g  dt= %g \n", All.NumCurrentTiStep, All.Time, All.TimeStep);
      fprintf(FdTimings, "Nf= %d%09d  total-Nf= %d%09d  ex-frac= %g  iter= %d\n",
	      (int) (ntot / 1000000000), (int) (ntot % 1000000000),
	      (int) (All.TotNumOfForces / 1000000000), (int) (All.TotNumOfForces % 1000000000),
	      nexport / ((double) ntot), iter);
      /* note: on Linux, the 8-byte integer could be printed with the format identifier "%qd", but doesn't work on AIX */

      fac = NTask / ((double) All.TotNumPart);

      for(i = 0, maxt = timetreelist[0], sumt = 0, plb_max = 0,
	  maxnumnodes = 0, costtotal = 0, sumcomm = 0, ewaldtot = 0; i < NTask; i++)
	{
	  costtotal += costtreelist[i];

	  sumcomm += timecommlist[i];

	  if(maxt < timetreelist[i])
	    maxt = timetreelist[i];
	  sumt += timetreelist[i];

	  plb = nrecv[i] * fac;

	  if(plb > plb_max)
	    plb_max = plb;

	  if(numnodeslist[i] > maxnumnodes)
	    maxnumnodes = numnodeslist[i];

	  ewaldtot += ewaldlist[i];
	}
      fprintf(FdTimings, "work-load balance: %g  max=%g avg=%g PE0=%g\n",
	      maxt / (sumt / NTask), maxt, sumt / NTask, timetreelist[0]);
      fprintf(FdTimings, "particle-load balance: %g\n", plb_max);
      fprintf(FdTimings, "max. nodes: %d, filled: %g\n", maxnumnodes,
	      maxnumnodes / (All.TreeAllocFactor * All.MaxPart));
      fprintf(FdTimings, "part/sec=%g | %g  ia/part=%g (%g)\n", ntot / (sumt + 1.0e-20),
	      ntot / (maxt * NTask), ((double) (costtotal)) / ntot, ((double) ewaldtot) / ntot);
      fprintf(FdTimings, "\n");

      fflush(FdTimings);

      All.CPU_TreeWalk += sumt / NTask;
      All.CPU_Imbalance += sumimbalance / NTask;
      All.CPU_CommSum += sumcomm / NTask;
    }


  // DY: external potential, following Volker's suggestion 
  // Pos: comoving coordinates in units h^-1 kpc
  // physical pos: Pos*a/h = Pos*All.Time/All.HubbleParam
  // printf("\n All.Time %g h: %g \n",All.Time,All.HubbleParam);


  // 2023.03.31 ZXY: Begin MW potential ......
  // 2024.03.12 ZXY: Add a switch of external potential ......
  // 2024.04.17 ZXY: Revise MW potential according to McMillan (2017) ...
  // 2024.04.18 ZXY: According to McMillan (2017), there are other two gas discs, named d3 and d4 here ...

  if(All.ExternalPotential == 1){

  double M_d1 = 3.5186,   a_d1 = 2.50,   b_d1 = 0.3;
  double M_d2 = 1.0487,   a_d2 = 3.02,   b_d2 = 0.9;
  double M_d3 = 1.1,      a_d3 = 7,      b_d3 = 0.085;
  double M_d4 = 0.12,     a_d4 = 1.5,    b_d4 = 0.045;

  double M_Hern = 0.923,   r_Hern = 1.3;
  double rhos_NFW = 8.54e-4,  rs_NFW = 19.6;

  for(i=0; i < NumPart; i++)
     if(P[i].Ti_endstep==All.Ti_Current)
     {
        //double rphys = All.Time/All.HubbleParam*sqrt((P[i].Pos[0]) * (P[i].Pos[0]) + (P[i].Pos[1]) * (P[i].Pos[1])
        //               + (P[i].Pos[2]) * (P[i].Pos[2]));
        double rad_ZXY   = sqrt((P[i].Pos[0]) * (P[i].Pos[0]) + (P[i].Pos[1]) * (P[i].Pos[1])
                  + (P[i].Pos[2]) * (P[i].Pos[2]));
        double Rad_ZXY = sqrt( (P[i].Pos[0]) * (P[i].Pos[0]) + (P[i].Pos[1]) * (P[i].Pos[1]) );
        // Use rphys to compute mass, use comoving distance to compute acceleration

        //   2023.03.31 ZXY: MW external potential  ----  two disk + Buldge + Halo ......

        double FR_d1 = FR_disk_ZXY(M_d1, a_d1, b_d1, Rad_ZXY, P[i].Pos[2]);
        double FR_d2 = FR_disk_ZXY(M_d2, a_d2, b_d2, Rad_ZXY, P[i].Pos[2]);
        double FR_d3 = FR_disk_ZXY(M_d3, a_d3, b_d3, Rad_ZXY, P[i].Pos[2]);
        double FR_d4 = FR_disk_ZXY(M_d4, a_d4, b_d4, Rad_ZXY, P[i].Pos[2]);

        double Fz_d1 = Fz_disk_ZXY(M_d1, a_d1, b_d1, Rad_ZXY, P[i].Pos[2]);
        double Fz_d2 = Fz_disk_ZXY(M_d2, a_d2, b_d2, Rad_ZXY, P[i].Pos[2]);
        double Fz_d3 = Fz_disk_ZXY(M_d3, a_d3, b_d3, Rad_ZXY, P[i].Pos[2]);
        double Fz_d4 = Fz_disk_ZXY(M_d4, a_d4, b_d4, Rad_ZXY, P[i].Pos[2]);

        double Fr_H = Fr_Hern_ZXY(M_Hern, r_Hern, rad_ZXY);
        double Fr_N = Fr_NFW_ZXY(rhos_NFW, rs_NFW, rad_ZXY);

        P[i].GravAccel[0] += (FR_d1 + FR_d2 + FR_d3 + FR_d4) * P[i].Pos[0]/Rad_ZXY + (Fr_H + Fr_N) * P[i].Pos[0]/rad_ZXY;
        P[i].GravAccel[1] += (FR_d1 + FR_d2) * P[i].Pos[1] / Rad_ZXY + (Fr_H + Fr_N) * P[i].Pos[1]/rad_ZXY;
        P[i].GravAccel[2] += Fz_d1 + Fz_d2 + (Fr_H + Fr_N) * P[i].Pos[2]/rad_ZXY;
/*        
        if(i%100000 == 0)
        {
          fprintf(stderr,"***********************************ZXY: test of gravtree.c *************************************\n");
          fprintf(stderr,"Particle %d,  x = %g,   y = %g,   z = %g \n", i, P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
          fprintf(stderr,"FR_d1 = %g,  FR_d2 = %g,  Fz_d1 = %g,  Fz_d2 = %g,  Fr_H = %g,  Fr_N = %g \n", FR_d1, FR_d2, Fz_d1, Fz_d2, Fr_H, Fr_N);
          fprintf(stderr,"Fx = %g,  Fy = %g,  Fz = %g  \n", (FR_d1 + FR_d2) * P[i].Pos[0] / Radi + (Fr_H + Fr_N) * P[i].Pos[0] / rad, (FR_d1 + FR_d2) * P[i].Pos[1] / Radi + (Fr_H + Fr_N) * P[i].Pos[1] / rad, Fz_d1 + Fz_d2 + (Fr_H + Fr_N) * P[i].Pos[2] / rad);
        }
      // 2023.04.18 ZXY: Validity of force calculation has been done ... Compared with MMA force function ......
*/

     }
  }

  free(nrecv);
  free(ewaldlist);
  free(numnodeslist);
  free(costtreelist);
  free(timecommlist);
  free(timetreelist);
}



/*! This function sets the (comoving) softening length of all particle
 *  types in the table All.SofteningTable[...].  We check that the physical
 *  softening length is bounded by the Softening-MaxPhys values.
 */
void set_softenings(void)
{
  int i;

  if(All.ComovingIntegrationOn)
    {
      if(All.SofteningGas * All.Time > All.SofteningGasMaxPhys)
        All.SofteningTable[0] = All.SofteningGasMaxPhys / All.Time;
      else
        All.SofteningTable[0] = All.SofteningGas;
      
      if(All.SofteningHalo * All.Time > All.SofteningHaloMaxPhys)
        All.SofteningTable[1] = All.SofteningHaloMaxPhys / All.Time;
      else
        All.SofteningTable[1] = All.SofteningHalo;
      
      if(All.SofteningDisk * All.Time > All.SofteningDiskMaxPhys)
        All.SofteningTable[2] = All.SofteningDiskMaxPhys / All.Time;
      else
        All.SofteningTable[2] = All.SofteningDisk;
      
      if(All.SofteningBulge * All.Time > All.SofteningBulgeMaxPhys)
        All.SofteningTable[3] = All.SofteningBulgeMaxPhys / All.Time;
      else
        All.SofteningTable[3] = All.SofteningBulge;
      
      if(All.SofteningStars * All.Time > All.SofteningStarsMaxPhys)
        All.SofteningTable[4] = All.SofteningStarsMaxPhys / All.Time;
      else
        All.SofteningTable[4] = All.SofteningStars;
      
      if(All.SofteningBndry * All.Time > All.SofteningBndryMaxPhys)
        All.SofteningTable[5] = All.SofteningBndryMaxPhys / All.Time;
      else
        All.SofteningTable[5] = All.SofteningBndry;
    }
  else
    {
      All.SofteningTable[0] = All.SofteningGas;
      All.SofteningTable[1] = All.SofteningHalo;
      All.SofteningTable[2] = All.SofteningDisk;
      All.SofteningTable[3] = All.SofteningBulge;
      All.SofteningTable[4] = All.SofteningStars;
      All.SofteningTable[5] = All.SofteningBndry;
    }

  for(i = 0; i < 6; i++)
    All.ForceSoftening[i] = 2.8 * All.SofteningTable[i];

  All.MinGasHsml = All.MinGasHsmlFractional * All.ForceSoftening[0];
}


/*! This function is used as a comparison kernel in a sort routine. It is
 *  used to group particles in the communication buffer that are going to
 *  be sent to the same CPU.
 */
int grav_tree_compare_key(const void *a, const void *b)
{
  if(((struct gravdata_index *) a)->Task < (((struct gravdata_index *) b)->Task))
    return -1;

  if(((struct gravdata_index *) a)->Task > (((struct gravdata_index *) b)->Task))
    return +1;

  return 0;
}
