#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"
//#define DYDEBUG

void sidm(double dt_drift, double dt_gravkick){

    long long ntot,nsumscatt,nsumngb;
    int i, j, ibullet, sendTask, source;
    int nexport, *numlist, *ndonelist;

    int *noffset, *nbuffer, *nsend, *nsend_local;
    long long ntotleft;
    int ndone, maxfill, ngrp;
    int k, place;
    int level, recvTask;
    MPI_Status status;
    double TypicalDist;

    // make it comoving...
    TypicalDist = All.ForceSoftening[1]; // as a radius of the search region... 
    //printf(">>>> Typical dist %g \n",TypicalDist); 

    // loop over all the particles on this PE, findout the range of velocities
    for(i = 0; i < NumPart; i++){ P[i].tscatt=0; }

    // Each PE contains a top-level tree, where missing data on other PEs  
    // are represented by pseudo-particles. 

    // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    //offsets of bunches in common list
    noffset = malloc(sizeof(int) * NTask); 
    nbuffer = malloc(sizeof(int) * NTask);
    nsend_local = malloc(sizeof(int) * NTask);
    nsend = malloc(sizeof(int) * NTask * NTask);
    ndonelist = malloc(sizeof(int) * NTask);

    numlist = malloc(NTask * sizeof(int) * NTask);
    MPI_Allgather(&NumPart, 1, MPI_INT, numlist, 1, MPI_INT, MPI_COMM_WORLD);

    for(i = 0, ntot = 0; i < NTask; i++)
	ntot += numlist[i];
    free(numlist);

    if(ThisTask == 0) printf(">>> Begin scattering with TypicalDist = %g.\n", TypicalDist); 

    ntotscat=0;
    ntotngb=0;
    ibullet = 0;                    /* begin with this index */
    /* Gas particles are always in the beginning */
    ntotleft = ntot;                    /* particles left for all tasks together */

    while(ntotleft > 0)
    {
	for(j = 0; j < NTask; j++)
	    nsend_local[j] = 0;

	/* do local particles and prepare export list */
	for(nexport = 0, ndone = 0; ibullet < NumPart && nexport < All.BunchSizeSIDM - NTask; ibullet++){
	    ndone++;

	    for(j = 0; j < NTask; j++) Exportflag[j] = 0;

	    doscatt(ibullet,0,dt_drift, dt_gravkick, TypicalDist);

	    for(j = 0; j < NTask; j++){
		if(Exportflag[j])
		{
		    for(k = 0; k < 3; k++){
			SIDMDataIn[nexport].Pos[k] = P[ibullet].Pos[k];
			SIDMDataIn[nexport].Vel[k] = P[ibullet].Vel[k];
		    }
		    SIDMDataIn[nexport].Type = P[ibullet].Type;
		    SIDMDataIn[nexport].Mass = P[ibullet].Mass;
		    //SIDMDataIn[nexport].OldAcc = P[ibullet].OldAcc;
		    SIDMDataIn[nexport].tscatt = P[ibullet].tscatt;
		    SIDMDataIn[nexport].sendTask = ThisTask; 
		    SIDMDataIn[nexport].Task = j; 
		    SIDMDataIn[nexport].Index = ibullet;
		    //SIDMDataIn[nexport].ID = P[ibullet].ID;
		    //SIDMDataIn[nexport].Pijmax = P[ibullet].Pijmax;
		    nexport++;
		    nsend_local[j]++;
		}
	    }
	}
	// Local particles done

	qsort(SIDMDataIn, nexport, sizeof(struct sidmdata_in), sidm_compare_key);

	for(j = 1, noffset[0] = 0; j < NTask; j++)
	    noffset[j] = noffset[j - 1] + nsend_local[j - 1];

	MPI_Allgather(nsend_local, NTask, MPI_INT, nsend, NTask, MPI_INT, MPI_COMM_WORLD);

#ifdef DYDEBUG
	printf(" \n >>> Starting communication phase >>> \n ");
#endif
	/* now do the particles that need to be exported */
	for(level = 1; level < (1 << PTask); level++) /* PTask: smallest integer such that NTask <= 2^PTask */
	{
	    for(j = 0; j < NTask; j++)
		nbuffer[j] = 0;
	    for(ngrp = level; ngrp < (1 << PTask); ngrp++)
	    {
#ifdef DYDEBUG
		if(level>1) printf("level = %d, ngrp= %d \n",level,ngrp);
#endif
		maxfill = 0;
		for(j = 0; j < NTask; j++)
		{
		    if((j ^ ngrp) < NTask)
			if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
			    maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
		}
		/*!< All.BunchSizeSIDM: number of particles fitting into the buffer in the parallel tree-force algorithm  */
		if(maxfill >= All.BunchSizeSIDM)
		{ printf(" maxfill >= All.BunchSizeSIDM "); break; }

		sendTask = ThisTask;
		recvTask = ThisTask ^ ngrp;

		if(recvTask < NTask)
		{
		    if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
			/* get the particles */
			// Send and receive from Task=recvTask
			MPI_Sendrecv(&SIDMDataIn[noffset[recvTask]],
				nsend_local[recvTask] * sizeof(struct sidmdata_in), MPI_BYTE,
				recvTask, TAG_SIDM_A,
				&SIDMDataGet[nbuffer[ThisTask]],
				nsend[recvTask * NTask + ThisTask] * sizeof(struct sidmdata_in), MPI_BYTE,
				recvTask, TAG_SIDM_A, MPI_COMM_WORLD, &status);
		    }
		}

		for(j = 0; j < NTask; j++)
		    if((j ^ ngrp) < NTask)
			nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	    }

	    if(ThisTask==0) printf("@@@@ Check # of particles in the SIDM buffer: %d, PE: %d \n",nbuffer[ThisTask],ThisTask);

	    for(j = 0; j < nbuffer[ThisTask]; j++)
	    {
		// Do scatterings for external particles  
		doscatt(j,1,dt_drift, dt_gravkick, TypicalDist);
	    }

	    /* do a block to explicitly measure imbalance */
	    //tstart = second();
	    //MPI_Barrier(MPI_COMM_WORLD);
	    //tend = second();
	    //timeimbalance += timediff(tstart, tend);

#ifdef DYDEBUG
	    if(level > 1) printf("############################################### Getting results ThisTask = %d\n",ThisTask);
#endif
	    /* get the result */
	    for(j = 0; j < NTask; j++)
		nbuffer[j] = 0;
	    for(ngrp = level; ngrp < (1 << PTask); ngrp++)
	    {
#ifdef DYDEBUG
		if(level>1) printf("level = %d, ngrp= %d \n",level,ngrp);
#endif
		maxfill = 0;
		for(j = 0; j < NTask; j++)
		{
		    if((j ^ ngrp) < NTask)
			if(maxfill < nbuffer[j] + nsend[(j ^ ngrp) * NTask + j])
			    maxfill = nbuffer[j] + nsend[(j ^ ngrp) * NTask + j];
		}
		if(maxfill >= All.BunchSizeSIDM)
		{ printf(" All.BunchSizeSIDM... "); break; }

		sendTask = ThisTask;
		recvTask = ThisTask ^ ngrp;

		if(recvTask < NTask) 
		{
		    if(nsend[ThisTask * NTask + recvTask] > 0 || nsend[recvTask * NTask + ThisTask] > 0)
		    {
			/* send the results */
			MPI_Sendrecv(&SIDMDataResult[nbuffer[ThisTask]],
				nsend[recvTask * NTask + ThisTask] * sizeof(struct sidmdata_in),
				MPI_BYTE, recvTask, TAG_SIDM_B,
				&SIDMDataOut[noffset[recvTask]],
				nsend_local[recvTask] * sizeof(struct sidmdata_in),
				MPI_BYTE, recvTask, TAG_SIDM_B, MPI_COMM_WORLD, &status);

			/* add the result to the particles */
			if( recvTask<sendTask ){ // receive updates only from PEs of lower task 
			    for(j = 0; j < nsend_local[recvTask]; j++)
			    {
				source = j + noffset[recvTask];
				place = SIDMDataIn[source].Index;

				if(SIDMDataOut[source].tscatt != 0){ // Need to be fixed
				    for(k = 0; k < 3; k++){
					P[place].Pos[k] = SIDMDataOut[source].Pos[k];
					P[place].Vel[k] = SIDMDataOut[source].Vel[k];
				    }
				    P[place].tscatt = SIDMDataOut[source].tscatt;
				    //P[place].ID = SIDMDataOut[source].ID;
				}
			    }
			}
		    }
		}

		for(j = 0; j < NTask; j++)
		    if((j ^ ngrp) < NTask)
			nbuffer[j] += nsend[(j ^ ngrp) * NTask + j];
	    }
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

    if(ThisTask == 0)
	printf(">>> End scattering.\n");

    MPI_Reduce(&ntotscat, &nsumscatt, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if(ThisTask == 0) printf("~~~~~~~~~~~~~~~NtotScatt=%d\n",nsumscatt);

#ifdef DYDEBUG
    MPI_Reduce(&ntotngb, &nsumngb, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if(ThisTask == 0) printf("~~~~~~~~~~~~~~~Ntot neighbors=%d\n",nsumngb);
#endif

}


/*! This routine is a comparison kernel used in a sort routine to group
 *  *  particles that are exported to the same processor.
 *   */
int sidm_compare_key(const void *a, const void *b)
{
    if(((struct sidmdata_in *) a)->Task < (((struct sidmdata_in *) b)->Task))
	return -1;

    if(((struct sidmdata_in *) a)->Task > (((struct sidmdata_in *) b)->Task))
	return +1;

    return 0;
}

