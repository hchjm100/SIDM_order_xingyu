#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"

//#define DYDEBUG

#ifdef PERIODIC
/*! Macro that maps a distance to the nearest periodic neighbour */
#define NEAREST(x) (((x)>boxhalf)?((x)-boxsize):(((x)<-boxhalf)?((x)+boxsize):(x)))
/*! Size of 3D lock-up table for Ewald correction force */
#define EN  64
#endif

/* ref: ngb_treefind_pairs */
/* identify all scattering pairs, then doscattering in doscatt.c */
// parsing an index, should use FLOAT, as defined for P[0].Pos
int ngb_sidm(FLOAT searchcenter[3], double TypicalDist, int *startnode){

    struct NODE *nop = 0;
    int no, i, j, p, numngb;
    FLOAT dx, dy, dz, eff_dist;
    FLOAT pos_x, pos_y, pos_z;

    numngb = 0;
    no = *startnode;

    pos_x = searchcenter[0];
    pos_y = searchcenter[1];
    pos_z = searchcenter[2];

    while(no >= 0)
    {

	if(no < All.MaxPart) // All.MaxPart: maximum number of particles that can be stored on one PE. 
	{
            p = no;
            no = Nextnode[no];

            if(P[p].Type!=1) continue;
	    //if(no==bullet){ no = Nextnode[no]; continue; } // should in doscatt  

	    /*In this case, the index of the node is the index of the particle */
	    dx = P[p].Pos[0] - pos_x;
	    dy = P[p].Pos[1] - pos_y;
	    dz = P[p].Pos[2] - pos_z;
#ifdef PERIODIC
	    dx = NEAREST(dx);
	    dy = NEAREST(dy);
	    dz = NEAREST(dz);
#endif

	    eff_dist = (FLOAT) TypicalDist;

	    if(dx < -eff_dist || dx > eff_dist) continue; 
	    if(dy < -eff_dist || dy > eff_dist) continue; 
	    if(dz < -eff_dist || dz > eff_dist) continue; 

	    // ############### Local neighboring particles ###################

	    Ngblist[numngb++] = p;

	    if(numngb == MAX_NGB)
	    {
		if(ThisTask==0) printf
		    ("ThisTask=%d: Need to do another neighbour loop in sidm for (%g|%g|%g) TypicalDist=%g no=%d\n",
		     ThisTask, searchcenter[0], searchcenter[1], searchcenter[2], TypicalDist, no);
		*startnode = no;
		return numngb;
	    }
	}
	else                      /* we have an  internal node */
	{

            // see which PE(s) should we export the bullet to...
	    if(no >= All.MaxPart + MaxNodes)      /* pseudo particle */
	    {
		Exportflag[DomainTask[no - (All.MaxPart + MaxNodes)]] = 1;
		no = Nextnode[no - MaxNodes];
		continue;
	    }

	    nop = &Nodes[no];
	    no = nop->u.d.sibling; /* in case the node can be discarded */

	    // Check if we need to open the node
	    eff_dist = TypicalDist + 0.5 * nop->len;

	    dx = nop->center[0] - pos_x; 
	    dy = nop->center[1] - pos_y; 
	    dz = nop->center[2] - pos_z;
#ifdef PERIODIC
	    dx = NEAREST(dx);
	    dy = NEAREST(dy);
	    dz = NEAREST(dz);
#endif

	    // If we are outside the node, check if it is nearby
	    // In most cases, the tree walk stop here...
	    if(dx < -eff_dist || dx > eff_dist) continue;
	    if(dy < -eff_dist || dy > eff_dist) continue;
	    if(dz < -eff_dist || dz > eff_dist) continue;
	    // The survived nodes should be very close to the bullet particle

	    no = nop->u.d.nextnode; /* ok, we need to open the node */
	}
    }

    *startnode = -1; // no==-1 once reached here;
    return numngb;
}

