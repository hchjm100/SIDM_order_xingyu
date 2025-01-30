#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"

//#define DYDEBUG
//#define DYDEBUGKIN
//#define DYCHECKTH

#ifdef PERIODIC
/*! Macro that maps a distance to the nearest periodic neighbour */
#define NEAREST(x) (((x)>boxhalf)?((x)-boxsize):(((x)<-boxhalf)?((x)+boxsize):(x)))
/*! Size of 3D lock-up table for Ewald correction force */
#define EN  64
#endif

int doscatt(int bullet, int mode,double dt_drift, double dt_gravkick, double TypicalDist){

#ifdef PERIODIC
    double boxsize, boxhalf;

    boxsize = All.BoxSize;
    boxhalf = 0.5 * All.BoxSize;
#endif

    struct NODE *nop = 0;
    int no, startnode, numngb, i, j, n, ptype, inttype, target;
    double vr2, vr, dx, dy, dz, mass, massH; 
    double dr2,dr,dt_scatt,b,costh,sinth,costhp,sinthp,phip,thp;
    double eff_dist;
    unsigned int idcheck = 1220152896;
    //unsigned int idcheck = 1226648464;

    double sigma_m_inel,sigma_m_el,sigma_m_tot,sigma0;
    double vw,vw2;
    double *pos;
    double pos_x, pos_y, pos_z;
    double vb_x, vb_y, vb_z;
    double vr_x, vr_y, vr_z;
    double vcmx, vcmy, vcmz;
    double vhatbp_x, vhatbp_y, vhatbp_z;
    double vrp_x, vrp_y, vrp_z;
    double ex,ey,ez;
    double cdx_x,cdx_y,cdx_z,ncd;
    double cdy_x,cdy_y,cdy_z;
    double cdz_x,cdz_y,cdz_z;
    double scattlength;
    double Probij,thisrnd,xrnd;
    //----------------

    //sigma0   = 24324.6;
    //vw=20.0*All.Time; // code vel variable u=a^2 dotx =a * vp
    //vw2=vw*vw;

    inttype = 0; // initialization, now as long as it is >0, we perform scattering. It is not used to distinguish elastic/inelastic...



    // ZXY 2022.02.04: deltaf the fractional mass difference, ratio_ds the ratio of down scatter ...... 
    // -------------------------------------------

    double deltaf = 0;
    double ratio_ds = 0;
    double Delta, m1, m2;

    // -------------------------------------------


    if(mode == 0) // The particle is on this PE
    {
	ptype = P[bullet].Type;

        // ZXY 2022.02.04: Now consider down scatter, not only type 1 can scatter......
        //                 Set lower level state to be type 2, which can still have elastic scatter ......
	if( (ptype!=1) && (ptype!=2) ) return 1;

    // predicted particle velocity for bullet
    double pred_vel_bu_x = P[bullet].Vel[0] + P[bullet].GravAccel[0] * dt_gravkick;
    double pred_vel_bu_y = P[bullet].Vel[1] + P[bullet].GravAccel[1] * dt_gravkick;
    double pred_vel_bu_z = P[bullet].Vel[2] + P[bullet].GravAccel[2] * dt_gravkick;

	pos = P[bullet].Pos;
	pos_x = P[bullet].Pos[0];
	pos_y = P[bullet].Pos[1];
	pos_z = P[bullet].Pos[2];
	vb_x = pred_vel_bu_x;
	vb_y = pred_vel_bu_y;
	vb_z = pred_vel_bu_z;
	mass = P[bullet].Mass;
    }
    else // mode==1, exported particles, working on the buffer
    {
	ptype = SIDMDataGet[bullet].Type;

        // ZXY 2022.02.04: The same as mode = 0 ......
	if( (ptype!=1) && (ptype!=2) ) return 1;

    // predicted exported particle velocity for bullet
    double pred_vel_ex_bu_x = SIDMDataGet[bullet].Vel[0] + GravDataGet[bullet].u.Acc[0] * dt_gravkick;
    double pred_vel_ex_bu_y = SIDMDataGet[bullet].Vel[1] + GravDataGet[bullet].u.Acc[1] * dt_gravkick;
    double pred_vel_ex_bu_z = SIDMDataGet[bullet].Vel[2] + GravDataGet[bullet].u.Acc[2] * dt_gravkick;

	// To avoid double counting, we consider pseudo particles only from PEs of higher task.
        // This operation necessarily remove some neighbors from considerations
	if(ThisTask > SIDMDataGet[bullet].sendTask) return 1;
	//if(SIDMDataGet[bullet].tscatt==1) return 1;
	pos   = SIDMDataGet[bullet].Pos;
	pos_x = SIDMDataGet[bullet].Pos[0];
	pos_y = SIDMDataGet[bullet].Pos[1];
	pos_z = SIDMDataGet[bullet].Pos[2];
	vb_x = pred_vel_ex_bu_x;
	vb_y = pred_vel_ex_bu_y;
	vb_z = pred_vel_ex_bu_z;
	mass = SIDMDataGet[bullet].Mass;
    }

    /* Now start the SIDM computation for this particle */
    startnode = All.MaxPart;

    do
    {
	numngb = ngb_sidm(&pos[0], TypicalDist, &startnode);
        
        randomize(Ngblist,numngb);

	for(n = 0; n < numngb; n++)
	{
	    target = Ngblist[n];

            #ifdef DYDEBUG
                if(mode==0 && P[bullet].ID==idcheck && P[target].ID!=idcheck){
                    ntotngb++;
                }
                // SIDMDataGet[bullet].ID cannot be idcheck, ThisTask!=Task of the event
                else if(mode==1 && P[target].ID == idcheck ){
                    ntotngb++;
                }
            #endif

            if(target==bullet) continue; 

            //-------- Reduce to spherical search region --------

            dx = P[target].Pos[0] - pos_x;
            dy = P[target].Pos[1] - pos_y;
            dz = P[target].Pos[2] - pos_z;
#ifdef PERIODIC
            dx = NEAREST(dx);
            dy = NEAREST(dy);
            dz = NEAREST(dz);
#endif
            dr2 = dx * dx + dy * dy + dz * dz;
            if(dr2 > TypicalDist*TypicalDist) continue;
               
            //-------- Setup coordinate system --------

            // predicted particle velocity for target
            double pred_vel_ta_x = P[target].Vel[0] + P[target].GravAccel[0] * dt_gravkick;
            double pred_vel_ta_y = P[target].Vel[1] + P[target].GravAccel[1] * dt_gravkick;
            double pred_vel_ta_z = P[target].Vel[2] + P[target].GravAccel[2] * dt_gravkick;

            vr_x = vb_x - pred_vel_ta_x;
            vr_y = vb_y - pred_vel_ta_y;
            vr_z = vb_z - pred_vel_ta_z;
            vr2 = vr_x * vr_x + vr_y * vr_y + vr_z * vr_z;
            vr = sqrt(vr2);

            dr = sqrt(dr2);
            costh = (dx*vr_x+dy*vr_y+dz*vr_z)/(dr*vr);
            sinth = sqrt(1.0-costh*costh);
  
            //sigma_m_el = sigma0/(1.0+vr2/vw2); // total cross section  
            sigma_m_tot = 50;

	    // In order that everything is in comoving coordinate, # of collisions should remain the same
            // n sigma/m (v dt) = (n a^3) (sigma/m/a^2) dr/a = (n a^3) (sigma/m/a^2) d_x
            // vr*dt_drift = d_x
            // sigma/m cm^2/g -> sigma/m *h (cm/h)^2/(g/h)
            // double Probij = (sigma_m_tot*All.HubbleParam/pow(All.Time,2))*mass*All.UnitMass_in_g*pow(All.UnitLength_in_cm,-2)*vr*dt_drift/(4*3.14159265358979*pow(TypicalDist,3)/3);
            
            // ZXY 2022.02.08: In non-cosmological simulation, we cannot use this probability again ......
            // -------------------------------------------------------------------------------------------
/*
            double Probij = (sigma_m_tot*All.HubbleParam/pow(All.Time,2))*mass*All.UnitMass_in_g*pow(All.UnitLength_in_cm,-2)*vr*dt_drift;
            Probij=Probij*get_Wij(TypicalDist,dr);
	    if(Probij>0.1) printf("XXXXXXXX Pij=%g >1, scattering prob too large! \n",Probij);

            // consider half the probability
            if(mode == 0){   if(gsl_rng_uniform(random_generator)>Probij/2) continue; }
            // consider the full probability
            else{            if(gsl_rng_uniform(random_generator)>Probij) continue; }
*/
            // ------------------------------------------------------------------------------------------


            b = dr * sinth;
            scattlength = sqrt(sigma_m_tot * mass*All.UnitMass_in_g / 3.14159265358979) / All.UnitLength_in_cm;
            dt_scatt = fabs(dr*costh/vr);

            // ZXY 2022.02.08: Determine whether scatter will happen ......
            if ( costh < 0 || b > scattlength || dt_scatt > dt_drift || P[target].Type > 2 ) continue;

	    cdz_x = vr_x/vr;
            cdz_y = vr_y/vr;
            cdz_z = vr_z/vr;

            vcmx = 0.5*(vb_x+pred_vel_ta_x);
            vcmy = 0.5*(vb_y+pred_vel_ta_y);
            vcmz = 0.5*(vb_z+pred_vel_ta_z);

            if(vcmx==0&&vcmy==0&&vcmz==0){
               printf("XXXXXXXXX DY: vcm==0, extremely rare... %g\n");
               double norm=sqrt(cdz_y*cdz_y+cdz_z*cdz_z);
               cdx_x = 0;
               cdx_y = -cdz_z/norm;
               cdx_z = cdz_y/norm;
            }
            else{
                 cdx_x = (cdz_y*vcmz-cdz_z*vcmy);
                 cdx_y = (cdz_z*vcmx-cdz_x*vcmz);
                 cdx_z = (cdz_x*vcmy-cdz_y*vcmx);
                 ncd = sqrt(cdx_x*cdx_x+cdx_y*cdx_y+cdx_z*cdx_z);
                 cdx_x/=ncd;
                 cdx_y/=ncd;
                 cdx_z/=ncd;
            }
            cdy_x = cdz_y*cdx_z-cdz_z*cdx_y;
            cdy_y = cdz_z*cdx_x-cdz_x*cdx_z;
            cdy_z = cdz_x*cdx_y-cdz_y*cdx_x;
      
            //xrnd=gsl_rng_uniform(random_generator);
            //costhp = (vr2+vw2-vr2*xrnd-2.0*vw2*xrnd)/(vr2+vw2-vr2*xrnd) ; // valued in (-1,1) 
            // Rejection sampling of the costhp distribution, 
            // Rutherford or Moller scattering 
            //double fmax=fcosthR(1,vr,vw); 
            //double fmax=fcosthM(1,vr,vw); 
            //double urnd=fmax*gsl_rng_uniform(random_generator);
            //double yrnd=2*gsl_rng_uniform(random_generator)-1;
            //while(urnd>fcosthR(yrnd,vr,vw)){
            //while(urnd>fcosthM(yrnd,vr,vw)){
            //    urnd=fmax*gsl_rng_uniform(random_generator);
            //    yrnd=2*gsl_rng_uniform(random_generator)-1;
            //}
            //costhp = yrnd;
            costhp = 2.0*gsl_rng_uniform(random_generator)-1.0;
            sinthp = sqrt(1.-costhp*costhp);
            phip = gsl_rng_uniform(random_generator) * 2.0* 3.14159265358979;
	    if(costhp<0){
		costhp = -costhp;
		phip+=3.14159265358979;
	    }

            #ifdef DYCHECKTH
                FILE * fp;
                fp = fopen("points_thetaR.txt","a");
                fprintf (fp, "%g %g %g\n",acos(costhp)/3.141593*180,costhp,vr);
                //fprintf (fp, "%g %g %g %g %g %g \n",acos(costhp)/3.141593*180,costhp,acos(costhH)/3.141593*180,costhH,acos(costhH2)/3.141593*180,costhH2);
                fclose (fp);
            #endif

            // direction in the CM frame
            vhatbp_x = sinthp*cos(phip)*cdx_x + sinthp*sin(phip)*cdy_x + costhp*cdz_x;
            vhatbp_y = sinthp*cos(phip)*cdx_y + sinthp*sin(phip)*cdy_y + costhp*cdz_y;
            vhatbp_z = sinthp*cos(phip)*cdx_z + sinthp*sin(phip)*cdy_z + costhp*cdz_z;

            //printf(" XXXXXXXXXXXXXXXXXXXXXX %g %g %g %g %g %g \n",cdy_x,cdy_y,cdy_z,cdz_x,cdz_y,cdz_z);
            //-------------------------------------------

	    if(mode == 0)
	    {

            // ZXY 2022.02.04: Here is after cross section judgement, means that scatter will certainly occur ......
            //
            // ZXY 2022.02.04: 1)  Determine particle type, then decide scatter type accordingly.
            //
            //                 2)  If P.Type are both 1, consider probability  prob_d  and  prob_e  for down scatter
            //                     and elastic scatter ( prob_d + prob_e = 1 ).
            //
            //                     3) If down scatter, revise: Delta, m2, P.Type
            //
            //                        Else, Delta = 0, m2 = m1, P.Type remain unchanged ......
            //
            //                 4)  Else, consider only elastic scatter.
            //
            // ZXY 2022.02.04: Warning: P.Type is used to jump out of scatter in previous part !!!!!! ( Revised )
            // ZXY 2022.02.04: Warning: In sidm.c, P.Type may need to be updated ...... ( Revised )

            //---------------------------------------------------------------------------------------------
            double prob_d = gsl_rng_uniform(random_generator);

              if( (P[bullet].Type == 1) && (P[target].Type == 1) && prob_d < ratio_ds ){

                  Delta = mass*deltaf; // deltaf the fractional mass difference, defined previously ......
                  m1 = mass;
                  m2 = m1 - Delta;

                  P[bullet].Type = 2;
                  P[bullet].Mass = m2;

                  P[target].Type = 2;
                  P[target].Mass = m2;

                }
              else{
                  Delta = 0;
                  m1 = mass;
                  m2 = mass;
              }
            //---------------------------------------------------------------------------------------------

              //  double m1=mass;
              //  double m2=mass;
              //  double Delta=0;

                // predicted particle velocity for bullet
                double pred_vel_bu_x = P[bullet].Vel[0] + P[bullet].GravAccel[0] * dt_gravkick;
                double pred_vel_bu_y = P[bullet].Vel[1] + P[bullet].GravAccel[1] * dt_gravkick;
                double pred_vel_bu_z = P[bullet].Vel[2] + P[bullet].GravAccel[2] * dt_gravkick;

                double clight = 299792.458;
                double vbx = pred_vel_bu_x;
                double vby = pred_vel_bu_y;
                double vbz = pred_vel_bu_z;
                double vtx = pred_vel_ta_x;
                double vty = pred_vel_ta_y;
                double vtz = pred_vel_ta_z;
                double vcmx = 0.5*(vbx+vtx);
                double vcmy = 0.5*(vby+vty);
                double vcmz = 0.5*(vbz+vtz);
                double vrp = sqrt((4.0*Delta*clight*clight+m1*vr*vr)/m2);

                #ifdef DYDEBUGKIN
                   printf("############     In the Halo frame    #####################################\n");
                   double v1i=sqrt(pow(P[bullet].Vel[0],2)+pow(P[bullet].Vel[1],2)+pow(P[bullet].Vel[2],2));
                   double v2i=sqrt(pow(P[target].Vel[0],2)+pow(P[target].Vel[1],2)+pow(P[target].Vel[2],2));
                   double vcmi=sqrt(pow(vcmx,2)+pow(vcmy,2)+pow(vcmz,2));
                   double vrinit=vr;
                   printf(">>>>> initial v1=%g v2=%g vr=%g vcmi=%g \n",v1i,v2i,vrinit,vcmi);
                   double vrxi=P[bullet].Vel[0]-P[target].Vel[0];
                   double vryi=P[bullet].Vel[1]-P[target].Vel[1];
                   double vrzi=P[bullet].Vel[2]-P[target].Vel[2];
                   printf("############     In the CM frame    #####################################\n");
                   double v1cmx = vbx - vcmx;
                   double v2cmx = vtx - vcmx;
                   double v1cmy = vby - vcmy;
                   double v2cmy = vty - vcmy;
                   double v1cmz = vbz - vcmz;
                   double v2cmz = vtz - vcmz;
                   printf("v1cmx = -v2cmx = vrx/2: %g, %g, %g  \n",v1cmx,-v2cmx,vrxi/2);
                   printf("v1cmy = -v2cmy = vry/2: %g, %g, %g  \n",v1cmy,-v2cmy,vryi/2);
                   printf("v1cmz = -v2cmz = vrz/2: %g, %g, %g  \n",v1cmz,-v2cmz,vrzi/2);
                #endif
          
                // if inelastic scattering
                //P[bullet].Type=2; 
                //P[target].Type=2; 
                //P[bullet].Mass=m2;
                //P[target].Mass=m2;
          
                P[bullet].Vel[0] = m1/m2*vcmx + 0.5*vrp*vhatbp_x;
                P[bullet].Vel[1] = m1/m2*vcmy + 0.5*vrp*vhatbp_y;
                P[bullet].Vel[2] = m1/m2*vcmz + 0.5*vrp*vhatbp_z;
                P[target].Vel[0] = m1/m2*vcmx - 0.5*vrp*vhatbp_x;
                P[target].Vel[1] = m1/m2*vcmy - 0.5*vrp*vhatbp_y;
                P[target].Vel[2] = m1/m2*vcmz - 0.5*vrp*vhatbp_z;
          
                #ifdef DYDEBUGKIN
                   double v1f=sqrt(pow(P[bullet].Vel[0],2)+pow(P[bullet].Vel[1],2)+pow(P[bullet].Vel[2],2));
                   double v2f=sqrt(pow(P[target].Vel[0],2)+pow(P[target].Vel[1],2)+pow(P[target].Vel[2],2));
                   double vcmf=0.5*sqrt(pow(P[bullet].Vel[0]+P[target].Vel[0],2)+pow(P[bullet].Vel[1]+P[target].Vel[1],2)+pow(P[bullet].Vel[2]+P[target].Vel[2],2));
                   double vrfinal=sqrt(pow(P[bullet].Vel[0]-P[target].Vel[0],2)+pow(P[bullet].Vel[1]-P[target].Vel[1],2)+pow(P[bullet].Vel[2]-P[target].Vel[2],2));
                   printf("############     In the Halo frame    #####################################\n");
                   printf(">>>>> scatted v1=%g v2=%g vr=%g vcm=%g \n",v1f,v2f,vrfinal,vcmf);
                   printf("***** energy   conservation: %g==%g \n",0.5*m1*v1i*v1i+0.5*m1*v2i*v2i+Delta*clight*clight,0.5*m2*v1f*v1f+0.5*m2*v2f*v2f);
                   printf("***** momentum conservation x: %g==%g \n",m1*vcmx,m2*(P[bullet].Vel[0]+P[target].Vel[0])/2);
                   printf("***** momentum conservation y: %g==%g \n",m1*vcmy,m2*(P[bullet].Vel[1]+P[target].Vel[1])/2);
                   printf("***** momentum conservation z: %g==%g \n",m1*vcmz,m2*(P[bullet].Vel[2]+P[target].Vel[2])/2);
          
                   double costhH=(vbx*P[bullet].Vel[0]+vby*P[bullet].Vel[1]+vbz*P[bullet].Vel[2])/(v1i*v1f);
                   printf("***** deflection angle: %g (degree), costhH=%g \n",acos(costhH)/3.141593*180,costhH);
                   printf("############     In the CM frame    #####################################\n");
                   printf("***** deflection angle CM frame: %g (degree), costhp=%g, costhp(reco)=%g \n",acos(costhp),costhp,vhatbp_x*cdz_x+vhatbp_y*cdz_y+vhatbp_z*cdz_z);
                   printf("***** energy   conservation: %g==%g \n",0.5*m1*vr*vr/4+0.5*m1*vr*vr/4+Delta*clight*clight,0.5*m2*vrp*vrp/4+0.5*m2*vrp*vrp/4);
                #endif
          
                P[bullet].tscatt+=1;
                P[target].tscatt+=1;
                ntotscat++;
                return 1;
	    }
	    else // mode1
	    {
            // ZXY 2022.02.04: Down scatter for mode = 1 ......
            //-----------------------------------------------------------------------------------------
            double prob_d = gsl_rng_uniform(random_generator);

              if( (SIDMDataResult[bullet].Type == 1) && (P[target].Type == 1) && prob_d < ratio_ds ){

                  Delta = mass*deltaf; // deltaf the fractional mass difference, defined previously ......
                  m1 = mass;
                  m2 = m1 - Delta;

                  SIDMDataResult[bullet].Type = 2;
                  SIDMDataResult[bullet].Mass = m2;
                //  P[bullet].Type = 2;
                //  P[bullet].Mass = m2;

                  P[target].Type = 2;
                  P[target].Mass = m2;

                }
              else{
                  Delta = 0;
                  m1 = mass;
                  m2 = mass;
              }
             //-----------------------------------------------------------------------------------------
             // double m1 = mass;
             // double Delta = deltaf*P[bullet].Mass;
             // double Delta=0.;
             // double m2=m1-0.5*Delta;
                double clight = 299792.458;

                // predicted exported particle velocity for bullet
                double pred_vel_dr_bu_x = SIDMDataResult[bullet].Vel[0] + GravDataResult[bullet].u.Acc[0] * dt_gravkick;
                double pred_vel_dr_bu_y = SIDMDataResult[bullet].Vel[1] + GravDataResult[bullet].u.Acc[1] * dt_gravkick;
                double pred_vel_dr_bu_z = SIDMDataResult[bullet].Vel[2] + GravDataResult[bullet].u.Acc[2] * dt_gravkick;

                double vbx = pred_vel_dr_bu_x;
                double vby = pred_vel_dr_bu_y;
                double vbz = pred_vel_dr_bu_z;
                double vtx = pred_vel_ta_x;
                double vty = pred_vel_ta_y;
                double vtz = pred_vel_ta_z;
                double vcmx = 0.5*(vbx+vtx);
                double vcmy = 0.5*(vby+vty);
                double vcmz = 0.5*(vbz+vtz);
                double vrp = sqrt((4.0*Delta*clight*clight+m1*vr*vr)/m2);
          
                #ifdef DYDEBUGKIN
                   printf("############     In the Halo frame (Comm Buffer)   ########################\n");
                   double v1i=sqrt(pow(vbx,2)+pow(vby,2)+pow(vbz,2));
                   double v2i=sqrt(pow(P[target].Vel[0],2)+pow(P[target].Vel[1],2)+pow(P[target].Vel[2],2));
                   double vcmi=sqrt(pow(vcmx,2)+pow(vcmy,2)+pow(vcmz,2));
                   double vrinit=vr;
                   printf(">>>>> initial v1=%g v2=%g vr=%g vcmi=%g \n",v1i,v2i,vrinit,vcmi);
                   double vrxi=SIDMDataResult[bullet].Vel[0]-P[target].Vel[0];
                   double vryi=SIDMDataResult[bullet].Vel[1]-P[target].Vel[1];
                   double vrzi=SIDMDataResult[bullet].Vel[2]-P[target].Vel[2];
                #endif
          
                SIDMDataResult[bullet].Vel[0] = m1/m2*vcmx + 0.5*vrp*vhatbp_x;
                SIDMDataResult[bullet].Vel[1] = m1/m2*vcmy + 0.5*vrp*vhatbp_y;
                SIDMDataResult[bullet].Vel[2] = m1/m2*vcmz + 0.5*vrp*vhatbp_z;
                P[target].Vel[0]              = m1/m2*vcmx - 0.5*vrp*vhatbp_x;
                P[target].Vel[1]              = m1/m2*vcmy - 0.5*vrp*vhatbp_y;
                P[target].Vel[2]              = m1/m2*vcmz - 0.5*vrp*vhatbp_z;
          
                #ifdef DYDEBUGKIN
                   double v1f=sqrt(pow(SIDMDataResult[bullet].Vel[0],2)+pow(SIDMDataResult[bullet].Vel[1],2)+pow(SIDMDataResult[bullet].Vel[2],2));
                   double v2f=sqrt(pow(P[target].Vel[0],2)+pow(P[target].Vel[1],2)+pow(P[target].Vel[2],2));
                   double vcmf=0.5*sqrt(pow(SIDMDataResult[bullet].Vel[0]+P[target].Vel[0],2)+pow(SIDMDataResult[bullet].Vel[1]+P[target].Vel[1],2)+pow(SIDMDataResult[bullet].Vel[2]+P[target].Vel[2],2));
                   double vrfinal=sqrt(pow(SIDMDataResult[bullet].Vel[0]-P[target].Vel[0],2)+pow(SIDMDataResult[bullet].Vel[1]-P[target].Vel[1],2)+pow(SIDMDataResult[bullet].Vel[2]-P[target].Vel[2],2));
                   printf(">>>>> scatted v1=%g v2=%g vr=%g vcm=%g \n",v1f,v2f,vrfinal,vcmf);
                   printf("***** energy   conservation: %g==%g \n",0.5*m1*v1i*v1i+0.5*m1*v2i*v2i+Delta*clight*clight,0.5*m2*v1f*v1f+0.5*m2*v2f*v2f);
                   printf("***** momentum conservation x: %g==%g \n",m1*vcmx,m2*(SIDMDataResult[bullet].Vel[0]+P[target].Vel[0])/2);
                   printf("***** momentum conservation y: %g==%g \n",m1*vcmy,m2*(SIDMDataResult[bullet].Vel[1]+P[target].Vel[1])/2);
                   printf("***** momentum conservation z: %g==%g \n",m1*vcmz,m2*(SIDMDataResult[bullet].Vel[2]+P[target].Vel[2])/2);
          
                   double costhH=(vbx*SIDMDataResult[bullet].Vel[0]+vby*SIDMDataResult[bullet].Vel[1]+vbz*SIDMDataResult[bullet].Vel[2])/(v1i*v1f);
                   printf("***** deflection angle: %g (degree) \n",acos(costhH)/3.141593*180);
                   printf("############     In the CM frame (Comm Buffer)   ########################\n");
                   printf("***** energy   conservation: %g==%g \n",0.5*m1*vr*vr/4+0.5*m1*vr*vr/4+Delta*clight*clight,0.5*m2*vrp*vrp/4+0.5*m2*vrp*vrp/4);
                #endif
          
                SIDMDataResult[bullet].tscatt += 1;
                SIDMDataResult[bullet].sendTask = SIDMDataGet[bullet].sendTask;
                SIDMDataResult[bullet].Task = SIDMDataGet[bullet].Task;
                P[target].tscatt += 1;
                ntotscat++;

	    }
	}
    }
    while(startnode >= 0);

    #ifdef DYDEBUG
    if( mode==0&&P[bullet].ID==idcheck ){
        printf("PID==idcheck, mode==%d, task:%d, pos_x=%g, total number of neighbors: %d \n",mode,ThisTask,pos_x,numngb);
    }
    #endif

}

void randomize(int list[], int n) {
    srand(time(0)); // it is changing roughly once a sec
    // printf("time: %d",time(0));
    int i;
    for(i = n-1; i > 0; i--) {
        // swap with elements in its front
	int j = rand() % (i+1);
	int temp = list[i];
	list[i] = list[j];
	list[j] = temp;
    }
}

double get_Wij(double h,double dr){

   double Wij=8.0/(3.14159265358979*pow(h,3));
   double q=dr/h;
   if(dr < h/2.) Wij*=(1-6*q*q+6*pow(q,3));
   else if(dr<h) Wij*=(2*pow(1-q,3));
   else Wij=0;

   return Wij;
}

double get_gij(double h,double dr){

   double gij;
   gij=0;
   // dr should always be larger than 0
   if(dr < h/2.){
      gij=(2*(4160*pow(dr,9) - 11520*pow(dr,8)*h + 8064*pow(dr,7)*pow(h,2) +
       2688*pow(dr,6)*pow(h,3) - 2016*pow(dr,5)*pow(h,4) - 5040*pow(dr,4)*pow(h,5) +
       3864*pow(dr,3)*pow(h,6) + 1116*pow(dr,2)*pow(h,7) - 1854*dr*pow(h,8) +
       491*pow(h,9)))/(315.*pow(h,12)*3.14159265358979);
   }
   else if(dr<h){
      gij=(-128*pow(dr - h,6)*(5*pow(dr,3) - 6*pow(dr,2)*h - 3*dr*pow(h,2) - 3*pow(h,3)))/
   (105.*pow(h,12)*3.14159265358979);
   }
   else gij=0;
   // normalize gij: int d^3 x gij(|x|) = 1
   return 12.0974*gij;
}

// costhp distribution of Rutherford scatterings
double fcosthR(double costh,double v,double vw){
    return (2*pow(vw,2)*(pow(v,2) + pow(vw,2)))/
   pow((-1 + costh)*pow(v,2) - 2*pow(vw,2),2);
}

// costhp distribution of Moller scatterings
double fcosthM(double costh,double v,double vw){

    return ((1 + 3*pow(costh,2))*pow(v,4) + 4*pow(v,2)*pow(vw,2) + 4*pow(vw,4))/
   (pow((-1 + pow(costh,2))*pow(v,4) - 4*pow(v,2)*pow(vw,2) - 4*pow(vw,4),2)*(1/(pow(v,2)*pow(vw,2) + pow(vw,4)) + log(pow(vw,2)/(pow(v,2) + pow(vw,2)))/(pow(v,4) + 2*pow(v,2)*pow(vw,2))));

}
