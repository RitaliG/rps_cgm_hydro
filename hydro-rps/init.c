/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Contains basic functions for problem initialization.

  The init.c file collects most of the user-supplied functions useful 
  for problem configuration.
  It is automatically searched for by the makefile.

  \author A. Mignone (mignone@ph.unito.it)
  \date   March 5, 2017
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"
#include "local_pluto.h"

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*! 
 * The Init() function can be used to assign initial conditions as
 * as a function of spatial position.
 *
 * \param [out] v   a pointer to a vector of primitive variables
 * \param [in] x1   coordinate point in the 1st dimension
 * \param [in] x2   coordinate point in the 2nd dimension
 * \param [in] x3   coordinate point in the 3rdt dimension
 *
 * The meaning of x1, x2 and x3 depends on the geometry:
 * \f[ \begin{array}{cccl}
 *    x_1  & x_2    & x_3  & \mathrm{Geometry}    \\ \noalign{\medskip}
 *     \hline
 *    x    &   y    &  z   & \mathrm{Cartesian}   \\ \noalign{\medskip}
 *    R    &   z    &  -   & \mathrm{cylindrical} \\ \noalign{\medskip}
 *    R    & \phi   &  z   & \mathrm{polar}       \\ \noalign{\medskip}
 *    r    & \theta & \phi & \mathrm{spherical} 
 *    \end{array}
 *  \f]
 *
 * Variable names are accessed by means of an index v[nv], where
 * nv = RHO is density, nv = PRS is pressure, nv = (VX1, VX2, VX3) are
 * the three components of velocity, and so forth.
 *
 *********************************************************************** */
{
}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*! 
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 *
 *
 *********************************************************************** */
{
  g_gamma      = 5/3.;
  double oth_mu[4];
  double mu    = MeanMolecularWeight((double*)d->Vc, oth_mu);
  printLog("mu = %.3f\tmue = %.3f\tmui = %.3f\n",mu, oth_mu[0], oth_mu[1]);

  g_minCoolingTemp = g_inputParam[T_FLOOR];
  g_rhoICM         = g_inputParam[N_ICM] * mu * CONST_mp / UNIT_DENSITY; // density (code units)

  /* set pointershortcuts */
  int i, j, k;
  double *x1 = grid->x[IDIR];
  double *x2 = grid->x[JDIR];
  double *x3 = grid->x[KDIR];

  double x, y, z, R;
  double z_offset = g_inputParam[OFFSET];

   double theta    = g_inputParam[INCL] * CONST_PI / 180.; // angle (in degrees -> radians) between direction of wind and the normal to the galactic disk
  double vx1_prime=0., vx2_prime=0., vx3_prime=0.;

  /* read from dbl file and interpolate for simulation domain */
  int id_rho, id_prs, id_vx1, id_vx2, id_vx3;
  int offset, size[3];

  #if WIND != NO && WIND_DOMAIN != NO
  double vWind = g_inputParam[V_WIND] * 1.e5 / UNIT_VELOCITY;  // (km/s -> code units)
  double v_smooth;
  #endif

  double Temp;
  double radCGM    = g_inputParam[R_CGM];
  double sigmaRCGM = g_inputParam[SIGMA_R_CGM];
  
  // over pressurised ICM around
  double rhoICM  = g_rhoICM; //g_inputParam[N_ICM] * mu * CONST_mp / UNIT_DENSITY; // density (code units)
  double TWind   = g_inputParam[T_ICM];
  double pICM    = ( (rhoICM*UNIT_DENSITY) * CONST_kB * TWind / (mu*CONST_mp) ) / (UNIT_DENSITY*pow(UNIT_VELOCITY,2));

  /* Interpolate density (1st record in the file) */
  id_rho = InputDataOpen("ic_file.dbl", "grid0.out", " ", 0, CENTER);
  TOT_LOOP(k,j,i){
    d->Vc[RHO][k][j][i]   = InputDataInterpolate(id_rho, x1[i], x2[j], x3[k]);
    if ((x3[k] - g_domBeg[KDIR])/g_domBeg[KDIR] < 1e-3) g_rhoICM = d->Vc[RHO][k][j][i];   
  }
 
  InputDataGridSize(id_rho, size);   // get grid size
  offset = DIM_EXPAND((long)size[0], *(long)size[1], *(long)size[2];);
  InputDataClose(id_rho);

  /* Interpolate pressure (5th record in the file) */
  id_prs = InputDataOpen("ic_file.dbl", "grid0.out", " ", 4*offset, CENTER);
  double pCGM;
  double pIsoBar = 1e24; // isobaric pressure for CGM and ISM 
  int count =0;
  
  TOT_LOOP(k,j,i){
    d->Vc[PRS][k][j][i]   = InputDataInterpolate(id_prs, x1[i], x2[j], x3[k]);
    d->Vc[TRC][k][j][i]   = 0.;
    d->Vc[TRC+1][k][j][i] = 0.;
    d->Vc[TRC+2][k][j][i] = 0.;
  }
  InputDataClose(id_prs);

 
  /* Interpolate velocities (2nd, 3rd and 4th record in the file) */
  DIM_EXPAND(
    id_vx1 = InputDataOpen("ic_file.dbl", "grid0.out", " ", 1*offset, CENTER);,
    id_vx2 = InputDataOpen("ic_file.dbl", "grid0.out", " ", 2*offset, CENTER);,
    // id_vx3 = InputDataOpen("ic_file.dbl", "gridRB06.out", " ", 3*offset, CENTER);
    )

  TOT_LOOP(k,j,i){
    /* convert to cartesian axes ------------------------------- */
    x = x1[i] * cos(x2[j]);
    y = x1[i] * sin(x2[j]);
    R = x1[i] ; // galactocentric radius
  
    /* get values from galaxy's frame of reference ------------- */
    DIM_EXPAND(
      vx1_prime = InputDataInterpolate(id_vx1, x1[i], x2[j], x3[k]);,
      vx2_prime = InputDataInterpolate(id_vx2, x1[i], x2[j], x3[k]);,
      // vx3_prime = InputDataInterpolate(id_vx3, x1[i], x2[j], x3[k]);
      )      

    /* convert into polar components --------------------------- */
    DIM_EXPAND(
      d->Vc[iVR][k][j][i]   = ( vx1_prime*x + vx2_prime*cos(theta)*y )/R ;,
      d->Vc[iVPHI][k][j][i] = ( vx2_prime*cos(theta)*x - vx1_prime*y )/R ;,
      d->Vc[iVZ][k][j][i]   =   vx2_prime * sin(theta) ;  // vx3_prime is zero for our setup initially
    )
    
    Temp = (d->Vc[PRS][k][j][i]/d->Vc[RHO][k][j][i])*pow(UNIT_VELOCITY,2)*(CONST_mp*mu)/CONST_kB; // Kelvin
    if (Temp>=TWind/3.) {
      #if WIND != NO && WIND_DOMAIN != NO
        v_smooth = 0.5 * (1 + tanh(( sqrt(pow(x1[i], 2.) + pow(x3[k]-z_offset, 2.)) - radCGM)/sigmaRCGM));
        d->Vc[iVZ][k][j][i] = vWind  * v_smooth;
      #endif
      d->Vc[RHO][k][j][i]   = rhoICM * v_smooth;
      d->Vc[PRS][k][j][i]   = pICM   * v_smooth;
    }
   
    d->Vc[TRC][k][j][i]   = ((Temp <= 3*g_inputParam[T_ISM]) ? 1.0 : 0. );
    d->Vc[TRC+1][k][j][i] = ((Temp > 3*g_inputParam[T_ISM] && Temp <= 3*g_inputParam[T_CGM]) ? 1.0 : 0.);
    d->Vc[TRC+2][k][j][i] = 1. - (d->Vc[TRC][k][j][i]+d->Vc[TRC+1][k][j][i]);
  }
  
  DIM_EXPAND(
    InputDataClose(id_vx2);,
    InputDataClose(id_vx1);,
    // InputDataClose(id_vx3);
  )
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/*! 
 *  Perform runtime data analysis.
 *
 * \param [in] d the PLUTO Data structure
 * \param [in] grid   pointer to array of Grid structures  
 *
 *********************************************************************** */
{
  int i, j, k;
  double *x1, *x2, *x3;
  double *dx1, *dx2, *dx3;
  double dV, rad=0., hei=0.;

  x1  = grid->x[IDIR];   x2 = grid->x[JDIR];   x3 = grid->x[KDIR];
  dx1 = grid->dx[IDIR]; dx2 = grid->dx[JDIR]; dx3 = grid->dx[KDIR];
  
  double dummy[4];
  double mu = MeanMolecularWeight((double*)d->Vc, dummy);
  double UNIT_MASS = (UNIT_DENSITY * pow(UNIT_LENGTH, 3.));

  double z_offset = g_inputParam[OFFSET];
  double theta    = g_inputParam[INCL] * CONST_PI / 180.; // angle (in degrees -> radians) between direction of wind and the normal to the galactic disk

  double x, y, z;
  double x_prime, y_prime, z_prime;

  double Temp;
  double vWind     = g_inputParam[V_WIND] * (1.e5 / UNIT_VELOCITY);  // (km/s -> code units)
  double TWind     = g_inputParam[T_ICM]; // Kelvin
  double TDisk     = g_inputParam[T_ISM]; // Kelvin
  double rad_limit = g_inputParam[R_DISK_CYL];
  double z_limit   = g_inputParam[Z_DISK_CYL]; // limit of cylinder to be considered (code units)

  /* --- variables used --- */
  double rho;
  double vx1, vx2, vx3;
  double vr, vphi, vz;
  double dense_mass = 0., dense_vol  = 0.;
  double disk_mass  = 0., disk_vol   = 0.;
  double tail_mass  = 0., tail_vol   = 0.;
  double cgm_mass   = 0., cgm_vol    = 0.;
  double cm_disk    = 0., cm_cgm     = 0.; // center of mass : component along vertical
  double vr_disk    = 0., vphi_disk  = 0.,    vz_disk = 0.;
  double vr_cgm     = 0., vphi_cgm   = 0.,    vz_cgm  = 0.;

  double mass_tot   = 0., vol_tot    = 0., mass_tot2 = 0.;
  double disk_tracer= 0., cgm_tracer = 0., wind_tracer= 0.; 
  double tail_mass_cold  = 0., tail_vol_cold   = 0.;
  double disk_mass_cold  = 0., disk_vol_cold   = 0.;

  // double vxTr = 0., vyTr = 0., vzTr = 0., vxWind = 0., vyWind = 0., vzWind = 0.;

  /* Check if restart of the code was done */
  static int restart = 1;
  static long int nstep = -1;
  static double tanl = 0.;

  /* ----- Main Loop ----- */
  if (g_stepNumber == 0) restart = 0;

  if (restart == 1){ /* means we have restarted */
    restart = 0;
    FILE *fp;
    char fname[512];
    int dummy;
    sprintf (fname, "%s/restart-analysis.out",RuntimeGet()->output_dir);
    fp = fopen(fname,"r");
    dummy = fscanf(fp, "%ld", &nstep);
    dummy = fscanf(fp, "%lf", &tanl);
    fclose(fp);
    //printLog("This is step %ld!\n", g_stepNumber);
    //printLog("Analysis should resume from step %d:\n", nstep);
    //printLog("Initial Tracer read as %e\n", iniDiskTrc);
  }
  
  if (g_stepNumber<=nstep && g_time<=(tanl+0.5*g_anldt)) return;

  /* ----- Main Loop ----- */
  DOM_LOOP(k,j,i){
    dV  = grid->dV[k][j][i];
    rho = d->Vc[RHO][k][j][i];
    vx1 = d->Vc[VX1][k][j][i];  // x-velocity 
    vx2 = d->Vc[VX2][k][j][i];  // y-velocity
    vx3 = d->Vc[VX3][k][j][i];  // z-velocity 

    #if GEOMETRY == POLAR
    /* convert POLAR (rad, phi, z) -> CARTESIAN (x, y, z) in lab frame ---- */
    x =  x1[i] * cos(x2[j]);
    y =  x1[i] * sin(x2[j]);
    z =  x3[k] - z_offset;
    /* rotated coordinates in galaxy's frame of reference ----------------- */
    x_prime =  x;
    y_prime =  cos(theta) * y + sin(theta) * z;
    z_prime = -sin(theta) * y + cos(theta) * z;
    rad     =  sqrt(pow(x_prime, 2.) + pow(y_prime, 2.));
    hei     =  z_prime;
    vr = vx1; vphi = vx2; vz = vx3;

    #elif GEOMETRY == CARTESIAN
    rad  = sqrt( pow(x1[i] - z_offset, 2.) + pow(x2[j] - z_offset, 2.) );
    hei  = x3[k] - z_offset;
    #endif

    Temp = (d->Vc[PRS][k][j][i]/d->Vc[RHO][k][j][i])*pow(UNIT_VELOCITY,2)*(CONST_mp*mu)/CONST_kB; // Kelvin
    mass_tot  += rho * dV;
    // cold mass evolution : temperature based : all volume
    if (Temp <= 3*g_inputParam[T_ISM]){
      dense_mass += rho * dV; // temperature based (all volume)
      dense_vol  += dV;
    }
    // warm mass evolution : temperature based : all volume
    else if (Temp>3*g_inputParam[T_ISM] && Temp<=3*g_inputParam[T_CGM]){
       cgm_mass += rho * dV;
       cgm_vol  += dV; 
    }
    
    // disk mass evolution : tracer based
    if (fabs(hei)<z_limit  && rad<=rad_limit){
      disk_mass += rho * d->Vc[TRC][k][j][i] *  dV;
      disk_vol  += d->Vc[TRC+1][k][j][i] * dV;
      if (Temp <= 3*g_inputParam[T_ISM]){
        disk_mass_cold += rho * dV;
        disk_vol_cold  += dV;
      }
    }
    // tail mass increase from stripping : tracer based
    if (hei >= z_limit && hei<80 && rad<=1.5*rad_limit){
      tail_mass += rho * d->Vc[TRC+1][k][j][i] * dV;
      tail_vol  += d->Vc[TRC+1][k][j][i] * dV;
      // how much of that is cold
      if (Temp <= 3*g_inputParam[T_ISM]){
        tail_mass_cold += rho * dV;
        tail_vol_cold  += dV;
      }
    }
    
    // check on tracer leaking out
    disk_tracer += rho*d->Vc[TRC][k][j][i]*dV;
    cgm_tracer  += rho*d->Vc[TRC+1][k][j][i]*dV;
    wind_tracer += rho*d->Vc[TRC+2][k][j][i]*dV;

    // center of mass of various components
    // for disk gas
    cm_disk =  z*rho*d->Vc[TRC][k][j][i]*dV; // along vertical 
    vr_disk = vr*rho*d->Vc[TRC][k][j][i]*dV;
    vz_disk = vz*rho*d->Vc[TRC][k][j][i]*dV;
    // for cgm gas
    cm_cgm  =  z * (rho*d->Vc[TRC+1][k][j][i]) * dV ; // along vertical
    vr_cgm  = vr * (rho*d->Vc[TRC+1][k][j][i]) * dV;
    vz_cgm  = vz * (rho*d->Vc[TRC+1][k][j][i]) * dV;
  }

  /* ------ Parallel data reduction ------ */
  #ifdef PARALLEL
  int count = 22;
  int a = 0;
  double sendArray[count], recvArray[count];

  sendArray[a++] = mass_tot;
  sendArray[a++] = dense_mass;      sendArray[a++] = dense_vol;
  sendArray[a++] = cgm_mass;        sendArray[a++] = cgm_vol;
  sendArray[a++] = disk_mass;       sendArray[a++] = disk_vol;

  sendArray[a++] = disk_mass_cold;  sendArray[a++] = disk_vol_cold;
  sendArray[a++] = tail_mass;       sendArray[a++] = tail_vol;
  sendArray[a++] = tail_mass_cold;  sendArray[a++] = tail_vol_cold;
  sendArray[a++] = disk_tracer;     sendArray[a++] = cgm_tracer; 

  sendArray[a++] = wind_tracer;
  sendArray[a++] = cm_disk;         sendArray[a++] = cm_cgm; 
  sendArray[a++] = vr_disk;         sendArray[a++] = vr_cgm; 
  sendArray[a++] = vz_disk;         sendArray[a++] = vz_cgm;  

  MPI_Reduce(sendArray, recvArray, count, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);

  if (prank == 0){
     a = 0;
     mass_tot       = recvArray[a++];
     dense_mass     = recvArray[a++];   dense_vol     = recvArray[a++];
     cgm_mass       = recvArray[a++];   cgm_vol       = recvArray[a++];
     disk_mass      = recvArray[a++];   disk_vol      = recvArray[a++];

     disk_mass_cold = recvArray[a++];   disk_vol_cold = recvArray[a++];
     tail_mass      = recvArray[a++];   tail_vol      = recvArray[a++];
     tail_mass_cold = recvArray[a++];   tail_vol_cold = recvArray[a++];
     disk_tracer    = recvArray[a++];   cgm_tracer    = recvArray[a++];

     wind_tracer    = recvArray[a++];
     cm_disk        = recvArray[a++];   cm_cgm        = recvArray[a++];
     vr_disk        = recvArray[a++];   vr_cgm        = recvArray[a++];
     vz_disk        = recvArray[a++];   vz_cgm        = recvArray[a++];
  }
  #endif
  /* --- end of parallel data reduction --- */

  /* --- Write ascii file "analysis.dat" to disk --- */
  if (prank == 0){
     cm_disk /= disk_tracer;  vr_disk /= disk_tracer; vz_disk /= disk_tracer;
     cm_cgm  /= cgm_tracer;   vr_cgm  /= cgm_tracer;  vz_cgm  /= cgm_tracer; 
     
     FILE *fp;
     char fname[512];
     static double tpos = -1.0;
     sprintf (fname, "%s/analysis.dat", RuntimeGet()->output_dir);
     if (g_stepNumber == 0){    /* Open for writing only when we're starting */
     fp = fopen(fname,"w");     /* from beginning */
     fprintf (fp,"# %s\t\t%s\t\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s \n",
                     "(0) t",            "(1) dt",        "(2) MTot [M.]",     
                     "(3) MDense [M.]",  "(4) denVol",
                     "(5) MCgm   [M.]",  "(6) cgmVol", 
                     "(7) MDisk  [M.]",  "(8) diskVol",

                     "(9) MDiskC [M.]",  "(10)diskVolC",
                     "(11)MTail  [M.]",  "(12)tailVol",
                     "(13)MTailC [M.]",  "(14)tailVolC", 
                     "(15)MDiskTr[M.]",  "(16)MCgmTr[M.]",

                     "(17)MWindTr[M.]",  
                     "(18)Z_Disk[kpc]",  "(19)Z_Cgm [kpc]",
                     "(20)vrDisk[kpc]",  "(21)vrCgm [kpc]",
                     "(22)vzDisk[kpc]",  "(23)vzCgm [kpc]"
             );
     }
     else{                        /* write if current time t>0*/
     if (tpos < 0.0){           /* obtain time coordinate of to last written row*/
        char   sline[512];
        fp = fopen(fname,"r");
        while (fgets(sline, 512, fp))  {}
        sscanf(sline, "%lf\n",&tpos);  /* tpos = time of the last written row */
        fclose(fp);
     }
     fp = fopen(fname,"a");
    }
    if (g_time > tpos){
       fprintf (fp, "%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e \n",
         g_time,                                g_dt,           (mass_tot*UNIT_MASS/CONST_Msun),
         (dense_mass*UNIT_MASS/CONST_Msun),     dense_vol,
         (cgm_mass*UNIT_MASS/CONST_Msun),       cgm_vol,
         (disk_mass*UNIT_MASS/CONST_Msun),      disk_vol, 

         (disk_mass_cold*UNIT_MASS/CONST_Msun), disk_vol_cold, 
         (tail_mass*UNIT_MASS/CONST_Msun),      tail_vol, 
         (tail_mass_cold*UNIT_MASS/CONST_Msun), tail_vol_cold, 
         (disk_tracer*UNIT_MASS/CONST_Msun),    (cgm_tracer*UNIT_MASS/CONST_Msun), 

         (wind_tracer*UNIT_MASS/CONST_Msun),
         (cm_disk*UNIT_LENGTH/(1.e3*CONST_pc)), (cm_cgm*UNIT_LENGTH/(1.e3*CONST_pc)),
         (vr_disk*UNIT_VELOCITY/1.e5),          (vr_cgm*UNIT_VELOCITY/1.e5),
         (vz_disk*UNIT_VELOCITY/1.e5),          (vz_cgm*UNIT_VELOCITY/1.e5)
       );
    }
    fclose(fp);

    /* Write restart file */
    //printLog("Step %d: Writing Analysis restart!\n", g_stepNumber);
    FILE *frestart;
    sprintf (fname, "%s/restart-analysis.out", RuntimeGet()->output_dir);
    frestart = fopen(fname,"w");
    fprintf(frestart,"%ld\n", g_stepNumber);
    fprintf(frestart,"%lf\n", g_time);
    fclose(frestart);
  }
  /* --- end of writing "analyis.dat" file --- */
}
#if PHYSICS == MHD
/* ********************************************************************* */
void BackgroundField (double x1, double x2, double x3, double *B0)
/*!
 * Define the component of a static, curl-free background 
 * magnetic field.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * \param [out] B0 array containing the vector componens of the background
 *                 magnetic field
 *********************************************************************** */
{
   B0[0] = 0.0;
   B0[1] = 0.0;
   B0[2] = 0.0;
}
#endif

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*! 
 *  Assign user-defined boundary conditions.
 *
 * \param [in,out] d  pointer to the PLUTO data structure containing
 *                    cell-centered primitive quantities (d->Vc) and 
 *                    staggered magnetic fields (d->Vs, when used) to 
 *                    be filled.
 * \param [in] box    pointer to a RBox structure containing the lower
 *                    and upper indices of the ghost zone-centers/nodes
 *                    or edges at which data values should be assigned.
 * \param [in] side   specifies the boundary side where ghost zones need
 *                    to be filled. It can assume the following 
 *                    pre-definite values: X1_BEG, X1_END,
 *                                         X2_BEG, X2_END, 
 *                                         X3_BEG, X3_END.
 *                    The special value side == 0 is used to control
 *                    a region inside the computational domain.
 * \param [in] grid  pointer to an array of Grid structures.
 *
 *********************************************************************** */
{
  int   i, j, k, nv;
  double  *x1, *x2, *x3;

  x1 = grid->x[IDIR];
  x2 = grid->x[JDIR];
  x3 = grid->x[KDIR];
  double Temp;

  double vWind = g_inputParam[V_WIND]*1.e5/UNIT_VELOCITY;  // (km/s -> code units)
  double TWind = g_inputParam[T_ICM];
  double dummy[4];
  double mu    = MeanMolecularWeight((double*)d->Vc, dummy);
  double rhoICM  = g_rhoICM; //g_inputParam[N_ICM] * mu * CONST_mp / UNIT_DENSITY; // density (code units)
  double pICM    = ( (rhoICM*UNIT_DENSITY) * CONST_kB * TWind / (mu*CONST_mp) ) / (UNIT_DENSITY*pow(UNIT_VELOCITY,2));

  //if (side == 0) {    /* -- check solution inside domain -- */
  /* TOT_LOOP(k,j,i){
      Temp = (d->Vc[PRS][k][j][i]/d->Vc[RHO][k][j][i])*pow(UNIT_VELOCITY,2)*(CONST_mp*mu)/CONST_kB; // Kelvin
      if(Temp<g_minCoolingTemp) d->Vc[PRS][k][j][i] = d->Vc[RHO][k][j][i] * g_minCoolingTemp / (KELVIN*mu);
    }
  }
  */
  int l; double rad, hei;
  double rad_limit = g_inputParam[R_DISK_CYL], z_limit = g_inputParam[Z_DISK_CYL]; // limit of cylinder to be considered (code units)
  // double offset = g_inputParam[OFFSET];
  double z_offset = g_inputParam[OFFSET];
  double theta    = g_inputParam[INCL] * CONST_PI / 180.; // angle (in degrees -> radians) between direction of wind and the normal to the galactic disk

  double x, y, z;
  double x_prime, y_prime, z_prime;

  if (side == X1_BEG){  /* -- X1_BEG boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X1_END){  /* -- X1_END boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X2_BEG){  /* -- X2_BEG boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X2_END){  /* -- X2_END boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X3_BEG){  /* -- X3_BEG boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){
        #if WIND != NO   /* --- wind conditions to be decided --- */
          #if WIND_DOMAIN != NO  /* --- outflow in the boundary --- */
          d->Vc[RHO][k][j][i] = d->Vc[RHO][KBEG][j][i];;
          d->Vc[PRS][k][j][i] = d->Vc[PRS][KBEG][j][i];
          DIM_EXPAND(
            d->Vc[VX1][k][j][i] = d->Vc[VX1][KBEG][j][i];,
            d->Vc[VX2][k][j][i] = d->Vc[VX2][KBEG][j][i];,
            d->Vc[VX3][k][j][i] = d->Vc[VX3][KBEG][j][i];
          )
          #else   /* --- put wind from lower boundary --- */
          d->Vc[RHO][k][j][i] = rhoICM;
          d->Vc[PRS][k][j][i] = pICM;
          DIM_EXPAND(
            d->Vc[VX1][k][j][i] = d->Vc[VX1][KBEG][j][i];,
            d->Vc[VX2][k][j][i] = d->Vc[VX2][KBEG][j][i];,
            d->Vc[VX3][k][j][i] = vWind;
          )
          #endif
        #else     /* --- no wind => outflow boundary --- */
        d->Vc[RHO][k][j][i] = d->Vc[RHO][KBEG][j][i];;
        d->Vc[PRS][k][j][i] = d->Vc[PRS][KBEG][j][i];
        DIM_EXPAND(
          d->Vc[VX1][k][j][i] = d->Vc[VX1][KBEG][j][i];,
          d->Vc[VX2][k][j][i] = d->Vc[VX2][KBEG][j][i];,
          d->Vc[VX3][k][j][i] = d->Vc[VX3][KBEG][j][i];
        )
        #endif     /* --- end of wind conditions to be decide --- */
        d->Vc[TRC][k][j][i]   = 0.;
        d->Vc[TRC+1][k][j][i] = 0.;
        d->Vc[TRC+2][k][j][i] = 1.;
      }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }

  if (side == X3_END){  /* -- X3_END boundary -- */
    if (box->vpos == CENTER) {
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X1FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X2FACE){
      BOX_LOOP(box,k,j,i){  }
    }else if (box->vpos == X3FACE){
      BOX_LOOP(box,k,j,i){  }
    }
  }
}

#if BODY_FORCE != NO
/* ********************************************************************* */
void BodyForceVector(double *v, double *g, double x1, double x2, double x3)
/*!
 * Prescribe the acceleration vector as a function of the coordinates
 * and the vector of primitive variables *v.
 *
 * \param [in] v  pointer to a cell-centered vector of primitive 
 *                variables
 * \param [out] g acceleration vector
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 *
 *********************************************************************** */
{
  g[IDIR] = 0.0;
  g[JDIR] = 0.0;
  g[KDIR] = 0.0;
}
/* ********************************************************************* */
double BodyForcePotential(double x1, double x2, double x3)
/*!
 * Return the gravitational potential as function of the coordinates.
 *
 * \param [in] x1  position in the 1st coordinate direction \f$x_1\f$
 * \param [in] x2  position in the 2nd coordinate direction \f$x_2\f$
 * \param [in] x3  position in the 3rd coordinate direction \f$x_3\f$
 * 
 * \return The body force potential \f$ \Phi(x_1,x_2,x_3) \f$.
 *
 *********************************************************************** */
{
  double x, y, z;
  double x_prime, y_prime, z_prime, rad;
  double z_offset = g_inputParam[OFFSET];
  double offset   = g_inputParam[OFFSET];
  double theta    = g_inputParam[INCL] * CONST_PI / 180.; // angle (in degrees -> radians) between direction of wind and the normal to the galactic disk

  #if GEOMETRY == CARTESIAN
  double R = sqrt(DIM_EXPAND((x1)*(x1), + (x2)*(x2), ));  // galactocentric radius
  return phi_modNFW(R, x3-offset) + phi_MN(R, x3-offset); // + phi_Bl(R, x3-offset));   // (code units)
  // return (phi_MB(R, x3-offset) + phi_MN(R, x3-offset) + phi_Bl(R, x3-offset));   // (code units)
  
  #elif GEOMETRY == POLAR
  /* convert POLAR (rad, phi, z) -> CARTESIAN (x, y, z) in lab frame ---- */
  x =  x1 * cos(x2);
  y =  x1 * sin(x2);
  z =  x3 - z_offset;
  /* rotated coordinates in galaxy's frame of reference ----------------- */
  x_prime =  x;
  y_prime =  cos(theta) * y + sin(theta) * z;
  z_prime = -sin(theta) * y + cos(theta) * z;
  rad     =  sqrt(pow(x_prime, 2.) + pow(y_prime, 2.));
  return  phi_modNFW(rad, z_prime) + phi_MN(rad, z_prime); // + phi_Bl(rad, z_prime));  // (code units)
  // return (phi_MB(rad, z_prime) + phi_MN(rad, z_prime) + phi_Bl(rad, z_prime));   // (code units)
  
  #endif
}
#endif // end of [if BODY_FORCE != NO]
