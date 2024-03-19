#include "pluto.h"
#include "local_pluto.h"

/* *************************************************************** */
void ComputeUserVar (const Data *d, Grid *grid)
/*
 *
 *  PURPOSE
 *
 *    Define user-defined output variables
 *
 *
 *
 ***************************************************************** */
{
  int i,j,k;
  double ***T, ***n, ***p, ***mach;
  double ***prs, ***rho;
  double *dx, *dy, *dz;
  #if COOLING==NO || COOLING==TABULATED || COOLING==TOWNSEND
  double dummy[4];
  double mu = MeanMolecularWeight((double*)d->Vc, dummy);
  #else
  double mu = MeanMolecularWeight((double*)d->Vc);
  #endif
  T    = GetUserVar("Temp");
  n    = GetUserVar("ndens");
  p    = GetUserVar("PbykB");
  mach = GetUserVar("mach");
  rho = d->Vc[RHO];    // pointer shortcut to density
  prs = d->Vc[PRS];    // pointer shortcut to pressure
  dx = grid->dx[IDIR]; // shortcut to dx 
  dy = grid->dx[JDIR]; // shortcut to dy
  dz = grid->dx[KDIR]; // shortcut to dy 
  TOT_LOOP(k,j,i){
   T[k][j][i] = (d->Vc[PRS][k][j][i]/d->Vc[RHO][k][j][i])*pow(UNIT_VELOCITY,2)*(CONST_mp*mu)/CONST_kB;
   n[k][j][i] = d->Vc[RHO][k][j][i] * UNIT_DENSITY / (mu * CONST_mp);
   p[k][j][i] = d->Vc[PRS][k][j][i] * UNIT_DENSITY * pow(UNIT_VELOCITY,2.) / CONST_kB;
   mach[k][j][i] = sqrt( DIM_EXPAND(d->Vc[iVR][k][j][i]*d->Vc[iVR][k][j][i], 
                                        + d->Vc[iVZ][k][j][i]*d->Vc[iVZ][k][j][i], 
                                        + d->Vc[iVPHI][k][j][i]*d->Vc[iVPHI][k][j][i]) )/sqrt(g_gamma*(d->Vc[PRS][k][j][i]/d->Vc[RHO][k][j][i]));
  }
}
/* ************************************************************* */
void ChangeOutputVar ()
/* 
 *
 * 
 *************************************************************** */
{ 
    
/*
  Image *image;

  //SetOutputVar("bx1", FLT_OUTPUT, NO);
  SetOutputVar("rho", PNG_OUTPUT, YES);
  SetOutputVar("T", PNG_OUTPUT, YES);
  SetOutputVar("tr1", PNG_OUTPUT, YES);
  SetOutputVar("prs", PNG_OUTPUT, YES);
  //SetOutputVar("vortz", PNG_OUTPUT, YES);

  image = GetImage ("rho");
  image->slice_plane = X12_PLANE;
  image->slice_coord = 0.;
  //image->max = image->min = 0.0;
  image->logscale = 1;
  image->colormap = "red";

  image = GetImage ("prs");
  image->slice_plane = X12_PLANE;
  image->slice_coord = 0.;
  //image->max = image->min = 0.0;
  image->logscale = 1;
  image->colormap = "red";

  image = GetImage ("T");
  image->slice_plane = X12_PLANE;
  image->slice_coord = 0.;
  //image->max = image->min = 0.0;
  image->logscale = 1;
  image->colormap = "red";

  image = GetImage ("tr1");
  image->slice_plane = X12_PLANE;
  image->slice_coord = 0.;
  //image->max = image->min = 0.0;
  image->logscale = 1;
  image->colormap = "red";
  */

#ifdef PARTICLES
  //SetOutputVar ("energy",PARTICLES_FLT_OUTPUT, NO);
//  SetOutputVar ("x1",    PARTICLES_FLT_OUTPUT, NO);
  //SetOutputVar ("vx1",   PARTICLES_FLT_OUTPUT, NO);
#endif

}





