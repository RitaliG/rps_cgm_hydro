/* ******************************************************************** */
/*  Various functions for isothermal components of Disk-Halo system 
 *  *  @author : ritali ghosh 
 * 
 *  References:
 *  1) Sarkar K.C., Nath B. B., Sharma P., Shchekinov Y., 2015 MNRAS 
 */
/* ******************************************************************** */
#include "pluto.h"
#include "local_pluto.h"


double phi_MN(double R_val, double Z_val){
    /* Returns Miyamoto-Nagai disc potential (ref: Miyamoto & Nagai (1975)) (code units) 
     * \params:  R_val - galactocentric radius (code units) 
     *           Z_val - height above the disc (code units) 
     */
    double MDisk    = g_inputParam[M_STAR];  // solar units
    double aDisk    = g_inputParam[A_STAR] * 1.e3 * CONST_pc / UNIT_LENGTH;  // (code units)
    double bDisk    = g_inputParam[B_STAR] * 1.e3 * CONST_pc / UNIT_LENGTH;  // (code units)
    
    double factor   = - CONST_G * MDisk * CONST_Msun / UNIT_LENGTH; // cgs
    double dist2    = pow(R_val, 2.) + pow(aDisk + sqrt(pow(Z_val, 2.) + pow(bDisk,2.)), 2.);
    double phiMN    = 1. / sqrt(dist2);
    phiMN          *= (factor / pow(UNIT_VELOCITY, 2.)); // converted to code units
    return phiMN;
}


double phi_MB(double R_val, double Z_val){
    /* Returns dark matter potential from Mori and Burkert model 2000 (code units)
     * \params:  R_val - galactocentric radius (code units) 
     *           Z_val - height above the disc (code units) 
     */
    double rDM      = g_inputParam[R_DM];    // (code units)
    double MrDM     = g_inputParam[M_RDM];   // solar units
    double rhoDM    = 3.e-24 * pow( (rDM*UNIT_LENGTH/(1.e3*CONST_pc)), -2/3.) / UNIT_DENSITY; //g_inputParam[RHO_DM];  // (code units)
    double factor   = -CONST_PI * (CONST_G*pow(UNIT_LENGTH,2)*UNIT_DENSITY/pow(UNIT_VELOCITY, 2.)) * rhoDM * pow(rDM, 2.); // code units
    double radDist  = sqrt( pow(R_val, 2.) + pow(Z_val, 2.) );
    double x        = (radDist/rDM ); // ratio
    double eps      = 1.e-6;
    double phiMB    = ( -CONST_PI/(x+eps) + 2*(1 + 1/(x+eps))*atan(1/(x+eps)) + 2*(1 + 1/(x+eps))*log(1 + (x+eps)) - (1 - 1/(x+eps))*log(1 + pow(x, 2.)) + eps );
    //double phiMB  = ( -CONST_PI/x + 2*(1 + 1/x)*atan(1/x) + 2*(1 + 1/x)*log(1 + x) - (1 - 1/x)*log(1 + pow(x, 2.)) );
    //    // double phiMB    = (1/x) * (-CONST_PI + 2.*(1+x)*atan(1/x) + 2*(1+x)*log(1+x) + (1-x)*(log(1+ pow(x,2.))));
    phiMB          *= factor;  // converted to code units
    return phiMB;
}


double phi_modNFW(double R_val, double Z_val){
    /* Returns modified NFW potential (code units)
     * \params:  R_val - galactocentric radius (code units) 
     *           Z_val - height above the disc (code units) 
     */
    double MVir     = g_inputParam[M_VIR];   // solar units
    double H        = g_inputParam[HUNIV] * (1.e5 / (1.e6 * CONST_pc));        // km/s/Mpc --> cgs
    double rhoCrit  = 3 * H * H / (8. * CONST_PI * CONST_G);  // critical density of universe in cgs
    double r200     = pow((MVir * CONST_Msun / (200 * rhoCrit * 4 * CONST_PI / 3.)), 1/3.) / UNIT_LENGTH; // (code units)
    double c200     = g_inputParam[C200];  
    double dHalo    = g_inputParam[D_HALO] * 1.e3 * CONST_pc / UNIT_LENGTH; // (code units)
    double fc       = log(1+c200) - c200/(1+c200);
    double rs       = r200 / c200;
    
    double factor   = -CONST_G * MVir * CONST_Msun / (fc * rs * UNIT_LENGTH); // cgs
    double radDist2 = pow(R_val, 2.) + pow(Z_val, 2.) + pow(dHalo, 2.);
    double phimNFW  = (log(1 + sqrt(radDist2)/rs) / (sqrt(radDist2)/rs));
    phimNFW        *= (factor / pow(UNIT_VELOCITY, 2.));  // converted to code units
    return phimNFW;
}


double phi_Bl(double R_val, double Z_val){
    /* Returns Hernquist (1993) potential for Bulge (code units)
 *      * \params:  R_val - galactocentric radius (code units) 
 *           *           Z_val - height above the disc (code units) 
 *                */
    double UNIT_MASS  = UNIT_DENSITY * pow(UNIT_LENGTH,3);
    double rad    = sqrt(pow(R_val, 2.) + pow(Z_val, 2.));
    double MBulge = g_inputParam[M_BULGE];
    double radBl  = g_inputParam[R_BULGE] * 1.e3 * CONST_pc / UNIT_LENGTH;
    double factor = (CONST_G*pow(UNIT_LENGTH,2)*UNIT_DENSITY/pow(UNIT_VELOCITY, 2.));
    double phiBl  = factor * (MBulge*CONST_Msun/UNIT_MASS) / (rad + radBl);
    return phiBl;
}


double delPhiDelZ(double R_val, double Z_val){
    /* Calculates gradient of potential along z (code units)
     *
     * \params:  R_val - galactocentric radius (code units) 
     *           Z_val - height above the disc (code units)
     */ 
    double MDisk    = g_inputParam[M_STAR];  // solar units
    double aDisk    = g_inputParam[A_STAR] * 1.e3 * CONST_pc / UNIT_LENGTH;
    double bDisk    = g_inputParam[B_STAR] * 1.e3 * CONST_pc / UNIT_LENGTH;

    double MVir     = g_inputParam[M_VIR];   // solar units
    double H        = g_inputParam[HUNIV] * (1.e5 / (1.e6 * CONST_pc));        // km/s/Mpc --> cgs
    double rhoCrit  = 3 * H * H / (8. * CONST_PI * CONST_G);  // critical density of universe in cgs
    double r200     = pow((MVir * CONST_Msun / (200 * rhoCrit * 4 * CONST_PI / 3.)), 1/3.) / UNIT_LENGTH; // (code units)
    double c200     = g_inputParam[C200];   // (code units)
    double dHalo    = g_inputParam[D_HALO] * 1.e3 * CONST_pc / UNIT_LENGTH; // (code units)
    double fc       = log(1+c200) - c200/(1+c200);
    double rs       = r200 / c200;

    double factorMN = CONST_G * MDisk * CONST_Msun / pow(UNIT_LENGTH, 2.);  // cgs
    double dist2    = pow(Z_val, 2.) + pow(bDisk, 2.);
    double dist1    = pow(R_val, 2.) + pow(aDisk + sqrt(dist2), 2.);
    double dPhiDzMN = Z_val * (aDisk + sqrt(dist2));
    dPhiDzMN       /= (pow(dist1, 3/2.) * pow(dist2, 1/2.));
    dPhiDzMN       *= (factorMN * UNIT_LENGTH / pow(UNIT_VELOCITY, 2.)); // converted to code units
   
    double factorDM = CONST_G * MVir * CONST_Msun / (fc * rs * pow(UNIT_LENGTH, 2.)); // cgs
    double radDist2 = pow(R_val, 2.) + pow(Z_val, 2.) + pow(dHalo, 2.);
    double term1    = rs * log(1. + sqrt(radDist2)/rs) / pow(radDist2, 3/2.);
    double term2    = 1. / (radDist2 * (1. + sqrt(radDist2)/rs));
    double dPhiDzDM = Z_val * (term1 - term2);
    dPhiDzDM       *= (factorDM * UNIT_LENGTH / pow(UNIT_VELOCITY, 2.)); // converted to code units

    return (dPhiDzMN + dPhiDzDM);
}


double delPhiDelR(double R_val, double Z_val){
    /* Calculates gradient of potential along R (code units)
     *
     * \params:  R_val - galactocentric radius (code units) 
     *           Z_val - height above the disc (code units)
     */
    double MDisk    = g_inputParam[M_STAR];  // solar units
    double aDisk    = g_inputParam[A_STAR] * 1.e3 * CONST_pc / UNIT_LENGTH;
    double bDisk    = g_inputParam[B_STAR] * 1.e3 * CONST_pc / UNIT_LENGTH;

    double MVir     = g_inputParam[M_VIR];   // solar units
    double H        = 70 * (1.e5 / (1.e6 * CONST_pc));        // km/s/Mpc --> cgs
    double rhoCrit  = 3 * H * H / (8. * CONST_PI * CONST_G);  // critical density of universe in cgs
    double r200     = pow((MVir * CONST_Msun / (200 * rhoCrit * 4 * CONST_PI / 3.)), 1/3.) / UNIT_LENGTH; // (code units)
    double c200     = g_inputParam[C200];   // (code units)
    double dHalo    = g_inputParam[D_HALO] * 1.e3 * CONST_pc / UNIT_LENGTH; // (code units)
    double fc       = log(1+c200) - c200/(1+c200);
    double rs       = r200 / c200;

    double factorMN = CONST_G * MDisk * CONST_Msun / pow(UNIT_LENGTH, 2.);  // cgs
    double dist2    = pow(Z_val, 2.) + pow(bDisk, 2.);
    double dist1    = pow(R_val, 2.) + pow(aDisk + sqrt(dist2), 2.);
    double dPhiDRMN = R_val / (pow(dist1, 3/2.));
    dPhiDRMN       *= (factorMN * UNIT_LENGTH / pow(UNIT_VELOCITY, 2.)); // converted to code units

    double factorDM = CONST_G * MVir * CONST_Msun / (fc * rs * pow(UNIT_LENGTH, 2.)); // cgs
    double radDist2 = pow(R_val, 2.) + pow(Z_val, 2.) + pow(dHalo, 2.);
    double term1    = rs * log(1 + sqrt(radDist2)/rs) / (pow(radDist2, 3/2.));
    double term2    = 1./ (radDist2 * (1. + sqrt(radDist2)/rs)); 
    double dPhiDRDM = R_val * (term1 - term2);
    dPhiDRDM       *= (factorDM * UNIT_LENGTH / pow(UNIT_VELOCITY, 2.)); // converted to code units
    
    return (dPhiDRMN + dPhiDRDM);
}
