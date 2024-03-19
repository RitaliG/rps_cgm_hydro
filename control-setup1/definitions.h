#define  PHYSICS                        HD
#define  DIMENSIONS                     3
#define  GEOMETRY                       POLAR
#define  BODY_FORCE                     NO
#define  COOLING                        NO
#define  RECONSTRUCTION                 LINEAR
#define  TIME_STEPPING                  RK2
#define  NTRACER                        3
#define  PARTICLES                      NO
#define  USER_DEF_PARAMETERS            24

/* -- physics dependent declarations -- */

#define  DUST_FLUID                     NO
#define  EOS                            IDEAL
#define  ENTROPY_SWITCH                 NO
#define  THERMAL_CONDUCTION             NO
#define  VISCOSITY                      NO
#define  ROTATING_FRAME                 NO
#define  VERBOSE                        NO
#define  ID_NZ_MAX                      2
#define  WIND                           YES
#define  WIND_DOMAIN                    YES

/* -- user-defined parameters (labels) -- */

#define  M_VIR                          0
#define  HUNIV                          1
#define  C200                           2
#define  D_HALO                         3
#define  M_RDM                          4
#define  R_DM                           5
#define  M_STAR                         6
#define  A_STAR                         7
#define  B_STAR                         8
#define  M_BULGE                        9
#define  R_BULGE                        10
#define  ZMET                           11
#define  OFFSET                         12
#define  N_ICM                          13
#define  T_ICM                          14
#define  V_WIND                         15
#define  INCL                           16
#define  T_FLOOR                        17
#define  T_CGM                          18
#define  T_ISM                          19
#define  R_CGM                          20
#define  SIGMA_R_CGM                    21
#define  R_DISK_CYL                     22
#define  Z_DISK_CYL                     23

/* [Beg] user-defined constants (do not change this line) */

#define  UNIT_LENGTH                    (1.e3*CONST_pc)
#define  UNIT_DENSITY                   (1.e-4*CONST_mp)
#define  UNIT_VELOCITY                  1.e7

/* [End] user-defined constants (do not change this line) */
