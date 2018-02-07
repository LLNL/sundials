/*---------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *---------------------------------------------------------------
 * LLNS/SMU Copyright Start
 * Copyright (c) 2015, Southern Methodist University and 
 * Lawrence Livermore National Security
 *
 * This work was performed under the auspices of the U.S. Department 
 * of Energy by Southern Methodist University and Lawrence Livermore 
 * National Laboratory under Contract DE-AC52-07NA27344.
 * Produced at Southern Methodist University and the Lawrence 
 * Livermore National Laboratory.
 *
 * All rights reserved.
 * For details, see the LICENSE file.
 * LLNS/SMU Copyright End
 *---------------------------------------------------------------
 * This is the implementation file for Butcher table structure
 * (and built-in tables) for the ARKODE solver.
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>

#include "arkode_impl.h"
#include <sundials/sundials_math.h>

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define RSYM ".32Lg"
#else
#define RSYM ".16g"
#endif


/*---------------------------------------------------------------
  Utility routine to allocate a table structure
  ---------------------------------------------------------------*/
ARKodeButcherTable AllocButcherTable(int stages, booleantype embedded) {

  int i;
  ARKodeButcherTable B;

  /* Check for legal 'stages' value */
  if (stages < 1)  return(NULL);

  /* Allocate Butcher table structure */
  B = (ARKodeButcherTable) malloc(sizeof(struct ARKodeButcherTableMem));

  /* set stages into B structure */
  B->stages = stages;
  
  /* Allocate fields within Butcher table structure */
  B->A = (realtype **) calloc( stages, sizeof(realtype*) );
  for (i=0; i<stages; i++)
    B->A[i] = (realtype *) calloc( stages, sizeof(realtype) );
  B->b = (realtype *) calloc( stages, sizeof(realtype) );
  B->c = (realtype *) calloc( stages, sizeof(realtype) );
  if (embedded)
    B->d = (realtype *) calloc( stages, sizeof(realtype) );

  /* initialize order parameters */
  B->q = 0;
  B->p = 0;

  return(B);
}


/*---------------------------------------------------------------
  Utility routine to query the size of a Butcher table structure
  ---------------------------------------------------------------*/
void ButcherTableSpace(ARKodeButcherTable B, sunindextype *liw,
                       sunindextype *lrw)
{
  /* initialize outputs and return if B is not allocated */
  *liw = 0;  *lrw = 0;
  if (B == NULL)  return;

  /* fill outputs based on B */
  *liw = 3;
  if (B->d != NULL) {
    *lrw = B->stages * (B->stages + 3);
  } else {
    *lrw = B->stages * (B->stages + 2);
  }
}


/*---------------------------------------------------------------
  Utility routine to free a Butcher table structure
  ---------------------------------------------------------------*/
void FreeButcherTable(ARKodeButcherTable B) {

  int i;
  
  /* Free each field within Butcher table structure, and then 
     free structure itself */
  if (B != NULL) {
    if (B->d != NULL)  free(B->d);
    if (B->c != NULL)  free(B->c);
    if (B->b != NULL)  free(B->b);
    if (B->A != NULL) {
      for (i=0; i<B->stages; i++)
        if (B->A[i] != NULL)  free(B->A[i]);
      free(B->A);
    }

    free(B);
  }
}


/*---------------------------------------------------------------
  Utility routine to print a Butcher table structure
  ---------------------------------------------------------------*/
void WriteButcherTable(ARKodeButcherTable B, FILE *outfile) {
  int i, j;
  if (B == NULL) return;
  fprintf(outfile, "  A = \n");
  for (i=0; i<B->stages; i++) {
    fprintf(outfile, "      ");
    for (j=0; j<B->stages; j++)
      fprintf(outfile, "%"RSYM"  ", B->A[i][j]);
    fprintf(outfile, "\n");
  }
  fprintf(outfile, "  c = ");
  for (i=0; i<B->stages; i++) 
    fprintf(outfile, "%"RSYM"  ", B->c[i]);
  fprintf(outfile, "\n");
  fprintf(outfile, "  b = ");
  for (i=0; i<B->stages; i++) 
    fprintf(outfile, "%"RSYM"  ", B->b[i]);
  fprintf(outfile, "\n");
  fprintf(outfile, "  d = ");
  for (i=0; i<B->stages; i++) 
    fprintf(outfile, "%"RSYM"  ", B->d[i]);
  fprintf(outfile, "\n");
}


/*---------------------------------------------------------------
 Returns Butcher table structure for pre-set Runge Kutta methods.  

 Input:  imeth -- integer key for the desired method (see below)

 Allowed 'method' names and properties (those in an ARK pair are marked 
 with a *).  All method names are of the form <name>_s_p_q.  The method 
 'type' is one of 
    ERK -- explicit Runge Kutta
    SDIRK -- singly-diagonally implicit Runge Kutta
    SDIRK -- explicit [1st stage] singly-diagonally implicit Runge Kutta
 The 'A-stable' and 'L-stable' columns are based on numerical estimates 
 of each property.  The 'QP' column denotes whether the coefficients 
 of the method are known precisely enough for use in 'long double' 
 (128-bit) calculations.

   imeth                       type  A-stable  L-stable  QP 
  ----------------------------------------------------------
   HEUN_EULER_2_1_2             ERK     N         N       Y
   BOGACKI_SHAMPINE_4_2_3       ERK     N         N       Y
   ARK324L2SA_ERK_4_2_3*        ERK     N         N       N
   ZONNEVELD_5_3_4              ERK     N         N       Y
   ARK436L2SA_ERK_6_3_4*        ERK     N         N       N
   SAYFY_ABURUB_6_3_4           ERK     N         N       N
   CASH_KARP_6_4_5              ERK     N         N       Y
   FEHLBERG_6_4_5               ERK     N         N       Y
   DORMAND_PRINCE_7_4_5         ERK     N         N       Y
   ARK548L2SA_ERK_8_4_5*        ERK     N         N       N
   VERNER_8_5_6                 ERK     N         N       Y
   FEHLBERG_13_7_8              ERK     N         N       Y
   SDIRK_2_1_2                SDIRK     Y         N       Y
   BILLINGTON_3_3_2           SDIRK     N         N       N
   TRBDF2_3_3_2              ESDIRK     N         N       Y
   KVAERNO_4_2_3             ESDIRK     Y         Y       N
   ARK324L2SA_DIRK_4_2_3*    ESDIRK     Y         Y       N
   CASH_5_2_4                 SDIRK     Y         Y       N
   CASH_5_3_4                 SDIRK     Y         Y       N
   SDIRK_5_3_4                SDIRK     Y         Y       Y
   KVAERNO_5_3_4             ESDIRK     Y         N       N
   ARK436L2SA_DIRK_6_3_4*    ESDIRK     Y         Y       N
   KVAERNO_7_4_5             ESDIRK     Y         Y       N
   ARK548L2SA_DIRK_8_4_5*    ESDIRK     Y         Y       N
  ----------------------------------------------------------

---------------------------------------------------------------*/
ARKodeButcherTable ARKodeLoadButcherTable(int imethod) 
{

  ARKodeButcherTable B;
  B = NULL;

  /* fill in coefficients based on method name */
  switch(imethod) {

  case(HEUN_EULER_2_1_2):    /* Heun-Euler-ERK */
    B = AllocButcherTable(2, SUNTRUE);
    B->q = 2;
    B->p = 1;
      
    B->A[1][0] = RCONST(1.0);

    B->b[0] = RCONST(1.0)/RCONST(2.0);
    B->b[1] = RCONST(1.0)/RCONST(2.0);

    B->d[0] = RCONST(1.0);

    B->c[1] = RCONST(1.0);
    break;

  case(BOGACKI_SHAMPINE_4_2_3):    /* Bogacki-Shampine-ERK */
    B = AllocButcherTable(4, SUNTRUE);
    B->q = 3;
    B->p = 2;
    B->A[1][0] = RCONST(1.0)/RCONST(2.0);
    B->A[2][1] = RCONST(3.0)/RCONST(4.0);
    B->A[3][0] = RCONST(2.0)/RCONST(9.0);
    B->A[3][1] = RCONST(1.0)/RCONST(3.0);
    B->A[3][2] = RCONST(4.0)/RCONST(9.0);

    B->b[0] = RCONST(2.0)/RCONST(9.0);
    B->b[1] = RCONST(1.0)/RCONST(3.0);
    B->b[2] = RCONST(4.0)/RCONST(9.0);

    B->d[0] = RCONST(7.0)/RCONST(24.0);
    B->d[1] = RCONST(1.0)/RCONST(4.0);
    B->d[2] = RCONST(1.0)/RCONST(3.0);
    B->d[3] = RCONST(1.0)/RCONST(8.0);

    B->c[1] = RCONST(1.0)/RCONST(2.0);
    B->c[2] = RCONST(3.0)/RCONST(4.0);
    B->c[3] = RCONST(1.0);
    break;

  case(ARK324L2SA_ERK_4_2_3):    /* ARK3(2)4L[2]SA-ERK */
    B = AllocButcherTable(4, SUNTRUE);
    B->q = 3;
    B->p = 2;
    B->A[1][0] = RCONST(1767732205903.0)/RCONST(2027836641118.0);
    B->A[2][0] = RCONST(5535828885825.0)/RCONST(10492691773637.0);
    B->A[2][1] = RCONST(788022342437.0)/RCONST(10882634858940.0);
    B->A[3][0] = RCONST(6485989280629.0)/RCONST(16251701735622.0);
    B->A[3][1] = RCONST(-4246266847089.0)/RCONST(9704473918619.0);
    B->A[3][2] = RCONST(10755448449292.0)/RCONST(10357097424841.0);

    B->b[0] = RCONST(1471266399579.0)/RCONST(7840856788654.0);
    B->b[1] = RCONST(-4482444167858.0)/RCONST(7529755066697.0);
    B->b[2] = RCONST(11266239266428.0)/RCONST(11593286722821.0);
    B->b[3] = RCONST(1767732205903.0)/RCONST(4055673282236.0);

    B->d[0] = RCONST(2756255671327.0)/RCONST(12835298489170.0);
    B->d[1] = RCONST(-10771552573575.0)/RCONST(22201958757719.0);
    B->d[2] = RCONST(9247589265047.0)/RCONST(10645013368117.0);
    B->d[3] = RCONST(2193209047091.0)/RCONST(5459859503100.0);

    B->c[1] = RCONST(1767732205903.0)/RCONST(2027836641118.0);
    B->c[2] = RCONST(3.0)/RCONST(5.0);
    B->c[3] = RCONST(1.0);
    break;

  case(ZONNEVELD_5_3_4):    /* Zonneveld */
    B = AllocButcherTable(5, SUNTRUE);
    B->q = 4;
    B->p = 3;
    B->A[1][0] = RCONST(0.5);
    B->A[2][1] = RCONST(0.5);
    B->A[3][2] = RCONST(1.0);
    B->A[4][0] = RCONST(5.0)/RCONST(32.0);
    B->A[4][1] = RCONST(7.0)/RCONST(32.0);
    B->A[4][2] = RCONST(13.0)/RCONST(32.0);
    B->A[4][3] = RCONST(-1.0)/RCONST(32.0);

    B->b[0] = RCONST(1.0)/RCONST(6.0);
    B->b[1] = RCONST(1.0)/RCONST(3.0);
    B->b[2] = RCONST(1.0)/RCONST(3.0);
    B->b[3] = RCONST(1.0)/RCONST(6.0);

    B->d[0] = RCONST(-1.0)/RCONST(2.0);
    B->d[1] = RCONST(7.0)/RCONST(3.0);
    B->d[2] = RCONST(7.0)/RCONST(3.0);
    B->d[3] = RCONST(13.0)/RCONST(6.0);
    B->d[4] = RCONST(-16.0)/RCONST(3.0);

    B->c[1] = RCONST(0.5);
    B->c[2] = RCONST(0.5);
    B->c[3] = RCONST(1.0);
    B->c[4] = RCONST(0.75);
    break;

  case(ARK436L2SA_ERK_6_3_4):    /* ARK4(3)6L[2]SA-ERK */
    B = AllocButcherTable(6, SUNTRUE);
    B->q = 4;
    B->p = 3;
    B->A[1][0] = RCONST(0.5);
    B->A[2][0] = RCONST(13861.0)/RCONST(62500.0);
    B->A[2][1] = RCONST(6889.0)/RCONST(62500.0);
    B->A[3][0] = RCONST(-116923316275.0)/RCONST(2393684061468.0);
    B->A[3][1] = RCONST(-2731218467317.0)/RCONST(15368042101831.0);
    B->A[3][2] = RCONST(9408046702089.0)/RCONST(11113171139209.0);
    B->A[4][0] = RCONST(-451086348788.0)/RCONST(2902428689909.0);
    B->A[4][1] = RCONST(-2682348792572.0)/RCONST(7519795681897.0);
    B->A[4][2] = RCONST(12662868775082.0)/RCONST(11960479115383.0);
    B->A[4][3] = RCONST(3355817975965.0)/RCONST(11060851509271.0);
    B->A[5][0] = RCONST(647845179188.0)/RCONST(3216320057751.0);
    B->A[5][1] = RCONST(73281519250.0)/RCONST(8382639484533.0);
    B->A[5][2] = RCONST(552539513391.0)/RCONST(3454668386233.0);
    B->A[5][3] = RCONST(3354512671639.0)/RCONST(8306763924573.0);
    B->A[5][4] = RCONST(4040.0)/RCONST(17871.0);

    B->b[0] = RCONST(82889.0)/RCONST(524892.0);
    B->b[2] = RCONST(15625.0)/RCONST(83664.0);
    B->b[3] = RCONST(69875.0)/RCONST(102672.0);
    B->b[4] = RCONST(-2260.0)/RCONST(8211.0);
    B->b[5] = RCONST(1.0)/RCONST(4.0);

    B->d[0] = RCONST(4586570599.0)/RCONST(29645900160.0);
    B->d[2] = RCONST(178811875.0)/RCONST(945068544.0);
    B->d[3] = RCONST(814220225.0)/RCONST(1159782912.0);
    B->d[4] = RCONST(-3700637.0)/RCONST(11593932.0);
    B->d[5] = RCONST(61727.0)/RCONST(225920.0);

    B->c[1] = RCONST(1.0)/RCONST(2.0);
    B->c[2] = RCONST(83.0)/RCONST(250.0);
    B->c[3] = RCONST(31.0)/RCONST(50.0);
    B->c[4] = RCONST(17.0)/RCONST(20.0);
    B->c[5] = RCONST(1.0);
    break;

  case(SAYFY_ABURUB_6_3_4):    /* Sayfy-Aburub-4-3-ERK */
    B = AllocButcherTable(6, SUNTRUE);
    B->q = 4;
    B->p = 3;
    B->A[1][0] = RCONST(1.0)/RCONST(2.0); 
    B->A[2][0] = RCONST(-1.0); 
    B->A[2][1] = RCONST(2.0); 
    B->A[3][0] = RCONST(1.0)/RCONST(6.0);
    B->A[3][1] = RCONST(2.0)/RCONST(3.0);
    B->A[3][2] = RCONST(1.0)/RCONST(6.0);
    B->A[4][0] = RCONST(0.137);
    B->A[4][1] = RCONST(0.226);
    B->A[4][2] = RCONST(0.137);
    B->A[5][0] = RCONST(0.452);
    B->A[5][1] = RCONST(-0.904);
    B->A[5][2] = RCONST(-0.548);
    B->A[5][4] = RCONST(2.0);

    B->b[0] = RCONST(1.0)/RCONST(6.0);
    B->b[1] = RCONST(1.0)/RCONST(3.0);
    B->b[2] = RCONST(1.0)/RCONST(12.0);
    B->b[3] = RCONST(0.0);
    B->b[4] = RCONST(1.0)/RCONST(3.0);
    B->b[5] = RCONST(1.0)/RCONST(12.0);

    B->d[0] = RCONST(1.0)/RCONST(6.0);
    B->d[1] = RCONST(2.0)/RCONST(3.0);
    B->d[2] = RCONST(1.0)/RCONST(6.0);

    B->c[1] = RCONST(1.0)/RCONST(2.0);
    B->c[2] = RCONST(1.0);
    B->c[3] = RCONST(1.0);
    B->c[4] = RCONST(1.0)/RCONST(2.0);
    B->c[5] = RCONST(1.0);
    break;

  case(CASH_KARP_6_4_5):    /* Cash-Karp-ERK */
    B = AllocButcherTable(6, SUNTRUE);
    B->q = 5;
    B->p = 4;
    B->A[1][0] = RCONST(1.0)/RCONST(5.0);
    B->A[2][0] = RCONST(3.0)/RCONST(40.0);
    B->A[2][1] = RCONST(9.0)/RCONST(40.0);
    B->A[3][0] = RCONST(3.0)/RCONST(10.0);
    B->A[3][1] = RCONST(-9.0)/RCONST(10.0);
    B->A[3][2] = RCONST(6.0)/RCONST(5.0);
    B->A[4][0] = RCONST(-11.0)/RCONST(54.0);
    B->A[4][1] = RCONST(5.0)/RCONST(2.0);
    B->A[4][2] = RCONST(-70.0)/RCONST(27.0);
    B->A[4][3] = RCONST(35.0)/RCONST(27.0);
    B->A[5][0] = RCONST(1631.0)/RCONST(55296.0);
    B->A[5][1] = RCONST(175.0)/RCONST(512.0);
    B->A[5][2] = RCONST(575.0)/RCONST(13824.0);
    B->A[5][3] = RCONST(44275.0)/RCONST(110592.0);
    B->A[5][4] = RCONST(253.0)/RCONST(4096.0);

    B->b[0] = RCONST(37.0)/RCONST(378.0);
    B->b[2] = RCONST(250.0)/RCONST(621.0);
    B->b[3] = RCONST(125.0)/RCONST(594.0);
    B->b[5] = RCONST(512.0)/RCONST(1771.0);

    B->d[0] = RCONST(2825.0)/RCONST(27648.0);
    B->d[2] = RCONST(18575.0)/RCONST(48384.0);
    B->d[3] = RCONST(13525.0)/RCONST(55296.0);
    B->d[4] = RCONST(277.0)/RCONST(14336.0);
    B->d[5] = RCONST(1.0)/RCONST(4.0);

    B->c[1] = RCONST(1.0)/RCONST(5.0);
    B->c[2] = RCONST(3.0)/RCONST(10.0);
    B->c[3] = RCONST(3.0)/RCONST(5.0);
    B->c[4] = RCONST(1.0);
    B->c[5] = RCONST(7.0)/RCONST(8.0);
    break;

  case(FEHLBERG_6_4_5):    /* Fehlberg-ERK */
    B = AllocButcherTable(6, SUNTRUE);
    B->q = 5;
    B->p = 4;
    B->A[1][0] = RCONST(1.0)/RCONST(4.0);
    B->A[2][0] = RCONST(3.0)/RCONST(32.0);
    B->A[2][1] = RCONST(9.0)/RCONST(32.0);
    B->A[3][0] = RCONST(1932.0)/RCONST(2197.0);
    B->A[3][1] = RCONST(-7200.0)/RCONST(2197.0);
    B->A[3][2] = RCONST(7296.0)/RCONST(2197.0);
    B->A[4][0] = RCONST(439.0)/RCONST(216.0);
    B->A[4][1] = RCONST(-8.0);
    B->A[4][2] = RCONST(3680.0)/RCONST(513.0);
    B->A[4][3] = RCONST(-845.0)/RCONST(4104.0);
    B->A[5][0] = RCONST(-8.0)/RCONST(27.0);
    B->A[5][1] = RCONST(2.0);
    B->A[5][2] = RCONST(-3544.0)/RCONST(2565.0);
    B->A[5][3] = RCONST(1859.0)/RCONST(4104.0);
    B->A[5][4] = RCONST(-11.0)/RCONST(40.0);

    B->b[0] = RCONST(16.0)/RCONST(135.0);
    B->b[2] = RCONST(6656.0)/RCONST(12825.0);
    B->b[3] = RCONST(28561.0)/RCONST(56430.0);
    B->b[4] = RCONST(-9.0)/RCONST(50.0);
    B->b[5] = RCONST(2.0)/RCONST(55.0);

    B->d[0] = RCONST(25.0)/RCONST(216.0);
    B->d[2] = RCONST(1408.0)/RCONST(2565.0);
    B->d[3] = RCONST(2197.0)/RCONST(4104.0);
    B->d[4] = RCONST(-1.0)/RCONST(5.0);
      
    B->c[1] = RCONST(1.0)/RCONST(4.0);
    B->c[2] = RCONST(3.0)/RCONST(8.0);
    B->c[3] = RCONST(12.0)/RCONST(13.0);
    B->c[4] = RCONST(1.0);
    B->c[5] = RCONST(1.0)/RCONST(2.0);
    break;

  case(DORMAND_PRINCE_7_4_5):    /* Dormand-Prince-ERK */
    B = AllocButcherTable(7, SUNTRUE);
    B->q = 5;
    B->p = 4;
    B->A[1][0] = RCONST(1.0)/RCONST(5.0);
    B->A[2][0] = RCONST(3.0)/RCONST(40.0);
    B->A[2][1] = RCONST(9.0)/RCONST(40.0);
    B->A[3][0] = RCONST(44.0)/RCONST(45.0);
    B->A[3][1] = RCONST(-56.0)/RCONST(15.0);
    B->A[3][2] = RCONST(32.0)/RCONST(9.0);
    B->A[4][0] = RCONST(19372.0)/RCONST(6561.0);
    B->A[4][1] = RCONST(-25360.0)/RCONST(2187.0);
    B->A[4][2] = RCONST(64448.0)/RCONST(6561.0);
    B->A[4][3] = RCONST(-212.0)/RCONST(729.0);
    B->A[5][0] = RCONST(9017.0)/RCONST(3168.0);
    B->A[5][1] = RCONST(-355.0)/RCONST(33.0);
    B->A[5][2] = RCONST(46732.0)/RCONST(5247.0);
    B->A[5][3] = RCONST(49.0)/RCONST(176.0);
    B->A[5][4] = RCONST(-5103.0)/RCONST(18656.0);
    B->A[6][0] = RCONST(35.0)/RCONST(384.0);
    B->A[6][2] = RCONST(500.0)/RCONST(1113.0);
    B->A[6][3] = RCONST(125.0)/RCONST(192.0);
    B->A[6][4] = RCONST(-2187.0)/RCONST(6784.0);
    B->A[6][5] = RCONST(11.0)/RCONST(84.0);

    B->b[0] = RCONST(35.0)/RCONST(384.0);
    B->b[2] = RCONST(500.0)/RCONST(1113.0);
    B->b[3] = RCONST(125.0)/RCONST(192.0);
    B->b[4] = RCONST(-2187.0)/RCONST(6784.0);
    B->b[5] = RCONST(11.0)/RCONST(84.0);

    B->d[0] = RCONST(5179.0)/RCONST(57600.0);
    B->d[2] = RCONST(7571.0)/RCONST(16695.0);
    B->d[3] = RCONST(393.0)/RCONST(640.0);
    B->d[4] = RCONST(-92097.0)/RCONST(339200.0);
    B->d[5] = RCONST(187.0)/RCONST(2100.0);
    B->d[6] = RCONST(1.0)/RCONST(40.0);

    B->c[1] = RCONST(1.0)/RCONST(5.0);
    B->c[2] = RCONST(3.0)/RCONST(10.0);
    B->c[3] = RCONST(4.0)/RCONST(5.0);
    B->c[4] = RCONST(8.0)/RCONST(9.0);
    B->c[5] = RCONST(1.0);
    B->c[6] = RCONST(1.0);
    break;

  case(ARK548L2SA_ERK_8_4_5):    /* ARK5(4)8L[2]SA-ERK */
    B = AllocButcherTable(8, SUNTRUE);
    B->q = 5;
    B->p = 4;
    B->A[1][0] = RCONST(41.0)/RCONST(100.0);
    B->A[2][0] = RCONST(367902744464.0)/RCONST(2072280473677.0);
    B->A[2][1] = RCONST(677623207551.0)/RCONST(8224143866563.0);
    B->A[3][0] = RCONST(1268023523408.0)/RCONST(10340822734521.0);
    B->A[3][2] = RCONST(1029933939417.0)/RCONST(13636558850479.0);
    B->A[4][0] = RCONST(14463281900351.0)/RCONST(6315353703477.0);
    B->A[4][2] = RCONST(66114435211212.0)/RCONST(5879490589093.0);
    B->A[4][3] = RCONST(-54053170152839.0)/RCONST(4284798021562.0);
    B->A[5][0] = RCONST(14090043504691.0)/RCONST(34967701212078.0);
    B->A[5][2] = RCONST(15191511035443.0)/RCONST(11219624916014.0);
    B->A[5][3] = RCONST(-18461159152457.0)/RCONST(12425892160975.0);
    B->A[5][4] = RCONST(-281667163811.0)/RCONST(9011619295870.0);
    B->A[6][0] = RCONST(19230459214898.0)/RCONST(13134317526959.0);
    B->A[6][2] = RCONST(21275331358303.0)/RCONST(2942455364971.0);
    B->A[6][3] = RCONST(-38145345988419.0)/RCONST(4862620318723.0);
    B->A[6][4] = RCONST(-1.0)/RCONST(8.0);
    B->A[6][5] = RCONST(-1.0)/RCONST(8.0);
    B->A[7][0] = RCONST(-19977161125411.0)/RCONST(11928030595625.0);
    B->A[7][2] = RCONST(-40795976796054.0)/RCONST(6384907823539.0);
    B->A[7][3] = RCONST(177454434618887.0)/RCONST(12078138498510.0);
    B->A[7][4] = RCONST(782672205425.0)/RCONST(8267701900261.0);
    B->A[7][5] = RCONST(-69563011059811.0)/RCONST(9646580694205.0);
    B->A[7][6] = RCONST(7356628210526.0)/RCONST(4942186776405.0);

    B->b[0] = RCONST(-872700587467.0)/RCONST(9133579230613.0);
    B->b[3] = RCONST(22348218063261.0)/RCONST(9555858737531.0);
    B->b[4] = RCONST(-1143369518992.0)/RCONST(8141816002931.0);
    B->b[5] = RCONST(-39379526789629.0)/RCONST(19018526304540.0);
    B->b[6] = RCONST(32727382324388.0)/RCONST(42900044865799.0);
    B->b[7] = RCONST(41.0)/RCONST(200.0);

    B->d[0] = RCONST(-975461918565.0)/RCONST(9796059967033.0);
    B->d[3] = RCONST(78070527104295.0)/RCONST(32432590147079.0);
    B->d[4] = RCONST(-548382580838.0)/RCONST(3424219808633.0);
    B->d[5] = RCONST(-33438840321285.0)/RCONST(15594753105479.0);
    B->d[6] = RCONST(3629800801594.0)/RCONST(4656183773603.0);
    B->d[7] = RCONST(4035322873751.0)/RCONST(18575991585200.0);

    B->c[1] = RCONST(41.0)/RCONST(100.0);
    B->c[2] = RCONST(2935347310677.0)/RCONST(11292855782101.0);
    B->c[3] = RCONST(1426016391358.0)/RCONST(7196633302097.0);
    B->c[4] = RCONST(92.0)/RCONST(100.0);
    B->c[5] = RCONST(24.0)/RCONST(100.0);
    B->c[6] = RCONST(3.0)/RCONST(5.0);
    B->c[7] = RCONST(1.0);
    break;

  case(VERNER_8_5_6):    /* Verner-6-5 */
    B = AllocButcherTable(8, SUNTRUE);
    B->q = 6;
    B->p = 5;
    B->A[1][0] = RCONST(1.0)/RCONST(6.0);
    B->A[2][0] = RCONST(4.0)/RCONST(75.0);
    B->A[2][1] = RCONST(16.0)/RCONST(75.0);
    B->A[3][0] = RCONST(5.0)/RCONST(6.0);
    B->A[3][1] = RCONST(-8.0)/RCONST(3.0);
    B->A[3][2] = RCONST(5.0)/RCONST(2.0);
    B->A[4][0] = RCONST(-165.0)/RCONST(64.0);
    B->A[4][1] = RCONST(55.0)/RCONST(6.0);
    B->A[4][2] = RCONST(-425.0)/RCONST(64.0);
    B->A[4][3] = RCONST(85.0)/RCONST(96.0);
    B->A[5][0] = RCONST(12.0)/RCONST(5.0);
    B->A[5][1] = RCONST(-8.0);
    B->A[5][2] = RCONST(4015.0)/RCONST(612.0);
    B->A[5][3] = RCONST(-11.0)/RCONST(36.0);
    B->A[5][4] = RCONST(88.0)/RCONST(255.0);
    B->A[6][0] = RCONST(-8263.0)/RCONST(15000.0);
    B->A[6][1] = RCONST(124.0)/RCONST(75.0);
    B->A[6][2] = RCONST(-643.0)/RCONST(680.0);
    B->A[6][3] = RCONST(-81.0)/RCONST(250.0);
    B->A[6][4] = RCONST(2484.0)/RCONST(10625.0);
    B->A[7][0] = RCONST(3501.0)/RCONST(1720.0);
    B->A[7][1] = RCONST(-300.0)/RCONST(43.0);
    B->A[7][2] = RCONST(297275.0)/RCONST(52632.0);
    B->A[7][3] = RCONST(-319.0)/RCONST(2322.0);
    B->A[7][4] = RCONST(24068.0)/RCONST(84065.0);
    B->A[7][6] = RCONST(3850.0)/RCONST(26703.0);

    B->b[0] = RCONST(3.0)/RCONST(40.0);
    B->b[2] = RCONST(875.0)/RCONST(2244.0);
    B->b[3] = RCONST(23.0)/RCONST(72.0);
    B->b[4] = RCONST(264.0)/RCONST(1955.0);
    B->b[6] = RCONST(125.0)/RCONST(11592.0);
    B->b[7] = RCONST(43.0)/RCONST(616.0);

    B->d[0] = RCONST(13.0)/RCONST(160.0);
    B->d[2] = RCONST(2375.0)/RCONST(5984.0);
    B->d[3] = RCONST(5.0)/RCONST(16.0);
    B->d[4] = RCONST(12.0)/RCONST(85.0);
    B->d[5] = RCONST(3.0)/RCONST(44.0);

    B->c[0] = RCONST(0.0);
    B->c[1] = RCONST(1.0)/RCONST(6.0);
    B->c[2] = RCONST(4.0)/RCONST(15.0);
    B->c[3] = RCONST(2.0)/RCONST(3.0);
    B->c[4] = RCONST(5.0)/RCONST(6.0);
    B->c[5] = RCONST(1.0);
    B->c[6] = RCONST(1.0)/RCONST(15.0);
    B->c[7] = RCONST(1.0);
    break;

  case(FEHLBERG_13_7_8):    /* Fehlberg-8-7 */
    B = AllocButcherTable(13, SUNTRUE);
    B->q = 8;
    B->p = 7;
    B->A[1][0] = RCONST(2.0)/RCONST(27.0);
    B->A[2][0] = RCONST(1.0)/RCONST(36.0);
    B->A[2][1] = RCONST(1.0)/RCONST(12.0);
    B->A[3][0] = RCONST(1.0)/RCONST(24.0);
    B->A[3][2] = RCONST(1.0)/RCONST(8.0);
    B->A[4][0] = RCONST(5.0)/RCONST(12.0);
    B->A[4][2] = RCONST(-25.0)/RCONST(16.0);
    B->A[4][3] = RCONST(25.0)/RCONST(16.0);
    B->A[5][0] = RCONST(1.0)/RCONST(20.0);
    B->A[5][3] = RCONST(1.0)/RCONST(4.0);
    B->A[5][4] = RCONST(1.0)/RCONST(5.0);
    B->A[6][0] = RCONST(-25.0)/RCONST(108.0);
    B->A[6][3] = RCONST(125.0)/RCONST(108.0);
    B->A[6][4] = RCONST(-65.0)/RCONST(27.0);
    B->A[6][5] = RCONST(125.0)/RCONST(54.0);
    B->A[7][0] = RCONST(31.0)/RCONST(300.0);
    B->A[7][4] = RCONST(61.0)/RCONST(225.0);
    B->A[7][5] = RCONST(-2.0)/RCONST(9.0);
    B->A[7][6] = RCONST(13.0)/RCONST(900.0);
    B->A[8][0] = RCONST(2.0);
    B->A[8][3] = RCONST(-53.0)/RCONST(6.0);
    B->A[8][4] = RCONST(704.0)/RCONST(45.0);
    B->A[8][5] = RCONST(-107.0)/RCONST(9.0);
    B->A[8][6] = RCONST(67.0)/RCONST(90.0);
    B->A[8][7] = RCONST(3.0);
    B->A[9][0] = RCONST(-91.0)/RCONST(108.0);
    B->A[9][3] = RCONST(23.0)/RCONST(108.0);
    B->A[9][4] = RCONST(-976.0)/RCONST(135.0);
    B->A[9][5] = RCONST(311.0)/RCONST(54.0);
    B->A[9][6] = RCONST(-19.0)/RCONST(60.0);
    B->A[9][7] = RCONST(17.0)/RCONST(6.0);
    B->A[9][8] = RCONST(-1.0)/RCONST(12.0);
    B->A[10][0] = RCONST(2383.0)/RCONST(4100.0);
    B->A[10][3] = RCONST(-341.0)/RCONST(164.0);
    B->A[10][4] = RCONST(4496.0)/RCONST(1025.0);
    B->A[10][5] = RCONST(-301.0)/RCONST(82.0);
    B->A[10][6] = RCONST(2133.0)/RCONST(4100.0);
    B->A[10][7] = RCONST(45.0)/RCONST(82.0);
    B->A[10][8] = RCONST(45.0)/RCONST(164.0);
    B->A[10][9] = RCONST(18.0)/RCONST(41.0);
    B->A[11][0] = RCONST(3.0)/RCONST(205.0);
    B->A[11][5] = RCONST(-6.0)/RCONST(41.0);
    B->A[11][6] = RCONST(-3.0)/RCONST(205.0);
    B->A[11][7] = RCONST(-3.0)/RCONST(41.0);
    B->A[11][8] = RCONST(3.0)/RCONST(41.0);
    B->A[11][9] = RCONST(6.0)/RCONST(41.0);
    B->A[12][0] = RCONST(-1777.0)/RCONST(4100.0);
    B->A[12][3] = RCONST(-341.0)/RCONST(164.0);
    B->A[12][4] = RCONST(4496.0)/RCONST(1025.0);
    B->A[12][5] = RCONST(-289.0)/RCONST(82.0);
    B->A[12][6] = RCONST(2193.0)/RCONST(4100.0);
    B->A[12][7] = RCONST(51.0)/RCONST(82.0);
    B->A[12][8] = RCONST(33.0)/RCONST(164.0);
    B->A[12][9] = RCONST(12.0)/RCONST(41.0);
    B->A[12][11] = RCONST(1.0);

    B->b[5]  = RCONST(34.0)/RCONST(105.0);
    B->b[6]  = RCONST(9.0)/RCONST(35.0);
    B->b[7]  = RCONST(9.0)/RCONST(35.0);
    B->b[8]  = RCONST(9.0)/RCONST(280.0);
    B->b[9]  = RCONST(9.0)/RCONST(280.0);
    B->b[11] = RCONST(41.0)/RCONST(840.0);
    B->b[12] = RCONST(41.0)/RCONST(840.0);

    B->d[0]  = RCONST(41.0)/RCONST(840.0);
    B->d[5]  = RCONST(34.0)/RCONST(105.0);
    B->d[6]  = RCONST(9.0)/RCONST(35.0);
    B->d[7]  = RCONST(9.0)/RCONST(35.0);
    B->d[8]  = RCONST(9.0)/RCONST(280.0);
    B->d[9]  = RCONST(9.0)/RCONST(280.0);
    B->d[10] = RCONST(41.0)/RCONST(840.0);

    B->c[1]  = RCONST(2.0)/RCONST(27.0);
    B->c[2]  = RCONST(1.0)/RCONST(9.0);
    B->c[3]  = RCONST(1.0)/RCONST(6.0);
    B->c[4]  = RCONST(5.0)/RCONST(12.0);
    B->c[5]  = RCONST(1.0)/RCONST(2.0);
    B->c[6]  = RCONST(5.0)/RCONST(6.0);
    B->c[7]  = RCONST(1.0)/RCONST(6.0);
    B->c[8]  = RCONST(2.0)/RCONST(3.0);
    B->c[9]  = RCONST(1.0)/RCONST(3.0);
    B->c[10] = RCONST(1.0);
    B->c[12] = RCONST(1.0);
    break;

  case(SDIRK_2_1_2):   /* SDIRK-2-1 (A,B stable) */
    B = AllocButcherTable(2, SUNTRUE);
    B->q = 2;
    B->p = 1;
      
    B->A[0][0] = RCONST(1.0);
    B->A[1][0] = RCONST(-1.0);
    B->A[1][1] = RCONST(1.0);

    B->b[0] = RCONST(0.5);
    B->b[1] = RCONST(0.5);

    B->d[0] = RCONST(1.0);

    B->c[0] = RCONST(1.0);
    B->c[1] = RCONST(0.0);
    break;

  case(BILLINGTON_3_3_2):    /* Billington-SDIRK */
    B = AllocButcherTable(3, SUNTRUE);
    B->q = 2;
    B->p = 3;

    B->A[0][0] = RCONST(0.292893218813);
    B->A[1][0] = RCONST(0.798989873223);
    B->A[1][1] = RCONST(0.292893218813);
    B->A[2][0] = RCONST(0.740789228841);
    B->A[2][1] = RCONST(0.259210771159);
    B->A[2][2] = RCONST(0.292893218813);

    B->d[0] = RCONST(0.691665115992);
    B->d[1] = RCONST(0.503597029883);
    B->d[2] = RCONST(-0.195262145876);

    B->b[0] = RCONST(0.740789228840);
    B->b[1] = RCONST(0.259210771159);

    B->c[0] = RCONST(0.292893218813);
    B->c[1] = RCONST(1.091883092037);
    B->c[2] = RCONST(1.292893218813);
    break;

  case(TRBDF2_3_3_2):    /* TRBDF2-ESDIRK */
    B = AllocButcherTable(3, SUNTRUE);
    B->q = 2;
    B->p = 3;

    B->A[1][0] = (RCONST(2.0)-SUNRsqrt(RCONST(2.0)))/RCONST(2.0);
    B->A[1][1] = (RCONST(2.0)-SUNRsqrt(RCONST(2.0)))/RCONST(2.0);
    B->A[2][0] = SUNRsqrt(RCONST(2.0))/RCONST(4.0);
    B->A[2][1] = SUNRsqrt(RCONST(2.0))/RCONST(4.0);
    B->A[2][2] = (RCONST(2.0)-SUNRsqrt(RCONST(2.0)))/RCONST(2.0);

    B->d[0] = (RCONST(1.0)-SUNRsqrt(RCONST(2.0))/RCONST(4.0))/RCONST(3.0);
    B->d[1] = (RCONST(3.0)*SUNRsqrt(RCONST(2.0))/RCONST(4.0)+RCONST(1.0))/RCONST(3.0);
    B->d[2] = (RCONST(2.0)-SUNRsqrt(RCONST(2.0)))/RCONST(6.0);

    B->b[0] = SUNRsqrt(RCONST(2.0))/RCONST(4.0);
    B->b[1] = SUNRsqrt(RCONST(2.0))/RCONST(4.0);
    B->b[2] = (RCONST(2.0)-SUNRsqrt(RCONST(2.0)))/RCONST(2.0);

    B->c[1] = RCONST(2.0)-SUNRsqrt(RCONST(2.0));
    B->c[2] = RCONST(1.0);
    break;

  case(KVAERNO_4_2_3):    /* Kvaerno(4,2,3)-ESDIRK */
    B = AllocButcherTable(4, SUNTRUE);
    B->q = 3;
    B->p = 2;
    B->A[1][0] = RCONST(0.4358665215);
    B->A[1][1] = RCONST(0.4358665215);
    B->A[2][0] = RCONST(0.490563388419108);
    B->A[2][1] = RCONST(0.073570090080892);
    B->A[2][2] = RCONST(0.4358665215);
    B->A[3][0] = RCONST(0.308809969973036);
    B->A[3][1] = RCONST(1.490563388254106);
    B->A[3][2] = RCONST(-1.235239879727145);
    B->A[3][3] = RCONST(0.4358665215);

    B->b[0] = RCONST(0.308809969973036);
    B->b[1] = RCONST(1.490563388254106);
    B->b[2] = RCONST(-1.235239879727145);
    B->b[3] = RCONST(0.4358665215);

    B->d[0] = RCONST(0.490563388419108);
    B->d[1] = RCONST(0.073570090080892);
    B->d[2] = RCONST(0.4358665215);

    B->c[1] = RCONST(0.871733043);
    B->c[2] = RCONST(1.0);
    B->c[3] = RCONST(1.0);
    break;

  case(ARK324L2SA_DIRK_4_2_3):    /* ARK3(2)4L[2]SA-ESDIRK */
    B = AllocButcherTable(4, SUNTRUE);
    B->q = 3;
    B->p = 2;
    B->A[1][0] = RCONST(1767732205903.0)/RCONST(4055673282236.0);
    B->A[1][1] = RCONST(1767732205903.0)/RCONST(4055673282236.0);
    B->A[2][0] = RCONST(2746238789719.0)/RCONST(10658868560708.0);
    B->A[2][1] = RCONST(-640167445237.0)/RCONST(6845629431997.0);
    B->A[2][2] = RCONST(1767732205903.0)/RCONST(4055673282236.0);
    B->A[3][0] = RCONST(1471266399579.0)/RCONST(7840856788654.0);
    B->A[3][1] = RCONST(-4482444167858.0)/RCONST(7529755066697.0);
    B->A[3][2] = RCONST(11266239266428.0)/RCONST(11593286722821.0);
    B->A[3][3] = RCONST(1767732205903.0)/RCONST(4055673282236.0);

    B->b[0] = RCONST(1471266399579.0)/RCONST(7840856788654.0);
    B->b[1] = RCONST(-4482444167858.0)/RCONST(7529755066697.0);
    B->b[2] = RCONST(11266239266428.0)/RCONST(11593286722821.0);
    B->b[3] = RCONST(1767732205903.0)/RCONST(4055673282236.0);

    B->d[0] = RCONST(2756255671327.0)/RCONST(12835298489170.0);
    B->d[1] = RCONST(-10771552573575.0)/RCONST(22201958757719.0);
    B->d[2] = RCONST(9247589265047.0)/RCONST(10645013368117.0);
    B->d[3] = RCONST(2193209047091.0)/RCONST(5459859503100.0);

    B->c[1] = RCONST(1767732205903.0)/RCONST(2027836641118.0);
    B->c[2] = RCONST(3.0)/RCONST(5.0);
    B->c[3] = RCONST(1.0);
    break;

  case(CASH_5_2_4):    /* Cash(5,2,4)-SDIRK */
    B = AllocButcherTable(5, SUNTRUE);
    B->q = 4;
    B->p = 2;
    B->A[0][0] = RCONST(0.435866521508);
    B->A[1][0] = RCONST(-1.13586652150);
    B->A[1][1] = RCONST(0.435866521508);
    B->A[2][0] = RCONST(1.08543330679);
    B->A[2][1] = RCONST(-0.721299828287);
    B->A[2][2] = RCONST(0.435866521508);
    B->A[3][0] = RCONST(0.416349501547);
    B->A[3][1] = RCONST(0.190984004184);
    B->A[3][2] = RCONST(-0.118643265417);
    B->A[3][3] = RCONST(0.435866521508);
    B->A[4][0] = RCONST(0.896869652944);
    B->A[4][1] = RCONST(0.0182725272734);
    B->A[4][2] = RCONST(-0.0845900310706);
    B->A[4][3] = RCONST(-0.266418670647);
    B->A[4][4] = RCONST(0.435866521508);

    B->b[0] = RCONST(0.896869652944);
    B->b[1] = RCONST(0.0182725272734);
    B->b[2] = RCONST(-0.0845900310706);
    B->b[3] = RCONST(-0.266418670647);
    B->b[4] = RCONST(0.435866521508);

    B->d[0] = (RCONST(-0.7)-RCONST(0.5))/(RCONST(-0.7)-RCONST(0.435866521508));
    B->d[1] = (RCONST(0.5)-RCONST(0.435866521508))/(RCONST(-0.7)-RCONST(0.435866521508));

    B->c[0] = RCONST(0.435866521508);
    B->c[1] = RCONST(-0.7);
    B->c[2] = RCONST(0.8);
    B->c[3] = RCONST(0.924556761814);
    B->c[4] = RCONST(1.0);
    break;

  case(CASH_5_3_4):    /* Cash(5,3,4)-SDIRK */
    B = AllocButcherTable(5, SUNTRUE);
    B->q = 4;
    B->p = 3;
    B->A[0][0] = RCONST(0.435866521508);
    B->A[1][0] = RCONST(-1.13586652150);
    B->A[1][1] = RCONST(0.435866521508);
    B->A[2][0] = RCONST(1.08543330679);
    B->A[2][1] = RCONST(-0.721299828287);
    B->A[2][2] = RCONST(0.435866521508);
    B->A[3][0] = RCONST(0.416349501547);
    B->A[3][1] = RCONST(0.190984004184);
    B->A[3][2] = RCONST(-0.118643265417);
    B->A[3][3] = RCONST(0.435866521508);
    B->A[4][0] = RCONST(0.896869652944);
    B->A[4][1] = RCONST(0.0182725272734);
    B->A[4][2] = RCONST(-0.0845900310706);
    B->A[4][3] = RCONST(-0.266418670647);
    B->A[4][4] = RCONST(0.435866521508);

    B->b[0] = RCONST(0.896869652944);
    B->b[1] = RCONST(0.0182725272734);
    B->b[2] = RCONST(-0.0845900310706);
    B->b[3] = RCONST(-0.266418670647);
    B->b[4] = RCONST(0.435866521508);

    B->d[0] = RCONST(0.776691932910);
    B->d[1] = RCONST(0.0297472791484);
    B->d[2] = RCONST(-0.0267440239074);
    B->d[3] = RCONST(0.220304811849);

    B->c[0] = RCONST(0.435866521508);
    B->c[1] = RCONST(-0.7);
    B->c[2] = RCONST(0.8);
    B->c[3] = RCONST(0.924556761814);
    B->c[4] = RCONST(1.0);
    break;

  case(SDIRK_5_3_4):    /* SDIRK-5-4 */
    B = AllocButcherTable(5, SUNTRUE);
    B->q = 4;
    B->p = 3;
    B->A[0][0] = RCONST(0.25);
    B->A[1][0] = RCONST(0.5);
    B->A[1][1] = RCONST(0.25);
    B->A[2][0] = RCONST(17.0)/RCONST(50.0);
    B->A[2][1] = RCONST(-1.0)/RCONST(25.0);
    B->A[2][2] = RCONST(0.25);
    B->A[3][0] = RCONST(371.0)/RCONST(1360.0);
    B->A[3][1] = RCONST(-137.0)/RCONST(2720.0);
    B->A[3][2] = RCONST(15.0)/RCONST(544.0);
    B->A[3][3] = RCONST(0.25);
    B->A[4][0] = RCONST(25.0)/RCONST(24.0);
    B->A[4][1] = RCONST(-49.0)/RCONST(48.0);
    B->A[4][2] = RCONST(125.0)/RCONST(16.0);
    B->A[4][3] = RCONST(-85.0)/RCONST(12.0);
    B->A[4][4] = RCONST(0.25);

    B->b[0] = RCONST(25.0)/RCONST(24.0);
    B->b[1] = RCONST(-49.0)/RCONST(48.0);
    B->b[2] = RCONST(125.0)/RCONST(16.0);
    B->b[3] = RCONST(-85.0)/RCONST(12.0);
    B->b[4] = RCONST(0.25);

    B->d[0] = RCONST(59.0)/RCONST(48.0);
    B->d[1] = RCONST(-17.0)/RCONST(96.0);
    B->d[2] = RCONST(225.0)/RCONST(32.0);
    B->d[3] = RCONST(-85.0)/RCONST(12.0);

    B->c[0] = RCONST(0.25);
    B->c[1] = RCONST(0.75);
    B->c[2] = RCONST(11.0)/RCONST(20.0);
    B->c[3] = RCONST(0.5);
    B->c[4] = RCONST(1.0);
    break;

  case(KVAERNO_5_3_4):    /* Kvaerno(5,3,4)-ESDIRK */
    B = AllocButcherTable(5, SUNTRUE);
    B->q = 4;
    B->p = 3;
    B->A[1][0] = RCONST(0.4358665215); 
    B->A[1][1] = RCONST(0.4358665215); 
    B->A[2][0] = RCONST(0.140737774731968);
    B->A[2][1] = RCONST(-0.108365551378832);
    B->A[2][2] = RCONST(0.4358665215);
    B->A[3][0] = RCONST(0.102399400616089);
    B->A[3][1] = RCONST(-0.376878452267324);
    B->A[3][2] = RCONST(0.838612530151233);
    B->A[3][3] = RCONST(0.4358665215);
    B->A[4][0] = RCONST(0.157024897860995);
    B->A[4][1] = RCONST(0.117330441357768);
    B->A[4][2] = RCONST(0.61667803039168);
    B->A[4][3] = RCONST(-0.326899891110444);
    B->A[4][4] = RCONST(0.4358665215);

    B->b[0] = RCONST(0.157024897860995);
    B->b[1] = RCONST(0.117330441357768);
    B->b[2] = RCONST(0.61667803039168);
    B->b[3] = RCONST(-0.326899891110444);
    B->b[4] = RCONST(0.4358665215);

    B->d[0] = RCONST(0.102399400616089);
    B->d[1] = RCONST(-0.376878452267324);
    B->d[2] = RCONST(0.838612530151233);
    B->d[3] = RCONST(0.4358665215);

    B->c[1] = RCONST(0.871733043);
    B->c[2] = RCONST(0.468238744853136);
    B->c[3] = RCONST(1.0);
    B->c[4] = RCONST(1.0);
    break;

  case(ARK436L2SA_DIRK_6_3_4):    /* ARK4(3)6L[2]SA-ESDIRK */
    B = AllocButcherTable(6, SUNTRUE);
    B->q = 4;
    B->p = 3;
    B->A[1][0] = RCONST(1.0)/RCONST(4.0);
    B->A[1][1] = RCONST(1.0)/RCONST(4.0);
    B->A[2][0] = RCONST(8611.0)/RCONST(62500.0);
    B->A[2][1] = RCONST(-1743.0)/RCONST(31250.0);
    B->A[2][2] = RCONST(1.0)/RCONST(4.0);
    B->A[3][0] = RCONST(5012029.0)/RCONST(34652500.0);
    B->A[3][1] = RCONST(-654441.0)/RCONST(2922500.0);
    B->A[3][2] = RCONST(174375.0)/RCONST(388108.0);
    B->A[3][3] = RCONST(1.0)/RCONST(4.0);
    B->A[4][0] = RCONST(15267082809.0)/RCONST(155376265600.0);
    B->A[4][1] = RCONST(-71443401.0)/RCONST(120774400.0);
    B->A[4][2] = RCONST(730878875.0)/RCONST(902184768.0);
    B->A[4][3] = RCONST(2285395.0)/RCONST(8070912.0);
    B->A[4][4] = RCONST(1.0)/RCONST(4.0);
    B->A[5][0] = RCONST(82889.0)/RCONST(524892.0);
    B->A[5][2] = RCONST(15625.0)/RCONST(83664.0);
    B->A[5][3] = RCONST(69875.0)/RCONST(102672.0);
    B->A[5][4] = RCONST(-2260.0)/RCONST(8211.0);
    B->A[5][5] = RCONST(1.0)/RCONST(4.0);

    B->b[0] = RCONST(82889.0)/RCONST(524892.0);
    B->b[2] = RCONST(15625.0)/RCONST(83664.0);
    B->b[3] = RCONST(69875.0)/RCONST(102672.0);
    B->b[4] = RCONST(-2260.0)/RCONST(8211.0);
    B->b[5] = RCONST(1.0)/RCONST(4.0);

    B->c[1] = RCONST(1.0)/RCONST(2.0);
    B->c[2] = RCONST(83.0)/RCONST(250.0);
    B->c[3] = RCONST(31.0)/RCONST(50.0);
    B->c[4] = RCONST(17.0)/RCONST(20.0);
    B->c[5] = RCONST(1.0);

    B->d[0] = RCONST(4586570599.0)/RCONST(29645900160.0);
    B->d[2] = RCONST(178811875.0)/RCONST(945068544.0);
    B->d[3] = RCONST(814220225.0)/RCONST(1159782912.0);
    B->d[4] = RCONST(-3700637.0)/RCONST(11593932.0);
    B->d[5] = RCONST(61727.0)/RCONST(225920.0);
    break;

  case(KVAERNO_7_4_5):    /* Kvaerno(7,4,5)-ESDIRK */
    B = AllocButcherTable(7, SUNTRUE);
    B->q = 5;
    B->p = 4;
    B->A[1][0] = RCONST(0.26);
    B->A[1][1] = RCONST(0.26);
    B->A[2][0] = RCONST(0.13);
    B->A[2][1] = RCONST(0.84033320996790809);
    B->A[2][2] = RCONST(0.26);
    B->A[3][0] = RCONST(0.22371961478320505);
    B->A[3][1] = RCONST(0.47675532319799699);
    B->A[3][2] = RCONST(-0.06470895363112615);
    B->A[3][3] = RCONST(0.26);
    B->A[4][0] = RCONST(0.16648564323248321);
    B->A[4][1] = RCONST(0.10450018841591720);
    B->A[4][2] = RCONST(0.03631482272098715);
    B->A[4][3] = RCONST(-0.13090704451073998);
    B->A[4][4] = RCONST(0.26);
    B->A[5][0] = RCONST(0.13855640231268224);
    B->A[5][2] = RCONST(-0.04245337201752043);
    B->A[5][3] = RCONST(0.02446657898003141);
    B->A[5][4] = RCONST(0.61943039072480676);
    B->A[5][5] = RCONST(0.26);
    B->A[6][0] = RCONST(0.13659751177640291);
    B->A[6][2] = RCONST(-0.05496908796538376);
    B->A[6][3] = RCONST(-0.04118626728321046);
    B->A[6][4] = RCONST(0.62993304899016403);
    B->A[6][5] = RCONST(0.06962479448202728);
    B->A[6][6] = RCONST(0.26);

    B->b[0] = RCONST(0.13659751177640291);
    B->b[2] = RCONST(-0.05496908796538376);
    B->b[3] = RCONST(-0.04118626728321046);
    B->b[4] = RCONST(0.62993304899016403);
    B->b[5] = RCONST(0.06962479448202728);
    B->b[6] = RCONST(0.26);

    B->d[0] = RCONST(0.13855640231268224);
    B->d[2] = RCONST(-0.04245337201752043);
    B->d[3] = RCONST(0.02446657898003141);
    B->d[4] = RCONST(0.61943039072480676);
    B->d[5] = RCONST(0.26);

    B->c[1] = RCONST(0.52);
    B->c[2] = RCONST(1.230333209967908);
    B->c[3] = RCONST(0.895765984350076);
    B->c[4] = RCONST(0.436393609858648);
    B->c[5] = RCONST(1.0);
    B->c[6] = RCONST(1.0);
    break;

  case(ARK548L2SA_DIRK_8_4_5):    /* ARK5(4)8L[2]SA-ESDIRK */
    B = AllocButcherTable(8, SUNTRUE);
    B->q = 5;
    B->p = 4;
    B->A[1][0] = RCONST(41.0)/RCONST(200.0);
    B->A[1][1] = RCONST(41.0)/RCONST(200.0);
    B->A[2][0] = RCONST(41.0)/RCONST(400.0);
    B->A[2][1] = RCONST(-567603406766.0)/RCONST(11931857230679.0);
    B->A[2][2] = RCONST(41.0)/RCONST(200.0);
    B->A[3][0] = RCONST(683785636431.0)/RCONST(9252920307686.0);
    B->A[3][2] = RCONST(-110385047103.0)/RCONST(1367015193373.0);
    B->A[3][3] = RCONST(41.0)/RCONST(200.0);
    B->A[4][0] = RCONST(3016520224154.0)/RCONST(10081342136671.0);
    B->A[4][2] = RCONST(30586259806659.0)/RCONST(12414158314087.0);
    B->A[4][3] = RCONST(-22760509404356.0)/RCONST(11113319521817.0);
    B->A[4][4] = RCONST(41.0)/RCONST(200.0);
    B->A[5][0] = RCONST(218866479029.0)/RCONST(1489978393911.0);
    B->A[5][2] = RCONST(638256894668.0)/RCONST(5436446318841.0);
    B->A[5][3] = RCONST(-1179710474555.0)/RCONST(5321154724896.0);
    B->A[5][4] = RCONST(-60928119172.0)/RCONST(8023461067671.0);
    B->A[5][5] = RCONST(41.0)/RCONST(200.0);
    B->A[6][0] = RCONST(1020004230633.0)/RCONST(5715676835656.0);
    B->A[6][2] = RCONST(25762820946817.0)/RCONST(25263940353407.0);
    B->A[6][3] = RCONST(-2161375909145.0)/RCONST(9755907335909.0);
    B->A[6][4] = RCONST(-211217309593.0)/RCONST(5846859502534.0);
    B->A[6][5] = RCONST(-4269925059573.0)/RCONST(7827059040749.0);
    B->A[6][6] = RCONST(41.0)/RCONST(200.0);
    B->A[7][0] = RCONST(-872700587467.0)/RCONST(9133579230613.0);
    B->A[7][3] = RCONST(22348218063261.0)/RCONST(9555858737531.0);
    B->A[7][4] = RCONST(-1143369518992.0)/RCONST(8141816002931.0);
    B->A[7][5] = RCONST(-39379526789629.0)/RCONST(19018526304540.0);
    B->A[7][6] = RCONST(32727382324388.0)/RCONST(42900044865799.0);
    B->A[7][7] = RCONST(41.0)/RCONST(200.0);

    B->b[0] = RCONST(-872700587467.0)/RCONST(9133579230613.0);
    B->b[3] = RCONST(22348218063261.0)/RCONST(9555858737531.0);
    B->b[4] = RCONST(-1143369518992.0)/RCONST(8141816002931.0);
    B->b[5] = RCONST(-39379526789629.0)/RCONST(19018526304540.0);
    B->b[6] = RCONST(32727382324388.0)/RCONST(42900044865799.0);
    B->b[7] = RCONST(41.0)/RCONST(200.0);

    B->d[0] = RCONST(-975461918565.0)/RCONST(9796059967033.0);
    B->d[3] = RCONST(78070527104295.0)/RCONST(32432590147079.0);
    B->d[4] = RCONST(-548382580838.0)/RCONST(3424219808633.0);
    B->d[5] = RCONST(-33438840321285.0)/RCONST(15594753105479.0);
    B->d[6] = RCONST(3629800801594.0)/RCONST(4656183773603.0);
    B->d[7] = RCONST(4035322873751.0)/RCONST(18575991585200.0);

    B->c[1] = RCONST(41.0)/RCONST(100.0);
    B->c[2] = RCONST(2935347310677.0)/RCONST(11292855782101.0);
    B->c[3] = RCONST(1426016391358.0)/RCONST(7196633302097.0);
    B->c[4] = RCONST(92.0)/RCONST(100.0);
    B->c[5] = RCONST(24.0)/RCONST(100.0);
    B->c[6] = RCONST(3.0)/RCONST(5.0);
    B->c[7] = RCONST(1.0);
    break;

  default:

    arkProcessError(NULL, ARK_ILL_INPUT, "ARKODE", 
                    "ARKodeLoadButcherTable", "Unknown Butcher table");
    return(NULL);

  }

  return(B);
}


/*---------------------------------------------------------------
  EOF
  ---------------------------------------------------------------*/
