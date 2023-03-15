/*---------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *---------------------------------------------------------------
 * 
 *--------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <arkode/arkode_sprk.h>
#include <arkode/arkode.h>
#include <arkode/arkode_butcher.h>
#include <sundials/sundials_types.h>

ARKodeSPRKMem ARKodeSymplecticEuler() {
  ARKodeSPRKMem sprk_mem = ARKodeSPRKMem_Alloc(1);
  sprk_mem->q = 1;
  sprk_mem->stages = 1;
  sprk_mem->b[0] = SUN_RCONST(1.0);
  sprk_mem->B[0] = SUN_RCONST(1.0);
  return sprk_mem;
}

ARKodeSPRKMem ARKodeSymplecticLeapfrog() {
  ARKodeSPRKMem sprk_mem = ARKodeSPRKMem_Alloc(2);
  sprk_mem->q = 2;
  sprk_mem->stages = 2;
  sprk_mem->b[0] = SUN_RCONST(0.5);
  sprk_mem->b[1] = SUN_RCONST(0.5);
  sprk_mem->B[0] = SUN_RCONST(1.0); 
  sprk_mem->B[1] = SUN_RCONST(0.0);
  return sprk_mem;
}

ARKodeSPRKMem ARKodeSymplecticRuth3() {
  ARKodeSPRKMem sprk_mem = ARKodeSPRKMem_Alloc(3);
  sprk_mem->q = 3;
  sprk_mem->stages = 3;
  sprk_mem->b[0] = SUN_RCONST(7.0)/SUN_RCONST(24.0);
  sprk_mem->b[1] = SUN_RCONST(3.0)/SUN_RCONST(4.0);
  sprk_mem->b[2] = -SUN_RCONST(1.0)/SUN_RCONST(24.0);
  sprk_mem->B[0] = SUN_RCONST(2.0)/SUN_RCONST(3.0);
  sprk_mem->B[1] = -SUN_RCONST(2.0)/SUN_RCONST(3.0);
  sprk_mem->B[2] = SUN_RCONST(1.0);
  return sprk_mem;
}

ARKodeSPRKMem ARKodeSymplecticMcLauchlan4() {
  ARKodeSPRKMem sprk_mem = ARKodeSPRKMem_Alloc(4);
  sprk_mem->q = 4;
  sprk_mem->stages = 4;
  sprk_mem->b[0] = SUN_RCONST(0.134496199277431089);
  sprk_mem->b[1] = -SUN_RCONST(0.224819803079420806);
  sprk_mem->b[2] = SUN_RCONST(0.756320000515668291);
  sprk_mem->b[3] = SUN_RCONST(0.33400360328632142);
  sprk_mem->B[0] = SUN_RCONST(0.515352837431122936); 
  sprk_mem->B[1] = -SUN_RCONST(0.085782019412973646);
  sprk_mem->B[2] = SUN_RCONST(0.441583023616466524);
  sprk_mem->B[3] = SUN_RCONST(0.128846158365384185);
  return sprk_mem;
}

ARKodeSPRKMem ARKodeSPRKMem_Alloc(int stages)
{
  ARKodeSPRKMem sprk_mem;

  sprk_mem = (ARKodeSPRKMem) malloc(sizeof(struct ARKodeSPRKMem_s));

  sprk_mem->q = 0;
  sprk_mem->stages = stages;
  sprk_mem->b = (sunrealtype*) malloc(stages*sizeof(sunrealtype));
  sprk_mem->B = (sunrealtype*) malloc(stages*sizeof(sunrealtype));

  return sprk_mem;
}

ARKodeSPRKMem ARKodeSPRKMem_Load(ARKODE_SPRKMethodID id)
{
  switch(id) {
    case ARKODE_SYMPLECTIC_EULER_1:
      return ARKodeSymplecticEuler();
    case ARKODE_SYMPLECTIC_LEAPFROG_2:
      return ARKodeSymplecticLeapfrog();
    case ARKODE_SYMPLECTIC_RUTH_3:
      return ARKodeSymplecticRuth3();
    case ARKODE_SYMPLECTIC_MCLAUCHLAN_4:
      return ARKodeSymplecticMcLauchlan4();
    default:
      return NULL;
  }
}

ARKodeSPRKMem ARKodeSPRKMem_Copy(ARKodeSPRKMem that_sprk_mem)
{
  int i;
  ARKodeSPRKMem sprk_mem;

  sprk_mem = ARKodeSPRKMem_Alloc(that_sprk_mem->stages);

  sprk_mem->q = that_sprk_mem->q;

  for (i = 0; i < sprk_mem->stages; ++i)
  {
    sprk_mem->b[i] = that_sprk_mem->b[i];
    sprk_mem->B[i] = that_sprk_mem->B[i];
  }

  return sprk_mem;
}

void ARKodeSPRKMem_Space(ARKodeSPRKMem sprk_mem, sunindextype *liw, sunindextype *lrw)
{
  *liw = 2;
  *lrw = sprk_mem->stages * 2;
  return;
}

void ARKodeSPRKMem_Free(ARKodeSPRKMem sprk_mem)
{
  if (sprk_mem)
  {
    free(sprk_mem->b);
    free(sprk_mem->B);
    free(sprk_mem);
  }
  return;
}

int ARKodeSPRKMem_ToButcher(ARKodeSPRKMem sprk_mem, ARKodeButcherTable* b_ptr, ARKodeButcherTable* B_ptr)
{
  int i, j;
  ARKodeButcherTable b, B;

  b = *b_ptr;
  B = *B_ptr;

  b = ARKodeButcherTable_Alloc(sprk_mem->stages, SUNFALSE);
  B = ARKodeButcherTable_Alloc(sprk_mem->stages, SUNFALSE);

  /* DIRK table */
  for (i = 0; i < sprk_mem->stages; ++i)
  {
    b->b[i] = sprk_mem->b[i];
    for (j = 0; j <= i; ++j)
    {
      b->A[i][j] = sprk_mem->b[j];
    }
  }

  /* Explicit table */
  for (i = 0; i < sprk_mem->stages; ++i)
  {
    B->b[i] = sprk_mem->B[i];
    for (j = 0; j < i; ++j)
    {
      B->A[i][j] = sprk_mem->B[j];
    }
  }

  /* Time weights: C_j = sum_{i=0}^{j-1} B_i */
  for (j = 0; j < sprk_mem->stages; ++j) {
    for (i = 0; i < j; ++i) {
      b->c[j] += sprk_mem->B[i];
    }
  }

  /* Set method order */
  b->q = sprk_mem->q;
  B->q = sprk_mem->q;

  /* No embedding, so set embedding order to 0 */
  b->p = 0;
  B->p = 0;

  return ARK_SUCCESS;
}
