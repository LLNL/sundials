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
#include <arkode/arkode.h>
#include <arkode/arkode_sprk.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

ARKodeSPRKMem ARKodeSymplecticEuler() {
  ARKodeSPRKMem sprk_mem = ARKodeSPRKMem_Alloc(1);
  sprk_mem->q = 1;
  sprk_mem->stages = 1;
  sprk_mem->b[0] = SUN_RCONST(1.0);
  sprk_mem->a[0] = SUN_RCONST(1.0);
  return sprk_mem;
}

/* 
  The following methods are from:

  J Candy, W Rozmus, A symplectic integration algorithm for separable Hamiltonian functions,
  Journal of Computational Physics, Volume 92, Issue 1, 1991, Pages 230-256, ISSN 0021-9991,
  https://doi.org/10.1016/0021-9991(91)90299-Z.
 */

ARKodeSPRKMem ARKodeSymplecticLeapfrog2() {
  ARKodeSPRKMem sprk_mem = ARKodeSPRKMem_Alloc(2);
  sprk_mem->q = 2;
  sprk_mem->stages = 2;
  sprk_mem->a[0] = SUN_RCONST(0.5); 
  sprk_mem->a[1] = SUN_RCONST(0.5);
  sprk_mem->b[0] = SUN_RCONST(0.0);
  sprk_mem->b[1] = SUN_RCONST(1.0);
  return sprk_mem;
}

ARKodeSPRKMem ARKodeSymplecticPseudoLeapfrog2() {
  ARKodeSPRKMem sprk_mem = ARKodeSPRKMem_Alloc(2);
  sprk_mem->q = 2;
  sprk_mem->stages = 2;
  sprk_mem->a[0] = SUN_RCONST(1.0); 
  sprk_mem->a[1] = SUN_RCONST(0.0);
  sprk_mem->b[0] = SUN_RCONST(0.5);
  sprk_mem->b[1] = SUN_RCONST(0.5);
  return sprk_mem;
}

ARKodeSPRKMem ARKodeSymplecticCandyRozmus4() {
  ARKodeSPRKMem sprk_mem = ARKodeSPRKMem_Alloc(4);
  sprk_mem->q = 4;
  sprk_mem->stages = 4;
  sprk_mem->a[0] = ( SUN_RCONST(2.0) + SUNRpowerR(SUN_RCONST(2.0), SUN_RCONST(1.0)/SUN_RCONST(3.0)) 
                     + SUNRpowerR(SUN_RCONST(2.0), -SUN_RCONST(1.0)/SUN_RCONST(3.0)) ) / SUN_RCONST(6.0);
  sprk_mem->a[1] = ( SUN_RCONST(1.0) - SUNRpowerR(SUN_RCONST(2.0), SUN_RCONST(1.0)/SUN_RCONST(3.0)) 
                     - SUNRpowerR(SUN_RCONST(2.0), -SUN_RCONST(1.0)/SUN_RCONST(3.0)) ) / SUN_RCONST(6.0);
  sprk_mem->a[2] = sprk_mem->a[1];
  sprk_mem->a[3] = sprk_mem->a[0];
  sprk_mem->b[0] = SUN_RCONST(0.0);
  sprk_mem->b[1] = SUN_RCONST(1.0) / 
                   ( SUN_RCONST(2.0) - SUNRpowerR(SUN_RCONST(2.0), SUN_RCONST(1.0)/SUN_RCONST(3.0)) );
  sprk_mem->b[2] = SUN_RCONST(1.0) / 
                   ( SUN_RCONST(1.0) - SUNRpowerR(SUN_RCONST(2.0), SUN_RCONST(2.0)/SUN_RCONST(3.0)) );
  sprk_mem->b[3] = sprk_mem->b[1];
  return sprk_mem;
}

/* 
  The following methods are from: 

  Ruth, R. D. (1983). A CANONICAL INTEGRATION TECHNIQUE.
  IEEE Transactions on Nuclear Science, 30(4).
  https://accelconf.web.cern.ch/p83/PDF/PAC1983_2669.PDF 
 */

ARKodeSPRKMem ARKodeSymplecticRuth3() {
  ARKodeSPRKMem sprk_mem = ARKodeSPRKMem_Alloc(3);
  sprk_mem->q = 3;
  sprk_mem->stages = 3;
  sprk_mem->a[0] = SUN_RCONST(2.0)/SUN_RCONST(3.0);
  sprk_mem->a[1] = -SUN_RCONST(2.0)/SUN_RCONST(3.0);
  sprk_mem->a[2] = SUN_RCONST(1.0);
  sprk_mem->b[0] = SUN_RCONST(7.0)/SUN_RCONST(24.0);
  sprk_mem->b[1] = SUN_RCONST(3.0)/SUN_RCONST(4.0);
  sprk_mem->b[2] = -SUN_RCONST(1.0)/SUN_RCONST(24.0);
  return sprk_mem;
}

/* 
  The following methods are from: 

  McLachlan, R.I., Atela, P.: The accuracy of symplectic integrators. 
  Nonlinearity. 5, 541–562 (1992). https://doi.org/10.1088/0951-7715/5/2/011
 */

ARKodeSPRKMem ARKodeSymplecticMcLachlan2() {
  ARKodeSPRKMem sprk_mem = ARKodeSPRKMem_Alloc(2);
  sprk_mem->q = 2;
  sprk_mem->stages = 2;
  sprk_mem->a[1] = SUN_RCONST(1.0) - (SUN_RCONST(1.0) / SUN_RCONST(2.0)) * SUNRsqrt(2.0);
  sprk_mem->a[0] = SUN_RCONST(1.0) - sprk_mem->a[1];
  sprk_mem->b[1] = SUN_RCONST(1.0) / (SUN_RCONST(2.0) * (SUN_RCONST(1.0) - sprk_mem->a[1]));
  sprk_mem->b[0] = SUN_RCONST(1.0) - sprk_mem->b[1];
  return sprk_mem;
}

ARKodeSPRKMem ARKodeSymplecticMcLachlan3() {
  ARKodeSPRKMem sprk_mem = ARKodeSPRKMem_Alloc(3);
  sprk_mem->q = 3;
  sprk_mem->stages = 3;
  sprk_mem->a[0] = SUN_RCONST(0.919661523017399857);
  sprk_mem->a[1] = SUN_RCONST(0.25)/sprk_mem->a[0] - sprk_mem->a[0]/SUN_RCONST(2.0);
  sprk_mem->a[2] = SUN_RCONST(1.0) - sprk_mem->a[0] - sprk_mem->a[1];
  sprk_mem->b[0] = sprk_mem->a[2];
  sprk_mem->b[1] = sprk_mem->a[1];
  sprk_mem->b[2] = sprk_mem->a[0];
  return sprk_mem;
}

ARKodeSPRKMem ARKodeSymplecticMcLachlan4() {
  ARKodeSPRKMem sprk_mem = ARKodeSPRKMem_Alloc(4);
  sprk_mem->q = 4;
  sprk_mem->stages = 4;
  sprk_mem->a[0] = SUN_RCONST(0.515352837431122936); 
  sprk_mem->a[1] = -SUN_RCONST(0.085782019412973646);
  sprk_mem->a[2] = SUN_RCONST(0.441583023616466524);
  sprk_mem->a[3] = SUN_RCONST(0.128846158365384185);
  sprk_mem->b[0] = SUN_RCONST(0.134496199277431089);
  sprk_mem->b[1] = -SUN_RCONST(0.224819803079420806);
  sprk_mem->b[2] = SUN_RCONST(0.756320000515668291);
  sprk_mem->b[3] = SUN_RCONST(0.33400360328632142);
  return sprk_mem;
}

ARKodeSPRKMem ARKodeSymplecticMcLachlan5() {
  ARKodeSPRKMem sprk_mem = ARKodeSPRKMem_Alloc(6);
  sprk_mem->q = 5;
  sprk_mem->stages = 6;
  sprk_mem->a[0] = SUN_RCONST(0.339839625839110000);
  sprk_mem->a[1] = -SUN_RCONST(0.088601336903027329);
  sprk_mem->a[2] = SUN_RCONST(0.5858564768259621188);
  sprk_mem->a[3] = -SUN_RCONST(0.603039356536491888);
  sprk_mem->a[4] = SUN_RCONST(0.3235807965546976394);
  sprk_mem->a[5] = SUN_RCONST(0.4423637942197494587);
  sprk_mem->b[0] = SUN_RCONST(0.1193900292875672758);
  sprk_mem->b[1] = SUN_RCONST(0.6989273703824752308);
  sprk_mem->b[2] = -SUN_RCONST(0.1713123582716007754);
  sprk_mem->b[3] = SUN_RCONST(0.4012695022513534480);
  sprk_mem->b[4] = SUN_RCONST(0.0107050818482359840);
  sprk_mem->b[5] = -SUN_RCONST(0.0589796254980311632);
  return sprk_mem;
}

/* 
  The following methods are from: 

  Yoshida, H.: Construction of higher order symplectic integrators. 
  Phys Lett A. 150, 262–268 (1990). 
  https://doi.org/10.1016/0375-9601(90)90092-3

 */

ARKodeSPRKMem ARKodeSymplecticYoshida6() {
  ARKodeSPRKMem sprk_mem = ARKodeSPRKMem_Alloc(8);
  sprk_mem->q = 6;
  sprk_mem->stages = 8;
  sprk_mem->a[0] = SUN_RCONST(0.78451361047755726382);
  sprk_mem->a[1] = SUN_RCONST(0.23557321335935813368);
  sprk_mem->a[2] = -SUN_RCONST(1.17767998417887100695);
  sprk_mem->a[3] = SUN_RCONST(1.3151863206839);
  sprk_mem->a[4] = sprk_mem->a[2];
  sprk_mem->a[5] = sprk_mem->a[1];
  sprk_mem->a[6] = sprk_mem->a[0];
  sprk_mem->a[7] = SUN_RCONST(0.0);
  sprk_mem->b[0] = sprk_mem->a[0] / SUN_RCONST(2.0);
  sprk_mem->b[1] = (sprk_mem->a[0] + sprk_mem->a[1]) / SUN_RCONST(2.0);
  sprk_mem->b[2] = (sprk_mem->a[1] + sprk_mem->a[2]) / SUN_RCONST(2.0);
  sprk_mem->b[3] = (sprk_mem->a[2] + sprk_mem->a[3]) / SUN_RCONST(2.0);
  sprk_mem->b[4] = sprk_mem->b[3];
  sprk_mem->b[5] = sprk_mem->b[2];
  sprk_mem->b[6] = sprk_mem->b[1];
  sprk_mem->b[7] = sprk_mem->b[0];
  return sprk_mem;
}

/* 
  The following methods are from: 

  McLachlan, R.I.: On the Numerical Integration of Ordinary Differential Equations
  by Symmetric Composition Methods. Siam J Sci Comput. 16, 151–168 (1995).
  https://doi.org/10.1137/0916010

 */

ARKodeSPRKMem ARKodeSymplecticMcLachlan8() {
  ARKodeSPRKMem sprk_mem = ARKodeSPRKMem_Alloc(16);
  sprk_mem->q = 8;
  sprk_mem->stages = 16;
  sprk_mem->a[0] = SUN_RCONST(0.74167036435061295344822780);
  sprk_mem->a[1] = -SUN_RCONST(0.40910082580003159399730010);
  sprk_mem->a[2] = SUN_RCONST(0.19075471029623837995387626);
  sprk_mem->a[3] = -SUN_RCONST(0.57386247111608226665638773);
  sprk_mem->a[4] = SUN_RCONST(0.29906418130365592384446354);
  sprk_mem->a[5] = SUN_RCONST(0.33462491824529818378495798);
  sprk_mem->a[6] = SUN_RCONST(0.31529309239676659663205666);
  sprk_mem->a[7] = -SUN_RCONST(0.79688793935291635401978884);
  sprk_mem->a[8] = sprk_mem->a[6];
  sprk_mem->a[9] = sprk_mem->a[5];
  sprk_mem->a[10] = sprk_mem->a[4];
  sprk_mem->a[11] = sprk_mem->a[3];
  sprk_mem->a[12] = sprk_mem->a[2];
  sprk_mem->a[13] = sprk_mem->a[1];
  sprk_mem->a[14] = sprk_mem->a[0];
  sprk_mem->a[15] = SUN_RCONST(0.0);
  sprk_mem->b[0] = sprk_mem->a[0] / SUN_RCONST(2.0);
  sprk_mem->b[1] = (sprk_mem->a[0] + sprk_mem->a[1]) / SUN_RCONST(2.0);
  sprk_mem->b[2] = (sprk_mem->a[1] + sprk_mem->a[2]) / SUN_RCONST(2.0);
  sprk_mem->b[3] = (sprk_mem->a[2] + sprk_mem->a[3]) / SUN_RCONST(2.0);
  sprk_mem->b[4] = (sprk_mem->a[3] + sprk_mem->a[4]) / SUN_RCONST(2.0);
  sprk_mem->b[5] = (sprk_mem->a[4] + sprk_mem->a[5]) / SUN_RCONST(2.0);
  sprk_mem->b[6] = (sprk_mem->a[5] + sprk_mem->a[6]) / SUN_RCONST(2.0);
  sprk_mem->b[7] = (sprk_mem->a[6] + sprk_mem->a[7]) / SUN_RCONST(2.0);
  sprk_mem->b[8] = sprk_mem->b[7];
  sprk_mem->b[9] = sprk_mem->b[6];
  sprk_mem->b[10] = sprk_mem->b[5];
  sprk_mem->b[11] = sprk_mem->b[4];
  sprk_mem->b[12] = sprk_mem->b[3];
  sprk_mem->b[13] = sprk_mem->b[2];
  sprk_mem->b[14] = sprk_mem->b[1];
  sprk_mem->b[15] = sprk_mem->b[0];
  return sprk_mem;
}

/* 
  The following methods are from: 

  Sofroniou, M., Spaletta, G.: Derivation of symmetric composition constants for 
  symmetric integrators. Optim Methods Softw. 20, 597–613 (2005). 
  https://doi.org/10.1080/10556780500140664
  
 */

ARKodeSPRKMem ARKodeSymplecticSofroniou10() {
  ARKodeSPRKMem sprk_mem = ARKodeSPRKMem_Alloc(36);
  sprk_mem->q = 10;
  sprk_mem->stages = 36;

  sprk_mem->a[0] = SUN_RCONST(0.078795722521686419263907679337684);
  sprk_mem->a[1] = SUN_RCONST(0.31309610341510852776481247192647);
  sprk_mem->a[2] = SUN_RCONST(0.027918383235078066109520273275299);
  sprk_mem->a[3] = -SUN_RCONST(0.22959284159390709415121339679655);
  sprk_mem->a[4] = SUN_RCONST(0.13096206107716486317465685927961);
  sprk_mem->a[5] = -SUN_RCONST(0.26973340565451071434460973222411);
  sprk_mem->a[6] = SUN_RCONST(0.074973343155891435666137105641410);
  sprk_mem->a[7] = SUN_RCONST(0.11199342399981020488957508073640);
  sprk_mem->a[8] = SUN_RCONST(0.36613344954622675119314812353150);
  sprk_mem->a[9] = -SUN_RCONST(0.39910563013603589787862981058340);
  sprk_mem->a[10] = SUN_RCONST(0.10308739852747107731580277001372);
  sprk_mem->a[11] = SUN_RCONST(0.41143087395589023782070411897608);
  sprk_mem->a[12] = -SUN_RCONST(0.0048663605831352617621956593099771);
  sprk_mem->a[13] = -SUN_RCONST(0.39203335370863990644808193642610);
  sprk_mem->a[14] = SUN_RCONST(0.051942502962449647037182904015976);
  sprk_mem->a[15] = SUN_RCONST(0.050665090759924496335874344156866);
  sprk_mem->a[16] = SUN_RCONST(0.049674370639729879054568800279461);
  sprk_mem->a[17] = SUN_RCONST(0.049317735759594537917680008339338);
  sprk_mem->a[18] = sprk_mem->a[16];
  sprk_mem->a[19] = sprk_mem->a[15];
  sprk_mem->a[20] = sprk_mem->a[14];
  sprk_mem->a[21] = sprk_mem->a[13];
  sprk_mem->a[22] = sprk_mem->a[12];
  sprk_mem->a[23] = sprk_mem->a[11];
  sprk_mem->a[24] = sprk_mem->a[10];
  sprk_mem->a[25] = sprk_mem->a[9];
  sprk_mem->a[26] = sprk_mem->a[8];
  sprk_mem->a[27] = sprk_mem->a[7];
  sprk_mem->a[28] = sprk_mem->a[6];
  sprk_mem->a[29] = sprk_mem->a[5];
  sprk_mem->a[30] = sprk_mem->a[4];
  sprk_mem->a[31] = sprk_mem->a[3];
  sprk_mem->a[32] = sprk_mem->a[2];
  sprk_mem->a[33] = sprk_mem->a[1];
  sprk_mem->a[34] = sprk_mem->a[0];
  sprk_mem->a[35] = SUN_RCONST(0.0);
  sprk_mem->b[0] = sprk_mem->a[0] / SUN_RCONST(2.0);
  sprk_mem->b[1] = (sprk_mem->a[0] + sprk_mem->a[1]) / SUN_RCONST(2.0);
  sprk_mem->b[2] = (sprk_mem->a[1] + sprk_mem->a[2]) / SUN_RCONST(2.0);
  sprk_mem->b[3] = (sprk_mem->a[2] + sprk_mem->a[3]) / SUN_RCONST(2.0);
  sprk_mem->b[4] = (sprk_mem->a[3] + sprk_mem->a[4]) / SUN_RCONST(2.0);
  sprk_mem->b[5] = (sprk_mem->a[4] + sprk_mem->a[5]) / SUN_RCONST(2.0);
  sprk_mem->b[6] = (sprk_mem->a[5] + sprk_mem->a[6]) / SUN_RCONST(2.0);
  sprk_mem->b[7] = (sprk_mem->a[6] + sprk_mem->a[7]) / SUN_RCONST(2.0);
  sprk_mem->b[8] = (sprk_mem->a[7] + sprk_mem->a[8]) / SUN_RCONST(2.0);
  sprk_mem->b[9] = (sprk_mem->a[8] + sprk_mem->a[9]) / SUN_RCONST(2.0);
  sprk_mem->b[10] = (sprk_mem->a[9] + sprk_mem->a[10]) / SUN_RCONST(2.0);
  sprk_mem->b[11] = (sprk_mem->a[10] + sprk_mem->a[11]) / SUN_RCONST(2.0);
  sprk_mem->b[12] = (sprk_mem->a[11] + sprk_mem->a[12]) / SUN_RCONST(2.0);
  sprk_mem->b[13] = (sprk_mem->a[12] + sprk_mem->a[13]) / SUN_RCONST(2.0);
  sprk_mem->b[14] = (sprk_mem->a[13] + sprk_mem->a[14]) / SUN_RCONST(2.0);
  sprk_mem->b[15] = (sprk_mem->a[14] + sprk_mem->a[15]) / SUN_RCONST(2.0);
  sprk_mem->b[16] = (sprk_mem->a[15] + sprk_mem->a[16]) / SUN_RCONST(2.0);
  sprk_mem->b[17] = (sprk_mem->a[16] + sprk_mem->a[17]) / SUN_RCONST(2.0);
  sprk_mem->b[18] = sprk_mem->b[17];
  sprk_mem->b[19] = sprk_mem->b[16];
  sprk_mem->b[20] = sprk_mem->b[15];
  sprk_mem->b[21] = sprk_mem->b[14];
  sprk_mem->b[22] = sprk_mem->b[13];
  sprk_mem->b[23] = sprk_mem->b[12];
  sprk_mem->b[24] = sprk_mem->b[11];
  sprk_mem->b[25] = sprk_mem->b[10];
  sprk_mem->b[26] = sprk_mem->b[9];
  sprk_mem->b[27] = sprk_mem->b[8];
  sprk_mem->b[28] = sprk_mem->b[7];
  sprk_mem->b[29] = sprk_mem->b[6];
  sprk_mem->b[30] = sprk_mem->b[5];
  sprk_mem->b[31] = sprk_mem->b[4];
  sprk_mem->b[32] = sprk_mem->b[3];
  sprk_mem->b[33] = sprk_mem->b[2];
  sprk_mem->b[34] = sprk_mem->b[1];
  sprk_mem->b[35] = sprk_mem->b[0];

  return sprk_mem;
}

ARKodeSPRKMem ARKodeSPRKMem_Alloc(int stages)
{
  ARKodeSPRKMem sprk_mem;

  sprk_mem = (ARKodeSPRKMem) malloc(sizeof(struct ARKodeSPRKMem_s));

  sprk_mem->q = 0;
  sprk_mem->stages = stages;
  sprk_mem->b = (sunrealtype*) malloc(stages*sizeof(sunrealtype));
  sprk_mem->a = (sunrealtype*) malloc(stages*sizeof(sunrealtype));

  return sprk_mem;
}

ARKodeSPRKMem ARKodeSPRKMem_Load(ARKODE_SPRKMethodID id)
{
  switch(id) {
    case ARKODE_SYMPLECTIC_EULER_1:
      return ARKodeSymplecticEuler();
    case ARKODE_SYMPLECTIC_LEAPFROG_2:
      return ARKodeSymplecticLeapfrog2();
    case ARKODE_SYMPLECTIC_RUTH_3:
      return ARKodeSymplecticRuth3();
    case ARKODE_SYMPLECTIC_MCLACHLAN_2:
      return ARKodeSymplecticMcLachlan2();
    case ARKODE_SYMPLECTIC_MCLACHLAN_3:
      return ARKodeSymplecticMcLachlan3();
    case ARKODE_SYMPLECTIC_MCLACHLAN_4:
      return ARKodeSymplecticMcLachlan4();
    case ARKODE_SYMPLECTIC_CANDY_ROZMUS_4:
      return ARKodeSymplecticCandyRozmus4();
    case ARKODE_SYMPLECTIC_MCLACHLAN_5:
      return ARKodeSymplecticMcLachlan5();
    case ARKODE_SYMPLECTIC_YOSHIDA_6:
      return ARKodeSymplecticYoshida6();
    case ARKODE_SYMPLECTIC_MCLACHLAN_8:
      return ARKodeSymplecticMcLachlan8();
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
    sprk_mem->a[i] = that_sprk_mem->a[i];
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
    free(sprk_mem->a);
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
    B->b[i] = sprk_mem->a[i];
    for (j = 0; j < i; ++j)
    {
      B->A[i][j] = sprk_mem->a[j];
    }
  }

  /* Time weights: C_j = sum_{i=0}^{j-1} B_i */
  for (j = 0; j < sprk_mem->stages; ++j) {
    for (i = 0; i < j; ++i) {
      b->c[j] += sprk_mem->a[i];
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
