/*---------------------------------------------------------------
 * Programmer(s): Cody J. Balos @ LLNL
 *---------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2023, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 *--------------------------------------------------------------*/

#include <arkode/arkode.h>
#include <arkode/arkode_sprk.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sundials/sundials_math.h>
#include <sundials/sundials_types.h>

ARKodeSPRKStorage ARKodeSymplecticEuler()
{
  ARKodeSPRKStorage sprk_storage = ARKodeSPRKStorage_Alloc(1);
  if (!sprk_storage) { return NULL; }
  sprk_storage->q       = 1;
  sprk_storage->stages  = 1;
  sprk_storage->a[0]    = SUN_RCONST(1.0);
  sprk_storage->ahat[0] = SUN_RCONST(1.0);
  return sprk_storage;
}

/*
  The following methods are from:

  J Candy, W Rozmus, A symplectic integration algorithm for separable
  Hamiltonian functions, Journal of Computational Physics, Volume 92, Issue 1,
  1991, Pages 230-256, ISSN 0021-9991,
  https://doi.org/10.1016/0021-9991(91)90299-Z.
 */

ARKodeSPRKStorage ARKodeSymplecticLeapfrog2()
{
  ARKodeSPRKStorage sprk_storage = ARKodeSPRKStorage_Alloc(2);
  if (!sprk_storage) { return NULL; }
  sprk_storage->q       = 2;
  sprk_storage->stages  = 2;
  sprk_storage->a[0]    = SUN_RCONST(0.5);
  sprk_storage->a[1]    = SUN_RCONST(0.5);
  sprk_storage->ahat[0] = SUN_RCONST(0.0);
  sprk_storage->ahat[1] = SUN_RCONST(1.0);
  return sprk_storage;
}

ARKodeSPRKStorage ARKodeSymplecticPseudoLeapfrog2()
{
  ARKodeSPRKStorage sprk_storage = ARKodeSPRKStorage_Alloc(2);
  if (!sprk_storage) { return NULL; }
  sprk_storage->q       = 2;
  sprk_storage->stages  = 2;
  sprk_storage->a[0]    = SUN_RCONST(1.0);
  sprk_storage->a[1]    = SUN_RCONST(0.0);
  sprk_storage->ahat[0] = SUN_RCONST(0.5);
  sprk_storage->ahat[1] = SUN_RCONST(0.5);
  return sprk_storage;
}

ARKodeSPRKStorage ARKodeSymplecticCandyRozmus4()
{
  ARKodeSPRKStorage sprk_storage = ARKodeSPRKStorage_Alloc(4);
  if (!sprk_storage) { return NULL; }
  sprk_storage->q      = 4;
  sprk_storage->stages = 4;
  sprk_storage->a[0] =
    (SUN_RCONST(2.0) +
     SUNRpowerR(SUN_RCONST(2.0), SUN_RCONST(1.0) / SUN_RCONST(3.0)) +
     SUNRpowerR(SUN_RCONST(2.0), -SUN_RCONST(1.0) / SUN_RCONST(3.0))) /
    SUN_RCONST(6.0);
  sprk_storage->a[1] =
    (SUN_RCONST(1.0) -
     SUNRpowerR(SUN_RCONST(2.0), SUN_RCONST(1.0) / SUN_RCONST(3.0)) -
     SUNRpowerR(SUN_RCONST(2.0), -SUN_RCONST(1.0) / SUN_RCONST(3.0))) /
    SUN_RCONST(6.0);
  sprk_storage->a[2]    = sprk_storage->a[1];
  sprk_storage->a[3]    = sprk_storage->a[0];
  sprk_storage->ahat[0] = SUN_RCONST(0.0);
  sprk_storage->ahat[1] =
    SUN_RCONST(1.0) /
    (SUN_RCONST(2.0) -
     SUNRpowerR(SUN_RCONST(2.0), SUN_RCONST(1.0) / SUN_RCONST(3.0)));
  sprk_storage->ahat[2] =
    SUN_RCONST(1.0) /
    (SUN_RCONST(1.0) -
     SUNRpowerR(SUN_RCONST(2.0), SUN_RCONST(2.0) / SUN_RCONST(3.0)));
  sprk_storage->ahat[3] = sprk_storage->ahat[1];
  return sprk_storage;
}

/*
  The following methods are from:

  Ruth, R. D. (1983). A CANONICAL INTEGRATION TECHNIQUE.
  IEEE Transactions on Nuclear Science, 30(4).
  https://accelconf.web.cern.ch/p83/PDF/PAC1983_2669.PDF
 */

ARKodeSPRKStorage ARKodeSymplecticRuth3()
{
  ARKodeSPRKStorage sprk_storage = ARKodeSPRKStorage_Alloc(3);
  if (!sprk_storage) { return NULL; }
  sprk_storage->q       = 3;
  sprk_storage->stages  = 3;
  sprk_storage->a[0]    = SUN_RCONST(2.0) / SUN_RCONST(3.0);
  sprk_storage->a[1]    = -SUN_RCONST(2.0) / SUN_RCONST(3.0);
  sprk_storage->a[2]    = SUN_RCONST(1.0);
  sprk_storage->ahat[0] = SUN_RCONST(7.0) / SUN_RCONST(24.0);
  sprk_storage->ahat[1] = SUN_RCONST(3.0) / SUN_RCONST(4.0);
  sprk_storage->ahat[2] = -SUN_RCONST(1.0) / SUN_RCONST(24.0);
  return sprk_storage;
}

/*
  The following methods are from:

  McLachlan, R.I., Atela, P.: The accuracy of symplectic integrators.
  Nonlinearity. 5, 541–562 (1992). https://doi.org/10.1088/0951-7715/5/2/011
 */

ARKodeSPRKStorage ARKodeSymplecticMcLachlan2()
{
  ARKodeSPRKStorage sprk_storage = ARKodeSPRKStorage_Alloc(2);
  if (!sprk_storage) { return NULL; }
  sprk_storage->q      = 2;
  sprk_storage->stages = 2;
  sprk_storage->a[1]   = SUN_RCONST(1.0) -
                       (SUN_RCONST(1.0) / SUN_RCONST(2.0)) * SUNRsqrt(2.0);
  sprk_storage->a[0] = SUN_RCONST(1.0) - sprk_storage->a[1];
  sprk_storage->ahat[1] =
    SUN_RCONST(1.0) / (SUN_RCONST(2.0) * (SUN_RCONST(1.0) - sprk_storage->a[1]));
  sprk_storage->ahat[0] = SUN_RCONST(1.0) - sprk_storage->ahat[1];
  return sprk_storage;
}

ARKodeSPRKStorage ARKodeSymplecticMcLachlan3()
{
  sunrealtype w                  = 0.0;
  sunrealtype y                  = 0.0;
  sunrealtype z                  = 0.0;
  ARKodeSPRKStorage sprk_storage = ARKodeSPRKStorage_Alloc(3);
  if (!sprk_storage) { return NULL; }

  sprk_storage->q      = 3;
  sprk_storage->stages = 3;

  z = -SUNRpowerR((SUN_RCONST(2.0) / SUN_RCONST(27.0)) -
                    SUN_RCONST(1.0) / (SUN_RCONST(9.0) * SUNRsqrt(3.0)),
                  SUN_RCONST(1.0) / SUN_RCONST(3.0));
  w = -SUN_RCONST(2.0) / SUN_RCONST(3.0) +
      SUN_RCONST(1.0) / (SUN_RCONST(9.0) * z) + z;
  y                  = (SUN_RCONST(1.0) + w * w) / SUN_RCONST(4.0);
  sprk_storage->a[0] = SUNRsqrt(SUN_RCONST(1.0) / (SUN_RCONST(9.0) * y) -
                                w / SUN_RCONST(2.0) + SUNRsqrt(y)) -
                       SUN_RCONST(1.0) / (SUN_RCONST(3.0) * SUNRsqrt(y));
  sprk_storage->a[1] = SUN_RCONST(0.25) / sprk_storage->a[0] -
                       sprk_storage->a[0] / SUN_RCONST(2.0);
  sprk_storage->a[2] = SUN_RCONST(1.0) - sprk_storage->a[0] - sprk_storage->a[1];
  sprk_storage->ahat[0] = sprk_storage->a[2];
  sprk_storage->ahat[1] = sprk_storage->a[1];
  sprk_storage->ahat[2] = sprk_storage->a[0];
  return sprk_storage;
}

ARKodeSPRKStorage ARKodeSymplecticMcLachlan4()
{
  ARKodeSPRKStorage sprk_storage = ARKodeSPRKStorage_Alloc(4);
  if (!sprk_storage) { return NULL; }
  sprk_storage->q       = 4;
  sprk_storage->stages  = 4;
  sprk_storage->a[0]    = SUN_RCONST(0.515352837431122936);
  sprk_storage->a[1]    = -SUN_RCONST(0.085782019412973646);
  sprk_storage->a[2]    = SUN_RCONST(0.441583023616466524);
  sprk_storage->a[3]    = SUN_RCONST(0.128846158365384185);
  sprk_storage->ahat[0] = SUN_RCONST(0.134496199277431089);
  sprk_storage->ahat[1] = -SUN_RCONST(0.224819803079420806);
  sprk_storage->ahat[2] = SUN_RCONST(0.756320000515668291);
  sprk_storage->ahat[3] = SUN_RCONST(0.33400360328632142);
  return sprk_storage;
}

ARKodeSPRKStorage ARKodeSymplecticMcLachlan5()
{
  ARKodeSPRKStorage sprk_storage = ARKodeSPRKStorage_Alloc(6);
  if (!sprk_storage) { return NULL; }
  sprk_storage->q       = 5;
  sprk_storage->stages  = 6;
  sprk_storage->a[0]    = SUN_RCONST(0.339839625839110000);
  sprk_storage->a[1]    = -SUN_RCONST(0.088601336903027329);
  sprk_storage->a[2]    = SUN_RCONST(0.5858564768259621188);
  sprk_storage->a[3]    = -SUN_RCONST(0.603039356536491888);
  sprk_storage->a[4]    = SUN_RCONST(0.3235807965546976394);
  sprk_storage->a[5]    = SUN_RCONST(0.4423637942197494587);
  sprk_storage->ahat[0] = SUN_RCONST(0.1193900292875672758);
  sprk_storage->ahat[1] = SUN_RCONST(0.6989273703824752308);
  sprk_storage->ahat[2] = -SUN_RCONST(0.1713123582716007754);
  sprk_storage->ahat[3] = SUN_RCONST(0.4012695022513534480);
  sprk_storage->ahat[4] = SUN_RCONST(0.0107050818482359840);
  sprk_storage->ahat[5] = -SUN_RCONST(0.0589796254980311632);
  return sprk_storage;
}

/*
  The following methods are from:

  Yoshida, H.: Construction of higher order symplectic integrators.
  Phys Lett A. 150, 262–268 (1990).
  https://doi.org/10.1016/0375-9601(90)90092-3

 */

ARKodeSPRKStorage ARKodeSymplecticYoshida6()
{
  ARKodeSPRKStorage sprk_storage = ARKodeSPRKStorage_Alloc(8);
  if (!sprk_storage) { return NULL; }
  sprk_storage->q      = 6;
  sprk_storage->stages = 8;
  sprk_storage->a[0]   = SUN_RCONST(0.7845136104775572638194976338663498757768);
  sprk_storage->a[1]   = SUN_RCONST(0.2355732133593581336847931829785346016865);
  sprk_storage->a[2]   = -SUN_RCONST(1.177679984178871006946415680964315734639);
  sprk_storage->a[3]   = SUN_RCONST(1.315186320683911218884249728238862514352);
  sprk_storage->a[4]   = sprk_storage->a[2];
  sprk_storage->a[5]   = sprk_storage->a[1];
  sprk_storage->a[6]   = sprk_storage->a[0];
  sprk_storage->a[7]   = SUN_RCONST(0.0);
  sprk_storage->ahat[0] = sprk_storage->a[0] / SUN_RCONST(2.0);
  sprk_storage->ahat[1] = (sprk_storage->a[0] + sprk_storage->a[1]) /
                          SUN_RCONST(2.0);
  sprk_storage->ahat[2] = (sprk_storage->a[1] + sprk_storage->a[2]) /
                          SUN_RCONST(2.0);
  sprk_storage->ahat[3] = (sprk_storage->a[2] + sprk_storage->a[3]) /
                          SUN_RCONST(2.0);
  sprk_storage->ahat[4] = sprk_storage->ahat[3];
  sprk_storage->ahat[5] = sprk_storage->ahat[2];
  sprk_storage->ahat[6] = sprk_storage->ahat[1];
  sprk_storage->ahat[7] = sprk_storage->ahat[0];
  return sprk_storage;
}

/*
  The following methods are from:

  (Original) Suzuki, M., & Umeno, K. (1993). Higher-order decomposition theory
  of exponential operators and its applications to QMC and nonlinear dynamics.
  Computer simulation studies in condensed-matter physics VI, 74-86.
  https://doi.org/10.1007/978-3-642-78448-4_7

  McLachlan, R.I.: On the Numerical Integration of Ordinary Differential
  Equations by Symmetric Composition Methods. Siam J Sci Comput. 16, 151–168
  (1995). https://doi.org/10.1137/0916010

 */

ARKodeSPRKStorage ARKodeSymplecticSuzukiUmeno816()
{
  ARKodeSPRKStorage sprk_storage = ARKodeSPRKStorage_Alloc(16);
  if (!sprk_storage) { return NULL; }
  sprk_storage->q      = 8;
  sprk_storage->stages = 16;
  sprk_storage->a[0]   = SUN_RCONST(0.7416703643506129534482278017838063156035);
  sprk_storage->a[1]  = -SUN_RCONST(0.4091008258000315939973000958935634173099);
  sprk_storage->a[2]  = SUN_RCONST(0.1907547102962383799538762564503716627355);
  sprk_storage->a[3]  = -SUN_RCONST(0.5738624711160822666563877266355357421595);
  sprk_storage->a[4]  = SUN_RCONST(0.2990641813036559238444635406886029882258);
  sprk_storage->a[5]  = SUN_RCONST(0.3346249182452981837849579798821822886337);
  sprk_storage->a[6]  = SUN_RCONST(0.3152930923967665966320566638110024309941);
  sprk_storage->a[7]  = -SUN_RCONST(0.7968879393529163540197888401737330534463);
  sprk_storage->a[8]  = sprk_storage->a[6];
  sprk_storage->a[9]  = sprk_storage->a[5];
  sprk_storage->a[10] = sprk_storage->a[4];
  sprk_storage->a[11] = sprk_storage->a[3];
  sprk_storage->a[12] = sprk_storage->a[2];
  sprk_storage->a[13] = sprk_storage->a[1];
  sprk_storage->a[14] = sprk_storage->a[0];
  sprk_storage->a[15] = SUN_RCONST(0.0);
  sprk_storage->ahat[0] = sprk_storage->a[0] / SUN_RCONST(2.0);
  sprk_storage->ahat[1] = (sprk_storage->a[0] + sprk_storage->a[1]) /
                          SUN_RCONST(2.0);
  sprk_storage->ahat[2] = (sprk_storage->a[1] + sprk_storage->a[2]) /
                          SUN_RCONST(2.0);
  sprk_storage->ahat[3] = (sprk_storage->a[2] + sprk_storage->a[3]) /
                          SUN_RCONST(2.0);
  sprk_storage->ahat[4] = (sprk_storage->a[3] + sprk_storage->a[4]) /
                          SUN_RCONST(2.0);
  sprk_storage->ahat[5] = (sprk_storage->a[4] + sprk_storage->a[5]) /
                          SUN_RCONST(2.0);
  sprk_storage->ahat[6] = (sprk_storage->a[5] + sprk_storage->a[6]) /
                          SUN_RCONST(2.0);
  sprk_storage->ahat[7] = (sprk_storage->a[6] + sprk_storage->a[7]) /
                          SUN_RCONST(2.0);
  sprk_storage->ahat[8]  = sprk_storage->ahat[7];
  sprk_storage->ahat[9]  = sprk_storage->ahat[6];
  sprk_storage->ahat[10] = sprk_storage->ahat[5];
  sprk_storage->ahat[11] = sprk_storage->ahat[4];
  sprk_storage->ahat[12] = sprk_storage->ahat[3];
  sprk_storage->ahat[13] = sprk_storage->ahat[2];
  sprk_storage->ahat[14] = sprk_storage->ahat[1];
  sprk_storage->ahat[15] = sprk_storage->ahat[0];
  return sprk_storage;
}

/*
  The following methods are from:

  Sofroniou, M., Spaletta, G.: Derivation of symmetric composition constants for
  symmetric integrators. Optim Methods Softw. 20, 597–613 (2005).
  https://doi.org/10.1080/10556780500140664

 */

ARKodeSPRKStorage ARKodeSymplecticSofroniou10()
{
  ARKodeSPRKStorage sprk_storage = ARKodeSPRKStorage_Alloc(36);
  if (!sprk_storage) { return NULL; }
  sprk_storage->q      = 10;
  sprk_storage->stages = 36;

  sprk_storage->a[0]    = SUN_RCONST(0.078795722521686419263907679337684);
  sprk_storage->a[1]    = SUN_RCONST(0.31309610341510852776481247192647);
  sprk_storage->a[2]    = SUN_RCONST(0.027918383235078066109520273275299);
  sprk_storage->a[3]    = -SUN_RCONST(0.22959284159390709415121339679655);
  sprk_storage->a[4]    = SUN_RCONST(0.13096206107716486317465685927961);
  sprk_storage->a[5]    = -SUN_RCONST(0.26973340565451071434460973222411);
  sprk_storage->a[6]    = SUN_RCONST(0.074973343155891435666137105641410);
  sprk_storage->a[7]    = SUN_RCONST(0.11199342399981020488957508073640);
  sprk_storage->a[8]    = SUN_RCONST(0.36613344954622675119314812353150);
  sprk_storage->a[9]    = -SUN_RCONST(0.39910563013603589787862981058340);
  sprk_storage->a[10]   = SUN_RCONST(0.10308739852747107731580277001372);
  sprk_storage->a[11]   = SUN_RCONST(0.41143087395589023782070411897608);
  sprk_storage->a[12]   = -SUN_RCONST(0.0048663605831352617621956593099771);
  sprk_storage->a[13]   = -SUN_RCONST(0.39203335370863990644808193642610);
  sprk_storage->a[14]   = SUN_RCONST(0.051942502962449647037182904015976);
  sprk_storage->a[15]   = SUN_RCONST(0.050665090759924496335874344156866);
  sprk_storage->a[16]   = SUN_RCONST(0.049674370639729879054568800279461);
  sprk_storage->a[17]   = SUN_RCONST(0.049317735759594537917680008339338);
  sprk_storage->a[18]   = sprk_storage->a[16];
  sprk_storage->a[19]   = sprk_storage->a[15];
  sprk_storage->a[20]   = sprk_storage->a[14];
  sprk_storage->a[21]   = sprk_storage->a[13];
  sprk_storage->a[22]   = sprk_storage->a[12];
  sprk_storage->a[23]   = sprk_storage->a[11];
  sprk_storage->a[24]   = sprk_storage->a[10];
  sprk_storage->a[25]   = sprk_storage->a[9];
  sprk_storage->a[26]   = sprk_storage->a[8];
  sprk_storage->a[27]   = sprk_storage->a[7];
  sprk_storage->a[28]   = sprk_storage->a[6];
  sprk_storage->a[29]   = sprk_storage->a[5];
  sprk_storage->a[30]   = sprk_storage->a[4];
  sprk_storage->a[31]   = sprk_storage->a[3];
  sprk_storage->a[32]   = sprk_storage->a[2];
  sprk_storage->a[33]   = sprk_storage->a[1];
  sprk_storage->a[34]   = sprk_storage->a[0];
  sprk_storage->a[35]   = SUN_RCONST(0.0);
  sprk_storage->ahat[0] = sprk_storage->a[0] / SUN_RCONST(2.0);
  sprk_storage->ahat[1] = (sprk_storage->a[0] + sprk_storage->a[1]) /
                          SUN_RCONST(2.0);
  sprk_storage->ahat[2] = (sprk_storage->a[1] + sprk_storage->a[2]) /
                          SUN_RCONST(2.0);
  sprk_storage->ahat[3] = (sprk_storage->a[2] + sprk_storage->a[3]) /
                          SUN_RCONST(2.0);
  sprk_storage->ahat[4] = (sprk_storage->a[3] + sprk_storage->a[4]) /
                          SUN_RCONST(2.0);
  sprk_storage->ahat[5] = (sprk_storage->a[4] + sprk_storage->a[5]) /
                          SUN_RCONST(2.0);
  sprk_storage->ahat[6] = (sprk_storage->a[5] + sprk_storage->a[6]) /
                          SUN_RCONST(2.0);
  sprk_storage->ahat[7] = (sprk_storage->a[6] + sprk_storage->a[7]) /
                          SUN_RCONST(2.0);
  sprk_storage->ahat[8] = (sprk_storage->a[7] + sprk_storage->a[8]) /
                          SUN_RCONST(2.0);
  sprk_storage->ahat[9] = (sprk_storage->a[8] + sprk_storage->a[9]) /
                          SUN_RCONST(2.0);
  sprk_storage->ahat[10] = (sprk_storage->a[9] + sprk_storage->a[10]) /
                           SUN_RCONST(2.0);
  sprk_storage->ahat[11] = (sprk_storage->a[10] + sprk_storage->a[11]) /
                           SUN_RCONST(2.0);
  sprk_storage->ahat[12] = (sprk_storage->a[11] + sprk_storage->a[12]) /
                           SUN_RCONST(2.0);
  sprk_storage->ahat[13] = (sprk_storage->a[12] + sprk_storage->a[13]) /
                           SUN_RCONST(2.0);
  sprk_storage->ahat[14] = (sprk_storage->a[13] + sprk_storage->a[14]) /
                           SUN_RCONST(2.0);
  sprk_storage->ahat[15] = (sprk_storage->a[14] + sprk_storage->a[15]) /
                           SUN_RCONST(2.0);
  sprk_storage->ahat[16] = (sprk_storage->a[15] + sprk_storage->a[16]) /
                           SUN_RCONST(2.0);
  sprk_storage->ahat[17] = (sprk_storage->a[16] + sprk_storage->a[17]) /
                           SUN_RCONST(2.0);
  sprk_storage->ahat[18] = sprk_storage->ahat[17];
  sprk_storage->ahat[19] = sprk_storage->ahat[16];
  sprk_storage->ahat[20] = sprk_storage->ahat[15];
  sprk_storage->ahat[21] = sprk_storage->ahat[14];
  sprk_storage->ahat[22] = sprk_storage->ahat[13];
  sprk_storage->ahat[23] = sprk_storage->ahat[12];
  sprk_storage->ahat[24] = sprk_storage->ahat[11];
  sprk_storage->ahat[25] = sprk_storage->ahat[10];
  sprk_storage->ahat[26] = sprk_storage->ahat[9];
  sprk_storage->ahat[27] = sprk_storage->ahat[8];
  sprk_storage->ahat[28] = sprk_storage->ahat[7];
  sprk_storage->ahat[29] = sprk_storage->ahat[6];
  sprk_storage->ahat[30] = sprk_storage->ahat[5];
  sprk_storage->ahat[31] = sprk_storage->ahat[4];
  sprk_storage->ahat[32] = sprk_storage->ahat[3];
  sprk_storage->ahat[33] = sprk_storage->ahat[2];
  sprk_storage->ahat[34] = sprk_storage->ahat[1];
  sprk_storage->ahat[35] = sprk_storage->ahat[0];

  return sprk_storage;
}

ARKodeSPRKStorage ARKodeSPRKStorage_Alloc(int stages)
{
  ARKodeSPRKStorage sprk_storage = NULL;

  sprk_storage = (ARKodeSPRKStorage)malloc(sizeof(struct ARKodeSPRKStorage_s));
  if (!sprk_storage) { return NULL; }

  memset(sprk_storage, 0, sizeof(struct ARKodeSPRKStorage_s));

  sprk_storage->ahat = (sunrealtype*)malloc(stages * sizeof(sunrealtype));
  if (!(sprk_storage->ahat))
  {
    ARKodeSPRKStorage_Free(sprk_storage);
    return NULL;
  }

  sprk_storage->a = (sunrealtype*)malloc(stages * sizeof(sunrealtype));
  if (!(sprk_storage->a))
  {
    ARKodeSPRKStorage_Free(sprk_storage);
    return NULL;
  }

  sprk_storage->stages = stages;
  sprk_storage->ahat   = (sunrealtype*)malloc(stages * sizeof(sunrealtype));
  sprk_storage->a      = (sunrealtype*)malloc(stages * sizeof(sunrealtype));

  return sprk_storage;
}

ARKodeSPRKStorage ARKodeSPRKStorage_Load(ARKODE_SPRKMethodID id)
{
  switch (id)
  {
  case ARKODE_SPRK_EULER_1_1: return ARKodeSymplecticEuler();
  case ARKODE_SPRK_LEAPFROG_2_2: return ARKodeSymplecticLeapfrog2();
  case ARKODE_SPRK_PSEUDO_LEAPFROG_2_2:
    return ARKodeSymplecticPseudoLeapfrog2();
  case ARKODE_SPRK_RUTH_3_3: return ARKodeSymplecticRuth3();
  case ARKODE_SPRK_MCLACHLAN_2_2: return ARKodeSymplecticMcLachlan2();
  case ARKODE_SPRK_MCLACHLAN_3_3: return ARKodeSymplecticMcLachlan3();
  case ARKODE_SPRK_MCLACHLAN_4_4: return ARKodeSymplecticMcLachlan4();
  case ARKODE_SPRK_CANDY_ROZMUS_4_4: return ARKodeSymplecticCandyRozmus4();
  case ARKODE_SPRK_MCLACHLAN_5_6: return ARKodeSymplecticMcLachlan5();
  case ARKODE_SPRK_YOSHIDA_6_8: return ARKodeSymplecticYoshida6();
  case ARKODE_SPRK_SUZUKI_UMENO_8_16: return ARKodeSymplecticSuzukiUmeno816();
  case ARKODE_SPRK_SOFRONIOU_10_36: return ARKodeSymplecticSofroniou10();
  default: return NULL;
  }
}

ARKodeSPRKStorage ARKodeSPRKStorage_LoadByName(const char* method)
{
  if (!strcmp(method, "ARKODE_SPRK_EULER_1_1"))
  {
    return ARKodeSymplecticEuler();
  }
  if (!strcmp(method, "ARKODE_SPRK_LEAPFROG_2_2"))
  {
    return ARKodeSymplecticLeapfrog2();
  }
  if (!strcmp(method, "ARKODE_SPRK_PSEUDO_LEAPFROG_2_2"))
  {
    return ARKodeSymplecticPseudoLeapfrog2();
  }
  if (!strcmp(method, "ARKODE_SPRK_RUTH_3_3"))
  {
    return ARKodeSymplecticRuth3();
  }
  if (!strcmp(method, "ARKODE_SPRK_MCLACHLAN_2_2"))
  {
    return ARKodeSymplecticMcLachlan2();
  }
  if (!strcmp(method, "ARKODE_SPRK_MCLACHLAN_3_3"))
  {
    return ARKodeSymplecticMcLachlan3();
  }
  if (!strcmp(method, "ARKODE_SPRK_MCLACHLAN_4_4"))
  {
    return ARKodeSymplecticMcLachlan4();
  }
  if (!strcmp(method, "ARKODE_SPRK_CANDY_ROZMUS_4_4"))
  {
    return ARKodeSymplecticCandyRozmus4();
  }
  if (!strcmp(method, "ARKODE_SPRK_MCLACHLAN_5_6"))
  {
    return ARKodeSymplecticMcLachlan5();
  }
  if (!strcmp(method, "ARKODE_SPRK_YOSHIDA_6_8"))
  {
    return ARKodeSymplecticYoshida6();
  }
  if (!strcmp(method, "ARKODE_SPRK_SUZUKI_UMENO_8_16"))
  {
    return ARKodeSymplecticSuzukiUmeno816();
  }
  if (!strcmp(method, "ARKODE_SPRK_SOFRONIOU_10_36"))
  {
    return ARKodeSymplecticSofroniou10();
  }
  return NULL;
}

ARKodeSPRKStorage ARKodeSPRKStorage_Copy(ARKodeSPRKStorage that_sprk_storage)
{
  int i                          = 0;
  ARKodeSPRKStorage sprk_storage = NULL;

  sprk_storage = ARKodeSPRKStorage_Alloc(that_sprk_storage->stages);

  sprk_storage->q = that_sprk_storage->q;

  for (i = 0; i < sprk_storage->stages; ++i)
  {
    sprk_storage->ahat[i] = that_sprk_storage->ahat[i];
    sprk_storage->a[i]    = that_sprk_storage->a[i];
  }

  return sprk_storage;
}

void ARKodeSPRKStorage_Space(ARKodeSPRKStorage sprk_storage, sunindextype* liw,
                             sunindextype* lrw)
{
  *liw = 2;
  *lrw = sprk_storage->stages * 2;
}

void ARKodeSPRKStorage_Free(ARKodeSPRKStorage sprk_storage)
{
  if (sprk_storage)
  {
    if (sprk_storage->ahat) { free(sprk_storage->ahat); }
    if (sprk_storage->a) { free(sprk_storage->a); }
    free(sprk_storage);
  }
}

int ARKodeSPRKStorage_ToButcher(ARKodeSPRKStorage sprk_storage,
                                ARKodeButcherTable* erk_ptr,
                                ARKodeButcherTable* dirk_ptr)
{
  int i                = 0;
  int j                = 0;
  ARKodeButcherTable a = NULL;
  ARKodeButcherTable b = NULL;

  a = ARKodeButcherTable_Alloc(sprk_storage->stages, SUNFALSE);
  if (!a) { return ARK_MEM_FAIL; }
  b = ARKodeButcherTable_Alloc(sprk_storage->stages, SUNFALSE);
  if (!b) { return ARK_MEM_FAIL; }

  /* DIRK table */
  for (i = 0; i < sprk_storage->stages; ++i)
  {
    b->b[i] = sprk_storage->ahat[i];
    for (j = 0; j <= i; ++j) { b->A[i][j] = sprk_storage->ahat[j]; }
    /* Time weights: C_j = sum_{i=0}^{j} b_i */

    /* Time weights: C_j = sum_{i=0}^{j-1} b_i */
    for (j = 0; j < sprk_storage->stages; ++j)
    {
      for (i = 0; i <= j; ++i) { b->c[j] += sprk_storage->ahat[i]; }
    }

    /* Explicit table */
    for (i = 0; i < sprk_storage->stages; ++i)
    {
      a->b[i] = sprk_storage->a[i];
      for (j = 0; j < i; ++j) { a->A[i][j] = sprk_storage->a[j]; }
    }

    /* Time weights: c_j = sum_{i=0}^{j-1} a_i */
    for (j = 0; j < sprk_storage->stages; ++j)
    {
      for (i = 0; i < j; ++i) { a->c[j] += sprk_storage->a[i]; }
    }

    /* Set method order */
    a->q = sprk_storage->q;
    b->q = sprk_storage->q;

    /* No embedding, so set embedding order to 0 */
    a->p = 0;
    b->p = 0;

  }

  *erk_ptr  = a;
  *dirk_ptr = b;

  return ARK_SUCCESS;
}
