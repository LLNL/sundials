/*******************************************************************
 *                                                                 *
 * File          : sundialsmath.c                                  *
 * Programmers   : Scott D. Cohen and Alan C. Hindmarsh @ LLNL     *
 * Version of    : 26 June 2002                                    *
 *-----------------------------------------------------------------*
 * Copyright (c) 2002, The Regents of the University of California *
 * Produced at the Lawrence Livermore National Laboratory          *
 * All rights reserved                                             *
 * For details, see sundials/shared/LICENSE                        *
 *-----------------------------------------------------------------*
 * This is the implementation file for a C math library.           *
 *                                                                 *
 *******************************************************************/


#include <stdio.h>
#include <math.h>
#include "sundialsmath.h"
#include "sundialstypes.h"


#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)
#define TWO  RCONST(2.0)


realtype UnitRoundoff(void)
{
  realtype u;
  volatile realtype one_plus_u;
  
  u = ONE;
  one_plus_u = ONE + u;
  while (one_plus_u != ONE) {
    u /=  TWO;
    one_plus_u = ONE + u;
  }
  u *=  TWO;
  
  return(u);
}


realtype RPowerI(realtype base, int exponent)
{
  int i, expt;
  realtype prod;

  prod = ONE;
  expt = ABS(exponent);
  for(i=1; i <= expt; i++) prod *= base;
  if (exponent < 0) prod = ONE/prod;
  return(prod);
}


realtype RPowerR(realtype base, realtype exponent)
{
 
  if (base <= ZERO) return(ZERO);

  return((realtype)pow((double)base,(double)exponent));
}


realtype RSqrt(realtype x)
{
  if (x <= ZERO) return(ZERO);

  return((realtype) sqrt((double) x));
}
