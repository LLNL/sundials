/* -----------------------------------------------------------------------------
 * Programmer(s): David J. Gardner @ LLNL
 * -----------------------------------------------------------------------------
 * SUNDIALS Copyright Start
 * Copyright (c) 2002-2022, Lawrence Livermore National Security
 * and Southern Methodist University.
 * All rights reserved.
 *
 * See the top-level LICENSE and NOTICE files for details.
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * SUNDIALS Copyright End
 * -----------------------------------------------------------------------------
 * Serial diffusion functions
 * ---------------------------------------------------------------------------*/

#include "diffusion_2D.hpp"

// Diffusion function
int laplacian(realtype t, N_Vector u, N_Vector f, void *user_data)
{
  int          flag;
  sunindextype i, j;

  // Access problem data
  UserData *udata = (UserData *) user_data;

  SUNDIALS_CXX_MARK_FUNCTION(udata->prof);

  // Start exchange
  flag = udata->start_exchange(u);
  if (check_flag(&flag, "SendData", 1)) return -1;

  // Shortcuts to local number of nodes
  sunindextype nx_loc = udata->nx_loc;
  sunindextype ny_loc = udata->ny_loc;

  // Determine iteration range excluding the overall domain boundary
  sunindextype istart = (udata->HaveNbrW) ? 0      : 1;
  sunindextype iend   = (udata->HaveNbrE) ? nx_loc : nx_loc - 1;
  sunindextype jstart = (udata->HaveNbrS) ? 0      : 1;
  sunindextype jend   = (udata->HaveNbrN) ? ny_loc : ny_loc - 1;

  // Constants for computing diffusion term
  realtype cx = udata->kx / (udata->dx * udata->dx);
  realtype cy = udata->ky / (udata->dy * udata->dy);
  realtype cc = -TWO * (cx + cy);

  // Access data arrays
  realtype *uarray = N_VGetArrayPointer(u);
  if (check_flag((void *) uarray, "N_VGetArrayPointer", 0)) return -1;

  realtype *farray = N_VGetArrayPointer(f);
  if (check_flag((void *) farray, "N_VGetArrayPointer", 0)) return -1;

  // Initialize rhs vector to zero (handles boundary conditions)
  N_VConst(ZERO, f);

  // Iterate over subdomain and compute rhs forcing term
  if (udata->forcing)
  {
    realtype x, y;
    realtype sin_sqr_x, sin_sqr_y;
    realtype cos_sqr_x, cos_sqr_y;

    realtype bx = (udata->kx) * TWO * PI * PI;
    realtype by = (udata->ky) * TWO * PI * PI;

    realtype sin_t_cos_t = sin(PI * t) * cos(PI * t);
    realtype cos_sqr_t   = cos(PI * t) * cos(PI * t);

    for (j = jstart; j < jend; j++)
    {
      for (i = istart; i < iend; i++)
      {
        x = (udata->is + i) * udata->dx;
        y = (udata->js + j) * udata->dy;

        sin_sqr_x = sin(PI * x) * sin(PI * x);
        sin_sqr_y = sin(PI * y) * sin(PI * y);

        cos_sqr_x = cos(PI * x) * cos(PI * x);
        cos_sqr_y = cos(PI * y) * cos(PI * y);

        farray[IDX(i,j,nx_loc)] =
          -TWO * PI * sin_sqr_x * sin_sqr_y * sin_t_cos_t
          -bx * (cos_sqr_x - sin_sqr_x) * sin_sqr_y * cos_sqr_t
          -by * (cos_sqr_y - sin_sqr_y) * sin_sqr_x * cos_sqr_t;
      }
    }
  }

  // Iterate over subdomain interior and add rhs diffusion term
  for (j = 1; j < ny_loc - 1; j++)
  {
    for (i = 1; i < nx_loc - 1; i++)
    {
      farray[IDX(i,j,nx_loc)] +=
        cc * uarray[IDX(i,j,nx_loc)]
        + cx * (uarray[IDX(i-1,j,nx_loc)] + uarray[IDX(i+1,j,nx_loc)])
        + cy * (uarray[IDX(i,j-1,nx_loc)] + uarray[IDX(i,j+1,nx_loc)]);
    }
  }

  // Wait for exchange receives
  flag = udata->end_exchange();
  if (check_flag(&flag, "UserData::end_excahnge", 1)) return -1;

  // Iterate over subdomain boundaries and add rhs diffusion term
  realtype *Warray = udata->Wrecv;
  realtype *Earray = udata->Erecv;
  realtype *Sarray = udata->Srecv;
  realtype *Narray = udata->Nrecv;

  // West face (updates south-west and north-west corners if necessary)
  if (udata->HaveNbrW)
  {
    i = 0;
    if (udata->HaveNbrS)  // South-West corner
    {
      j = 0;
      farray[IDX(i,j,nx_loc)] +=
        cc * uarray[IDX(i,j,nx_loc)]
        + cx * (Warray[j] + uarray[IDX(i+1,j,nx_loc)])
        + cy * (Sarray[i] + uarray[IDX(i,j+1,nx_loc)]);
    }

    for (j = 1; j < ny_loc - 1; j++)
    {
      farray[IDX(i,j,nx_loc)] +=
        cc * uarray[IDX(i,j,nx_loc)]
        + cx * (Warray[j] + uarray[IDX(i+1,j,nx_loc)])
        + cy * (uarray[IDX(i,j-1,nx_loc)] + uarray[IDX(i,j+1,nx_loc)]);
    }

    if (udata->HaveNbrN)  // North-West corner
    {
      j = ny_loc - 1;
      farray[IDX(i,j,nx_loc)] +=
        cc * uarray[IDX(i,j,nx_loc)]
        + cx * (Warray[j] + uarray[IDX(i+1,j,nx_loc)])
        + cy * (uarray[IDX(i,j-1,nx_loc)] + Narray[i]);
    }
  }

  // East face (updates south-east and north-east corners if necessary)
  if (udata->HaveNbrE)
  {
    i = nx_loc - 1;
    if (udata->HaveNbrS)  // South-East corner
    {
      j = 0;
      farray[IDX(i,j,nx_loc)] +=
        cc * uarray[IDX(i,j,nx_loc)]
        + cx * (uarray[IDX(i-1,j,nx_loc)] + Earray[j])
        + cy * (Sarray[i] + uarray[IDX(i,j+1,nx_loc)]);
    }

    for (j = 1; j < ny_loc - 1; j++)
    {
      farray[IDX(i,j,nx_loc)] +=
        cc * uarray[IDX(i,j,nx_loc)]
        + cx * (uarray[IDX(i-1,j,nx_loc)] + Earray[j])
        + cy * (uarray[IDX(i,j-1,nx_loc)] + uarray[IDX(i,j+1,nx_loc)]);
    }

    if (udata->HaveNbrN)  // North-East corner
    {
      j = ny_loc - 1;
      farray[IDX(i,j,nx_loc)] +=
        cc * uarray[IDX(i,j,nx_loc)]
        + cx * (uarray[IDX(i-1,j,nx_loc)] + Earray[j])
        + cy * (uarray[IDX(i,j-1,nx_loc)] + Narray[i]);
    }
  }

  // South face (excludes corners)
  if (udata->HaveNbrS)
  {
    j = 0;
    for (i = 1; i < nx_loc - 1; i++)
    {
      farray[IDX(i,j,nx_loc)] +=
        cc * uarray[IDX(i,j,nx_loc)]
        + cx * (uarray[IDX(i-1,j,nx_loc)] + uarray[IDX(i+1,j,nx_loc)])
        + cy * (Sarray[i] + uarray[IDX(i,j+1,nx_loc)]);
    }
  }

  // North face (excludes corners)
  if (udata->HaveNbrN)
  {
    j = udata->ny_loc - 1;
    for (i = 1; i < nx_loc - 1; i++)
    {
      farray[IDX(i,j,nx_loc)] +=
        cc * uarray[IDX(i,j,nx_loc)]
        + cx * (uarray[IDX(i-1,j,nx_loc)] + uarray[IDX(i+1,j,nx_loc)])
        + cy * (uarray[IDX(i,j-1,nx_loc)] + Narray[i]);
    }
  }

  // Return success
  return 0;
}

#if defined(USE_SUPERLU_DIST)
// Compute the global index of a node
//   i -- local x index
//   j -- local y index
//   x -- x processor coordinate
//   y -- y processor coordinate
//   qx -- base number of nodes per x process
//   qy -- base number of nodes per y process
//   rx -- number of extra nodes for x processes
//   ry -- number of extra nodes for y processes
//   nx -- global number of nodes in the x-direction
sunindextype global_index(sunindextype i, sunindextype j, int x, int y,
                          sunindextype qx, sunindextype qy, sunindextype rx,
                          sunindextype ry, sunindextype nx)
{
  // offset from previous process rows
  sunindextype offset_p = nx * ((qy + 1) * min(y, ry) + qy * max(y - ry, 0));

  // offset within current process row
  sunindextype offset_c;
  if (y < ry)
  {
    offset_c = (qy + 1) * ((qx + 1) * min(x, rx) + qx * max(x - rx, 0));
  }
  else
  {
    offset_c = qy * ((qx + 1) * min(x, rx) + qx * max(x - rx, 0));
  }

  // base node index for this process
  sunindextype base = offset_p + offset_c;

  // number of local x nodes
  sunindextype nx_loc;
  if (x < rx) nx_loc = qx + 1;
  else nx_loc = qx;

  // global index for the requested local node
  sunindextype idx = base + j * nx_loc + i;

  return idx;
}

// Compute the global column indices for a row
//   i -- local x index
//   j -- local y index
//   x -- x processor coordinate
//   y -- y processor coordinate
//   qx -- base number of nodes per x process
//   qy -- base number of nodes per y process
//   rx -- number of extra nodes for x processes
//   ry -- number of extra nodes for y processes
//   nx -- global number of nodes in the x-direction
int matrix_columns(sunindextype i, sunindextype j, int x, int y, int npx,
                   int npy, sunindextype qx, sunindextype qy, sunindextype rx,
                   sunindextype ry, sunindextype nx, sunrealtype* vals,
                   sunindextype* col_idx, sunindextype* row_nnz)
{
  // number of local x nodes
  if (x < rx) nx_loc = qx + 1;
  else nx_loc = qx;

  // number of local y nodes
  if (y < ry) ny_loc = qy + 1;
  else ny_loc = qy;

  // -----------------
  // global boundaries
  // -----------------

  // south boundary
  if (y == 0 && j == 0)
  {
    c_idx = global_index(i, j, x, y, qx, qy, rx, ry, nx);
    return [0], [c_idx], 1;
  }

  // west boundary
  if (x == 0 && i == 0)
  {
    c_idx = global_index(i, j, x, y, qx, qy, rx, ry, nx);
    return [0], [c_idx], 1;
  }

  // east boundary
  if (x == npx - 1 && i == nx_loc - 1)
  {
    c_idx = global_index(i, j, x, y, qx, qy, rx, ry, nx);
    return [0], [c_idx], 1;
  }

  // north boundary
  if (y == npx - 1 && j == ny_loc - 1)
  {
    c_idx = global_index(i, j, x, y, qx, qy, rx, ry, nx);
    return [0], [c_idx], 1;
  }

  // --------
  // interior
  // --------

  // ordering of columns for output
  //   0. South Neighbor
  //   1. West Neighbor
  //   Interior:
  //      2. South
  //      3. West
  //      4. Center
  //      5. East
  //      6. North
  //   7. East Neighbor
  //   8. North Neighbor

  // List of bools for where values come from
  values  = 5 * [None];
  columns = 5 * [None];
  idx     = 0;

  // south neighbor
  if (j == 0)
  {
    if (y - 1 < ry) ny_loc = qy + 1;
    else ny_loc = qy;
    values[idx]  = cy;
    columns[idx] = global_index(i, ny_loc - 1, x, y - 1, qx, qy, rx, ry, nx);
    idx += 1;
  }

  // west neighbor
  if (i == 0)
  {
    if (x - 1 < rx) nx_loc = qx + 1;
    else:
      nx_loc = qx;
    values[idx]  = cx;
    columns[idx] = global_index(nx_loc - 1, j, x - 1, y, qx, qy, rx, ry, nx);
    idx += 1;
  }

  // south interior
  if (j > 0)
  {
    values[idx]  = cy;
    columns[idx] = global_index(i, j - 1, x, y, qx, qy, rx, ry, nx);
    idx += 1;
  }

  // west interior
  if (i > 0)
  {
    values[idx]  = cx;
    columns[idx] = global_index(i - 1, j, x, y, qx, qy, rx, ry, nx);
    idx += 1;
  }

  // center
  values[idx]  = cc;
  columns[idx] = global_index(i, j, x, y, qx, qy, rx, ry, nx);
  idx += 1;

  // east interior
  if (i < nx_loc - 1)
  {
    values[idx]  = cx;
    columns[idx] = global_index(i + 1, j, x, y, qx, qy, rx, ry, nx);
    idx += 1;
  }

  // north interior
  if (j < ny_loc - 1)
  {
    values[idx]  = cy;
    columns[idx] = global_index(i, j + 1, x, y, qx, qy, rx, ry, nx);
    idx += 1;
  }

  // east neighbor
  if (i == nx_loc - 1)
  {
    values[idx]  = cx;
    columns[idx] = global_index(0, j, x + 1, y, qx, qy, rx, ry, nx);
    idx += 1;
  }

  // north neighbor
  if (j == ny_loc - 1)
  {
    values[idx]  = cy;
    columns[idx] = global_index(i, 0, x, y + 1, qx, qy, rx, ry, nx);
    idx += 1;
  }

  return values, columns, 5;
}

int laplacian_matrix(N_Vector u, SUNMatrix L, void* user_data)
{
  // Access problem data
  UserData* udata = (UserData*)user_data;

  // Set shortcuts
  SuperMatrix* Lsuper   = SUNMatrix_SLUNRloc_SuperMatrix(L);
  NRformat_loc* Lstore  = (NRformat_loc*)Lsuper->Store;
  sunindextype* rowptrs = Lstore->rowptr;
  sunindextype* colinds = Lstore->colind;
  realtype* Ldata       = (realtype*)Lstore->nzval;

  sunindextype rowptr = 0;
  for (sunindextype j = 0; j < udata->ny_loc; j++)
  {
    for (sunindextype i = 0; i < udata->nx_loc; i++)
    {
      matrix_columns(i, j, x, y, npx, npy, qx, qy, rx, ry, nx);
      rowptr += rptr;
    }
  }

  // Return success
  return 0;
}
#endif
