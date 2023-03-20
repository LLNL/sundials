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
 * CUDA/HIP diffusion functions
 * ---------------------------------------------------------------------------*/

#include "diffusion_2D.hpp"

#if defined(USE_HIP)
#include <hip/hip_runtime.h>
#define BLOCK_SIZE   256
#define BLOCK_SIZE_X 16
#define BLOCK_SIZE_Y 16
#elif defined(USE_CUDA)
#include <cuda_runtime.h>
#define BLOCK_SIZE   256
#define BLOCK_SIZE_X 16
#define BLOCK_SIZE_Y 16
#else
#error Define USE_CUDA or USE_HIP
#endif


// Forcing device function
__device__ void add_forcing(const realtype t,
                            const realtype x, const realtype y,
                            const realtype kx, const realtype ky,
                            const sunindextype c, realtype *udot)
{
  realtype sin_sqr_x = sin(PI * x) * sin(PI * x);
  realtype sin_sqr_y = sin(PI * y) * sin(PI * y);

  realtype cos_sqr_x = cos(PI * x) * cos(PI * x);
  realtype cos_sqr_y = cos(PI * y) * cos(PI * y);

  realtype sin_t_cos_t = sin(PI * t) * cos(PI * t);
  realtype cos_sqr_t   = cos(PI * t) * cos(PI * t);

  realtype bx = kx * TWO * PI * PI;
  realtype by = ky * TWO * PI * PI;

  udot[c] += -TWO * PI * sin_sqr_x * sin_sqr_y * sin_t_cos_t
    -bx * (cos_sqr_x - sin_sqr_x) * sin_sqr_y * cos_sqr_t
    -by * (cos_sqr_y - sin_sqr_y) * sin_sqr_x * cos_sqr_t;
}


// Interior diffusion kernel
__global__ void diffusion_interior_kernel(const realtype t,
                                          const realtype *u,
                                          realtype *udot,
                                          const sunindextype is,
                                          const sunindextype js,
                                          const sunindextype nx_loc,
                                          const sunindextype ny_loc,
                                          const realtype dx,
                                          const realtype dy,
                                          const realtype kx,
                                          const realtype ky,
                                          const bool forcing)
{
  // Thread location in the local grid
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;

  // Only update the interior points
  bool interior = (i > 0 && i < nx_loc - 1 &&
                   j > 0 && j < ny_loc - 1);
  if (interior)
  {
    // 1D array index for center, west, east, south, and north nodes
    int c = i + j * nx_loc;
    int w = c - 1;
    int e = c + 1;
    int s = c - nx_loc;
    int n = c + nx_loc;

    // Set diffusion term
    realtype cx = kx / (dx * dx);
    realtype cy = ky / (dy * dy);
    realtype cc = -TWO * (cx + cy);

    udot[c] = cc * u[c] + cx * (u[w] + u[e]) + cy * (u[s] + u[n]);

    if (forcing)
    {
      realtype x = (is + i) * dx;
      realtype y = (js + j) * dy;

      add_forcing(t, x, y, kx, ky, c, udot);
    }
  }
}


// Interior boundary kernel
__global__ void diffusion_boundary_kernel(const realtype t,
                                          const realtype *u,
                                          realtype *udot,
                                          const sunindextype is,
                                          const sunindextype js,
                                          const sunindextype nx_loc,
                                          const sunindextype ny_loc,
                                          const realtype dx,
                                          const realtype dy,
                                          const realtype kx,
                                          const realtype ky,
                                          const bool forcing,
                                          const realtype *wbuf,
                                          const realtype *ebuf,
                                          const realtype *sbuf,
                                          const realtype *nbuf)
{
  // Thread ID
  int i = blockIdx.x * blockDim.x + threadIdx.x;

  // Set diffusion term
  realtype cx = kx / (dx * dx);
  realtype cy = ky / (dy * dy);
  realtype cc = -TWO * (cx + cy);

  // West and East faces excluding corners
  if (i > 0 && i < ny_loc - 1)
  {
    // West face: 1D array index for center, west, east, south, and north nodes
    int c = i * nx_loc;
    int w = i;
    int e = c + 1;
    int s = c - nx_loc;
    int n = c + nx_loc;

    if (wbuf)
    {
      // West processor boundary
      udot[c] = cc * u[c] + cx * (wbuf[w] + u[e]) + cy * (u[s] + u[n]);

      if (forcing)
      {
        realtype x = is * dx;
        realtype y = (js + i) * dy;

        add_forcing(t, x, y, kx, ky, c, udot);
      }
    }
    else
    {
      // West physical boundary
      udot[c] = ZERO;
    }

    // East face: 1D array index for center, west, east, south, and north nodes
    c = (i + 1) * nx_loc - 1;
    w = c - 1;
    e = i;
    s = c - nx_loc;
    n = c + nx_loc;

    if (ebuf)
    {
      // East processor boundary
      udot[c] = cc * u[c] + cx * (u[w] + ebuf[e]) + cy * (u[s] + u[n]);

      if (forcing)
      {
        realtype x = (is + nx_loc - 1) * dx;
        realtype y = (js + i) * dy;

        add_forcing(t, x, y, kx, ky, c, udot);
      }
    }
    else
    {
      // East physical boundary
      udot[c] = ZERO;
    }
  }

  __syncthreads();

  // South and North faces excluding corners
  if (i > 0 && i < nx_loc - 1)
  {
    // South face: 1D array index for center, west, east, south, and north nodes
    int c = i;
    int w = c - 1;
    int e = c + 1;
    int s = i;
    int n = c + nx_loc;

    if (sbuf)
    {
      // South processor boundary
      udot[c] = cc * u[c] + cx * (u[w] + u[e]) + cy * (sbuf[s] + u[n]);

      if (forcing)
      {
        realtype x = (is + i) * dx;
        realtype y = js * dy;

        add_forcing(t, x, y, kx, ky, c, udot);
      }
    }
    else
    {
      // South physical boundary
      udot[c] = ZERO;
    }

    // North face: 1D array index for center, west, east, south, and north nodes
    c = i + (ny_loc - 1) * nx_loc;
    w = c - 1;
    e = c + 1;
    s = c - nx_loc;
    n = i;

    if (nbuf)
    {
      // North processor boundary
      udot[c] = cc * u[c] + cx * (u[w] + u[e]) + cy * (u[s] + nbuf[n]);

      if (forcing)
      {
        realtype x = (is + i) * dx;
        realtype y = (js + ny_loc - 1) * dy;

        add_forcing(t, x, y, kx, ky, c, udot);
      }
    }
    else
    {
      // North physical boundary
      udot[c] = ZERO;
    }
  }

  __syncthreads();

  // Corners
  if (i == 0)
  {
    // South-West
    int c = 0;
    int w = 0;
    int e = c + 1;
    int s = 0;
    int n = c + nx_loc;

    if (sbuf && wbuf)
    {
      // South-West processor boundary
      udot[c] = cc * u[c] + cx * (wbuf[w] + u[e]) + cy * (sbuf[s] + u[n]);

      if (forcing)
      {
        realtype x = is * dx;
        realtype y = js * dy;

        add_forcing(t, x, y, kx, ky, c, udot);
      }
    }
    else
    {
      // South-West physical boundary
      udot[c] = ZERO;
    }

    // South-East
    c = nx_loc - 1;
    w = c - 1;
    e = 0;
    s = nx_loc - 1;
    n = c + nx_loc;

    if (sbuf && ebuf)
    {
      // South-East processor boundary
      udot[c] = cc * u[c] + cx * (u[w] + ebuf[e]) + cy * (sbuf[s] + u[n]);

      if (forcing)
      {
        realtype x = (is + nx_loc - 1) * dx;
        realtype y = js * dy;

        add_forcing(t, x, y, kx, ky, c, udot);
      }
    }
    else
    {
      // South-East physical boundary
      udot[c] = ZERO;
    }

    // North-West
    c = (ny_loc - 1) * nx_loc;
    w = ny_loc - 1;
    e = c + 1;
    s = c - nx_loc;
    n = 0;

    if (nbuf && wbuf)
    {
      // North-West processor boundary
      udot[c] = cc * u[c] + cx * (wbuf[w] + u[e]) + cy * (u[s] + nbuf[n]);

      if (forcing)
      {
        realtype x = is * dx;
        realtype y = (js + ny_loc - 1) * dy;

        add_forcing(t, x, y, kx, ky, c, udot);
      }
    }
    else
    {
      // North-West physical boundary
      udot[c] = ZERO;
    }

    // North-East
    c = nx_loc * ny_loc - 1;
    w = c - 1;
    e = ny_loc - 1;
    s = c - nx_loc;
    n = nx_loc - 1;

    if (nbuf && ebuf)
    {
      // North-East processor boundary
      udot[c] = cc * u[c] + cx * (u[w] + ebuf[e]) + cy * (u[s] + nbuf[n]);

      if (forcing)
      {
        realtype x = (is + nx_loc - 1) * dx;
        realtype y = (js + ny_loc - 1) * dy;

        add_forcing(t, x, y, kx, ky, c, udot);
      }
    }
    else
    {
      // North-East physical boundary
      udot[c] = ZERO;
    }
  }
}


// Diffusion function
int laplacian(realtype t, N_Vector u, N_Vector f, UserData* udata)
{
  SUNDIALS_CXX_MARK_FUNCTION(udata->prof);

  int flag;

  // Start exchange
  flag = udata->start_exchange(u);
  if (check_flag(&flag, "UserData::start_exchange", 1)) return -1;

  // Extract needed constants from user data
  const sunindextype is      = udata->is;
  const sunindextype js      = udata->js;
  const sunindextype nx_loc  = udata->nx_loc;
  const sunindextype ny_loc  = udata->ny_loc;
  const realtype     dx      = udata->dx;
  const realtype     dy      = udata->dy;
  const realtype     kx      = udata->kx;
  const realtype     ky      = udata->ky;
  const bool         forcing = udata->forcing;

  // Access data arrays
  const realtype *uarray = N_VGetDeviceArrayPointer(N_VGetLocalVector_MPIPlusX(u));
  if (check_flag((void *) uarray, "N_VGetDeviceArrayPointer", 0)) return -1;

  realtype *farray = N_VGetDeviceArrayPointer(N_VGetLocalVector_MPIPlusX(f));
  if (check_flag((void *) farray, "N_VGetDeviceArrayPointer", 0)) return -1;

  // Update subdomain interior
  dim3 iblock(BLOCK_SIZE_X, BLOCK_SIZE_Y);
  dim3 igrid(ICEIL(nx_loc, BLOCK_SIZE_X), ICEIL(ny_loc, BLOCK_SIZE_Y));

  diffusion_interior_kernel<<<igrid,iblock>>>(t, uarray, farray, is, js,
                                              nx_loc, ny_loc, dx, dy, kx, ky,
                                              forcing);

  // Wait for exchange receives
  flag = udata->end_exchange();
  if (check_flag(&flag, "UserData::end_exchagne", 1)) return -1;

  // Update subdomain boundary
  const realtype *Warray = (udata->HaveNbrW) ? udata->Wrecv : NULL;
  const realtype *Earray = (udata->HaveNbrE) ? udata->Erecv : NULL;
  const realtype *Sarray = (udata->HaveNbrS) ? udata->Srecv : NULL;
  const realtype *Narray = (udata->HaveNbrN) ? udata->Nrecv : NULL;

  sunindextype maxdim = max(nx_loc, ny_loc);
  dim3 bblock(BLOCK_SIZE);
  dim3 bgrid(ICEIL(maxdim, BLOCK_SIZE));

  diffusion_boundary_kernel<<<bgrid,bblock>>>(t, uarray, farray, is, js,
                                              nx_loc, ny_loc, dx, dy, kx, ky,
                                              forcing, Warray, Earray, Sarray,
                                              Narray);

  // Return success
  return 0;
}

#if defined(USE_SUPERLU_DIST)
// Compute the global index of a node
//   i -- local x index
//   j -- local y index
//   x -- x processor coordinate
//   y -- y processor coordinate
sunindextype global_index(sunindextype i, sunindextype j, int x, int y,
                          UserData* udata)
{
  // Unpack values (same for all ranks)
  sunindextype qx = udata->qx;
  sunindextype qy = udata->qy;
  sunindextype rx = udata->rx;
  sunindextype ry = udata->ry;
  sunindextype nx = udata->nx;

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

  // nx_loc for the input (x,y) process (not necessarily the same as nx_loc in
  // user data because this could be the index for a neighboring rank)
  sunindextype nx_loc;
  if (x < rx) nx_loc = qx + 1;
  else nx_loc = qx;

  // global index for the requested local node
  return offset_p + offset_c + j * nx_loc + i;
}

// Compute the global column indices for a row
//   i -- local x index
//   j -- local y index
//   x -- x processor coordinate
//   y -- y processor coordinate
int matrix_columns(sunindextype i, sunindextype j, int x, int y,
                   UserData* udata, sunrealtype* vals, sunindextype* col_idx,
                   sunindextype* row_nnz)
{
  // Unpack values
  int npx = udata->npx;
  int npy = udata->npy;

  sunindextype qx = udata->qx;
  sunindextype qy = udata->qy;
  sunindextype rx = udata->rx;
  sunindextype ry = udata->ry;

  // -----------------
  // global boundaries
  // -----------------

  if ((y == 0 && j == 0) ||                       // south
      (x == 0 && i == 0) ||                       // west
      (x == npx - 1 && i == udata->nx_loc - 1) || // east
      (y == npy - 1 && j == udata->ny_loc - 1))   // north
  {
    vals[0]    = 0;
    col_idx[0] = global_index(i, j, x, y, udata);
    *row_nnz   = 1;
    return 0;
  }

  // --------
  // interior
  // --------

  // Ordering of columns for output
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

  // Constants for computing Jacobian
  sunrealtype cx = udata->kx / (udata->dx * udata->dx);
  sunrealtype cy = udata->ky / (udata->dy * udata->dy);
  sunrealtype cc = -TWO * (cx + cy);

  // Value and Columns index
  sunindextype idx = 0;

  // south neighbor
  if (j == 0)
  {
    // neighbor ny_loc
    sunindextype ny_loc;
    if (y - 1 < ry) ny_loc = qy + 1;
    else ny_loc = qy;

    vals[idx]    = cy;
    col_idx[idx] = global_index(i, ny_loc - 1, x, y - 1, udata);
    idx += 1;
  }

  // west neighbor
  if (i == 0)
  {
    // neighbor nx_loc
    sunindextype nx_loc;
    if (x - 1 < rx) nx_loc = qx + 1;
    else nx_loc = qx;

    vals[idx]    = cx;
    col_idx[idx] = global_index(nx_loc - 1, j, x - 1, y, udata);
    idx += 1;
  }

  // south interior
  if (j > 0)
  {
    vals[idx]    = cy;
    col_idx[idx] = global_index(i, j - 1, x, y, udata);
    idx += 1;
  }

  // west interior
  if (i > 0)
  {
    vals[idx]    = cx;
    col_idx[idx] = global_index(i - 1, j, x, y, udata);
    idx += 1;
  }

  // center
  vals[idx]    = cc;
  col_idx[idx] = global_index(i, j, x, y, udata);
  idx += 1;

  // east interior
  if (i < udata->nx_loc - 1)
  {
    vals[idx]    = cx;
    col_idx[idx] = global_index(i + 1, j, x, y, udata);
    idx += 1;
  }

  // north interior
  if (j < udata->ny_loc - 1)
  {
    vals[idx]    = cy;
    col_idx[idx] = global_index(i, j + 1, x, y, udata);
    idx += 1;
  }

  // east neighbor
  if (i == udata->nx_loc - 1)
  {
    vals[idx]    = cx;
    col_idx[idx] = global_index(0, j, x + 1, y, udata);
    idx += 1;
  }

  // north neighbor
  if (j == udata->ny_loc - 1)
  {
    vals[idx]    = cy;
    col_idx[idx] = global_index(i, 0, x, y + 1, udata);
    idx += 1;
  }

  // non-zeros in this row
  *row_nnz = 5;

  return 0;
}

int laplacian_matrix_sludist(N_Vector u, SUNMatrix L, UserData* udata)
{
  // Set shortcuts
  SuperMatrix*  Lsuper   = SUNMatrix_SLUNRloc_SuperMatrix(L);
  NRformat_loc* Lstore   = (NRformat_loc*)Lsuper->Store;
  sunindextype* row_ptrs = Lstore->rowptr;
  sunindextype* col_idxs = Lstore->colind;
  sunrealtype*  data     = (sunrealtype*)Lstore->nzval;

  int x = udata->idx;
  int y = udata->idy;

  // Initial row pointer value
  row_ptrs[0] = 0;

  sunindextype idx     = 0;
  sunindextype row_nnz = 0;

  for (sunindextype j = 0; j < udata->ny_loc; j++)
  {
    for (sunindextype i = 0; i < udata->nx_loc; i++)
    {
      matrix_columns(i, j, x, y, udata, data + row_ptrs[idx],
                     col_idxs + row_ptrs[idx], &row_nnz);
      row_ptrs[idx + 1] = row_ptrs[idx] + row_nnz;
      idx++;
    }
  }

  Lstore->nnz_loc = row_ptrs[idx];
  Lstore->m_loc   = udata->nx_loc * udata->ny_loc;
  Lstore->fst_row = global_index(0, 0, x, y, udata);

  for (int ip = 0; ip < udata->np; ip++)
  {
    if (ip == udata->myid_c)
    {
      std::cout << "My ID   = " << udata->myid_c << "\n";
      std::cout << "x coord = " << udata->idx << "\n";
      std::cout << "y coord = " << udata->idy << "\n";

      std::cout << "nx_loc  = " << udata->nx_loc << "\n";
      std::cout << "ny_loc  = " << udata->ny_loc << "\n";
      std::cout << "nnz_loc = " << Lstore->nnz_loc << "\n";
      std::cout << "m_loc   = " << Lstore->m_loc << "\n";
      std::cout << "fst_row = " << Lstore->fst_row << "\n";

      std::cout << "row_ptrs = ";
      for (int i = 0; i <= udata->nx_loc * udata->ny_loc; i++)
      {
        std::cout << row_ptrs[i] << " ";
      }
      std::cout << std::endl << std::endl;

      std::cout << "col_idxs = ";
      for (int i = 0; i < Lstore->nnz_loc; i++)
      {
        std::cout << col_idxs[i] << " ";
      }
      std::cout << std::endl << std::endl;

      std::cout << "data = ";
      for (int i = 0; i < Lstore->nnz_loc; i++)
      {
        std::cout << data[i] << " ";
      }
      std::cout << std::endl << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  // Return success
  return 0;
}
#endif
