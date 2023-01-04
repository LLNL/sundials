/* -----------------------------------------------------------------------------
 * Programmer(s): Daniel R. Reynolds @ SMU
 *                Cody J. Balos @ LLNL
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
 * A simple implementation of a parallel structured Cartesian mesh class that
 * supports up to 3 spatial dimensions and an arbitrary number of degrees of
 * freedom, and that uses Kokkos views to store communication buffer data.
 * ----------------------------------------------------------------------------*/

#ifndef _KOKKOSPARGRID_H
#define _KOKKOSPARGRID_H

#include <iomanip>
#include <iostream>
#include <fstream>
#include <mpi.h>
#include <Kokkos_Core.hpp>

namespace sundials_tools
{

// Types of boundaries supported.
enum class BoundaryType
{
  PERIODIC
};

// Types of stencils supported.
enum class StencilType
{
  UPWIND
};

template<typename REAL, typename GLOBALINT, int NDIMS>
class ParallelGrid
{
public:
  // Constructor that creates a new ParallelGrid object.
  // [in] - the memory helper to use for allocating the MPI buffers
  // [in,out] comm - on input, the overal MPI communicator, on output, the cartesian communicator
  // [in] a[] - an array of length NDIMS which defines the domain [a,b]
  // [in] b[] - an array of length NDIMS which defines the domain [a,b]
  // [in] npts[] - an array of length NDIMS which defines the number of mesh points in each dimension
  // [in] dof - the number of degrees of freedom in each dimension
  // [in] bc - the type of boundary conditions (see BoundaryType)
  // [in] st - the stencil to use (see StencilType)
  // [in] width - the stencil width; defaults to 1
  // [in] npxyz - the number of processors in each dimension; defaults to 0 which means MPI will choose
  // [in] reorder - should MPI_Cart_create do process reordering to optimize or not; defaults to false (some MPI implementations ignore this)
  ParallelGrid(MPI_Comm* comm, const REAL a[], const REAL b[], const GLOBALINT npts[],
               int dof, BoundaryType bc, StencilType st, const REAL c, int width = 1,
               const int npxyz[] = nullptr, bool reorder = false)
    : nx(1), ny(1), nz(1),
      nxl(1), nyl(1), nzl(1),
      npx(1), npy(1), npz(1),
      dx(0.0), dy(0.0), dz(0.0),
      ax(0.0), ay(0.0), az(0.0),
      bx(0.0), by(0.0), bz(0.0),
      dof(dof), dims{0,0,0}, coords{0,0,0},
      bc(bc), st(st), width(width),
      upwindRight(true)
  {
    static_assert((NDIMS >= 1 && NDIMS <= 3), "ParallelGrid NDIMS must be 1, 2 or 3");
    assert(st == StencilType::UPWIND);

    /* Set up MPI Cartesian communicator */
    if (npxyz)
    {
      dims[0] = npxyz[0];
      if (NDIMS >= 2) dims[1] = npxyz[1];
      if (NDIMS == 3) dims[2] = npxyz[2];
    }

    int retval, nprocs;
    MPI_Comm_size(*comm, &nprocs);
    retval = MPI_Dims_create(nprocs, NDIMS, dims);
    assert(retval == MPI_SUCCESS);

    int periods[] = { bc == BoundaryType::PERIODIC,
                      bc == BoundaryType::PERIODIC,
                      bc == BoundaryType::PERIODIC };
    retval = MPI_Cart_create(*comm, NDIMS, dims, periods, reorder, comm);
    assert(retval == MPI_SUCCESS);

    retval = MPI_Cart_get(*comm, NDIMS, dims, periods, coords);
    assert(retval == MPI_SUCCESS);

    cart_comm = *comm;

    /* Set upwinding direction */
    upwindRight = (c > 0.0);

    /* Set up information for the first spatial dimension */
    npx = dims[0];
    nx  = npts[0];
    ax  = a[0];
    bx  = b[0];
    dx  = (bx-ax) / (REAL) nx;
    int is = nx*(coords[0])/npx;
    int ie = nx*(coords[0]+1)/npx-1;
    nxl = ie-is+1;
    neq = dof * nxl;

    if (NDIMS >= 2)
    {
      /* Set up information for the second spatial dimension (if applicable) */
      npy = dims[1];
      ny  = npts[1];
      ay  = a[1];
      by  = b[1];
      dy  = (by-ay) / (REAL) ny;
      int js = ny*(coords[1])/npy;
      int je = ny*(coords[1]+1)/npy-1;
      nyl = je-js+1;
      neq *= nyl;
    }

    if (NDIMS == 3)
    {
      /* Set up information for the third spatial dimension (if applicable) */
      npz = dims[2];
      nz  = npts[2];
      az  = a[2];
      bz  = b[2];
      dz  = (bz-az) / (REAL) nz;
      int ks = nz*(coords[2])/npz;
      int ke = nz*(coords[2]+1)/npz-1;
      nzl = ke-ks+1;
      neq *= nzl;
    }

    /* Allocate buffers for nearest-neighbor exchange */
    if (st == StencilType::UPWIND)
      AllocateBuffersUpwind();

  }

  // TODO:
  //  - support non-periodic boundary conditions
  // For all faces where neighbors exist: determine neighbor process indices.
  // For all faces: allocate upwind exchange buffers.
  void AllocateBuffersUpwind()
  {

    /* Allocate send/receive buffers and determine ID for communication West */
    if (upwindRight)
      Wrecv_ = Kokkos::View<REAL*>("Wrecv", dof*width*nyl*nzl);
    else
      Wsend_ = Kokkos::View<REAL*>("Wsend", dof*width*nyl*nzl);
    ipW = MPI_PROC_NULL;
    if ((coords[0] > 0) || (bc == BoundaryType::PERIODIC)) {
      int nbcoords[] = {coords[0]-1, coords[1], coords[2]};
      int retval = MPI_Cart_rank(cart_comm, nbcoords, &ipW);
      assert(retval == MPI_SUCCESS);
    }

    /* Allocate send/receive buffers and determine ID for communication East */
    if (upwindRight)
      Esend_ = Kokkos::View<REAL*>("Esend", dof*width*nyl*nzl);
    else
      Erecv_ = Kokkos::View<REAL*>("Erecv", dof*width*nyl*nzl);
    ipE = MPI_PROC_NULL;
    if ((coords[0] < dims[0]-1) || (bc == BoundaryType::PERIODIC)) {
      int nbcoords[] = {coords[0]+1, coords[1], coords[2]};
      int retval = MPI_Cart_rank(cart_comm, nbcoords, &ipE);
      assert(retval == MPI_SUCCESS);
    }

    if (NDIMS >= 2)
    {
      /* Allocate send/receive buffers and determine ID for communication South */
      if (upwindRight)
        Srecv_ = Kokkos::View<REAL*>("Srecv", dof*width*nxl*nzl);
      else
        Ssend_ = Kokkos::View<REAL*>("Ssend", dof*width*nxl*nzl);
      ipS = MPI_PROC_NULL;
      if ((coords[1] > 0) || (bc == BoundaryType::PERIODIC)) {
        int nbcoords[] = {coords[0], coords[1]-1, coords[2]};
        int retval = MPI_Cart_rank(cart_comm, nbcoords, &ipS);
        assert(retval == MPI_SUCCESS);
      }

      /* Allocate send/receive buffers and determine ID for communication North */
      if (upwindRight)
        Nsend_ = Kokkos::View<REAL*>("Nsend", dof*width*nxl*nzl);
      else
        Nrecv_ = Kokkos::View<REAL*>("Nrecv", dof*width*nxl*nzl);
      ipN = MPI_PROC_NULL;
      if ((coords[1] < dims[1]-1) || (bc == BoundaryType::PERIODIC)) {
        int nbcoords[] = {coords[0], coords[1]+1, coords[2]};
        int retval = MPI_Cart_rank(cart_comm, nbcoords, &ipN);
        assert(retval == MPI_SUCCESS);
      }
    }

    if (NDIMS == 3)
    {
      /* Allocate send/receive buffers and determine ID for communication Back */
      if (upwindRight)
        Brecv_ = Kokkos::View<REAL*>("Brecv", dof*width*nxl*nyl);
      else
        Bsend_ = Kokkos::View<REAL*>("Bsend", dof*width*nxl*nyl);
      ipB = MPI_PROC_NULL;
      if ((coords[2] > 0) || (bc == BoundaryType::PERIODIC)) {
        int nbcoords[] = {coords[0], coords[1], coords[2]-1};
        int retval = MPI_Cart_rank(cart_comm, nbcoords, &ipB);
        assert(retval == MPI_SUCCESS);
      }

      /* Allocate send/receive buffers and determine ID for communication Front */
      if (upwindRight)
        Fsend_ = Kokkos::View<REAL*>("Fsend", dof*width*nxl*nyl);
      else
        Frecv_ = Kokkos::View<REAL*>("Frecv", dof*width*nxl*nyl);
      ipF = MPI_PROC_NULL;
      if ((coords[2] < dims[2]-1) || (bc == BoundaryType::PERIODIC)) {
        int nbcoords[] = {coords[0], coords[1], coords[2]+1};
        int retval = MPI_Cart_rank(cart_comm, nbcoords, &ipF);
        assert(retval == MPI_SUCCESS);
      }
    }

  }

  // Initiate non-blocking neighbor communication
  int ExchangeStart()
  {
    int retval = 0;
    nreq = 0;

    // Initialize all requests in array
    for (int i=0; i<12; i++)
      req[i] = MPI_REQUEST_NULL;

    // Open an Irecv buffer for each neighbor
    if ((ipW != MPI_PROC_NULL) && (upwindRight))
    {
      retval = MPI_Irecv(Wrecv_.data(), int(Wrecv_.size()), MPI_SUNREALTYPE, ipW,
                         1, cart_comm, req+nreq);
      assert(retval == MPI_SUCCESS);
      nreq++;
    }

    if ((ipE != MPI_PROC_NULL) && (!upwindRight))
    {
      retval = MPI_Irecv(Erecv_.data(), int(Erecv_.size()), MPI_SUNREALTYPE, ipE,
                         0, cart_comm, req+nreq);
      assert(retval == MPI_SUCCESS);
      nreq++;
    }

    if (NDIMS >= 2)
    {
      if ((ipS != MPI_PROC_NULL) && (upwindRight))
      {
        retval = MPI_Irecv(Srecv_.data(), int(Srecv_.size()), MPI_SUNREALTYPE, ipS,
                           3, cart_comm, req+nreq);
        assert(retval == MPI_SUCCESS);
        nreq++;
      }

      if ((ipN != MPI_PROC_NULL) && (!upwindRight))
      {
        retval = MPI_Irecv(Nrecv_.data(), int(Nrecv_.size()), MPI_SUNREALTYPE, ipN,
                           2, cart_comm, req+nreq);
        assert(retval == MPI_SUCCESS);
        nreq++;
      }
    }

    if (NDIMS >= 3)
    {
      if ((ipB != MPI_PROC_NULL) && (upwindRight))
      {
        retval = MPI_Irecv(Brecv_.data(), int(Brecv_.size()), MPI_SUNREALTYPE, ipB,
                           5, cart_comm, req+nreq);
        assert(retval == MPI_SUCCESS);
        nreq++;
      }

      if ((ipF != MPI_PROC_NULL) && (!upwindRight))
      {
        retval = MPI_Irecv(Frecv_.data(), int(Frecv_.size()), MPI_SUNREALTYPE, ipF,
                           4, cart_comm, req+nreq);
        assert(retval == MPI_SUCCESS);
        nreq++;
      }
    }

    // Send data to neighbors
    if ((ipW != MPI_PROC_NULL) && (!upwindRight))
    {
      retval = MPI_Isend(Wsend_.data(), int(Wsend_.size()), MPI_SUNREALTYPE, ipW, 0,
                         cart_comm, req+nreq);
      assert(retval == MPI_SUCCESS);
      nreq++;
    }

    if ((ipE != MPI_PROC_NULL) && (upwindRight))
    {
      retval = MPI_Isend(Esend_.data(), int(Esend_.size()), MPI_SUNREALTYPE, ipE, 1,
                         cart_comm, req+nreq);
      assert(retval == MPI_SUCCESS);
      nreq++;
    }

    if (NDIMS >= 2)
    {
      if ((ipS != MPI_PROC_NULL) && (!upwindRight))
      {
        retval = MPI_Isend(Ssend_.data(), int(Ssend_.size()), MPI_SUNREALTYPE, ipS, 2,
                           cart_comm, req+nreq);
        assert(retval == MPI_SUCCESS);
        nreq++;
      }

      if ((ipN != MPI_PROC_NULL) && (upwindRight))
      {
        retval = MPI_Isend(Nsend_.data(), int(Nsend_.size()), MPI_SUNREALTYPE, ipN, 3,
                           cart_comm, req+nreq);
        assert(retval == MPI_SUCCESS);
        nreq++;
      }
    }

    if (NDIMS == 3)
    {
      if ((ipB != MPI_PROC_NULL) && (!upwindRight))
      {
        retval = MPI_Isend(Bsend_.data(), int(Bsend_.size()), MPI_SUNREALTYPE, ipB, 4,
                           cart_comm, req+nreq);
        assert(retval == MPI_SUCCESS);
        nreq++;
      }

      if ((ipF != MPI_PROC_NULL) && (upwindRight))
      {
        retval = MPI_Isend(Fsend_.data(), int(Fsend_.size()), MPI_SUNREALTYPE, ipF, 5,
                           cart_comm, req+nreq);
        assert(retval == MPI_SUCCESS);
        nreq++;
      }
    }

    return retval;
  }

  // Waits for neighbor exchange to finish.
  int ExchangeEnd()
  {
    MPI_Status stat[12];
    int retval;

    // return automatically with success if there are no outstanding requests
    if (nreq == 0)
      return(0);

    // Wait for messages to finish send/receive
    retval = MPI_Waitall(nreq, req, stat);
    assert(retval == MPI_SUCCESS);

    return retval;
  }

  // Prints out information about the ParallelGrid to stdout.
  void PrintInfo()
  {
    printf("ParallelGrid Info:\n");
    printf("    dimensions = %d\n", NDIMS);
    printf("    processors = {%d, %d, %d}\n", npx, npy, npz);
    printf("        domain = {[%g,%g], [%g,%g], [%g,%g]}\n", ax, bx, ay, by, az, bz);
    printf("   global npts = {%li, %li, %li}\n", (long int) nx, (long int) ny, (long int) nz);
    printf("    local npts = {%d, %d, %d}\n", nxl, nyl, nzl);
    printf("  mesh spacing = {%g, %g, %g}\n", dx, dy, dz);
    if (upwindRight)
      printf("    upwind dir = right\n");
    else
      printf("    upwind dir = left\n");
  }

  // Saves the mesh to a file.
  //    First row is x. Second row is y. Third row is z.
  //    Can be loaded into MATLAB like so:
  //      mesh = loadtxt('mesh.txt');
  //      [X,Y,Z] = meshgrid(mesh(1,:),mesh(2,:),mesh(3,:));
  void MeshToFile(const std::string& fname)
  {
    std::ofstream mesh_file;
    mesh_file.open(fname);
    mesh_file << std::setprecision(16);
    for (GLOBALINT i = 0; i < nx; i++)
      mesh_file << " " << dx*i;
    mesh_file << std::endl;
    for (GLOBALINT i = 0; i < ny; i++)
      mesh_file << " " << dy*i;
    mesh_file << std::endl;
    for (GLOBALINT i = 0; i < nz; i++)
      mesh_file << " " << dz*i;
    mesh_file << std::endl;
    mesh_file.close();
  }

  int nprocs() const
  {
    return npx*npy*npz;
  }

  GLOBALINT npts() const
  {
    if (NDIMS == 1) return nx;
    if (NDIMS == 2) return nx*ny;
    if (NDIMS == 3) return nx*ny*nz;
  }

  GLOBALINT nptsl() const
  {
    if (NDIMS == 1) return nxl;
    if (NDIMS == 2) return nxl*nyl;
    if (NDIMS == 3) return nxl*nyl*nzl;
  }

  GLOBALINT neql() const
  {
    return dof*nptsl();
  }

  Kokkos::View<REAL****> GetRecvView(const std::string& direction)
  {
    if (direction == "WEST")
    {
      return Kokkos::View<REAL****>(Wrecv_.data(), 1, nyl, nzl, dof);
    }
    else if (direction == "EAST")
    {
      return Kokkos::View<REAL****>(Erecv_.data(), 1, nyl, nzl, dof);
    }
    else if (direction == "NORTH")
    {
      return Kokkos::View<REAL****>(Nrecv_.data(), nxl, 1, nzl, dof);
    }
    else if (direction == "SOUTH")
    {
      return Kokkos::View<REAL****>(Srecv_.data(), nxl, 1, nzl, dof);
    }
    else if (direction == "FRONT")
    {
      return Kokkos::View<REAL****>(Frecv_.data(), nxl, nyl, 1, dof);
    }
    else if (direction == "BACK")
    {
      return Kokkos::View<REAL****>(Brecv_.data(), nxl, nyl, 1, dof);
    }
    else
    {
      assert(direction == "ILLEGAL");
      return Kokkos::View<REAL****>();
    }
  }

  Kokkos::View<REAL****> GetSendView(const std::string& direction)
  {
    if (direction == "WEST")
    {
      return Kokkos::View<REAL****>(Wsend_.data(), 1, nyl, nzl, dof);
    }
    else if (direction == "EAST")
    {
      return Kokkos::View<REAL****>(Esend_.data(), 1, nyl, nzl, dof);
    }
    else if (direction == "NORTH")
    {
      return Kokkos::View<REAL****>(Nsend_.data(), nxl, 1, nzl, dof);
    }
    else if (direction == "SOUTH")
    {
      return Kokkos::View<REAL****>(Ssend_.data(), nxl, 1, nzl, dof);
    }
    else if (direction == "FRONT")
    {
      return Kokkos::View<REAL****>(Fsend_.data(), nxl, nyl, 1, dof);
    }
    else if (direction == "BACK")
    {
      return Kokkos::View<REAL****>(Bsend_.data(), nxl, nyl, 1, dof);
    }
    else
    {
      assert(direction == "ILLEGAL");
      return Kokkos::View<REAL****>();
    }
  }

  GLOBALINT nx, ny, nz;    /* number of intervals globally       */
  int       nxl, nyl, nzl; /* number of intervals locally        */
  int       npx, npy, npz; /* numner of processes                */
  REAL      dx, dy, dz;    /* mesh spacing                       */
  REAL      ax, ay, az;    /* domain in [a, b]                   */
  REAL      bx, by, bz;
  int       dof;           /* degrees of freedom per node        */
  int       neq;           /* total number of equations locally  */

  int       ipW, ipE;      /* MPI ranks for neighbor procs       */
  int       ipS, ipN;
  int       ipB, ipF;
  bool      upwindRight;   /* Upwind dir: true/false == R/L      */

  int       dims[3];
  int       coords[3];


private:
  MPI_Comm     cart_comm;  /* MPI cartesian communicator         */
  MPI_Request  req[12];
  int          nreq;

  BoundaryType bc;
  StencilType  st;
  int          width;

  Kokkos::View<REAL*> Wsend_;            /* MPI send/recv buffers              */
  Kokkos::View<REAL*> Esend_;
  Kokkos::View<REAL*> Ssend_;
  Kokkos::View<REAL*> Nsend_;
  Kokkos::View<REAL*> Bsend_;
  Kokkos::View<REAL*> Fsend_;
  Kokkos::View<REAL*> Wrecv_;
  Kokkos::View<REAL*> Erecv_;
  Kokkos::View<REAL*> Srecv_;
  Kokkos::View<REAL*> Nrecv_;
  Kokkos::View<REAL*> Brecv_;
  Kokkos::View<REAL*> Frecv_;

};

}

#endif
