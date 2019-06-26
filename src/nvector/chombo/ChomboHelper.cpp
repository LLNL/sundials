#include "ChomboHelper.H"


void export_fableveldata(const char* file, const FArrayBoxLD* data)
{
  const DisjointBoxLayout& grids( data->disjointBoxLayout() );
  const ProblemDomain& domain( grids.physDomain() );
  WriteUGHDF5( file, grids, *data, domain.domainBox() );
}

void box_extents(const Box& region, int* ilower, int* iupper)
{
  for (int idir=0; idir<CH_SPACEDIM; idir++)
  {
    ilower[idir] = region.smallEnd(idir);
    iupper[idir] = region.bigEnd(idir);
  }
}

void add_term(const FArrayBoxLD& R, 
        const Real coeff, 
        FArrayBoxLD& G)
{
  const DisjointBoxLayout& dbl(R.disjointBoxLayout());
  DataIterator dit(dbl.dataIterator());
  for (dit.reset(); dit.ok(); ++dit)
  {
    const FArrayBox& R_f(R[dit()]);
    FArrayBox& G_f(G[dit()]);
    G_f.plus(R_f, coeff);
  }
}

Box interiorbox(const Domain& D, 
        const DataIndex& index) 
{
  const DisjointBoxLayout& dbl( D.disjointBoxLayout() );

  Box i_box( dbl[index] );

  SideIterator side;
  for (int idir=0; idir<CH_SPACEDIM; idir++)
  {
    for (side.begin(); side.ok(); ++side)
    {
      if ( D.disjoint_box( adjCellBox(dbl[index], idir, side(), 1) ) )
        i_box.growDir(idir, side(), -1);
    } // side
  } // idir

  return i_box;
}

Box interiorfacebox(const Domain& D, 
        const DataIndex& index, 
        const int idir)
{
  const DisjointBoxLayout& dbl( D.disjointBoxLayout() );

  Box flux_box( dbl[index] );

  flux_box.growHi(idir, 1);
  SideIterator side;
  for (side.begin(); side.ok(); ++side)
  {
    if ( D.disjoint_box(adjCellBox(dbl[index], idir, side(), 1)) )
      flux_box.growDir(idir, side(), -1);
  }

  return flux_box;
}

Box boundarybox(const FArrayBox& F, 
        const IntVect& ghost, 
        const int idir, 
        const Side::LoHiSide& iside)
{
  Box face( F.box() );
  face.grow( -ghost );
  if (iside == Side::Lo)
    face.setBig(idir, face.smallEnd(idir));
  else
    face.setSmall(idir, face.bigEnd(idir));
  return face;
}

Box ghostbox(const FArrayBox& F,
        const IntVect& ghost,
        const int idir,
        const Side::LoHiSide& iside)
{
  // NOTE: only works with one layer of ghost cells
  Box gbox( F.box() );
  if (iside == Side::Lo)
    gbox.setBig(idir, gbox.smallEnd(idir));
  else
    gbox.setSmall(idir, gbox.bigEnd(idir));
  return gbox;
}

