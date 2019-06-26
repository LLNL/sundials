#ifndef _CHOMBO_HELPER_H
#define _CHOMBO_HELPER_H

#include "ChomboBoxTools.H"
#include "Domain.H"

void export_fableveldata(const char* file, const FArrayBoxLD* data);

void box_extents(const Box& region, int* ilower, int* iupper);

void add_term(const FArrayBoxLD& R, 
        const Real coeff, 
        FArrayBoxLD& G);


Box interiorbox(const Domain& D, 
        const DataIndex& index);

Box interiorfacebox(const Domain& D, 
        const DataIndex& index, 
        const int idir);

Box boundarybox(const FArrayBox& F, 
        const IntVect& ghost, 
        const int idir, 
        const Side::LoHiSide& iside);

Box ghostbox(const FArrayBox& F,
        const IntVect& ghost,
        const int idir,
        const Side::LoHiSide& iside);

#endif
