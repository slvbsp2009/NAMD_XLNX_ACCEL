/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

#include "InfoStream.h"
#include "LJTable.h"
#include "Node.h"
#include "Parameters.h"
// #define DEBUGM
#include "Debug.h"
#include "Molecule.h"

//----------------------------------------------------------------------  
LJTable::LJTable()
{
  Bool soluteScalingOn = Node::Object()->simParameters->soluteScalingOn;

  if (!soluteScalingOn) {
  table_dim = Node::Object()->parameters->get_num_vdw_params();
  } else {
  int ss_dim = Node::Object()->molecule->ss_num_vdw_params;
  table_dim = ss_dim + Node::Object()->parameters->get_num_vdw_params();
  }
  table_alloc = new char[2*table_dim*table_dim*sizeof(TableEntry) + 31];
  char *table_align = table_alloc;
  while ( (long)table_align % 32 ) table_align++;
  table = (TableEntry *) table_align;

  for (register int i=0; i < table_dim; i++)
    for (register int j=i; j < table_dim; j++)
    {
      TableEntry *curij = &(table[2*(i*table_dim+j)]);
      TableEntry *curji = &(table[2*(j*table_dim+i)]);
      compute_vdw_params(i,j,curij,curij+1);

      // Copy to transpose entry
      *curji = *curij;
      *(curji + 1) = *(curij + 1);
    }

}

//----------------------------------------------------------------------  
LJTable::~LJTable()
{
  delete [] table_alloc;
}

//----------------------------------------------------------------------
void LJTable::compute_vdw_params(int i, int j,
				 LJTable::TableEntry *cur, 
				 LJTable::TableEntry *cur_scaled)
{
  Parameters *params = Node::Object()->parameters;
  SimParameters *simParams = Node::Object()->simParameters;
  int useGeom = simParams->vdwGeometricSigma;
  Bool tabulatedEnergies = simParams->tabulatedEnergies;
  Bool soluteScalingOn = simParams->soluteScalingOn;
  BigReal soluteScalingFactor = simParams->soluteScalingFactorVdw;
  unsigned int table_dim_org = params->get_num_vdw_params();
  int ss_dim = Node::Object()->molecule->ss_num_vdw_params;
  Real A, B, A14, B14;
  int K = -1;
  int *ss_vdw_type = Node::Object()->molecule->ss_vdw_type;
  // BigReal sigma_max;
  //  We need the A and B parameters for the Van der Waals.  These can
  //  be explicitly be specified for this pair or calculated from the
  //  sigma and epsilon values for the two atom types
//  printf("Looking at interaction of  %i with %i\n", i, j);
  if ( tabulatedEnergies && params->get_table_pair_params(i,j,&K)) {
//    printf("Making this interaction tabulated. %i %i %i\n", i, j, K);
#ifdef NAMD_CUDA
    NAMD_die("Tabulated energies are not supported in CUDA-enabled NAMD");
#endif
    if ( K < 0 ) NAMD_bug(
        "LJTable::compute_vdw_params: energy table index is negative");

    cur->A = -1 - K;
    cur->B = 0;
    cur_scaled->A = -1 - K;
    cur_scaled->B = 0;
  }
  else if (params->get_vdw_pair_params(i,j, &A, &B, &A14, &B14))
  {
    cur->A = A;
    cur->B = B;
    cur_scaled->A = A14;
    cur_scaled->B = B14;

    if ( tabulatedEnergies && ( cur->A < 0 || cur_scaled->A < 0 ) )
      NAMD_die("LJ A is negative with tabulatedEnergies enabled");

    // BigReal sigma_ij, sigma_ij14;

    // if ((B == 0) || (A/B < 0)) sigma_ij = 0;
    // else sigma_ij = pow((BigReal)(A/B),(BigReal)(1./6.));

    // if ((B14 == 0) || (A14/B14 < 0)) sigma_ij14 = 0;
    // else sigma_ij14 = pow((BigReal)(A14/B14),(BigReal)(1./6.));

    // sigma_max = ( sigma_ij > sigma_ij14 ? sigma_ij : sigma_ij14 );
  }
  else
  {
    //  We didn't find explicit parameters for this pair. So instead,
    //  get the parameters for each atom type separately and use them
    //  to calculate the values we need
    Real sigma_i, sigma_i14, epsilon_i, epsilon_i14;
    Real sigma_j, sigma_j14, epsilon_j, epsilon_j14;

   if (!soluteScalingOn) {

    params->get_vdw_params(&sigma_i, &epsilon_i, &sigma_i14,
				       &epsilon_i14,i);
    params->get_vdw_params(&sigma_j, &epsilon_j, &sigma_j14, 
				       &epsilon_j14,j);
   } else {
   int i_type = (i >= table_dim_org)? ss_vdw_type[i-table_dim_org]:i;
   int j_type = (j >= table_dim_org)? ss_vdw_type[j-table_dim_org]:j;
   params->get_vdw_params(&sigma_i, &epsilon_i, &sigma_i14,
                                       &epsilon_i14,i_type);
   params->get_vdw_params(&sigma_j, &epsilon_j, &sigma_j14,
                                       &epsilon_j14,j_type);
   }  	
    BigReal sigma_ij =
       useGeom ? sqrt(sigma_i*sigma_j) : 0.5*(sigma_i+sigma_j);
    BigReal sigma_ij14 =
       useGeom ? sqrt(sigma_i14*sigma_j14) : 0.5 * (sigma_i14+sigma_j14);
    BigReal epsilon_ij = sqrt(epsilon_i*epsilon_j);
    BigReal epsilon_ij14 = sqrt(epsilon_i14*epsilon_j14);

    // sigma_max = ( sigma_ij > sigma_ij14 ? sigma_ij : sigma_ij14 );

    //  Calculate sigma^6
    sigma_ij *= sigma_ij*sigma_ij;
    sigma_ij *= sigma_ij;
    sigma_ij14 *= sigma_ij14*sigma_ij14;
    sigma_ij14 *= sigma_ij14;
    
    //  Calculate LJ constants A & B
    cur->B = 4.0 * sigma_ij * epsilon_ij;
    cur->A = cur->B * sigma_ij;
    cur_scaled->B = 4.0 * sigma_ij14 * epsilon_ij14;
    cur_scaled->A = cur_scaled->B * sigma_ij14;

    if (soluteScalingOn) {
     if (i >= table_dim_org && i < (table_dim_org+ss_dim) && j < table_dim_org) {
        cur->A *= sqrt(soluteScalingFactor);
        cur->B *= sqrt(soluteScalingFactor);
        cur_scaled->A *= sqrt(soluteScalingFactor);
        cur_scaled->B *= sqrt(soluteScalingFactor);
     }
     if (i < table_dim_org && j >= table_dim_org && j < (table_dim_org+ss_dim)) {
        cur->A *= sqrt(soluteScalingFactor);
        cur->B *= sqrt(soluteScalingFactor);
        cur_scaled->A *= sqrt(soluteScalingFactor);
        cur_scaled->B *= sqrt(soluteScalingFactor);
      }
     if (i >=table_dim_org && i < (table_dim_org+ss_dim) && j >= table_dim_org && j < (table_dim_org+ss_dim)) {
        cur->A *= soluteScalingFactor;
        cur->B *= soluteScalingFactor;
        cur_scaled->A *= soluteScalingFactor;
        cur_scaled->B *= soluteScalingFactor;
     }
    }

    if ( tabulatedEnergies && ( cur->A < 0 || cur_scaled->A < 0 ) )
      NAMD_die("LJ A is negative with tabulatedEnergies enabled");
  }
  //  Calculate exclcut2
  // cur_scaled->exclcut2 = cur->exclcut2 = 0.64 * sigma_max * sigma_max;

}

