/*!
 * \file numerics_direct_maxwell.cpp
 * \brief This file contains the Maxwell curl (vorticity) term discretization by utility of flux-splitting method.
 * \author Wenyin Wei 
 * \version 6.2.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#include "../include/numerics_structure.hpp"
#include <limits>

const unsigned short MAXW_EM_DIM = 6;
const bool BOOL_POSITIVE = true;
const bool BOOL_NEGATIVE = false;

/*!
  * \brief Compute the Aofn, short for $\tilde{A}(n)^{+}$ or $\tilde{A}(n)^{-}$. 
  * The dimension is 6 because Maxwell curl operator only works in 3D. One may improve it to be compatible with both simplified 2D and 3D. 
  * \param[in] eps - Electric permittivity of the material. 
  * \param[in] mu - Magnetic peameability of the material. 
  * \param[in] n - The unit direction vector of the cell face.
  * \param[in] is_positive - Whether the user wants a $\tilde{A}(n)^{+}$ or $\tilde{A}(n)^{-}$.
  * \param[out] Aofn - Pointer to the result.
*/
void Compute_Aofn(su2double& eps, su2double& mu, su2double*& n, bool is_positive, su2double** Aofn)
{
  
  su2double v = 1/sqrt(eps*mu);

  if (is_positive)
  {
    Aofn[0][0]=(pow(n[2],2)+pow(n[1],2))*v; Aofn[0][1]=-n[0]*n[1]*v;                Aofn[0][2]=-n[0]*n[2]*v;                Aofn[0][3]=0;                           Aofn[0][4]=n[2]/eps;                    Aofn[0][5]=-n[1]/eps;
    Aofn[1][0]=-n[0]*n[1]*v ;               Aofn[1][1]=(pow(n[0],2)+pow(n[2],2))*v; Aofn[1][2]=-n[1]*n[2]*v ;               Aofn[1][3]=-n[2]/eps ;                  Aofn[1][4]=0;                           Aofn[1][5]=n[0]/eps;
    Aofn[2][0]=-n[0]*n[2]*v ;               Aofn[2][1]=-n[1]*n[2]*v ;               Aofn[2][2]=(pow(n[0],2)+pow(n[1],2))*v; Aofn[2][3]=n[1]/eps ;                   Aofn[2][4]=-n[0]/eps;                   Aofn[2][5]=0;
    Aofn[3][0]=0 ;                          Aofn[3][1]=-n[2]/mu ;                   Aofn[3][2]=n[1]/mu;                     Aofn[3][3]=(pow(n[2],2)+pow(n[1],2))*v; Aofn[3][4]=-n[0]*n[1]*v;                Aofn[3][5]=-n[0]*n[2]*v;
    Aofn[4][0]=n[2]/mu ;                    Aofn[4][1]=0 ;                          Aofn[4][2]=-n[0]/mu;                    Aofn[4][3]=-n[0]*n[1]*v;                Aofn[4][4]=(pow(n[0],2)+pow(n[2],2))*v; Aofn[4][5]=-n[1]*n[2]*v;
    Aofn[5][0]=-n[1]/mu ;                   Aofn[5][1]=n[0]/mu ;                    Aofn[5][2]=0;                           Aofn[5][3]=-n[0]*n[2]*v;                Aofn[5][4]=-n[1]*n[2]*v;                Aofn[5][5]=(pow(n[0],2)+pow(n[1],2))*v;
  }
  else
  {
    Aofn[0][0]=-(pow(n[2],2)+pow(n[1],2))*v;Aofn[0][1]=n[0]*n[1]*v;                 Aofn[0][2]=n[0]*n[2]*v;                 Aofn[0][3]=0;                           Aofn[0][4]=n[2]/eps;                    Aofn[0][5]=-n[1]/eps;
    Aofn[1][0]=n[0]*n[1]*v ;                Aofn[1][1]=-(pow(n[0],2)+pow(n[2],2))*v;Aofn[1][2]=n[1]*n[2]*v ;                Aofn[1][3]=-n[2]/eps ;                  Aofn[1][4]=0;                           Aofn[1][5]=n[0]/eps;
    Aofn[2][0]=n[0]*n[2]*v ;                Aofn[2][1]=n[1]*n[2]*v ;                Aofn[2][2]=-(pow(n[0],2)+pow(n[1],2))*v;Aofn[2][3]=n[1]/eps ;                   Aofn[2][4]=-n[0]/eps;                   Aofn[2][5]=0;
    Aofn[3][0]=0 ;                          Aofn[3][1]=-n[2]/mu ;                   Aofn[3][2]=n[1]/mu;                     Aofn[3][3]=-(pow(n[2],2)+pow(n[1],2))*v; Aofn[3][4]=n[0]*n[1]*v;                Aofn[3][5]=n[0]*n[2]*v;
    Aofn[4][0]=n[2]/mu ;                    Aofn[4][1]=0 ;                          Aofn[4][2]=-n[0]/mu;                    Aofn[4][3]=n[0]*n[1]*v;                Aofn[4][4]=-(pow(n[0],2)+pow(n[2],2))*v; Aofn[4][5]=n[1]*n[2]*v;
    Aofn[5][0]=-n[1]/mu ;                   Aofn[5][1]=n[0]/mu ;                    Aofn[5][2]=0;                           Aofn[5][3]=n[0]*n[2]*v;                Aofn[5][4]=n[1]*n[2]*v;                  Aofn[5][5]=-(pow(n[0],2)+pow(n[1],2))*v;
  }
}


CSourceFluxSplit_Maxwell::CSourceFluxSplit_Maxwell(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  implicit        = (config->GetKind_TimeIntScheme_Maxwell() == EULER_IMPLICIT); //TODO, Check it works 

  Edge_Vector = new su2double [nDim];

  // Aofn_i, short for $\tilde{A}(n)^{+}$ or $\tilde{A}(n)^{-}$. 
  su2double** Aofn_i = new su2double* [MAXW_EM_DIM];
  su2double** Aofn_j = new su2double* [MAXW_EM_DIM];
  for (unsigned short iVar = 0; iVar < MAXW_EM_DIM; iVar++) 
  {
    Aofn_i[iVar] = new su2double [MAXW_EM_DIM]();
    Aofn_j[iVar] = new su2double [MAXW_EM_DIM]();
  }
}

CSourceFluxSplit_Maxwell::~CSourceFluxSplit_Maxwell(void) {

  delete [] Edge_Vector;
  for (iVar = 0; iVar < MAXW_EM_DIM; iVar++)
    delete [] Aofn_i[iVar]; 
  delete [] Aofn_i; 

}

void CSourceFluxSplit_Maxwell::ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config) {

  AD::StartPreacc();
  AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(Temp_i); AD::SetPreaccIn(Temp_j);
  AD::SetPreaccIn(Maxwell_U_i, nVar); AD::SetPreaccIn(Maxwell_U_j, nVar);
  AD::SetPreaccIn(ConsVar_Grad_i, nVar, nDim); AD::SetPreaccIn(ConsVar_Grad_j, nVar, nDim);
  AD::SetPreaccIn(Thermal_Diffusivity_i); AD::SetPreaccIn(Thermal_Conductivity_j);

  Thermal_Diffusivity_Mean = 0.5*(Thermal_Diffusivity_i + Thermal_Diffusivity_j);

  su2double* Temp_Six_Vector_i = new su2double [MAXW_EM_DIM]();
  su2double* Temp_Six_Vector_j = new su2double [MAXW_EM_DIM]();


  su2double eps_i, eps_j, mu_i,  mu_j;// permittivity, permeability
  su2double *n = Normal;
  Compute_Aofn(eps_i, mu_i, n, BOOL_POSITIVE, Aofn_i);
  Compute_Aofn(eps_j, mu_j, n, BOOL_NEGATIVE, Aofn_j);
  su2double Y_i = sqrt(eps_i/mu_i), Z_i =1/Y_i;
  su2double Y_j = sqrt(eps_j/mu_j), Z_j =1/Y_j;

  /*--- Compute vector going from iPoint to jPoint ---*/
  dist_ij_2 = 0; proj_vector_ij = 0; area_face = 0.0; 
  for (iDim = 0; iDim < nDim; iDim++) {
    area_face += Normal[iDim]*Normal[iDim];
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  };
  if (dist_ij_2 == 0.0) {proj_vector_ij = 0.0;}
  else proj_vector_ij = proj_vector_ij/dist_ij_2;
  if (area_face != 0) = sqrt(area_face); 

  /*--- Compute A*U matrix-vector product ---*/
  for (iVar = 0; iVar < MAXW_EM_DIM; iVar++) 
    for (jVar = 0; jVar < MAXW_EM_DIM; iVar++) 
    {
      Temp_Six_Vector_i[iVar] += Aofn_i[iVar][jVar] * Maxwell_U_i[jVar];
      Temp_Six_Vector_j[iVar] += Aofn_j[iVar][jVar] * Maxwell_U_j[jVar];
    }


  for (iVar = 0; iVar < MAXW_EM_DIM; iVar++) 
    val_residual[iVar] = *(Temp_Six_Vector_i[iVar]+Temp_Six_Vector_j[iVar]);
  

  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
  // TODO: Maxwell equation evolution have not yet implemented the flux contribution to Jacobian
  // if (implicit) {
  //   Jacobian_i[0][0] = -Thermal_Diffusivity_Mean*proj_vector_ij;
  //   Jacobian_j[0][0] = Thermal_Diffusivity_Mean*proj_vector_ij;
  // }

  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();

}


CSourceFluxSplitCorrected_Maxwell::CSourceFluxSplitCorrected_Maxwell(unsigned short val_nDim, unsigned short val_nVar,
                                                   CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  implicit        = (config->GetKind_TimeIntScheme_Maxwell() == EULER_IMPLICIT); //TODO, Check the config function is modified correctly
  // TODO, Heat module characters, many heat variables
  Edge_Vector = new su2double [nDim];
  Proj_Mean_GradHeatVar_Edge = new su2double [nVar];
  Proj_Mean_GradHeatVar_Kappa = new su2double [nVar];
  Proj_Mean_GradHeatVar_Corrected = new su2double [nVar];
  Mean_GradHeatVar = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradHeatVar[iVar] = new su2double [nDim];

}

CSourceFluxSplitCorrected_Maxwell::~CSourceFluxSplitCorrected_Maxwell(void) {

  delete [] Edge_Vector;
  delete [] Proj_Mean_GradHeatVar_Edge;
  delete [] Proj_Mean_GradHeatVar_Kappa;
  delete [] Proj_Mean_GradHeatVar_Corrected;
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradHeatVar[iVar];
  delete [] Mean_GradHeatVar;

}

void CSourceFluxSplitCorrected_Maxwell::ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config) {

  // TODO: No idea how to set an appropriate preprocess of adaptive mesh for Maxwell
  // AD::StartPreacc();
  // AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
  // AD::SetPreaccIn(Normal, nDim);
  // AD::SetPreaccIn(Temp_i); AD::SetPreaccIn(Temp_j);
  // AD::SetPreaccIn(ConsVar_Grad_i[0],nDim); AD::SetPreaccIn(ConsVar_Grad_j[0],nDim);
  // AD::SetPreaccIn(Thermal_Diffusivity_i); AD::SetPreaccIn(Thermal_Diffusivity_j);

  Thermal_Diffusivity_Mean = 0.5*(Thermal_Diffusivity_i + Thermal_Diffusivity_j);

  /*--- Compute vector going from iPoint to jPoint ---*/

  dist_ij_2 = 0; proj_vector_ij = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  }
  if (dist_ij_2 == 0.0) proj_vector_ij = 0.0;
  else proj_vector_ij = proj_vector_ij/dist_ij_2;

  /*--- Mean gradient approximation. Projection of the mean gradient
   in the direction of the edge ---*/
  // TODO, Heat module characters
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradHeatVar_Edge[iVar] = 0.0;
    Proj_Mean_GradHeatVar_Kappa[iVar] = 0.0;

    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradHeatVar[iVar][iDim] = 0.5*(ConsVar_Grad_i[iVar][iDim] + ConsVar_Grad_j[iVar][iDim]);
      Proj_Mean_GradHeatVar_Kappa[iVar] += Mean_GradHeatVar[iVar][iDim]*Normal[iDim];
      Proj_Mean_GradHeatVar_Edge[iVar] += Mean_GradHeatVar[iVar][iDim]*Edge_Vector[iDim];
    }
    Proj_Mean_GradHeatVar_Corrected[iVar] = Proj_Mean_GradHeatVar_Kappa[iVar];
    Proj_Mean_GradHeatVar_Corrected[iVar] -= Proj_Mean_GradHeatVar_Edge[iVar]*proj_vector_ij -
    (Temp_j-Temp_i)*proj_vector_ij;
  }

  val_residual[0] = Thermal_Diffusivity_Mean*Proj_Mean_GradHeatVar_Corrected[0];

  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/

  if (implicit) {
    Jacobian_i[0][0] = -Thermal_Diffusivity_Mean*proj_vector_ij;
    Jacobian_j[0][0] = Thermal_Diffusivity_Mean*proj_vector_ij;
  }

  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();
}
