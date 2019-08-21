/*!
 * \file numerics_direct_maxwell.cpp
 * \brief This file contains all the convective term discretization.
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

CUpwSca_Maxwell::CUpwSca_Maxwell(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  implicit        = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  grid_movement   = config->GetGrid_Movement();

  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];

  Laminar_Viscosity_i = config->GetViscosity_FreeStreamND();
  Laminar_Viscosity_j = config->GetViscosity_FreeStreamND();
}

CUpwSca_Maxwell::~CUpwSca_Maxwell(void) {

  delete [] Velocity_i;
  delete [] Velocity_j;

}

void CUpwSca_Maxwell::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {

  q_ij = 0.0;

  AD::StartPreacc();
  AD::SetPreaccIn(V_i, nDim+1); AD::SetPreaccIn(V_j, nDim+1);
  AD::SetPreaccIn(Temp_i); AD::SetPreaccIn(Temp_j);
  AD::SetPreaccIn(Normal, nDim);
  if (grid_movement) {
    AD::SetPreaccIn(GridVel_i, nDim); AD::SetPreaccIn(GridVel_j, nDim);
  }

  for (iDim = 0; iDim < nDim; iDim++) {
    Velocity_i[iDim] = V_i[iDim+1];
    Velocity_j[iDim] = V_j[iDim+1];
    q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
  }

  a0 = 0.5*(q_ij+fabs(q_ij));
  a1 = 0.5*(q_ij-fabs(q_ij));
  val_residual[0] = a0*Temp_i+a1*Temp_j;

  if (implicit) {
    val_Jacobian_i[0][0] = a0;
    val_Jacobian_j[0][0] = a1;
  }

  AD::SetPreaccOut(val_residual[0]);
  AD::EndPreacc();

}

CCentSca_Maxwell::CCentSca_Maxwell(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  implicit        = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  grid_movement   = config->GetGrid_Movement();

  MeanVelocity = new su2double [nDim];

  Laminar_Viscosity_i = config->GetViscosity_FreeStreamND(); //TODO, Heat module characters
  Laminar_Viscosity_j = config->GetViscosity_FreeStreamND(); //TODO, Heat module characters

  Param_Kappa_4 = config->GetKappa_4th_Heat(); //TODO, Heat module characters
  Diff_Lapl = new su2double [nVar];
}

CCentSca_Maxwell::~CCentSca_Maxwell(void) {

  delete [] MeanVelocity;
  delete [] Diff_Lapl;

}

void CCentSca_Maxwell::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {

  AD::StartPreacc();
  AD::SetPreaccIn(V_i, nDim+3); AD::SetPreaccIn(V_j, nDim+3);
  AD::SetPreaccIn(Temp_i); AD::SetPreaccIn(Temp_j);
  AD::SetPreaccIn(Und_Lapl_i, nVar); AD::SetPreaccIn(Und_Lapl_j, nVar);
  AD::SetPreaccIn(Normal, nDim);
  if (grid_movement) {
    AD::SetPreaccIn(GridVel_i, nDim); AD::SetPreaccIn(GridVel_j, nDim);
  }

  /*--- Primitive variables at point i and j ---*/

  Pressure_i =    V_i[0];       Pressure_j = V_j[0];
  DensityInc_i =  V_i[nDim+2];  DensityInc_j = V_j[nDim+2];
  BetaInc2_i =    V_i[nDim+3];  BetaInc2_j = V_j[nDim+3];

  /*--- Projected velocities at the current edge ---*/

  ProjVelocity = 0.0; ProjVelocity_i = 0.0; ProjVelocity_j = 0.0; Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    ProjVelocity_i += V_i[iDim+1]*Normal[iDim];
    ProjVelocity_j += V_j[iDim+1]*Normal[iDim];
    Area += Normal[iDim]*Normal[iDim];
  }
  Area = sqrt(Area);

  /*--- Computing the second order centered scheme part ---*/

  ProjVelocity = 0.5*(ProjVelocity_i+ProjVelocity_j);

  val_residual[0] = 0.5*(Temp_i + Temp_j)*ProjVelocity;

  if (implicit) {
    val_Jacobian_i[0][0] = 0.5*ProjVelocity;
    val_Jacobian_j[0][0] = 0.5*ProjVelocity;
  }

  /*--- Adding artificial dissipation to stabilize the centered scheme ---*/

  Diff_Lapl[0] = Und_Lapl_i[0]-Und_Lapl_j[0];

  SoundSpeed_i = sqrt(ProjVelocity_i*ProjVelocity_i + (BetaInc2_i/DensityInc_i)*Area*Area);
  SoundSpeed_j = sqrt(ProjVelocity_j*ProjVelocity_j + (BetaInc2_j/DensityInc_j)*Area*Area);

  Local_Lambda_i = fabs(ProjVelocity_i)+SoundSpeed_i;
  Local_Lambda_j = fabs(ProjVelocity_j)+SoundSpeed_j;
  MeanLambda = 0.5*(Local_Lambda_i+Local_Lambda_j);

  val_residual[0] += -Param_Kappa_4*Diff_Lapl[0]*MeanLambda;

  if (implicit) {
    cte_0 = Param_Kappa_4*su2double(Neighbor_i+1)*MeanLambda;
    cte_1 = Param_Kappa_4*su2double(Neighbor_j+1)*MeanLambda;

    val_Jacobian_i[0][0] += cte_0;
    val_Jacobian_j[0][0] -= cte_1;
  }

  AD::SetPreaccOut(val_residual[0]);
  AD::EndPreacc();

}

const unsigned short MAXW_EM_DIM = 6;
CAvgGrad_Maxwell::CAvgGrad_Maxwell(unsigned short val_nDim, unsigned short val_nVar, CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  implicit        = (config->GetKind_TimeIntScheme_Heat() == EULER_IMPLICIT); //TODO, Heat module characters

  Edge_Vector = new su2double [nDim];
  Proj_Mean_GradHeatVar_Normal = new su2double [nVar]; //TODO, Heat module characters
  Proj_Mean_GradHeatVar_Corrected = new su2double [nVar]; //TODO, Heat module characters
  Mean_GradHeatVar = new su2double* [nVar]; //TODO, Heat module characters
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradHeatVar[iVar] = new su2double [nDim]; //TODO, Heat module characters

  Aofn = new su2double* [MAXW_EM_DIM]; // Aofn, short for \tilde{A}(n)^{+}, The dimension is 6 because Maxwell curl operator only works in 3D
  for (iVar = 0; iVar < MAXW_EM_DIM; iVar++)
    Aofn[iVar] = new su2double [MAXW_EM_DIM]; 
  

}

CAvgGrad_Maxwell::~CAvgGrad_Maxwell(void) {

  delete [] Edge_Vector;
  delete [] Proj_Mean_GradHeatVar_Normal; //TODO, Heat module characters
  delete [] Proj_Mean_GradHeatVar_Corrected; //TODO, Heat module characters
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradHeatVar[iVar]; //TODO, Heat module characters
  delete [] Mean_GradHeatVar; //TODO, Heat module characters
  for (iVar = 0; iVar < MAXW_EM_DIM; iVar++)
    delete [] Aofn[iVar]; 
  delete [] Aofn; 

}

void CAvgGrad_Maxwell::ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config) {

  AD::StartPreacc();
  AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(Temp_i); AD::SetPreaccIn(Temp_j);
  AD::SetPreaccIn(ConsVar_Grad_i[0],nDim); AD::SetPreaccIn(ConsVar_Grad_j[0],nDim);
  AD::SetPreaccIn(Thermal_Diffusivity_i); AD::SetPreaccIn(Thermal_Conductivity_j);

  Thermal_Diffusivity_Mean = 0.5*(Thermal_Diffusivity_i + Thermal_Diffusivity_j);

  

  su2double eps ;// permittivity
  su2double mu  ;// permeability
  su2double v = 1/sqrt(eps*mu);
  Aofn[0][0]=(pow(n[2],2)+pow(n[1],2))*v; Aofn[0][1]=-n[0]*n[1]*v;                Aofn[0][2]=-n[0]*n[2]*v;                Aofn[0][3]=0;                           Aofn[0][4]=n[2]/eps;                    Aofn[0][5]=-n[1]/eps;
  Aofn[1][0]=-n[0]*n[1]*v ;               Aofn[1][1]=(pow(n[0],2)+pow(n[2],2))*v; Aofn[1][2]=-n[1]*n[2]*v ;               Aofn[1][3]=-n[2]/eps ;                  Aofn[1][4]=0;                           Aofn[1][5]=n[0]/eps;
  Aofn[2][0]=-n[0]*n[2]*v ;               Aofn[2][1]=-n[1]*n[2]*v ;               Aofn[2][2]=(pow(n[0],2)+pow(n[1],2))*v; Aofn[2][3]=n[1]/eps ;                   Aofn[2][4]=-n[0]/eps;                   Aofn[2][5]=0;
  Aofn[3][0]=0 ;                          Aofn[3][1]=-n[2]/mu ;                   Aofn[3][2]=n[1]/mu;                     Aofn[3][3]=(pow(n[2],2)+pow(n[1],2))*v; Aofn[3][4]=-n[0]*n[1]*v;                Aofn[3][5]=-n[0]*n[2]*v;
  Aofn[4][0]=n[2]/mu ;                    Aofn[4][1]=0 ;                          Aofn[4][2]=-n[0]/mu;                    Aofn[4][3]=-n[0]*n[1]*v;                Aofn[4][4]=(pow(n[0],2)+pow(n[2],2))*v; Aofn[4][5]=-n[1]*n[2]*v;
  Aofn[5][0]=-n[1]/mu ;                   Aofn[5][1]=n[0]/mu ;                    Aofn[5][2]=0;                           Aofn[5][3]=-n[0]*n[2]*v;                Aofn[5][4]=-n[1]*n[2]*v;                Aofn[5][5]=(pow(n[0],2)+pow(n[1],2))*v;

  /*--- Compute vector going from iPoint to jPoint ---*/

  dist_ij_2 = 0; proj_vector_ij = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  }
  if (dist_ij_2 == 0.0) proj_vector_ij = 0.0;
  else proj_vector_ij = proj_vector_ij/dist_ij_2;

  /*--- Mean gradient approximation. Projection of the mean gradient in the direction of the edge ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_GradHeatVar_Normal[iVar] = 0.0; //TODO, Heat module characters
    Proj_Mean_GradHeatVar_Corrected[iVar] = 0.0; //TODO, Heat module characters
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradHeatVar[iVar][iDim] = 0.5*(ConsVar_Grad_i[iVar][iDim] + ConsVar_Grad_j[iVar][iDim]); //TODO, Heat module characters
      Proj_Mean_GradHeatVar_Normal[iVar] += Mean_GradHeatVar[iVar][iDim]*Normal[iDim]; //TODO, Heat module characters
    }
    Proj_Mean_GradHeatVar_Corrected[iVar] = Proj_Mean_GradHeatVar_Normal[iVar]; //TODO, Heat module characters
  }

  val_residual[0] = Thermal_Diffusivity_Mean*Proj_Mean_GradHeatVar_Corrected[0]; //TODO, Heat module characters


  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/
  if (implicit) {
    Jacobian_i[0][0] = -Thermal_Diffusivity_Mean*proj_vector_ij;
    Jacobian_j[0][0] = Thermal_Diffusivity_Mean*proj_vector_ij;
  }

  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();

}


CAvgGradCorrected_Maxwell::CAvgGradCorrected_Maxwell(unsigned short val_nDim, unsigned short val_nVar,
                                                   CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

  implicit        = (config->GetKind_TimeIntScheme_Heat() == EULER_IMPLICIT); //TODO, Heat module characters
  // TODO, Heat module characters, many heat variables
  Edge_Vector = new su2double [nDim];
  Proj_Mean_GradHeatVar_Edge = new su2double [nVar];
  Proj_Mean_GradHeatVar_Kappa = new su2double [nVar];
  Proj_Mean_GradHeatVar_Corrected = new su2double [nVar];
  Mean_GradHeatVar = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_GradHeatVar[iVar] = new su2double [nDim];

}

CAvgGradCorrected_Maxwell::~CAvgGradCorrected_Maxwell(void) {

  delete [] Edge_Vector;
  delete [] Proj_Mean_GradHeatVar_Edge;
  delete [] Proj_Mean_GradHeatVar_Kappa;
  delete [] Proj_Mean_GradHeatVar_Corrected;
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_GradHeatVar[iVar];
  delete [] Mean_GradHeatVar;

}

void CAvgGradCorrected_Maxwell::ComputeResidual(su2double *val_residual, su2double **Jacobian_i, su2double **Jacobian_j, CConfig *config) {

  AD::StartPreacc();
  AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(Temp_i); AD::SetPreaccIn(Temp_j);
  AD::SetPreaccIn(ConsVar_Grad_i[0],nDim); AD::SetPreaccIn(ConsVar_Grad_j[0],nDim);
  AD::SetPreaccIn(Thermal_Diffusivity_i); AD::SetPreaccIn(Thermal_Diffusivity_j);

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
