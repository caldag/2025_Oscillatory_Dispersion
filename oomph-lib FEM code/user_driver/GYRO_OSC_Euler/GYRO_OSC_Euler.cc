/*
//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented, 
//LIC// multi-physics finite-element library, available 
//LIC// at http://www.oomph-lib.org.
//LIC// 
//LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
//LIC// 
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC// 
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC// 
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC// 
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC// 
//LIC//====================================================================

// Part of the supplementary material for Caldag & Bees, 
Fine-tuning the dispersion of active suspensions with oscillatory flows.

// Driver code for an adaptive 2D advection diffusion problem with flux boundary 
// conditions using two separate meshes for the bulk and surface meshes.*/

// Code mostly adapted from the code for:

// Bearon, R. N., Bees, M. A., & Croze, O. A. (2012). 
// Biased swimming cells do not disperse in pipes as tracers: 
// a population model based on microscale behaviour. Physics of fluids, 24(12).
///////////////////////////////////////////////////////////////////////////////
// The code runs with oomph-lib, an open-source finite-element solver.
// To learn how to run this driver code, visit:

// https://oomph-lib.github.io/oomph-lib/doc/the_distribution/html/index.html

// This code and Makefile.am in the same directory are the files needed to
// be placed in a folder in user_drivers in the oomph-lib distribution.

// The code utilizes custom-built libraries. These need to be installed before
// running the code. For instructions, see the link above. The library files are
// placed in user_src subdirectory.

// The custom-built library is built upon the generalised advection-diffusion
// library available in the oomph-lib distribution.

// The name of the code is kept joes_poisson_code.cc because the example code in
// the oomph-lib website uses it. One needs to change other parts of the Makefile.am
// to change the file name. Encountered errors while doing it, so keeping it like this for now.

//Generic routines
#include "generic.h"
#include <math.h>
#include <cmath>
#include <sys/stat.h>

#include <cstring> 
#include <string> 
#include <chrono>

#include "adv_diff_time_wind.h"

// The mesh
#include "meshes/rectangular_quadmesh.h"
#include <complex>

using namespace std; // for outputs
using namespace std::chrono; // For timers
using namespace oomph;

////===========================ParameterSpace=========================//
/// Namespace for the problem, includes main physical and geometric parameters,
/// flow, swimming and diffusivity functions.
//====================================================================
namespace ParameterSpace
{
 /// Peclet number
 double Peclet=2.0;    /// Peclet number, determines flow amplitude (based on root-mean-square velocity)
 double Sc=16.8;       /// Schmidt number, enters into the problem through Peclet-Strouhal number
 double Beta=1.0;  	   /// Swimming Peclet number (ie non-dimensional swimming speed)
 double flow_sign=1.0; /// Flow direction sign, changes both the flow and shear profiles
					   /// Upwards (positive flow_sign) flow is in the direction of gravity
					   /// (during the first half of the period)
 extern double Wo=5.0; /// Womersley number default, changed during runtime.
 extern double Peclet_St = Wo*Wo*Sc;  /// Peclet-Strouhal number, equivalent to Wo^2*Sc

extern double velsc=1.0;
/// velsc ensures the root-mean-square velocity is equal to the Pe defined above.
/// The value has to be set for each different Wo manually since there is not a good curve
/// fit for all Wo. velsc coefficients for the majority of Wo tested in the article are available under
/// the function setWo().

double Length=700.0;   /// Channel length (periodic conditions apply), change depending on simulation duration and Wo
 

//=========================== wind_function =========================//
// The wind function applies the fluid flow to the system.
// Different from the conventional wind function used in the advection-diffusion equations in oomph-lib,
// here the wind function is a function of time as we have a flow oscillating in time.
//====================================================================
void wind_function(const Vector<double>& x, double& tval, Vector<double>& wind)
 {
  const   complex<double> I(0.0,1.0);
  wind[0]=0.0; /// No flow in the cross-stream direction (x in simulation frame)

  /// In the simulation model, the cross-stream coordinates vary from 0 to 2 instead of
  /// -1 to 1, so we subtract 1 from the x[0] in these expressions
  if (Wo!=0.0) /// Oscillating flow if Wo is non-zero
  {
    wind[1]=flow_sign*Peclet*velsc*real(1.0/(I*(Wo*Wo))*(cosh(sqrt(I)*Wo*(x[0]-1))/(cosh(sqrt(I)*Wo))-1.0)*exp(I*tval));
  }
  else /// If Wo is zero we will have simple Poiseuille field
  {
    wind[1]=flow_sign*Peclet*(1-(x[0]-1)*(x[0]-1));
  }
 }

//=========================== swimming =========================//
// swimming function represents the swimming component of concentration field.
// It acts as a conservative wind function.
// Different from the conventional wind function used in the advection-diffusion equations in oomph-lib,
// here the wind function is a function of time as we have a flow oscillating in time.
//==============================================================
void swimming(const Vector<double> &x, double &tval, Vector<double> &swim)
{
  double sigma; /// the shear

   if (Wo!=0.0) /// shear profile for non-zero Womersley number
   {
    const   complex<double> I(0.0,1.0);  /// complex number i
    const   complex<double> c2(2.0,0.0); /// This value is non-complex, 2.0, but was required to be written in complex form

    sigma=-0.5*flow_sign*Peclet*velsc*real( 1.0*(pow(c2, 0.5)*sinh(pow(c2, 0.5)*Wo*(x[0]-1)*(0.5 + I/c2))*(cos(tval)+sin(tval)*I)*(0.5 - I/c2))/(Wo*cosh(pow(c2, 0.5)*Wo*(0.5 + I/c2))))/(Beta*Beta*(Wo*Wo*Sc));
   }
  else
  {
    sigma=0.5*flow_sign*(1.0-x[0])*2.0*Peclet; /// Non-oscillating case
  }

  /// Components of swimming expressed as functions of sigma, the fluid shear.
  /// Curves are fit to obtain coefficients a0, a2, c2 and c4. See Appendix (b). 
  
  // 2D swimming, y-direction (lateral), Chlamydomonas augustae (our case)
    double a0=3.814e-001;
    double a2= 4.489e-002;
    double c2= 2.959e-001;
    double c4= 4.405e-002;

  swim[0] = -Beta*sigma*(a0+a2*sigma*sigma)/(1+c2*sigma*sigma+c4*sigma*sigma*sigma*sigma);
  /// Lateral swimming component

  // 2D swimming, x- direction (vertical), Chlamydomonas augustae (our case)
  a0=7.281e-001;
  a2=7.52e-002;
  c2=3.584e-001;
  c4=6.06e-002;

  swim[1] = -Beta*(a0+a2*sigma*sigma)/(1+c2*sigma*sigma+c4*sigma*sigma*sigma*sigma);
  /// Vertical swimming component
 }
 
//=========================== source_function =========================//
// The system has no source term so it is set to zero. The function needs
// to exist in the model because the problem constructor expects it as an input.
//=====================================================================
  void source_function(const Vector<double>& x_vect, double& source)
 {
  source=0.0;
 }
 
//=========================== diff_function =========================//
// diff_function represents the diffusive component of concentration field.
// Different from the conventional diffusion function used in the advection-diffusion equations in oomph-lib,
// here the diffusion function is a function of time as we have a flow oscillating in time.
//===================================================================
 void diff_function(const Vector<double> &x, double& tval, DenseMatrix<double> &D)
{
  double sigma; /// the shear

   if (Wo!=0.0) /// shear profile for non-zero Womersley number
   {
    const   complex<double> I(0.0,1.0);  /// complex number i
    const   complex<double> c2(2.0,0.0); /// This value is non-complex, 2.0, but was required to be written in complex form

    sigma=-0.5*flow_sign*Peclet*velsc*real( 1.0*(pow(c2, 0.5)*sinh(pow(c2, 0.5)*Wo*(x[0]-1)*(0.5 + I/c2))*(cos(tval)+sin(tval)*I)*(0.5 - I/c2))/(Wo*cosh(pow(c2, 0.5)*Wo*(0.5 + I/c2))))/(Beta*Beta*(Wo*Wo*Sc));
   }
  else
  {
    sigma=0.5*flow_sign*(1.0-x[0])*2.0*Peclet; /// Non-oscillating case
  }

  /// Components of diffusion expressed as functions of sigma, the fluid shear.
  /// Curves are fit to obtain coefficients a0, a2, c2 and c4. See Appendix (b). 

    // 2D swimming, D^{yy}, Chlamydomonas augustae
	double a0= 1.7672e-001;
	double a2= 1.269;
	double c2= 7.959;
	double c4= 1.483; 

	D(0,0) = (a0+a2*sigma*sigma)/(1+c2*sigma*sigma+c4*sigma*sigma*sigma*sigma);
	/// Lateral diffusion term
  
    // 2D swimming, cross-diffusion (D^{xy}), Chlamydomonas augustae
	a0= 1.694e-001;
	a2= 5.182e-002;
	c2= 2.166e-001;
	c4= 9.147e-002; 

  D(0,1) =  -sigma*(a0+a2*sigma*sigma)/(1+c2*sigma*sigma+c4*sigma*sigma*sigma*sigma);
  /// Lateral-axial diffusion term

  //D^{xy}=D^{yx}
  D(1,0) =  D(0,1);
  
  // 2D swimming, axial diffusion (D^{xx}), Chlamydomonas augustae
  a0=6.9693e-002;
  a2= 3.046e-001; 
  c2= 1.627e-001;
  c4= 8.72e-002 ;

  D(1,1) =  (a0+a2*sigma*sigma)/(1+c2*sigma*sigma+c4*sigma*sigma*sigma*sigma);
  /// Axial diffusion term
 } 

 void setWo(double Woval); // Function that sets the Womersley and Pe-St numbers

} //end of namespace

namespace ParameterSpace
{
//=========================== setWo =========================//
// setWo is used to access and modify the Womersley number and Peclet-Strouhal number
// from main(). setWo also adjusts the velocity scaling parameter, velsc, as it
// is a function of Womersley number. Values for Wo values tested in the paper are
// provided below, for other Wo velsc has to be computed and entered below.
//===========================================================
  void setWo(double Woval)
  {
    Wo=Woval;
    Peclet_St = Wo*Wo*Sc;

	// This part of the code may be updated depending on the Wo simulated.
	// You can use velsc_calculator.m to compute the velsc coefficient for a given Wo.
	// velsc for several values are provided below.

    if (Wo==0.2)
    {
      velsc=3.8735;
    }
	
	if (Wo==0.3)
    {
      velsc=3.8756;
    }

    if (Wo==0.4)
    {
      velsc=3.8811;
    }

    if (Wo==0.5)
    {
      velsc=3.8928;
    }
	
	if (Wo==0.6)
    {
      velsc=3.9139;
    }
	
	if (Wo==0.7)
    {
      velsc=3.9485;
    }
	
	if (Wo==0.8)
    {
      velsc=4.0010;
    }
	
	if (Wo==1.0)
    {
      velsc=4.1785;
    }
	
	if (Wo==1.5)
    {
      velsc=5.2384;
    }
	
	if (Wo==0.0)
    {
      velsc=1.0;
    }
  }
}

//=========================== start_of_problem_class=========================//
// 2D AdvectionDiffusion problem on rectangular domain, discretised with
// 2D QAdvectionDiffusion elements. The specific type of 
// element is specified via the template parameter.
//===========================================================================
template<class ELEMENT> 
class BiasedActiveMatterDispersionProblem : public Problem
{

public:

 /// Constructor: Pass pointer to source, wind, conservative wind and diffusion functions
 BiasedActiveMatterDispersionProblem(RefineableGeneralisedAdvectionDiffusionEquationsTimeWind<2>::GeneralisedAdvectionDiffusionSourceFctPt source_fct_pt,
  RefineableGeneralisedAdvectionDiffusionEquationsTimeWind<2>::GeneralisedAdvectionDiffusionWindFctPt wind_fct_pt,
  RefineableGeneralisedAdvectionDiffusionEquationsTimeWind<2>::GeneralisedAdvectionDiffusionWindFctPt conserved_wind_fct_pt,
  RefineableGeneralisedAdvectionDiffusionEquationsTimeWind<2>::GeneralisedAdvectionDiffusionDiffFctPt diff_fct_pt);

 /// Destructor (empty)
 ~BiasedActiveMatterDispersionProblem(){}

 /// Doc the solution. DocInfo object stores flags/labels for where the
 /// output gets written to
 void doc_solution(DocInfo& doc_info, double curtime);
 
 void set_initial_condition(); /// sets the initial concentration distribution
 
double global_temporal_error_norm(); /// computes the error

double Target_error_safety_factor = 0.5; /// The factor allows some relaxation in time stepping

private:

 /// No action before solving (empty)
 void actions_before_newton_solve();

 /// Update the problem specs after solve
 void actions_after_newton_solve(){};

 /// Actions before adapt (empty)
 void actions_before_adapt();

 /// Actions after adapt: Rebuild the mesh
 void actions_after_adapt();
 
 /// Actions before implicit timestep (redefining the wind function as it is a function of time)
 void actions_before_implicit_timestep(){};
 
 /// Pointer to source function (set to zero, no effect)
 RefineableGeneralisedAdvectionDiffusionEquationsTimeWind<2>::GeneralisedAdvectionDiffusionSourceFctPt Source_fct_pt;

 /// Pointer to the "bulk" mesh
 //RefineableRectangularQuadMesh<ELEMENT>* Bulk_mesh_pt;
 
 RefineableRectangularQuadMesh<ELEMENT>* mesh_pt() /// Mesh pointer
{
 return dynamic_cast<RefineableRectangularQuadMesh<ELEMENT>*>(
  Problem::mesh_pt());
}
 
 /// Pointer to wind function
 RefineableGeneralisedAdvectionDiffusionEquationsTimeWind<2>::GeneralisedAdvectionDiffusionWindFctPt Wind_fct_pt;
 
 /// Pointer to swimming function
 RefineableGeneralisedAdvectionDiffusionEquationsTimeWind<2>::GeneralisedAdvectionDiffusionWindFctPt Conserved_wind_fct_pt;
 
 /// Pointer to diffusivity
 RefineableGeneralisedAdvectionDiffusionEquationsTimeWind<2>::GeneralisedAdvectionDiffusionDiffFctPt Diff_fct_pt;

 void possibly_disable_ALE() // Disabling arbitrary Lagrangian-Eulerian framework for speed
  {
   bool Use_ALE=0;
   // Loop over the elements 
   unsigned n_element = mesh_pt()->nelement();
   for(unsigned i=0;i<n_element;i++)
    {
     // Upcast from FiniteElement to the present element
     ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i));
     
     if (Use_ALE)
      {
       el_pt->enable_ALE();
      }
     else
      {
       el_pt->disable_ALE();
      }
    }
  }
}; // end of problem class


//===========================start_of_constructor=========================//
// Constructor for the problem: Pass pointer to source, wind and
// diffusion functions. 
//========================================================================
template<class ELEMENT>
BiasedActiveMatterDispersionProblem<ELEMENT>::
BiasedActiveMatterDispersionProblem(
 RefineableGeneralisedAdvectionDiffusionEquationsTimeWind<2>::GeneralisedAdvectionDiffusionSourceFctPt source_fct_pt,
 RefineableGeneralisedAdvectionDiffusionEquationsTimeWind<2>::GeneralisedAdvectionDiffusionWindFctPt wind_fct_pt, 
 RefineableGeneralisedAdvectionDiffusionEquationsTimeWind<2>::GeneralisedAdvectionDiffusionWindFctPt conserved_wind_fct_pt,
 RefineableGeneralisedAdvectionDiffusionEquationsTimeWind<2>::GeneralisedAdvectionDiffusionDiffFctPt diff_fct_pt)
 : Source_fct_pt(source_fct_pt), Wind_fct_pt(wind_fct_pt), Conserved_wind_fct_pt(conserved_wind_fct_pt), Diff_fct_pt(diff_fct_pt)
{ 
 /// Add the time stepper
 this->add_time_stepper_pt(new BDF<2>(true));
 
 // Setup "bulk" mesh

 /// # of elements in x-direction (cross-stream)
 unsigned n_x=20;

 /// # of elements in the y-direction (axial)
 unsigned n_y=ParameterSpace::Length;
 
 /// Domain length in x-direction
 /// Setting this l_x=2.0 so that x can vary between -1 and 1
 /// when 1 is subtracted from the value
 double l_x=2.0;

 /// Domain length in the y-direction
 double l_y=ParameterSpace::Length;

 /// Build "bulk" mesh
 Problem::mesh_pt()=new RefineableRectangularQuadMesh<ELEMENT>(n_x,n_y,l_x,l_y,this->time_stepper_pt());

 /// Create/set error estimator
 mesh_pt()->spatial_error_estimator_pt()=new Z2ErrorEstimator;
 
 mesh_pt()->max_permitted_error() = 1.0e-3;
 mesh_pt()->min_permitted_error() = 1.0e-5;
 mesh_pt()->max_refinement_level() = 10; /// Error preferences

  // Setting top and bottom boundaries as periodic.
  // Taken from oomph-lib documentation.

  unsigned nx=n_x;
  unsigned ny=n_y;
 
  unsigned n_node = mesh_pt()->nboundary_node(0);
  for(unsigned n=0;n<n_node;n++)
  {
   mesh_pt()->boundary_node_pt(0,n)
    ->make_periodic(mesh_pt()->boundary_node_pt(2,n));
  }
  
  /// Get pointers to tree roots associated with elements on the 
  /// left and right boundaries
  Vector<TreeRoot*> left_root_pt(ny);
  Vector<TreeRoot*> right_root_pt(ny);
  for(unsigned i=0;i<ny;i++) 
   {
    left_root_pt[i] = 
     dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(i*nx))->
     tree_pt()->root_pt();
    right_root_pt[i] = 
     dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(nx-1+i*nx))->
     tree_pt()->root_pt();
   }
   
  /// Switch on QuadTreeNames for enumeration of directions
   using namespace QuadTreeNames;
 
  /// Set the neighbour and periodicity
  for(unsigned i=0;i<ny;i++) 
   {
    /// The western neighbours of the elements on the left
    /// boundary are those on the right
    left_root_pt[i]->neighbour_pt(W) = right_root_pt[i];
    left_root_pt[i]->set_neighbour_periodic(W); 
    
    /// The eastern neighbours of the elements on the right
    /// boundary are those on the left
    right_root_pt[i]->neighbour_pt(E) = left_root_pt[i];
    right_root_pt[i]->set_neighbour_periodic(E);     
   } /// done

   // End of assigning periodic boundaries

  // Building the mesh
  //add_sub_mesh(mesh_pt());
  //build_global_mesh();

 // Loop over the bulk elements to set up element-specific 
 // things that cannot be handled by constructor: Pass pointer to 
 // wind and diffusion function
  unsigned n_element = mesh_pt()->nelement();
  for(unsigned e=0;e<n_element;e++)
  {
    /// Upcast from GeneralisedElement to current bulk element
    ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));

    /// Set the source function pointer
    el_pt->source_fct_pt() = Source_fct_pt;

    /// Set the wind function pointer
    el_pt->wind_fct_pt() = &ParameterSpace::wind_function;
   
    /// Swimming component
    el_pt->conserved_wind_fct_pt() =  &ParameterSpace::swimming;
   
    /// Diffusivity
    el_pt->diff_fct_pt() = &ParameterSpace::diff_function;

    /// Set the Peclet number
    el_pt->pe_pt() = &ParameterSpace::Peclet;
   
    /// Set the Peclet Strouhal number
    el_pt->pe_st_pt() = &ParameterSpace::Peclet_St;

    /// Set the Peclet number
    el_pt->pe_pt() = &ParameterSpace::Peclet;

	/// Time pointer
    el_pt->time_fct_pt()= time_pt();
  }
  // Set up equation numbering scheme
  cout <<"Number of equations: " << assign_eqn_numbers() << std::endl;
  possibly_disable_ALE();  
} // end of constructor

//===========================Temporal error norm function=========================//
// This function computes the temporal error to decide if the time step taken is
// appropriate. The error should be less than the parameter epsilon defined in main().
//================================================================================

template<class ELEMENT>
double BiasedActiveMatterDispersionProblem<ELEMENT>::global_temporal_error_norm()
{
 double global_error = 0.0;
   
 /// Find out how many nodes there are in the problem
 unsigned n_node = mesh_pt()->nnode();
 
 /// Loop over the nodes and calculate the estimated error in the values
 for(unsigned i=0;i<n_node;i++)
  {
   /// Get error in solution: Difference between predicted and actual
   /// value for nodal value 0
   double error = mesh_pt()->node_pt(i)->time_stepper_pt()->
    temporal_error_in_value(mesh_pt()->node_pt(i),0);
   
   /// Add the square of the individual error to the global error
   global_error += error*error;
  }
    
 /// Divide by the number of nodes
 global_error /= double(n_node);
 
 /// Return square root...
 return sqrt(global_error);
 
} // end of global_temporal_error_norm

//===========================start_of_actions_before_newton_solve=========================//
/// Update the problem specs before solve (empty)
//========================================================================================
template<class ELEMENT>
void BiasedActiveMatterDispersionProblem<ELEMENT>::actions_before_newton_solve()
{
} // end of actions before solve

//===========================set_initial_condition=========================//
// set_initial_condition initializes the concentration field.
// Any function of position can be provided.
// For numerical stability, the concentration field should not have sharp gradients.
//========================================================================
 template<class ELEMENT>
void BiasedActiveMatterDispersionProblem<ELEMENT>::
set_initial_condition()
{
 unsigned n_node = mesh_pt()->nnode(); /// Will go through nodes to assign the concentration to each
 for(unsigned n=0;n<n_node;n++) 	   /// Sweeping through all nodes
  {
    Node* nod_pt = mesh_pt()->node_pt(n); /// Take the current node
    double yv = nod_pt->x(0); 			  /// The cross-stream coordinate
    double xv = nod_pt->x(1); 			  /// The axial coordinate

    double u = 1.0*exp(-0.01*abs(xv-400.0)*abs(xv-400.0)); /// Concentrated axially at a certain position (400 here)
	/// The decay coefficient (-0.01) should be kept small (in amplitude) to prevent instabilities
	
    nod_pt->set_value(0,u); /// Set the concentration value for the current node
  }
} 

//===========================start_of_doc=========================//
/// Doc the solution: doc_info contains labels/output directory etc.
//================================================================
template<class ELEMENT>
void BiasedActiveMatterDispersionProblem<ELEMENT>::
doc_solution(DocInfo& doc_info, double curtime)
{ 
  ofstream some_file;
  char filename[100];

  /// Number of plot points
  unsigned npts;
  npts=5; 

  /// Output solution 
  sprintf(filename,"%s/soln%i_t=%lf.dat",doc_info.directory().c_str(),
         doc_info.number(),curtime);
  some_file.open(filename);
  mesh_pt()->output(some_file,npts);
  some_file.close();
} // end of doc
 

//===========================start_of_actions_before_adapt=========================//
// Remove the flux elements from the mesh (empty).
//=================================================================================
template<class ELEMENT>
void BiasedActiveMatterDispersionProblem<ELEMENT>::actions_before_adapt()
{
} // end of actions_before_adapt

//===========================start_of_actions_after_adapt=========================//
// Attach flux elements to the mesh.
//================================================================================
template<class ELEMENT>
void BiasedActiveMatterDispersionProblem<ELEMENT>::actions_after_adapt()
{
  unsigned n_element = mesh_pt()->nelement(); /// Assignments need to be redone if mesh is adapted

  for(unsigned e=0;e<n_element;e++)
  {
    // Upcast from GeneralisedElement to current bulk element
    ELEMENT *el_pt = dynamic_cast<ELEMENT*>(mesh_pt()->element_pt(e));
    
    //Set the source function pointer
    el_pt->source_fct_pt() = Source_fct_pt;

    //Set the wind function pointer
    el_pt->wind_fct_pt() = &ParameterSpace::wind_function;
   
    // Swimming component
    el_pt->conserved_wind_fct_pt() =  &ParameterSpace::swimming;
   
    // Diffusivity
    el_pt->diff_fct_pt() = &ParameterSpace::diff_function;

    // Set the Peclet number
    el_pt->pe_pt() = &ParameterSpace::Peclet;
   
    //Set the Peclet Strouhal number
    el_pt->pe_st_pt() = &ParameterSpace::Peclet_St;

    // Set the Peclet number
    el_pt->pe_pt() = &ParameterSpace::Peclet;

    el_pt->time_fct_pt()= time_pt();
    possibly_disable_ALE();
  }
} // end of actions_after_adapt

//===========================start_of_main=========================//
// Defines the problem. Parameter sweeps and export related settings adjusted.
// Time stepping carried out in a for loop.
// The outermost for loop sweeps through different Wo.
//=================================================================
int main(int argc, char* argv[])
{
  CommandLineArgs::setup(argc,argv);

  double Wos[7]={0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0}; /// Running for multiple Wo
  double Sc=16.8; /// The Schmidt number

 /// Set up the problem with 2D nine-node elements from the
 /// RefineableQGeneralisedAdvectionDiffusionElementTimeWind family.
 /// This is quite similar to RefineableQGeneralisedAdvectionDiffusionElement
 /// except for wind and diffusivity functions also getting current time as inputs.
  BiasedActiveMatterDispersionProblem<RefineableQGeneralisedAdvectionDiffusionElementTimeWind<2,
    3> > problem(&ParameterSpace::source_function,
          &ParameterSpace::wind_function,
          &ParameterSpace::swimming,
          &ParameterSpace::diff_function);
		  
  problem.target_error_safety_factor()=0.5; /// This allows some relaxation in time stepping
  
  for (int ii=0; ii<7; ii++) /// Sweeping through Womersley numbers
  {
    ParameterSpace::setWo(Wos[ii]); /// First set the Womersley number	
	
	double timeabs=300.0*Wos[ii]*Wos[ii]*Sc;/// The absolute minimum of simulation duration
  /// in terms of oscillation timescale. The value 300 is in terms of diffusive timescale, this is
  /// multiplied with Wo*Wo*Sc to obtain the corresponding value non-dimensionalized with the 
  /// oscillation time scale. Exact duration varies from Wo to Wo 
  /// because the simulation always finishes at the end of a complete period.
  
    double period=2*M_PI; /// Period of oscillation
    double ncycles=ceil(timeabs/period); /// Derive a number of cycles such that the simulation duration is larger than timeabs
    double t_max = ncycles*period; /// The actual simulation duration t_max is an integer number of periods.

    /// Create label for output
    DocInfo doc_info;
    /// Set output directory, numbering according to loop index
    int folder_no=ii;
    string out_dir="RESLT"+std::to_string(folder_no);
    int folder_length = out_dir.length();
    char* char_folder = new char[folder_length + 1];
    strcpy(char_folder, out_dir.c_str());
    mkdir(char_folder,0777);
    doc_info.set_directory(out_dir);
    /// Step number
    doc_info.number()=0;

    bool PAUSEREC=false;

  /// Self-testing before attempting the solution
  /// This is a built-in routine.
    cout << "\n\n\nProblem self-test ";
    if (problem.self_test()==0) 
    {
    cout << "passed: Problem can be solved." << std::endl;
    }
    else 
    {
    throw OomphLibError("Self test failed",
                        OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
    }

    double dt = period/200; /// Time step (modified by adaptive routines)
    double t_rec=period/8; /// Time interval for the results to be recorded
	/// We record the first and last 10 periods of oscillation only to save space.

    std::cout << "Period is  " << period << " seconds.\n";
    std::cout << "Record interval is  " << t_rec << " seconds.\n"; /// General info

    problem.set_initial_condition(); 			 /// Set initial condition
	problem.assign_initial_values_impulsive(dt); /// This ensures all the history values before t=0 are equal to the initial values.
    problem.doc_solution(doc_info,0.0); 		 /// Output solution (the initial condition)
    doc_info.number()++; 						 /// Increase the count to save with a new file name
    double epsilon_t=1.0e-2; 					 /// Error threshold between time steps

	bool first = true; 						     /// Checks whether the current step was the first time step of the problem

    unsigned istep=0; 							 /// Counter for time steps
    double last_time=0.0; 						 /// Keeps when a solution was last recorded
    problem.time_pt()->time()=0.0; 				 /// Set the timer to zero
    
    auto problemstart = high_resolution_clock::now(); /// used for timing the code, the start point of solving a case

  while (problem.time_pt()->time()<t_max) 		 /// Going through the time steps
    {
      auto start = high_resolution_clock::now(); /// Timing related, used to evaluate how long it takes to simulate a time step
      istep=istep+1; 							 /// Increase the step counter
      double dt_next=problem.adaptive_unsteady_newton_solve(dt,epsilon_t); /// Adaptive solver, reduces dt if error criteria not met
      
      auto problemstep=high_resolution_clock::now(); /// Timing related, used to evaluate how long it takes to simulate a time step
      auto simtime=duration_cast<minutes>(problemstep - problemstart);
    
      std::cout << "A time step of " << dt << "has been taken.\n"; /// Output to the terminal
      std::cout << "Simulation has been running for " << simtime.count() << "minutes. \n";
      std::cout << "The duration simulated so far: " << problem.time_pt()->time() << " out of " << t_max << ". \n";

      if(problem.time_pt()->time()==0.0) /// If this is the first time step
      {
        first=false;							 /// Future steps won't be the first any more.
        problem.doc_solution(doc_info,0.0); 	 /// Output solution 
        doc_info.number()++; 					 /// Increment counter for solutions
      }
  
      if (istep>0) 								 /// If number of steps are larger than 0
      {
        std::cout << "Current time: " << problem.time_pt()->time() << "\n"; /// Report simulation time
        if (problem.time_pt()->time() == last_time+t_rec && !PAUSEREC) 		/// If the simulation has advanced in time pi/4 amount
        {   
          std::cout << "----------------------------------------------------\n";
          std::cout << "----------------------------------------------------\n";
          std::cout << "----------------------------------------------------\n";
          std::cout << "Solving for simulation no " << ii << "\n";
          std::cout << "----------------------------------------------------\n";
          std::cout << "----------------------------------------------------\n";
          std::cout << "----------------------------------------------------\n";
		  
          last_time=problem.time_pt()->time(); 	 /// Set the last_time a solution was recorded to the current time
          std::cout << "Recording time: " << problem.time_pt()->time() << "\n";
          std::cout << "----------------------------------------------------\n";
          std::cout << "----------------------------------------------------\n";
          std::cout << "----------------------------------------------------\n";

          problem.doc_solution(doc_info,last_time); /// Record the solution 
          doc_info.number()++;						/// Increment the count used in file naming
          dt=dt_next;								/// Set the next time step

          if (dt>t_rec)								/// If the time step is larger than the interval of recording
          {
            dt=t_rec;								/// Then limit the time step with the recording interval
            cout<<"Upper limit on dt! Will retain t_rec...\n";
          }
        }
        else 										/// Determining the time stepping if data is not recorded
        {
          if (problem.time_pt()->time()+dt_next > last_time+t_rec && !PAUSEREC)
          /// If the code will pass the instance that I need to record with the current time step,
          /// set the time step such that I land on that next instance to record
          {
            dt=(last_time+t_rec-problem.time_pt()->time());

            if (dt<0) /// This is some error handling for negative dt, the code is never supposed to go in here
            {
              cout << "Negative dt found, may miss a recording point in time!!!\n";
              cout << "Using the suggested time step instead to prevent simulation crashing...\n";
              dt=dt_next;
            }

          }
          else /// Otherwise use whatever derived by the software
          {
            dt=dt_next; /// Record the time step from previous stepping, will be the best candidate
          }
        }

		/// We only record data from the first and last 10 periods of simulation to save space.
        if (problem.time_pt()->time()+dt > 10.0*period && problem.time_pt()->time()+dt < t_max-10.0*period)
        {
          PAUSEREC=true; /// Setting PAUSEREC to true prevents recording of solution data
        }
        if (problem.time_pt()->time()+dt >= t_max-10.0*period && PAUSEREC)
        {
          cout << "Approaching re-recording phase, time step should decrease here...\n"; 
		  /// If we are approaching the time to start recording, we should see it before taking the next time step.
          PAUSEREC=false;
          dt = t_max-10.0*period-problem.time_pt()->time();
          last_time=t_max-10.0*period-t_rec;
        }
      }
    } 
  }
} // end of main