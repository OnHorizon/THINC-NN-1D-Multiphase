/*
 * hype.h
 *      Author: sunder
 */

#ifndef HYPE_H_
#define HYPE_H_

#include<iostream>
#include<fstream>
#include<string>
#include <cmath>
#include <sstream>
#include <chrono>
#include <stdio.h>
#include <string.h> 
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#define BOOST_DISABLE_ASSERTS
#include "boost/multi_array.hpp"

using namespace boost;


enum recon{THINC, superbee, minmod, van_leer, mc, upwind};

const double small_num = 1.0e-12;

//----------------------------------------------------------------------------
// MLP Class
//----------------------------------------------------------------------------

class MultiLayerPerceptron {
    int n_input_neurons; 
    int n_hidden_layers; 
    int* n_hidden_neurons; 
    int n_output_neurons; 
    int* rows; 
    int* cols; 
    gsl_matrix**W;              // Weight matrix for each layer 
    gsl_vector**b;              // bias vector for each layer
    gsl_vector**z;              // Net sum (W*a^{l-1} + b) 
    gsl_vector**a;              // Activated output (sigma(z))
    inline double sigmoid(double x) {return 1./(1.+std::exp(-x));}  
    inline double reLU(double x) {return std::max(0., x);}
    void softMax(const gsl_vector*, gsl_vector*); 
public:
    MultiLayerPerceptron(); 
    ~MultiLayerPerceptron(); 
    double predict(gsl_vector*); 
}; 

//----------------------------------------------------------------------------
// Six-Equation Model (1D)
//----------------------------------------------------------------------------

// Various types of initial conditions

enum init_cond{smooth_problem,air_helium, air_water, interface_advection, water_cavitation, dodecane_vapor_mixture, epoxy_spinel, multi_component_shu_osher};


class Six_Equation {

    const int N_ph = 3;               // Number of ghost cells at each end
    const double CFL = 0.2;           // CFL Number
    const int nVar = 6;               // Number of variables in the system
    const double eps_vol = 1.0e-8;    // Small number for initializing zero volume fraction 

    multi_array<double,1> x;        // Cell centers
    multi_array<double,1> beta_vf;        // beta values for volume fraction
    multi_array<double,1> eta_vf;        // beta values for volume fraction
    multi_array<double,2> qh;       // Cell average conserved variable
    multi_array<double,2> qhi;       // Initial Condition    
    multi_array<double,2> qh0;      // Cell average conserved variable (for time stepping)
    multi_array<double,2> vh;       // Cell average primitive variable
    multi_array<double,2> vh0;       // IC of Cell average primitive variable
    multi_array<double,2> dqh;      // RHS value for time update
    multi_array<double,3> vbnd;     // Value of primitive variable at cell faces 
    multi_array<double,2> Dm;       // Upwind flux at each face
    multi_array<double,2> Dp;       // Non-Conservative flux at each face

    int IMAX;                   // Number of cells in the domain
    recon reconstruction;       // Reconstruction Method             
    double xmin, xmax;          // Domain boundaries
    double time;                // Simulation time
    double tend;                // Final time of simulation
    double dt;                  // Time step
    int time_step;              // Time step size
    double dx;                  // Mesh size
    double s_max;               // Maximum eigenvalue  (for time step calculation)
    init_cond icond;            // Initial condition
    std::string filename;       // Output file name

    double gamma_1, pi_1, gamma_2, pi_2;
    double GAMMA_1, PI_1, GAMMA_2, PI_2;


    // PDE related functions

    void cons2Prim(const multi_array<double,1>&, multi_array<double,1>&) const;
    void prim2Cons(const multi_array<double,1>&, multi_array<double,1>&) const;
    double speedOfSound(const multi_array<double,1>&);
    void calcStarStateHLLC(const multi_array<double,1>&, const multi_array<double,1>&, double&, double&, double&, multi_array<double,1>&, multi_array<double,1>&, multi_array<double,1>&);
    double calcFluctuationsHLLC(const multi_array<double,1>&, const multi_array<double,1>&, multi_array<double,1>&, multi_array<double,1>&);
    void calcFlucSrcHLLC(const multi_array<double,1>&, const multi_array<double,1>&, multi_array<double,1>&);

    // Main solver 

    void initCond(double,multi_array<double,1>&);  
    void applyBoundaryConditions(); 
    void computeRHS();
    void stepSSPRK22(); 
    void relaxPressure(); 

    // Limiters
    double minmod_lim(double,double) const;
    double vanleer_lim(double,double) const;
    double superbee_lim(double,double) const;
    double mc_lim(double,double) const;
    double upwind_lim(double,double) const;   
    void THINC_lim(double,double,double,double,double&,double&);  
    
    // Initial conditions

    void smoothProblem(double, multi_array<double,1>&);
    void waterCavitation(double, multi_array<double,1>&); 
    void airHelium(double, multi_array<double,1>&);
    void multiComponentShuOsher(double, multi_array<double,1>&);
    void interfaceAdvection(double, multi_array<double,1>&);
    void airWater(double, multi_array<double,1>&);
    void dodecaneVaporMixture(double, multi_array<double,1>&);
    void epoxySpinel(double, multi_array<double,1>&);
public:

    Six_Equation(int, recon, init_cond);
    ~Six_Equation();
    void run(); 
    void plot() const; 
    void beta_plot();
    // Errors
    void errors() const; 

};

#endif /* HYPE_H_ */
