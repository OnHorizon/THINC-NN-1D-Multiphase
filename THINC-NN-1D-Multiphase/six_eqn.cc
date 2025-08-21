/*
 * six_eqn.cc
 *      Author: sunder
 */

#include "hype.h"

//----------------------------------------------------------------------------
// 1D Tanh Reconstruction 
//----------------------------------------------------------------------------

double tanh1D(double x, double x0, double eps, double in, double out) {

    return 0.5*((in+out) + (out-in)*std::tanh((x-x0)/(eps)));
}

//----------------------------------------------------------------------------
// Constructor
//----------------------------------------------------------------------------

Six_Equation::Six_Equation(int _IMAX, recon _reconstruction, init_cond _icond) :
    x(extents[_IMAX]),
    beta_vf(extents[_IMAX]),
    eta_vf(extents[_IMAX]),
    qh(extents[_IMAX][nVar]),
    qhi(extents[_IMAX][nVar]),
    qh0(extents[_IMAX][nVar]),
    vh(extents[multi_array_types::extent_range(-N_ph,_IMAX+N_ph)][nVar]),
    vh0(extents[_IMAX][nVar]),
    dqh(extents[multi_array_types::extent_range(-N_ph,_IMAX+N_ph)][nVar]),
    vbnd(extents[multi_array_types::extent_range(-N_ph,_IMAX+N_ph)][2][nVar]), // 2 is for two faces per cell
    Dm(extents[_IMAX+1][nVar]),
    Dp(extents[_IMAX+1][nVar]),
    IMAX(_IMAX), 
    reconstruction(_reconstruction),
    time(0.0),
    icond(_icond)

{
    int i, c; 
    multi_array<double,1> Q0(extents[nVar]); 

    //delta = gsl_vector_alloc(4);

    xmin = -5.0;
    xmax =  5.0;

    if (icond == smooth_problem) {

        xmin = -1.0;
        xmax =  1.0;
        tend = 2.0;
        gamma_1 = 1.4;
        pi_1 = 0.;
        gamma_2 = 1.6;
        pi_2 = 7.5;
        filename = "Smooth-Problem.csv";

    }

    if (icond == air_helium) {
        tend = 1.5;
        gamma_1 = 1.4;
        pi_1 = 0.0;
        gamma_2 = 1.648;
        pi_2 = 0.0;
        filename = "Air-Helium-Shock-Tube.csv";
    }

    if (icond == multi_component_shu_osher) {
        tend = 1.8;
        gamma_1 = 1.4;
        pi_1 = 0.0;
        gamma_2 = 1.9;
        pi_2 = 0.0;
        filename = "Multi-Component-Shu-Osher.csv";
    }

    if (icond == interface_advection) {

        xmin = -5.0;
        xmax =  5.0;

        //double cycles = 2.0;
        
        tend = 0.279;
        gamma_1 = 4.4;
        pi_1 = 6000.0;
        gamma_2 = 1.4;
        pi_2 = 0.0;
        filename = "Interface-Advection.csv";
    }

    if (icond == air_water) {
        tend = 2.4e-2;
        gamma_1 = 4.4;
        pi_1 = 6000.0;
        gamma_2 = 1.4;
        pi_2 = 0.0;
        filename = "Air-Water-Shock-Tube.csv";
    }

    if (icond == water_cavitation) {
        tend = 0.32;
        gamma_1 = 4.4;
        pi_1 = 6000.0;
        gamma_2 = 1.4;
        pi_2 = 0.0;
        filename = "Water-Cavitation-Test.csv";
    }


    if (icond == dodecane_vapor_mixture) {
        tend = 0.35e-1;
        gamma_1 = 2.35;
        pi_1 = 4000.0;
        gamma_2 = 1.025;
        pi_2 = 0.0;
        filename = "Dodecane-Vapor-Mixture-Shock-Tube.csv";
    }

    if (icond == epoxy_spinel) {

        xmin = 0.0;
        xmax = 1.0;
        tend = 80.0e-6;
        gamma_1 = 2.43;
        pi_1 = 5.3e9;
        gamma_2 = 1.62;
        pi_2 = 141.0e9;
        filename = "E-S.csv";
    }
    
    

    dx = (xmax-xmin)/static_cast<double>(IMAX);
    s_max = 0.0; dt = 0.0; time_step = 0;

    GAMMA_1 = 1.0/(gamma_1-1.0); PI_1 = (gamma_1*pi_1)/(gamma_1-1.0);
    GAMMA_2 = 1.0/(gamma_2-1.0); PI_2 = (gamma_2*pi_2)/(gamma_2-1.0);
    
    for (i = 0; i < IMAX; ++i) {
        x[i] = xmin + (static_cast<double>(i)+0.5)*dx; 
        initCond(x[i],Q0); 
        for (c = 0; c < nVar; ++c) {
            qh[i][c] = Q0[c]; 
            qhi[i][c] = Q0[c];        
        }
    }
    

    
}

//----------------------------------------------------------------------------
// Destructor 
//----------------------------------------------------------------------------

Six_Equation::~Six_Equation() {
    //gsl_vector_free(delta);
}

//----------------------------------------------------------------------------
// Apply boundary conditions
//----------------------------------------------------------------------------

void Six_Equation::applyBoundaryConditions() {
    
    int oned_begin, oned_end, ilhs, irhs;

    // ---------------------- Left boundary ----------------------

    oned_begin = 0; oned_end = IMAX-1;

    for (int i = 0; i < N_ph; ++i) {

            ilhs = oned_begin - i - 1;
            irhs = oned_begin + i; // Transmissive 
            //irhs = oned_end - i;  //Perioidic 


            for (int c = 0; c < nVar; ++c) 
                vh[ilhs][c] = vh[irhs][c];
            
    }

    // ---------------------- Right boundary ----------------------

    for (int i = 0; i < N_ph; ++i) {

        ilhs = oned_end + i + 1;
        irhs = oned_end - i; // Transmissive
        //irhs = oned_begin + i; // Periodic 

        for (int c = 0; c < nVar; ++c) 
            vh[ilhs][c] = vh[irhs][c];
     
    }
}

//----------------------------------------------------------------------------
// Minmod limiter
//----------------------------------------------------------------------------

double Six_Equation::minmod_lim(double a, double b) const {

    if (a*b < 0.0)
        return 0.0; 
    else {
        if (std::abs(a) < std::abs(b)) 
            return a;
        else 
            return b;    
    }  
}

//----------------------------------------------------------------------------
// van-Leer limiter
//----------------------------------------------------------------------------

double Six_Equation::vanleer_lim(double a, double b) const {

    if (a*b < small_num) {
        return 0.0;     
    }
    
    else {
        if (a < 0.0)
            return -2.0*a*b/(std::abs(a) + std::abs(b)); 
        else
            return  2.0*a*b/(std::abs(a) + std::abs(b));
    }
}

//----------------------------------------------------------------------------
// super-bee limiter
//----------------------------------------------------------------------------

double Six_Equation::superbee_lim(double a, double b) const {

    if (a*b < small_num) {
        return 0.0;     
    }
    
    else {

        double beta = 2.0;         
        double value = std::max(std::min(std::abs(a), beta*std::abs(b)), std::min(beta*std::abs(a), std::abs(b)));

        if (a < 0.0)
            return -value; 
        else
            return  value;
    }
}

//----------------------------------------------------------------------------
// Monotized-Central limiter
//----------------------------------------------------------------------------

double Six_Equation::mc_lim(double a, double b) const {

    if (a*b < small_num) 
        return 0.0;     
   
    else {

        double beta = 2.0;         
        double value = std::min(0.5*std::abs(a+b), std::min(beta*std::abs(a), beta*std::abs(b)));

        if (a < 0.0)
            return -value; 
        else
            return  value;
    }
}

//----------------------------------------------------------------------------
// Upwind
//----------------------------------------------------------------------------

double Six_Equation::upwind_lim(double a, double b) const {

    return 0.5*(a+b);
}

double sigmoid(double x) {
    return 1./(1. + std::exp(-x));
}


double beta_NN(double q_im1, double q_i, double q_ip1) {
    
    double b0[] = {-6.991438269615173340e-01,-3.014139413833618164e+00,-5.259067416191101074e-01};
    double b1[] = {-6.548673152923583984e+00, -4.569244384765625000e+00};
    double b2[] = {5.281073093414306641e+00};
    
    double W0[3][1] = {
        {7.860679626464843750e+00},
        {3.703494787216186523e+00},
        {2.847473335266113281e+01}
    };
    
    double W1[2][3] = {
        {6.084172725677490234e+00, 4.970069885253906250e+00, 7.743915081024169922e+00},
        {2.945849180221557617e+00, 2.825659275054931641e+00, 8.736223578453063965e-01}
        
    };
    
    double W2[1][2] = {
        {-4.840313434600830078e+00, -6.451254367828369141e+00}
    };
    
    double r;
    
    r = (q_i-q_im1)/(q_ip1-q_i + 1.0e-12);
    
    double input[1] = {4.*std::pow(r,4.)/std::pow((1.0 +  std::pow(r,4)),2 ) };
    
    double y[3],z[2], beta_[1];
    
    for (int i = 0; i < 3; ++i) {
        y[i] = 0.0;
        for (int j = 0; j < 1; ++j) {
            y[i] += (W0[i][j] * input[j]);
            
        }
        y[i] += b0[i];
    }
    
    for (int i = 0; i < 3; ++i) {
        y[i] = sigmoid(y[i]); 
        
        
    }
    
    
    for (int i = 0; i < 2; ++i) {
        z[i] = 0.0;
        for (int j = 0; j < 3; ++j) {
            z[i] = z[i] + (W1[i][j] * y[j]);
            
        }
        z[i] = z[i] + b1[i];
        
    }
    
    for (int i = 0; i < 2; ++i) {
        z[i] = sigmoid(z[i]);
        
    }
     
   for (int i = 0; i < 1; ++i) {
        beta_[i] = 0.0;
        for (int j = 0; j < 2; ++j) {
            beta_[i] += (W2[i][j] * z[j]);
        }
        
        beta_[i] += b2[i];
    }
    
    for (int i = 0; i < 1; ++i)
        beta_[i] = (1.6 - std::log(3.))*sigmoid(beta_[i]) + std::log(3.);
    
    return beta_[0];
}

double beta_ETA(double q_im1, double q_i, double q_ip1) {

    double r;
    
    r = (q_i-q_im1)/(q_ip1-q_i + 1.0e-12);
    
    double eta = 4.*std::pow(r,4.)/std::pow(1.0 + std::pow(r,4.), 2.);
    double beta_s, beta_d;
    beta_s = std::log(3.);
    beta_d = 1.6;
    
    return eta*beta_s + (1.-eta)*beta_d;
}

//----------------------------------------------------------------------------
// THINC Reconstruction 
//----------------------------------------------------------------------------

void Six_Equation::THINC_lim(double q_im1, double q_i, double q_ip1, double beta, double& L, double& R) {

    if (  (q_ip1-q_i)*(q_i-q_im1) < 0.0 ) {
        L = q_i;
        R = q_i;     
    }

    else {
        
        //double beta = 1.6; 

        double min, max, gamma; 
        
        min = std::min(q_im1,q_ip1);
        max = std::max(q_im1,q_ip1); 
        
        if (q_im1 < q_ip1)
            gamma = 1.0; 
        else
            gamma = -1.0; 

        double C = (q_i - min + 1.0e-20)/(max-min+1.0e-20); 
        double B = std::exp(gamma*beta*(2.*C-1.)); 
        double A = (B/std::cosh(beta)  - 1.)/std::tanh(beta);
        double D = (std::tanh(beta) + A)/(1. + A*std::tanh(beta));


        L = min + 0.5*(max-min)*(1. + gamma*A);
        R = min + 0.5*(max-min)*(1. + gamma*D); 
    
    }
}

//----------------------------------------------------------------------------
// Calculate delta layer for the given stencil
//----------------------------------------------------------------------------
/*
void Six_Equation::calcDelta(double q_im1, double q_i, double q_ip1) {
    
    double d1 = std::abs(q_i-q_im1);
    double d2 = std::abs(q_ip1-q_i);
    double d3 = 0.5*std::abs(q_ip1-q_im1);
    double d4 = std::abs(q_ip1 - 2.*q_i + q_im1);

    double M = std::max(d1,d2) + 1.0e-15;
     
    gsl_vector_set(delta,0,d1/M);
    gsl_vector_set(delta,1,d2/M);
    gsl_vector_set(delta,2,d3/M);
    gsl_vector_set(delta,3,d4/M);
}
*/



void Six_Equation::beta_plot() {
    double r;
    
    for (int i = 0; i < IMAX; ++i){
    
        beta_vf[i] = beta_NN(vh[i-1][4], vh[i][4], vh[i+1][4]);
        r = (vh[i][4]-vh[i-1][4])/(vh[i+1][4]-vh[i][4] + 1.0e-12);
    
        eta_vf[i] = {4.*std::pow(r,4.)/std::pow((1.0 +  std::pow(r,4)),2 ) };
          
    
    }

    std::ofstream out_data;
    out_data.open ("Beta-plot.csv");
    out_data.flags( std::ios::dec | std::ios::scientific );
    out_data.precision(16);

    out_data << "x,beta_alpha_1\n";

    for (int i = 0; i < IMAX; ++i)
        out_data << x[i] << "," << beta_vf[i] << ","<<eta_vf[i]<<std::endl;

    out_data.close();

    

}







double beta_curve(double vim1, double vi, double vip1){
    double d1, d2, M, delta, beta_c;
    d1 = std::abs(vi-vim1);
    d2 = std::abs(vip1-vi);
    M  = std::max(d1,d2) + 1.0e-15;
    delta = std::min(d1,d2)/M;
    beta_c =  1.6 + (1.11-1.6)/(1. + std::exp(-30.*(delta - 0.7)));
    return beta_c;

}










//----------------------------------------------------------------------------
// Compute the RHS in each cell
//----------------------------------------------------------------------------

void Six_Equation::computeRHS() {

    int i, c; 
    const double r1_dx = 1./dx;
    double s, slope, beta; 
    multi_array<double,1> Q(extents[nVar]), V(extents[nVar]), V_x(extents[nVar]), BgradQ(extents[nVar]);
    multi_array<double,1> Flux(extents[nVar]), Fluc(extents[nVar]);
    multi_array<double,1> VL(extents[nVar]), VR(extents[nVar]);

    s_max = 0.0; 

    // Convert to primitive variables 

    for (i = 0; i < IMAX; ++i) {
        for (c = 0; c < nVar; ++c)
            Q[c] = qh[i][c]; 
        cons2Prim(Q,V);
        for (c = 0; c < nVar; ++c)
            vh[i][c] = V[c]; 
    }

    applyBoundaryConditions();

    for (i = -1; i < IMAX+1; ++i) {

        for (c = 0; c < nVar; ++c) {

            if (reconstruction == THINC) {


                beta = beta_NN(vh[i-1][c],vh[i][c],vh[i+1][c]); //beta_curve or beta_ETA or beta_NN
                
                
                //beta = 1.5;
                
                THINC_lim(vh[i-1][c],vh[i][c],vh[i+1][c],beta,vbnd[i][0][c],vbnd[i][1][c]);
                slope = vbnd[i][1][c] - vbnd[i][0][c];
            }

            else {

                if (reconstruction == mc)
                    slope = mc_lim(vh[i][c]-vh[i-1][c], vh[i+1][c]-vh[i][c]);
                else if (reconstruction == van_leer)
                    slope = vanleer_lim(vh[i][c]-vh[i-1][c], vh[i+1][c]-vh[i][c]);
                else if (reconstruction == superbee)
                    slope = superbee_lim(vh[i][c]-vh[i-1][c], vh[i+1][c]-vh[i][c]);
                else if (reconstruction == minmod)
                    slope = minmod_lim(vh[i][c]-vh[i-1][c], vh[i+1][c]-vh[i][c]);
                else if (reconstruction == upwind)
                    slope = upwind_lim(vh[i][c]-vh[i-1][c], vh[i+1][c]-vh[i][c]);
                else
                    slope = 0.0; // First order 

                vbnd[i][0][c] = vh[i][c] - 0.5*slope;   
                vbnd[i][1][c] = vh[i][c] + 0.5*slope;
            }            


            V[c]   = vh[i][c]; 
            V_x[c] = r1_dx*slope;  

            VL[c] = vbnd[i][0][c];
            VR[c] = vbnd[i][1][c];
            
        }

        calcFlucSrcHLLC(VL,VR,BgradQ);

        for (c = 0; c < nVar; ++c)
            dqh[i][c] = -r1_dx*BgradQ[c];
        
    }

    // Find upwind flux

    for (i = 0; i < IMAX + 1; ++i) {
      
        for (c = 0; c < nVar; ++c) {
            VL[c] = vbnd[i-1][1][c];
            VR[c] = vbnd[i][0][c];
        }
        
        s = calcFluctuationsHLLC(VL,VR,Flux,Fluc); if (s > s_max) s_max = s; 

        for (c = 0; c < nVar; ++c) {
            Dm[i][c] = Flux[c]; Dp[i][c] = Fluc[c];
        }
        
    }

    // Find RHS

    for (i = 0; i < IMAX; ++i) {       
        
        for (c = 0; c < nVar; ++c) 
            dqh[i][c] += -r1_dx*( Dm[i+1][c] + Dp[i][c]  ); 
        
    }
}

//----------------------------------------------------------------------------
// One step of SSPRK(2,2)
//----------------------------------------------------------------------------

void Six_Equation::stepSSPRK22() {

    int i, c; 

    computeRHS(); 
    
    dt = CFL*dx/s_max; 

    if (time+dt > tend) 
        dt = tend-time; 

    for (i = 0; i < IMAX; ++i) {
        for (c = 0; c < nVar; ++c) {
            qh0[i][c] = qh[i][c]; 
            qh[i][c] +=  dt*dqh[i][c]; 
        }
    }

    relaxPressure(); 
   
    computeRHS(); 

    for (i = 0; i < IMAX; ++i) {
        for (c = 0; c < nVar; ++c) { 
            qh[i][c] =  0.5*(qh0[i][c] + qh[i][c] + dt*dqh[i][c]); 
        }
    }

    relaxPressure(); 
}

//----------------------------------------------------------------------------
// Plot solution 
//----------------------------------------------------------------------------

void Six_Equation::plot() const {
    
    std::ofstream out_data;
    out_data.open (filename);
    out_data.flags( std::ios::dec | std::ios::scientific );
    out_data.precision(16);

    out_data << "x,rho,u,p,alpha_1\n";

    for (int i = 0; i < IMAX; ++i)
        out_data << x[i] << "," << vh[i][0]+vh[i][1] << "," << vh[i][2] << "," << vh[i][3] << "," << vh[i][4] << std::endl;

    out_data.close();
}

//----------------------------------------------------------------------------
// Put everything together and run the problem 
//----------------------------------------------------------------------------

void Six_Equation::run() {

    computeRHS(); 

    
    std::cout << time_step << " " << time << std::endl;

    
    while(time < tend) {
        stepSSPRK22(); 
        time += dt; 
        time_step ++; 
        std::cout << time_step << " " << time << std::endl;  
    }
   
}

//----------------------------------------------------------------------------
// Convert conservative variable to primitive variables
//----------------------------------------------------------------------------

void Six_Equation::cons2Prim(const multi_array<double,1>& Q, multi_array<double,1>& V) const {

    double prho_1  = Q[0]; // Partial Density of phase 1 
	double prho_2  = Q[1]; // Partial Density of phase 2 
	double rhou    = Q[2]; // x-momentum of the mixture 
	double pE_1    = Q[3]; // Partial total energy of phase 1 
	double pE_2    = Q[4]; // Partial total energy of phase 2 
	double alpha_1 = Q[5]; // Volume fraction of phase 1 

	double alpha_2 = 1.0-alpha_1;
	double rho     = prho_1 + prho_2;
	double E       = pE_1 + pE_2;
	double u       = rhou/rho;  
	double ke      = 0.5*rhou*u;
	double e       = E - ke;

	double GAMMA   = alpha_1*GAMMA_1 + alpha_2*GAMMA_2;
	double PI      = alpha_1*PI_1 + alpha_2*PI_2;

	double prs = (e - PI)/GAMMA;

    V[0] = prho_1;     // Partial density of phase 1
    V[1] = prho_2;     // Partial density of phase 2
    V[2] = u;          // Mixture x-velocity
    V[3] = prs;        // Mixture pressure
    V[4] = alpha_1;    // Volume fraction of phase 1 
    V[5] = rho;        // Mixture density
}

//----------------------------------------------------------------------------
// Convert primitive variable to conservative variables
//----------------------------------------------------------------------------

void Six_Equation::prim2Cons(const multi_array<double,1>& V, multi_array<double,1>& Q) const {

    double prho_1  = V[0];
    double prho_2  = V[1];
    double u       = V[2];
    double prs     = V[3];
    double alpha_1 = V[4];

    double alpha_2 = 1.0 - alpha_1;
    double rho     = prho_1 + prho_2;
    double rhou    = rho*u; 
    double pE_1    = alpha_1*(prs*GAMMA_1 + PI_1) + 0.5*prho_1*u*u;
    double pE_2    = alpha_2*(prs*GAMMA_2 + PI_2) + 0.5*prho_2*u*u;

    Q[0] = prho_1;
    Q[1] = prho_2;
    Q[2] = rhou;
    Q[3] = pE_1;
    Q[4] = pE_2;
    Q[5] = alpha_1;
}

//----------------------------------------------------------------------------
// Speed of Sound
//----------------------------------------------------------------------------

double Six_Equation::speedOfSound(const multi_array<double,1>& V) {
    
    double rho     = V[0] + V[1];
	double prs     = V[3];
	double alpha_1 = V[4];
	double alpha_2 = 1-alpha_1;

    double a_sq = alpha_1*(1.0/GAMMA_1 + 1.0)*(prs + PI_1/(GAMMA_1 + 1.0)) +
                     alpha_2*(1.0/GAMMA_2 + 1.0)*(prs + PI_2/(GAMMA_2 + 1.0)) ;

    a_sq = a_sq/rho;

    if (a_sq < 0.0) {
        printf("Negative square speed of sound = %f\n", a_sq);
        printf("Volume fraction (phase 1) = %f\n", alpha_1); 
        printf("Density = %f\n", rho);
        printf("Pressure = %f\n", prs);
                
        std::exit(EXIT_FAILURE);
    }

    return std::sqrt(a_sq);
}


//----------------------------------------------------------------------------
// Calculate Star State
//----------------------------------------------------------------------------

void Six_Equation::calcStarStateHLLC(const multi_array<double,1>& V_L, const multi_array<double,1>& V_R, 
                      double& s_L, double& s_Star, double& s_R, 
                      multi_array<double,1>& W_L, multi_array<double,1>& W_Star, multi_array<double,1>& W_R) {


    int c; 

    double Q_L[nVar], Q_R[nVar], Q_L_star[nVar], Q_R_star[nVar]; 


    // Extract states

	double prho_1_L  = V_L[0]; double prho_1_R  = V_R[0];
	double prho_2_L  = V_L[1]; double prho_2_R  = V_R[1];
	double u_L       = V_L[2]; double u_R       = V_R[2];
	double prs_L     = V_L[3]; double prs_R     = V_R[3];
	double alpha_1_L = V_L[4]; double alpha_1_R = V_R[4];

	double alpha_2_L = 1. - alpha_1_L;      double alpha_2_R = 1. - alpha_1_R;
	double rho_L     = prho_1_L + prho_2_L; double rho_R     = prho_1_R + prho_2_R;

	double pE_1_L = alpha_1_L*(prs_L*GAMMA_1 + PI_1) + 0.5*prho_1_L*(u_L*u_L);
	double pE_2_L = alpha_2_L*(prs_L*GAMMA_2 + PI_2) + 0.5*prho_2_L*(u_L*u_L);

	double pE_1_R = alpha_1_R*(prs_R*GAMMA_1 + PI_1) + 0.5*prho_1_R*(u_R*u_R);
	double pE_2_R = alpha_2_R*(prs_R*GAMMA_2 + PI_2) + 0.5*prho_2_R*(u_R*u_R);

	double rho_1_L = prho_1_L/alpha_1_L; double rho_2_L = prho_2_L/alpha_2_L;
	double rho_1_R = prho_1_R/alpha_1_R; double rho_2_R = prho_2_R/alpha_2_R;

	double E_1_L = pE_1_L/alpha_1_L; double E_2_L = pE_2_L/alpha_2_L;
	double E_1_R = pE_1_R/alpha_1_R; double E_2_R = pE_2_R/alpha_2_R;

	Q_L[0] = prho_1_L;  Q_R[0] = prho_1_R;
    Q_L[1] = prho_2_L;  Q_R[1] = prho_2_R;
    Q_L[2] = rho_L*u_L; Q_R[2] = rho_R*u_R;
    Q_L[3] = pE_1_L;    Q_R[3] = pE_1_R;
    Q_L[4] = pE_2_L;    Q_R[4] = pE_2_R;
    Q_L[5] = alpha_1_L; Q_R[5] = alpha_1_R;

	double a_L = speedOfSound(V_L); 
	double a_R = speedOfSound(V_R); 


	double S_L    = std::min(u_L-a_L, u_R-a_R);
	double S_R    = std::max(u_L+a_L, u_R+a_R);
	double S_Star = (prs_R-prs_L + rho_L*u_L*(S_L-u_L) - rho_R*u_R*(S_R-u_R)  )/(rho_L*(S_L-u_L) - rho_R*(S_R-u_R));

    double xi_L = (S_L - u_L)/(S_L - S_Star);
	double xi_R = (S_R - u_R)/(S_R - S_Star);

	Q_L_star[0] = prho_1_L*xi_L;
	Q_L_star[1] = prho_2_L*xi_L;
	Q_L_star[2] = rho_L*xi_L*S_Star;
	Q_L_star[3] = prho_1_L*xi_L*(E_1_L/rho_1_L + (S_Star - u_L)*(S_Star + prs_L/(rho_1_L*(S_L - u_L))));
	Q_L_star[4] = prho_2_L*xi_L*(E_2_L/rho_2_L + (S_Star - u_L)*(S_Star + prs_L/(rho_2_L*(S_L - u_L))));
	Q_L_star[5] = alpha_1_L;

	Q_R_star[0] = prho_1_R*xi_R;
	Q_R_star[1] = prho_2_R*xi_R;
	Q_R_star[2] = rho_R*xi_R*S_Star;
	Q_R_star[3] = prho_1_R*xi_R*(E_1_R/rho_1_R + (S_Star - u_R)*(S_Star + prs_R/(rho_1_R*(S_R - u_R))));
	Q_R_star[4] = prho_2_R*xi_R*(E_2_R/rho_2_R + (S_Star - u_R)*(S_Star + prs_R/(rho_2_R*(S_R - u_R))));
	Q_R_star[5] = alpha_1_R;

	for(c = 0; c < nVar; ++c) {
		W_L[c]     = Q_L_star[c] - Q_L[c];
		W_Star[c]  = Q_R_star[c] - Q_L_star[c];
		W_R[c]     = Q_R[c]      - Q_R_star[c];
	}


    s_L = S_L; s_Star = S_Star; s_R = S_R;  
}


//----------------------------------------------------------------------------
// Calculate fluctuation terms for the given left and right states (HLLC)
//----------------------------------------------------------------------------

double Six_Equation::calcFluctuationsHLLC(const multi_array<double,1>&  V_L, const multi_array<double,1>&  V_R, multi_array<double,1>& Dm, multi_array<double,1>& Dp) {
    
    int c; 
    multi_array<double,1> W_L(extents[nVar]), W_R(extents[nVar]), W_Star(extents[nVar]); 
    double s_L, s_R, s_Star;
    
    calcStarStateHLLC(V_L, V_R, s_L, s_Star, s_R, W_L, W_Star, W_R); 

	double s_L_m    = std::min(0.,s_L);    double s_L_p    = std::max(0.,s_L);
	double s_Star_m = std::min(0.,s_Star); double s_Star_p = std::max(0.,s_Star);
	double s_R_m    = std::min(0.,s_R);    double s_R_p    = std::max(0.,s_R);

	for(c = 0; c < nVar; ++c) {
		Dm[c] = s_L_m*W_L[c] +  s_Star_m*W_Star[c] + s_R_m*W_R[c];
		Dp[c] = s_L_p*W_L[c] +  s_Star_p*W_Star[c] + s_R_p*W_R[c];
	}


    return std::max(std::abs(s_L), std::max(std::abs(s_Star), std::abs(s_R)));
} 

//----------------------------------------------------------------------------
// Calculate fluctuation source term (HLLC)
//----------------------------------------------------------------------------

void Six_Equation::calcFlucSrcHLLC(const multi_array<double,1>&  V_L, const multi_array<double,1>&  V_R, multi_array<double,1>&  Flux_Src) {

    int c; 
    multi_array<double,1> W_L(extents[nVar]), W_R(extents[nVar]), W_Star(extents[nVar]); 
    double s_L, s_R, s_Star;
    
    calcStarStateHLLC(V_L, V_R, s_L, s_Star, s_R, W_L, W_Star, W_R); 

	for(c = 0; c < nVar; ++c)
		Flux_Src[c] = s_L*W_L[c] +  s_Star*W_Star[c] + s_R*W_R[c];
} 

//----------------------------------------------------------------------------
// Pressure relaxation procedure 
//----------------------------------------------------------------------------

void Six_Equation::relaxPressure() {

    int i, c;


    multi_array<double,1> Q(extents[nVar]);

    double prho_1, prho_2, rhou, pE_1, pE_2, alpha_1, alpha_2, rho, u, E_1, E_2, rho_1, rho_2, v_sq, prs_1, prs_2, aa, bb, dd, D;
    double e_1, e_2, prs;

    for (i = 0; i < IMAX; ++i) {
        
        for (c = 0; c < nVar; ++c) {
            Q[c] = qh[i][c];
        }


	    prho_1  = Q[0];
	    prho_2  = Q[1];
	    rhou    = Q[2];
	    pE_1    = Q[3];
	    pE_2    = Q[4];
	    alpha_1 = Q[5];
	    alpha_2 = 1.0 - alpha_1;

        
	    rho = prho_1 + prho_2;
	    u = rhou/rho;

	    E_1 = pE_1/alpha_1;
	    E_2 = pE_2/alpha_2;

	    rho_1 = prho_1/alpha_1;
	    rho_2 = prho_2/alpha_2;

        v_sq = u*u; 

	    prs_1 = (gamma_1-1.)*(E_1 - 0.5*rho_1*v_sq) - gamma_1*pi_1;
	    prs_2 = (gamma_2-1.)*(E_2 - 0.5*rho_2*v_sq) - gamma_2*pi_2;


        aa =  gamma_1*alpha_2 + gamma_2*alpha_1;
        bb = -gamma_1*alpha_2*(prs_2-pi_1) - gamma_2*alpha_1*(prs_1-pi_2);
        dd = -gamma_1*alpha_2*prs_2*pi_1 - gamma_2*alpha_1*prs_1*pi_2;

        D = bb*bb - 4.*aa*dd;
        
        if (D >= 0.0) {
            
            D = std::sqrt(D);

            if (bb > 0.0)
	            prs = 2.*dd/(-bb-D);
            else
	            prs = (-bb+D)/(2.*aa);            

            alpha_1 = alpha_1*((gamma_1-1.)*prs + prs_1 + gamma_1*pi_1)/((gamma_1-1.)*prs + prs + gamma_1*pi_1);
            alpha_2 = 1. - alpha_1;

            e_1 = (prs + gamma_1*pi_1)/(gamma_1-1.);
            e_2 = (prs + gamma_2*pi_2)/(gamma_2-1.);

            qh[i][3] = alpha_1*e_1 + 0.5*prho_1*v_sq;
            qh[i][4] = alpha_2*e_2 + 0.5*prho_2*v_sq;
            qh[i][5] = alpha_1;
        }

    }
}

//----------------------------------------------------------------------------
// Errors 
//----------------------------------------------------------------------------

void Six_Equation::errors() const {
    
    double err_rho_l2 = 0., err_rho_max = 0.; 
    double err; 

    for (int i = 0; i < IMAX; ++i) {
        double rho   =  qh[i][0] + qh[i][1];
        double rho_e =  qhi[i][0] + qhi[i][1];

        err = std::abs(rho_e-rho);

        if (err > err_rho_max)
            err_rho_max = err; 

        err_rho_l2 += err*err; 
            
    }

    err_rho_l2 = std::sqrt(err_rho_l2/static_cast<double>(IMAX));

    printf("Density: L2  Error = %.5e\n", err_rho_l2);
    printf("Density: Max Error = %.5e\n", err_rho_max);    

}


//----------------------------------------------------------------------------
// Initial condition function 
//----------------------------------------------------------------------------

void Six_Equation::smoothProblem(double xx, multi_array<double,1>& Q0) {

    multi_array<double,1> V0(extents[nVar]);

    double rho_1 = 1.0;
    double rho_2 = 1.5 + 0.4*std::cos(M_PI*xx);
    double u = 1.0;
    double p = 1.0;    
    double alpha_1 = 0.5 + 0.4*std::sin(M_PI*xx);
    double alpha_2 = 1.0-alpha_1; 

    V0[0] = alpha_1*rho_1;
    V0[1] = alpha_2*rho_2;
    V0[2] = u;
    V0[3] = p; 
    V0[4] = alpha_1;
    
    prim2Cons(V0,Q0);
}

void Six_Equation::airHelium(double xx, multi_array<double,1>& Q0) {

    multi_array<double,1> V0(extents[nVar]);

    if (xx < 0.0) {

        V0[4] = 1.0-eps_vol;

        V0[0] = 1.0*V0[4];
        V0[1] = 1.0*(1.-V0[4]);
        V0[2] = 0.0;
        V0[3] = 1.0;        
    }

    else {
    
        V0[4] = eps_vol;
        
        V0[0] = 1.0*V0[4];
        V0[1] = 0.125*(1.-V0[4]);
        V0[2] = 0.0;
        V0[3] = 0.1;        
    }

    prim2Cons(V0,Q0);
}

void Six_Equation::multiComponentShuOsher(double xx, multi_array<double,1>& Q0) {

    multi_array<double,1> V0(extents[nVar]);

    double eps_he = 1.0e-8;
    double eps_air = 1.0-eps_he;
    

    if (xx <= -4.) {
    
        V0[0] = 3.857143*eps_air;
        V0[1] = (1 + 0.2*std::sin(5.*xx))*eps_he;
        V0[2] = 2.629369;
        V0[3] = 31./3.;
        V0[4] = eps_air; 
    }

    else {

        V0[0] = 3.857143*(1.-eps_air);
        V0[1] = (1 + 0.2*std::sin(5.*xx))*(1.-eps_he);
        V0[2] = 0.;
        V0[3] = 1.;
        V0[4] = eps_he; 
    }

    prim2Cons(V0,Q0);
}



void Six_Equation::interfaceAdvection(double xx, multi_array<double,1>& Q0) {

    multi_array<double,1> V0(extents[nVar]);

    double eps = 1.0e-8; 

    double rho_1 = 1.0; 
    double rho_2 = 0.01; 
    double u = 10.0; 
    double p = 1.0;   
    
    V0[2] = u;  
    V0[3] = p;
    
    if(xx<0){
        
        V0[0] = (1.0-eps)*rho_1;
        V0[1] = eps*rho_2;
        V0[4] = 1.0-eps;
    }
    
    else
    {
        V0[0] = eps*rho_1;
        V0[1] = (1.0-eps)*rho_2;
        V0[4] = eps;
    
    }
   

    prim2Cons(V0,Q0);
}

void Six_Equation::airWater(double xx, multi_array<double,1>& Q0) {

    multi_array<double,1> V0(extents[nVar]);

    double rho_water = 1.0, rho_air = 0.001;

    if (xx < 2.5) {
    
        V0[2] = 0.0;
        V0[3] = 1.0e4;
        V0[4] = 1.0-1.0e-6; 

        V0[0] = V0[4]*rho_water;
        V0[1] = (1.0-V0[4])*rho_air;
    }

    else {

        V0[2] = 0.0;
        V0[3] = 1.0;
        V0[4] = 1.0e-6;

        V0[0] = V0[4]*rho_water;
        V0[1] = (1.0-V0[4])*rho_air;
    }

    prim2Cons(V0,Q0);
}

void Six_Equation::waterCavitation(double xx, multi_array<double,1>& Q0) {

    multi_array<double,1> V0(extents[nVar]);
    
    double alpha_v = 1.0e-2;
    double alpha_l = 1.0 - alpha_v;
    double p = 1.0; 
    double rho_v = 0.00063;
    double rho_l = 1.150;  
    double u = 0.2; 
     
    V0[0] = alpha_l*rho_l;
    V0[1] = alpha_v*rho_v;

    if (xx < 0.0) 
        V0[2] = -1.0*u;
    else   
        V0[2] = u; 

    V0[3] = p;
    V0[4] = alpha_l;
    
    prim2Cons(V0,Q0);
}

void Six_Equation::dodecaneVaporMixture(double xx, multi_array<double,1>& Q0) {

    multi_array<double,1> V0(extents[nVar]);
    
    double alpha_1; 

    if (xx < 2.5) {

        alpha_1 = 0.7; 
        V0[0] = alpha_1*0.5; 
        V0[1] = (1.0-alpha_1)*0.002;
        V0[2] = 0.0; 
        V0[3] = 1000.0;
        V0[4] = alpha_1; 
     }
    else   {

        alpha_1 = 0.3; 
        V0[0] = alpha_1*0.5; 
        V0[1] = (1.0-alpha_1)*0.002;
        V0[2] = 0.0; 
        V0[3] = 1.0;
        V0[4] = alpha_1; 
    } 
    

    prim2Cons(V0,Q0);
}

void Six_Equation::epoxySpinel(double xx, multi_array<double,1>& Q0) {

    multi_array<double,1> V0(extents[nVar]);
    
	double rho_epoxy = 1185.0;
	double rho_spinel = 3622.0;
	double alpha_epoxy = 0.5954;
	double alpha_spinel = 1.0-alpha_epoxy;


	V0[0] = alpha_epoxy*rho_epoxy;
	V0[1] = alpha_spinel*rho_spinel;
	V0[2] = 0.0;
	V0[4] = alpha_epoxy;

	if (xx < 0.6)
		V0[3] = 1.0e10;
	else
		V0[3] = 1.0e5;
    

    prim2Cons(V0,Q0);
}

void Six_Equation::initCond(double xx, multi_array<double,1>& Q0) {
    
    if (icond == smooth_problem)
        smoothProblem(xx,Q0); 
    else if (icond == air_helium)
        airHelium(xx,Q0); 
    else if (icond == multi_component_shu_osher)
        multiComponentShuOsher(xx,Q0);
    else if (icond == water_cavitation)
        waterCavitation(xx,Q0);
    else if (icond == interface_advection)
        interfaceAdvection(xx,Q0);
    else if (icond == dodecane_vapor_mixture)
        dodecaneVaporMixture(xx,Q0);
    else if (icond == epoxy_spinel)
        epoxySpinel(xx,Q0);
    else
        airWater(xx,Q0);
}
