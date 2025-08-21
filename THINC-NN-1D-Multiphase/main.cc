/*
 * main.cc
 *      Author: sunder
 */

#include "hype.h"

int main() {

    /*
     MultiLayerPerceptron MLP;

     gsl_vector* delta = gsl_vector_alloc(4);

     gsl_vector_set(delta,0,3.961779076302204272e-01);
     gsl_vector_set(delta,1,9.999999999999967804e-01);
     gsl_vector_set(delta,2,6.980889538151086038e-01);
     gsl_vector_set(delta,3,6.038220923697766862e-01);
    
    std::cout << MLP.predict(delta)  << std::endl;
    
    gsl_vector_free(delta);
    */
   
    
    int IMAX = 256;
    Six_Equation Sol(IMAX, THINC, interface_advection);

    Sol.run();   
    Sol.errors();
    Sol.plot();
    //Sol.beta_plot();
    
    return 0;
}
