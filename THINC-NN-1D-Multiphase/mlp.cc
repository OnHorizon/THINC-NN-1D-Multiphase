/*
 * mlp.cc
 *      Author: sunder
 */

#include "hype.h"

//------------------------------------------------------------------------------------
// Soft-Max Activation function 
//----------------------------------------------------------------------------------

void MultiLayerPerceptron::softMax(const gsl_vector* x, gsl_vector* out) {
    double total = 0.0; 
    for (int i = 0; i < n_output_neurons; ++i) {
        gsl_vector_set(out,i,std::exp(gsl_vector_get(x,i))); 
        total += gsl_vector_get(out,i);    
    }

    gsl_vector_scale(out,1./total); 
}


//------------------------------------------------------------------------------------
// Constructor. Read and store the weights matrices and bias vectors 
//------------------------------------------------------------------------------------

MultiLayerPerceptron::MultiLayerPerceptron() {

    int i, j, l;
    std::string line;   
	std::string buffer;
    std::string filename; 
    double temp; 
     
    // 1) Get the size of the neural network     

    std::ifstream size_file("Neural-Network/sizes.csv");

	if ( !(size_file.is_open()) ) {
		std::cerr << "Error. Unable to open file" << std::endl;
		std::exit(EXIT_FAILURE);
	}

    std::getline(size_file, line);
	std::istringstream instr(line);
    instr >> n_input_neurons;
    instr.clear();
	
    std::getline(size_file, line);
    instr.str(line);
    instr >> n_hidden_layers;
    instr.clear();

    n_hidden_neurons = new int(n_hidden_layers); 

    for (l = 0; l < n_hidden_layers; ++l) {
        std::getline(size_file, line);
        instr.str(line);
        instr >> n_hidden_neurons[l];
        instr.clear();
    }

    std::getline(size_file, line);
    instr.str(line);
    instr >> n_output_neurons;
    instr.clear();

    size_file.close(); 

    rows = new int[n_hidden_layers+1];
    cols = new int[n_hidden_layers+1]; 
    
    cols[0] = n_input_neurons; 

    for (l = 1; l < n_hidden_layers+1; ++l) 
        cols[l] = n_hidden_neurons[l-1]; 

    for (l = 0; l < n_hidden_layers; ++l)
        rows[l] = n_hidden_neurons[l]; 

    rows[n_hidden_layers] = n_output_neurons;

    // 2) Assign memory for the weight matrices and bias vectors 

    W = new gsl_matrix*[n_hidden_layers+1];
    b = new gsl_vector*[n_hidden_layers+1];
    z = new gsl_vector*[n_hidden_layers+1];
    a = new gsl_vector*[n_hidden_layers+1];

    for (l = 0; l <= n_hidden_layers; ++l) {
         W[l] = gsl_matrix_alloc (rows[l],cols[l]);
         b[l] = gsl_vector_alloc (rows[l]);
         z[l] = gsl_vector_alloc (rows[l]); 
         a[l] = gsl_vector_alloc (rows[l]);    
    }
   
    // 3) Read weights from the file and fill the weight matrices  
    
    for (l = 0; l <= n_hidden_layers; ++l) {
            
        filename = "Neural-Network/weight-" + std::to_string(l) + ".csv";

        std::ifstream file(filename);

    	if ( !(file.is_open()) ) {
		    std::cerr << "Error. Unable to open file " << filename << std::endl;
		    std::exit(EXIT_FAILURE);
	    }

        for (i = 0; i < rows[l]; ++i) {
            std::getline(file, line);
            instr.str(line);           
            for (j = 0; j < cols[l]; ++j) {
                instr >> temp;
                gsl_matrix_set(W[l],i,j,temp); 
            }
            instr.clear();
        }

        file.close();         
    }
        
    // 4) Read intercepts from the file and fill the bias vectors
    
    for (l = 0; l <= n_hidden_layers; ++l) {
            
        filename = "Neural-Network/bias-" + std::to_string(l) + ".csv";

        std::ifstream file(filename);

    	if ( !(file.is_open()) ) {
		    std::cerr << "Error. Unable to open file " << filename << std::endl;
		    std::exit(EXIT_FAILURE);
	    }

        for (j = 0; j < rows[l]; ++j) {
            std::getline(file, line);
            instr.str(line);
            instr >> temp;
            gsl_vector_set(b[l],j,temp);
            instr.clear();
        }

        file.close();         
    }
}

//------------------------------------------------------------------------------------
// Destructor 
//------------------------------------------------------------------------------------

MultiLayerPerceptron::~MultiLayerPerceptron() {
    
    for (int l = 0; l <= n_hidden_layers; ++l){
        gsl_matrix_free(W[l]); gsl_vector_free(b[l]);
        gsl_vector_free(z[l]); gsl_vector_free(a[l]);     
    }

    delete[] n_hidden_neurons; delete[] rows; delete[] cols;
    delete[] W; delete[] b; delete[] z; delete[] a; 
     
}

//------------------------------------------------------------------------------------
// Classify the given data 
//------------------------------------------------------------------------------------

double MultiLayerPerceptron::predict(gsl_vector* in) {

    int l, j; 

    for (l = 0; l <= n_hidden_layers; ++l) {
    
        // 1) Find the weighted sum

        gsl_blas_dcopy(b[l], z[l]); // Copy b into z 

        if (l == 0)
            gsl_blas_dgemv(CblasNoTrans, 1.0, W[l], in, 1.0, z[l]);
        else
            gsl_blas_dgemv(CblasNoTrans, 1.0, W[l], a[l-1], 1.0, z[l]);


        // 2) Apply activation 
        
        if (l == n_hidden_layers) {
            for (j = 0; j < n_output_neurons; ++j) 
                gsl_vector_set( a[l], j, sigmoid( gsl_vector_get(z[l], j) ) );            
        }
    
        else {
            for (j = 0; j < n_hidden_neurons[l]; ++j) 
                gsl_vector_set( a[l], j, sigmoid( gsl_vector_get(z[l], j) ) ); 
        }
    }

    return (1.6-std::log(3.))*gsl_vector_get(a[n_hidden_layers],0) + std::log(3.);
}



