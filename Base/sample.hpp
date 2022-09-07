

#ifndef __SAMPLE_HPP__
#define __SAMPLE_HPP__

#include <math.h>
#include <random>
#include <Eigen/Dense>

#include "vector3d.hpp"


// random number generator
static std::mt19937 mtgen{0};

// python programmer is required to initialise generator
int init_generator(long pyseed);

double r2();


Vector3d modified_vmf(double kappa, double L, double R);


// be careful with defining the random number generator, 
// should be done once during program execution
//cout << "Using random seed: " << seed << endl;
//double what_seed();

// random seed
//srand(time(0));

// number in curly brackets is the seed
// determined
//static void expose_srand(long pyseed) {
  //srand(pyseed);
//}


//https://stackoverflow.com/questions/5008804/generating-random-integer-from-a-range

std::uniform_int_distribution<int> init_uniform_rint(int min, int max);

Eigen::RowVectorXd phi_row(int size);

#endif
