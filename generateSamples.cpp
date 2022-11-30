#include <iostream>
#include <math.h>
#include <cmath>
#include <cstdlib>
#include <stdio.h>
#include <random>
#include <iomanip>
#include <map>

using namespace std;

int main(int argc, char *argv[]){

    // number of samples per experiment matrix per uncertain parameter
    int n = atoi(argv[1]); 
    // fileout filenames for A and B matrices
    char *fileoutA = argv[2];
    char *fileoutB = argv[3];
    FILE * fout;
    // number of uncertain parameters (dimension of parameter space)
    static const int d=4;

    // upper/lower bounds (uniform) and mean/variance (Gaussian) for UQ parameters
    static const double lbTairIndoors    = 293.150;
    static const double ubTairIndoors    = 298.706;
    static const double lbTsrmIndoors    = 293.150;
    static const double ubTsrmIndoors    = 298.706;
    static const double meanTairOutdoors = 285.928;
    static const double varTairOutdoors  = 10.000 ;
    static const double meanTsrmOutdoors = 285.928;
    static const double varTsrmOutdoors  = 10.000 ;

    // creating first experiment matrix A
    double A[d][n];

    // creating second experiment matrix B
    double B[d][n];
    for (int id=0; id<d; id++) {
        for (int in=0; in<n; in++) {
            A[id][in]=0;
            B[id][in]=0;
        }
    }

    // initializing random seed for rng
    std::random_device rd; // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd());  // Standard mersenne_twister_engine seeded with rd()

    // drawing uniform random samples for indoor temperatures
    std::uniform_real_distribution<double> unif(lbTairIndoors, ubTairIndoors);
    std::normal_distribution<double> normal(meanTairOutdoors,sqrt(varTairOutdoors));
    for (int in=0; in<n; in++) {
        A[0][in] = unif(gen); B[0][in] = unif(gen);
        A[1][in] = unif(gen); B[1][in] = unif(gen);
        A[2][in] = normal(gen); B[2][in] = normal(gen);
        A[3][in] = normal(gen); B[3][in] = normal(gen);
    }

    // writing matrix A to vecfileoutA
    fout = fopen(fileoutA, "w");
    if (fout == NULL) {
        printf("Error opening output file !");
        exit(EXIT_FAILURE);
    }
    for (int in=0; in<n; in++) {
        fprintf(fout,"%.17g,",A[0][in]);
        fprintf(fout,"%.17g,",A[1][in]);
        fprintf(fout,"%.17g,",A[2][in]);
        fprintf(fout,"%.17g\n",A[3][in]);
    }
    // writing matrix B to vecfileoutB   
    fout = fopen(fileoutB, "w");
    if (fout == NULL) {
        printf("Error opening output file !");
        exit(EXIT_FAILURE);
    }
    for (int in=0; in<n; in++) {
        fprintf(fout,"%.17g,",B[0][in]);
        fprintf(fout,"%.17g,",B[1][in]);
        fprintf(fout,"%.17g,",B[2][in]);
        fprintf(fout,"%.17g\n",B[3][in]);
    }

}