#include <iostream>
#include <math.h>
#include <cmath>
#include <cstdlib>
#include <stdio.h>
#include <random>

using namespace std;

int main(int argc, char *argv[]){

    // number of samples per experiment matrix per uncertain parameter
    int n = atoi(argv[1]); 
    // fileout filenames for A, B, and C matrices
    char *fileoutA = argv[2];
    char *fileoutB = argv[3];
    char *fileoutC = argv[4];
    // number of cores to partition samples onto
    int nCores = atoi(argv[5]);
    // fileout filename for partitioning index vector
    char *fileoutP = argv[6];
    FILE * fout;
    // number of uncertain parameters (dimension of parameter space)
    static const int d=4;

    // distribution constants for UQ parameters  
    // lines below test case of analytically known indices
    static const double meanTsrmIndoors  = 1;
    static const double stdTsrmIndoors   = 1;
    static const double meanTsrmOutdoors = 3;
    static const double stdTsrmOutdoors  = 1.5;
    static const double meanTToRecovery  = 5;
    static const double stdTToRecovery   = 2;
    static const double meanTMachine     = 7;
    static const double stdTMachine      = 2.5;
    // uncomment the lines before for equal indices:
/*    static const double meanTsrmIndoors  = 1;
    static const double stdTsrmIndoors   = 1;
    static const double meanTsrmOutdoors = 1;
    static const double stdTsrmOutdoors  = 1;
    static const double meanTToRecovery  = 1;
    static const double stdTToRecovery   = 1;
    static const double meanTMachine     = 1;
    static const double stdTMachine      = 1;*/

    // creating first experiment matrix A
    double A[d][n];

    // creating second experiment matrix B
    double B[d][n];

    // initializing to 0 to be safe
    for (int id=0; id<d; id++) {
        for (int in=0; in<n; in++) {
            A[id][in]=0;
            B[id][in]=0;
        }
    }

    // initializing rng
    std::random_device rd; 
    std::mt19937 gen(rd());  // Standard mersenne_twister_engine seeded with rd()

    // drawing uniform random samples for temperatures and saving to arrays
    std::normal_distribution<double> norm1(meanTsrmIndoors,stdTsrmIndoors);
    std::normal_distribution<double> norm2(meanTsrmOutdoors,stdTsrmOutdoors);
    std::normal_distribution<double> norm3(meanTToRecovery,stdTToRecovery);
    std::normal_distribution<double> norm4(meanTMachine,stdTMachine);
    for (int in=0; in<n; in++) {
        A[0][in] = norm1(gen); B[0][in] = norm1(gen);
        A[1][in] = norm2(gen); B[1][in] = norm2(gen);
        A[2][in] = norm3(gen); B[2][in] = norm3(gen);
        A[3][in] = norm4(gen); B[3][in] = norm4(gen);
    }

    // writing matrix A to fileoutA
    fout = fopen(fileoutA, "w");
    if (fout == NULL) {
        printf("Error opening output file A !");
        exit(EXIT_FAILURE);
    }
    for (int in=0; in<n; in++) {
        fprintf(fout,"%.17g,",A[0][in]);
        fprintf(fout,"%.17g,",A[1][in]);
        fprintf(fout,"%.17g,",A[2][in]);
        fprintf(fout,"%.17g\n",A[3][in]);
    }
    // writing matrix B to fileoutB   
    fout = fopen(fileoutB, "w");
    if (fout == NULL) {
        printf("Error opening output file B !");
        exit(EXIT_FAILURE);
    }
    for (int in=0; in<n; in++) {
        fprintf(fout,"%.17g,",B[0][in]);
        fprintf(fout,"%.17g,",B[1][in]);
        fprintf(fout,"%.17g,",B[2][in]);
        fprintf(fout,"%.17g\n",B[3][in]);
    }

    // writing matrix C to fileoutC
    fout = fopen(fileoutC, "w");
    if (fout == NULL) {
        printf("Error opening output file C !");
        exit(EXIT_FAILURE);
    }
    for (int in=0; in<n; in++) {
        fprintf(fout,"%.17g,",A[0][in]);
        fprintf(fout,"%.17g,",A[1][in]);
        fprintf(fout,"%.17g,",A[2][in]);
        fprintf(fout,"%.17g\n",A[3][in]);
    }
    for (int in=0; in<n; in++) {
        fprintf(fout,"%.17g,",B[0][in]);
        fprintf(fout,"%.17g,",B[1][in]);
        fprintf(fout,"%.17g,",B[2][in]);
        fprintf(fout,"%.17g\n",B[3][in]);
    }

    // partitioning samples and writing to partition file
    int nPerProc = n/nCores;
    int remainder = n % nCores;
    fout = fopen(fileoutP,"w");
    if (fout == NULL) {
        printf("Error opening output file P !");
        exit(EXIT_FAILURE);
    }
    int row=0; int iStart=0; int iStop=0;
    int nThisProc=nPerProc;
    for (int ic=0; ic<nCores; ic++) {
        if (ic<remainder) {
            nThisProc=nPerProc+1;
        }
        else {
            nThisProc=nPerProc;
        }
        iStart=row;
        row=row+nThisProc;
        iStop=row-1;
        //fprintf(fout,"%d ",ic);
        fprintf(fout,"%d, ",iStart);
        fprintf(fout,"%d \n",iStop);
    }
}
