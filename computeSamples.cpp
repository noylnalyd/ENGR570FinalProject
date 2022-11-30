#include <iostream>
#include <math.h>
#include <cmath>
#include <cstdlib>
#include <stdio.h>
#include <fstream>
#include <string>
#include <sstream>

/// @brief Script that runs a given set of rows from A and B matrices and computes model outputs.
int main(int argc, char *argv[]){

    // start and end row indices (inclusive, 0-indexed) to consider in .txt files
    int iStart = atoi(argv[1]); 
    int iStop = atoi(argv[2]); 
    // total number of rows to read
    int nrows = iStop-iStart+1; 
    // fileout filenames for A and B matrices (to read and write from)
    char *fileoutA = argv[3];
    char *fileoutB = argv[4];
    // UQ variables
    int d=4; // dimension of parameter space
    // parallelization variables
    int pID = atoi(argv[5]);
    // variables for reading and writing
    FILE * fp, * fpB;
    size_t len = 0;
    size_t lenB = 0;
    ssize_t read;
    int linenum=0;
    char * line = NULL; char * lineB = NULL;
    char * token, * end, * endB;
    static const char separator[] = ",";

    // checks
    if (iStart>iStop) {
        std::cout << "iStart must be <= iStop!\n";
        exit(EXIT_FAILURE);
    }

    // creating A matrix
    double* A[nrows];
    for (int i=0; i<nrows; i++) {
        A[i] = (double*)calloc(d, sizeof(double));
    }    
    // initializing to 0 to be safe
    for (int id=0; id<d; id++) {
        for (int in=0; in<nrows; in++) {
            A[in][id]=0;
        }
    }
    // reading in A matrix
    fp = fopen(fileoutA,"r");
    if (fp == NULL) {
        exit(EXIT_FAILURE);
    }
    int in=0; //counts linenum w.r.t. local A matrix
    std::cout << "A matrix:\n";
    while ((read = getline(&line, &len, fp)) != -1) {
        // only read the lines of interest, skip the rest
        if ((linenum<=iStop) && (linenum>=iStart)) {
            int id=0;
            token=strtok(line,separator);
            A[in][id]=strtod(token, &end); id++; std::cout << A[in][id-1] << ' ';
            token=strtok(NULL,separator);
            A[in][id]=strtod(token, &end); id++; std::cout << A[in][id-1] << ' ';
            token=strtok(NULL,separator);
            A[in][id]=strtod(token, &end); id++; std::cout << A[in][id-1] << ' ';
            token=strtok(NULL,separator);
            A[in][id]=strtod(token, &end); id++;
            std::cout << A[in][id-1] << std::endl;
            in++; 
        }
        linenum++;
    }

    // creating B matrix
    double* B[nrows];
    for (int i=0; i<nrows; i++) {
        B[i] = (double*)calloc(d, sizeof(double));
    }    
    // initializing to 0 to be safe
    for (int id=0; id<d; id++) {
        for (int in=0; in<nrows; in++) {
            B[in][id]=0;
        }
    }
    // reading in B matrix
    fpB = fopen(fileoutB,"r");
    if (fpB == NULL) {
        std::cout << "Can't read from b file\n";
        exit(EXIT_FAILURE);
    }
    in=0; //counts linenum w.r.t. local A matrix
    linenum=0;
    std::cout << "B matrix:\n";
    while ((read = getline(&lineB, &lenB, fpB)) != -1) {
        // only read the lines of interest, skip the rest
        if ((linenum<=iStop) && (linenum>=iStart)) {
            int id=0;
            token=strtok(lineB,separator);
            B[in][id]=strtod(token, &endB); id++; std::cout << B[in][id-1] << ' ';
            token=strtok(NULL,separator);
            B[in][id]=strtod(token, &endB); id++; std::cout << B[in][id-1] << ' ';
            token=strtok(NULL,separator);
            B[in][id]=strtod(token, &endB); id++; std::cout << B[in][id-1] << ' ';
            token=strtok(NULL,separator);
            B[in][id]=strtod(token, &endB); id++;
            std::cout << B[in][id-1] << std::endl;
            in++; 
        }
        linenum++;
    }

    // creating A_B matrices, d sets of A_B[nrows][d] (3rd dimension is i index per A_B^i notation)
    double *** AB = new double**[nrows];
    for (int in = 0; in < nrows; in++) {
        AB[in] = new double*[d];
        for (int id = 0; id < d; id++) {
            AB[in][id] = new double[d];
        }
    }
    // setting all columns but column=kd to same as A
    for (int kd=0; kd<d; kd++) {
        for (int id=0; id<d; id++) {
            for (int in=0; in<nrows; in++) {
                if (id==kd) {
                    AB[in][id][kd]=B[in][id];
                }
                else {
                    AB[in][id][kd]=A[in][id];
                }
            }
        }
    }
    for (int kd=0; kd<d; kd++) {
        std::cout << "AB_" << kd <<std::endl;
        for (int in=0; in<nrows; in++) {
            std::cout << AB[in][0][kd] << ", ";
            std::cout << AB[in][1][kd] << ", ";
            std::cout << AB[in][2][kd] << ", ";
            std::cout << AB[in][3][kd] << std::endl;
        }
    }


    // initialize Simulator, BodyModel, and SimModel structures
    //Simulator* mySim = new Simulator();
    //BodyModel* myBody = defaultBody();
    //SimModel* mySimModel = new SimModel();
    //mySim->initializer();

    // run each case

    for (int i=0; i<nrows; i++) {
        free(A[i]); 
        free(B[i]);
        for (int id=0; id<d; id++) {
            free(AB[i][id]);
        }
    }
}