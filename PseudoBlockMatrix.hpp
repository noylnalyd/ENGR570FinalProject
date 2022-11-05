#ifndef _PSEUDOBLOCKMATRIX_
#define _PSEUDOBLOCKMATRIX_

#include <iostream>
#include <cstddef>
#include <math.h>

using namespace std;

namespace PSEUDOBLOCKMATRIX
{
    class PseudoBlockMatrix
    {
    private:
        double ***blocks; // 2D dense blocks
        double *denseRow; // 1D dense final row
        double *denseCol; // 1D dense final column
        double denseCorner; // 0D final corner value
        double *rowStart, *colStart; // Indices of first row/col per block. For binsearch of get/add/set

    public:
        int N,M; // Number of rows and cols
        int *blockN, *blockM; // Number of rows and cols per block
        int Nblocks; // Number of blocks

        PseudoBlockMatrix(
            int tNblocks,
            int* tblockN,
            int* tblockM)
        {
            Nblocks = tNblocks;
            blockN = tblockN;
            blockM = tblockM;

            N = 0;
            M = 0;

            // Allocate blocks, compute rowstart / colStart
            blocks = new double**[Nblocks];
            rowStart = new double[Nblocks];
            colStart = new double[Nblocks];
            for(int b=0;b<Nblocks;b++){
                blocks[b] = new double*[blockN[b]];
                for(int i=0;i<blockN[b];i++){
                    blocks[b][i] = new double[blockM[b]];
                }
                rowStart[b] = N;
                N = N + blockN[b];
                colStart[b] = M;
                M = M + blockM[b];
            }
            // Allocate dense col, row, corner
            denseCol = new double[N-1];
            denseRow = new double[M-1];
            denseCorner = 0;

            // Adjust for dense row/col!
            N = N+1;
            M = M+1;

        }

        ~PseudoBlockMatrix()
        {
            for(int b=0;b<Nblocks;b++){
                for(int i=0;i<blockN[b];i++){
                    delete [] blocks[b][i];
                }
            }
            for(int b=0;b<Nblocks;b++){
                delete [] blocks[b];
            }
            delete [] blocks;

            delete [] denseRow;
            delete [] denseCol;
            delete [] rowStart;
            delete [] colStart;
            delete [] blockN;
            delete [] blockM;
        }

    };
}

#endif