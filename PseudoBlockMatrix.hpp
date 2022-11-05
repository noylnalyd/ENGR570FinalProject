#ifndef _PSEUDOBLOCKMATRIX_
#define _PSEUDOBLOCKMATRIX_

#include <iostream>
#include <cstddef>
#include <assert.h>
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

        int binSearchBlock( int rowIdx );
        double getOverBlocks( int i, int j );
        double getFromBlock( int b, int i, int j );
        double getOverDense( int i, int j );
        double get( int i , int j );
        bool setOverBlocks( int i, int j, double val );
        bool setFromBlock( int b, int i, int j, double val );
        bool setOverDense( int i, int j, double val );
        bool set( int i , int j, double val );
        bool addOverBlocks( int i, int j, double val );
        bool addFromBlock( int b, int i, int j, double val );
        bool addOverDense( int i, int j, double val );
        bool add( int i , int j, double val );
        void fill( double val );
        PseudoBlockMatrix shallowCopy( PseudoBlockMatrix );
        void PBMVM( const double *x, double *b );
        // bool GaussElim( const double *rhs, double *x, PseudoBlockMatrix tmp );
        bool GaussSeidel( const double *rhs, const double *x0, double resTol, double convTol, double *x );


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

        int PseudoBlockMatrix::binSearchBlock( int rowIdx ){
            int lo = 0, hi = Nblocks-1,mid=-1;
            while(lo != hi)
            {
                mid = (lo+hi+1)/2;
                if(rowStart[mid]<=rowIdx){
                    lo = mid;
                }
                else{
                    hi = mid-1;
                }
            }
            return mid;
        }


        double PseudoBlockMatrix::getOverBlocks( int i, int j ){
            assert(i<N && i>=0 && j<M && j>=0);
            int tblock = binSearchBlock(i);
            return getFromBlock(tblock,i-rowStart[tblock],j-colStart[tblock]);
        }

        double PseudoBlockMatrix::getFromBlock( int b, int i, int j ){
            assert(i<blockN[b] && j<blockM[b] && i>=0 && j>=0);
            return blocks[b][i][j];
        }
        
        double PseudoBlockMatrix::getOverDense( int i, int j ){
            assert(i<N && i>=0 && j<M && j>=0);
            // Last row or corner!
            if(i==N-1){
                // Corner!
                if(j==M-1){
                    return denseCorner;
                }
                // Dense row!
                return denseRow[j];
            }
            // Last col!
            if(j==M-1){
                return denseCol[i];
            }
            return NULL;
        }

        double PseudoBlockMatrix::get( int i , int j )
        {
            assert(i<N && i>=0 && j<M && j>=0);
            double tmp = getOverDense(i,j);
            if(tmp!=NULL)
                return tmp;
            return getOverBlocks(i,j);
        }
        bool PseudoBlockMatrix::setOverBlocks( int i, int j, double val ){
            assert(i<N && i>=0 && j<M && j>=0);
            int tblock = binSearchBlock(i);
            return setFromBlock(tblock,i-rowStart[tblock],j-colStart[tblock],val);
        }

        bool PseudoBlockMatrix::setFromBlock( int b, int i, int j, double val ){
            assert(i<blockN[b] && j<blockM[b] && i>=0 && j>=0);
            blocks[b][i][j] = val;
            return true;
        }
        
        bool PseudoBlockMatrix::setOverDense( int i, int j, double val ){
            assert(i<N && i>=0 && j<M && j>=0);
            // Last row or corner!
            if(i==N-1){
                // Corner!
                if(j==M-1){
                    denseCorner = val;
                    return true;
                }
                // Dense row!
                denseRow[j] = val;
                return true;
            }
            // Last col!
            if(j==M-1){
                denseCol[i] = val;
                return true;
            }
            return false;
        }

        bool PseudoBlockMatrix::set( int i , int j, double val )
        {
            assert(i<N && i>=0 && j<M && j>=0);
            if(setOverDense(i,j,val))
                return true;
            return setOverBlocks(i,j,val);
        }
        bool PseudoBlockMatrix::addOverBlocks( int i, int j, double val ){
            assert(i<N && i>=0 && j<M && j>=0);
            int tblock = binSearchBlock(i);
            return addFromBlock(tblock,i-rowStart[tblock],j-colStart[tblock],val);
        }

        bool PseudoBlockMatrix::addFromBlock( int b, int i, int j, double val ){
            assert(i<blockN[b] && j<blockM[b] && i>=0 && j>=0);
            blocks[b][i][j] += val;
            return true;
        }
        
        bool PseudoBlockMatrix::addOverDense( int i, int j, double val ){
            assert(i<N && i>=0 && j<M && j>=0);
            // Last row or corner!
            if(i==N-1){
                // Corner!
                if(j==M-1){
                    denseCorner += val;
                    return true;
                }
                // Dense row!
                denseRow[j] += val;
                return true;
            }
            // Last col!
            if(j==M-1){
                denseCol[i] += val;
                return true;
            }
            return false;
        }

        bool PseudoBlockMatrix::add( int i , int j, double val )
        {
            assert(i<N && i>=0 && j<M && j>=0);
            if(addOverDense(i,j,val))
                return true;
            return addOverBlocks(i,j,val);
        }

}

#endif