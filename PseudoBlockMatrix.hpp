#ifndef _PSEUDOBLOCKMATRIX_
#define _PSEUDOBLOCKMATRIX_

#include <iostream>
#include <cstddef>
#include <assert.h>
#include <math.h>

namespace PSEUDOBLOCKMATRIX
{
    enum PBMState { undefined, initialized, allocated };

    class PseudoBlockMatrix
    {
    private:
        PBMState _state = undefined;
        double ***blocks; // 2D dense blocks
        double *denseRow; // 1D dense final row
        double *denseCol; // 1D dense final column
        double denseCorner; // 0D final corner value
        double *rowStart, *colStart; // Indices of first row/col per block. For binsearch of get/add/set

    public:
        int N,M; // Number of rows and cols
        int *blockN, *blockM; // Number of rows and cols per block
        int Nblocks; // Number of blocks

        void allocate(int tNblocks, int* tblockN, int* tblockM);
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

        PseudoBlockMatrix()
        {
            _state = initialized;
        }

        PseudoBlockMatrix(
            int tNblocks,
            int* tblockN,
            int* tblockM) : PseudoBlockMatrix()
        {
            allocate( tNblocks, tblockN, tblockM);
        }

        ~PseudoBlockMatrix()
        {
            if(_state==allocated){
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
        }

        


    };

    void PseudoBlockMatrix::allocate(
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
        
        _state= allocated;
    }

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
        if(!(i<blockN[b] && j<blockM[b] && i>=0 && j>=0))
            return NULL;
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

    void PseudoBlockMatrix::fill( double val )
    {
        for(int b=0;b<Nblocks;b++){
            for(int i=0;i<blockN[b];i++){
                for(int j=0;j<blockM[b];j++){
                    blocks[b][i][j] = val;
                }
            }
        }
        for(int i=0;i<N-1;i++)
            denseCol[i] = val;
        for(int j=0;j<M-1;j++)
            denseRow[j] = val;
        denseCorner = val;
    }
    PseudoBlockMatrix PseudoBlockMatrix::shallowCopy( PseudoBlockMatrix )
    {
        return PseudoBlockMatrix(Nblocks,blockN,blockM);
    }

    void PseudoBlockMatrix::PBMVM( const double *x, double *B )
    {
        int rs,cs,b,i,j;
        for(i=0;i<N-1;i++)
            B[i] = 0;
        for(b=0;b<Nblocks;b++){
            rs = rowStart[b];
            cs = colStart[b];
            for(i=0;i<blockN[b];i++){
                for(j=0;j<blockM[b];j++){
                    B[i+rs] += blocks[b][i][j]*x[j+cs];
                }
            }
        }
        for(i=0;i<N-1;i++)
            B[i] += denseCol[i]*x[M-1];
        for(j=0;j<N-1;j++)
            B[N-1] += denseRow[j]*x[j];
        B[N-1] += denseCorner*x[M-1];
    }
    // bool GaussElim( const double *rhs, double *x, PseudoBlockMatrix tmp );
    bool PseudoBlockMatrix::GaussSeidel( const double *rhs, const double *x0, double resTol, double convTol, double *x )
    {
        // Iterators
        int rs,cs,b,i,j,row,col,iter = 0,maxIter = 100;
        // Tmp
        double tmpx,tmpAii,tmpres;

        // Assert square!
        assert(N==M);

        // Assert GS conditions (no zero entries on diag).
        //  No pivoting. Mostly bc permutation would slow computation and is not necessary for body model.
        for(i=0;i<N;i++){
            int tmp = get(i,i);
            assert(tmp!=NULL);
            assert(tmp!=0);
        }

        // Copy x0 to x
        for(j=0;j<M-1;j++)
            x[j] = x0[j];

        // Residuals and convergences
        double res = 2*resTol;
        double conv = 2*convTol;
        while(res>resTol && conv>convTol && iter<maxIter){
            res = 0;
            conv = 0;
            ++iter;

            // Iterate over block rows
            for(b=0;b<Nblocks;b++){
                rs = rowStart[b];
                cs = colStart[b];
                
                for(i=0;i<blockN[b];i++){
                    row = rs+i;
                    // Temporary value of x
                    tmpx = 0;
                    // Add rhs
                    tmpx += rhs[row];
                    // Subtract last col*last x
                    tmpx -= denseCol[row]*x[M-1];
                    // Subtract jth col * jth x
                    for(j=0;j<blockM[b];j++){
                        tmpx -= blocks[b][i][j]*x[j+cs];
                    }
                    // Add back Aii * xii
                    tmpAii = getOverBlocks(row,row);
                    tmpx += tmpAii*x[row];
                    // Divide by Aii
                    tmpx /= tmpAii;
                    // compute conv
                    conv = max(conv,abs(x[row]-tmpx));
                    // Update x
                    x[row] = tmpx;
                }
            }
            // Iterate over last row
            // Temporary value of x
            tmpx = 0;
            // Add rhs
            tmpx += rhs[N-1];
            // Subtract jth col * jth x
            for(j=0;j<M-1;j++){
                tmpx -= denseRow[j]*x[j];
            }
            // Divide by corner
            tmpx /= denseCorner;
            // compute conv
            conv = max(conv,abs(x[M-1]-tmpx));
            // Update x
            x[M-1] = tmpx;


            // Check residual
            res = 0;


        }
        return iter<maxIter && res<resTol && conv<convTol;
    }
}

#endif