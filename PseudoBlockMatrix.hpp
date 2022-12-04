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
        double*** blocks; // 2D dense blocks
        double *denseRow; // 1D dense final row
        double *denseCol; // 1D dense final column
        double denseCorner; // 0D final corner value
        int *rowStart, *colStart; // Indices of first row/col per block. For binsearch of get/add/set

    public:
        int N,M; // Number of rows and cols
        int *blockN, *blockM; // Number of rows and cols per block
        int *blockPicker; // Finds block given row address
        int Nblocks; // Number of blocks

        void Print();
        void allocate(int tNblocks, int* tblockN, int* tblockM);
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
        double checkRes( const double *x, const double *B );
        PseudoBlockMatrix shallowCopy( PseudoBlockMatrix );
        void PBMVM( const double *x, double *b );
        // bool GaussElim( const double *rhs, double *x, PseudoBlockMatrix tmp );
        bool directsolve( const double *rhs, double *x, PseudoBlockMatrix tmp );
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
            allocate( tNblocks, tblockN, tblockM );
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
        blocks = new double**[Nblocks]();
        rowStart = new int[Nblocks];
        colStart = new int[Nblocks];
        for(int b=0;b<Nblocks;b++){
            blocks[b] = new double*[blockN[b]];
            rowStart[b] = N;
            colStart[b] = M;

            for(int i=0;i<blockN[b];i++){
                blocks[b][i] = new double[blockM[b]];
            }
            N = N + blockN[b];
            M = M + blockM[b];
        }

        blockPicker = new int[N];
        for(int b=0;b<Nblocks;b++){
            for(int i=0;i<blockN[b];i++){
                blockPicker[rowStart[b]+i] = b;
            }
        }
        // Adjust for dense row/col!
        N = N+1;
        M = M+1;
        // Allocate dense col, row, corner
        denseCol = new double[N-1];
        denseRow = new double[M-1];
        denseCorner = 0;        
        
        fill(NAN);

        _state= allocated;
    }

    double PseudoBlockMatrix::getOverBlocks( int i, int j ){
        assert(i<N && i>=0 && j<M && j>=0);
        int tblock = blockPicker[i];
        return getFromBlock(tblock,i-rowStart[tblock],j-colStart[tblock]);
    }

    double PseudoBlockMatrix::getFromBlock( int b, int i, int j ){
        if(!(i<blockN[b] && j<blockM[b] && i>=0 && j>=0))
            assert(1==0);
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
        return NAN;
    }

    double PseudoBlockMatrix::get( int i , int j )
    {
        assert(i<N && i>=0 && j<M && j>=0);
        double tmp = getOverDense(i,j);
        if(!isnan(tmp))
            return tmp;
        return getOverBlocks(i,j);
    }
    bool PseudoBlockMatrix::setOverBlocks( int i, int j, double val ){
        assert(i<N && i>=0 && j<M && j>=0);
        int tblock = blockPicker[i];
        return setFromBlock(tblock,i-rowStart[tblock],j-colStart[tblock],val);
    }

    bool PseudoBlockMatrix::setFromBlock( int b, int i, int j, double val ){
        assert(i<blockN[b] && j<blockM[b] && i>=0 && j>=0);
        blocks[b][i][j] = val;
        return true;
    }
    
    bool PseudoBlockMatrix::setOverDense( int i, int j, double val ){
        assert(i<N && i>=0 && j<M && j>=0);
        assert(!isnan(val));
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
        assert(!isnan(val));
        if(setOverDense(i,j,val))
            return true;
        return setOverBlocks(i,j,val);
    }
    bool PseudoBlockMatrix::addOverBlocks( int i, int j, double val ){
        assert(i<N && i>=0 && j<M && j>=0);
        assert(!isnan(val));
        int tblock = blockPicker[i];
        return addFromBlock(tblock,i-rowStart[tblock],j-colStart[tblock],val);
    }

    bool PseudoBlockMatrix::addFromBlock( int b, int i, int j, double val ){
        assert(i<blockN[b] && j<blockM[b] && i>=0 && j>=0);
        assert(!isnan(val));
        blocks[b][i][j] += val;
        return true;
    }
    
    bool PseudoBlockMatrix::addOverDense( int i, int j, double val ){
        assert(i<N && i>=0 && j<M && j>=0);
        assert(!isnan(val));
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
        assert(!isnan(val));
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

    double PseudoBlockMatrix::checkRes(const double *x, const double *B){
        double res = 0,tres=0;
        int rs,cs,b,i,j;
        for(b=0;b<Nblocks;b++){
            rs = rowStart[b];
            cs = colStart[b];
            for(i=0;i<blockN[b];i++){
                tres = -B[i+rs];
                for(j=0;j<blockM[b];j++){
                    tres += blocks[b][i][j]*x[j+cs];
                }
                tres += denseCol[i]*x[M-1];
                res += tres*tres;
            }
        }
        tres = -B[N-1];
        for(j=0;j<M-1;j++)
            tres += denseRow[j]*x[j];
        tres += denseCorner*x[M-1];
        res += tres*tres;
        return sqrt(res);
    }
    void PseudoBlockMatrix::Print(){
        int rs,cs,b,i,j;
        for(b=1;b<2;b++){
            rs = rowStart[b];
            cs = colStart[b];
            cout << "B "<<b <<endl;
            for(i=0;i<blockN[b];i++){
                for(j=0;j<blockM[b];j++){
                    cout<<blocks[b][i][j] << "\t";
                }
                cout << endl;
            }
        }

        for(i=rowStart[1];i<rowStart[2];i++)
            cout<<denseCol[i] << endl;
        cout << endl << endl;
        for(j=rowStart[1];j<rowStart[2];j++)
            cout << denseRow[j] << endl;
        cout << endl;
        cout << denseCorner << endl;
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
    bool PseudoBlockMatrix::directsolve( double *rhs, const double *x0, double *x )
    {
        // Iterators
        int rs,cs,b,i,j,row,row2,col;
        // Tmp
        double res=NAN,norm=NAN;

        // Assert square!
        assert(N==M);

        // Assert GS conditions (no zero entries on diag).
        //  No pivoting. Mostly bc permutation would slow computation and is not necessary for body model.
        for(i=0;i<N;i++){
            double tmp = get(i,i);
            assert(tmp!=NAN);
            assert(tmp!=0);
        }
        for(i=0;i<N;i++){
            double tmp = rhs[i];
            assert(tmp!=NAN);
        }
        for(i=0;i<M;i++){
            double tmp = x0[i];
            assert(tmp!=NAN);
        }
        // Precondition by norming values
        for(b=0;b<Nblocks;b++){
            rs = rowStart[b];
            cs = colStart[b];
            
            for(i=0;i<blockN[b];i++){
                row = rs+i;

                norm = 0;
                // Subtract last col*last x
                norm += abs(denseCol[row]);
                // Subtract jth col * jth x
                for(j=0;j<blockM[b];j++){
                    norm += abs(blocks[b][i][j]);
                }
                denseCol[row] /= norm;
                for(j=0;j<blockM[b];j++){
                    blocks[b][i][j]/=norm;
                }
                rhs[row]/=norm;
            }
        }
        norm = abs(denseCorner);
        for(j=0;j<M-1;j++){
            norm += abs(denseRow[j]);
        }
        // Divide by corner
        denseCorner /= norm;
        for(j=0;j<M-1;j++){
            denseRow[j] /= norm;
        }

        // Copy x0 to x
        for(j=0;j<M;j++)
            x[j] = x0[j];
        double A[N][M];
        for(b=0;b<Nblocks;b++){
            rs = rowStart[b];
            cs = colStart[b];
            
            for(i=0;i<blockN[b];i++){
                row = rs+i;
                for(j=0;j<blockM[b];j++){
                    col = cs+j;
                    A[row][col] = blocks[b][i][j]*x[j+cs];
                }
                A[row][M-1] = denseCol[row];
            }
        }
        for(j=0;j<M-1;j++){
            A[N-1][j] = denseRow[j];
        }
        A[N-1][M-1] = denseCorner;
        for(row=0;row<N;row++){
            double Aiiinv = 1.0/A[row][row];
            for(row2=row+1;row2<N;row2++){
                double srat = Aiiinv*A[row2][row];
                A[row2][row] = 0;
                for(col=row+1;col<M;col++){
                    A[row2][col] -= srat*A[row][col];
                }
                rhs[row2] -= srat*rhs[row];
            }
        }
        for(row=N-1;row>=0;row--){
            x[row] = rhs[row]/A[row][row];
            for(col=row+1;col<M;col++){
                x[row] -= A[row][col]/A[row][row];
            }
        }
        res = checkRes(x,rhs);
        cout << res << endl;
        return res<1e-1;
    }
    // bool GaussElim( const double *rhs, double *x, PseudoBlockMatrix tmp );
    bool PseudoBlockMatrix::GaussSeidel( const double *rhs, const double *x0, double resTol, double convTol, double *x )
    {
        // Iterators
        int rs,cs,b,i,j,row,col,iter = 0,maxIter = 100;
        // Tmp
        double tmpx=NAN,tmpAii=NAN,tmpres=NAN,norm=NAN;

        // Assert square!
        assert(N==M);

        // Assert GS conditions (no zero entries on diag).
        //  No pivoting. Mostly bc permutation would slow computation and is not necessary for body model.
        for(i=0;i<N;i++){
            double tmp = get(i,i);
            assert(tmp!=NAN);
            assert(tmp!=0);
        }
        for(i=0;i<N;i++){
            double tmp = rhs[i];
            assert(tmp!=NAN);
        }
        for(i=0;i<M;i++){
            double tmp = x0[i];
            assert(tmp!=NAN);
        }
        // Precondition by norming values
        for(b=0;b<Nblocks;b++){
            rs = rowStart[b];
            cs = colStart[b];
            
            for(i=0;i<blockN[b];i++){
                row = rs+i;

                norm = 0;
                // Subtract last col*last x
                norm += abs(denseCol[row]);
                // Subtract jth col * jth x
                for(j=0;j<blockM[b];j++){
                    norm += abs(blocks[b][i][j]);
                }
                denseCol[row] /= norm;
                for(j=0;j<blockM[b];j++){
                    blocks[b][i][j]/=norm;
                }
                rhs[row]/=norm;
            }
        }
        // Iterate over last row
        // Temporary value of x
        norm = abs(denseCorner);
        for(j=0;j<M-1;j++){
            norm += abs(denseRow[j]);
        }
        // Divide by corner
        denseCorner /= norm;
        for(j=0;j<M-1;j++){
            denseRow[j] /= norm;
        }

        // Copy x0 to x
        for(j=0;j<M;j++)
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
                    tmpAii = blocks[b][i][i];
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
            res = checkRes(x,rhs);
            cout << res << endl;
        }
        assert(iter<maxIter-1);
        return iter<maxIter && res<resTol && conv<convTol;
    }
}

#endif