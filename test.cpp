#include "Simulator.cpp"
#ifdef _WIN32
#include <Windows.h>
#else
#include <unistd.h>
#endif
int main(void){
    double *args = new double[4];
    args[0] = 300;
    args[1] = 305;
    args[2] = 300;
    args[3] = 300;

    double *outs = new double[1];
    
    //BODYMODEL::BodyModel* body = new BODYMODEL::BodyModel(10);
    int nb = 3;
    int ns[3] = {1, 2, 1};
    int ms[3] = {1, 2, 1};

    PSEUDOBLOCKMATRIX::PseudoBlockMatrix* A = new PSEUDOBLOCKMATRIX::PseudoBlockMatrix(nb,ns,ms);
    A->fill(0);
    A->set(1,1,7);
    A->set(2,1,3);
    A->set(2,2,4);
    A->set(3,3,.5);
    A->set(0,0,.2);
    A->set(4,4,.6);
    A->set(4,3,-.3);
    A->set(2,4,-.7);
    
    double rhs[4] = {1,1,1,1};
    double x0[4] = {1,1,1,1};
    double x[4];
    A->Print();
    A->directsolve(rhs,x);
    //A->GaussSeidel(rhs,x0,1e-3,1e-3,x);
    cout << x[0] << endl;
    cout << x[1] << endl;
    cout << x[2] << endl;
    cout << x[3] << endl;
    cout << x[4] << endl;
    /*
    BODYMODEL::BodyModel* body = BODYMODEL::defaultBody();
    SIMMODEL::SimModel* sim = new SIMMODEL::KonstasCase();

    SIMULATOR::Simulator* simulator = new SIMULATOR::Simulator(body,sim);
    simulator->initializer();
    simulator->runSim(args,outs);
    
    cout << outs[0] << endl;
    
    delete simulator;
    delete sim;
    delete body;
    delete [] args;
    delete [] outs;
 */   
}