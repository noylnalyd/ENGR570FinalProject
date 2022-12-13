#include "Simulator.cpp"
#ifdef _WIN32
#include <Windows.h>
#else
#include <unistd.h>
#endif
#include <chrono>
int main(void){
    double *args = new double[4];
    args[0] = 20+273.15;
    args[1] = 20+273.15;
    args[2] = 300;
    args[3] = 7+273.15;

    double *outs = new double[1];
    
    // int nb = 3;
    // int ns[3] = {1, 2, 1};
    // int ms[3] = {1, 2, 1};


    // PSEUDOBLOCKMATRIX::PseudoBlockMatrix* A = new PSEUDOBLOCKMATRIX::PseudoBlockMatrix(nb,ns,ms);
    // A->fill(0);
    // A->set(1,1,7);
    // A->set(2,1,3);
    // A->set(2,2,4);
    // A->set(3,3,.5);
    // A->set(0,0,.2);
    // A->set(4,4,.6);
    // A->set(4,3,-.3);
    // A->set(2,4,-.7);
    

    // double rhs[5] = {1,1,1,1,1};
    // double x0[5] = {1,1,1,1,1};
    // double x[5];
    
   
    // int nb = 10;
    // int ns[10] = {1+2*12, 1+1*6, 1+2*10,1+1*6,1+3*12,1+3*12,1+2*15,1+2*8,1+3*18,1+2*10};
    // int ms[10];
    // for(int i=0;i<10;i++)
    //     ms[i] = ns[i];


    // PSEUDOBLOCKMATRIX::PseudoBlockMatrix* A = new PSEUDOBLOCKMATRIX::PseudoBlockMatrix(nb,ns,ms);
    // double rhs[310] = {1,1,1,1,1};
    // double x0[310] = {1,1,1,1,1};
    // double x[310];
    // for(int i=0;i<310;i++){
    //     rhs[i] = 1;
    //     x0[i] = 1;
    // }
    // A->fillRandDD(1000);
    // // time_t begin, end;
    


    // A->Print();
    // using std::chrono::high_resolution_clock;
    // using std::chrono::duration_cast;
    // using std::chrono::duration;
    // using std::chrono::milliseconds;

    
    // auto t1 = high_resolution_clock::now();
    // using namespace std::chrono_literals;
    // //A->directsolve(rhs,x);
    // //A->GaussSeidel(rhs,x0,1e-5,1e-5,x);
    // A->denseGaussSeidel(rhs,x0,1e-5,1e-5,x);
    // auto t2 = high_resolution_clock::now();

    // /* Getting number of milliseconds as an integer. */
    // auto ms_int = duration_cast<milliseconds>(t2 - t1);

    // /* Getting number of milliseconds as a double. */
    // duration<double, std::milli> ms_double = t2 - t1;

    // cout << "Run time = " << ms_double.count() << " ms\n";
    // //A->GaussSeidel(rhs,x0,1e-3,1e-3,x);
    // cout << x[0] << endl;
    // cout << x[1] << endl;
    // cout << x[2] << endl;
    // cout << x[3] << endl;
    // cout << x[4] << endl;
    

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
}