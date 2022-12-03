#include "Simulator.cpp"
int main(void){
    double *args = new double[4];
    args[0] = 300;
    args[1] = 305;
    args[2] = 300;
    args[3] = 300;

    double *outs = new double[1];

    BODYMODEL::BodyModel* body = BODYMODEL::defaultBody();
    SIMMODEL::SimModel* sim = new SIMMODEL::InjuryCase();

    SIMULATOR::Simulator* simulator = new SIMULATOR::Simulator(body,sim);
    simulator->initializer();
    simulator->runSim(args,outs);
}