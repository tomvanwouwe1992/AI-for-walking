#include <OpenSim/OpenSim.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <ctime>
#include <windows.h>
#include "matrix.h"
#include "math.h"
#include <tuple>
#include "mex.h"
#include "SimTKcommon.h"

#include "SimTKmath.h"
#include "simmath/TimeStepper.h"
using namespace OpenSim;
using namespace SimTK;
using namespace std;

int find(Vector arr, double lengte, double seek)
{   int return_value;
    for (int i = 0; i < lengte; ++i)
    {
        if (arr[i] == seek) {return_value = i;}
    }
    return return_value;
}

Vector get_states(Model* osimModel, State& s) {
    
    const Vector& observation_original = s.getY();
    Vector observation(31, 0.0);
    
    observation[1] = observation_original.get(0);
    observation[2] = observation_original.get(1);
    observation[3] = observation_original.get(2);
    
    observation[4] = observation_original.get(10);
    observation[5] = observation_original.get(11);
    observation[6] = observation_original.get(12);
    
    observation[7]  = observation_original.get(3);
    observation[8]  = observation_original.get(4);
    observation[9]  = observation_original.get(5);
    observation[10] = observation_original.get(6);
    observation[11] = observation_original.get(7);
    observation[12] = observation_original.get(8);
    
    observation[13] = observation_original.get(13);
    observation[14] = observation_original.get(14);
    observation[15] = observation_original.get(15);
    observation[16] = observation_original.get(16);
    observation[17] = observation_original.get(17);
    observation[18] = observation_original.get(18);
    
    Vec3 massCenterPos = osimModel->calcMassCenterPosition(s);
    Vec3 massCenterVel = osimModel->calcMassCenterVelocity(s);
    
    observation[19] = massCenterPos.get(0);
    observation[20] = massCenterPos.get(1);
    
    observation[21] = massCenterVel.get(0);
    observation[22] = massCenterVel.get(1);
    
    const SimbodyEngine& SBE = osimModel->getSimbodyEngine();               // Get simbody engine - subclass that contains methods to provide information/calculation on the model (e.g.: calculate distance between two bodies)
    const BodySet& BS = osimModel->getBodySet();
    Vec3 &PosHEAD = Vec3(0.0, 0.0, 0.0);
    Vec3 &RelPosHEAD = Vec3(0.0, 0.0, 0.0);
    SBE.getPosition(s, BS.get("head"), RelPosHEAD, PosHEAD);
    
    observation[23] = PosHEAD.get(0);
    observation[24] = PosHEAD.get(1);
    
    Vec3 &PosPELVIS = Vec3(0.0, 0.0, 0.0);
    SBE.getPosition(s, BS.get("pelvis"), RelPosHEAD, PosPELVIS);
    observation[25] = PosPELVIS.get(0);
    observation[26] = PosPELVIS.get(1);
    
    Vec3 &PosTOESR = Vec3(0.0, 0.0, 0.0);
    SBE.getPosition(s, BS.get("toes_r"), RelPosHEAD, PosTOESR);
    observation[27] = PosTOESR.get(0);
    observation[28] = PosTOESR.get(1);
    
    Vec3 &PosTOESL = Vec3(0.0, 0.0, 0.0);
    SBE.getPosition(s, BS.get("toes_l"), RelPosHEAD, PosTOESL);
    observation[29] = PosTOESL.get(0);
    observation[30] = PosTOESL.get(1);
    
    return observation;
}

    Vector NN_controller(Vector inputs,  std::tuple<Vector,Vector,Vector> nodegenes, std::tuple<Vector,Vector,Vector,Vector>  connectiongenes, Vector NN_info, Model* osimModel, State& s, Vector ContactInfo){
    // We simplify the observationvector (inputs) to a limited number of states
    Vector adapted_inputs(8, 0.0);
    
    const SimbodyEngine& SBE = osimModel->getSimbodyEngine();               // Get simbody engine - subclass that contains methods to provide information/calculation on the model (e.g.: calculate distance between two bodies)
    const BodySet& BS = osimModel->getBodySet();
    
    Vec3 &RelPos = Vec3(0.0, 0.0, 0.0);
    
    Vec3 &PosCALCNL = Vec3(0.0, 0.0, 0.0);
    SBE.getPosition(s, BS.get("calcn_l"), RelPos, PosCALCNL);
    
    Vec3 &PosCALCNR = Vec3(0.0, 0.0, 0.0);
    SBE.getPosition(s, BS.get("calcn_r"), RelPos, PosCALCNR);
    
    Vec3 &VelHead = Vec3(0.0, 0.0, 0.0);
    SBE.getVelocity(s, BS.get("head"), RelPos, VelHead);
        
    
    // Lengths of the legs, seen as one segment.
    adapted_inputs[0] =sqrt((inputs[1] - PosCALCNL[0]) * (inputs[1] - PosCALCNL[0]) +  (inputs[2] - PosCALCNL[1]) * (inputs[2] - PosCALCNL[1]));
    adapted_inputs[1] =sqrt((inputs[1] - PosCALCNR[0]) * (inputs[1] - PosCALCNR[0]) +  (inputs[2] - PosCALCNR[1]) * (inputs[2] - PosCALCNR[1]));
    // Angle of the legs, seen as one segment.
    adapted_inputs[2] = atan2(inputs[28] - inputs[2], inputs[27] - inputs[1]);
    adapted_inputs[3] = atan2(inputs[30] - inputs[2], inputs[29] - inputs[1]);
    
    // Pelvis speed
    adapted_inputs[4] = inputs[5]; // horizontal
    adapted_inputs[5] = inputs[6]; // vertical
    
    // Trunk angle and angular velocity
    double y = inputs[24] - inputs[3];
    double x = inputs[23] - inputs[2];
    double dy = VelHead[1] - inputs[6];
    double dx = VelHead[0] - inputs[5];
            
    adapted_inputs[6] = atan2(y , x);
    adapted_inputs[7] = (-y / (x * x + y * y)) * dx + (x / (x * x + y * y)) * dy; //Derivative of atan2 -->  [-y/(x²+y²)]*dx + [x/(x²+y²)]*dy
        
    //Ankle angle
    adapted_inputs[8] = inputs[9];
    adapted_inputs[9] = inputs[12];
    
    //Contact Inputs
    if (ContactInfo[0] == 0){
        adapted_inputs[10] = 0;
        adapted_inputs[11] = 0;
    }
    else{
        Vec3 COM = osimModel->calcMassCenterPosition(s);
        adapted_inputs[10] = ContactInfo[0] - COM[0];
        adapted_inputs[11] = ContactInfo[1] - COM[0];
    }
    adapted_inputs[12] = ContactInfo[2];
    adapted_inputs[13] = ContactInfo[3];
    
    //////////////////////////////////////////////////////////////////////
    double change_treshold = 0.001;
    int number_of_inputnodes = NN_info[0];
    int number_of_outputnodes = NN_info[1];
    
    int number_of_hiddennodes = NN_info[2];
    int number_of_connections = NN_info[3];
    int no_change_count=0;
    int index_loop=0;
    
    Vector output_signal(number_of_outputnodes, 0.0);
    
    Vector nodegenes_ID = get<0>(nodegenes);
    Vector nodegenes_input = get<1>(nodegenes);
    Vector nodegenes_output = get<2>(nodegenes);
    Vector connectiongenes_from = get<0>(connectiongenes);
    Vector connectiongenes_to = get<1>(connectiongenes);
    Vector connectiongenes_weight = get<2>(connectiongenes);
    Vector connectiongenes_enabled = get<3>(connectiongenes);
    
    
    for (int i = number_of_inputnodes + 1; i < number_of_inputnodes + number_of_outputnodes + number_of_hiddennodes + 1; i++){nodegenes_input[i] = 0;}
    nodegenes_input[number_of_inputnodes] = 1;
    
    for (int i = 0; i < number_of_inputnodes; i++){nodegenes_input[i] = adapted_inputs[i];}
    
    for(int i = 0; i < number_of_inputnodes + 1; i++){nodegenes_output[i] = nodegenes_input[i];}
    
    for (int i = number_of_inputnodes + 1; i < number_of_inputnodes + number_of_outputnodes + number_of_hiddennodes + 1; i++){nodegenes_output[i] = 1 / (1 + exp(-4.9 * nodegenes_input[i]));}
    
    Vector vector_node_state = nodegenes_output;
    double ID_connection_from_node;
    double ID_connection_to_node;
    double connection_weight;
    int index_connection_from_node;
    int index_connection_to_node;
    
    
    int checkchanges = number_of_inputnodes + number_of_outputnodes + number_of_hiddennodes + 1;
    
    while ((no_change_count<checkchanges) && index_loop<3*number_of_connections ) {
        
        index_loop++;
        vector_node_state = nodegenes_output;
        
        for (int i = 0; i < number_of_connections; i++){
            
            ID_connection_from_node = connectiongenes_from[i];
            ID_connection_to_node = connectiongenes_to[i];
            connection_weight = connectiongenes_weight[i];
            
            index_connection_from_node = find(nodegenes_ID, number_of_inputnodes+number_of_outputnodes+number_of_hiddennodes+1, ID_connection_from_node);
            
            index_connection_to_node   = find(nodegenes_ID, number_of_inputnodes+number_of_outputnodes+number_of_hiddennodes+1, ID_connection_to_node);
 

            if (connectiongenes_enabled[i] == 1){
                nodegenes_input[index_connection_to_node] = nodegenes_input[index_connection_to_node] + nodegenes_output[index_connection_from_node] * connection_weight;
            }
        }
        for (int j = number_of_inputnodes + 1; j < number_of_inputnodes+ 1 + number_of_hiddennodes + number_of_outputnodes; j++){ nodegenes_output[j] = 1 / (1 + exp(-4.9 * nodegenes_input[j]));}
        
        for (int j = number_of_inputnodes + 1; j < number_of_inputnodes+ 1 + number_of_hiddennodes + number_of_outputnodes; j++){nodegenes_input[j] = 0;}
        
        no_change_count = 0;
        for (int j = 0; j < number_of_inputnodes+ 1 + number_of_hiddennodes + number_of_outputnodes; j++)
        {
            if (nodegenes_output[j] - vector_node_state[j] < change_treshold){no_change_count++;}
        }
        
    }
    
    for (int i = 0; i < number_of_outputnodes; i++){output_signal[i] = nodegenes_output[number_of_inputnodes + 1 + i];}
    
    return output_signal;
}
    
Vector Calculate_ContactInfo(Vector Right_Foot_Contact){
    
    double F_x = Right_Foot_Contact[0];
    
    double F_y = Right_Foot_Contact[1];
    
    double F_z = Right_Foot_Contact[2];
   
    double M_x = Right_Foot_Contact[3];
    
    double M_y = Right_Foot_Contact[4]; 
    
    double M_z = Right_Foot_Contact[5];
    
    
    double COP_y = 0; //(M_y * F_y - M_z * F_z - M_x*F_x) / (2*F_y*F_x);
    double COP_x = (M_z - COP_y * F_x) / F_y;
    double COP_z = (M_y - COP_x * F_z) / F_x;
    
        
    Vector COP(3, 0.0);
    COP[0] = COP_x;
    COP[1] = COP_y;
    COP[2] = COP_z;
    
    return COP;
}

tuple<double,bool> Penalty_Calc(Vector observation, Model* osimModel, State& s) {
    
    
//     const SimTK::SimbodyMatterSubsystem& SMS = osimModel->getMatterSubsystem();
//     SpatialVec momentum = SMS.calcSystemCentralMomentum(s);
//     Vec3 angularmomentum = momentum[0];
//     
//     double momentum_Z = angularmomentum[2]*angularmomentum[2];
     double penalty = 0.0;
    bool stop_integration;
    stop_integration = false;
    if (observation[3] < 0.70) {stop_integration = true;mexPrintf("Pelvis low");}
    if (observation[5] < -0.05) {stop_integration = true;mexPrintf("Velocity negative");}
   // if (observation[26] < -0.06){stop_integration = true;}
   // if (observation[28] < -0.06){stop_integration = true;}
    
    
    std::tuple <double,bool> penalty_break;
    penalty_break = std::make_tuple(penalty,stop_integration);
    
    return penalty_break;
}
/////////////////////////////////////////

class MyIntegrationStopper : public TriggeredEventHandler {
public:
    MyIntegrationStopper(Model* osimModel) : TriggeredEventHandler(Stage::Dynamics),model(osimModel) {
    }
    Real SimTK::TriggeredEventHandler::getValue	(const State& s)	const{
    // model->getMultibodySystem().realize(state, Stage::Dynamics);
    Real return_value = 1;   
    const SimbodyEngine& SBE = model->getSimbodyEngine();               // Get simbody engine - subclass that contains methods to provide information/calculation on the model (e.g.: calculate distance between two bodies)
    const BodySet& BS = model->getBodySet();
    Vec3 &RelPosHEAD = Vec3(0.03, 0.04, 0.0);
    Vec3 &PosTOESL = Vec3(0.0, 0.0, 0.0);
    SBE.getPosition(s, BS.get("toes_l"), RelPosHEAD, PosTOESL);
    
    if (PosTOESL.get(1) < -0.00){return_value = 0;}
    
    Vec3 &PosTOESR = Vec3(0.0, 0.0, 0.0);
    SBE.getPosition(s, BS.get("toes_r"), RelPosHEAD, PosTOESR);
   
    double sphere_R = PosTOESR.get(1);
    if (PosTOESR.get(1) < -0.00){return_value = 0;}
    
     return return_value;
//         mexPrintf("check");
//         Storage forces = Freporter->getForceStorage();
//         const Real& int_time = state.getTime();
//         Real return_value = 1.0;
//         double contact_forceR;
//         double contact_forceL;
//         int size = forces.getSize()-1;
//         mexPrintf("%i", size );
//         forces.getData(forces.getSize()-1, 19, &contact_forceR);
//         forces.getData(forces.getSize()-1, 37, &contact_forceL);
//         mexPrintf("check3");
//         mexPrintf("%f",contact_forceL );
//         mexPrintf("%f",contact_forceR );
//         if (abs(contact_forceR) > 3500 || abs(contact_forceL > 3500)){return_value = 0.0;}
//         return return_value;
        }
    void handleEvent(State& state, Real accuracy, bool& shouldTerminate) const {
        mexPrintf("HUNTCROSSLEY");
        shouldTerminate = true;
    }
    
private:
    Model* model;
    
};


////////////////////////

class MyPeriodicController : public PeriodicEventHandler {
public:
    MyPeriodicController(Real eventinterval, Model* osimModel, tuple<tuple<Vector,Vector,Vector>,tuple<Vector,Vector,Vector,Vector>,Vector> NN_STRUCTURE, double* cost_function, ForceReporter* freporter) : PeriodicEventHandler(eventinterval), model(osimModel), NN_structure(NN_STRUCTURE), cost(cost_function), force_reporter(freporter){
    }
    void handleEvent(State& state, Real accuracy, bool& shouldTerminate) const {
        
        Vector states = get_states(model,  state);
        const Real& int_time = state.getTime();
        
		Vector forces(67, 0.0);
        force_reporter->getForceStorage().getDataAtTime(force_reporter->getForceStorage().getLastTime(),42,forces);
        
        
        Vector Right_Foot_Contact(6, 0.0);
        for (int l = 18; l < 24; l++ ){
            Right_Foot_Contact[l-18] = forces[l];
        }
        Vector Left_Foot_Contact(6, 0.0);
        for (int k = 36; k < 42; k++){
            Left_Foot_Contact[k-36] = forces[k];
        }
        double contact_Rx;
        double contact_Lx;
        if (Right_Foot_Contact[0] + Right_Foot_Contact[1] + Right_Foot_Contact[2] == 0){
            contact_Rx = 0;
            contact_Lx = 0;
        }
        else{
        //const Storage* force_storage = &(force_reporter->getForceStorage());
        //force_storage->print("GaitProblem_forces.mot");
        Vector contact_R = Calculate_ContactInfo(Right_Foot_Contact);
        contact_Rx = contact_R[0];
        //cout << contact_R << "   ";
        Vector contact_L = Calculate_ContactInfo(Left_Foot_Contact);
        contact_Lx = contact_L[0];
        //cout << contact_L << endl;
        }
        double Fy_R = Right_Foot_Contact[1];
        double Fy_L = Left_Foot_Contact[1];
        
        Vector ContactInfo(4, 0.0);   // contains [COPx right foot, COPx left foot, Fy right foot, Fy left foot]
        ContactInfo[0] = contact_Rx;
        ContactInfo[0] = contact_Lx;
        ContactInfo[0] = Fy_R;
        ContactInfo[0] = Fy_L;
                
        std::tuple<double, bool> result_penalty = Penalty_Calc(states, model, state);
        shouldTerminate = get<1>(result_penalty);
        //cost[0] = cost[0] + get<0>(result_penalty);
               
        Vector activation_vector = NN_controller(states, get<0>(NN_structure), get<1>(NN_structure), get<2>(NN_structure), model, state, ContactInfo);
            
        Vector& activations = state.updZ();
        
        for (int i = 0; i < 18; i++){
            if (activation_vector[i] > 0.99)  { activation_vector[i] = 0.99;}
            if (activation_vector[i] < 0.01)  { activation_vector[i] = 0.01;}
            
            activations[2*i] = activation_vector[i];}
    }
private:
    Model* model;
	ForceReporter* force_reporter;
    tuple<tuple<Vector,Vector,Vector>,tuple<Vector,Vector,Vector,Vector>,Vector> NN_structure;
    double* cost;
};

void main(double* TimeVec, Model* osimModel, double* cost_function, std::tuple<Vector,Vector,Vector> nodegenes,std::tuple<Vector,Vector,Vector,Vector> connectiongenes, Vector NN_info)
{
    std::streambuf* cout_sbuf = std::cout.rdbuf(); // save original sbuf
    std::ofstream   fout("cout.txt");
    std::cout.rdbuf(fout.rdbuf()); // redirect 'cout' to a 'fout'
    
    std::streambuf* cerr_sbuf = std::cerr.rdbuf(); // save original sbuf
    std::ofstream   ferr("cerr.txt");
    std::cerr.rdbuf(ferr.rdbuf()); // redirect 'cout' to a 'fout'
    cost_function[0] = 0.0;
    
    std::tuple<tuple<Vector,Vector,Vector>,tuple<Vector,Vector,Vector,Vector>,Vector> NN_structure = make_tuple(nodegenes, connectiongenes, NN_info);
    // ******************** SimBodyLand***********************//
    osimModel->buildSystem();
    Real eventinterval = 0.01;
    ForceReporter* reporter = new ForceReporter(osimModel);
    osimModel->addAnalysis(reporter);
    
    
    MultibodySystem& system_ref = osimModel->updMultibodySystem();
    
    MyPeriodicController* myPeriodicController = new MyPeriodicController(eventinterval, osimModel, NN_structure, cost_function, reporter);
    
    MyIntegrationStopper* myIntegrationStopper = new MyIntegrationStopper(osimModel);
    
    system_ref.addEventHandler(myPeriodicController);
    //system_ref.addEventHandler(myIntegrationStopper);
    
    // ******************** OpenSimLand***********************//
    State& s = osimModel->initializeState();
    
    double initialTime;
    double finalTime;
    initialTime = TimeVec[0]; finalTime = TimeVec[1];
    
    osimModel->equilibrateMuscles(s);
    osimModel->getMultibodySystem().realize(s, Stage::Dynamics);
    
   //  BodyKinematics* BodyKinematics_reporter = new BodyKinematics(osimModel, true);
   // BodyKinematics_reporter->setModel(*osimModel);
    //BodyKinematics_reporter->setRecordCenterOfMass(true);
    //osimModel->addAnalysis(BodyKinematics_reporter);
    
    RungeKuttaMersonIntegrator integrator(osimModel->getMultibodySystem());     //Initialize integrator (RungeKuttaMerson-type)
    integrator.setAccuracy(1e-6);
    Manager manager(*osimModel, integrator);
    manager.setInitialTime(initialTime);                                        //Set initial and final times of the manager
    manager.setFinalTime(finalTime);
    manager.integrate(s);
    Storage output_storage = manager.getStateStorage();                     //Store&Output of the states as .sto-file
    Array<double> end_state = output_storage.getLastStateVector()->getData();
    output_storage.print("output_1.sto");   
    Vec3 massCenterPos = osimModel->calcMassCenterPosition(s);
    
    const Storage* force_storage = &(reporter->getForceStorage());

    
//     Array<double> X_com_LToes; Array<double>& X_com_LToes_ref = X_com_LToes; 
//     Body_position->getDataColumn("toes_l_X",X_com_LToes_ref, 0.0);
    
    //BodyKinematics_reporter->printResults("BodyKinematics");
    //BodyK_storage->print("BodyKinematics.mot");
    
    
   
    double cost = massCenterPos[0];//+0.1*end_state[0];
    double time = output_storage.getLastTime();
    //double angular_cost = - 0.02 * (cost_function[0] / (output_storage.getLastTime() / 0.01));
    mexPrintf("  DISTANCE WALKED  =  ");
    mexPrintf("%f",cost);
    mexPrintf("\n");
    
    mexPrintf("  TIME WALKED  =  ");
    mexPrintf("%f",time);
    mexPrintf("\n");
   
    //force_storage->print("GaitProblem_forces.mot");
    cost_function[0] = cost*sqrt(time);

    return;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //	OPEN Model
    
    char *osimFileName_ptr;
    osimFileName_ptr = mxArrayToString(prhs[0]);
    std::string osimFileName = std::string(&osimFileName_ptr[0]);
    
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //	IMPORT timevector
    double *TimeVec; double time_instants; const int *dimTime;
    TimeVec = mxGetPr(prhs[1]);
    dimTime = mxGetDimensions(prhs[1]);
    time_instants = *(dimTime + 1);
    
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //	IMPORT neural-network == CONTROLLER
    double *nodegenes_ptr;
    double *connectiongenes_ptr;
    const int *dim_nodegenes; const int *dim_connectiongenes;
    
    dim_nodegenes = mxGetDimensions(prhs[2]);
    int n_nodes = *(dim_nodegenes+1);
    int n_rows = *(dim_nodegenes+0);
    nodegenes_ptr = mxGetPr(prhs[2]);
    
    mexPrintf("%i",n_nodes);
   
    
    Vector IDs    = Vector(n_nodes, 0.0);
    Vector input  = Vector(n_nodes, 0.0);
    Vector output = Vector(n_nodes, 0.0);
    
    for (int j = 0; j < n_nodes; j++)
    {
        IDs[j] = nodegenes_ptr[j*4];
        input[j] = nodegenes_ptr[j*4+2 ];
        output[j] = nodegenes_ptr[j*4+3 ];
    }
    std::tuple<Vector,Vector,Vector> nodegenes = make_tuple(IDs, input, output);
    
    
    dim_connectiongenes = mxGetDimensions(prhs[3]);
    int n_connections = *(dim_connectiongenes+1);
    connectiongenes_ptr = mxGetPr(prhs[3]);
    
    Vector from = Vector(n_connections, 0.0);
    Vector to = Vector(n_connections, 0.0);
    Vector weight = Vector(n_connections, 0.0);
    Vector enabled = Vector(n_connections, 0.0);
    mexPrintf("%i",n_connections);
    for (int j = 0; j < n_connections; j++)
    {
        from[j]    = connectiongenes_ptr[j*5+1];
        to[j]      = connectiongenes_ptr[j*5 + 2];
        weight[j]  = connectiongenes_ptr[j*5 + 3];
        enabled[j] = connectiongenes_ptr[j*5 + 4];
    }
    std::tuple<Vector,Vector,Vector,Vector> connectiongenes = make_tuple(from, to, weight, enabled);
    
    double* NN_info_ptr;
    NN_info_ptr = mxGetPr(prhs[4]);
    Vector NN_info = Vector(4, 0.0);
    NN_info[0] = NN_info_ptr[0];
    NN_info[1] = NN_info_ptr[1];
    NN_info[2] = NN_info_ptr[2];
    NN_info[3] = NN_info_ptr[3];
    
    Model model = Model(osimFileName);
    Model *osimModel = &model;
    
    
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    double *cost_function;
    cost_function	= mxGetPr(plhs[0]);
    
    main(TimeVec, osimModel, cost_function, nodegenes, connectiongenes, NN_info);
    
    
}



//
    // COP calculation
//     Storage* Body_position = BodyKinematics_reporter->getPositionStorage();
//     Array<double> X_com; Array<double>& X_com_ref = X_com; 
//     Body_position->getDataColumn("center_of_mass_X", X_com_ref, 0.0);
//     
//     Array<double> X_com_RToes; Array<double>& X_com_RToes_ref = X_com_RToes; 
//     Body_position->getDataColumn("toes_r_X", X_com_RToes_ref, 0.0);
//     Array<double> Y_com_RToes; Array<double>& Y_com_RToes_ref = Y_com_RToes; 
//     Body_position->getDataColumn("toes_r_Y", Y_com_RToes_ref, 0.0);
//     Array<double> Fx_RToes; Array<double>& Fx_RToes_ref = Fx_RToes; 
//     force_storage->getDataColumn(30, Fx_RToes_ref);
//     Array<double> Fy_RToes; Array<double>& Fy_RToes_ref = Fy_RToes; 
//     force_storage->getDataColumn(31, Fy_RToes_ref);
//     Array<double> M; Array<double>& M_ref = M; 
//     force_storage->getDataColumn(35, M_ref);
//     Array<double> x_cop(0.0, M.getSize());
//     for (int i = 0; i < M.getSize(); i++){
//         x_cop[i] =((M[i]-Fx_RToes[i] * Y_com_RToes[i]) / Fy_RToes[i]) + X_com_RToes[i];
//     }
//     cout << x_cop;