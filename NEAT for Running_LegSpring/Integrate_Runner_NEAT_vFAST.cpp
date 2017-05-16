
//Including the necessary header files
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

//Define used namespaces
using namespace OpenSim;
using namespace SimTK;
using namespace std;


int find(Vector arr, double length, double seek)
////Function that searches a vector for a value and returns the index. 
//INPUT: 
//arr = vector to search through
//length = length of the vector being searched
//seek = double searching for
//OUTPUT: index as integer
{
	int return_value;
	for (int i = 0; i < length; ++i)
	{
		if (arr[i] == seek) { return_value = i; }
	}
	return return_value;
}

Vector Calculate_ContactInfo(Vector Right_Foot_Contact) {

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

Vector get_NN_inputs(Model* osimModel, State& s, ForceReporter* forcereporter, Vector NN_info)
////Function that returns the inputs to the neural network. These depend naturally on the model and the state it is in at the time instant.
//Inputs: 
//osimModel = pointer to the opensim model
//s = reference to state of the model
//forcereporter = pointer to the forcereporter analysis of the system
//Output: Vector containing all inputinformation
{
	Vector inputs(NN_info[0], 0.0);  //Size of the inputs vector is the number of inputs to the neural network

	const SimbodyEngine& SBE = osimModel->getSimbodyEngine();  // Get simbody engine - subclass that contains methods to provide information/calculation on the model (e.g.: calculate distance between two bodies)
	const BodySet& BS = osimModel->getBodySet();

	//Calculate some useful cody positions and velocities.
	Vec3 &RelPos = Vec3(0.0, 0.0, 0.0);
	// Left Calcaneus
	Vec3 &PosCALCNL = Vec3(0.0, 0.0, 0.0);
	SBE.getPosition(s, BS.get("calcn_l"), RelPos, PosCALCNL);
	// Right Calcaneus
	Vec3 &PosCALCNR = Vec3(0.0, 0.0, 0.0);
	SBE.getPosition(s, BS.get("calcn_r"), RelPos, PosCALCNR);
	// Head position
	Vec3 &PosHead = Vec3(0.0, 0.0, 0.0);
	SBE.getPosition(s, BS.get("head"), RelPos, PosHead);
	// Head velocity
	Vec3 &VelHead = Vec3(0.0, 0.0, 0.0);
	SBE.getVelocity(s, BS.get("head"), RelPos, VelHead);

	//Get the vector with bodyforces on the first 42 bodies (to be sure we include forces acting on the ground)
	Vector forces(42, 0.0);
	forcereporter->getForceStorage().getDataAtTime(forcereporter->getForceStorage().getLastTime(), 42, forces);
	//Get the forcec and moments from the left and right foot on the ground.
	Vector Right_Foot_Contact(6, 0.0);
	for (int l = 18; l < 24; l++) { Right_Foot_Contact[l - 18] = forces[l];	}
	Vector Left_Foot_Contact(6, 0.0);
	for (int k = 36; k < 42; k++) { Left_Foot_Contact[k - 36] = forces[k];	}
	//Calculate the COP of the contacts w.r.t. the origin
	double contact_Rx; double contact_Lx;
	if (Right_Foot_Contact[0] + Right_Foot_Contact[1] + Right_Foot_Contact[2] == 0) { contact_Rx = 0; } //Foot not in contact
	else {
		Vector contact_R = Calculate_ContactInfo(Right_Foot_Contact);
		contact_Rx = contact_R[0];
	}
    if (Left_Foot_Contact[0] + Left_Foot_Contact[1] + Left_Foot_Contact[2] == 0) { contact_Lx = 0; } //Foot not in contact
	else {
		Vector contact_L = Calculate_ContactInfo(Left_Foot_Contact);
		contact_Lx = contact_L[0];
	}
	const Vector& kinematics = s.getY();
	/* Kinematics:
	POSITION
	0. pelvis_tilt
	1. pelvis_x
	2. pelvis_y
	3. hip_flexion_r
	4. knee_angle_r
	5. ankle_angle_r
	6. hip_flexion_l
	7. knee_angle_l
	8. ankle_angle_l
	9. lumbar extension
	VELOCITY
	10. pelvis_tilt
	11. pelvis_x
	12. pelvis_y
	13. hip_flexion_r
	14. knee_angle_r
	15. ankle_angle_r
	16. hip_flexion_l
	17. knee_angle_l
	18. ankle_angle_l
	19. lumbar extension
	*/


	// Lengths of the legs, seen as one segment (pelvis to calcaneus).
	inputs[0] = sqrt((kinematics[1] - PosCALCNL[0]) * (kinematics[1] - PosCALCNL[0]) + (kinematics[2] - PosCALCNL[1]) * (kinematics[2] - PosCALCNL[1]));
	inputs[1] = sqrt((kinematics[1] - PosCALCNR[0]) * (kinematics[1] - PosCALCNR[0]) + (kinematics[2] - PosCALCNR[1]) * (kinematics[2] - PosCALCNR[1]));

	// Angle of the legs, seen as one segment (pelvis to calcaneus). atan2(y,x) 
	inputs[2] = atan2(PosCALCNL[1] - kinematics[2], PosCALCNL[0] - kinematics[1]) /  3.1415;
	inputs[3] = atan2(PosCALCNR[1] - kinematics[2], PosCALCNR[0] - kinematics[1]) /  3.1415;
	
	// Pelvis y-position
	inputs[4] = kinematics[2];

	// Pelvis speed (x & y)
	inputs[5] = kinematics[11]; // horizontal
	inputs[6] = kinematics[12]; // vertical
	
	// Trunk angle and angular velocity
	double x = PosHead[0] - kinematics[1];
	double y = PosHead[1] - kinematics[2];
	double dx = VelHead[0] - kinematics[11];
	double dy = VelHead[1] - kinematics[12];
	//Angle
	inputs[7] = atan2(y, x) / 3.1415;
	//Angular Velocity - based on derivative of atan2 =  [-y/(x²+y²)]*dx + [x/(x²+y²)]*dy
	inputs[8] = (-y / (x * x + y * y)) * dx + (x / (x * x + y * y)) * dy; 
    cout << inputs[7] << endl;
	//Ankle angle
	inputs[9] = kinematics[5];
	inputs[10] = kinematics[8];
	
	//Contact Information - COPvsCOM
    Vec3 COM = osimModel->calcMassCenterPosition(s);
	if (contact_Lx == 0.0) { } //No contact - COP is zero
	else { inputs[11] = contact_Lx - COM[0]; }
    if (contact_Rx == 0.0) { } //No contact - COP is zero
	else { inputs[12] = contact_Rx - COM[0]; }
    
	//Contact Information - vertical component of the contact force normalized by total body weight
	inputs[13] = Left_Foot_Contact[1]/(9.81*73);
	inputs[14] = Right_Foot_Contact[1]/(9.81*73);
	
	return inputs;
    /*INPUTS:
     0. Left leg length
     1. Right leg length
     2. Left leg angle wrt worldframe
     3. Right leg angle wrt worldframe
     4. Pelvis y-position
     5. Pelvis x-velocity
     6. Pelvis y-velocity
     7. Trunk angle
     8. Trunk angular velocity
     9. Left ankle angle
     10. Right ankle angle
     11. Left COP to COM
     12. Right COP to COM
     13. Left vertical contact force
     14. Right vertical contact force
     */
}

Vector NN_controller(Vector inputs, std::tuple<Vector, Vector, Vector> nodegenes, std::tuple<Vector, Vector, Vector, Vector>  connectiongenes, Vector NN_info) {
	/*Implementation of the NEAT neural network encoding into an input - output relation*/
	double change_treshold = 0.001;
	int number_of_inputnodes = NN_info[0];
	int number_of_outputnodes = NN_info[1];

	int number_of_hiddennodes = NN_info[2];
	int number_of_connections = NN_info[3];
	int no_change_count = 0;
	int index_loop = 0;

	Vector output_signal(number_of_outputnodes, 0.0);
	Vector nodegenes_ID = get<0>(nodegenes);
	Vector nodegenes_input = get<1>(nodegenes);
	Vector nodegenes_output = get<2>(nodegenes);
	Vector connectiongenes_from = get<0>(connectiongenes);
	Vector connectiongenes_to = get<1>(connectiongenes);
	Vector connectiongenes_weight = get<2>(connectiongenes);
	Vector connectiongenes_enabled = get<3>(connectiongenes);
	;

	for (int i = number_of_inputnodes + 1; i < number_of_inputnodes + number_of_outputnodes + number_of_hiddennodes + 1; i++) { nodegenes_input[i] = 0; }
	nodegenes_input[number_of_inputnodes] = 1;

	for (int i = 0; i < number_of_inputnodes; i++) { nodegenes_input[i] = inputs[i]; }
	for (int i = 0; i < number_of_inputnodes + 1; i++) { nodegenes_output[i] = nodegenes_input[i]; }

	for (int i = number_of_inputnodes + 1; i < number_of_inputnodes + number_of_outputnodes + number_of_hiddennodes + 1; i++) { nodegenes_output[i] = 1 / (1 + exp(-4.9 * nodegenes_input[i])); }

	Vector vector_node_state = nodegenes_output;
	double ID_connection_from_node;
	double ID_connection_to_node;
	double connection_weight;
	int index_connection_from_node;
	int index_connection_to_node;
	
	int checkchanges = number_of_inputnodes + number_of_outputnodes + number_of_hiddennodes + 1;

	while ((no_change_count < checkchanges) && index_loop < 3 * number_of_connections) {

		index_loop++;
		vector_node_state = nodegenes_output;

		for (int i = 0; i < number_of_connections; i++) {

			ID_connection_from_node = connectiongenes_from[i];
			ID_connection_to_node = connectiongenes_to[i];
			connection_weight = connectiongenes_weight[i];

			index_connection_from_node = find(nodegenes_ID, number_of_inputnodes + number_of_outputnodes + number_of_hiddennodes + 1, ID_connection_from_node);

			index_connection_to_node = find(nodegenes_ID, number_of_inputnodes + number_of_outputnodes + number_of_hiddennodes + 1, ID_connection_to_node);


			if (connectiongenes_enabled[i] == 1) {
				nodegenes_input[index_connection_to_node] = nodegenes_input[index_connection_to_node] + nodegenes_output[index_connection_from_node] * connection_weight;
			}
		}
		for (int j = number_of_inputnodes + 1; j < number_of_inputnodes + 1 + number_of_hiddennodes + number_of_outputnodes; j++) { nodegenes_output[j] = 1 / (1 + exp(-4.9 * nodegenes_input[j])); }

		for (int j = number_of_inputnodes + 1; j < number_of_inputnodes + 1 + number_of_hiddennodes + number_of_outputnodes; j++) { nodegenes_input[j] = 0; }

		no_change_count = 0;
		for (int j = 0; j < number_of_inputnodes + 1 + number_of_hiddennodes + number_of_outputnodes; j++)
		{
			if (nodegenes_output[j] - vector_node_state[j] < change_treshold) { no_change_count++; }
		}

	}

	for (int i = 0; i < number_of_outputnodes; i++) { output_signal[i] = nodegenes_output[number_of_inputnodes + 1 + i]; }

	return output_signal;
}



class MyPeriodicController : public PeriodicEventHandler {
	/*Implementation of the abstract class periodiceventhandler: 
	- The neural network controller to recalculate activations every 10ms based on the current state.
	- The integration is interrupted if model is going into "bad" condition (non-minimal negative velocity or too low pelvis position)
	*/
public:
	MyPeriodicController(Real eventinterval, Model* osimModel, tuple<tuple<Vector, Vector, Vector>, tuple<Vector, Vector, Vector, Vector>, Vector> NN_STRUCTURE, double* cost_function, ForceReporter* freporter) : PeriodicEventHandler(eventinterval), model(osimModel), NN_structure(NN_STRUCTURE), cost(cost_function), force_reporter(freporter) {
	}
	void handleEvent(State& state, Real accuracy, bool& shouldTerminate) const {
		
		//Get the inputs to the neural network controller - based on the model and its state.
		Vector inputs = get_NN_inputs(model, state, force_reporter, get<2>(NN_structure));
		// Part of our costfunction is the integral of the trunk angular velocity squared. We add it every period to the cost function.
		cost[0] = cost[0] + inputs[8] *  inputs[8];
		cout << "cost = " << cost[0] << endl;
		// For reasons of efficiency we want to stop the integration early if the pelvis goes to low or the velocity becomes negative.		
		if (inputs[4] < 0.70) { shouldTerminate = true; mexPrintf("Pelvis low \n"); }
		if (inputs[5] < -0.10) { shouldTerminate = true; mexPrintf("Velocity negative \n"); }
		
		// Calculate the output of the controller --> muscle activations
		Vector activation_vector = NN_controller(inputs, get<0>(NN_structure), get<1>(NN_structure), get<2>(NN_structure));

		// Change the muscle activations in the model
		Vector& activations = state.updZ();

		for (int i = 0; i < 18; i++) {
			// Cap the activations between 0.01 and 0.99
			if (activation_vector[i] > 0.99) { activation_vector[i] = 0.99; }
			if (activation_vector[i] < 0.01) { activation_vector[i] = 0.01; }

			activations[2 * i] = activation_vector[i];
		}
	}
private:
	Model* model;
	ForceReporter* force_reporter;
	tuple<tuple<Vector, Vector, Vector>, tuple<Vector, Vector, Vector, Vector>, Vector> NN_structure;
	double* cost;
};

void main(double* TimeVec, Model* osimModel, double* cost_function, std::tuple<Vector, Vector, Vector> nodegenes, std::tuple<Vector, Vector, Vector, Vector> connectiongenes, Vector NN_info, double print)
{
	std::streambuf* cout_sbuf = std::cout.rdbuf(); // save original sbuf
	std::ofstream   fout("cout.txt");
	std::cout.rdbuf(fout.rdbuf()); // redirect 'cout' to a 'fout'

	std::streambuf* cerr_sbuf = std::cerr.rdbuf(); // save original sbuf
	std::ofstream   ferr("cerr.txt");
	std::cerr.rdbuf(ferr.rdbuf()); // redirect 'cout' to a 'fout'
	cost_function[0] = 0.0;

	// Make a tuple that contains all necessary information on the neural network.
	std::tuple<tuple<Vector, Vector, Vector>, tuple<Vector, Vector, Vector, Vector>, Vector> NN_structure = make_tuple(nodegenes, connectiongenes, NN_info);
	
	// ******************** SimBodyLand***********************//
	osimModel->buildSystem();   // Build the simbody system of the model
	Real controlinterval = 0.01;  // Control interval is set to 10ms.
	
	ForceReporter* reporter = new ForceReporter(osimModel); //Add a forcereporter to the model
	osimModel->addAnalysis(reporter);
	
	// Add the implemented class MyPeriodicController to the multibodysystem
	MultibodySystem& system_ref = osimModel->updMultibodySystem();
	MyPeriodicController* myPeriodicController = new MyPeriodicController(controlinterval, osimModel, NN_structure, cost_function, reporter);
	system_ref.addEventHandler(myPeriodicController);


	// ******************** OpenSimLand***********************//
	State& s = osimModel->initializeState(); // Initialize the state of the model

	double initialTime; double finalTime; // Set boundaries on the time of the integration
	initialTime = TimeVec[0]; finalTime = TimeVec[1];

	osimModel->equilibrateMuscles(s);    // Equilibrate the muscles in the current state
	osimModel->getMultibodySystem().realize(s, Stage::Dynamics);

	//Perform the controlled integration
	RungeKuttaMersonIntegrator integrator(osimModel->getMultibodySystem());     // Initialize integrator (RungeKuttaMerson-type)
	integrator.setAccuracy(1e-6);
	Manager manager(*osimModel, integrator);
	manager.setInitialTime(initialTime);                                        // Set initial and final times of the manager
	manager.setFinalTime(finalTime);
	manager.integrate(s);

	// Get the states and forces from the simulation
	Storage output_storage = manager.getStateStorage();                         // Store&Output of the states as .sto-file
	Array<double> end_state = output_storage.getLastStateVector()->getData();
	const Storage* force_storage = &(reporter->getForceStorage());
		
	// We already added a term integrating the trunk angular velocity wrt the vertical. We need to divide through the integrationtime (number of steps at which a sample is taken :: 10ms)
	cost_function[0] = cost_function[0] / (force_storage->getLastTime() / 0.01);
	cost_function[0] = 1 / (0.25 + cost_function[0]);  // Make it a regularization term (the NEAT algorithm will maximize the cost function)
			mexPrintf("   Trunk Regularization Term = %f", cost_function[0]);
			mexPrintf("\n");
    // Calculate the distance moved of the COM in x direction
	Vec3 massCenterPos = osimModel->calcMassCenterPosition(s);
	double COM_x = massCenterPos[0];
	// Time the integration has run is part of the cost function -- measure for stability of the model
	double time = output_storage.getLastTime();

			mexPrintf("  DISTANCE WALKED  =  ");
			mexPrintf("%f", COM_x);
			mexPrintf("\n");
			mexPrintf("  TIME WALKED  =  ");
			mexPrintf("%f", time);
			mexPrintf("\n");

	double time_distance = COM_x;
	mexPrintf("Time-distance term =  %f", time_distance);
	mexPrintf("\n");

	cost_function[0] = 0.2*cost_function[0] + time_distance;

	// If print option is turned 'on' the states and forces are printed into a file.
	if (print == 1) { output_storage.print("output_states.sto");  force_storage->print("output_forces.mot"); }
	return;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	OPEN Model
	char *osimFileName_ptr;
	osimFileName_ptr = mxArrayToString(prhs[0]);
	std::string osimFileName = std::string(&osimFileName_ptr[0]);
	Model model = Model(osimFileName);
	Model *osimModel = &model;

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
	int n_nodes = *(dim_nodegenes + 1);
	int n_rows = *(dim_nodegenes + 0);
	nodegenes_ptr = mxGetPr(prhs[2]);

	Vector IDs = Vector(n_nodes, 0.0);
	Vector input = Vector(n_nodes, 0.0);
	Vector output = Vector(n_nodes, 0.0);

	for (int j = 0; j < n_nodes; j++)
	{
		IDs[j] = nodegenes_ptr[j * 4];
		input[j] = nodegenes_ptr[j * 4 + 2];
		output[j] = nodegenes_ptr[j * 4 + 3];
	}
	// Make tuple that holds all necessary info of the nodegenes (ID - inputvalue - outputvalue)
	std::tuple<Vector, Vector, Vector> nodegenes = make_tuple(IDs, input, output);

	dim_connectiongenes = mxGetDimensions(prhs[3]);
	int n_connections = *(dim_connectiongenes + 1);
	connectiongenes_ptr = mxGetPr(prhs[3]);

	Vector from = Vector(n_connections, 0.0);
	Vector to = Vector(n_connections, 0.0);
	Vector weight = Vector(n_connections, 0.0);
	Vector enabled = Vector(n_connections, 0.0);

	for (int j = 0; j < n_connections; j++)
	{
		from[j] = connectiongenes_ptr[j * 5 + 1];
		to[j] = connectiongenes_ptr[j * 5 + 2];
		weight[j] = connectiongenes_ptr[j * 5 + 3];
		enabled[j] = connectiongenes_ptr[j * 5 + 4];
	}
	// Make tuple that holds all necessary info of the connections (from node to node - weight - connection enabled?)
	std::tuple<Vector, Vector, Vector, Vector> connectiongenes = make_tuple(from, to, weight, enabled);

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	IMPORT information on the size of the neural network
	double* NN_info_ptr;
	NN_info_ptr = mxGetPr(prhs[4]);
	Vector NN_info = Vector(4, 0.0);
	NN_info[0] = NN_info_ptr[0];
	NN_info[1] = NN_info_ptr[1];
	NN_info[2] = NN_info_ptr[2];
	NN_info[3] = NN_info_ptr[3];

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	Print Information (check whether output will be printed)
	double* print_ptr;
	print_ptr = mxGetPr(prhs[5]);
	double print = print_ptr[0];

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	Prepare output to matlab
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	double *cost_function;
	cost_function = mxGetPr(plhs[0]);

	//Call main function
	main(TimeVec, osimModel, cost_function, nodegenes, connectiongenes, NN_info, print);
}


