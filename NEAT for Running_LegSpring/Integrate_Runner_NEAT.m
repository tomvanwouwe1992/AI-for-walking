#include <OpenSim/OpenSim.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <ctime>  
#include <windows.h>
#include "matrix.h"
#include "mex.h"
using namespace OpenSim;
using namespace SimTK;
using namespace std;

int main()
{
	return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// INPUT:
//		  0 - model name
//		  1 - position coordinates
//		  2 - velocity coordinates
//		  3 - acceleration coordinates
//		  4 - name of the two bodies that we want to coincide


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	OPEN Model
       
    char *osimFileName_ptr;
	osimFileName_ptr = mxArrayToString(prhs[0]);
	std::string osimFileName = std::string(&osimFileName_ptr[0]);
	Model osimModel = Model(osimFileName);
   	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	IMPORT position coordinates
	double *posiz;
	const int *dimPosiz;
	int Cposiz, Rposiz;
	posiz = mxGetPr(prhs[1]);
	dimPosiz = mxGetDimensions(prhs[1]);
	Cposiz = *(dimPosiz + 1);  // Columns in the positionvector
	Rposiz = *(dimPosiz + 0);  // Rows in the positionvector
 
   	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	IMPORT velocity coordinates
	double *veloc;
	veloc = mxGetPr(prhs[2]);
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	IMPORT acceleration coordinates
	double *accel;
	accel = mxGetPr(prhs[3]);   
 	   
       
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//	USA OPENSIM
	// ¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤  //
	//		---		QUESTA PARTE USATA PER DEFINIRE LE VARIABILI CHE SERVONO DOPO
	SimTK::State &sss = osimModel.initSystem();                  // Get state	
	const SimTK::MultibodySystem* mbs = &osimModel.getMultibodySystem();   // Get multibodysystem
	const SimbodyEngine& SBE = osimModel.getSimbodyEngine();               // Get simbody engine - subclass that contains methods to provide information/calculation on the model (e.g.: calculate distance between two bodies)
	SimTK::Vector posizvector = Vector(Cposiz, 0.0);   //Create empty positionvector. //SimTK::Vector posizvector = sss.getQ();
	SimTK::Vector velocvector = Vector(Cposiz, 0.0);   //Create empty velocityvector
	SimTK::Vector accelvector = Vector(Cposiz, 0.0);   //Create empty acceleration vector
	int nb = osimModel.getNumBodies();     
	SimTK::Vec3 g = osimModel.getGravity();
	SimTK::Array_<SimTK::Vector> jointTorques(Rposiz, SimTK::Vector(Cposiz, 0.0));  //Create empty torquematrix - size:rows = timepoints, columns =  number of coordinates (Cposiz)
	SimTK::Vector q(Cposiz, 0.0);
	SimTK::Vector u(Cposiz, 0.0);
	SimTK::Vector udot(Cposiz, 0.0);
	const SimTK::SimbodyMatterSubsystem& SMS = osimModel.getMatterSubsystem();
	SimTK::MobilizedBody mobod = SMS.getGround();;
	SimTK::Vector_<SimTK::SpatialVec> totCorForce_bodyFrame; //coriolis and gyroscopic wrenches of each body expressed in body origin
	SimTK::Vector_<SimTK::SpatialVec> gravForces_bodyFrame; //gravitational forces in the body frame
	SimTK::Vector_<SimTK::SpatialVec> GRF_D_bodyFrame; //modeled GRF at heel in body frame
	SimTK::Vector_<SimTK::SpatialVec> GRF_S_bodyFrame; //modeled GRF at heel in body frame
	SimTK::Vector MobilityForces_bodyFrame; //moobility forces always 0
	totCorForce_bodyFrame.resize(nb);
	gravForces_bodyFrame.resize(nb);
	GRF_D_bodyFrame.resize(nb);
	GRF_S_bodyFrame.resize(nb);
	MobilityForces_bodyFrame.resize(Cposiz);
	SimTK::Vector_<SimTK::SpatialVec> totF_frame(SMS.getNumBodies(), SimTK::SpatialVec());
	SimTK::Vector totF_system(Cposiz, 0.0);
	SimTK::Vector m_udot(Cposiz, 0.0);
	const BodySet& BS = osimModel.getBodySet();
	// ¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤-¤ //
	// Iterate over all collocation points
	for (int i = 0; i < Rposiz; i++) // Rposiz sono le righe della matrice delle coordinate che corrispondono ai varii istanti di tempo
	{
		GRF_S_bodyFrame.setToZero();
		GRF_D_bodyFrame.setToZero();
		MobilityForces_bodyFrame.setToZero();
		totCorForce_bodyFrame.setToZero();
		gravForces_bodyFrame.setToZero();

		for (int j = 0; j < Cposiz; j++)
		{
			posizvector[j] = posiz[i + j*Rposiz]; // In every iteration over a collocation point we change the value of the position/vel/Acc vector to the correct one
            velocvector[j] = veloc[i + j*Rposiz];
			accelvector[j] = accel[i + j*Rposiz];
  		}
		sss.setQ(posizvector);   // We change the state of our system to the correct one
		sss.setU(velocvector);
		mbs->realize(sss, SimTK::Stage::Dynamics);
		
        
	
		// Calculate ID torques
		for (SimTK::MobilizedBodyIndex mbi(1); mbi < SMS.getNumBodies(); ++mbi)   // We loop over all bodies (except mbi(0) - that is the ground)
		{
			mobod = SMS.getMobilizedBody(mbi);          
			SMS.addInStationForce(sss, mbi, mobod.getBodyMassCenterStation(sss), mobod.getBodyMass(sss)*g, gravForces_bodyFrame);  // We add to the gravForces_bodyframe, a force (gravity) on body 'mbi' at the COM.
		}

		// totF_frame --> all forces acting on body
		totF_frame = GRF_D_bodyFrame + GRF_S_bodyFrame + gravForces_bodyFrame;

		// MobilityForces_bodyFrame --> e.g. a perturbation
		SMS.calcResidualForceIgnoringConstraints(sss, MobilityForces_bodyFrame, totF_frame, accelvector, jointTorques[i]); //actual calculation of the jointTorques as residuals that push the mulitbodysystem in the correct state.
	}
// 	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 	//	OUTPUT
// 	// 0 - matrice ID
// 	// 1 - coordinate in Ground dei punti da seguire -> definito sopra
//  	plhs[0] = mxCreateDoubleMatrix(1, Cposiz, mxREAL);
//  	double*qqqID = mxGetPr(plhs[0]);
// 
// 		for (int c = 0; c < Cposiz; c++)
// 		{
// 			qqqID[c] = posizvector[c];
// 		}

    plhs[0] = mxCreateDoubleMatrix(Rposiz, Cposiz, mxREAL);
 	double*qqqID = mxGetPr(plhs[0]);
 	for (int r = 0; r < Rposiz; r++)
	{
		for (int c = 0; c < Cposiz; c++)
		{
			qqqID[r + c*(Rposiz)] = jointTorques[r][c];
		}
	}
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	return;
}