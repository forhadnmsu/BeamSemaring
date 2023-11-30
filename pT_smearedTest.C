#include <ROOT/RDataFrame.hxx>
#include "ana_const.h"
using namespace std;
void pT_smearedTest(){
	ROOT::EnableImplicitMT();
	gStyle->SetOptFit(1);
	gStyle->SetOptStat(0);
	TH2::SetDefaultSumw2();
	ROOT::RDataFrame df("result_mc", "data/mc_test_beampT.root"); 
	TRandom3 randomGenerator;

	auto df1 = df.Define("beampT", [&randomGenerator]() { 
			return randomGenerator.Uniform(0.0, 0.5);
			}).Define("beamPhi", [&randomGenerator]() { 
				return randomGenerator.Uniform(-3.1415, 3.1415);
				}).Define("randomE", [&randomGenerator]() {
					// Generate a random number between 0 and 1
					double randNum = randomGenerator.Rndm();
					// Assign 1.0 to 10% of the events
					if (randNum <= 0.1) {
					return 1.0;
					} else {
					return 0.0;
					}
					});

	auto filtered_df1 = df1.Define("pass2", cutRecoMC).Define("Cut", [&](bool pass_selection)->int { return pass_selection ? 1 : 0; }, {"pass2"});
	auto df2= filtered_df1.Filter("Cut==0").Define("costh2", "-999.0").Define("phi2", "-999.0");
	df2.Snapshot("result_mc", "data_beampT/data_beampT_test/result_mc_not_passing_cut.root", {"true_costh", "true_phi", "mass","costh", "phi", "pT", "xF","xT","xB", "pass", "weight", "costh2", "phi2", "beampT", "beamPhi", "randomE"});

	auto df_passed = filtered_df1.Filter("Cut==1");

	auto df3 = df_passed.Define("costh2", [&randomGenerator]( float true_pT, float pT, float xB, float costh, float phi,float mass, float px1, float py1, float pz1, float px2, float py2, float pz2, float dpx, float dpy, float dpz, double beampT, double beamPhi, double randomE ) {

			float E1 = sqrt(px1 * px1 + py1 * py1 + pz1 * pz1 + 0.10566 * 0.10566);  //mu+ E
			float E2 = sqrt(px2 * px2 + py2 * py2 + pz2 * pz2 + 0.10566 * 0.10566);  //mu- E
			TVector3 dimu3Vec(dpx, dpy, dpz);  // dimu vector
			float p_T_beam = beampT; // gettting the random beam pT
			float phi_ran = beamPhi; // getting the random phi
						 // Calculate smeared pT
			float x_component = xB * p_T_beam * std::cos(phi_ran);
			float y_component = xB * p_T_beam * std::sin(phi_ran);
			TVector3 newVector(x_component, y_component, 0.0);
			TVector3 sumVector = dimu3Vec + newVector; // getting the new vector for smeared dimu PT
			float pT_sm = sumVector.Perp();
			double mod_costh;
			if(randomE==1.0){ //10%  random reco events have  randomE==1. defined in line 19. 
			mod_costh = 2. * (E2 * pz1 - E1 * pz2) / mass / sqrt(mass * mass + pT_sm * pT_sm);
			} else{mod_costh = costh;}
			return mod_costh;
			}, {"true_pT", "pT", "xB", "costh","phi", "mass", "px1", "py1", "pz1", "px2", "py2", "pz2", "dpx", "dpy", "dpz", "beampT", "beamPhi", "randomE"});

	auto df4 = df3.Define("phi2", [&randomGenerator](float true_pT, float pT, float xB, float costh, float phi,float mass, float px1, float py1, float pz1, float px2, float py2, float pz2, float dpx, float dpy, float dpz, double beampT, double beamPhi, double randomE) {
			float E1 = sqrt(px1 * px1 + py1 * py1 + pz1 * pz1 + 0.10566 * 0.10566);
			float E2 = sqrt(px2 * px2 + py2 * py2 + pz2 * pz2 + 0.10566 * 0.10566);
			TVector3 dimu3Vec(dpx, dpy, dpz);

			float p_T_beam = beampT;
			float phi_ran = beamPhi;
			//cout << "phi_ran: " << phi_ran << endl;
			// Calculate smeared pT
			float x_component = xB * p_T_beam * std::cos(phi_ran);
			float y_component = xB * p_T_beam * std::sin(phi_ran);
			TVector3 newVector(x_component, y_component, 0.0);
			TVector3 sumVector = dimu3Vec + newVector;
			float pT_sm = sumVector.Perp();

			double mod_phi;
			if(randomE==1.0){
			mod_phi =  atan2(2.*sqrt(mass*mass + pT_sm*pT_sm)*(px2*py1 - px1*py2), mass*(px1*px1- px2*px2 + py1*py1- py2*py2));
			} else{mod_phi = phi;}

			cout << "Smeared phi===: " << mod_phi << " Default phi===: " << phi << endl;
			cout << "===============================" << endl;

			return mod_phi;
	}, { "true_pT", "pT", "xB", "costh","phi", "mass", "px1", "py1", "pz1", "px2", "py2", "pz2", "dpx", "dpy", "dpz", "beampT", "beamPhi", "randomE"});

	df4.Snapshot("result_mc", "data_beampT/data_beampT_test/mc_lh2_test_beampT_passed.root", {"true_costh", "true_phi", "mass","costh", "phi", "pT", "xF","xT","xB", "pass", "weight", "costh2", "phi2", "beampT", "beamPhi","randomE", "D1" });	

}
