//C++
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
//ROOT
#include <TH1F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TSystem.h>

#if !defined(__CINT__) || defined(__MAKECINT__)
//WCSim
#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
//BONSAI
#include "WCSimBonsai.hh"
#endif

// low energy reconstruction
int energetic_bonsai(char *filename="../wcsim.root", bool verbose=false) {
	//setPlotStyle(); // defined below

	// set up histograms
	//TH1F *actX = new TH1F("Actual X", "Actual X", 200, -1500, 1500);
	TH1F *recTheta = new TH1F("Reconstructed Theta", "Reconstructed Theta", 35, -3.500, 3.500);
	TH1F *recPhi = new TH1F("Reconstructed Phi", "Reconstructed Phi", 35, -3.500, 3.500);
	TH1F *recAlpha = new TH1F("Reconstructed Alpha", "Reconstructed Alpha", 35, -3.500, 3.500);
	TH1F *recCherenkov = new TH1F("Reconstructed Cherenkov angle", "Reconstructed Cherenkov angle", 100, -1, 1);
	//TH1F *recEllipticity = new TH1F("Reconstructed ellipticity", "Reconstructed ellipticity", 100, -1, 1);
	//TH1F *recLikelihood = new TH1F("Reconstructed likelihood", "Reconstructed likelihood", 100, -100, 900);
	//TH1F *recR = new TH1F("Reconstructed R", "Reconstructed R", 200, 0, 1000);

#if !defined(__MAKECINT__)
	// Load the library with class dictionary info (create with "gmake shared")
	if (getenv ("WCSIMDIR") !=  NULL) {
		gSystem->Load("${WCSIMDIR}/libWCSimRoot.so");
	} else {
		gSystem->Load("../hk-BONSAI/libWCSimRoot.so");
	}
	if (getenv ("BONSAIDIR") !=  NULL) {
		gSystem->Load("${BONSAIDIR}/libWCSimBonsai.so");
	} else {
		gSystem->Load("../hk-BONSAI/libWCSimBonsai.so");
	}
#endif

	WCSimBonsai* bonsai = new WCSimBonsai();

	// Open the ROOT file
	TFile *file = new TFile(filename,"read");
	if (!file->IsOpen()) {
		std::cout << "ERROR: Could not open input file " << filename << std::endl;
		return -1;
	}

	// Read geometry from ROOT file and initialize Bonsai
	TTree *geotree = (TTree*)file->Get("wcsimGeoT");
	WCSimRootGeom *geo = 0;
	geotree->SetBranchAddress("wcsimrootgeom", &geo);
	if (geotree->GetEntries() == 0) exit(9); // exit if no geometry is defined in the ROOT file
	geotree->GetEntry(0);
	bonsai->Init(geo);

	TTree *tree = (TTree*)file->Get("wcsimT"); // Get a pointer to the tree from the file
	WCSimRootEvent* event = new WCSimRootEvent(); // Create WCSimRootEvent to put stuff from the tree in
	tree->SetBranchAddress("wcsimrootevent", &event); // Set branch address for reading from tree
	tree->GetBranch("wcsimrootevent")->SetAutoDelete(kTRUE); // Force deletion to prevent memory leak
	WCSimRootTrigger* trigger; // will contain triggers of event later (0: initial particle; 1..n: decay products)

	// Now loop over events
	for (int ev=0; ev < tree->GetEntries(); ev++) {
		// Read the event from the tree into the WCSimRootEvent instance
		tree->GetEntry(ev);
		trigger = event->GetTrigger(0);

		// See chapter 5 of doc/DetectorDocumentation.pdf in the WCSim repository
		// for more information on the structure of the root file.

		// Loop over triggers in the event
		float bsT[500],bsQ[500];
		float bsVertex[4],bsResult[6];
		float bsGood[3];
		float distance[500];
		float tCorrected[500];
		float PMTX[500],PMTY[500],PMTZ[500];
		int  n50Array[500];
		int bsCAB[500];
		int *bsNhit;
		int bsNsel[2];

		for (int index = 0 ; index < event->GetNumberOfEvents(); index++) {
			trigger = event->GetTrigger(index);
			int ncherenkovdigihits = trigger->GetNcherenkovdigihits();
			bsNhit = & ncherenkovdigihits;

			// get time, charge and PMT number for each WCSimRootCherenkovDigiHit in the trigger
			for (int i=0;i<ncherenkovdigihits;i++) {
				TObject *element = (trigger->GetCherenkovDigiHits())->At(i);
				WCSimRootCherenkovDigiHit *cherenkovdigihit = dynamic_cast<WCSimRootCherenkovDigiHit*>(element);

				// this will crash for high-E events, where ncherenkovdigihits > 500
				bsT[i] = cherenkovdigihit->GetT();
				bsQ[i] = cherenkovdigihit->GetQ();
				bsCAB[i] = cherenkovdigihit->GetTubeId();

				WCSimRootPMT pmt = geo->GetPMT(cherenkovdigihit->GetTubeId());
				PMTX[i] = pmt.GetPosition(0);
				PMTY[i] = pmt.GetPosition(1);
				PMTZ[i] = pmt.GetPosition(2);
			} // End of loop over Cherenkov digihits

			// fit vertex position and direction using BONSAI
			bonsai->BonsaiFit( bsVertex, bsResult, bsGood, bsNsel, bsNhit, bsCAB, bsT, bsQ);

			// fill histograms with reconstructed vertex position (bsVertex[i]) and direction (bsResult[i])
			recTheta->Fill(bsResult[0]);
			recPhi->Fill(bsResult[1]);
			recAlpha->Fill(bsResult[2]);
			recCherenkov->Fill(bsResult[3]);
			//recEllipticity->Fill(bsResult[4]);
			//recLikelihood->Fill(bsResult[5]);
			//recR->Fill(sqrt(pow(bsVertex[0], 2) + pow(bsVertex[1], 2) + pow(bsVertex[2], 2)));
			//actX->Fill(trigger->GetVtx(0)); // x component of the true vertex position, for comparison
			//std::cout << "reconstructed direction: " << bsResult[0] << " " << bsResult[1] << " " << bsResult[2] << " " << bsResult[3] << " " << bsResult[4] << " " << bsResult[5] << std::endl;


			// ****************************************
			// ** Energy estimation for this trigger **
			// ****************************************
			// For a detailed description of the energy estimation formulas in SK-IV, see ch. 4.3 of
			// http://www-sk.icrr.u-tokyo.ac.jp/sk/_pdf/articles/2016/doc_thesis_naknao.pdf

			for (int i=0;i<ncherenkovdigihits;i++) {// Loop through all WCSimRootCherenkovDigiHits in this trigger
				// get distance of hit (=PMT position) to reconstructed vertex (bsVertex[i])
				distance[i] = sqrt(pow((PMTX[i]-bsVertex[0]), 2) + pow((PMTY[i]-bsVertex[1]), 2) + pow((PMTZ[i]-bsVertex[2]), 2));
				// substract time-of-flight from measured arrival time bsT[i] --> tCorrected[i]
				tCorrected[i] = bsT[i] - (distance[i]/0.225407); // speed of light in water (refraction index n=1.33)
			
			}
			
			// sort tCorrected array into ascending order
			int sort(void) {

				int tmp, j;

	                        for (int i=0;i<500;i++) { //Loop through tCorrected array

        	                        scanf("%d"; &tCorrected[i]);

                	                for (int i=0;i<499;i++) {

                        	                for (int j=i+1;j<(500-i-1);j++) {

                                	                if tCorrected[j] > tCorrected[j+1] {
							// swap
                                        	        tmp = tCorrected[j];
	                                               	tCorrected[j] = tCorrected[j+1];
                                             		tCorrected[j+1] = tmp;

                                                	}
                                        	}
                                	}
				}//end of loop through tCorrected array

			}

			int n50[500];//dummy number of elements for n50 array TODO work out sensible number for array
			int n100[500];//dummy number of elements for n100 array
			int tMin[500];
			int tMax50[500];
			float distance50[500];
			int bsCab50[500];
			int tMax100[500];

			// look for the 50 ns interval with the maximum total number of hits --> start/end time: tMin, tMax50		

			for (int i=0;i<ncherenkovdigihits;i++) { // Loop through values of tCorrected for all hits to find the number of hits for each 50 ns window for sliding minimum value from 0 to 500 ms
				
				tMin[i]=i;
				tMax50[i]=i+50;
				tMax100[i]=i+100;				

				if (tMin[i] < tCorrected[i] && tCorrected[i] < tMax50[i]) {
				
					n50[i]++;

				}
					

				if (tMin[i] < tCorrected[i] && tCorrected[i] < tMax100) { 
	
						n100[i]++;
				
				}
				
			} // end of loop over ncherenkovdigihits
			
			// find the maximum value in the n50 array

			int n50Max = n50[i];
			
			for (int i=0;i<500;i++) { //Loop through elements in the n50 array dummy value = 500
				
				if (n50[i] > n50Max) {
					
					n50Max = n50[i];		
					int iValue = i;
				}
		
			} // end of loop over n50 array
			
			// find the number of hits in the 100 ns interval corresponding to the maximal 50 ns window
			int n100Max = n100[iValue];

			// calculate distance from vertex in cm and tubeID and send to separate arrays for all hits in n50Max 

			int iMaxValue = iValue+n50Max;

			for (int i=iValue;i<(iMaxValue);i++) { // loop over hits in n50Max, (i still corresponding to range 0<i<cherenkovdigihits)

				distance50[i] = sqrt(pow((PMTX[i]-bsVertex[0]), 2) + pow((PMTY[i]-bsVertex[1]), 2) + pow((PMTZ[i]-bsVertex[2]), 2));
                                bsCAB50[i] = bsCAB[i];

                        } // end of loop over hits in n50Max
 

			int nPMTs = 1; // total number of PMTs (dummy value)
			int nWorkingPMTs = 1; // number of working PMTs (dummy value)
			int darkRate = 1; // dark noise rate of the PMT (dummy value)
			float lambdaEff = 100*100; // scattering length in cm (dummy value, based on Design Report II.2.E.1)
			float nEff = 0; // effective number of hits
			for (i=0;i<n50;i++) { // loop over hits in 50 ns interval and calculate nEff
				// correct for multiple hits on a single PMT
				float occupancy = occupancy(bsCAB50[i], n50, bsCAB50);

				// correct for delayed hits (e.g. due to scattering)
				float lateHits = (n100 - n50 - (nWorkingPMTs * darkRate * 50)) / n50;

				// substract dark noise hits
				float darkNoise = (nWorkingPMTs * darkRate * 50) / n50;

				// correct for photoCathodeCoverage
				float photoCathodeCoverage = 1 / effCoverage(bsCAB50[i], bsVertex, r50[i]);

				// correct for scattering in water
				float waterTransparency = exp(r50[i] / lambdaEff);

				// correct for quantum efficiency of PMT
				// TODO: Get this from root file, once https://github.com/WCSim/WCSim/pull/198 is merged
				float quantumEfficiency = 1/0.315; // dummy value

				float nEffHit = (occupancy + lateHits - darkNoise) * photoCathodeCoverage * waterTransparency * quantumEfficiency;
				nEff += nEffHit;
			}

			nEff *= nPMTs / nWorkingPMTs; // correct for dead PMTs


		} // End of loop over triggers in event

		// reinitialize event between loops
		event->ReInitialize();

	} // End of loop over events

	// display histograms
	float winScale = 0.75;
	int nWide = 2;
	int nHigh = 2;
	TCanvas* c1 = new TCanvas("c1", "First canvas", 500*nWide*winScale, 500*nHigh*winScale);
	c1->Draw();
	c1->Divide(nWide, nHigh);
	c1->cd(1); recTheta->Draw();
	c1->cd(2); recPhi->Draw();
	c1->cd(3); recAlpha->Draw();
	c1->cd(4); recCherenkov->Draw();
	//c1->cd(5); recEllipticity->Draw();
	//c1->cd(6); recLikelihood->Draw();

	return 0;
}

int occupancy(int tubeID, int n50, int *tubeIDs) {
	// In a 3x3 grid around PMT 'tubeID', what proportion of PMTs has seen a hit?
	// TODO: Treat PMTs at the edge (that have fewer neighbors) differently!

	WCSimRootPMT p = geo->GetPMT(tubeID);
	float x = p.GetPosition(0);
	float y = p.GetPosition(1);
	float z = p.GetPosition(2);

	int nearbyHits = 0;
	WCSimRootPMT pmt;
	for (int i=0; i<n50; i++) {
		pmt = geo->GetPMT(tubeIDs[i]);
		if (sqrt(pow(x - pmt.GetPosition(0), 2) + pow(y - pmt.GetPosition(1), 2) + pow(z - pmt.GetPosition(2), 2)) < 101) {
			// distance to neighboring PMTs is 70.71 cm (100 cm diagonally)
			nearbyHits++;
		}
	}

	float ratio = nearbyHits / 9;

	if (ratio < 1) {
		return log(1 / (1-ratio)) / ratio;
	} else {
		return 3.0;
	}
}

float effCoverage (int tubeID, float *bsVertex, float distance) {
	// dependent on angle of incidence
	WCSimRootPMT pmt = geo->GetPMT(tubeID);

	// calculate theta, phi in Fig. 4.5 (left) of http://www-sk.icrr.u-tokyo.ac.jp/sk/_pdf/articles/2016/doc_thesis_naknao.pdf
	float dotProduct = pmt.GetOrientation(0)*(bsVertex[0] - pmt.GetPosition(0)) + pmt.GetOrientation(1)*(bsVertex[1] - pmt.GetPosition(1)) + pmt.GetOrientation(2)*(bsVertex[2] - pmt.GetPosition(2));
	float incidentAngle = acos( dotProduct / distance);
	float azimuthAngle = 0; // dummy value

	// TODO: return S(theta, phi) as show in Fig. 4.5 (right)

	return 0.4; // dummy value (equivalent to incidentAngle = 0)
}

int setPlotStyle() {
	gStyle->SetOptStat(0);
	gStyle->SetCanvasColor(0);
	gStyle->SetTitleColor(1);
	gStyle->SetStatColor(0);
	gStyle->SetFrameFillColor(0);
	gStyle->SetPadColor(0);
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	gStyle->SetTitleSize(0.04);
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetFrameBorderMode(0);
	gStyle->SetFrameLineWidth(2);
	gStyle->SetPadBorderMode(0);
	gStyle->SetPalette(1);
	gStyle->SetTitleAlign(23);
	gStyle->SetTitleX(.5);
	gStyle->SetTitleY(0.99);
	gStyle->SetTitleBorderSize(0);
	gStyle->SetTitleFillColor(0);
	gStyle->SetHatchesLineWidth(2);
	gStyle->SetLineWidth(1.5);
	gStyle->SetTitleFontSize(0.07);
	gStyle->SetLabelSize(0.05,"X");
	gStyle->SetLabelSize(0.05,"Y");
	gStyle->SetTitleSize(0.04,"X");
	gStyle->SetTitleSize(0.04,"Y");
	gStyle->SetTitleBorderSize(0);
	gStyle->SetCanvasBorderMode(0);

	return 0;
}
