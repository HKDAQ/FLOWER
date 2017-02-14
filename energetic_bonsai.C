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
	TH1F *recEnergy = new TH1F("Reconstructed Energy", "Reconstructed Energy", 25, 0, 500);
	//TH1F *recEllipticity = new TH1F("Reconstructed ellipticity", "Reconstructed ellipticity", 100, -1, 1);
	//TH1F *recLikelihood = new TH1F("Reconstructed likelihood", "Reconstructed likelihood", 100, -100, 900);
	//TH1F *recR = new TH1F("Reconstructed R", "Reconstructed R", 200, 0, 1000);

#if !defined(__MAKECINT__)
	// Load the library with class dictionary info (create with "gmake shared")
	if (getenv ("WCSIMDIR") != NULL) {
		gSystem->Load("${WCSIMDIR}/libWCSimRoot.so");
	} else {
		gSystem->Load("../hk-BONSAI/libWCSimRoot.so");
	}
	if (getenv ("BONSAIDIR") != NULL) {
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
		int n50Array[500];
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
				tCorrected[i] = bsT[i] - (distance[i]/21.58333); // speed of light in water, value from https://github.com/hyperk/hk-BONSAI/blob/d9b227dad26fb63f2bfe80f60f7f58b5a703250a/bonsai/hits.h#L5
			}

			// create a copy of tCorrected
			float tCorrectedSort[500];
			for (int i=0;i<500;i++) {
				tCorrectedSort[i] = tCorrected[i];
			}

			// sort tCorrectedSort array into ascending order
			int tmp;
			for (int i=0;i<500;i++) {
				for (int j=i+1;j<500;j++) {
					if (tCorrectedSort[i] > tCorrectedSort[j]) {
						tmp = tCorrectedSort[i];
						tCorrectedSort[i] = tCorrectedSort[j];
						tCorrectedSort[j] = tmp;
					}
				}
			}

			int n50[500];
			int n100[500];

			// look for the 50 ns interval with the maximum total number of hits --> start time: tMin
			float tMin;
			int n50tmp;
			int n100tmp;

			for (int i=0;i<500;i++) { // loop through tCorrectedSort array: take each element as tMin and find 50 and 100 ns window
				tMin = tCorrectedSort[i];
				for (int j=0;j<500;j++) { // loop over tCorrected array to find number of hits in each window
					if (tMin <tCorrectedSort[j] && tCorrectedSort[j] < tMin + 100) {
						n100tmp++;
						if (tCorrectedSort[j] < tMin + 50) n50tmp++;
					}
				} // end of loop over tCorrectedSort array to allocate hits to 50 ns window

				//create arrays giving number of hits in each 50ns and 100ns interval
				n50[i] = n50tmp;
				n50tmp=0;
				n100[i] = n100tmp;
				n100tmp = 0;
			}

			// find the maximum value in the n50 array
			int n50Max = 0;
			int iValue;

			for (int i=0;i<500;i++) { //Loop through elements in the n50 array
				if (n50[i] > n50Max) {
					n50Max = n50[i];
					iValue = i;
				}
			} // end of loop over n50 array

			// find the number of hits in the 100 ns interval corresponding to the maximal 50 ns window
			int n100Max = n100[iValue];

			float distance50[500];
			int tubeID[500];
			int j=0;
			// create arrays of distance from vertex in cm and tubeID for each hit in maximal interval
			// NB arrays have 500 elements but only nMax50 elements are required
			// TODO either change array length to n50Max, or will need to specify length when using arrays (esp tubeID)
			for (int i=0;i<ncherenkovdigihits;i++) { // loop each hit in ncherenkovdigihits
				tMin = tCorrectedSort[iValue];
				if (tMin <tCorrected[i] && tCorrected[i] < tMin + 50) {
					distance50[j] = sqrt(pow((PMTX[i]-bsVertex[0]), 2) + pow((PMTY[i]-bsVertex[1]), 2) + pow((PMTZ[i]-bsVertex[2]), 2));
					tubeID[j] = bsCAB[i];
					j++;
				}
			} // end of loop over hits

			if (verbose) {
				std::cout << "tMin: " << tMin << "\n";
				std::cout << "n50Max: " << n50Max << "\n";
				std::cout << "n100Max: " << n100Max << "\n";
			}

			int nPMTs = 11146; // total number of PMTs (dummy value)
			int nWorkingPMTs = 11146; // number of working PMTs (dummy value)
			float darkRate = 4.2/1000000; // dark noise rate (per ns) of the PMT (dummy value, based on 8.4kHz for B&L PMT)
			float lambdaEff = 100*100; // scattering length in cm (dummy value, based on Design Report II.2.E.1)
			float nEff = 0; // effective number of hits
			float occupancy;
			float eRecArray[500];
			for (i=0;i<n50Max;i++) { // loop over hits in 50 ns interval and calculate nEff

				// calculate occupancy to correct for multiple hits on a single PMT
				// In a 3x3 grid around PMT 'tubeID', what proportion of PMTs has seen a hit?
				// TODO: Treat PMTs at the edge (that have fewer neighbors) differently!

				WCSimRootPMT p = geo->GetPMT(tubeID[i]);
				float x = p.GetPosition(0);
				float y = p.GetPosition(1);
				float z = p.GetPosition(2);

				int nearbyHits = 0;
				WCSimRootPMT pmt;
				for (int j=0; j<n50Max; j++) { // loop through all hit PMTs, count number of hits in nearby PMTs
					pmt = geo->GetPMT(tubeID[j]);
					if (sqrt(pow(x - pmt.GetPosition(0), 2) + pow(y - pmt.GetPosition(1), 2) + pow(z - pmt.GetPosition(2), 2)) < 101) {
						// distance to neighboring PMTs is 70.71 cm (100 cm diagonally)
						nearbyHits++;
					}
				}

				double ratio = float(nearbyHits) / 9;
				if (ratio < 1) {
					occupancy= log(1 / (1-ratio)) / ratio; // from Poisson statistics
				} else {
					occupancy= 3.0;
				}


				// correct for delayed hits (e.g. due to scattering)
				float lateHits = (n100Max - n50Max - (nWorkingPMTs * darkRate * 50)) / float(n50Max);

				// substract dark noise hits
				float darkNoise = (nWorkingPMTs * darkRate * 50) / float(n50Max);

				// calculate effective coverage to correct for photoCathodeCoverage
				// this depends on angle of incidence, see Fig. 4.5 (left) of http://www-sk.icrr.u-tokyo.ac.jp/sk/_pdf/articles/2016/doc_thesis_naknao.pdf
				WCSimRootPMT pmt = geo->GetPMT(tubeID[i]);
//				float dotProduct = pmt.GetOrientation(0)*(bsVertex[0] - pmt.GetPosition(0)) + pmt.GetOrientation(1)*(bsVertex[1] - pmt.GetPosition(1)) + pmt.GetOrientation(2)*(bsVertex[2] - pmt.GetPosition(2));
//				float incidentAngle = acos( dotProduct / distance50[i]);
//				float azimuthAngle = 0; // dummy value

				// TODO: return S(incidentAngle, azimuthAngle) as show in Fig. 4.5 (right)
				float effCoverage= 0.4; // dummy value (equivalent to incidentAngle = 0)

				float photoCathodeCoverage = 1 / effCoverage;

				// correct for scattering in water
				float waterTransparency = exp(distance50[i] / lambdaEff);

				float nEffHit = (occupancy + lateHits - darkNoise) * photoCathodeCoverage * waterTransparency;
				nEff += nEffHit;

				if (verbose) {
					std::cout << "\n*** i = " << i << " ***************************************\n";
					std::cout << "occupancy (ratio of hits in 3x3 grid): " << occupancy << " (" << ratio << ")\n";
					std::cout << "lateHits: " << lateHits << ")\n";
					std::cout << "darkNoise: " << darkNoise << ")\n";
					std::cout << "photoCathodeCoverage: " << photoCathodeCoverage << ")\n";
					std::cout << "waterTransparency: " << waterTransparency << ")\n";
					std::cout << "nEff for this 1 hit: " << nEffHit << ")\n";
				}
			} // end of loop over hits in 50 ns interval


			nEff *= nPMTs / float(nWorkingPMTs); // correct for dead PMTs; convert nWorkingPMTs to float because integer division is inaccurate

			// reconstruct energy from nEff; this is approximately linear, except at very low energies
			// TODO: determine fit parameters
			float eRec = 0;
			float a[5]= {0.82, 0.13, -1.11*pow(10, -4), 1.25*pow(10, -6), -3.42*pow(10, -9)};
			if (nEff<189.8) {
				for (int n=0;n<5;n++) {
					eRec += a[n]*pow(nEff, n);
				}
			} else {
				eRec=25.00 + 0.138*(nEff-189.8);
			}
			recEnergy->Fill(eRec);

			if (verbose) {
				std::cout << "Neff = " << nEff << std::endl;
				std::cout << "Reconstructed energy = " << eRec << std::endl;
			}
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
	c1->cd(3); recEnergy->Draw();
	c1->cd(4); recCherenkov->Draw();
	//c1->cd(5); recEllipticity->Draw();
	//c1->cd(6); recLikelihood->Draw();

	return 0;
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
