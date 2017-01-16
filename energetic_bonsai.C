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
int reconstruct_energy(char *filename="../wcsim.root", bool verbose=false) {
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
		gSystem->Load("../libWCSimRoot.so");
	}
	if (getenv ("BONSAIDIR") !=  NULL) {
		gSystem->Load("${BONSAIDIR}/libWCSimBonsai.so");
	} else {
		gSystem->Load("../libWCSimBonsai.so");
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
		int bsCAB[500];
		int *bsNhit;
		int bsNsel[2];

		for (int index = 0 ; index < event->GetNumberOfEvents(); index++) {
			trigger = event->GetTrigger(index);
			int ncherenkovdigihits = trigger->GetNcherenkovdigihits();
			bsNhit = & ncherenkovdigihits;

			// get time, charge and PMT number for each WCSimRootCherenkovDigiHit in the trigger
			for (i=0;i<ncherenkovdigihits;i++) {
				TObject *element = (trigger->GetCherenkovDigiHits())->At(i);
				WCSimRootCherenkovDigiHit *cherenkovdigihit = dynamic_cast<WCSimRootCherenkovDigiHit*>(element);

				// this will crash for high-E events, where ncherenkovdigihits > 500
				bsT[i] = cherenkovdigihit->GetT();
				bsQ[i] = cherenkovdigihit->GetQ();
				bsCAB[i] = cherenkovdigihit->GetTubeId();
			} // End of loop over Cherenkov digihits

			// fit vertex position and direction using BONSAI
			bonsai->BonsaiFit( bsVertex, bsResult, bsGood, bsNsel, bsNhit, bsCAB, bsT, bsQ);

			// fill histograms with reconstructed vertex position (bsVertex[i]) and direction (bsResult[i])
			recTheta->Fill(bsResult[0]);
			recPhi->Fill(bsResult[1]);
			recAlpha->Fill(bsResult[2]);
			recCherenkov->Fill(bsResult[3]);
			recEllipticity->Fill(bsResult[4]);
			recLikelihood->Fill(bsResult[5]);
			//recR->Fill(sqrt(pow(bsVertex[0], 2) + pow(bsVertex[1], 2) + pow(bsVertex[2], 2)));
			//actX->Fill(trigger->GetVtx(0)); // x component of the true vertex position, for comparison
			//std::cout << "reconstructed direction: " << bsResult[0] << " " << bsResult[1] << " " << bsResult[2] << " " << bsResult[3] << " " << bsResult[4] << " " << bsResult[5] << std::endl;


			// ****************************************
			// ** Energy estimation for this trigger **
			// ****************************************
			// For a detailed description of the energy estimation formulas in SK-IV, see ch. 4.3 of
			// http://www-sk.icrr.u-tokyo.ac.jp/sk/_pdf/articles/2016/doc_thesis_naknao.pdf

			for (i=0;i<ncherenkovdigihits;i++) {// Loop through all WCSimRootCherenkovDigiHits in this trigger
				// get distance of hit (=PMT position) to reconstructed vertex (bsVertex[i])

				// substract time-of-flight from measured arrival time bsT[i]
			}

			// look for 50/100 ms interval with maximal number of hits --> start/end time: tMin50, tMax50, …


			int n50 = 0; // number of hits in 50 ms interval
			for (i=0;i<ncherenkovdigihits;i++) {// Loop through elements in the TClonesArray of WCSimRootCherenkovDigiHits
				if (tMin50 < tCorrected[i] && tCorrected[i] < tMax50) {
					n50++;
					// for convenience: save distance from vertex, incident angles on PMT, tubeID to separate arrays
				}
			}

			int nEff = 0; // effective number of hits
			for (i=0;i<n50;i++) { // loop over hits in 50 ms interval and calculate nEff
				//nEff += …;
			}

			int deadPMTCorrection = 1; // total number of PMTs / number of working PMTs
			nEff *= deadPMTCorrection;


		} // End of loop over triggers in event

		// reinitialize event between loops
		event->ReInitialize();

	} // End of loop over events

	// display histograms
	float winScale = 0.75;
	int nWide = 3;
	int nHigh = 2;
	TCanvas* c1 = new TCanvas("c1", "First canvas", 500*nWide*winScale, 500*nHigh*winScale);
	c1->Draw();
	c1->Divide(nWide, nHigh);
	c1->cd(1); recTheta->Draw();
	c1->cd(2); recPhi->Draw();
	c1->cd(3); recAlpha->Draw();
	c1->cd(4); recCherenkov->Draw();
	c1->cd(5); recEllipticity->Draw();
	c1->cd(6); recLikelihood->Draw();

	return 0;
}

int occupancy(int tubeID) { // TODO
	// in 3x3 grid around PMT 'tubeID', what proportion x of PMTs has seen a hit?
	// return log(1/(1-x))/x, if x<1
	// return 3, if x=1
	return 1;
}

float photocathodeCoverage (int tubeID) { // TODO
	// dependent on angle of incidence
	return 0.4;
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
