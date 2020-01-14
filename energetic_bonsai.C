//C++
#include <algorithm> // std::sort
#include <iostream>
#include <math.h> // pow
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
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
#include "WCSimEBonsai.h"
#endif

// low energy reconstruction
int energetic_bonsai(const char *filename="../wcsim.root", int verbose=1, const char *detector="SuperK") {
	// set up histogram
	TH1F *recEnergy = new TH1F("Reconstructed Energy", "Reconstructed Energy", 50, 0, 100);

#if !defined(__MAKECINT__)
	// Load the library with class dictionary info (create with "gmake shared")
	if (getenv ("WCSIMDIR") != NULL) {
		gSystem->Load("${WCSIMDIR}/libWCSimRoot.so");
	} else {
		gSystem->Load("../WCSim/libWCSimRoot.so");
	}
	if (getenv ("BONSAIDIR") != NULL) {
		gSystem->Load("${BONSAIDIR}/libWCSimBonsai.so");
	} else {
		gSystem->Load("../hk-BONSAI/libWCSimBonsai.so");
	}
	if (getenv ("EBONSAIDIR") != NULL) {
		gSystem->Load("${EBONSAIDIR}/libWCSimEBonsai.so");
	} else {
		gSystem->Load("../energetic-BONSAI/libWCSimEBonsai.so");
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

	WCSimEBonsai * energetic_bonsai = new WCSimEBonsai(detector, geo, verbose);

	TTree *tree = (TTree*)file->Get("wcsimT"); // Get a pointer to the tree from the file
	WCSimRootEvent* event = new WCSimRootEvent(); // Create WCSimRootEvent to put stuff from the tree in
	tree->SetBranchAddress("wcsimrootevent", &event); // Set branch address for reading from tree
	tree->GetBranch("wcsimrootevent")->SetAutoDelete(kTRUE); // Force deletion to prevent memory leak
	WCSimRootTrigger* trigger; // will contain triggers of event later (0: initial particle; 1..n: decay products)

	int ncherenkovdigihits, i;
	double x, y, z, eRec;	

	// used by bonsai->BonsaiFit()
	float bsVertex[4],bsResult[6], bsGood[3];
	int *bsNhit;
	int bsNsel[2];

	// Loop over events
	for (int ev=0; ev < tree->GetEntries(); ev++) {
		if (verbose) std::cout << "event number: " << ev << std::endl;

		// Read the event from the tree into the WCSimRootEvent instance
		tree->GetEntry(ev);
		trigger = event->GetTrigger(0);

		// See chapter 5 of doc/DetectorDocumentation.pdf in the WCSim repository
		// for more information on the structure of the root file.

		// Loop over triggers in the event
		for (int index = 0 ; index < event->GetNumberOfEvents(); index++) {
			// We don't have secondary decays so shouldn't observe more than 1 trigger per event.
			// If we do, it's noise (because the NDigits trigger in WCSim is not optimized for HK).
			// Those would get filtered later in energy cuts, but for simplicity, we ignore them here.
			if (index == 1) continue;

			trigger = event->GetTrigger(index);
			ncherenkovdigihits = trigger->GetNcherenkovdigihits();
			if (verbose) std::cout << "ncherenkovdigihits: " << ncherenkovdigihits << std::endl;
			if (ncherenkovdigihits == 0) {
				std::cout << "t, PID, 0, 0, 0, 0" << std::endl;
				continue;
			}
			std::vector<float> bsT (ncherenkovdigihits,0);
			std::vector<float> bsQ (ncherenkovdigihits,0);
			std::vector<int> bsCAB (ncherenkovdigihits,0);

			bsNhit = & ncherenkovdigihits;

			// Get time, charge and PMT number for each hit
			for (i=0; i<ncherenkovdigihits; i++) {
				TObject *element = (trigger->GetCherenkovDigiHits())->At(i);
				WCSimRootCherenkovDigiHit *cherenkovdigihit = dynamic_cast<WCSimRootCherenkovDigiHit*>(element);

				bsT[i] = cherenkovdigihit->GetT();
				bsQ[i] = cherenkovdigihit->GetQ();
				bsCAB[i] = cherenkovdigihit->GetTubeId();

				if(bsCAB[i] == 0)
				  cout << "WTF!?!? Digit has tube ID 0" << endl;
			}

			// fit vertex position and direction using BONSAI
			int* bsCAB_a = &bsCAB[0]; // turn std::vector into a pointer for BONSAI
			float* bsT_a = &bsT[0];
			float* bsQ_a = &bsQ[0];

			if (verbose) std::cout << "Fitting event vertex with hk-BONSAI ..." << std::endl;
			bonsai->BonsaiFit(bsVertex, bsResult, bsGood, bsNsel, bsNhit, bsCAB_a, bsT_a, bsQ_a);
			if (verbose) std::cout << "BONSAI done!" << std::endl;


			// Energy estimation for this trigger
			eRec = energetic_bonsai->GetEnergy(bsT, bsCAB, &bsVertex[0]);

			recEnergy->Fill(eRec);

			// reconstructed lepton direction (bsResult[0] is theta, bsResult[1] is phi)
			x = sin(bsResult[0]) * cos(bsResult[1]);
			y = sin(bsResult[0]) * sin(bsResult[1]);
			z = cos(bsResult[0]);

			// Print out reconstruction results in the format `time, PID, energy, direction (xyz), vertex (xyz)`
			// PID (i.e. electron vs. positron) and absolute time (as opposed to time relative to the trigger window) are not available.
			std::cout << "t, PID, " << eRec << ", " << x << ", " << y << ", " << z << ", " << bsVertex[0] << ", " << bsVertex[1] << ", " << bsVertex[2] << std::endl;

			// Free the memory used by these vectors
			std::vector<int>().swap(bsCAB);
			std::vector<float>().swap(bsT);
			std::vector<float>().swap(bsQ);
		} // End of loop over triggers in event

		// reinitialize event between loops
		event->ReInitialize();
	} // End of loop over events

	// display histograms
	float winScale = 0.75;
	int nWide = 1;
	int nHigh = 1;
	TCanvas* c1 = new TCanvas("c1", "First canvas", 500*nWide*winScale, 500*nHigh*winScale);
	c1->Draw();
	c1->Divide(nWide, nHigh);
	c1->cd(1); recEnergy->Draw();

	return 0;
}
