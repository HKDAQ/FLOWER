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
#include <TStopwatch.h>

//WCSim
#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include "WCSimRootOptions.hh"
//BONSAI
#include "WCSimBonsai.hh"
//FLOWER
#include "WCSimFLOWER.h"

// low energy reconstruction
int flower_with_bonsai(const char *detector, //sets the default nearest neighbour distances etc. Check DetectorEnumFromString inside WCSimFLOWER.cpp for allowed values
		       const char *filename="../wcsim.root",
		       const int  verbose=1,
		       const bool overwrite_nearest = false, //if true, will overwrite the cached nearest neighbours file
		       const double override_dark_rate = -99 //if positive, will override the default dark rate (taken from WCSimRootOptions)
		       )
{
	// set up histogram
	TH1F *recEnergy = new TH1F("hErec", "Reconstructed Energy", 50, 0, 100);

	WCSimBonsai* bonsai = new WCSimBonsai();

	// Open the ROOT file
	TFile *file = new TFile(filename,"read");
	if (!file->IsOpen()) {
		std::cout << "ERROR: Could not open input file " << filename << std::endl;
		return -1;
	}

	// Read geometry from WCSim file
	TTree *geotree = (TTree*)file->Get("wcsimGeoT");
	WCSimRootGeom *geo = 0;
	geotree->SetBranchAddress("wcsimrootgeom", &geo);
	if (geotree->GetEntries() == 0) exit(9); // exit if no geometry is defined in the WCSim file
	geotree->GetEntry(0);

	//Initalise BONSAI
	bonsai->Init(geo);

	// Read options tree from WCSim file
	TTree * opttree;
	file->GetObject("wcsimRootOptionsT", opttree);
	WCSimRootOptions * opt = 0;
	opttree->SetBranchAddress("wcsimrootoptions", &opt);
	if(opttree->GetEntries() == 0) exit(9); // exit if no options are defined in the WCSim file
	opttree->GetEntry(0);
	
	//Initialise FLOWER
	bool get_npmts_from_wcsimrootgeom = true;
	WCSimFLOWER * flower = new WCSimFLOWER(detector, geo, get_npmts_from_wcsimrootgeom, overwrite_nearest, verbose);
	if(override_dark_rate >= 0)
	  flower->SetDarkRate(override_dark_rate);
	else
	  flower->SetDarkRate(opt->GetPMTDarkRate("tank"));
	if(!flower->CheckNearestNeighbours()) {
	  cerr << "Nearest neighbour distribution failure. Set nearest neighbour and try again" << endl;
	  exit(8);
	}
	  
	//Read event tree from WCSim file
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
	cout << "Starting timer for the event loop" << endl;
	TStopwatch timer;
	for (int ev=0; ev < tree->GetEntries(); ev++) {
		if (verbose) std::cout << "event number: " << ev << std::endl;

		// Read the event from the tree into the WCSimRootEvent instance
		tree->GetEntry(ev);
		trigger = event->GetTrigger(0);

		// See chapter 5 of doc/DetectorDocumentation.pdf in the WCSim repository
		// for more information on the structure of the root file.

		// Loop over triggers in the event
		for (int index = 0 ; index < event->GetNumberOfEvents(); index++) {
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
				  cout << "Digit has tube ID 0. This shouldn't happen. WCSim TubeIds run from 1 to N. Has WCSim changed?" << endl;
			}

			// fit vertex position and direction using BONSAI
			int* bsCAB_a = &bsCAB[0]; // turn std::vector into a pointer for BONSAI
			float* bsT_a = &bsT[0];
			float* bsQ_a = &bsQ[0];

			if (verbose) std::cout << "Fitting event vertex with hk-BONSAI ..." << std::endl;
			bonsai->BonsaiFit(bsVertex, bsResult, bsGood, bsNsel, bsNhit, bsCAB_a, bsT_a, bsQ_a);
			if (verbose){
			  std::cout << "Vertex found at:";
			  for(int iv = 0; iv < 4; iv++)
			    std::cout << " " << bsVertex[iv];
			  std::cout << std::endl << "BONSAI done!" << std::endl;
			}

			// Energy estimation for this trigger
			eRec = flower->GetEnergy(bsT, bsCAB, &bsVertex[0]);

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
	timer.Print();

	// display histogram(s)
	float winScale = 0.75;
	int nWide = 1;
	int nHigh = 1;
	TCanvas* c1 = new TCanvas("c1", "First canvas", 500*nWide*winScale, 500*nHigh*winScale);
	c1->Draw();
	c1->Divide(nWide, nHigh);
	c1->cd(1); recEnergy->Draw();
	c1->SaveAs("Erec.pdf");

	return 0;
}
