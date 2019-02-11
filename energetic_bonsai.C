//C++
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
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
#endif

// low energy reconstruction
int energetic_bonsai(char *filename="../wcsim.root", bool verbose=false, bool is_HK=false) {
	// set up histogram
	TH1F *recEnergy = new TH1F("Reconstructed Energy", "Reconstructed Energy", 50, 0, 100);

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

	// used by bonsai->BonsaiFit()
	float bsVertex[4],bsResult[6], bsGood[3];
	int *bsNhit;
	int bsNsel[2];

	// used for energy reconstruction
	int i, j, ncherenkovdigihits, n50tmp, n100tmp, n50Max, n100Max, nPMTs, nWorkingPMTs, nearbyHits;
	float x,y,z, ratio, occupancy, ttmp, tStart, lateHits, darkRate, darkNoise, dotProduct, theta, photoCoverage, waterTransparency, nEffHit, nEff, eRec;
	float lambdaEff = 100*100; // scattering length in cm (based on Design Report II.2.E.1)
	const float effCoverages[] = {0.4, 0.4, 0.4, 0.4, 0.4068, 0.4244, 0.4968, 0.5956, 0.67}; // from MC: coverage at theta = 5, 15, ..., 85 degree
	WCSimRootPMT pmt, otherpmt;

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
			// TODO: NDigits trigger in WCSim is not optimized for HK, sometimes adds a second trigger for noise. We ignore those here.
			if (index == 1) continue;

			trigger = event->GetTrigger(index);
			ncherenkovdigihits = trigger->GetNcherenkovdigihits();
			if (ncherenkovdigihits == 0) {
				std::cout << "t, PID, 0, 0, 0, 0\n";
				continue;
			}
			std::vector<float> bsT (ncherenkovdigihits,0);
			std::vector<float> bsQ (ncherenkovdigihits,0);
			std::vector<float> distance (ncherenkovdigihits,0);
			std::vector<float> tCorrected (ncherenkovdigihits,0);
			std::vector<int> bsCAB (ncherenkovdigihits,0);

			bsNhit = & ncherenkovdigihits;

			// Get time, charge and PMT number for each hit
			for (i=0; i<ncherenkovdigihits; i++) {
				TObject *element = (trigger->GetCherenkovDigiHits())->At(i);
				WCSimRootCherenkovDigiHit *cherenkovdigihit = dynamic_cast<WCSimRootCherenkovDigiHit*>(element);

				bsT[i] = cherenkovdigihit->GetT();
				bsQ[i] = cherenkovdigihit->GetQ();
				bsCAB[i] = cherenkovdigihit->GetTubeId();
			}

			// fit vertex position and direction using BONSAI
			int* bsCAB_a = &bsCAB[0]; // turn std::vector into a pointer for BONSAI
			float* bsT_a = &bsT[0];
			float* bsQ_a = &bsQ[0];

			bonsai->BonsaiFit(bsVertex, bsResult, bsGood, bsNsel, bsNhit, bsCAB_a, bsT_a, bsQ_a);


			// ****************************************
			// ** Energy estimation for this trigger **
			// ****************************************
			// For a detailed description of the energy estimation formulas in SK-IV, see ch. 4.3 of
			// http://www-sk.icrr.u-tokyo.ac.jp/sk/_pdf/articles/2016/doc_thesis_naknao.pdf

			// Correct all hit times for time-of-flight between vertex and PMT
			for (i=0; i<ncherenkovdigihits; i++) {
				pmt = geo->GetPMT(bsCAB[i]);
				x = pmt.GetPosition(0);
				y = pmt.GetPosition(1);
				z = pmt.GetPosition(2);
				distance[i] = sqrt(pow((x-bsVertex[0]), 2) + pow((y-bsVertex[1]), 2) + pow((z-bsVertex[2]), 2));
				tCorrected[i] = bsT[i] - (distance[i]/21.58333); // speed of light in water, value from https://github.com/hyperk/hk-BONSAI/blob/d9b227dad26fb63f2bfe80f60f7f58b5a703250a/bonsai/hits.h#L5
			}

			// Sort a copy of tCorrected, so tCorrected stays in the same order as other vectors
			std::vector<float> tCorrected_sorted (tCorrected);
			std::sort(tCorrected_sorted.begin(), tCorrected_sorted.end());

			// Find the 50 ns interval with the highest number of hits and its start time
			n50Max = 0;
			for (i=0; i<ncherenkovdigihits; i++) { // count hits within 50 and 100 ns after each entry
				n50tmp = 0;
				n100tmp = 0;
				ttmp = tCorrected_sorted[i];
				for (j=i; j<ncherenkovdigihits; j++) {
					if (tCorrected_sorted[j] < ttmp + 100) {
						n100tmp++;
						if (tCorrected_sorted[j] < ttmp + 50) n50tmp++;
					}
				}

				if (n50tmp > n50Max) {
					n50Max = n50tmp;
					n100Max = n100tmp;
					tStart = ttmp;
				}
			}


			// Create arrays of distance from vertex (in cm) and tubeID for each hit in 50 ns interval
			std::vector<float> distance50 (n50Max,0);
			std::vector<int> tubeID (n50Max,0);
			j=0;
			for (i=0; i<ncherenkovdigihits; i++) { // loop over tCorrected, which is in the same order as the vectors distance and bsCAB
				if (tStart <= tCorrected[i] && tCorrected[i] < tStart + 50) {
					distance50[j] = distance[i];
					tubeID[j] = bsCAB[i];
					j++;
				}
			} // end of loop over hits

			if (verbose) {
				std::cout << "tStart: " << tStart << "\n";
				std::cout << "n50Max: " << n50Max << "\n";
				std::cout << "n100Max: " << n100Max << "\n";
			}

			// TODO: these are detector-dependent
			nPMTs = 11146; // total number of PMTs
			darkRate = 4.2/1000000; // dark noise rate (per ns) of the PMT (8.4kHz for B&L PMT)
			if (is_HK) {
				nPMTs = 38448;
				darkRate = 8.4/1000000;
			}

			nWorkingPMTs = nPMTs; // number of working PMTs (may be lower in real detector due to defects)
			nEff = 0; // effective number of hits

			for (i=0; i<n50Max; i++) { // loop over hits in 50 ns interval and calculate nEff
				pmt = geo->GetPMT(tubeID[i]);
				x = pmt.GetPosition(0);
				y = pmt.GetPosition(1);
				z = pmt.GetPosition(2);

				// Calculate occupancy to correct for multiple hits on a single PMT:
				// In a 3x3 grid around PMT 'tubeID', what proportion of PMTs has seen a hit?
				nearbyHits = 0;
				for (j=0; j<n50Max; j++) { // loop through all hit PMTs, count number of hits in nearby PMTs
					otherPMT = geo->GetPMT(tubeID[j]);
					if (sqrt(pow(x - otherPMT.GetPosition(0), 2) + pow(y - otherPMT.GetPosition(1), 2) + pow(z - otherPMT.GetPosition(2), 2)) < 101) {
						// distance to neighboring PMTs is 70.71 cm (100 cm diagonally)
						nearbyHits++;
					}
				}

				ratio = nearbyHits / 8.0;
				// PMTs at the top/bottom edge of the barrel only have five neighbours
				if ((z < 2690)&&(z>2670) || (z > -2690)&&(z < -2670)) {
					ratio = nearbyHits / 5.0;
				} // TODO: deal with PMTs at the edge of the top/bottom

				if (ratio < 1) {
					occupancy = log(1 / (1-ratio)) / ratio; // from Poisson statistics
				} else {
					occupancy = 3.0;
				}


				// correct for delayed hits (e.g. due to scattering)
				lateHits = (n100Max - n50Max - (nWorkingPMTs * darkRate * 50)) / n50Max;

				// substract dark noise hits
				darkNoise = (nWorkingPMTs * darkRate * 50) / n50Max;

				// calculate effective coverage to correct for photoCoverage
				// this depends on angle of incidence, see Fig. 4.5 (left) of http://www-sk.icrr.u-tokyo.ac.jp/sk/_pdf/articles/2016/doc_thesis_naknao.pdf
				// TODO: take into account azimuthal dependence as show in Fig. 4.5 (right); currently assumes phi=0
				dotProduct = pmt.GetOrientation(0)*(bsVertex[0] - x) + pmt.GetOrientation(1)*(bsVertex[1] - y) + pmt.GetOrientation(2)*(bsVertex[2] - z);
				theta = acos( dotProduct / distance50[i]) * 180 / 3.14159265;
				if (theta > 89.99) theta = 0; // we have apparently mis-reconstructed the vertex, so let's set ...
				if (theta < 0) theta = 0; // ... the coverage to the most likely value of 0.4 (i.e. theta < 40 degrees)
				photoCoverage = 1 / effCoverages[int(theta/10)];

				// correct for scattering in water
				waterTransparency = exp(distance50[i] / lambdaEff);

				nEffHit = (occupancy + lateHits - darkNoise) * photoCoverage * waterTransparency;
				nEff += nEffHit;

				if (verbose) {
					std::cout << "\n*** event #" << ev << ", PMT hit #" << i << " ***************************************\n";
					std::cout << "occupancy (ratio of hits in 3x3 grid): " << occupancy << " (" << ratio << ")\n";
					std::cout << "lateHits:  " << lateHits << "\n";
					std::cout << "darkNoise: " << darkNoise << "\n";
					std::cout << "photoCoverage: " << photoCoverage << "\n";
					std::cout << "waterTransparency: " << waterTransparency << "\n";
					std::cout << "nEff for this 1 hit: " << nEffHit << "\n";
				}
			} // end of loop over hits in 50 ns interval


			nEff *= nPMTs / float(nWorkingPMTs); // correct for dead PMTs; convert nWorkingPMTs to float because integer division is inaccurate

			// reconstruct energy from nEff; this is approximately linear, except at very low energies
			// TODO: determine fit parameters
// 			float a[5]= {0.82, 0.13, -1.11*pow(10, -4), 1.25*pow(10, -6), -3.42*pow(10, -9)};
// 			if (nEff<189.8) {
// 				for (int n=0;n<5;n++) {
// 					eRec += a[n]*pow(nEff, n);
// 				}
// 			} else {
				eRec=(25.00 + 0.138*(nEff-189.8))*0.378; // TODO: dummy value; needs to be determined/tested much more precisely!
// 			}
			if (is_HK) {
				eRec = eRec /2.584; // TODO: HK observes more photons for a given energy due to different PMTs etc.
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
	int nWide = 1;
	int nHigh = 1;
	TCanvas* c1 = new TCanvas("c1", "First canvas", 500*nWide*winScale, 500*nHigh*winScale);
	c1->Draw();
	c1->Divide(nWide, nHigh);
	c1->cd(1); recEnergy->Draw();

	return 0;
}
