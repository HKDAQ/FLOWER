/* vim:set noexpandtab tabstop=4 wrap */
//C++
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iterator>
#include <algorithm>
//ROOT
#include <TH1F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TMath.h>
#include <TStyle.h>

#if !defined(__CINT__) || defined(__MAKECINT__)
//WCSim
#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
//BONSAI
#include "WCSimBonsai.hh"
#endif

template <typename T>
std::vector<T> valsinrange(typename std::vector<T>::iterator rangebegin, typename std::vector<T>::iterator rangeend, T lowerlimit, T upperlimit){
  auto lowerit = std::lower_bound(rangebegin, rangeend, lowerlimit);
  auto upperit = std::upper_bound(rangebegin, rangeend, upperlimit);
  return typename std::vector<T> (lowerit, upperit);
}

// anonymous namespace to declare constants used within this translation unit.
namespace {
	constexpr float LIGHT_SPEED = 21.58333; // speed of light in water, value from https://github.com/hyperk/hk-BONSAI/blob/d9b227dad26fb63f2bfe80f60f7f58b5a703250a/bonsai/hits.h#L5
	constexpr int MAXDIGITS=500; // maximum number of digits we allow, due to use of preallocated vector sizes.
	constexpr int nWorkingPMTs = -1; // number of working PMTs If negative, = nPMTs from WCSimRootGeom
	constexpr float darkRate = 4.2/1000000; // dark noise rate (per ns) of the PMT (dummy value, based on 8.4kHz for B&L PMT)
	constexpr float lambdaEff = 100*100; // scattering length in cm (dummy value, based on Design Report II.2.E.1)
	constexpr float INTER_PMT_SEPARATION_CM = 106.f; // distance to neighboring PMTs is 106cm azimuthally, slightly less for ring separation, except at the end rings, which are further
	constexpr float effCoverages[] = {0.4, 0.4, 0.4, 0.4, 0.4068, 0.4244, 0.4968, 0.5956, 0.67};
	// from MC: effective photodetector coverage at an incidence angle of theta = 5, 15, ..., 85 degree
	constexpr double triggeroffset = 950.0; // WCSimRootTrigger::offset. TODO put in WCSimRootOptions
	constexpr int BOGUS_INT = std::numeric_limits<int>::max();
	constexpr float BOGUS_FLOAT = std::numeric_limits<float>::max();
	
	constexpr bool isANNIE = true; // enable modifications specific for ANNIE (i.e. z is the BEAM AXIS!)
	
	// function to calculate the reconstructed energy from the number of 'effective' hits:
	// The form and values in this function should be pre-determined your detector!
	// note a constexpr has many restrictions: it is possible this qualifier may need to be lifted if you change the form.
	constexpr float EnergyFromnEff(float nEff){
		return ((25.00 + 0.138*(nEff-189.8))*0.378);
//		float a[5]= {0.82, 0.13, -1.11*pow(10, -4), 1.25*pow(10, -6), -3.42*pow(10, -9)};
//		if (nEff<189.8) {
//			for (int n=0;n<5;n++) {
//				eRec += a[n]*pow(nEff, n);
//			}
//		} else {
//			eRec=(25.00 + 0.138*(nEff-189.8))*0.378; // TODO: dummy value; needs to be determined/tested much more precisely!
//		}
	}
}

// low energy reconstruction
int energetic_bonsai(const char *filename="../wcsim.root", bool verbose=false) {

	//setPlotStyle(); // defined below

	// set up histograms
	//TH1F *actX = new TH1F("Actual X", "Actual X", 200, -1500, 1500);
	TH1F *recTheta = new TH1F("Reconstructed Theta", "Reconstructed Theta; Theta [rads]; Num Events", 35, -3.500, 3.500);
	TH1F *recPhi = new TH1F("Reconstructed Phi", "Reconstructed Phi; Phi [rads]; Num Events", 35, -3.500, 3.500);
	TH1F *recAlpha = new TH1F("Reconstructed Alpha", "Reconstructed Alpha; Alpha [rads]; Num Events", 35, -3.500, 3.500);
	TH1F *recCherenkov = new TH1F("Reconstructed Cherenkov angle; Cherenkov Angle [rads]; Num Events", "Reconstructed Cherenkov angle", 100, -1, 1);
	TH1F *recEnergy = new TH1F("Reconstructed Energy", "Reconstructed Energy; Energy [MeV]; Num Events", 50, 0, 100);
	//TH1F *recEllipticity = new TH1F("Reconstructed ellipticity", "Reconstructed ellipticity", 100, -1, 1);
	//TH1F *recLikelihood = new TH1F("Reconstructed likelihood", "Reconstructed likelihood", 100, -100, 900);
	//TH1F *recR = new TH1F("Reconstructed R", "Reconstructed R", 200, 0, 1000);
	TH1F *hVtxDiff_X = new TH1F("Vertex Error X", "Vertex Error X; X Error [cm]; Num Events", 100, -200, 200);
	TH1F *hVtxDiff_Y = new TH1F("Vertex Error Y", "Vertex Error Y; Y Error [cm]; Num Events", 100, -200, 200);
	TH1F *hVtxDiff_Z = new TH1F("Vertex Error Z", "Vertex Error Z; Z Error [cm]; Num Events", 100, -200, 200);
	TH1F *hVtxDiff_T = new TH1F("Vertex Error T", "Vertex Error T; T Error [ns]; Num Events", 100, -200, 200);
	TH1F *hVtxDiff_Tot = new TH1F("Vertex Error Total", "Vertex Error Total; Tot Error [cm]; Num Events", 100, 0, 300);

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
	const int nPMTs = geo->GetWCNumPMT(); // total number of PMTs
	const int nWorkingPMTsToUse = (nWorkingPMTs<0) ? nPMTs : nWorkingPMTs;
	//bonsai->Init(geo, isANNIE); // WCSimBonsai modded to also introduce geo->GetWCOffset's into PMT posn's
	bonsai->Init(geo);

	TTree *tree = (TTree*)file->Get("wcsimT"); // Get a pointer to the tree from the file
	WCSimRootEvent* event = new WCSimRootEvent(); // Create WCSimRootEvent to put stuff from the tree in
	tree->SetBranchAddress("wcsimrootevent", &event); // Set branch address for reading from tree
	tree->GetBranch("wcsimrootevent")->SetAutoDelete(kTRUE); // Force deletion to prevent memory leak
	WCSimRootTrigger* trigger; // will contain triggers of event later (0: initial particle; 1..n: decay products)

	// static arrays for the fit
	float bsVertex[4],bsResult[6];
	float bsGood[3];
	std::vector<float> bsT(MAXDIGITS,0), bsQ(MAXDIGITS,0), distance(MAXDIGITS,0), PMTX(MAXDIGITS,0), PMTY(MAXDIGITS,0), PMTZ(MAXDIGITS,0), tCorrected(MAXDIGITS,0);
	std::vector<int> n50Array(MAXDIGITS,0), bsCAB(MAXDIGITS,0);
	int *bsNhit;
	int bsNsel[2];

	// Now loop over events
	for (int ev=0; ev < tree->GetEntries(); ev++) {
		// Read the event from the tree into the WCSimRootEvent instance
		tree->GetEntry(ev);
		trigger = event->GetTrigger(0);

		// See chapter 5 of doc/DetectorDocumentation.pdf in the WCSim repository
		// for more information on the structure of the root file.

		// Loop over triggers in the event
		for (int index = 0 ; index < event->GetNumberOfEvents(); index++) {
			trigger = event->GetTrigger(index);
			
			WCSimRootEventHeader* header = trigger->GetHeader();
			int triggertime = header->GetDate();
			// if (verbose) cout<<"adding offset to digits: + "<<triggertime<<" - "<<triggeroffset<<endl;
			int ncherenkovdigihits = trigger->GetNcherenkovdigihits();
			if(ncherenkovdigihits==0) continue;
			bsNhit = & ncherenkovdigihits;

			// get time, charge and PMT number for each WCSimRootCherenkovDigiHit in the trigger
			
			for (int i=0;i<ncherenkovdigihits;i++) {
				TObject *element = (trigger->GetCherenkovDigiHits())->At(i);
				WCSimRootCherenkovDigiHit *cherenkovdigihit = dynamic_cast<WCSimRootCherenkovDigiHit*>(element);

				// this will crash for high-E events, where ncherenkovdigihits > MAXDIGITS
				bsT.at(i) = cherenkovdigihit->GetT(); //+ triggertime - triggeroffset;
				bsQ.at(i) = cherenkovdigihit->GetQ();
				bsCAB.at(i) = cherenkovdigihit->GetTubeId();

				WCSimRootPMT pmt = geo->GetPMT(cherenkovdigihit->GetTubeId()-1);
				if(not isANNIE){
					PMTX[i] = pmt.GetPosition(0)-geo->GetWCOffset(0);
					PMTY[i] = pmt.GetPosition(1)-geo->GetWCOffset(1);
					PMTZ[i] = pmt.GetPosition(2)-geo->GetWCOffset(2);
				} else {
					PMTX[i] = pmt.GetPosition(0)-geo->GetWCOffset(0);
					PMTY[i] = pmt.GetPosition(2)-geo->GetWCOffset(2);
					PMTZ[i] = pmt.GetPosition(1)-geo->GetWCOffset(1);
				}
			} // End of loop over Cherenkov digihits
			

			// fit vertex position and direction using BONSAI
			cout<<". ";
			int vertexfound = bonsai->BonsaiFit( bsVertex, bsResult, bsGood, bsNsel, bsNhit, bsCAB.data(), bsT.data(), bsQ.data());
			if(vertexfound==0) continue;

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
			if(not isANNIE){
				hVtxDiff_X->Fill(bsVertex[0]-(trigger->GetVtx(0)-geo->GetWCOffset(0)));
				hVtxDiff_Y->Fill(bsVertex[1]-(trigger->GetVtx(1)-geo->GetWCOffset(1)));
				hVtxDiff_Z->Fill(bsVertex[2]-(trigger->GetVtx(2)-geo->GetWCOffset(2)));
				hVtxDiff_T->Fill(bsVertex[3]-triggertime);
				hVtxDiff_Tot->Fill(sqrt(pow(bsVertex[0]-(trigger->GetVtx(0)-geo->GetWCOffset(0)),2.)+
										pow(bsVertex[1]-(trigger->GetVtx(1)-geo->GetWCOffset(1)),2.)+
										pow(bsVertex[2]-(trigger->GetVtx(2)-geo->GetWCOffset(2)),2.)));
			} else {
				hVtxDiff_X->Fill(bsVertex[0]-(trigger->GetVtx(0)-geo->GetWCOffset(0)));
				hVtxDiff_Y->Fill(bsVertex[2]-(trigger->GetVtx(1)-geo->GetWCOffset(1)));
				hVtxDiff_Z->Fill(bsVertex[1]-(trigger->GetVtx(2)-geo->GetWCOffset(2)));
				hVtxDiff_T->Fill(bsVertex[3]-triggertime);
				hVtxDiff_Tot->Fill(sqrt(pow(bsVertex[0]-(trigger->GetVtx(0)-geo->GetWCOffset(0)),2.)+
										pow(bsVertex[2]-(trigger->GetVtx(1)-geo->GetWCOffset(1)),2.)+
										pow(bsVertex[1]-(trigger->GetVtx(2)-geo->GetWCOffset(2)),2.)));
			}


			// ****************************************
			// ** Energy estimation for this trigger **
			// ****************************************
			// For a detailed description of the energy estimation formulas in SK-IV, see ch. 4.3 of
			// http://www-sk.icrr.u-tokyo.ac.jp/sk/_pdf/articles/2016/doc_thesis_naknao.pdf

			for (int i=0;i<ncherenkovdigihits;i++) {// Loop through all WCSimRootCherenkovDigiHits in this trigger
				// get distance of hit (=PMT position) to reconstructed vertex (bsVertex[i])
				distance.at(i) = sqrt(pow((PMTX.at(i)-bsVertex[0]), 2) + pow((PMTY.at(i)-bsVertex[1]), 2) + pow((PMTZ.at(i)-bsVertex[2]), 2));

				// substract time-of-flight from measured arrival time bsT[i] --> tCorrected.at(i)
				tCorrected.at(i) = bsT.at(i) - (distance.at(i)/LIGHT_SPEED);
			}
			auto tCorrectedEnd = tCorrected.begin()+ncherenkovdigihits;

			// sort tCorrected array into ascending order
			std::sort(tCorrected.begin(),tCorrectedEnd);

			// look for the 50 ns interval with the maximum total number of hits --> start time: tMin
			float tMin;
			int n50Max=0;
			int iValue=0;

			for (int i=0;i<ncherenkovdigihits;i++) { // loop through tCorrected array: take each element as tMin and find 50ns window
				tMin = tCorrected.at(i);
				auto hitsintimerange = valsinrange(tCorrected.begin(),tCorrectedEnd,tMin,tMin+50.f);
				if(hitsintimerange.size() > n50Max){
					n50Max = hitsintimerange.size();
					iValue = i;
				}
			}
			tMin = tCorrected.at(iValue);

			// find the number of hits in the 100 ns interval corresponding to the maximal 50 ns window
			auto hitsintimerange = valsinrange(tCorrected.begin(),tCorrectedEnd,tMin,tMin+100.f);
			int n100Max = hitsintimerange.size();

			int j=0;
			std::vector<float> distance50(n50Max);
			std::vector<int> tubeID(n50Max);
			// fill arrays of distance from vertex in cm and tubeID for each hit in maximal interval
			for (int i=0;i<ncherenkovdigihits;i++) { // loop each hit in ncherenkovdigihits
				if (tMin <= tCorrected.at(i) && tCorrected.at(i) < tMin + 50.f) {
					distance50.at(j) = sqrt(pow((PMTX.at(i)-bsVertex[0]), 2) + pow((PMTY.at(i)-bsVertex[1]), 2) + pow((PMTZ.at(i)-bsVertex[2]), 2));
					tubeID.at(j) = bsCAB.at(i);
					j++;
				}
			} // end of loop over hits

			if (verbose) {
				std::cout << "tMin: " << tMin << "\n";
				std::cout << "n50Max: " << n50Max << "\n";
				std::cout << "n100Max: " << n100Max << "\n";
			}

			float nEff = 0; // effective number of hits
			float occupancy;
			for (int i=0;i<n50Max;i++) { // loop over hits in 50 ns interval and calculate nEff
				WCSimRootPMT pmt = geo->GetPMT(tubeID.at(i)-1);
				float x,y,z;
				if(not isANNIE){
					x = pmt.GetPosition(0)-geo->GetWCOffset(0);
					y = pmt.GetPosition(1)-geo->GetWCOffset(1);
					z = pmt.GetPosition(2)-geo->GetWCOffset(2);
				} else {
					x = pmt.GetPosition(0)-geo->GetWCOffset(0);
					y = pmt.GetPosition(2)-geo->GetWCOffset(2);
					z = pmt.GetPosition(1)-geo->GetWCOffset(1);
				}

				// calculate occupancy to correct for multiple hits on a single PMT
				// In a 3x3 grid around PMT 'tubeID', what proportion of PMTs has seen a hit?
				// TODO: Treat PMTs at the edge (that have fewer neighbors) differently!
				int nearbyHits = 0;
				WCSimRootPMT otherPMT;
				for (int j=0; j<n50Max; j++) { // loop through all hit PMTs, count number of hits in nearby PMTs
					otherPMT = geo->GetPMT(tubeID.at(j)-1);
					double distancetoother;
					if(not isANNIE){
						distancetoother = sqrt(pow(x - (otherPMT.GetPosition(0)-geo->GetWCOffset(0)), 2) + pow(y - (otherPMT.GetPosition(1)-geo->GetWCOffset(1)), 2) + pow(z - (otherPMT.GetPosition(2)-geo->GetWCOffset(2)), 2));
					} else {
						distancetoother = sqrt(pow(x - (otherPMT.GetPosition(0)-geo->GetWCOffset(0)), 2) + pow(y - (otherPMT.GetPosition(2)-geo->GetWCOffset(2)), 2) + pow(z - (otherPMT.GetPosition(1)-geo->GetWCOffset(1)), 2));
					}
					if (distancetoother < INTER_PMT_SEPARATION_CM) {
						nearbyHits++;
					}
				}

				double ratio = float(nearbyHits) / 9.;
				if (ratio < 1.) {
					occupancy= log(1. / (1.-ratio)) / ratio; // from Poisson statistics
				} else {
					occupancy= 3.0;
				}

				// correct for delayed hits (e.g. due to scattering)
				float lateHits = (n100Max - n50Max - (nWorkingPMTsToUse * darkRate * 50.f)) / float(n50Max);

				// substract dark noise hits
				float darkNoise = (nWorkingPMTsToUse * darkRate * 50.f) / float(n50Max);

				// calculate effective coverage to correct for photoCathodeCoverage
				// this depends on angle of incidence, see Fig. 4.5 (left) of http://www-sk.icrr.u-tokyo.ac.jp/sk/_pdf/articles/2016/doc_thesis_naknao.pdf
				// TODO: take into account azimuthal dependence as show in Fig. 4.5 (right); currently assumes phi=0
				float dotProduct;
				if(not isANNIE){
					dotProduct = pmt.GetOrientation(0)*(bsVertex[0] - x) + pmt.GetOrientation(1)*(bsVertex[1] - y) + pmt.GetOrientation(2)*(bsVertex[2] - z);
				} else {
					dotProduct = pmt.GetOrientation(0)*(bsVertex[0] - x) + pmt.GetOrientation(2)*(bsVertex[1] - y) + pmt.GetOrientation(1)*(bsVertex[2] - z);
				}
				float theta = acos( dotProduct / distance50.at(i)) * 180.f / TMath::Pi();
				if (theta > 89.99) theta = 0.f; // we have apparently mis-reconstructed the vertex, so let's set ...
				if (theta < 0.f) theta = 0.f; // ... the coverage to the most likely value of 0.4 (i.e. theta < 40 degrees)
//				float phi = 0; // dummy value
				float photoCathodeCoverage = 1 / effCoverages[int(theta/10)];

				// correct for scattering in water
				float waterTransparency = exp(distance50.at(i) / lambdaEff);

				float nEffHit = (occupancy + lateHits - darkNoise) * photoCathodeCoverage * waterTransparency;
				nEff += nEffHit;

				if (verbose) {
					std::cout << "\n*** event #" << ev << ", PMT hit #" << i << " ***************************************\n";
					std::cout << "occupancy (ratio of hits in 3x3 grid): " << occupancy << " (" << ratio << ")\n";
					std::cout << "lateHits:  " << lateHits << "\n";
					std::cout << "darkNoise: " << darkNoise << "\n";
					std::cout << "photoCathodeCoverage: " << photoCathodeCoverage << "\n";
					std::cout << "waterTransparency: " << waterTransparency << "\n";
					std::cout << "nEff for this 1 hit: " << nEffHit << "\n";
				}
			} // end of loop over hits in 50 ns interval


			nEff *= nPMTs / float(nWorkingPMTsToUse); // correct for dead PMTs; convert nWorkingPMTs to float because integer division is inaccurate

			// reconstruct energy from nEff; this is approximately linear, except at very low energies
			// TODO: determine fit parameters
			float eRec = EnergyFromnEff(nEff);
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
	
	TCanvas* c2 = new TCanvas("c2", "Second canvas", 500*nWide*winScale, 500*nHigh*winScale);
	c2->Draw();
	c2->Divide(nWide, nHigh);
	c2->cd(1); hVtxDiff_X->Draw();
	c2->cd(2); hVtxDiff_Y->Draw();
	c2->cd(3); hVtxDiff_Z->Draw();
	//c2->cd(4); hVtxDiff_T->Draw();
	c2->cd(4); hVtxDiff_Tot->Draw();

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
