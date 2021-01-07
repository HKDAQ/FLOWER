#include "WCSimFLOWER.h"

#include <iostream>
#include <cassert>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"

using std::cout;
using std::endl;
using std::cerr;

const float WCSimFLOWER::fEffCoverages[9] = {0.4, 0.4, 0.4, 0.4, 0.4068, 0.4244, 0.4968, 0.5956, 0.67}; // from MC: coverage at theta = 5, 15, ..., 85 degree
const int FIRST_PMT2 = 100000;

WCSimFLOWER::WCSimFLOWER(const char * detectorname, WCSimRootGeom * geom, bool overwrite_nearest, int verbose)
  : fLightGroupSpeed(21.58333), // speed of light in water (in cm/ns), value from https://github.com/hyperk/hk-BONSAI/blob/d9b227dad26fb63f2bfe80f60f7f58b5a703250a/bonsai/hits.h#L5
    fLambdaEff(100*100), // scattering length in cm (based on Design Report II.2.E.1)
    fDetectorName(detectorname),
    fDetector(DetectorEnumFromString(detectorname)),
    fVerbose(verbose),
    fGeom(geom),
    fLongDuration(100),
    fShortDuration(50)
{
  switch (fDetector) {
  case kSuperK:
    fDarkRate = 4.2;
    fDarkRate2 = 0;
    fNPMTs = 11146;
    fNPMTs2 = 0;
    fNeighbourDistance = 102;
    fTopBottomDistanceLo = 1767;
    fTopBottomDistanceHi = 1768;
    break;
  case kHyperK40:
    fDarkRate = 8.4;
    fDarkRate2 = 0;
    fNPMTs = 38448;
    fNPMTs2 = 0;
    fNeighbourDistance = 102;
    fTopBottomDistanceLo = 2670;
    fTopBottomDistanceHi = 2690;
    break;
  case kHyperK20:
    fDarkRate = 8.4;
    fDarkRate2 = 0;
    fNPMTs = 19462;
    fNPMTs2 = 0;
    fNeighbourDistance = 145;
    fTopBottomDistanceLo = 2670;
    fTopBottomDistanceHi = 2690;
    break;
  case kHyperK20BnL10mPMT:
    fDarkRate = 8.4;  // dark rate of 50cm Box&Line PMTs
    fDarkRate2 = 0.3; // dark rate of 7.5cm PMTs in mPMT modules
    fNPMTs = 19208;   // total number of 50cm Box&Line PMTs
    fNPMTs2 = 182704; // total number of 7.5cm PMTs (10k mPMT modules, 19 PMTs each)
    fNeighbourDistance = 145;
    fTopBottomDistanceLo = 2670;
    fTopBottomDistanceHi = 2690;
    break;
  default:
    cerr << "Unknown detector name " << fDetectorName << endl;
    exit(-1);
    break;
  }

  fNallPMTs = fNPMTs;
  if (fNPMTs2 > 0)
    fNallPMTs = FIRST_PMT2 + fNPMTs2;

  cout << "Using default values for detector " << fDetectorName << " " << fDetector << endl;

  cout << "\tUsing default dark rate " << fDarkRate << " kHz" << endl;
  fDarkRate /= 1000000; //convert to per ns

  cout << "\tUsing default NPMTs " << fNPMTs << endl;
  
  cout << "\tUsing default dark rate for 2nd PMT type " << fDarkRate2 << " kHz" << endl;
  fDarkRate2 /= 1000000; //convert to per ns

  cout << "\tUsing default NPMTs for 2nd PMT type " << fNPMTs2 << endl;
  
  cout << "\tAssuming all PMTs are working" << endl;
  fNWorkingPMTs = fNPMTs;
  fNWorkingPMTs2 = fNPMTs2;

  cout << "\tUsing default neighbour distance " << fNeighbourDistance << endl;

  cout << "\tUsing default short duration to look for max number of hits "
       << fShortDuration << " ns" << endl;

  cout << "\tUsing default long duration to look for late hits "
       << fLongDuration << " ns" << endl;

  cout << "\tUsing default abs(z)-cut to get barrel PMTs near caps that have fewer neighbours "
       << fTopBottomDistanceLo << " to " << fTopBottomDistanceHi << endl;
  cout << "\tA similar r-cut for caps with fewer neighbours is NOT YET implemented" << endl;

  //pre-assign vector memory
  fTimes               .reserve(fNallPMTs);
  fTubeIds             .reserve(fNallPMTs);
  fDistance            .reserve(fNallPMTs);
  fTimesCorrected      .reserve(fNallPMTs);
  fTimesCorrectedSorted.reserve(fNallPMTs);

  GetNearestNeighbours(overwrite_nearest);
}

// Note: Changing these values without also changing the calibration in CorrectEnergy() may not give completely accurate results.
void WCSimFLOWER::SetDarkRate(float darkrate)
{
  fDarkRate = darkrate;
  if (fVerbose >= 0) cout << "Setting DarkRate to " << fDarkRate << " kHz" << endl;
  fDarkRate /= 1000000; //convert to per ns
}

void WCSimFLOWER::SetDarkRate2(float darkrate)
{
  fDarkRate2 = darkrate;
  if (fVerbose >= 0) cout << "Setting DarkRate to " << fDarkRate2 << " kHz" << endl;
  fDarkRate2 /= 1000000; //convert to per ns
}

void WCSimFLOWER::SetNPMTs(int npmts)
{
  fNPMTs = npmts;
  if (fVerbose >= 0) cout << "Setting NPMTs to " << fNPMTs << endl;
}

void WCSimFLOWER::SetNPMTs2(int npmts)
{
  fNPMTs2 = npmts;
  if (fVerbose >= 0) cout << "Setting NPMTs2 to " << fNPMTs2 << endl;
}

void WCSimFLOWER::SetNWorkingPMTs(int nworkingpmts)
{
  fNWorkingPMTs = nworkingpmts;
  if (fVerbose >= 0) cout << "Setting NWorkingPMTs to " << fNWorkingPMTs << endl;
}

void WCSimFLOWER::SetNWorkingPMTs2(int nworkingpmts)
{
  fNWorkingPMTs2 = nworkingpmts;
  if (fVerbose >= 0) cout << "Setting NWorkingPMTs2 to " << fNWorkingPMTs2 << endl;
}

void WCSimFLOWER::SetNeighbourDistance(float neighbour_distance, bool overwrite_nearest)
{
  fNeighbourDistance = neighbour_distance;
  if (fVerbose >= 0) cout << "Setting NeighbourDistance to " << fNeighbourDistance << endl;
  GetNearestNeighbours(overwrite_nearest);
}

void WCSimFLOWER::SetShortDuration(float shortduration)
{
  fShortDuration = shortduration;
  if (fVerbose >= 0) cout << "Setting ShortDuration to " << fShortDuration << " ns" << endl;
}

void WCSimFLOWER::SetLongDuration(float longduration)
{
  fLongDuration = longduration;
  if (fVerbose >= 0) cout << "Setting LongDuration to " << fLongDuration << " ns" << endl;
}

void WCSimFLOWER::SetTopBottomDistance(float hi, float lo)
{
  if(abs(lo) <= abs(hi)) {
    fTopBottomDistanceLo = abs(lo);
    fTopBottomDistanceHi = abs(hi);
  }
  else {
    fTopBottomDistanceHi = abs(lo);
    fTopBottomDistanceLo = abs(hi);
    cerr << "Second argument should have the smallest absolute value. This has been handled, but swap your arguments to supress this warning" << endl;
  }
  cout << "Set abs(z)-cut to get barrel PMTs near caps that have fewer neighbours "
       << fTopBottomDistanceLo << " to " << fTopBottomDistanceHi << endl;
}


WCSimFLOWER::kDetector_t WCSimFLOWER::DetectorEnumFromString(std::string name)
{
  if(!name.compare("HyperK") || !name.compare("HyperK_40perCent"))
    return kHyperK40;
  else if(!name.compare("SuperK"))
    return kSuperK;
  else if (!name.compare("HyperK_20perCent"))
    return kHyperK20;
  else if (!name.compare("HyperK_20BnL10mPMT")) // TODO: name should be consistent with https://github.com/bquilain/WCSim/blob/hybridPMT/src/WCSimDetectorConfigs.cc
    return kHyperK20BnL10mPMT;
  cerr << "DetectorEnumFromString() Unknown detector name: " << name << endl;
  exit(-1);
}

float WCSimFLOWER::GetEnergy(std::vector<float> times, std::vector<int> tubeIds, float * vertex)
{
  assert(times.size() == tubeIds.size());
  fNDigits = times.size();
  //resize the vector (shouldn't reduce capacity, but increases capacity if required)
  fDistance           .resize(fNDigits);
  fTimesCorrected      .resize(fNDigits);
  fTimesCorrectedSorted.resize(fNDigits);
  fTimes              .resize(fNDigits);
  fTubeIds            .resize(fNDigits);

  fTimes   = times;
  fTubeIds = tubeIds;
  for(int i = 0; i < 3; i++)
    fVertex[i]  = vertex[i];

  CorrectHitTimes();
  FindMaxTimeInterval();
  GetNEff();
  CorrectEnergy();

  return fERec;
}

void WCSimFLOWER::CorrectHitTimes()
{
  // Correct all hit times for time-of-flight between vertex and PMT
  float x, y, z;
  WCSimRootPMT pmt;
  for (unsigned int i = 0; i < fNDigits; i++) {
    //WCSim tube IDs run from 1 to N
    //But they're stored in a TClonesArray which runs from 0 to N-1
    //Hence a -1 here
    pmt = fGeom->GetPMT(fTubeIds[i] - 1);
    x = pmt.GetPosition(0);
    y = pmt.GetPosition(1);
    z = pmt.GetPosition(2);
    fDistance[i] = sqrt(pow((x - fVertex[0]), 2) +
			pow((y - fVertex[1]), 2) +
			pow((z - fVertex[2]), 2));
    fTimesCorrected[i] = fTimes[i] - (fDistance[i] / fLightGroupSpeed);
  }

  // Sort a copy of tCorrected, so tCorrected stays in the same order as other vectors
  fTimesCorrectedSorted = fTimesCorrected;
  std::sort(fTimesCorrectedSorted.begin(), fTimesCorrectedSorted.end());
}

void WCSimFLOWER::FindMaxTimeInterval()
{
  // Find the x ns interval with the highest number of hits, and its start time
  fNMaxShort = 0;
  fNMaxLong  = 0;
  int nshort, nlong;
  float time;
  for (unsigned int i = 0; i < fNDigits; i++) {
    // count hits within 50 and 100 ns after each entry
    nshort = 0;
    nlong  = 0;
    time   = fTimesCorrectedSorted[i];
    for (unsigned int j = i; j < fNDigits; j++) {
      if (fTimesCorrectedSorted[j] < time + fLongDuration) {
        nlong++;
        if (fTimesCorrectedSorted[j] < time + fShortDuration)
          nshort++;
      } else
        	break;
    }//j

    if (nshort > fNMaxShort) {
      fNMaxShort = nshort;
      fNMaxLong  = nlong;
      fStartTime = time;
    }
  }//i

  unsigned int j = 0;
  fDistanceShort.resize(fNMaxShort);
  fTubeIdsShort .resize(fNMaxShort);
  for (unsigned int i = 0; i < fNDigits; i++) {
    if(fTimesCorrected[i] >= fStartTime && fTimesCorrected[i] < (fStartTime + fShortDuration)) {
      fDistanceShort[j] = fDistance[i];
      fTubeIdsShort [j] = fTubeIds [i];
      if (fTubeIdsShort[j] <= 0 || fTubeIdsShort[j] > fNallPMTs || (fTubeIdsShort[j] > fNPMTs && fTubeIdsShort[j] < FIRST_PMT2))
        cerr << "fTubeIdsShort has picked up an invalid ID: " << fTubeIdsShort[j] << endl
             << "Are you using the correct input file / geometry option combination" << endl;
      j++;
    }
  }//i

  if(fVerbose > 0)
    std::cout << "Maximum of " << fNMaxShort << " hits in " 
	      << fShortDuration << " ns window starting at "
	      << fStartTime << " ns" << std::endl
	      << "In extended " << fLongDuration
	      << " ns window there are " << fNMaxLong
	      << " hits" << std::endl;
}

void WCSimFLOWER::GetNearestNeighbours(bool overwrite_root_file)
{
  TString fname = TString::Format("%s/%s_%.5f.root", GetFLOWERDataDir().Data(), fDetectorName.c_str(), fNeighbourDistance);
  TFile f(fname);
  TTree * t;
  int tubeID;
  if(!overwrite_root_file && !f.IsZombie()) {
    if(fVerbose > 0)
      cout << "GetNearestNeighbours(): Getting information from file: " << fname << endl;
    //setup tree
    f.GetObject("tNeighbours", t);
    t->SetBranchAddress("tubeID", &tubeID);
    std::vector<int> * neighbours = 0;
    t->SetBranchAddress("neighbours", &neighbours);
    //loop over tree
    for(long ipmt = 0; ipmt < t->GetEntries(); ipmt++) {
      t->GetEntry(ipmt);
      fNeighbours[tubeID] = *neighbours;
      if(fVerbose > 2)
	      cout << "\tTube " << tubeID << " has " << neighbours->size() << " neighbours" << endl;
    }//ipmt
  }//reading values
  else {
    if(fVerbose > 0)
      cout << "GetNearestNeighbours(): calculating; starting loop over " << fNPMTs << " PMTs" << endl;
    //setup tree
    std::vector<int> neighbours;
    f.Open(fname, "RECREATE");
    t = new TTree("tNeighbours", "PMT nearest neighbours");
    t->Branch("tubeID", &tubeID);
    t->Branch("neighbours", &neighbours);
    //variables for the loop
    WCSimRootPMT pmt, otherPMT;
    int tubeID_j;
    float x, y, z;
    // Loop over 50cm PMTs to find their neighbours (don't search for neighbours of mPMTs, see https://github.com/HKDAQ/FLOWER/issues/12)
    for (unsigned int ipmt = 0; ipmt < fNPMTs; ipmt++) {
      if(fVerbose > 0 && ipmt % (fNPMTs / 10) == 0)
      	cout << "Finding nearest neighbours for PMT " << ipmt << " of " << fNPMTs << endl;
      pmt = fGeom->GetPMT(ipmt);
      tubeID = pmt.GetTubeNo();
      x = pmt.GetPosition(0);
      y = pmt.GetPosition(1);
      z = pmt.GetPosition(2);
      if(fVerbose > 3)
	      cout << "Tube " << tubeID << " x,y,z " << x << "," << y << "," << z << endl;

      // loop over all PMTs and get the IDs of each ones that are closer that fNeighbourDistance
      neighbours.clear();
      for (unsigned int jpmt = 0; jpmt < fNallPMTs; jpmt++) {
        if(jpmt == ipmt) continue; // don't count the current PMT itself
        if (jpmt < FIRST_PMT2)
          otherPMT = fGeom->GetPMT(jpmt);
        else
          otherPMT = fGeom->GetPMT(jpmt - FIRST_PMT2, true);

        tubeID_j = otherPMT.GetTubeNo();
        if (sqrt(pow(x - otherPMT.GetPosition(0), 2) + 
                 pow(y - otherPMT.GetPosition(1), 2) +
                 pow(z - otherPMT.GetPosition(2), 2)) < fNeighbourDistance) {
          neighbours.push_back(tubeID_j);
        }
      }//jpmt
      fNeighbours[tubeID] = neighbours;
      if(fVerbose > 2)
      	cout << "\tTube " << tubeID << " has " << neighbours.size() << " neighbours" << endl;
      //and save in the tree
      t->Fill();
    }//ipmt
    t->Write();
    delete t;
  }//calculating values
  if(fVerbose > 1)
    cout << "GetNearestNeighbours() finished" << endl;
}

void WCSimFLOWER::GetNEff()
{
  int nearbyHits, tubeID;
  float ratio, occupancy, lateHits, darkNoise, photoCoverage, waterTransparency, liveAreaFraction;
  float nEffHit, nEffHitMod, x, y, z;
  vector<int> neighbours;
  WCSimRootPMT pmt;
  fNEff = 0;
  fNEffMod = 0;

  for (unsigned int i = 0; i < fNMaxShort; i++) {
    tubeID = fTubeIdsShort[i];
    //WCSim tube IDs run from 1 to N
    //But they're stored in a TClonesArray which runs from 0 to N-1
    //Hence a -1 here
    pmt = fGeom->GetPMT(tubeID - 1);
    x = pmt.GetPosition(0);
    y = pmt.GetPosition(1);
    z = pmt.GetPosition(2);

    // Calculate occupancy to correct for multiple hits on a single PMT
    if (tubeID <= fNPMTs) {
      nearbyHits = 0;
      neighbours = fNeighbours[tubeID];
      for(unsigned int j = 0; j < fNMaxShort; j++) {
        if(std::find(neighbours.begin(), neighbours.end(), fTubeIdsShort[j]) != neighbours.end())
          nearbyHits++;
      }//j

      // PMTs at the top/bottom edge of the barrel only have five neighbours
      if(abs(z) > fTopBottomDistanceLo && abs(z) < fTopBottomDistanceHi)
        // TODO: deal with PMTs at the edge of the top/bottom
        ratio = nearbyHits / 5.0;
      else
        ratio = nearbyHits / 8.0;

      if (ratio == 0) {
        occupancy = 1.0;
      } else if (ratio < 1) {
        occupancy = log(1 / (1-ratio)) / ratio; // from Poisson statistics
      } else {
        occupancy = 3.0;
      }
    } else {
      // Ignore occupancy correction for 3" PMTs in mPMTs, since they are highly unlikely to receive multiple hits due to their small size.
      occupancy = 1.0;
    }

    // correct for delayed hits (e.g. due to scattering)
    lateHits = (fNMaxLong - fNMaxShort - (fNWorkingPMTs * fDarkRate * fShortDuration)
                                       - (fNWorkingPMTs2 * fDarkRate2 * fShortDuration)) / fNMaxShort;

    // substract dark noise hits
    darkNoise = (fNWorkingPMTs * fDarkRate * fShortDuration + fNWorkingPMTs2 * fDarkRate2 * fShortDuration) / fNMaxShort;

    // calculate effective coverage to correct for photoCoverage
    // this depends on angle of incidence, see Fig. 4.5 (left) of http://www-sk.icrr.u-tokyo.ac.jp/sk/_pdf/articles/2016/doc_thesis_naknao.pdf
    // TODO: take into account azimuthal dependence as show in Fig. 4.5 (right); currently assumes phi=0
    float dotProduct = pmt.GetOrientation(0) * (fVertex[0] - x) +
                       pmt.GetOrientation(1) * (fVertex[1] - y) +
                       pmt.GetOrientation(2) * (fVertex[2] - z);
    float theta = acos( dotProduct / fDistanceShort[i]) * 180 / 3.14159265;
    if (!((theta >= 0.) and (theta <= 89.99))) theta = 0; // we have apparently mis-reconstructed the vertex, so let's set ...
                                                          // ... the coverage to the most likely value of 0.4 (i.e. theta < 40 degrees)
    
    photoCoverage = 1 / WCSimFLOWER::fEffCoverages[int(theta/10)]; // 1 over photocoverage, assuming SuperK or HyperK_40perCent
    if (fDetector == kHyperK20)
      photoCoverage *= 38448/float(19462); // ratio of number of PMTs
    if (fDetector == kHyperK20BnL10mPMT)
      photoCoverage *= 38448/float(19462); // ratio of number of B&L PMTs
      photoCoverage /= 1.22; // 10k mPMT modules (19 PMTs, 3" diameter) have 22% the area of 19462 20-inch PMTs

    // correct for scattering in water
    waterTransparency = exp(fDistanceShort[i] / fLambdaEff);

    nEffHit = (occupancy + lateHits - darkNoise) * photoCoverage * waterTransparency;
    fNEff += nEffHit;
    // ad-hoc modification for low photocoverage
    nEffHitMod = (pow(occupancy, 1.4) + lateHits - darkNoise) * photoCoverage * waterTransparency;
    fNEffMod += nEffHitMod;
    
    if (fVerbose > 1) {
      std::cout << "\n*** PMT hit #" << i << " on tube " << tubeID << " ***************************************\n";
      std::cout << "occupancy (ratio of hits in 3x3 grid): " << occupancy << " (" << ratio << ")\n";
      std::cout << "lateHits:  " << lateHits << "\n";
      std::cout << "darkNoise: " << darkNoise << "\n";
      std::cout << "photoCoverage: " << photoCoverage << "\n";
      std::cout << "waterTransparency: " << waterTransparency << "\n";
      std::cout << "nEff for this 1 hit: " << nEffHit << "\n";
    }
  } // i //end of loop over hits in 50 ns interval


  // correct for dead PMTs
  liveAreaFraction = (fNPMTs + 0.4275 * fNPMTs2) / (fNWorkingPMTs + 0.4275 * fNWorkingPMTs2); // photosensitive area of 19x3" PMTs is 0.4275 times that of 1x20" PMT
  fNEff *= liveAreaFraction;
  fNEffMod *= liveAreaFraction;

  if(fVerbose > 0) {
    std::cout << endl << "***************************************" << endl
	      << "nEff for this event: " << fNEff << endl
	      << " (nEffMod for low photo-coverage is: " << fNEffMod << ")" << endl;
  }
}
void WCSimFLOWER::CorrectEnergy()
{
  // reconstruct energy from fNEff; this should be approximately linear
  switch(fDetector) {
  case kSuperK:
    if (fNEff<392) 
      fERec = 0.00002*pow(fNEff, 2) + 0.039*fNEff + 1.67;
    else
      fERec = 0.0522*fNEff - 0.46;
    break;
  case kHyperK40:
    if (fNEff<1320)
      fERec = 0.02360*fNEff + 0.082;
    else
      fERec = 0.02524*fNEff - 2.081;
    break;
  case kHyperK20:
    if (fNEff<701)
      // use fNEff, as normal
      fERec = 0.00000255*pow(fNEff, 2) + 0.0215*fNEff + 0.429;
    else
      // use fNEffMod (with occupancy to power of 1.4)
      fERec = 0.000001148*pow(fNEffMod, 2) + 0.02032*fNEffMod + 1.94;
    break;
  case kHyperK20BnL10mPMT:
    fERec = 1.0*fNEff + 0.0; // TODO: calibrate relation for this detector geometry
    break;
  }
  if(fVerbose > 0) {
    std::cout << "Reconstructed energy = " << fERec << " MeV" << std::endl;
  }
}

TString WCSimFLOWER::GetFLOWERDataDir()
{
  TString dir(gSystem->Getenv("FLOWERDATADIR"));
  if(!dir.Length() || !gSystem->OpenDirectory(dir.Data())) {
    cout << "$FLOWERDATADIR not set, or points to directory that doesn't exist. Using $FLOWERDIR/data/" << endl;
    dir = TString::Format("%s/data/", gSystem->Getenv("FLOWERDIR"));
  }
  return dir;
}
