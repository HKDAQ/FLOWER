#include "WCSimEBonsai.h"

#include <iostream>
#include <cassert>
#include <algorithm>

using std::cout;
using std::endl;
using std::cerr;

const float WCSimEBonsai::fEffCoverages[9] = {0.4, 0.4, 0.4, 0.4, 0.4068, 0.4244, 0.4968, 0.5956, 0.67}; // from MC: coverage at theta = 5, 15, ..., 85 degree

WCSimEBonsai::WCSimEBonsai(const char * detectorname, WCSimRootGeom * geom, int verbose)
  : fLightGroupSpeed(21.58333), // speed of light in water, value from https://github.com/hyperk/hk-BONSAI/blob/d9b227dad26fb63f2bfe80f60f7f58b5a703250a/bonsai/hits.h#L5
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
    fNPMTs = 11146;
    fNeighbourDistance = 102;
    fTopBottomDistanceLo = 1767;
    fTopBottomDistanceHi = 1768;
    break;
  case kHyperK40:
    fDarkRate = 8.4;
    fNPMTs = 38448;
    fNeighbourDistance = 102;
    fTopBottomDistanceLo = 2670;
    fTopBottomDistanceHi = 2690;
    break;
  case kHyperK20:
    fDarkRate = 8.4;
    fNPMTs = 19462;
    fNeighbourDistance = 145;
    fTopBottomDistanceLo = 2670;
    fTopBottomDistanceHi = 2690;
    break;
  default:
    cerr << "Unknown detector name " << fDetectorName << endl;
    exit(-1);
    break;
  }

  cout << "Using default values for detector " << fDetectorName << " " << fDetector << endl;

  cout << "\tUsing default dark rate " << fDarkRate << " kHz" << endl;
  fDarkRate /= 1000000; //convert to per ns

  cout << "\tUsing default NPMTs " << fNPMTs << endl;
  
  cout << "\tAssumming all PMTs are working" << endl;
  fNWorkingPMTs = fNPMTs;

  cout << "\tUsing default neighbour distance " << fNeighbourDistance << endl;

  cout << "\tUsing default short duration to look for max number of hits "
       << fShortDuration << " ns" << endl;

  cout << "\tUsing default long duration to look for late hits "
       << fLongDuration << " ns" << endl;

  cout << "\tUsing default abs(z)-cut to get barrel PMTs near caps that have fewer neighbours "
       << fTopBottomDistanceLo << " to " << fTopBottomDistanceHi << endl;
  cout << "\tA similar r-cut for caps with fewer neighbours is NOT YET implemented" << endl;

  //pre-assign vector memory
  fTimes               .reserve(fNPMTs);
  fTubeIds             .reserve(fNPMTs);
  fDistance            .reserve(fNPMTs);
  fTimesCorrected      .reserve(fNPMTs);
  fTimesCorrectedSorted.reserve(fNPMTs);

  GetNearestNeighbours();
}

void WCSimEBonsai::SetDarkRate(float darkrate)
{
  fDarkRate = darkrate;
  cout << "Setting DarkRate to " << fDarkRate << " kHz" << endl;
  fDarkRate /= 1000000; //convert to per ns
}

void WCSimEBonsai::SetNPMTs(int npmts)
{
  fNPMTs = npmts;
  cout << "Setting NPMTs to " << fNPMTs << endl;
}

void WCSimEBonsai::SetNWorkingPMTs(int nworkingpmts)
{
  fNWorkingPMTs = nworkingpmts;
  cout << "Setting NWorkingPMTs to " << fNWorkingPMTs << endl;
}

void WCSimEBonsai::SetNeighbourDistance(float neighbourdistance)
{
  fNeighbourDistance = neighbourdistance;
  cout << "Setting NeighbourDistance to " << fNeighbourDistance << endl;
}

void WCSimEBonsai::SetShortDuration(float shortduration)
{
  fShortDuration = shortduration;
  cout << "Setting ShortDuration to " << fShortDuration << " ns" << endl;
}

void WCSimEBonsai::SetLongDuration(float longduration)
{
  fLongDuration = longduration;
  cout << "Setting LongDuration to " << fLongDuration << " ns" << endl;
}

void WCSimEBonsai::SetTopBottomDistance(float hi, float lo)
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


WCSimEBonsai::kDetector_t WCSimEBonsai::DetectorEnumFromString(std::string name)
{
  if(!name.compare("HyperK") || !name.compare("HyperK_40perCent"))
    return kHyperK40;
  else if(!name.compare("SuperK"))
    return kSuperK;
  else if (!name.compare("HyperK_20perCent"))
    return kHyperK20;
  cerr << "FromString() Unknown detector name " << name << endl;
  exit(-1);
}

float WCSimEBonsai::GetEnergy(std::vector<float> times, std::vector<int> tubeIds, float * vertex)
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
  for(int i = 0; i < 4; i++)
    fVertex[i]  = vertex[i];

  CorrectHitTimes();
  FindMaxTimeInterval();
  GetNEff();
  CorrectEnergy();

  return fERec;
}

void WCSimEBonsai::CorrectHitTimes()
{
  // Correct all hit times for time-of-flight between vertex and PMT
  float x, y, z;
  WCSimRootPMT pmt;
  for (unsigned int i = 0; i < fNDigits; i++) {
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

void WCSimEBonsai::FindMaxTimeInterval()
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
      }
      else
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
      if(fTubeIdsShort[j] <= 0 || fTubeIdsShort[j] > fNPMTs)
	cerr << "WTF?! TubeId short has picked up an invalid ID" << fTubeIdsShort[j] << endl;
      j++;
    }
  }//i

  if(fVerbose)
    std::cout << "Maximum of " << fNMaxShort << " hits in " 
	      << fShortDuration << " ns window starting at "
	      << fStartTime << " ns" << std::endl
	      << "In extended " << fLongDuration
	      << " ns window there are " << fNMaxLong
	      << " hits" << std::endl;
}

void WCSimEBonsai::GetNearestNeighbours()
{
  WCSimRootPMT pmt, otherPMT;
  int tubeID_i, tubeID_j;
  float x, y, z;
  for (unsigned int i = 0; i < fNPMTs; i++) {
    pmt = fGeom->GetPMT(i);
    tubeID_i = pmt.GetTubeNo();
    x = pmt.GetPosition(0);
    y = pmt.GetPosition(1);
    z = pmt.GetPosition(2);
    if(fVerbose > 2)
      cout << "Tube " << tubeID_i << " x,y,z " << x << "," << y << "," << z << endl;

    // loop over all PMTs and get the IDs of each ones that are "close"
    // In a 3x3 grid around PMT 'tubeID'
    // explitly uses distance to neighboring PMTs 
    //  (for 40% coverage with 50cm PMTs, this is 70.71 cm (100 cm diagonally))
    vector<int> thisNeighbours;

    for (unsigned int j = 0; j < fNPMTs; j++) {
      if (j == i) continue; // don't count the current PMT itself
      otherPMT = fGeom->GetPMT(j);
      tubeID_j = otherPMT.GetTubeNo();
      if (sqrt(pow(x - otherPMT.GetPosition(0), 2) + 
	       pow(y - otherPMT.GetPosition(1), 2) +
	       pow(z - otherPMT.GetPosition(2), 2)) < fNeighbourDistance) {
	thisNeighbours.push_back(tubeID_j);
      }
    }//j
    fNeighbours[tubeID_i] = thisNeighbours;
    if(fVerbose > 2)
      cout << "\tTube " << tubeID_i << " has " << thisNeighbours.size() << " neighbours" << endl;
  }//i
}

void WCSimEBonsai::GetNEff()
{
  int nearbyHits, tubeID;
  float ratio, occupancy, lateHits, darkNoise, photoCoverage, waterTransparency;
  float nEffHit, nEffHit2, x, y, z;
  vector<int> neighbours;
  WCSimRootPMT pmt;
  fNEff = 0;
  fNEff2 = 0;

  for (unsigned int i = 0; i < fNMaxShort; i++) {
    tubeID = fTubeIdsShort[i];
    pmt = fGeom->GetPMT(tubeID - 1);
    x = pmt.GetPosition(0);
    y = pmt.GetPosition(1);
    z = pmt.GetPosition(2);

    // Calculate occupancy to correct for multiple hits on a single PMT
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

    // correct for delayed hits (e.g. due to scattering)
    lateHits = (fNMaxLong - fNMaxShort - (fNWorkingPMTs * fDarkRate * fShortDuration)) / fNMaxShort;

    // substract dark noise hits
    darkNoise = (fNWorkingPMTs * fDarkRate * fShortDuration) / fNMaxShort;

    // calculate effective coverage to correct for photoCoverage
    // this depends on angle of incidence, see Fig. 4.5 (left) of http://www-sk.icrr.u-tokyo.ac.jp/sk/_pdf/articles/2016/doc_thesis_naknao.pdf
    // TODO: take into account azimuthal dependence as show in Fig. 4.5 (right); currently assumes phi=0
    float dotProduct = pmt.GetOrientation(0)*(fVertex[0] - x) +
      pmt.GetOrientation(1)*(fVertex[1] - y) +
      pmt.GetOrientation(2)*(fVertex[2] - z);
    float theta = acos( dotProduct / fDistanceShort[i]) * 180 / 3.14159265;
    if (theta > 89.99) theta = 0; // we have apparently mis-reconstructed the vertex, so let's set ...
    if (theta < 0) theta = 0; // ... the coverage to the most likely value of 0.4 (i.e. theta < 40 degrees)
    
    photoCoverage = 1 / WCSimEBonsai::fEffCoverages[int(theta/10)];
    if (fDetector == kHyperK20)
      photoCoverage *= 38448/float(19462); // ratio of number of PMTs is not exactly 2

    // correct for scattering in water
    waterTransparency = exp(fDistanceShort[i] / fLambdaEff);

    nEffHit = (occupancy + lateHits - darkNoise) * photoCoverage * waterTransparency;
    fNEff += nEffHit;
    // ad-hoc modification for low photocoverage
    nEffHit2 = (pow(occupancy, 1.4) + lateHits - darkNoise) * photoCoverage * waterTransparency;
    fNEff2 += nEffHit2;
    
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


 // correct for dead PMTs; convert nWorkingPMTs to float because integer division is inaccurate
  fNEff  *= fNPMTs / float(fNWorkingPMTs);
  fNEff2 *= fNPMTs / float(fNWorkingPMTs);

  if(fVerbose) {
    std::cout << endl << "***************************************" << endl
	      << "Average nEff for all hits: " << fNEff << endl
	      << " (nEff2 for low photo-coverage is: " << fNEff2 << ")" << endl;
  }
}
void WCSimEBonsai::CorrectEnergy()
{
  // reconstruct energy from fNEff; this is approximately linear, except at very low energies
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
      // use fNEff2 (with occupancy to power of 1.4)
      fERec = 0.000001148*pow(fNEff2, 2) + 0.02032*fNEff2 + 1.94;
    break;
  }
  if(fVerbose) {
    std::cout << "Reconstructed energy = " << fERec << " MeV" << std::endl;
  }
}
