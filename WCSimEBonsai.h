#ifndef __WCSIMEBONSAI_H__
#define __WCSIMEBONSAI_H__

#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "WCSimRootGeom.hh"

using std::vector;

// ****************************************
// ** Energy estimation for a trigger    **
// ****************************************
// For a detailed description of the energy estimation formulas in SK-IV, see ch. 4.3 of
// http://www-sk.icrr.u-tokyo.ac.jp/sk/_pdf/articles/2016/doc_thesis_naknao.pdf

class WCSimEBonsai {

 public:
  WCSimEBonsai(const char * detectorname, WCSimRootGeom * geom, int verbose);
  ~WCSimEBonsai() {};

  double GetEnergy(std::vector<float> times, std::vector<int> tubeIds, float * vertex);

  //override default values with these methods
  void SetDarkRate(double darkrate);
  void SetNPMTs   (int npmts);
  void SetNWorkingPMTs(int nworkingpmts);
  void SetNeighbourDistance(double neighbourdistance);
  void SetShortDuration(double shortduration);
  void SetLongDuration(double longduration);
  void SetTopBottomDistance(double hi, double lo);

 private:
  enum kDetector_t {kSuperK = 0, kHyperK40, kHyperK20};

  kDetector_t DetectorEnumFromString(std::string name);
  void CorrectHitTimes();
  void FindMaxTimeInterval();
  void GetNearestNeighbours();
  void GetNEff();
  void CorrectEnergy();

  const double fLightGroupSpeed;
  const double fLambdaEff;
  static const double fEffCoverages[9];

  std::string fDetectorName;
  kDetector_t fDetector;
  double    fDarkRate;
  int       fNPMTs;
  int       fNWorkingPMTs;
  double    fNeighbourDistance;
  double    fTopBottomDistanceHi;
  double    fTopBottomDistanceLo;
  int       fVerbose;

  WCSimRootGeom * fGeom;

  int    fNDigits;
  int    fNMaxShort;
  int    fNMaxLong;
  double fLongDuration;
  double fShortDuration;
  double fStartTime;
  float  fVertex[4];

  double fNEff;
  double fNEff2;
  double fERec;

  vector<int>    fTubeIds;
  vector<double> fDistance;
  vector<float>  fTimes;
  vector<float>  fTimesCorrected;
  vector<float>  fTimesCorrectedSorted;

  vector<double> fDistanceShort;
  vector<double> fTubeIdsShort;

  std::map<int, std::vector<int> > fNeighbours;
};

#endif //__WCSIMEBONSAI_H__
