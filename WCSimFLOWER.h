#ifndef __WCSIMFLOWER_H__
#define __WCSIMFLOWER_H__

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

class WCSimFLOWER {

 public:
  WCSimFLOWER(const char * detectorname, WCSimRootGeom * geom, bool overwrite_nearest, int verbose);
  ~WCSimFLOWER() {};

  float GetEnergy(std::vector<float> times, std::vector<int> tubeIds, float * vertex);

  //override default values with these methods
  void SetDarkRate(float darkrate);
  void SetDarkRate2(float darkrate);
  void SetNPMTs   (int npmts);
  void SetNPMTs2  (int npmts);
  void SetNWorkingPMTs(int nworkingpmts);
  void SetNWorkingPMTs2(int nworkingpmts);
  void SetNeighbourDistance(float neighbourdistance, bool overwrite_nearest);
  void SetShortDuration(float shortduration);
  void SetLongDuration(float longduration);
  void SetTopBottomDistance(float hi, float lo);

  TString GetFLOWERDataDir();

 private:
  enum kDetector_t {kSuperK = 0, kHyperK40, kHyperK20, kHyperK20BnL10mPMT};

  kDetector_t DetectorEnumFromString(std::string name);
  void CorrectHitTimes();
  void FindMaxTimeInterval();
  void GetNearestNeighbours(bool overwrite_root_file);
  void GetNEff();
  void CorrectEnergy();

  const float fLightGroupSpeed;
  const float fLambdaEff;
  static const float fEffCoverages[9];

  std::string fDetectorName;
  kDetector_t fDetector;
  float    fDarkRate;
  float    fDarkRate2;
  int      fNallPMTs;
  int       fNPMTs;
  int       fNPMTs2;
  int       fNWorkingPMTs;
  int       fNWorkingPMTs2;
  float    fNeighbourDistance;
  float    fTopBottomDistanceHi;
  float    fTopBottomDistanceLo;
  int       fVerbose;

  WCSimRootGeom * fGeom;

  int    fNDigits;
  int    fNMaxShort;
  int    fNMaxLong;
  float fLongDuration;
  float fShortDuration;
  float fStartTime;
  float  fVertex[3];

  float fNEff;
  float fNEffMod;
  float fERec;

  vector<int>    fTubeIds;
  vector<float> fDistance;
  vector<float>  fTimes;
  vector<float>  fTimesCorrected;
  vector<float>  fTimesCorrectedSorted;

  vector<float> fDistanceShort;
  vector<float> fTubeIdsShort;

  std::map<int, std::vector<int> > fNeighbours;
};

#endif //__WCSIMFLOWER_H__
