# energetic-bonsai
Low-energy reconstruction script that works with (hk-)BONSAI

It also works with other low-energy reconstruction tools - it just requires a positional vertex

## Compiling
Make sure the following are sourced
* WCSim (`$WCSIMDIR` set)

Then
```bash
export EBONSAIDIR=/path/to/energetic-bonsai
make
```

## Running the WCSim/BONSAI script bundled here
Make sure the following are sourced
* WCSim (`$WCSIMDIR` set)
* BONSAI (`$BONSAIDIR` set)
* energetic_bonsai (`$EBONSAIDIR` set)

Then
```bash
export EBONSAIDIR=/path/to/energetic-bonsai
export PATH=$EBONSAIDIR:$PATH
rootebonsai -b -q energetic_bonsai.C+'("/path/to/wcsim/file.root",1,"SuperK")'
```

## Running in other code
* Construct a class member, giving it a detectorname (it sets default values and determines the exact energy correction factors) and detector geometry
```
WCSimEBonsai(const char * detectorname, WCSimRootGeom * geom, int verbose);
```
* Override default values of parameters with these methods
```
void SetDarkRate(double darkrate);
void SetNPMTs   (int npmts);
void SetNWorkingPMTs(int nworkingpmts);
void SetNeighbourDistance(double neighbourdistance);
void SetShortDuration(double shortduration);
void SetLongDuration(double longduration);
void SetTopBottomDistance(double hi, double lo);
```
* Every event/trigger, send it a list of hit times, hit tube IDs, and reconstructed vertex position (x,y,z)
```
double GetEnergy(std::vector<float> times, std::vector<int> tubeIds, float * vertex);
```