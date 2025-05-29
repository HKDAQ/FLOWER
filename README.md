> [!CAUTION]
> This repo has moved to HK gitlab at https://git.hyperk.org/hyperk/recon/FLOWER/
>
> Please clone this package from there.
>
> If you already have cloned this package, please see the migration instructions [here](https://git.hyperk.org/hyperk/wiki/-/wikis/FromGitHub#users--developers)
>
> Please also put new issues & merge requests in only on the HK gitlab


# FLOWER: soFtware for LOW-Energy Reconstruction
An energy reconstruction script for low-energy events in Hyper-Kamiokande.

FLOWER requires a trigger with a reconstructed event vertex from an external tool, such as (hk-)BONSAI. It then iterates over detected hits in a short time window and applies corrections for detector effects to estimate the total number of Cherenkov photons emitted by the particle. Finally, it matches this number to a reconstructed particle energy based on MC simulations.
A detailed description of this method (which is based on the one used in Super-Kamiokande) and the individual corrections can be found in chapter 3.4.3 of [Jost Migenda’s PhD thesis](https://inspirehep.net/literature/1778697).

## Compiling `WCSimFlower`
Make sure the following are sourced
* WCSim (`$WCSIM_BUILD_DIR` set)
  Note that the WCSim dependence of the `WCSimFlower` class is confined to geometry information.

Then
```bash
export FLOWERDIR=/path/to/FLOWER
make
```

## Running the FLOWER/BONSAI script bundled here
Make sure the following software is setup correctly
* WCSim (`$WCSIM_BUILD_DIR` set)
* BONSAI (`$BONSAIDIR` set)
* FLOWER (`$FLOWERDIR` set)

Then
```bash
export PATH=$FLOWERDIR/rootflower:$PATH
export FLOWERDATADIR=/path/to/writable/directory/ #e.g. $FLOWERDIR/data/
rootflower -b -q flower_with_bonsai.C+'("/path/to/wcsim/file.root",1,"SuperK")'
```

Note that `$FLOWERDATADIR` is used to store information that takes a while to calculate, and that only needs to be done once per geometry (e.g. the nearest neighbours of each PMT). If `$FLOWERDATADIR` is not set, or set to a directory that doesn't exist, `$FLOWERDIR/data/` will be used instead

## Running in other code
* Construct a class member, giving it a `detectorname` (it sets default values and determines the exact energy correction factors) and detector geometry
```
WCSimFlower(const char * detectorname, WCSimRootGeom * geom, int verbose);
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
