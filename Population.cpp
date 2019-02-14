/*
 * Population.cpp
 *
 * A network of "cells".  Each cell can have multiple Communities.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "GridCells.h"
#include "Population.h"
#include "R0Community.h"

using namespace std;

Population::Population() {
  _community = NULL;
  _nCommunityStart = NULL;
  _nNumCommunities = 0;
  _grid = NULL;
  _nGridSizeX=0;
  _nGridSizeY=0;
  _riverDist = NULL;
  _riverNeighbors = NULL;
  _newVibrio = NULL;
  _riverVibrio = _riverHyperVibrio = _riverVibrio2 = _riverHyperVibrio2 = NULL;
  _szLabels = _szLabels2 = NULL;
  _nNumUniqueLabels = 0;
  _fWorkingFraction = 0.0;
  _fXOrigin = _fYOrigin = 0;
  _nRiverShedCycles = 0;
  _fRiverShedFraction = 0.0;
  _fRiverFlowDelta = 0.0;
  _fTravelProb = 0.0;
  _fVaccinationTarget = 0.0;
  _nVaccinationThreshold = 0;
  _nVaccinationDelay = -1;
  _bGridVaccinated = NULL;
  _nVaccinationDay = NULL;
  _nNumVaccinesUsed = 0;
  _nNumVaccinesAvailable = 100000000;
  _fHygieneTarget = 0.0;
  _nTargetCommunitySize = 2000;
  _nMaxCommunitiesOnRiver = 6;
  _fDriveProb = 0.00001;
  _bHighway = NULL;
  _nHighwayDests = NULL;
  _nNumHighwayDests = NULL;
  _nHighwayDepth = NULL;
  _nDay = -1;
  _nDayStartOffset = 0;
  _nDayToAge = -1;
  _nNumRainfallLabels = 0;
  _fRainfallData=NULL;
  _fRainSheddingMultipliers[0] = _fRainSheddingMultipliers[1] = _fRainSheddingMultipliers[2] = _fRainSheddingMultipliers[3] = 1.0;
  _nDaysToWane = 3*365;
}

Population::~Population() {
  if (_nHighwayDests) {
    for (int i=0; i<_grid->getSize(); i++)
      delete [] _nHighwayDests[i];
    delete [] _nHighwayDests;
  }
  if (_community)
    delete [] _community;
  if (_nCommunityStart)
    delete [] _nCommunityStart;
  if (_grid)
    delete _grid;
  if (_riverDist)
    delete [] _riverDist;
  if (_riverVibrio)
    delete [] _riverVibrio;
  if (_riverHyperVibrio)
    delete [] _riverHyperVibrio;
  if (_riverVibrio2)
    delete [] _riverVibrio2;
  if (_riverHyperVibrio2)
    delete [] _riverHyperVibrio2;
  if (_riverNeighbors)
    delete [] _riverNeighbors;
  if (_newVibrio)
    delete [] _newVibrio;
  if (_bHighway)
    delete [] _bHighway;
  if (_nHighwayDepth)
    delete [] _nHighwayDepth;
  if (_nNumHighwayDests)
    delete [] _nNumHighwayDests;
  if (_bGridVaccinated)
    delete [] _bGridVaccinated;
  if (_nVaccinationDay)
    delete [] _nVaccinationDay;
  if (_szLabels)
    delete [] _szLabels;
  if (_szLabels2)
    delete [] _szLabels2;
  if (_fRainfallData)
    delete [] _fRainfallData;
}

int Population::loadGrid(gsl_rng *rng, int nGridSize, int nCommunitySize) {
  //  _nGridSize = nGridSize;
  _nGridSizeX = _nGridSizeY = nGridSize;
  _grid = new HexGridCells(_nGridSizeX,_nGridSizeY,false,false);
  _community = new Community[_grid->getSize()]; // populations of people
  _newVibrio = new double[_grid->getSize()];

  // initialize grid and non-working residents
  initGridPopulation(rng, nCommunitySize);
  return 1;
}

// Initialize a population using line list file of ages and households (popfilename).
// The population will be on a rectangular grid based on the coordinates in dimfilename.
int Population::loadPopulation(gsl_rng *rng,
			       double workingfraction,
			       double tau1, double tau2, double rho,
			       const char *dimfilename,
			       const char *popfilename,
			       bool bR0run) {
  const int MAXCELLPOP=50000;
  const int MAXCOMMUNITIES=41000;
  // get size of population grid
  ifstream issdim(dimfilename);
  if (!issdim) {
    cerr << "ERROR: " << dimfilename << " not found." << endl;
    return -1;
  }
  issdim >> _fXOrigin >> _fXMax >> _nGridSizeX;
  issdim >> _fYOrigin >> _fYMax >> _nGridSizeY;

  cerr << "Loading " << popfilename << endl;
  _fWorkingFraction = workingfraction;
  ostringstream oss;
  if (popfilename)
    oss.str(popfilename);
  else
    return -1;

  _grid = new SquareGridCells(_nGridSizeX,_nGridSizeY,false,false);
  _nCommunityStart = new int[_grid->getSize()+1]; // index of first community per raster cell
  _nCommunities = new int[_grid->getSize()+1]; // number of communities per raster cell
  _szLabels = new string[_grid->getSize()];
  _szLabels2 = new string[_grid->getSize()];
  _newVibrio = new double[_grid->getSize()];

  for (int i=0; i<=_grid->getSize(); i++) {
    _nCommunityStart[i] = -1;
    _nCommunities[i] = 0;
  }
  _nNumCommunities = 0;

  // read population file
  ifstream iss(popfilename);
  if (!iss) {
    cerr << "ERROR: " << popfilename << " not found." << endl;
    return -1;
  }
  // load population data and populate communities
  if (bR0run)
    _community = new R0Community[MAXCOMMUNITIES];
  else 
    _community = new Community[MAXCOMMUNITIES]; // well-mixed groups of people
  string line;
  getline(iss, line); // throw away header line
  int ages[MAXCELLPOP]; // ages read from raster cell
  int householdids[MAXCELLPOP]; // households read from raster cell
  int cellpop = 0; // population in this raster cell
  int totalpop = 0; // total people read so far
  double lastlongitude=-1,lastlatitude=-1; // assumes all individuals in a given raster cell are adjacent in the file
  string lastdept;
  do {
    int personid;
    string dept;
    double longitude,latitude;
    string department;
    // read in people until there is a new lat/long
    do {
      if (iss >> personid >> householdids[cellpop] >> longitude >> latitude >> dept >> ages[cellpop]) { // read data for one individual
	cellpop++;
	totalpop++;
	if (cellpop>=MAXCELLPOP) {
	  cerr << "ERROR: Too many people in a cell " << longitude << "," << latitude << "," << dept << " : " << cellpop << endl;
	  exit(-1);
	}
      } else
	break; // end of file? bad input?
    } while (lastlongitude==longitude && lastlatitude==latitude && !iss.eof()); // grid cell is done

    {
      if (!(lastlongitude==longitude && lastlatitude==latitude))
	cellpop--; // last person does not belong in this cell
      if (cellpop>0) {
	int xloc = (int)(round((_nGridSizeX-1)*(lastlongitude-_fXOrigin)/(_fXMax-_fXOrigin)));
	int yloc = (int)(round((_nGridSizeY-1)*(lastlatitude-_fYOrigin)/(_fYMax-_fYOrigin)));
	if (xloc>=0 && xloc<_nGridSizeX &&
	    yloc>=0 && yloc<_nGridSizeY) {
	  int gridindex = (int)(round(yloc*_nGridSizeX+xloc));
	  assert (gridindex<_grid->getSize());
	  _szLabels[gridindex] = lastdept;
	  bool bFound=false;
	  for (int j=0; j<_nNumUniqueLabels; j++)
	    if (lastdept.compare(_szUniqueLabels[j])==0) {
	      bFound=true;
	      break;
	    }
	  if (!bFound)
	    _szUniqueLabels[_nNumUniqueLabels++] = lastdept;

	  // create new communities by dividing cell population among them
	  int numcomms = round(cellpop/_nTargetCommunitySize);
	  int popleft = cellpop; // how many people are left to allocate
	  if (numcomms==0)
	    numcomms=1;
	  _nCommunityStart[gridindex]=_nNumCommunities;
	  _nCommunities[gridindex]=numcomms;
	  for (int i=0; i<numcomms; i++) {
	    //	    cout << "grid id = " << gridindex << " : " << _nCommunityStart[gridindex]+i << "(" << i << ") : " << round(popleft/(numcomms-i)) << " , " << numcomms << "/" << cellpop << " : " << popleft << endl;
	    Community &comm = _community[_nNumCommunities];
	    //	    cerr << " grid " << gridindex << " ; " << " comm " << i << " " << _nNumCommunities << ", " << comm.getid() << endl;
	    comm.populate(round(popleft/(numcomms-i)), popleft, householdids+cellpop-popleft, ages+cellpop-popleft);
	    popleft-=comm.getNumResidents();
	    _nNumCommunities++;
	    //	    cout << "comm " << gridindex << "," << _nNumCommunities << "," << cellpop << "," << popleft << endl;
	    if (_nNumCommunities>=MAXCOMMUNITIES) {
	      cerr << "ERROR: Too many communities" << _nNumCommunities << endl;
	      assert(-1);
	    }
	  }
	} else {
	  cerr << "ERROR: out of bounds: " << longitude <<"," << latitude <<"," << dept << " : " << xloc << "," << yloc << endl;
	}
	//      }
      //      cout << xloc << "," << yloc <<"," << ages[index] << "," << dept << endl;

      }
      if (!(lastlongitude==longitude && lastlatitude==latitude)) {
        // put last person read in the first slot for the next cell
	ages[0] = ages[cellpop];
	householdids[0] = householdids[cellpop];
	cellpop = 1;
      } else {
	cellpop = 0;
      }
      lastlongitude=longitude;
      lastlatitude=latitude;
      lastdept=dept;
    }
  } while (!iss.eof());
  iss.close();
  cerr << "DONE reading and populating " << totalpop << " people and " << _nNumCommunities << " communities" << endl;
  _bGridVaccinated = new bool[_grid->getSize()];
  for (int i=0; i<_grid->getSize(); i++)
    _bGridVaccinated[i] = false;
  _nVaccinationDay = new int[_grid->getSize()];
  for (int i=0; i<_grid->getSize(); i++)
    _nVaccinationDay[i] = 10000000;

  // assign people workplaces
  for (int sourceindex=0; sourceindex<_grid->getSize(); sourceindex++)
    assignWorkPlaces(rng, sourceindex, _fWorkingFraction,
		     tau1, tau2, rho);
  cerr << "DONE assigning workplaces (" << (100*_fWorkingFraction) << "% employment of working age)" << endl;
  return 1;
}

// assignWorkPlaces - assigns workplaces for residents of cell nSourceID
// returns number of employed people
int Population::assignWorkPlaces(gsl_rng *rng, int nSourceID, 
				 double fWorkFraction,
				 double tau1, double tau2, double rho) {
  if (_nCommunityStart[nSourceID]<0 || _nCommunities[nSourceID]<=0)
    return 0;
  const int maxradius = 40;
  if (fWorkFraction<=0.0)
    return 0;
  int P1 = getNumResidents(nSourceID); // population in this cell
  if (P1==0)
    return 0;
  double sourcex = _grid->getX(nSourceID);
  double sourcey = _grid->getY(nSourceID);

  // calculate probabilities of working in other neighborhoods
  int neighbors[maxradius*maxradius*4]; // id of neighbors
  double probs[maxradius*maxradius*4];  // probabilities based on gravity model
  int totalneighbors = 1;
  for (int destindex=0; destindex<_grid->getSize(); destindex++) {
    double destx = _grid->getX(destindex);
    double desty = _grid->getY(destindex);
    if (destindex!=nSourceID && fabs(destx-sourcex)<maxradius && fabs(desty-sourcey)<maxradius) {
      double dist = sqrt((sourcex-destx)*(sourcex-destx) + (sourcey-desty)*(sourcey-desty))*fabs((_fXMax-_fXOrigin)/_nGridSizeX*111.0); // in kilometers, assuming 1 degree=111km
      int P2 = getNumResidents(destindex);
      if (P2>0 && dist<maxradius) {
	neighbors[totalneighbors] = destindex;
	probs[totalneighbors] = pow(P1, tau1) * pow(P2, tau2) / pow(dist, rho); // do we need theta * ???
	//	      cerr << "dist=" << dist << "," << probs[totalneighbors] <<  "," << 		pow(P1, tau1)  << ","  << pow(P2, tau2) << ","  << (pow(dist, rho)) << endl;
	totalneighbors++;
	assert(totalneighbors<maxradius*maxradius*4);
      }
    }
  }
  if (totalneighbors==1) {
    return 0;
  } else {
    int count = 0;
    neighbors[0] = nSourceID;
    probs[0] = 0.3; // 30% chance of working in home community
    double sum=0.0;
    for (int i=1; i<totalneighbors; i++)
      sum += probs[i];
    sum/=(1.0-probs[0]); // normalize over neighbors >=1 cell away
    for (int i=1; i<totalneighbors; i++)
      probs[i] = probs[i-1]+probs[i]/sum;
    assert(fabs(1.0-probs[totalneighbors-1])<0.000001);
	   
    for (int commindex=_nCommunityStart[nSourceID]; commindex<_nCommunityStart[nSourceID]+_nCommunities[nSourceID]; commindex++) {
      Community &sourcecomm = _community[commindex];
      for (int i=0; i<sourcecomm.getNumResidents(); i++) {
	int id = sourcecomm.getResidentID(i);
	Person &person = Person::personArray[id];
	if (person.getHomeCommunity()!=commindex) {
	  cerr << person.getHomeCommunity() << " != " << commindex << " in " << nSourceID << " " << _nCommunityStart[nSourceID] << "+" << _nCommunities[nSourceID] << endl;
	  cerr << " " << _nCommunityStart[nSourceID-1] << "+" << _nCommunities[nSourceID-1] << endl;
	}
	assert(person.getHomeCommunity()==commindex);
	if (person.getAge()>=15 && gsl_rng_uniform(rng)<fWorkFraction) { // this person works
	  double r = gsl_rng_uniform(rng);
	  for (int j=0; j<totalneighbors; j++)
	    if (probs[j]>=r) {
	      int dest = _nCommunityStart[neighbors[j]];
	      int diff = _nCommunities[neighbors[j]];
	      if (diff>1)
		dest += gsl_rng_uniform_int(rng, diff);
	      person.setWorkCommunity(dest);
	      _community[dest].addVisitor(id);
	      count++;
	      break;
	    }
	}
      }
    }
    return count;
  }
}

int Population::loadRivers(const char *riverfilename) {
  ostringstream oss;
  if (riverfilename)
    oss.str(riverfilename);
  else
    return -1;
  ifstream iss(oss.str().c_str());
  if (!iss) {
    cerr << "ERROR: " << riverfilename << " not found." << endl;
    return -1;
  }
  _riverDist = new int[_grid->getSize()];
  _riverNeighbors = new int[_grid->getSize()];
  _riverVibrio = new double[_grid->getSize()];
  _riverHyperVibrio = new double[_grid->getSize()];
  _riverVibrio2 = new double[_grid->getSize()];
  _riverHyperVibrio2 = new double[_grid->getSize()];
  for (int i=0; i<_grid->getSize(); i++) {
    _riverDist[i] = -1;
    _riverVibrio[i] = 0.0;
    _riverHyperVibrio[i] = 0.0;
  }

  cerr << "loading rivers " << _nGridSizeX << " x " << _nGridSizeY << endl;
  string line;
  getline(iss, line); // throw away header line
  int nNumLinesRead = 0;
  while (getline(iss, line)) {
    istringstream iss;
    iss.str(line);
    int xdeg,xmin,xsec,ydeg,ymin,ysec,dist;
    iss >> ydeg >> ymin >> ysec >> xdeg >> xmin >> xsec >> dist;
    nNumLinesRead++;
    double longitude=-(xdeg+xmin/60.0+xsec/60.0/60.0); // river file assumes W longitude instead of negative
    double latitude=ydeg+ymin/60.0+ysec/60.0/60.0;
    int xloc = (int)(round((_nGridSizeX-1)*(longitude-_fXOrigin)/(_fXMax-_fXOrigin)));
    int yloc = (int)(round((_nGridSizeY-1)*(latitude-_fYOrigin)/(_fYMax-_fYOrigin)));
    if (dist>0 &&
	xloc>=0 && xloc<_nGridSizeX &&
	yloc>=0 && yloc<_nGridSizeY) {
      //      cerr << "river x,y loc: " << xloc << "," << yloc << "," << dist << endl;
      int gridindex = xloc + yloc*_nGridSizeX;
      assert (gridindex<_grid->getSize());
      _riverDist[gridindex] = dist;
      for (int commindex=_nCommunityStart[gridindex]; commindex<_nCommunityStart[gridindex]+_nCommunities[gridindex]; commindex++)
	_community[commindex].setRiver(true);
    }
  }
  // count cells that are downstream from each cell
  for (int i=0; i < _grid->getSize(); i++) {
    _riverNeighbors[i] = 0;
    if (_riverDist[i]>0) {
      for (int x=-1; x<=1; x++)
	for (int y=-1; y<=1; y++)
	  if ((x!=0 || y!=0) &&
	      _riverDist[i+x+y*_nGridSizeX] > _riverDist[i])
	    _riverNeighbors[i]++;
    }
  }

  return nNumLinesRead;
}

int Population::loadHighways(const char *highwayfilename) {
  ostringstream oss;
  if (highwayfilename)
    oss.str(highwayfilename);
  else
    return -1;
  ifstream iss(oss.str().c_str());
  if (!iss) {
    cerr << "ERROR: " << highwayfilename << " not found." << endl;
    return -1;
  }
  _bHighway = new bool[_grid->getSize()];
  for (int i=0; i<_grid->getSize(); i++)
    _bHighway[i] = false;

  cerr << "loading highways " << _nGridSizeX << " x " << _nGridSizeY << endl;
  string line;
  getline(iss, line); // throw away header line
  int nNumLinesRead = 0;
  while (getline(iss, line)) {
    istringstream iss;
    iss.str(line);
    int xdeg,xmin,xsec,ydeg,ymin,ysec,dist;
    iss >> ydeg >> ymin >> ysec >> xdeg >> xmin >> xsec >> dist;
    nNumLinesRead++;
    double longitude=-(xdeg+xmin/60.0+xsec/60.0/60.0); // highway file assumes W longitude instead of negative
    double latitude=ydeg+ymin/60.0+ysec/60.0/60.0;
    int xloc = (int)(round((_nGridSizeX-1)*(longitude-_fXOrigin)/(_fXMax-_fXOrigin)));
    int yloc = (int)(round((_nGridSizeY-1)*(latitude-_fYOrigin)/(_fYMax-_fYOrigin)));
    if (dist>0 &&
	xloc>=0 && xloc<_nGridSizeX &&
	yloc>=0 && yloc<_nGridSizeY) {
      int gridindex = xloc + yloc*_nGridSizeX; // assumes rectangular lattice
      assert (gridindex<_grid->getSize());
      _bHighway[gridindex] = true;
    }
  }
  return nNumLinesRead;
}

int Population::loadRainfall(const char *rainfallfilename) {
  ostringstream oss;
  if (rainfallfilename)
    oss.str(rainfallfilename);
  else
    return -1;
  ifstream iss(oss.str().c_str());
  if (!iss) {
    cerr << "ERROR: " << rainfallfilename << " not found." << endl;
    return -1;
  }
  cerr << "Loading rainfall from " << rainfallfilename << endl;
  string line;
  getline(iss, line); // header line
  istringstream linestream(line);
  _nNumRainfallLabels = 0;
  string subs;
  while (linestream >> subs) {
    if (_nNumRainfallLabels>0) { // ignore first label
      _szRainfallLabels[_nNumRainfallLabels-1] = subs;
    }
    _nNumRainfallLabels++;
    if (_nNumRainfallLabels>=MAXRAINFALLLABELS) {
      cerr << "ERROR: Too many rainfall columns" << endl;
      exit(-1);
    }
  }
  _nNumRainfallLabels--;
  cerr << _nNumRainfallLabels << " rainfall locations" << endl;
  const int MAXRAINFALLDAYS=8000; // maximum number of days after Jan 1 2010 for rainfall data
  _fRainfallData = new double[_nNumRainfallLabels*MAXRAINFALLDAYS];
  memset(_fRainfallData, 0, _nNumRainfallLabels*MAXRAINFALLDAYS*sizeof(double));
  _nNumRainfallDays = 0;
  while (getline(iss, line)) {
    istringstream linestream(line);
    int rainday;
    linestream >> rainday;
    if (rainday>=MAXRAINFALLDAYS) {
      cerr << "ERROR: too many rainfall days: " << rainday << endl;
    }
    if (rainday>_nNumRainfallDays)
      _nNumRainfallDays=rainday;
    int count=0; // which location
    double d;
    while (linestream >> d) {
      _fRainfallData[rainday*_nNumRainfallLabels+count] = d;
      count++;
    }
  }
  return 1;
}

// population is on a uniform lattice
int Population::initGridPopulation(gsl_rng *rng, int nCommunitySize) {
  // initialize grid and non-working residents
  int maxradius = min(15,_nGridSizeX);
  double nnorm = 0.0;
  double fade = 0.75;
  for (int i=0; i<maxradius; i++)
    nnorm += pow(fade,i);

  for (int i=0; i<_grid->getSize(); i++) {
    int neighbors[1000];
    double probs[1000];
    int pop = nCommunitySize;
    // figure out where people in this cell work
    if (pop>0 && _fWorkingFraction>0.0) {
      neighbors[0] = i;
      probs[0] = 1.0/nnorm;
      int totalneighbors = 1;
      for (int radius=1; radius<15; radius++) {
	int numneighbors = getRadius(_grid, i, radius, neighbors+totalneighbors, 1000);
	// remove neighbors with no residents
	for (int j=numneighbors-1; j>=0; j--) {
	  if (getNumResidents(neighbors[totalneighbors+j])<=0) {
	    neighbors[totalneighbors+j] = neighbors[totalneighbors+numneighbors-1]; // copy last valid neighbor to this spot
	    numneighbors--;
	  }
	}
	for (int j=0; j<numneighbors; j++)
	  probs[totalneighbors+j] = pow(fade,radius)/nnorm/numneighbors;
	totalneighbors+=numneighbors;
	assert(totalneighbors<1000);
      }
      for (int j=1; j<totalneighbors; j++)
	probs[j] = probs[j-1]+probs[j];
      _community[i].populate(rng, pop);
    } else {
      _community[i].populate(rng, pop);
    }
  }

  for (int sourceindex=0; sourceindex<_grid->getSize(); sourceindex++)
    assignWorkPlaces(rng, sourceindex, _fWorkingFraction,
		     0.0,0.0,2); //tau1, tau2, rho); // fix this!!!!!!!

  return 1;
}

int Population::getRadius(GridCells *g, int center, int radius, int *buf, int bufsize) {
  const int maxneighbors = 6;
  int searched[4000];
  int searchedsize=1;
  int oldsearchedsize=1;
  searched[0] = center;
  // breadth first search
  for (int r=0; r<radius; r++) {
    oldsearchedsize=searchedsize;
    for (int i=searchedsize-1; i>=0; i--) {
      int neighbors[maxneighbors];
      int numneighbors = _grid->getNeighbors(searched[i], neighbors, maxneighbors);
      for (int j=0; j<numneighbors; j++) {
	bool bFound=false;
	for (int k=0; k<searchedsize && !bFound; k++) {
	  if (searched[k]==neighbors[j])
	    bFound=true;
	}
	if (!bFound)
	  searched[searchedsize++] = neighbors[j];
	assert(searchedsize<1200);
      }
    }
  }
  for (int i=0; i<searchedsize-oldsearchedsize; i++)
    buf[i] = searched[oldsearchedsize+i];
  return searchedsize-oldsearchedsize;
}

// each degree is broken into 120 intervals
//int Population::convertDegToInt(int deg, int min, int sec) {
//  return deg*120+min*2+(sec>30?1:0);
//}

// infect - infect frac of the people in the cell gridindex
int Population::infect(gsl_rng *rng, int gridindex, double frac, bool bMakeSymptomatic) {
  int count = 0;
  for (int commindex=_nCommunityStart[gridindex];
       commindex < _nCommunityStart[gridindex]+_nCommunities[gridindex];
       commindex++) {
	 if (_community[commindex].getNumResidents()>0) {
	int nNumInfect = gsl_ran_binomial(rng, 
					   frac,
					  _community[commindex].getNumResidents());
	count += _community[commindex].infect(rng, nNumInfect, bMakeSymptomatic);
	 }
       }
  return count;
}

// infect - infect num people in the cell closest to x,y
int Population::infect(gsl_rng *rng, int num, double longitude, double latitude, bool bMakeSymptomatic) { 
  //	int xloc = (int)(round((_nGridSizeX-1)*(lastlongitude-_fXOrigin)/(_fXMax-_fXOrigin)));
  //	int yloc = (int)(round((_nGridSizeY-1)*(lastlatitude-_fYOrigin)/(_fYMax-_fYOrigin)));
  int index=0;
  double xindex, yindex;
  double distindex=1000000.0;
  for (int i=0; i<_grid->getSize(); i++) {
    double newx = _grid->getX(i);
    double newy = _grid->getY(i);
    double newlongitude = _fXOrigin+newx*(_fXMax-_fXOrigin)/(_nGridSizeX-1);
    double newlatitude = _fYOrigin+newy*(_fYMax-_fYOrigin)/(_nGridSizeY-1);
    //    double newdist = (newx-x)*(newx-x)+(newy-y)*(newy-y);
    double newdist = (newlongitude-longitude)*(newlongitude-longitude)+(newlatitude-latitude)*(newlatitude-latitude);
    if (newdist<distindex && getNumResidents(i)>0) {
      index = i;
      xindex = newx;
      yindex = newy;
      distindex = newdist;
    }
  }
  cerr << "seed infection at : " << latitude << "," << longitude << " = cell " << xindex << "," <<  yindex << endl;
  if (_community[_nCommunityStart[index]].getNumResidents()>0)
    _community[_nCommunityStart[index]].infect(rng, num, bMakeSymptomatic);
  return index;
}

// infectOne - infects a single person in a random community
int Population::infectOne(gsl_rng *rng, bool bMakeSymptomatic) {
  int commindex = gsl_rng_uniform_int(rng, _nNumCommunities);
  _community[commindex].infect(rng, 1, bMakeSymptomatic);
  return -1;
}

// prevaccinate
// we do not check how much vaccine is available
void Population::prevaccinate(gsl_rng *rng) {
  if (_fVaccinationTarget>0.0 && _bGridVaccinated) {
    for (int i=0; i < _grid->getSize(); i++)
      if (!_bGridVaccinated[i]) {
	for (int j=_nCommunityStart[i]; j<_nCommunityStart[i+1]; j++)
	  _nNumVaccinesUsed += _community[j].prevaccinate(rng, _fVaccinationTarget);
	_bGridVaccinated[i] = true;
	_nVaccinationDay[i] = -1;
      }
  }
}

void Population::prioritizeCell(int gridindex) {
  if (!_bGridVaccinated[gridindex])
    _nVaccinationDay[gridindex] = _nDay+_nDayStartOffset;
}

int Population::computedrivelength(gsl_rng *rng, int start) {
  const int MAXDIST=200;
  const int MAXQUEUESIZE=MAXDIST*30;
  int gridqueue[MAXQUEUESIZE+1];
  if (_nHighwayDepth==NULL)
    _nHighwayDepth = new int[_grid->getSize()];
  memset(_nHighwayDepth, 0, sizeof(int)*_grid->getSize());
  int head = 0;
  int tail = 1;
  gridqueue[head] = start;
  _nHighwayDepth[start] = 1;
  //  cerr << "start at " << start << ": " << (start%_nGridSizeX) << "," << (start/_nGridSizeX) << endl;
  while(head<tail && tail<MAXQUEUESIZE-1) {
    int loc = gridqueue[head++];
    for (int x=-1; x<=1; x++)
      for (int y=-1; y<=1; y++) {
	int pos = loc+_nGridSizeX*y+x;
	if (pos>=0 && pos<_grid->getSize())
	  if ((x!=0 || y!=0) &&
	      _bHighway[pos] &&
	      _nHighwayDepth[pos]<=0) {
	    //	    	  cerr << "  " << loc << ": " << (loc%_nGridSizeX)+x << "," << (loc/_nGridSizeX)+y << "," << (_nHighwayDepth[loc]+1) << endl;
	    if (_nHighwayDepth[loc]<MAXDIST) {
	      assert(tail<MAXQUEUESIZE);
	      if (tail>=MAXQUEUESIZE)
		break;
	      gridqueue[tail++] = loc+_nGridSizeX*y+x;
	      _nHighwayDepth[pos] = _nHighwayDepth[loc]+1;
	    }
	  }
      }
  }
  return start;
}

int Population::drive(gsl_rng *rng) {
  if (_fDriveProb>0.0)
    for (int i=0; i < _grid->getSize(); i++) {
      if (_bHighway[i] && (_nNumHighwayDests==NULL || _nNumHighwayDests[i]!=0) && getNumResidents(i)>0) {
	int nNumDrivers = gsl_ran_binomial(rng, 
					   _fDriveProb,
					   getNumResidents(i));
	if (nNumDrivers>0) {
	  if (!_nNumHighwayDests) {
	    _nNumHighwayDests = new int[_grid->getSize()];
	    for (int j=0; j < _grid->getSize(); j++)
	      _nNumHighwayDests[j] = -1;
	    _nHighwayDests = new int *[_grid->getSize()];
	  }
	  if (_nNumHighwayDests[i]<0) {
	    // store possible destinations
	    _nNumHighwayDests[i] = 0;
	    computedrivelength(rng, i);
	    _nHighwayDests[i] = new int[500];
	    for (int j=0; j<_grid->getSize(); j++)
	      if (_nHighwayDepth[j]>30 && getNumResidents(j)>2000) {
		_nHighwayDests[i][_nNumHighwayDests[i]++] = j;
		assert(_nNumHighwayDests[i]<500);
	      }
	  }
	  if (_nNumHighwayDests[i]>0 && _nHighwayDests[i])
	    while (--nNumDrivers >= 0) {
	      // swap non-symptomatic person at source with dest
	      int destnum = gsl_rng_uniform_int(rng, _nNumHighwayDests[i]);
	      int commsource = _nCommunityStart[i];
	      if (_nCommunities[i]>0)
		commsource += gsl_rng_uniform_int(rng, _nCommunities[i]);
	      int commdest = _nCommunityStart[_nHighwayDests[i][destnum]];
	      if (_nCommunities[_nHighwayDests[i][destnum]]>0)
		commdest += gsl_rng_uniform_int(rng, _nCommunities[_nHighwayDests[i][destnum]]);
	      int person1id = _community[commsource].getResidentID(gsl_rng_uniform_int(rng, _community[commsource].getNumResidents()));
	      int person2id = _community[commdest].getResidentID(gsl_rng_uniform_int(rng, _community[commdest].getNumResidents()));
	      Person &p1 = Person::personArray[person1id];
	      Person &p2 = Person::personArray[person2id];
	      if ((!p1.isSymptomatic() || p1.getDaysInfectious()==0) &&
		  (!p2.isSymptomatic() || p2.getDaysInfectious()==0)) { // don't travel if symptomatic for more than one day
		Person temp;
		temp.copyInfectionStatus(p1);
		p1.copyInfectionStatus(p2);
		p2.copyInfectionStatus(temp);
	      }
	      //	    cout << (start%_nGridSizeX) << "," << (start/_nGridSizeX) << "," << _nHighwayDepth[start] << endl;
	    }
	}
      }
    }
  return 0;
}

// age population by one year
// draw immunity status from earlier age cohort from the same community (aging)
// make those <1 year old fully susceptible (birth)
void Population::agePopulation(gsl_rng *rng) {
  const int AGEMAX = 100; // maximum age in years
  const int COHORTMAX = 3000; // largest number of people in one age bin
  cerr << "Aging population: simulation day " << _nDay << ", calendar day " << (_nDay+_nDayStartOffset) << endl;
  for (int gridnum=0; gridnum < _grid->getSize(); gridnum++) {
    for (int commnum=0; commnum<_nCommunities[gridnum]; commnum++) {
      // sort one-year age cohorts for each community
      Person *cohort[AGEMAX+1][COHORTMAX];
      int cohortsize[AGEMAX+1];
      memset(cohortsize, 0, (AGEMAX+1)*sizeof(int));
      for (int personnum=0; personnum<_community[_nCommunityStart[gridnum]+commnum].getNumResidents(); personnum++) {
	Person &p = (Person::personArray[_community[_nCommunityStart[gridnum]+commnum].getResidentID(personnum)]);
	assert(p.getAge()<=AGEMAX);
	assert(cohortsize[p.getAge()]+1<COHORTMAX);
	cohort[p.getAge()][cohortsize[p.getAge()]++]=&p;
      }
      // copy vaccine and immune status from younger age cohorts
      // aging doesn't happen if there are not enough people in the grid cell
      for (int agecohort=AGEMAX; agecohort>0; agecohort--) {
	if (cohortsize[agecohort]>0 && cohortsize[agecohort-1]>0) {
	  if (cohortsize[agecohort-1]==1) { // younger cohort is just one person
	    for (int pnum=0; pnum<cohortsize[agecohort]; pnum++)
	      cohort[agecohort][pnum]->copyImmuneStatus(*cohort[agecohort-1][0]); // copy status from younger cohort
	  } else  { // draw statuses from younger cohort
	    for (int pnum=0; pnum<cohortsize[agecohort]; pnum++)
	      cohort[agecohort][pnum]->copyImmuneStatus(*cohort[agecohort-1][gsl_rng_uniform_int(rng, cohortsize[agecohort-1])]); // copy status from younger cohort
	  }
	}
      }
      // reset youngest cohort to naive
      if (cohortsize[0]>0) {
	for (int pnum=0; pnum<cohortsize[0]; pnum++) {
	  cohort[0][pnum]->resetInfectionStatus();
	  cohort[0][pnum]->setTemporaryImmunityCountdown(gsl_rng_uniform_int(rng,365)); // newborn will not be exposed to any infection for 0 to 364 days.
	  // this prevents a synchronous introduction of susceptibles
	}
      }
    }
  }
}

int Population::step(gsl_rng *rng) {
  _nDay++;
  if ((_nDayToAge>=0) && ((_nDay+_nDayStartOffset)%365==_nDayToAge))
    agePopulation(rng);
  // shed vibrio and record amount shed today
  string lastlabel="----------"; // which location's rainfall data?
  int lastlabelmatch=-1; // which location's rainfall data?
  for (int i=0; i < _grid->getSize(); i++) {
    double rainmultiplier = 1.0;
    // daily rainfall data for this grid cell used to modify environmental shedding
    if (_fRainfallData!=NULL && (_nDay+_nDayStartOffset)<=_nNumRainfallDays) {
      string s = _szLabels[i];
      if (lastlabelmatch<0 || s.compare(lastlabel)!=0) { // different from previous cell's label
	lastlabelmatch=-1;
	lastlabel = _szLabels[i];
	for (int j=0; j<_nNumRainfallLabels; j++) {
	  if (_szLabels[i].compare(_szRainfallLabels[j])==0) {
	    lastlabelmatch = j;
	    break;
	  }
	}
      }
      if (lastlabelmatch>0) {
	double raintoday = _fRainfallData[(_nDay+_nDayStartOffset)*_nNumRainfallLabels+lastlabelmatch];
	if (raintoday==0.0) // no rain
	  rainmultiplier = _fRainSheddingMultipliers[0];
	else if (raintoday<1.0) // 0.1 mm
	  rainmultiplier = _fRainSheddingMultipliers[1];
	else if (raintoday<10.0) // 1 mm
	  rainmultiplier = _fRainSheddingMultipliers[2];
	else if (raintoday<50.0) // 5 mm
	  rainmultiplier = _fRainSheddingMultipliers[3];
	else // even more rain
	  rainmultiplier = _fRainSheddingMultipliers[4];
      } else
	  rainmultiplier = _fRainSheddingMultipliers[0];
    }
    _newVibrio[i] = 0;
    for (int j=0; j<_nCommunities[i]; j++)
      if (j<_nMaxCommunitiesOnRiver)
	_newVibrio[i] += _community[_nCommunityStart[i]+j].poop(rng, rainmultiplier);
      else 
	_community[_nCommunityStart[i]+j].poop(rng, rainmultiplier);
  }
  if (_riverDist && _fRiverShedFraction>0.0 && _nRiverShedCycles>0) {
    // decay of vibrio in rivers
    for (int i=0; i < _grid->getSize(); i++)
      if (_riverDist[i]>0) {
	_riverVibrio[i]*=Community::getVibrioDecay();     // level of Vibrio goes down
	_riverHyperVibrio[i] = 0.0;          // hyperinfectious vibrio disappears each day
      }
    // send vibrio down river (diffusion)
    for (int cycle=0; cycle<_nRiverShedCycles; cycle++) {
      for (int i=0; i < _grid->getSize(); i++)
	if (_riverDist[i]>0) {
	  _riverVibrio2[i] = 0.0;
	  _riverHyperVibrio2[i] = 0.0;
	}
      // copy vibrio to one cell downstream
      for (int i=0; i < _grid->getSize(); i++) 
	if (_riverDist[i]>0 &&
	    (_riverVibrio[i]>0.0 || _riverHyperVibrio[i]>0.0) &&
	    _riverNeighbors[i]>0) 
	  for (int x=-1; x<=1; x++)
	    for (int y=-1; y<=1; y++) {
	      int index = i+x+y*_nGridSizeX;
	      if ((x!=0 || y!=0) &&
		  _riverDist[index] > _riverDist[i]) {
		_riverVibrio2[index] += (1.0-_fRiverFlowDelta) * _riverVibrio[i]/_riverNeighbors[i];
		_riverHyperVibrio2[index] += (1.0-_fRiverFlowDelta) * _riverHyperVibrio[i]/_riverNeighbors[i];
	      }
	    }
      // add new vibrio from original sources
      for (int i=0; i < _grid->getSize(); i++)
	if (_riverDist[i]>0) {
	  double temp = _newVibrio[i]*_fRiverShedFraction/(_nRiverShedCycles+1);
	  _riverVibrio2[i] += temp;
	  _riverHyperVibrio2[i] += temp;
	}
      double *temp = _riverVibrio;
      _riverVibrio = _riverVibrio2;
      _riverVibrio2 = temp;
      temp = _riverHyperVibrio;
      _riverHyperVibrio = _riverHyperVibrio2;
      _riverHyperVibrio2 = temp;
    }
    // update river vibrio levels in communities
    for (int i=0; i < _grid->getSize(); i++) {
      if (_riverDist[i]>0) {
	int diff = _nCommunities[i];
	if (diff>_nMaxCommunitiesOnRiver)
	  diff=_nMaxCommunitiesOnRiver;
	for (int j=0; j<diff; j++) {
	  _community[_nCommunityStart[i]+j].setRiverVibrioLevel(_riverVibrio[i]);
	  _community[_nCommunityStart[i]+j].setRiverHyperVibrioLevel(_riverHyperVibrio[i]);
	}
      }
    }
  }

  // waning immunity from natural infection
  if (_nDaysToWane>0) {
    double probperday=1.0-exp(-1.0/(_nDaysToWane)); // daily waning probability
    int skip = ceil(log(gsl_rng_uniform(rng))/log(1-probperday)); // how many people in a row do not wane? draw from geometric distribution
    for (int personnum=0; personnum<Person::getLastPersonID(); personnum+=skip) {
      // make this person fully susceptible
      if (Person::personArray[personnum].getBaseSusceptibility()<1.0)
	Person::personArray[personnum].setBaseSusceptibility(1.0);
      skip = ceil(log(gsl_rng_uniform(rng))/log(1-probperday)); // choose next person to make susceptible
    }
  }

  // infect susceptibles
  for (int i=0; i < _nNumCommunities; i++)
    _community[i].drink(rng);
  for (int i=0; i < _nNumCommunities; i++)
    _community[i].tick(rng);

  // reactive vaccination if appropriate
  if (_fVaccinationTarget>0.0) {
    if (_nVaccinationThreshold>=0 && _nVaccinationDelay>=0) {
      // trigger (schedule) reactive vaccination
      for (int gridnum=0; gridnum < _grid->getSize(); gridnum++) {
	if (!_bGridVaccinated[gridnum] &&
	    _nVaccinationDay[gridnum]>_nDay+_nDayStartOffset+_nVaccinationDelay && 
	    getCumulativeSymptomatic(gridnum)>=_nVaccinationThreshold) {
	  _nVaccinationDay[gridnum] = _nDay+_nDayStartOffset+_nVaccinationDelay;
	}
      }
    }

    // vaccinate if ready
    if (_nNumVaccinesAvailable>_nNumVaccinesUsed) {
      //      cerr << "_nNumVaccinesAvailable>_nNumVaccinesUsed : " << _nNumVaccinesAvailable << ">" << _nNumVaccinesUsed << endl; ////////////
      int nNumPeopleWant = 0; // about how many people want vaccine
      int nNumCellsWant = 0;  // how many grid cells want vaccine
      for (int gridnum=0; gridnum < _grid->getSize(); gridnum++)
	if (!_bGridVaccinated[gridnum] &&
	    _nVaccinationDay[gridnum]<=_nDay+_nDayStartOffset &&
	    getNumResidents(gridnum)>0) {
	  nNumPeopleWant+=getNumResidents(gridnum);
	  nNumCellsWant++;
	}
      nNumPeopleWant *= _fVaccinationTarget;
      //      cerr << "Want vaccine: " << nNumPeopleWant << " people" <<  ", " << nNumCellsWant << " cells," << _fVaccinationTarget << endl;
      //      cerr << "Have: " << _nNumVaccinesAvailable-_nNumVaccinesUsed << " vaccines" << endl;
      if (nNumPeopleWant<=_nNumVaccinesAvailable-_nNumVaccinesUsed) {
	// enough vaccine for everyone who wants it
	for (int gridnum=0; gridnum < _grid->getSize(); gridnum++) {
	  if (!_bGridVaccinated[gridnum] &&
	      _nVaccinationDay[gridnum]<=_nDay+_nDayStartOffset &&
	      _nNumVaccinesAvailable>_nNumVaccinesUsed) {
	    for (int j=0; j<_nCommunities[gridnum]; j++) {
	      _nNumVaccinesUsed += _community[_nCommunityStart[gridnum]+j].vaccinate(rng, _fVaccinationTarget);
	      if (_fHygieneTarget>0.0)
		_community[_nCommunityStart[gridnum]+j].setHygiene(_fHygieneTarget);
	    }
	    _bGridVaccinated[gridnum] = true;
	  }
	}
      } else if (nNumPeopleWant>0) {
	// not enough vaccine - randomly choose grid cells to vaccinate until we run out
	int miss = 0; // number of times we chose grid cells that demand more vaccine than we have
	while (_nNumVaccinesAvailable>_nNumVaccinesUsed && nNumCellsWant>0 && miss<5) {
	  int gridcount = gsl_rng_uniform_int(rng, nNumCellsWant); // choose the nth unvaccinated cell
	  int gridid = 0;
	  for (gridid=0; gridid<_grid->getSize() && gridcount>=0; gridid++)
	    if (!_bGridVaccinated[gridid] &&
		_nVaccinationDay[gridid]<=_nDay+_nDayStartOffset &&
		getNumResidents(gridid)>0) {
	      gridcount--;
	    }
	  gridid--;

	  if (getNumResidents(gridid)*_fVaccinationTarget>_nNumVaccinesAvailable-_nNumVaccinesUsed) {
	    // not enough vaccine for this cell. try again
	    miss++;
	  } else {
	    // vaccinate this cell
	    miss=0;
	    for (int j=_nCommunityStart[gridid]; j<_nCommunityStart[gridid]+_nCommunities[gridid]; j++)
	      _nNumVaccinesUsed += _community[j].vaccinate(rng, _fVaccinationTarget);
	    _bGridVaccinated[gridid] = true;
	    nNumCellsWant--;
	  }
	}
      }
    }
  }
 
  // some people drive on the highways
  drive(rng);

  // "travel" by swapping infection and vaccination status of 2 random (healthy) people
  if (_fTravelProb>0.0) {
    int nNumTravelers = gsl_ran_binomial(rng, 
					 _fTravelProb,
					 Person::getLastPersonID()+1);
    for (int i=0; i<nNumTravelers; i++) {
      // should we check for symptomatics?????
      int person1id = gsl_rng_uniform_int(rng, Person::getLastPersonID()+1);
      int person2id = gsl_rng_uniform_int(rng, Person::getLastPersonID()+1);
      Person &p1 = Person::personArray[person1id];
      Person &p2 = Person::personArray[person2id];
      if ((!p1.isSymptomatic() || p1.getDaysInfectious()==0) &&
	  (!p2.isSymptomatic() || p2.getDaysInfectious()==0)) { // don't travel if symptomatic for more than one day
	Person temp;
	temp.copyInfectionStatus(p1);
	p1.copyInfectionStatus(p2);
	p2.copyInfectionStatus(temp);
      }
    }
  }
  return 1;
}
