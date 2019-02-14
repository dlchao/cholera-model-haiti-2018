/*
 * driver.cpp
 *
 * Each grid cell has a single cholera model.
 * Each pair of grid cells (i,j) for which residents of i work in j also has a cholera model.
 */

#include <math.h>
#include <assert.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "GridCells.h"
#include "Community.h"
#include "Population.h"

using namespace std;

int getRadius(GridCells *g, int center, int radius, int *buf, int bufsize) {
  const int maxneighbors = 6;
  int searched[1200];
  int searchedsize=1;
  int oldsearchedsize=1;
  searched[0] = center;
  // breadth first search
  for (int r=0; r<radius; r++) {
    oldsearchedsize=searchedsize;
    for (int i=searchedsize-1; i>=0; i--) {
      int neighbors[maxneighbors];
      int numneighbors = g->getNeighbors(searched[i], neighbors, maxneighbors);
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

int main(int argc, char *argv[]) {
  gsl_rng * rng = gsl_rng_alloc(gsl_rng_taus2);
  int randomseed = 5489;
  GridCells *g=NULL;
  Community *model=NULL;
  Population *pop=NULL;
  int communitysize = -1;
  bool bNoTransmit = false; // is cholera transmissible?
  int nStartCell = -1;   // where was the first case?
  int nGridSize = -1;    // use grid grid?
  double kappa = 70;     // vibrio environment half saturation
  double hyper = 100.00; // one-day hyperinfectiousness of shed vibrio (0.0)
  double beta = 0.12;    // 1.0; might be 0.25 according to paper
  double rho = 3.8;      // gravity model
  int runlength = 250;   // in days
  double fVaccinateFraction = 0.0;
  int nVaccinateThreshold = 0; // set to -1 for pre-vaccination
  int nVaccinateDelay = -1;    // set to >=0 for reactive vaccination
  int nVaccineStockpile = 1000000000;  // amount of vaccine available (after first day)
  int nVaccineCap = 1000000000;  // amount of vaccine ever available (maximum supply)
  int nVaccineFirstDay = -1;// when vaccines arrive
  int nNumVaccineLocationFirstDays=0;     // how many vaccination locations?
  int *nVaccineLocationFirstDay=NULL;     // when vaccination starts in one "location"
  string *szVaccineLocationLabel=NULL;    // names of locations for vaccination
  int nVaccinePerDay = 0;  // amount of vaccine that arrives per day after first day
  double fHygiene = 0.0;
  double fWorkingFraction=0.6; // fraction of working-age people (15y+) who work, from https://fred.stlouisfed.org/series/SLEMPTOTLSPZSHTI
  double fRiverShedFraction=0.1;
  double fRiverFlowDelta = 0.05;
  int nRiverShedCycles=9;
  bool bPrioritizeRiver = false;
  double fTravelProb = 0.000005;
  double fDriveProb =  0.00005;
  double fRainSheddingMultiplier1 = 1.0;
  double fRainSheddingMultiplier2 = 1.0;
  double fRainSheddingMultiplier3 = 1.0;
  double fRainSheddingMultiplier4 = 1.0;
  double fRainSheddingMultiplier5 = 1.0;
  double fVibrioHalflife = 14.0; // in days
  string szRainFile="haiti-rainfall.txt";
  string szTractFile=""; //grid-cells.csv";
  string szDailyOutputFile="location-infections.csv"; // daily output by location (labels)
  string szGridOutputFile=""; //="grid-infections.csv"; // daily output by grid square (huge raster)
  string szPeopleFile="";
  string szSummaryFile="";
  double fHouseholdContactProbability=0.01;
  double fSymptomaticFraction=0.2;
  double fAsymptomaticInfectiousness=0.1;
  int nTotalSymptomatic=0; // count of symptomatic cholera
  int nStartCalendarDay = 31+28+31+30+31+30+31+31+30+9; // first calendar day of simulation (October 9)
  int nWaningDays = -1;  // average number of days for immunity to wane
  
  if (argc>1) {
    string configfilename = argv[1];
    ifstream iss(configfilename.c_str());
    if (!iss) {
      cerr << "ERROR: " << configfilename << " not found." << endl;
      return -1;
    } else {
      cerr << "Reading parameters from " << configfilename << endl;
    }
    string line;
    while (getline(iss, line)) {
      istringstream linestream(line);
      string argname;
      if (linestream >> argname) {
      if (argname[0]=='#') {
	cerr << "Comment: " << line << endl;
      } else if (argname.compare("p")==0) {
	double p;
	linestream >> p;
	cerr << "p = " << p << endl;
	fHouseholdContactProbability = p;
      } else if (argname.compare("communitysize")==0) {
	linestream >> communitysize;
	cerr << "community size = " << communitysize << endl;
      } else if (argname.compare("k")==0 || argname.compare("kappa")==0) {
	linestream >> kappa;
	cerr << "kappa = " << kappa << endl;
      } else if (argname.compare("hyper")==0) {
	linestream >> hyper;
	cerr << "hyperinfectiousness = " << hyper << endl;
      } else if (argname.compare("beta")==0) {
	linestream >> beta;
	cerr << "beta = " << beta << endl;
      } else if (argname.compare("rho")==0) {
	linestream >> rho;
	cerr << "rho = " << rho << endl;
      } else if (argname.compare("m")==0 || argname.compare("asymptomaticmultiplier")==0) {
	linestream >> fAsymptomaticInfectiousness;
	cerr << "asymptomatic infectiousness multiplier = " << fAsymptomaticInfectiousness << endl;
      } else if (argname.compare("symptomaticfraction")==0) {
	linestream >> fSymptomaticFraction;
	cerr << "symptomatic fraction = " << fSymptomaticFraction << endl;
      } else if (argname.compare("randomseed")==0) {
	linestream >> randomseed;
	cerr << "random seed = " << randomseed << endl;
      } else if (argname.compare("vaccinatefraction")==0) {
	linestream >> fVaccinateFraction;
	cerr << "vaccinate fraction = " << fVaccinateFraction << endl;
      } else if (argname.compare("vaccinatethreshold")==0) {
	linestream >> nVaccinateThreshold;
	cerr << "vaccination threshold = " << nVaccinateThreshold << endl;
      } else if (argname.compare("vaccinatedelay")==0) {
	linestream >> nVaccinateDelay;
	cerr << "vaccination delay = " << nVaccinateDelay << endl;
      } else if (argname.compare("vaccinefirstday")==0) {
	linestream >> nVaccineFirstDay;
	cerr << "first (calendar) day of vaccination = " << nVaccineFirstDay << endl;
      } else if (argname.compare("vaccinelocationfirstday")==0) {
	if (nNumVaccineLocationFirstDays<=0) {
	  nVaccineLocationFirstDay = new int[100];
	  szVaccineLocationLabel = new string[100];
	  nNumVaccineLocationFirstDays = 0;
	}
	linestream >> szVaccineLocationLabel[nNumVaccineLocationFirstDays]; // vaccinate cells with this label
	linestream >> nVaccineLocationFirstDay[nNumVaccineLocationFirstDays]; // on this calendar day
	nNumVaccineLocationFirstDays++;
      } else if (argname.compare("vaccineperday")==0) {
	linestream >> nVaccinePerDay;
	cerr << "new vaccine available per day = " << nVaccinePerDay << endl;
      } else if (argname.compare("wanedays")==0) {
	linestream >> nWaningDays;
	cerr << "average days for waning immunity = " << nWaningDays << endl;
      } else if (argname.compare("workingfraction")==0) {
	linestream >> fWorkingFraction;
	cerr << "working fraction = " << fWorkingFraction << endl;
      } else if (argname.compare("rivershedfraction")==0) {
	linestream >> fRiverShedFraction;
	cerr << "fraction shed into the river = " << fRiverShedFraction << endl;
      } else if (argname.compare("rivershedcycles")==0) {
	linestream >> nRiverShedCycles;
	cerr << "river shed cycles per day = " << nRiverShedCycles << endl;
      } else if (argname.compare("riverflowdelta")==0) {
	linestream >> fRiverFlowDelta;
	cerr << "fraction that disappears per distance travelled in river = " << fRiverFlowDelta << endl;
      } else if (argname.compare("rainshedmultipliers")==0) {
	linestream >> fRainSheddingMultiplier1 >> fRainSheddingMultiplier2 >> fRainSheddingMultiplier3 >> fRainSheddingMultiplier4 >> fRainSheddingMultiplier5;
	cerr << "no/low/medium/high/extreme rainfall shedding multipliers = " << fRainSheddingMultiplier1 << " / " << fRainSheddingMultiplier2 << " / " << fRainSheddingMultiplier3 << " / " << fRainSheddingMultiplier4 << " / " << fRainSheddingMultiplier5 << endl;
      } else if (argname.compare("vibriohalflife")==0) {
	linestream >> fVibrioHalflife;
	cerr << "environmental vibrio half life (days) = " << fVibrioHalflife << endl;
      } else if (argname.compare("prioritizeriver")==0) {
	bPrioritizeRiver = true;
	cerr << "Prioritizing rivers" << endl;
      } else if (argname.compare("travelprob")==0) {
	linestream >> fTravelProb;
	cerr << "daily travel probability = " << fTravelProb << endl;
      } else if (argname.compare("driveprob")==0) {
	linestream >> fDriveProb;
	cerr << "daily highway driving probability = " << fDriveProb << endl;
      } else if (argname.compare("ves")==0 || argname.compare("VES")==0) {
	double temp;
	linestream >> temp;
	cerr << "VE_S = " << temp << endl;
	Community::setVES(temp);
      } else if (argname.compare("vep")==0 || argname.compare("VEP")==0) {
	double temp;
	linestream >> temp;
	cerr << "VE_P = " << temp << endl;
	Community::setVEP(temp);
      } else if (argname.compare("vei")==0 || argname.compare("VEI")==0) {
	double temp;
	linestream >> temp;
	cerr << "VE_I = " << temp << endl;
	Community::setVEI(temp);
      } else if (argname.compare("veu5reduction")==0) {
	double temp;
	linestream >> temp;
	cerr << "VE u5 reduction = " << temp << endl;
	Community::setVE_u5reduction(temp);
      } else if (argname.compare("vewaning")==0) {
	double temp;
	linestream >> temp;
	cerr << "VE linear waning per day = " << temp << endl;
	Community::setVaccineEfficacyWaning(temp);
      } else if (argname.compare("mininfectiousperiod")==0) {
	int temp;
	linestream >> temp;
	cerr << "minimum infectious period = " << temp << endl;
	Community::setMinInfectiousDays(temp);
      } else if (argname.compare("maxinfectiousperiod")==0) {
	int temp;
	linestream >> temp;
	cerr << "maximum infectious period = " << temp << endl;
	Community::setMaxInfectiousDays(temp);
      } else if (argname.compare("days")==0 || argname.compare("runlength")==0) {
	linestream >> runlength;
	cerr << "runlength (simulation days) = " << runlength << endl;
      } else if (argname.compare("vaccinestockpile")==0) {
	linestream >> nVaccineStockpile;
	cerr << "vaccine stockpile size = " << nVaccineStockpile << endl;
      } else if (argname.compare("vaccinecap")==0) {
	linestream >> nVaccineCap;
	cerr << "total vaccines = " << nVaccineCap << endl;
      } else if (argname.compare("hygiene")==0) {
	linestream >> fHygiene;
	cerr << "hygiene = " << fHygiene << endl;
      } else if (argname.compare("cellfile")==0) {
	linestream >> szTractFile;
	cerr << "cell file = " << szTractFile << endl;
      } else if (argname.compare("rainfile")==0) {
	linestream >> szRainFile;
	cerr << "rainfall file = " << szRainFile << endl;
      } else if (argname.compare("outputfile")==0) {
	if (linestream >> szDailyOutputFile) {
	  cerr << "daily location file = " << szDailyOutputFile << endl;
	} else {
	  cerr << "no daily location output file" << endl;
	  szDailyOutputFile="";
	}
      } else if (argname.compare("summaryfile")==0) {
	if (linestream >> szSummaryFile) {
	  cerr << "summary file = " << szSummaryFile << endl;
	} else {
	  cerr << "no summary file" << endl;
	  szSummaryFile="";
	}
      } else if (argname.compare("peoplefile")==0) {
	if (linestream >> szPeopleFile) {
	  cerr << "people file = " << szPeopleFile << endl;
	} else {
	  cerr << "no people file" << endl;
	  szPeopleFile="";
	}
      } else if (argname.compare("dailycellfile")==0) {
	if (linestream >> szGridOutputFile) {
	  cerr << "daily cell file = " << szGridOutputFile << endl;
	} else {
	  cerr << "no daily cell file" << endl;
	  szGridOutputFile="";
	}
      } else if (argname.compare("grid")==0) {
	linestream >> nGridSize;
	cerr << "grid dimension = " << nGridSize << endl;
      } else if (argname.compare("notransmit")==0) {
	bNoTransmit=true;
	cerr << "no secondary transmission" << endl;
      } else {
	cerr << "****Unknown parameter: " << argname << endl;
	exit(-1);
      }
    }
    }
  }  
  gsl_rng_set(rng, randomseed);
  Community::setVibrioDecay(1.0-1.0/fVibrioHalflife);
  Community::setVibrio50(kappa);
  Community::setHyperVibrioMultiplier(hyper);
  Community::setBeta(beta);
  if (fHouseholdContactProbability>=0.0)
    Community::setHouseholdContactProbability(fHouseholdContactProbability);
  if (fSymptomaticFraction>=0.0)
    Community::setSymptomaticFraction(fSymptomaticFraction);
  if (fAsymptomaticInfectiousness>=0.0)
    Community::setAsymptomaticInfectiousnessMultiplier(fAsymptomaticInfectiousness);

  if (nNumVaccineLocationFirstDays>0) {
    for (int i=0; i<nNumVaccineLocationFirstDays; i++) {
      cerr << "Vaccinate location " << szVaccineLocationLabel[i] << " starting on calendar day " << nVaccineLocationFirstDay[i] << endl;
    }
  }

    pop=new Population();
    pop->loadPopulation(rng, 
		      fWorkingFraction,
		      0.3, 0.64, rho,
		      "haitipopulation-dim.txt",
		      "haitipopulation.txt",
		      bNoTransmit);
    pop->loadRainfall(szRainFile.c_str());
    pop->loadRivers("rivers.txt");
    pop->loadHighways("highways.txt");
    pop->setRiverShedFraction(fRiverShedFraction);
    pop->setRiverShedCycles(nRiverShedCycles);
    pop->setRiverFlowDelta(fRiverFlowDelta);
    pop->setTravelProb(fTravelProb);
    pop->setDriveProb(fDriveProb);
    if (communitysize>0) {
      pop->setTargetCommunitySize(communitysize);
    }
    pop->setRainSheddingMultipliers(fRainSheddingMultiplier1,fRainSheddingMultiplier2,fRainSheddingMultiplier3,fRainSheddingMultiplier4,fRainSheddingMultiplier5);
    pop->setVaccinationTarget(fVaccinateFraction);
    pop->setDayToStart(nStartCalendarDay); // October 9
    pop->setWaningDays(nWaningDays);
    if (runlength>365) {
      pop->setDayToAge(31+28+15); // March 15
    }
    if (nVaccineFirstDay<0)
      pop->setNumVaccinesAvailable(nVaccineStockpile);
    else
      pop->setNumVaccinesAvailable(0);
    pop->setHygieneTarget(fHygiene);
    if (bPrioritizeRiver) {
      // send vaccine to river when available
      for (int i=0; i < pop->getNumCells(); i++)
	if (pop->getRiver(i)>0)
	  pop->prioritizeCell(i);
      pop->setVaccinationThreshold(100000000);
      pop->setVaccinationDelay(-1);
    } else {
      pop->setVaccinationThreshold(nVaccinateThreshold);
      pop->setVaccinationDelay(nVaccinateDelay);
    }

    // start epidemic in Mirebalais 18°50′0″N 72°6′19″W
    //    int source = pop->infect(rng, 100, 57, 98);
    if (!bNoTransmit) {
      // Petite Riviere de l'Artibonite: 19.121733, -72.477046
      //      pop->infect(rng, 50, double longitude, double latitude, bool bMakeSymptomatic      
      // seed in PETITE_RIVIERE_DE_L'ARTI, SAINT-MARC, VERRETTES, LASCAHOBAS, MIREBALAIS
      // pops of 144467, 215390, 112986, 28507, 83375
      // The first days of cas vu in Artibonite: 1111 1840 1833 2138 1675
      // First days of cas vu in Centre: 61  92 244 165
      cerr << "seeding epidemic" << endl;
      for (int i=0; i < pop->getNumCells(); i++) {
	if (pop->getRiver(i)>0) {
	  double newlongitude = pop->getlongitude(i);
	  double newlatitude = pop->getlatitude(i);
	  if ((newlongitude+72.477046)*(newlongitude+72.477046)+(newlatitude-19.121733)*(newlatitude-19.121733) < 0.0004) { // Petite Riviere de l'Artibonite
	    //	    cerr << "Petite Riviere vibrio at cell " << i << " : " << pop->getX(i) << "," << pop->getY(i) << " ; " << newlongitude << "," << newlatitude << endl;
	    pop->setRiverVibrioLevel(i, 100);
	  } else if ((newlongitude+72.695976)*(newlongitude+72.695976)+(newlatitude-19.108722)*(newlatitude-19.108722) < 0.0004) { // Saint-Marc
	    //	    cerr << "Saint-Marc vibrio at cell " << i << " : " << pop->getX(i) << "," << pop->getY(i) << " ; " << newlongitude << "," << newlatitude << endl;
	    pop->setRiverVibrioLevel(i, 100);
	  } else if ((newlongitude+72.465456)*(newlongitude+72.465456)+(newlatitude-19.050562)*(newlatitude-19.050562) < 0.0004) { // Verrettes
	    //	    cerr << "Verrettes vibrio at cell " << i << " : " << pop->getX(i) << "," << pop->getY(i) << " ; " << newlongitude << "," << newlatitude << endl;
	    pop->setRiverVibrioLevel(i, 100);
	  } else if ((newlongitude+72.104514)*(newlongitude+72.104514)+(newlatitude-18.835532)*(newlatitude-18.835532) < 0.0004) { // Mirbalais
	    //	    cerr << "Mirbalais vibrio at cell " << i << " : " << pop->getX(i) << "," << pop->getY(i) << " ; " << newlongitude << "," << newlatitude << endl;
	    //	    pop->setRiverVibrioLevel(i, 1);
	    pop->setRiverVibrioLevel(i, 100);
	  }
	}
      }
    } else {
      // R0 test run
      nStartCell = pop->infectOne(rng);
      cerr << "infected one person" << endl;
    }

  // pre-vaccinate?
  if (fVaccinateFraction>0.0 && nVaccinateThreshold<0)
    pop->prevaccinate(rng);

  // output tract file
  if (pop && szTractFile.length()>0) {
    ofstream tractFile;
    tractFile.open(szTractFile.c_str());
    if(tractFile.fail()) {
      cerr << "ERROR: Cell file '" << szTractFile << "' cannot be open for writing." << endl;
      return false;
    }
    cerr << "outputing cell information to " << szTractFile << endl;
    tractFile << "id,x,y,long,lat,pop,river,highway,label" << endl;
    for (int i=0; i < pop->getNumCells(); i++) {
      double popsize=pop->getNumResidents(i);
      int river = pop->getRiver(i);
      int highway = (pop->isHighway(i)?1:0);
      if (popsize>0 || river>0) {
	double x = pop->getX(i);
	double y = pop->getY(i);
	tractFile << i << "," << x << "," << y << "," << (pop->getlongitude(i)) <<  "," << (pop->getlatitude(i)) <<  "," << popsize << "," << river << "," << highway << "," << pop->getLabel(i) << endl;
      }
    }
    tractFile.close();
  }

  // output people file
  if (szPeopleFile.length()>0) {
    cerr << "outputing people information to " << szPeopleFile << endl;
    ofstream peopleFile;
    peopleFile.open(szPeopleFile.c_str());
    if(peopleFile.fail()) {
      cerr << "ERROR: People file '" << szPeopleFile << "' cannot be open for writing." << endl;
      return false;
    }
    peopleFile << "id,age,home,work" << endl;
    for (int i=0; i<Person::getLastPersonID(); i++) { // IDs start at 0
      peopleFile << Person::personArray[i].getID() << "," << Person::personArray[i].getAge() << "," << Person::personArray[i].getHomeCommunity() << "," << Person::personArray[i].getWorkCommunity() << endl;
    }
    peopleFile.close();
  }

  // set up main output file (location-level)
  ofstream outputFile;
  if (szDailyOutputFile.length()>0) {
    outputFile.open(szDailyOutputFile.c_str());
    if(outputFile.fail()) {
      cerr << "ERROR: Output file '" << szDailyOutputFile << "' cannot be open for writing." << endl;
      return false;
    }
    cerr << "outputing daily information to " << szDailyOutputFile << endl;
    outputFile << "time,location,residents,susceptible,symptomatic,new symptomatic,new symptomatic U5,infectious,vaccinated,vaccinesused"<< endl;;
  } else {
    cerr << "not outputing information" << endl;
  }
  // set up grid output file
  ofstream gridOutputFile;
  if (szGridOutputFile.length()>0) {
    gridOutputFile.open(szGridOutputFile.c_str());
    if(gridOutputFile.fail()) {
      cerr << "ERROR: Grid output file '" << szGridOutputFile << "' cannot be open for writing." << endl;
      return false;
    }
    cerr << "outputing daily grid information to " << szGridOutputFile << endl;
    gridOutputFile << "time,id,susceptible,infectious,symptomatic,newsymptomatic,cumulativesymptomatic,vibrio,hypervibrio,rivervibrio,riverhypervibrio,vaccines" << endl;
  } else {
    cerr << "not outputing grid information" << endl;
  }
  
  // main loop
  for (int t=0; t<runlength; t++) {
    if (t%100==0)
      cerr << "simulation day " << t << endl;
    if (nNumVaccineLocationFirstDays>0) {
      for (int i=0; i<nNumVaccineLocationFirstDays; i++) {
	if (nVaccineLocationFirstDay[i]==nStartCalendarDay+t) {
	  cerr << "Start vaccinating " << szVaccineLocationLabel[i] << " on calendar day " << nVaccineLocationFirstDay[i] << endl;
	  if ((nVaccineFirstDay<0) || (nVaccineFirstDay>nStartCalendarDay+t))
	    nVaccineFirstDay = nStartCalendarDay+t; // trigger "global" vaccination flag
	  for (int j=0; j < pop->getNumCells(); j++)
	    if (szVaccineLocationLabel[i].compare(pop->getLabel(j))==0) {
	      pop->prioritizeCell(j);
	    }
	}
      }
    }
    //    cerr << t << " vaccines available: " << pop->getNumVaccinesAvailable() << ", starting on day " << nVaccineFirstDay << endl;
    if (nVaccineFirstDay==nStartCalendarDay+t) {
      pop->setNumVaccinesAvailable(nVaccineStockpile);
      pop->setVaccinationThreshold(nVaccinateThreshold);
      pop->setVaccinationDelay(nVaccinateDelay);
      cerr << "Vaccination starts. Vaccines available: " << nVaccineStockpile << endl;
    }
    if (nVaccinePerDay>0 && nVaccineFirstDay<=nStartCalendarDay+t)
      pop->setNumVaccinesAvailable(pop->getNumVaccinesAvailable()+nVaccinePerDay);
    if (pop->getNumVaccinesAvailable()>nVaccineCap)
      pop->setNumVaccinesAvailable(nVaccineCap);

    if (pop) {
      pop->step(rng);
      if (bPrioritizeRiver) {
	// check to see if river is all vaccinated
	int count = 0;
	for (int i=0; i < pop->getNumCells(); i++)
	  if (pop->wantVaccine(i) && !pop->isVaccinated(i))
	    count++;
	if (count==0) {
	  // done vaccinating river
	  bPrioritizeRiver = false;
	  pop->setVaccinationThreshold(nVaccinateThreshold);
	  pop->setVaccinationDelay(nVaccinateDelay);
	}
      }
    } else {
      cerr << " no pop!" << endl;
      for (int i=0; i < g->getSize(); i++)
	model[i].poop(rng);
      for (int i=0; i < g->getSize(); i++)
	model[i].drink(rng);
      for (int i=0; i < g->getSize(); i++)
	model[i].tick(rng);
    }
    // print location-level results
    if (pop && outputFile.is_open()) {
      const int MAXLABELS=50;
      int residents[MAXLABELS]; // can only handle 50 unique "labels" for cells
      int susceptibles[MAXLABELS];
      int newsymptomatic[MAXLABELS];
      int newsymptomaticu5[MAXLABELS]; // new symptomatics under 5 years old
      int symptomatic[MAXLABELS];
      int infectious[MAXLABELS];
      int vaccinated[MAXLABELS];
      int vaccused[MAXLABELS];
      memset(residents, 0, sizeof(int)*MAXLABELS);
      memset(susceptibles, 0, sizeof(int)*MAXLABELS);
      memset(newsymptomatic, 0, sizeof(int)*MAXLABELS);
      memset(newsymptomaticu5, 0, sizeof(int)*MAXLABELS);
      memset(symptomatic, 0, sizeof(int)*MAXLABELS);
      memset(infectious, 0, sizeof(int)*MAXLABELS);
      memset(vaccinated, 0, sizeof(int)*MAXLABELS);
      memset(vaccused, 0, sizeof(int)*MAXLABELS);
      int labelnum = -1;
      string lastdept = "----";
      for (int cellnum=0; cellnum < pop->getNumCells(); cellnum++) {
	if (lastdept.compare(pop->getLabel(cellnum))!=0) { // this grid cell probably has same label as previous one
	  lastdept = pop->getLabel(cellnum);
	  for (labelnum=0; labelnum<pop->getNumUniqueLabels(); labelnum++) {
	    if (lastdept.compare(pop->getUniqueLabel(labelnum))==0) { // found it
	      break;
	    }
	  }
	}
	if (labelnum<MAXLABELS) {
	  residents[labelnum] += pop->getNumResidents(cellnum);
	  newsymptomatic[labelnum]+=pop->getNumNewSymptomatic(cellnum);
	  newsymptomaticu5[labelnum]+=pop->getNumNewSymptomatic(cellnum,0,4);
	  symptomatic[labelnum]+=pop->getNumSymptomatic(cellnum);
	  susceptibles[labelnum]+=pop->getNumSusceptible(cellnum);
	  infectious[labelnum]+=pop->getNumInfectious(cellnum);
	  vaccinated[labelnum]+=pop->getNumVaccinated(cellnum);
	  vaccused[labelnum]+=pop->getNumVaccinesUsed(cellnum);
	}
      }
      for (int labelnum=0; labelnum<pop->getNumUniqueLabels(); labelnum++) {
	outputFile << pop->getDay() << "," << pop->getUniqueLabel(labelnum) << "," << residents[labelnum] << "," << susceptibles[labelnum] << "," << symptomatic[labelnum] << "," << newsymptomatic[labelnum] <<  "," << newsymptomaticu5[labelnum] <<  "," << infectious[labelnum] << "," << vaccinated[labelnum] <<  "," << vaccused[labelnum] << endl;
	if (!gridOutputFile.is_open())
	  nTotalSymptomatic += newsymptomatic[labelnum];
      }
    }
    
    // print day's grid square results
    if (gridOutputFile.is_open()) {
      for (int i=0; i < (pop?pop->getNumCells():g->getSize()); i++) {
	int infectious, s,news,sus,cases,res,vac;
	double vlevel, hvlevel, rvlevel, rhvlevel;
	if (pop) {
	  res=pop->getNumResidents(i);
	  infectious=pop->getNumInfectious(i);
	  s=pop->getNumSymptomatic(i);
	  news=pop->getNumNewSymptomatic(i);
	  sus=pop->getNumSusceptible(i);
	  cases=pop->getCumulativeSymptomatic(i);
	  vac=pop->getNumVaccinated(i);
	  vlevel = pop->getVibrioLevel(i);
	  hvlevel = pop->getHyperVibrioLevel(i);
	  rvlevel = pop->getRiverVibrioLevel(i);
	  rhvlevel = pop->getRiverHyperVibrioLevel(i);
	  nTotalSymptomatic += news;
	} else {
	  res=model[i].getNumResidents();
	  infectious=model[i].getNumInfectious();
	  s=model[i].getNumSymptomatic();
	  news=model[i].getNumNewSymptomatic();
	  sus=model[i].getNumSusceptible();
	  cases=model[i].getCumulativeSymptomatic();
	  vac=pop->getNumVaccinated(i);
	  vlevel = model[i].getVibrioLevel();
	  hvlevel = model[i].getHyperVibrioLevel();
	  rvlevel = 0.0;
	  rhvlevel = 0.0;
	}
	if (res>0)
	  gridOutputFile << pop->getDay() << "," << i << "," << sus << "," << infectious << "," << s << "," << news << "," <<  cases << "," << vlevel << "," << hvlevel << "," << rvlevel << "," << rhvlevel << "," << vac << endl;
      }
    }
  }
  if (outputFile.is_open())
    outputFile.close();
  if (gridOutputFile.is_open())
    gridOutputFile.close();

  cerr << "Total symptomatic = " << nTotalSymptomatic << endl;
  
  // final output
  // family attack rates
  int family = -1;
  int numcase = 0;
  int numinfected = 0;
  int numpeople = 0;
  int numcasetotal = 0;
  int numinfectedtotal = 0;
  int numpeopletotal = 0;
  int numcasetotalminus = 0;
  int numinfectedtotalminus = 0;
  int numpeopletotalminus = 0;
  int count = 0;
  int countcase = 0;
  for (int personid=0; personid<Person::getLastPersonID(); personid++) {
    Person &person = Person::personArray[personid];
    if (person.getFamily()!=family) {
      if (numcase>0) {
	numcasetotal+=numcase;
	numinfectedtotal+=numinfected;
	numpeopletotal+=numpeople;
	numcasetotalminus+=numcase-1;
	numinfectedtotalminus+=numinfected-1;
	numpeopletotalminus+=numpeople-1;
      }
      family = person.getFamily();
      numcase = numinfected = numpeople = 0;
    }
    int temp=person.getInfectionCount();
    numinfected+=temp;
    count+=temp;
    temp=person.getSymptomaticCount();
    numcase+=temp;
    countcase+=temp;
    numpeople++;
  }
  cerr << "total cases/infected/pop: " << countcase << "/" << count << "/" << Person::getLastPersonID() << endl;
  //  cerr << "total family infected with case: " << numcasetotal << "," << numinfectedtotal << "/" << numpeopletotal << endl;
  //  cerr << "total family infected with case (minus one infected person): " << numinfectedtotalminus << "/" << numpeopletotalminus << endl;

  if (szSummaryFile.length()>0) {
    ofstream summaryFile;
    summaryFile.open(szSummaryFile.c_str());
    if(summaryFile.fail()) {
      cerr << "ERROR: Summary file '" << szSummaryFile << "' cannot be open for writing." << endl;
    } else {
      int nSumCases = 0;
      int nSumRecovered = 0;
      summaryFile << "Started in cell, " << nStartCell << endl;
      summaryFile << "Immune at origin, " << pop->getNumImmune(nStartCell) << endl;
      for (int i=0; i < pop->getNumCells(); i++) {
	nSumCases += pop->getCumulativeSymptomatic(i);
	nSumRecovered += pop->getNumImmune(i);
      }
      summaryFile << "Cases total (est), " << nSumCases << endl;
      summaryFile << "Cases total (real), " << nTotalSymptomatic << endl;
      summaryFile << "Immune total, " << nSumRecovered << endl;
      summaryFile.close();
    }
  }
  
  // clean up and exit
  delete [] model;
  delete g;
  delete pop;
  gsl_rng_free(rng);
  return 0;
}
