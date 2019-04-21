/*
 * Community.cpp
 * 11/2010 (updated 2/2019)
 * Contains classes for Person and Community (well-mixed collection of people)
 */

#include <iostream>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "Community.h"

using namespace std;

const int Person::MAXPERSONS = 11000000; // total people in simulation
Person Person::personArray[MAXPERSONS];
int Person::_nNextID = 0;
int Person::_nUsedID = 0;
const int Community::MAXFAMILYSIZE = 30;
double Community::_fFamiltySizeCDF[Community::MAXFAMILYSIZE] = 
  {1.0,1.0,1.0,1.0,1.0,
   1.0,1.0,1.0,1.0,1.0,
   1.0,1.0,1.0,1.0,1.0,
   1.0,1.0,1.0,1.0,1.0,
   1.0,1.0,1.0,1.0,1.0,
   1.0,1.0,1.0,1.0,1.0}; // cumulative fraction of families of size n+1
int Community::familyArray[Person::MAXPERSONS];
double Community::_fWorkTimeFraction = 0.3;
double Community::_fHouseholdContactProbability = 0.0;
double Community::_fVibrioDecay = 1.0-1.0/30.0;
double Community::_fVibrio50 = 650000;
double Community::_fHyperVibrioMultiplier=700.0;
double Community::_fVibrioBeta=1.0;
int Community::_nNextID = 0;
double Community::_fSymptomaticFraction = 0.1;
double Community::_fAsymptomaticInfectiousnessMultiplier = 0.1;
double Community::_fRiverSheddingMultiplier=1.0;

int Community::_nMinInfectiousDays = 7;
int Community::_nMaxInfectiousDays = 13;
double fWithdrawFraction = 0.75;
double Community::_fVES = 0.0;
double Community::_fVEI = 0.5;
double Community::_fVEP = 0.64;
double Community::_fVE_u5_reduction = 0.0;
double Community::_fVE_linearwaning = 0.0; // default is 0 waning of VE
const int _nMaxIncubationDays = 5;
double _fIncubationCDF[_nMaxIncubationDays] = {0.4,0.4,0.07,0.07,0.06};
int _fIncubationCDF100[100] = 
  {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
   1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
   2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
   2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
   3,3,3,3,3,3,3,4,4,4,4,4,4,4,5,5,5,5,5,5};

void Person::vaccinate() { 
  if (_nVaccinationDay<0) {
    _nVaccinationDay=0;
    _nAgeAtVaccination=_nAge;
  }
}

void Person::prevaccinate() { 
  if (_nVaccinationDay<0) {
    _nVaccinationDay=0;
    _nAgeAtVaccination=_nAge;
    if (Community::_fVE_u5_reduction==0.0 || _nAgeAtVaccination>=5) {
      _fInfectiousness = _fBaseInfectiousness * (1.0-Community::_fVEI);
      _fSusceptibility = _fBaseSusceptibility * (1.0-Community::_fVES);
    } else {
      _fInfectiousness = _fBaseInfectiousness * (1.0-Community::_fVEI*(1.0-Community::_fVE_u5_reduction));
      _fSusceptibility = _fBaseSusceptibility * (1.0-Community::_fVES*(1.0-Community::_fVE_u5_reduction));
    }
  }
}

void Person::infect(gsl_rng *rng) {
  if (_nTemporaryImmunityCountdown<0 && getSusceptibility()>0.0) {
    int firstday = _fIncubationCDF100[gsl_rng_uniform_int(rng, 100)];
    _nIncubationCountdown=firstday;
    _fBaseSusceptibility = _fSusceptibility = 0.0;
    _nInfectionCount++;
  }
}

// makesymptomatic - turns this (possibly uninfected) person into a symptomatic case
int Person::makesymptomatic(gsl_rng *rng) {
  if (_bSymptomatic)
    return 0;
  _nInfectiousCountdown = Community::_nMinInfectiousDays+gsl_rng_uniform_int(rng, 1+Community::_nMaxInfectiousDays-Community::_nMinInfectiousDays); // note that infectiousness can be up to Community::_nMaxInfectiousDays
  _nInfectiousCountup=0;
  _fBaseSusceptibility=_fSusceptibility=0.0;
  _bSymptomatic = true;
  _nSymptomaticCount++;
  _fBaseInfectiousness=_fInfectiousness=1.0;
  return 1;
}

// makeImmune - makes this person immune to infection
void Person::makeImmune() {
  _fBaseSusceptibility=_fSusceptibility=0.0;
}

// copyInfectionStatus - copies infection and vaccination data from p
void Person::copyInfectionStatus(const Person &p) {
  _nVaccinationDay = p.getVaccinationDay();
  _nAgeAtVaccination = p.getAgeAtVaccination();
  _nIncubationCountdown = p.getIncubationCountdown();
  _nInfectiousCountdown = p.getInfectiousCountdown();
  _nInfectiousCountup = p.getInfectiousCountup();
  _nTemporaryImmunityCountdown = p.getTemporaryImmunityCountdown();
  _fSusceptibility = p.getSusceptibility();
  _fInfectiousness = p.getInfectiousness();
  _fBaseSusceptibility = p.getBaseSusceptibility();
  _fBaseInfectiousness = p.getBaseInfectiousness();
  _bSymptomatic = p.isSymptomatic();
  _nSymptomaticCount = p.getSymptomaticCount();
  _nInfectionCount = p.getInfectionCount();
}

// copyImmuneStatus - copies immunity and vaccination data from p (but not active infections)
void Person::copyImmuneStatus(const Person &p) {
  _nVaccinationDay = p.getVaccinationDay();
  _nAgeAtVaccination = p.getAgeAtVaccination();
  //  _nIncubationCountdown = p.getIncubationCountdown();
  //  _nInfectiousCountdown = p.getInfectiousCountdown();
  //  _nInfectiousCountup = p.getInfectiousCountup();
  _nTemporaryImmunityCountdown = p.getTemporaryImmunityCountdown();
  _fSusceptibility = p.getSusceptibility();
  //  _fInfectiousness = p.getInfectiousness();
  _fBaseSusceptibility = p.getBaseSusceptibility();
  //  _fBaseInfectiousness = p.getBaseInfectiousness();
  //  _bSymptomatic = p.isSymptomatic();
  _nSymptomaticCount = p.getSymptomaticCount();
  _nInfectionCount = p.getInfectionCount();
}

// resetInfectionStatus - makes person naive and unvaccinated
void Person::resetInfectionStatus() {
  _nVaccinationDay=-1;
  _nAgeAtVaccination=-1;
  _nIncubationCountdown=-1;
  _nInfectiousCountdown=-1;
  _nInfectiousCountup=-1;
  _nTemporaryImmunityCountdown=-1;
  _fInfectiousness=0.0;
  _fSusceptibility=1.0;
  _fBaseInfectiousness=0.0;
  _fBaseSusceptibility=1.0;
  _bSymptomatic=false;
  //  _nSymptomaticCount=0;
  //    _fHygiene = 0.0;
}

void Person::setBaseSusceptibility(double d) {
  _fBaseSusceptibility = d;
  if (_nVaccinationDay>=0) {
    if (_nInfectiousCountdown>0) {
      //      _fInfectiousness = _fBaseInfectiousness * (1.0-Community::_fVEI*Community::_fVaccineBuildup[_nVaccinationDay]);
    } else if (_fBaseSusceptibility>0.0) {
      if (Community::_fVE_u5_reduction==0.0 || _nAgeAtVaccination>=5) 
	_fSusceptibility = _fBaseSusceptibility * (1.0-Community::_fVES*Community::getVaccineEfficacyWaning(_nVaccinationDay));
      else
	_fSusceptibility = _fBaseSusceptibility * (1.0-Community::_fVES*(1.0-Community::_fVE_u5_reduction)*Community::getVaccineEfficacyWaning(_nVaccinationDay));
    }
  } else
    _fSusceptibility = _fBaseSusceptibility;
}

// step - update timers and internal state for this person
bool Person::step(gsl_rng *rng) {
  if (_nTemporaryImmunityCountdown>=0)
    _nTemporaryImmunityCountdown--;
  if (_nInfectiousCountdown>0) {
    _nInfectiousCountdown--;
    _nInfectiousCountup++;
  } else if (_nInfectiousCountdown==0) {
    _fBaseInfectiousness=_fInfectiousness=0.0;
    _bSymptomatic = false;
  }
  if (_nIncubationCountdown>0)
    _nIncubationCountdown--;
  else if (_nIncubationCountdown==0) {
    // done incubating; now infectious
    _nIncubationCountdown=-1;
    if ((_nVaccinationDay<0 && 
	 gsl_rng_uniform(rng)<Community::getSymptomaticFraction()) ||
	(_nVaccinationDay>=0 && 
	 gsl_rng_uniform(rng)<Community::getSymptomaticFraction()*(1.0-Community::_fVEP*Community::getVaccineEfficacyWaning(_nVaccinationDay)*((Community::_fVE_u5_reduction==0.0 || _nAgeAtVaccination>=5)?1.0:(1.0-Community::_fVE_u5_reduction))))) {
      _bSymptomatic = true;
      _nSymptomaticCount++;
      _fBaseInfectiousness=_fInfectiousness=1.0;
    } else {
      _bSymptomatic = false;
      _fBaseInfectiousness=_fInfectiousness=Community::getAsymptomaticInfectiousnessMultiplier();
    }
    _nInfectiousCountdown = Community::_nMinInfectiousDays+gsl_rng_uniform_int(rng, 1+Community::_nMaxInfectiousDays-Community::_nMinInfectiousDays); // note that infectiousness can be up to Community::_nMaxInfectiousDays
    _nInfectiousCountup=0;
  }
  if (_nVaccinationDay>=0) {
    if (_nInfectiousCountdown>0) {
      if (Community::_fVE_u5_reduction==0.0 || _nAgeAtVaccination>=5)
	_fInfectiousness = _fBaseInfectiousness * (1.0-Community::_fVEI*Community::getVaccineEfficacyWaning(_nVaccinationDay));
      else
	_fInfectiousness = _fBaseInfectiousness * (1.0-Community::_fVEI*Community::getVaccineEfficacyWaning(_nVaccinationDay)*(1.0-Community::_fVE_u5_reduction));
    } else if (_fBaseSusceptibility>0.0) {
      if (Community::_fVE_u5_reduction==0.0 || _nAgeAtVaccination>=5)
	_fSusceptibility = _fBaseSusceptibility * (1.0-Community::_fVES*Community::getVaccineEfficacyWaning(_nVaccinationDay));
      else
	_fSusceptibility = _fBaseSusceptibility * (1.0-Community::_fVES*Community::getVaccineEfficacyWaning(_nVaccinationDay)*(1.0-Community::_fVE_u5_reduction));
    }
  }  
  if (_nVaccinationDay>=0) // && _nVaccinationDay<Community::getMaxVaccineDays()-1)
    _nVaccinationDay++; // days since vaccination
  return true;
}

// populate - add at least nMinPeople but no more than nMaxPeople to a new community
// will add families until at least nMinPeople are used
int Community::populate(int nMinPeople,int nMaxPeople,int *familyids,int *ages) {
  int familycount=0; // number of families in this community
  _nNumResidents=0;  // number of people in this community
  _nNumVisitors=0;
  _residents = new int[nMinPeople+Community::MAXFAMILYSIZE];
  _nMaxVisitors=nMinPeople/2+1;
  _visitors = new int[_nMaxVisitors];
  
  while (_nNumResidents<nMinPeople) {
    // new family - how large is this family?
    int familysize=1;
    for (; (_nNumResidents+familysize<nMaxPeople && familyids[_nNumResidents]==familyids[_nNumResidents+familysize]); familysize++)
      ;
    assert(familysize<=Community::MAXFAMILYSIZE);
    for (int i=0; i<familysize; i++) {
      int id = Person::getNextPersonID();
      _residents[_nNumResidents] = id;
      Person &person = Person::personArray[id];
      person.setAge(ages[_nNumResidents]); // bugbugbug????????
      person.setHomeCommunity(_nID);
      person.setWorkCommunity(_nID);
      person.setFamilySize(familysize); // we needed to know the family size in advance, so that's why we are doing this loop
      person.setFamily(familycount);
      familyArray[id] = familycount;
      _nNumResidents++;
    }
    familycount++;
  }
  return 1;
}

int Community::populate(gsl_rng *rng, int nNumPeople) {
  int familycount=0;
  int leftinfamily=0;
  int familysize=0;

  _nNumResidents=0;
  _nNumVisitors=0;
  _residents = new int[nNumPeople];
  _nMaxVisitors=nNumPeople/2+1;
  _visitors = new int[_nMaxVisitors];
  
  for (int i=nNumPeople; i>0; i--) {
    int id = Person::getNextPersonID();
    _residents[_nNumResidents++] = id;
    Person &person = Person::personArray[id];
    person.setHomeCommunity(_nID);
    person.setWorkCommunity(_nID);
    if (leftinfamily<=0) {
      double r = gsl_rng_uniform(rng);
      for (; leftinfamily < MAXFAMILYSIZE && _fFamiltySizeCDF[leftinfamily]<r; leftinfamily++)
	;
      familysize = leftinfamily+1;
      familycount++;
    } else
      leftinfamily--;
    familyArray[id] = familycount;
    person.setFamilySize(familysize);
    person.setFamily(familycount);
  }
  return 1;
}

void Community::addVisitor(int id) {
  if (_nNumVisitors+1>=_nMaxVisitors) {
    int *temp = new int[_nMaxVisitors*2];
    memcpy(temp, _visitors, sizeof(int)*_nNumVisitors);
    _nMaxVisitors*=2;
    delete [] _visitors;
    _visitors=temp;
  }
  _visitors[_nNumVisitors++] = id;
}

int Community::getNumSusceptible() {
  int total = 0;
  for (int i=0; i<_nNumResidents; i++) {
    Person &p = Person::personArray[_residents[i]];
    if (p.getSusceptibility()>0.0)
      total++;
  }
  return total;
}

int Community::getNumImmune() {
  int total = 0;
  for (int i=0; i<_nNumResidents; i++) {
    Person &p = Person::personArray[_residents[i]];
    if (p.getSusceptibility()==0.0)
      total++;
  }
  return total;
}

int Community::getNumInfectious() {
  int total = 0;
  for (int i=0; i<_nNumResidents; i++) {
    Person &p = Person::personArray[_residents[i]];
    if (p.getInfectiousness()>0.0)
      total++;
  }
  return total;
}

int Community::getNumSymptomatic() {
  int total = 0;
  for (int i=0; i<_nNumResidents; i++) {
    Person &p = Person::personArray[_residents[i]];
    if (p.isSymptomatic())
      total++;
  }
  return total;
}

// get number of newly symptomatic people from ages [agemin,agemax] inclusive
int Community::getNumNewSymptomatic(int agemin, int agemax) {
  int total = 0;
  if (agemin<0 && agemax<0) {
    for (int i=0; i<_nNumResidents; i++) {
      Person &p = Person::personArray[_residents[i]];
      if (p.isSymptomatic() && p.getInfectiousCountup()==0)
	total++;
    }
  } else {
    for (int i=0; i<_nNumResidents; i++) {
      Person &p = Person::personArray[_residents[i]];
      if (p.isSymptomatic() && p.getInfectiousCountup()==0 && p.getAge()>=agemin && p.getAge()<=agemax)
	total++;
    }
  }
  return total;
}

int Community::getCumulativeSymptomatic() {
  int total = 0;
  for (int i=0; i<_nNumResidents; i++) {
    Person &p = Person::personArray[_residents[i]];
    total += p.getSymptomaticCount();
  }
  return total;
}

int Community::getNumVaccinated() {
  int total = 0;
  for (int i=0; i<_nNumResidents; i++) {
    Person &p = Person::personArray[_residents[i]];
    if (p.isVaccinated())
      total++;
  }
  return total;
}

// vaccinate - vaccinates fraction of the population
// returns number vaccinated
int Community::vaccinate(gsl_rng *rng, double f) {
  int nNumVaccinate=0;
  for (int i=0; i<_nNumResidents; i++) {
    if (gsl_rng_uniform(rng)<f) {
      Person &p = Person::personArray[_residents[i]];
      if (!p.isVaccinated()) {
	p.vaccinate();
	nNumVaccinate++;
	_nNumVaccinesUsed++;
      }
    }
  }
  return nNumVaccinate;
}

int Community::prevaccinate(gsl_rng *rng, double f) {
  int nNumVaccinate=0;
  for (int i=0; i<_nNumResidents; i++) {
    if (gsl_rng_uniform(rng)<f) {
      Person &p = Person::personArray[_residents[i]];
      if (!p.isVaccinated()) {
	p.prevaccinate();
	nNumVaccinate++;
	_nNumVaccinesUsed++;
      }
    }
  }
  return nNumVaccinate;
}

// getVaccineEfficacyWaning - returns the relative efficacy of vaccine "days" after vaccination
double Community::getVaccineEfficacyWaning(int days) {
  if (_fVE_linearwaning<=0.0) {
    return 1.0;
  } else {
    double temp=1.0-_fVE_linearwaning*days;
    if (temp<0.0)
      temp=0.0;
    return temp;
  }
} 

// setHygiene - sets the hygiene levels of the entire population
// returns population size
int Community::setHygiene(double f) {
  for (int i=0; i<_nNumResidents; i++) {
    Person &p = Person::personArray[_residents[i]];
    p.setHygiene(f);
  }
  return _nNumResidents;
}

// infect
// n is the total number of new infections in an unvaccinated population.
// vaccinated people have lower chance of being infected.
// returns actual number of people infected.
int Community::infect(gsl_rng *rng, int n, bool bMakeSymptomatic) {
  int count=0;
  if (_nNumResidents>0) {
    for (int i=0; i<n; i++) {
      int id = _residents[gsl_rng_uniform_int(rng, _nNumResidents)]; // person to try infecting
      Person &p = Person::personArray[id]; 
      if (bMakeSymptomatic) {
	p.makesymptomatic(rng);
	count++;
      } else {
	double s = p.getSusceptibility();
	if (s==1.0 || (s>0.0 && gsl_rng_uniform(rng)<s)) {
	  p.infect(rng);
	  count++;
	}
      }
    }
  }
  return count;
}

// poop
// update amount of Vibrio in the environment
// returns the amount of new poop
double Community::poop(gsl_rng *rng, double rainmultiplier) {
  //  cerr << "poop " << _nNumResidents << "," << _nID << endl;
  double hf = (1.0-_fWorkTimeFraction); // time spent at home
  _fVibrio*=_fVibrioDecay;              // level of Vibrio goes down
  _fHyperVibrio=0.0;                    // hyperinfectious vibrio disappears each day
  double newvibrio = 0.0;
  for (int i=0; i<_nNumResidents; i++) {
    assert(_residents[i]<Person::getLastPersonID());
    Person &p = Person::personArray[_residents[i]];
    //    cerr << "comm " << _nID << ", person " << _residents[i] << ": " << i <<  " / " << _nNumResidents << endl;
    if (p.getInfectiousness()>0.0) {
      if (p.getWorkCommunity()==_nID) { // are they here all day?
	newvibrio += p.getInfectiousness() * rainmultiplier;
	_fHyperVibrio += p.getInfectiousness() * rainmultiplier;
      } else {
	newvibrio += hf*p.getInfectiousness() * rainmultiplier;
	_fHyperVibrio += hf*p.getInfectiousness() * rainmultiplier;
      }
    }
  }

  for (int i=0; i<_nNumVisitors; i++) {
    assert(_visitors[i]<Person::getLastPersonID());
    Person &p = Person::personArray[_visitors[i]];
    if (p.getWorkCommunity()==_nID) { // make sure they are still working here
      newvibrio += _fWorkTimeFraction*p.getInfectiousness() * rainmultiplier;
      _fHyperVibrio += _fWorkTimeFraction*p.getInfectiousness() * rainmultiplier;
    }
  }
  if (_bRiver) {
    newvibrio *= _fRiverSheddingMultiplier;
    _fHyperVibrio *= _fRiverSheddingMultiplier;
  }
  _fVibrio += newvibrio;
  return newvibrio;
}

// drink
// infect people (from environment and households)
void Community::drink(gsl_rng *rng) {
  double v = (_fVibrio + _fHyperVibrioMultiplier*_fHyperVibrio)/_nNumResidents + _fRiverVibrio + _fHyperVibrioMultiplier*_fRiverHyperVibrio;
  assert(_fRiverVibrio>=0.0);
  assert(_fRiverHyperVibrio>=0.0);
  double p = _fVibrioBeta * v/(_fVibrio50 + v);

  // transmission from environment
  //  if (p>0.0) {
  if (p>(1.0/gsl_rng_max(rng))) { // probability must be greater than lower limit of random number generator
    double ph = p*(1.0-_fWorkTimeFraction);
    double pw = p*_fWorkTimeFraction;
    for (int i=0; i<_nNumResidents; i++) {
      Person &person = Person::personArray[_residents[i]];
      if (person.getSusceptibility()>0.0 &&
	  person.getTemporaryImmunityCountdown()<0 && 
	  gsl_rng_uniform(rng)<(person.getWorkCommunity()==_nID?p:ph)*person.getSusceptibility()*(1.0-person.getHygiene()))
	person.infect(rng);
    }
    for (int i=0; i<_nNumVisitors; i++) {
      Person &person = Person::personArray[_visitors[i]];
      if (person.getWorkCommunity()==_nID &&
	  person.getSusceptibility()>0.0 && 
	  person.getTemporaryImmunityCountdown()<0 && 
	  gsl_rng_uniform(rng)<pw*person.getSusceptibility()*(1.0-person.getHygiene()))
	person.infect(rng);
    }
  }

  // transmission within families
  if (Community::_fHouseholdContactProbability>0.0) {
    for (int i=0; i<_nNumResidents; i++) {
      Person &person = Person::personArray[_residents[i]];
      // find infected person
      if (person.getInfectiousness()>0.0) {
	// find other infected people in the family
	double escapeprob = 1.0;
	int family = person.getFamily();
	int familylast; // id+1 of last person in family
	int familyfirst;// id of first person in family
	for (familylast=i; familylast<_nNumResidents && familylast<i+person.getFamilySize(); familylast++) {
	  Person &person2 = Person::personArray[_residents[familylast]];
	  if (person2.getFamily()!=family)
	    break;
	  if (person2.getInfectiousness()>0.0)
	    escapeprob *= (1.0-Community::_fHouseholdContactProbability*person2.getInfectiousness());
	}
	for (familyfirst=i-1; familyfirst>=0 && familyfirst>i-person.getFamilySize(); familyfirst--) {
	    Person &person2 = Person::personArray[_residents[familyfirst]];
	    if (person2.getFamily()!=family)
	      break;
	}
	familyfirst++;
	// expose all susceptibles in family
	double infectprob = 1.0-escapeprob;
	//	cerr << " person: " << person.getID() << ", family size: " << person.getFamilySize() << ",  p=" << Community::_fHouseholdContactProbability << ", infprob=" << infectprob << endl;
	if (infectprob>0.0) {
	  for (int j=familyfirst; j<familylast; j++) {
	    Person &person2 = Person::personArray[_residents[j]];
	    if (person2.getSusceptibility()>0.0 && 
		person.getTemporaryImmunityCountdown()<0 && 
		gsl_rng_uniform(rng)<infectprob*person2.getSusceptibility()*(1.0-person2.getHygiene()))
	      person2.infect(rng);
	  }
	}
	// go to next family
	i = familylast;
      }
    }
  }
}

// tick - advance counters for each resident
int Community::tick(gsl_rng *rng) {
  //  cerr << "tick" << endl;
  int count=0;
  for (int i=0; i<_nNumResidents; i++) {
    Person &p = Person::personArray[_residents[i]];
    p.step(rng);
    if (p.getInfectiousCountup()==0) {
      count++;
    } else if (p.getDaysInfectious()==0 && p.isSymptomatic() &&
	       p.getHomeCommunity()!=p.getWorkCommunity() &&
	       gsl_rng_uniform(rng)<fWithdrawFraction) {
      // this worker withdraws
      // equivalent to turning the worker into a non-worker
      p.setWorkCommunity(p.getHomeCommunity());
    }
  }
  return count;
}
