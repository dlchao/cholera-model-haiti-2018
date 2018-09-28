/*
 * CholeraModel.cpp
 */

#include <iostream>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "CholeraModel.h"

using namespace std;

int CholeraModel::_nMinInfectiousDays = 7;
int CholeraModel::_nMaxInfectiousDays = 14;
double CholeraModel::_fSymptomaticFraction = 0.1;
double CholeraModel::_fSymptomaticInfectiousnessMultiplier = 10.0;
double CholeraModel::_VES = 0.5;
const int CholeraModel::_nMaxIncubationDays = 5;
double CholeraModel::_fIncubationCDF[CholeraModel::_nMaxIncubationDays] = {.4,.4,.07,.07,.06};
const int CholeraModel::_nMaxVaccineDays = 21;
double CholeraModel::_fVaccineBuildup[CholeraModel::_nMaxVaccineDays] = 
  {0.01, 0.02, 0.03, 0.04, 0.05, 
   0.06, 0.07, 0.08, 0.09, 0.1,
   0.5,0.5,0.5,0.5,0.5,
   0.5,0.5,0.5,0.5,1.0,
   1.0};
int CholeraModel::_fIncubationCDF100[100] = 
  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
   1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
   1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
   2,2,2,2,2,2,2,3,3,3,3,3,3,3,4,4,4,4,4,4};

CholeraModel::CholeraModel() { 
  _nNumSusceptible = 0;
  _nNumInfectedSymptomatic = new int[_nMaxInfectiousDays];
  _nNumInfectedAsymptomatic = new int[_nMaxInfectiousDays];
  memset(_nNumInfectedSymptomatic, 0, sizeof(int)*_nMaxInfectiousDays);
  memset(_nNumInfectedAsymptomatic, 0, sizeof(int)*_nMaxInfectiousDays);
  memset(_nNumIncubating, 0, sizeof(int)*_nMaxIncubationDays);
  _nNumRecovered = 0;
  _p=0.0;
  _nNumSusceptibleVaccinatedTotal=0;
  memset(_nNumSusceptibleVaccinated, 0, sizeof(int)*_nMaxVaccineDays);
}

CholeraModel::~CholeraModel() { 
  delete _nNumInfectedSymptomatic;
  delete _nNumInfectedAsymptomatic;
}

int CholeraModel::getNumInfected() {
  int total=0;
  for (int i=0; i<_nMaxInfectiousDays; i++)
    total+=_nNumInfectedSymptomatic[i]+_nNumInfectedAsymptomatic[i];
  return total;
}

int CholeraModel::getNumInfectedSymptomatic() {
  int total=0;
  for (int i=0; i<_nMaxInfectiousDays; i++)
    total+=_nNumInfectedSymptomatic[i];
  return total;
}

int CholeraModel::getNumInfectedAsymptomatic() {
  int total=0;
  for (int i=0; i<_nMaxInfectiousDays; i++)
    total+=_nNumInfectedAsymptomatic[i];
  return total;
}

// vaccinate
// vaccinate n susceptibles
// should change this to susceptibles + exposed!!!!!!!!!!
int CholeraModel::vaccinate(gsl_rng *rng, int n) {
  if (_nNumSusceptibleVaccinatedTotal+n>_nNumSusceptible)
    n=_nNumSusceptible-_nNumSusceptibleVaccinatedTotal;
  _nNumSusceptibleVaccinatedTotal+=n;
  _nNumSusceptibleVaccinated[0]+=n;
  return 0;
}

// infect
// n is the total number of new infections in an unvaccinated population.
// vaccinated people have lower chance of being infected.
// returns actual number of people infected.
int CholeraModel::infect(gsl_rng *rng, int n) {
  int infectunvaccinated = n;
  int oldsusceptible = _nNumSusceptible;
  if (_nNumSusceptibleVaccinatedTotal>0) {
    infectunvaccinated = gsl_ran_binomial(rng, 
					  1.0-_nNumSusceptibleVaccinatedTotal/(double)_nNumSusceptible,
					  n);
  }
  // infect unvaccinated
  for (int i=0; i<infectunvaccinated; i++) {
    int firstday = _fIncubationCDF100[gsl_rng_uniform_int(rng, 100)];
    _nNumIncubating[firstday]++;
  }
  _nNumSusceptible-=infectunvaccinated;

  // try to infect vaccinated
  int infectvaccinated = n-infectunvaccinated;
  for (int i=0; i<infectvaccinated; i++) {
    int person = gsl_rng_uniform_int(rng, _nNumSusceptibleVaccinatedTotal);
    int sum=0;
    int day=0;
    for (; day<_nMaxVaccineDays; day++) {
      sum+=_nNumSusceptibleVaccinated[day];
      if (sum>=person)
	break;
    }
    if (gsl_rng_uniform(rng)<(1.0-_VES*_fVaccineBuildup[day])) { // vaccinated person infected
      int firstday = _fIncubationCDF100[gsl_rng_uniform_int(rng, 100)];
      _nNumIncubating[firstday]++;
      _nNumSusceptible--;
      _nNumSusceptibleVaccinatedTotal--;
    }
  }
  
  return oldsusceptible-_nNumSusceptible;
}

double CholeraModel::getFOI() {
  double p = 1.0-
    pow(1.0-_p*_fSymptomaticInfectiousnessMultiplier,getNumInfectedSymptomatic())*
    pow(1.0-_p,getNumInfectedAsymptomatic());
  return p;
}

int CholeraModel::step(gsl_rng *rng) {
  return step(rng, 0.0);
}

int CholeraModel::step(gsl_rng *rng, double externalFOI) {
  // calculate infection probability
  double p = getFOI();
  if (externalFOI>0.0)
    p = 1.0-(1.0-p)*(1.0-externalFOI);

  // advance infection compartments
  _nNumRecovered+=_nNumInfectedSymptomatic[_nMaxInfectiousDays-1]+
    _nNumInfectedAsymptomatic[_nMaxInfectiousDays-1];
  for (int i=_nMaxInfectiousDays-1; i>=1; i--) {
    _nNumInfectedSymptomatic[i]=_nNumInfectedSymptomatic[i-1];
    _nNumInfectedAsymptomatic[i]=_nNumInfectedAsymptomatic[i-1];
  }
  _nNumInfectedAsymptomatic[0]=0;
  _nNumInfectedSymptomatic[0]=0;

  // incubation -> infected
  if (_nNumIncubating[_nMaxIncubationDays-1]>0) {
    int asymptomatic = 0; // number of new asymptomatic infections
    if (_fSymptomaticFraction<1.0)
      asymptomatic = gsl_ran_binomial(rng, 
				      1.0-_fSymptomaticFraction,
				      _nNumIncubating[_nMaxIncubationDays-1]);
    if (_nMinInfectiousDays!=_nMaxInfectiousDays) {
      int firstdayrange = _nMaxInfectiousDays-_nMinInfectiousDays+1;
      for (int i=0; i<_nNumIncubating[_nMaxIncubationDays-1]; i++) {
	int firstday = gsl_rng_uniform_int(rng, firstdayrange);
	if (asymptomatic>0) {
	  _nNumInfectedAsymptomatic[firstday]++;
	  asymptomatic--;
	} else 
	  _nNumInfectedSymptomatic[firstday]++;
      }
    } else {
      _nNumInfectedSymptomatic[0]+=_nNumIncubating[_nMaxIncubationDays-1]-asymptomatic;
      _nNumInfectedAsymptomatic[0]+=asymptomatic;
    }
  }

  // advance incubation compartments
  for (int i=_nMaxIncubationDays-1; i>=1; i--)
    _nNumIncubating[i]=_nNumIncubating[i-1];
  _nNumIncubating[0]=0;

  // infect susceptibles
  if (p>0.0) {
    if (_nNumSusceptible-_nNumSusceptibleVaccinatedTotal>0) {
      int infected = gsl_ran_binomial(rng, 
				      p,
				      _nNumSusceptible-_nNumSusceptibleVaccinatedTotal);
      infect(rng, infected);
    }
  }

  // advance vaccination compartments
  _nNumSusceptibleVaccinated[_nMaxVaccineDays-1] += _nNumSusceptibleVaccinated[_nMaxVaccineDays-2];
  for (int i=_nMaxVaccineDays-2; i>=1; i--)
    _nNumSusceptibleVaccinated[i]=_nNumSusceptibleVaccinated[i-1];
  _nNumSusceptibleVaccinated[0]=0;

  return 0;
}
