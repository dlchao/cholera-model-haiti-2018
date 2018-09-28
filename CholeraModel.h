/*
 * CholeraModel.h
 * 10/2010
 */
#ifndef __CHOLERAMODEL_H
#define __CHOLERAMODEL_H

#include <math.h> 
#include <string.h> 
#include <gsl/gsl_rng.h>
#include "EpiModel.h"

class CholeraModel: public EpiModel {
 public:
  CholeraModel();
  ~CholeraModel();
  void setNumPeople(int n) { _nNumPeople=n; _nNumSusceptible=n;}
  int getNumSusceptible() { return _nNumSusceptible; } // residents only
  int getNumRecovered() { return _nNumRecovered; }
  int getNumInfected();  // actually, returns number infectious
  int getNumInfectedSymptomatic();
  int getNumInfectedAsymptomatic();
  int vaccinate(gsl_rng *rng, int n);
  int infect(gsl_rng *rng, int n);
  int step(gsl_rng *rng);
  int step(gsl_rng *rng, double externalFOI);
  double getFOI();  // force of infection in this community
  void setSymptomaticFraction(double f) { _fSymptomaticFraction=f; }
  double getSymptomaticFraction() { return _fSymptomaticFraction; }
  void setSymptomaticInfectiousnessMultiplier(double f) { _fSymptomaticInfectiousnessMultiplier=f; }
  double getSymptomaticInfectiousnessMultiplier() { return _fSymptomaticInfectiousnessMultiplier; }
  void setInfectionProbability(double f) { _p=f; }
 protected:
  int _nNumSusceptible;
  int _nNumRecovered;
  int *_nNumInfectedSymptomatic;
  int *_nNumInfectedAsymptomatic;
  double _p;
  int _nNumIncubating[20];

  // vaccinated states
  int _nNumSusceptibleVaccinated[50]; // # of susceptibles vaccinated n days ago
  int _nNumSusceptibleVaccinatedTotal;// total number of susceptibles vaccinated
  int _nNumIncubatingVaccinated[20];
  int *_nNumInfectedSymptomaticVaccinated;
  int *_nNumInfectedAsymptomaticVaccinated;

  static double _fWorkerInfectiousnessMultiplier; // fraction of time that a commuting working person's time is in the work community
  static double _fSymptomaticFraction;  // fraction of infected people who become symptomatic
  static double _fSymptomaticInfectiousnessMultiplier; // symptomatic people are more infectious
  static int _nMinInfectiousDays;  // minimum number of days that an infected person is infectious
  static int _nMaxInfectiousDays;  // maximum number of days that an infected person is infectious
  const static int _nMaxVaccineDays; // number of days for vaccine to reach full efficacy
  static double _fVaccineBuildup[]; // relative vaccine efficacy after n days
  const static int _nMaxIncubationDays;
  static double _fIncubationCDF[];
  static int _fIncubationCDF100[];
  static double _VES;
};
#endif
