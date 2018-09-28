/*
 * EpiModel.cpp
 */

#include <iostream>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "EpiModel.h"

using namespace std;

int SIRModel::getNumInfected() {
  int total=0;
  for (int i=0; i<_nMaxInfectiousDays; i++)
    total+=_nNumInfected[i];
  return total;
}

int SIRModel::getNumInfected(int day) {
  return _nNumInfected[day];
}

int SIRModel::infect(gsl_rng *rng, int n) {
  if (_nMaxInfectiousDays>0)
    _nNumInfected[0] += n;
  _nNumSusceptible-=n;
  return n;
}

int SIRModel::step(gsl_rng *rng) {
  // calculate infection probability
  double p = 1.0-pow(1.0-_p,getNumInfected());

  // advance infection compartments
  _nNumRecovered+=_nNumInfected[_nMaxInfectiousDays-1];
  for (int i=_nMaxInfectiousDays-1; i>=1; i--)
    _nNumInfected[i]=_nNumInfected[i-1];
  _nNumInfected[0]=0;

  // infect susceptibles
  if (p>0.0 && _nNumSusceptible>0) {
    int infected = gsl_ran_binomial(rng, 
				    p,
				    _nNumSusceptible);
    _nNumSusceptible -= infected;
    _nNumInfected[0] += infected;
  }

  return 0;
}
