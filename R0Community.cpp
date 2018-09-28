/*
 * R0Community.cpp
 * 12/2010
 * A community where there is no transmission.
 * Derived from Community.
 * The "drink" method is copied from Community, but rather
 * than infect people, people become immediately immune.
 */

#include <iostream>
#include <math.h>
#include <assert.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "R0Community.h"

using namespace std;

// drink
// infect people (from environment and households)
void R0Community::drink(gsl_rng *rng) {
  double v = (_fVibrio + _fHyperVibrioMultiplier*_fHyperVibrio)/_nNumResidents + _fRiverVibrio + _fHyperVibrioMultiplier*_fRiverHyperVibrio;
  assert(_fRiverVibrio>=0.0);
  assert(_fRiverHyperVibrio>=0.0);
  double p = _fVibrioBeta * v/(_fVibrio50 + v);

  // transmission from environment
  if (p>0.0) {
    double ph = p*(1.0-_fWorkTimeFraction);
    double pw = p*_fWorkTimeFraction;
    for (int i=0; i<_nNumResidents; i++) {
      Person &person = Person::personArray[_residents[i]];
      if (person.getSusceptibility()>0.0 && gsl_rng_uniform(rng)<(person.getWorkCommunity()==_nID?p:ph)*person.getSusceptibility()) {
	
	person.makeImmune();/////////////////
      }
    }
    for (int i=0; i<_nNumVisitors; i++) {
      Person &person = Person::personArray[_visitors[i]];
      if (person.getWorkCommunity()==_nID &&
	person.getSusceptibility()>0.0 && gsl_rng_uniform(rng)<pw*person.getSusceptibility())
	person.makeImmune();//////////
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
		gsl_rng_uniform(rng)<infectprob*person2.getSusceptibility())
	      person2.makeImmune();///////////
	  }
	}
	// go to next family
	i = familylast;
      }
    }
  }
}
