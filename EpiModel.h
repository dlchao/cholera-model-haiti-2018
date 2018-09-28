/*
 * EpiModel.h
 * 10/2010
 */
#ifndef __EPIMODEL_H
#define __EPIMODEL_H

#include <string.h> 
#include <gsl/gsl_rng.h>

class EpiModel {
 public:
  EpiModel() {  }
  virtual ~EpiModel() {  }
  virtual void setNumPeople(int n) { _nNumPeople=n; }
  int getNumPeople() { return _nNumPeople; }
  virtual int getNumInfected()=0;
  virtual int getNumSusceptible()=0;
  virtual int getNumRecovered()=0;
  virtual int infect(gsl_rng *rng, int n)=0; // returns the number infected
  virtual int step(gsl_rng *rng)=0; // returns force of infection
 protected:
  int _nNumPeople;
};

class SIRModel: public EpiModel {
 public:
  SIRModel() {
    _nNumPeople=0;
    _nNumSusceptible=0;
    _nNumInfected=NULL;
    _nNumRecovered=0;
    _nMaxInfectiousDays=0;
    _p=0;
  }
  SIRModel(int nNumPeople, int nNumDaysInfected) {
    _nNumPeople=nNumPeople;
    _nNumSusceptible=0;
    _nNumRecovered=0;
    _nMaxInfectiousDays=0;
    _nMaxInfectiousDays = nNumDaysInfected;
    _nNumInfected=new int[_nMaxInfectiousDays];
    _p = 0.01;
  }
  ~SIRModel() { delete _nNumInfected; }
  void setNumPeople(int n) { _nNumPeople=n; _nNumSusceptible=n;}
  void setDaysInfected(int n) { _nMaxInfectiousDays=n;     _nNumInfected=new int[_nMaxInfectiousDays];
}
  void setInfectionProbability(double f) { _p=f; }
  int getNumInfected();
  int getNumInfected(int day);
  int getNumSusceptible() { return _nNumSusceptible; }
  int getNumRecovered() { return _nNumRecovered; }
  int infect(gsl_rng *rng, int n);
  int step(gsl_rng *rng);
 protected:
  int _nNumSusceptible;
  int *_nNumInfected;
  int _nMaxInfectiousDays;
  int _nNumRecovered;
  double _p;
};
#endif
