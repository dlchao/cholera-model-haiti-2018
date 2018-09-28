/*
 * R0Community.h
 * 12/2010
 */
#ifndef __R0COMMUNITY_H
#define __R0COMMUNITY_H

#include <math.h> 
#include <string.h> 
#include <gsl/gsl_rng.h>
#include "Community.h"

using namespace std;

class R0Community : public Community {
 public:
  R0Community() {
  }
  virtual ~R0Community() {
  }
  virtual void drink(gsl_rng *rng);
};
#endif
