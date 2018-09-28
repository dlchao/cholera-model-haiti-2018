/*
 * GridCells.cpp
 */

#include <iostream>
#include <assert.h>
#include "GridCells.h"

using namespace std;

const double HexGridCells::HEXGRIDVERTICALSCALINGFACTOR = 0.75;

int SquareGridCells::getNeighbors(int nID, int *buf, int bufsize) {
  int x = nID%_nXSize;
  int y = nID/_nXSize;
  int count = 0;
  for (int i=-1; i<=1; i++) {
    for (int j=-1; j<=1; j++) {
      if ((i!=0 || j!=0) && 
	  (_bWrapX || (i+x>=0 && i+x<_nXSize)) && 
	  (_bWrapY || (j+y>=0 && j+y<_nYSize))) {
	int a = (x+i+_nXSize)%_nXSize + ((y+j+_nYSize)%_nYSize) *_nXSize;
	    //	    if (count+1>bufsize)
	    //	      return false;
	buf[count++] = a;
      }
    }
  }
  return count;
}

int HexGridCells::getNeighbors(int nID, int *buf, int bufsize) {
  int x = nID%_nXSize;
  int y = nID/_nXSize;
  int count = 0;
  if (!_bWrapX && !_bWrapY) {
    if (y>0) {
      buf[count++] = x + (y-1)*_nXSize;
      if (y%2==0) {
	if (x>0)
	  buf[count++] = (x-1) + (y-1)*_nXSize;
      } else if (y%2==1 && x+1<_nXSize) 
	buf[count++] = (x+1) + (y-1)*_nXSize;
    }
    if (x>0)
      buf[count++] = (x-1) + y*_nXSize;
    if (x+1<_nXSize)
      buf[count++] = (x+1) + y*_nXSize;
    if (y+1<_nYSize) {
      buf[count++] = x + (y+1)*_nXSize;
      if (y%2==0) {
	if (x>0)
	  buf[count++] = (x-1) + (y+1)*_nXSize;
      } else if (y%2==1 && x+1<_nXSize)
	buf[count++] = (x+1) + (y+1)*_nXSize;
    }
  } else if (_bWrapX && !_bWrapY) {
    int nextx = (x+1)%_nXSize;
    int prevx = ((x-1)+_nXSize)%_nXSize;
    int nexty = (y+1)%_nYSize;
    int prevy = ((y-1)+_nYSize)%_nYSize;
    if (y>0) {
      buf[count++] = x + prevy*_nXSize;
      if (y%2==0) {
	buf[count++] = prevx + prevy*_nXSize;
      } else if (y%2==1)
	buf[count++] = nextx + prevy*_nXSize;
    }
    buf[count++] = prevx + y*_nXSize;
    buf[count++] = nextx + y*_nXSize;
    if (nexty<_nYSize) {
      buf[count++] = x + nexty*_nXSize;
      if (y%2==0) {
	buf[count++] = prevx + nexty*_nXSize;
      } else if (y%2==1)
	buf[count++] = nextx + nexty*_nXSize;
    }
  } else if (_bWrapX && _bWrapY) {
    int nextx = (x+1)%_nXSize;
    int prevx = ((x-1)+_nXSize)%_nXSize;
    int nexty = (y+1)%_nYSize;
    int prevy = ((y-1)+_nYSize)%_nYSize;
    buf[count++] = x + prevy*_nXSize;
    if (y%2==0) {
      buf[count++] = prevx + prevy*_nXSize;
    } else if (y%2==1)
      buf[count++] = nextx + prevy*_nXSize;
    buf[count++] = prevx + y*_nXSize;
    buf[count++] = nextx + y*_nXSize;
    buf[count++] = x + nexty*_nXSize;
    if (y%2==0) {
      buf[count++] = prevx + nexty*_nXSize;
    } else if (y%2==1)
      buf[count++] = nextx + nexty*_nXSize;
  }
  return count;
}
