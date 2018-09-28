/* class GridCells
 *
 * Dennis Chao
 * 10/2010
 * Based on Dennis Chao's CryptSim classes, available at:
 * http://www.cs.unm.edu/~dlchao/cancer/index.html
 */

#ifndef __GRIDCELLS_H
#define __GRIDCELLS_H

class GridCells {
 public:
  GridCells() {
    _nNumCells=0;
  }
  virtual ~GridCells() {}
  int getSize() { return _nNumCells; }  // returns the # of cells
  virtual double getMaxX()=0;           // width of 2d plane
  virtual double getMaxY()=0;           // height of 2d plane
  virtual double getX(int n)=0;         // returns the x-coordinate of cell n on a plane
  virtual double getY(int n)=0;         // returns the y-coordinate of cell n on a plane
  virtual int getNeighbors(int nID, int *buf, int bufsize)=0;
 protected:
  int _nNumCells;
  bool _bWrapX;  // cyclic boundary on X axis
  bool _bWrapY;  // cyclic boundary on Y axis
};

class SquareGridCells: public GridCells {
 public:
  SquareGridCells(int nXSize, int nYSize, bool bWrapX=false, bool bWrapY=false) {
    _nXSize = nXSize;
    _nYSize = nYSize;
    _bWrapX = bWrapX;
    _bWrapY = bWrapY;
    _nNumCells=_nXSize*_nYSize;
  }
  ~SquareGridCells() {}
  double getMaxX()   { return _nXSize-1; }       // width of 2d plane
  double getMaxY()   { return _nYSize-1; }       // height of 2d plane
  double getX(int n) { return n%_nXSize; }       // returns the x-coordinate of cell n on a plane
  double getY(int n) { return n/_nXSize; }       // returns the y-coordinate of cell n on a plane
  int _nXSize,_nYSize;

  int getNeighbors(int nID, int *buf, int bufsize);
};

/* HexGridCells - hexagonal lattice

An example: of a 9xn hexagonal lattice:

columns:  1   2   3   4   5   6   7   8   9
       / \ / \ / \ / \ / \ / \ / \ / \ / \ / 
      |   |   |   |   |   |   |   |   |   |  
       \ / \ / \ / \ / \ / \ / \ / \ / \ / \ 
        |   |   |   |   |   |   |   |   |   |
       / \ / \ / \ / \ / \ / \ / \ / \ / \ / 
row 3 |   |   |   |   |   |   |   |   |   |  
       \ / \ / \ / \ / \ / \ / \ / \ / \ / \ 
row 2   |   |   |   |   |   |   |   |   |   |
       / \ / \ / \ / \ / \ / \ / \ / \ / \ / 
row 1 |   |   |   |   |   |   |   |   |   |  
       \ / \ / \ / \ / \ / \ / \ / \ / \ / \ 
columns:1   2   3   4   5   6   7   8   9
*/
class HexGridCells: public GridCells {
 public:
  HexGridCells(int nXSize, int nYSize, bool bWrapX=false, bool bWrapY=false) {
    _nXSize = nXSize;
    _nYSize = nYSize;
    _bWrapX = bWrapX;
    _bWrapY = bWrapY;
    _nNumCells=_nXSize*_nYSize;
  }
  ~HexGridCells() {}
  double getMaxX()   { return (_nXSize*2.0 + 1.0)/2.0; }       // width of 2d plane
  double getMaxY()   { return _nYSize*HEXGRIDVERTICALSCALINGFACTOR; }       // height of 2d plane
  double getX(int n) {
    double retval = n%_nXSize;
    if ((n/_nXSize)%2==1)
      retval += 0.5; 
    return  retval;
  }
  double getY(int n) { 
    int y = n/_nXSize;
    return y*HEXGRIDVERTICALSCALINGFACTOR;
  }
  int _nXSize,_nYSize;
  const static double HEXGRIDVERTICALSCALINGFACTOR;
  
  int getNeighbors(int nID, int *buf, int bufsize);
};
#endif
