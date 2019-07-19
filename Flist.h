#pragma once
#ifndef _FLIST_H_
#define _FLIST_H_

#include <string>
#include <cmath>
#include <iostream> 
#include <list>
#include "nr.h"
using namespace std;

struct Fibroblast{
    Fibroblast(DP x_,DP y_, DP theta_):
      xx(x_), yy(y_), theta(theta_){}

    DP xx,yy; 
    DP theta;
    DP speed;
};

class ECM;
class Chemokine;

class Flist{
public:
    Flist(int xstep, int ystep);
    void initialize(int FDensity, double a, double b);
    void initialize(int FDensity, int wound_radius);
 
    void Flist_move(Chemokine& PDGF,
                    ECM& extraCellularMatrix, DP time_step); 
    void Fcell_move(std::list<Fibroblast>::iterator curPtr, DP gradx, DP grady, DP fdensity, DP cdensity, 
                    Mat_DP& collagen_angle, DP time_step);
    void output_fibroblast();
    void output_flistFile();

protected:
    std::list<Fibroblast> mFCells;
    DP speed;
    int mXstep, mYstep;

    friend class ECM;
};

#endif
