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
    Flist(int xstep, int ystep, int FNinit);
        
    void Flist_move(Chemokine& PDGF,
                    ECM& extraCellularMatrix); 
    void Fcell_move(std::list<Fibroblast>::iterator curPtr, DP gradx, DP grady, DP fdensity, DP cdensity, 
                    Mat_DP& collagen_angle);
    void output_fibroblast();
    void output_flistFile();

protected:
    std::list<Fibroblast> mFCells;
    DP speed;
    int mXstep, mYstep;

    friend class ECM;
};

#endif
