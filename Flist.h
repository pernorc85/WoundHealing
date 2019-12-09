#pragma once
#ifndef _FLIST_H_
#define _FLIST_H_

#include <string>
#include <cmath>
#include <iostream> 
#include <list>
#include <map>
#include "nr.h"
using namespace std;

struct Fibroblast{
    Fibroblast(DP x_,DP y_, DP theta_):
      xx(x_), yy(y_), theta(theta_), activity(0.5) {}

    DP xx,yy; 
    DP theta;
    DP activity;
    DP speed;
};

class ECM;
class Chemokine;

class Flist{
public:
    Flist(int xstep, int ystep);
    void initialize(int FDensity, double a, double b);
    void initialize_rectangle(int FDensity, double a, double b);
    void initialize(int FDensity, int wound_radius);
    void initialize_oval(int FDensity);
 
    void Flist_move(Chemokine& PDGF,
                    ECM& extraCellularMatrix, 
                    const Mat_DP& tissue_displacement_x,
                    const Mat_DP& tissue_displacement_y, DP time_step); 
    void Fcell_activity_dynamics(Fibroblast* curPtr, DP conc, DP time_step);
    void Fcell_move(Fibroblast* curPtr, DP gradx, DP grady, DP fdensity, DP cdensity, 
                    Mat_DP& collagen_angle, Mat_DP& F, DP time_step);
    void output_fibroblast();
    void output_flistFile();

protected:
    std::map<int, Fibroblast> mFCellMap;
    DP speed;
    int mXstep, mYstep;

    friend class ECM;
};

#endif
