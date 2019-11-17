//#pragma once
#ifndef _ECM_H_
#define _ECM_H_

#include <string>
#include <cmath>
#include <iostream> 
#include <vector>
#include "nr.h"
using namespace std;

class Flist;

class ECM{
public:
    ECM(int ystep,int xstep):
            ECM_ystep(ystep),ECM_xstep(xstep),
            collagen(ystep/10, xstep/10),
            tension_theta_field(ystep/10, xstep/10),
            collagen_density(ystep/5, xstep/5),
            fibronectin_density(ystep/5, xstep/5),
            stretch_history(ystep/10, xstep/10){
                for(size_t i = 0; i < 18; i++) {
                    tensiontheta_bins_history.emplace_back(Mat_DP(ystep/10, xstep/10));
                }
            };
    void initialize(int wound_radius);
    void initialize(double a, double b);
    void initialize_rectangle(double a, double b);
    void initialize_oval();
    void collagen_orientation(Flist& fibroblastList, const Mat_DP& tdx, const Mat_DP& tdy, double time_step);
    void collagen_orientation_with_fibroblast(Flist& fList, double time_step);
    void output_collagen(Flist& fibroblastList);
    void output_tension();
    void detect_tension_field_singularity();

public:
    Mat_DP collagen;
    Mat_DP collagen_density;
    Mat_DP tension_theta_field;

protected:    
    Mat_DP fibronectin_density;
    Mat_DP stretch_history;
    vector<Mat_DP> tensiontheta_bins_history;//data intensive
    vector<pair<int, int>> mTensionFieldSingularPoints;
    int ECM_ystep;
    int ECM_xstep;
    friend class Flist;
    friend class Elist;

private: 
    void collagen_orientation_under_tension(const Mat_DP& tdx, const Mat_DP& tdy, double time_step);
    void collagen_production(Flist& fList, double time_step);
};


#endif
