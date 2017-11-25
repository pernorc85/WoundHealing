//#pragma once
#ifndef _ECM_H_
#define _ECM_H_

#include <string>
#include <cmath>
#include <iostream> 
#include "nr.h"
using namespace std;

const DP tlength = 0.15; //unit is hr

class Flist;

class ECM{
    public:
        ECM(int ystep,int xstep):
		ECM_ystep(ystep),ECM_xstep(xstep),
		collagen(ystep/10,xstep/10),collagen_density(ystep/5,xstep/5),
		fibronectin_density(ystep/5,xstep/5)
        {};
        void initiate(int wound_radius);
        void collagen_orientation(Flist& fibroblastList);
        void output_collagen(Flist& fibroblastList);

    protected:
        int ECM_ystep;
        int ECM_xstep;
        Mat_DP collagen;
        Mat_DP collagen_density;
        Mat_DP fibronectin_density;
        
    friend class Flist;
};



#endif
