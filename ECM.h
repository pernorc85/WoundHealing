//#pragma once
#ifndef _ECM_H_
#define _ECM_H_

#include <string>
#include <cmath>
#include <iostream> 
#include "nr.h"
#include "Flist.h"
using namespace std;

const DP tlength = 0.15; //unit is hr
const double PI = 3.1415926535897932;

typedef struct Fibroblast Fcell;
typedef Fcell* FcellPtr;
struct Fibroblast{
    DP xx,yy; 
    DP theta; 
    DP speed;
    FcellPtr next;
};

class ECM{
    public:
        ECM(int ystep,int xstep):ECM_ystep(ystep),ECM_xstep(xstep),collagen(ystep/10,xstep/10),collagen_density(ystep,xstep),fibronectin_density(ystep,xstep)
        {};
        void initiate();
        void collagen_orientation(Mat_DP& PDGF, Flist& fibroblastList);
        void output_collagen(Flist& fibroblastList);

    protected:
        int ECM_ystep;
        int ECM_xstep;
        Mat_DP collagen;
        Mat_DP collagen_density;
        Mat_DP fibronectin_density;
        
};
