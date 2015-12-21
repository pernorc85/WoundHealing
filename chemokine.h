#ifndef _CHEMOKINE_H_
#define _CHEMOKINE_H_

#include <string>
#include <cmath>
#include <iostream>
#include "nr.h"

class chemokine{
public:
    chemokine(int nx_, int ny_, double dx_, double dt_);
    ~chemokine();
    
    void initialize(double wound radius);
    void calculate_gradient();
    void diffusion(double D, Mat_DP& cellmatrix, double tlength);
    void output();
    void output_chemoFile();
    
    bool isGradientCalculated;

protected:
    int xstep, ystep;
    double dx;
    double dt;
    Mat_DP *conc;
    
public:
    Mat_DP *gradx;
    Mat_DP *grady;
};

#endif
