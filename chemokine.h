#ifndef _CHEMOKINE_H_
#define _CHEMOKINE_H_

#include <string>
#include <cmath>
#include <iostream>
#include "nr.h"

class Chemokine{
public:
    Chemokine(int nx_, int ny_, double dx_, double dt_);
    ~Chemokine();
    
    void initialize(double wound_radius);
    void initialize(double wound_a, double wound_b);
    void calculate_gradient();
    void diffusion(double D, double tlength);
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
