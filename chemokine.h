#ifndef _CHEMOKINE_H_
#define _CHEMOKINE_H_

#include <string>
#include <cmath>
#include <iostream>
#include "nr.h"

class Chemokine{
public:
    Chemokine(int xstep, int ystep, double dx_, double D);
    ~Chemokine();
    
    void initialize(double wound_radius);
    void initialize(double wound_a, double wound_b);
    void initialize_rectangle(double a, double b);
    void initialize_oval();
    void initialize2();
    void calculate_gradient();
    void diffusion(Mat_DP& vas, double tlength);
    void output(int out_i, string GF_name);
    void output_chemoFile();
    
    bool isGradientCalculated;
   
protected:
    int xstep, ystep;
    double dx;
    double diffusion_coeff;
    
public:
    Mat_DP *conc;
    Mat_DP *gradx;
    Mat_DP *grady;
};

#endif
