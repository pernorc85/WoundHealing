#include "util.h"
#include <cstdlib>
#include <iostream>
#include <cmath>

extern Mat_DP *fibroblast_density_ptr;
double fb_density(double y, double x){
    return (*fibroblast_density_ptr)[(int)y][(int)x];
}

bool isOnWoundEdge(double x, double y) {
    //return abs(abs(x-500.0)-500.0)<1.5 || abs(abs(y-500.0)-500.0) < 1.5;
    //return abs(sqrt(pow(x-500.0,2)+pow(y-500.0,2)) - 400.0) < 1.5;
    return abs(sqrt(pow(x-500.0,2)+pow((y-500.0)/0.8,2)) - 400.0) < 1.5;
}

void elasticity_coefficient(Mat_DP F, Vec_DP M, DP c_density, Mat4D_DP& A);

void matrix_inverse(Mat_DP& F_inverse,Mat_DP F) {
     double det_F=F[0][0]*F[1][1]-F[0][1]*F[1][0];
     if(det_F !=0){
         F_inverse[0][0] = F[1][1]/det_F;
         F_inverse[0][1] = -F[0][1]/det_F;
         F_inverse[1][0] = -F[1][0]/det_F;
         F_inverse[1][1] = F[0][0]/det_F;
     } else {
         F_inverse[0][0] = 0;
         F_inverse[0][1] = 0;
         F_inverse[1][0] = 0;
         F_inverse[1][1] = 0;
     }

     return;
}

