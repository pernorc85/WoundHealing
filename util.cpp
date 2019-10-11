#include "util.h"
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <cassert>

extern Mat_DP *fibroblast_density_ptr;
double fb_density(double y, double x){
    return (*fibroblast_density_ptr)[(int)y][(int)x];
}

bool isOnWoundEdge(double x, double y, int mXstep, int mYstep) {
    double wound_length = 800, wound_width = 600;
    return (abs(abs(x-mXstep/2)- wound_length/2)<1.5 and abs(y-mYstep/2) < wound_width/2) 
           || (abs(abs(y-mYstep/2)-wound_width/2) < 1.5 and abs(x-mXstep/2) < wound_length/2);
    //return abs(sqrt(pow(x-mXstep/2,2)+pow(y-mYstep/2,2)) - 350.0) < 2.0;
    //double a = 400, b = 300;
    //return abs(sqrt(pow(x-mXstep/2,2)+pow(a/b*(y-mYstep/2),2)) - a) < 2.0;
}

void elasticity_coefficient(Mat_DP F, Vec_DP M, DP c_density, Mat4D_DP& A);

double get_theta(double x, double y){
    double theta = -1;
    if(x != 0){
        if(x > 0 && y >= 0)theta = atan(y/x);//1st phase
        else if(x < 0 && y >= 0)theta = M_PI + atan(y/x);//2nd phase
        else if(x < 0 && y < 0)theta = M_PI + atan(y/x);//3nd phase
        else if(x > 0 && y < 0)theta = 2*M_PI + atan(y/x);//4th phase
    } else{
         if(y > 0)theta =  M_PI/2;
         else if(y < 0)theta = M_PI*3/2;
         else theta = -1;
    }
    return theta;
}

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

Mat_DP transpose(Mat_DP F){
    assert(F.nrows() == F.ncols());
    Mat_DP Ft(F.nrows(),F.nrows());
    for(size_t i = 0; i < F.nrows(); i++) {
        for(size_t j = 0; j < F.ncols(); j++) {
            Ft[i][j] = F[j][i];
        }
    }
    return Ft;
}

Mat_DP multiply(Mat_DP A, Mat_DP B){
    assert(A.ncols() == B.nrows());
    Mat_DP C(A.nrows(), B.ncols());
    for(size_t i = 0; i< A.nrows(); i++) {
        for(size_t j = 0; j < B.ncols(); j++) {
            C[i][j] = 0.0;
            for(size_t k = 0; k < A.ncols(); k++) {
                C[i][j] += A[i][k]*B[k][j];
            }
        }
    }
    return C;
}

Vec_DP eigen_decomposition(Mat_DP G, Mat_DP &eigen_vecs){
    Vec_DP eig(2);
    double det = sqrt(4*G[0][1]*G[0][1] + pow(G[0][0]-G[1][1], 2));
    eig[0] = 0.5 * (G[0][0] + G[1][1] + det);
    eig[1] = 0.5 * (G[0][0] + G[1][1] - det);

    double b = G[0][1];
    if (b) {
        double v0 = eig[0] - G[1][1];
        double v1 = b;
        double vn = sqrt(v0*v0 + v1*v1);
        eigen_vecs[0][0] = v0/vn;
        eigen_vecs[1][0] = v1/vn;

        v0 = eig[1] - G[1][1];
        vn = sqrt(v0*v0 + v1*v1);
        eigen_vecs[0][1] = v0/vn;
        eigen_vecs[1][1] = v1/vn;
    } else if (G[0][0] >= G[1][1]) {
        eigen_vecs[0][0] = 1;
        eigen_vecs[0][1] = 0;
        eigen_vecs[1][0] = 0;
        eigen_vecs[1][1] = 1;
    } else {
        eigen_vecs[0][0] = 0;
        eigen_vecs[0][1] = 1;
        eigen_vecs[1][0] = 1;
        eigen_vecs[1][1] = 0;
    } 
    return eig;
}

