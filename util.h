#include <cstdlib>
#include <iostream>
#include "nr.h"

double fb_density(double y, double x);
bool isOnWoundEdge(double x, double y, int mXstep, int mYstep);

void matrix_inverse(Mat_DP& F_inverse,Mat_DP F);
Mat_DP transpose(Mat_DP F);
Mat_DP multiply(Mat_DP A, Mat_DP B);
Vec_DP eigen_decomposition(Mat_DP G, Mat_DP &eigen_vecs);
