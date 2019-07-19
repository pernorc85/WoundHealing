#include <cstdlib>
#include <iostream>
#include "nr.h"

double fb_density(double y, double x);
bool isOnWoundEdge(double x, double y, int mXstep, int mYstep);

void matrix_inverse(Mat_DP& F_inverse,Mat_DP F);

