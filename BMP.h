#ifndef _BMP_H_
#define _BMP_H_

#include <map>
#include "nr.h"
using namespace std; 

namespace BMP{

void Prepare_BMP_Format(FILE *fp_BMP, int width, int height); 
void output_BMP (char* file_name,int n, const Mat_DP& data, int xstep,int ystep);
void output_BMP2(char* file_name,int n, const Mat_DP& data, const Mat_DP& cellmatrix, int xstep, int ystep);
void output_BMP3(char* file_name,int n, const map<pair<int, int>, int>& data, int xstep,int ystep);

}//end of namespace BMP

#endif
