//#pragma once
#ifndef _ECM_FLIST_BMP_H_
#define _ECM_FLIST_BMP_H_

#include <string>
#include <cmath>
#include <iostream> 
#include "nr.h"
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

class Flist;

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
        
    friend class Flist;
};




class Flist
{
    public:
        Flist(); 
        Flist(int FNinit);
        
        void buildFlist(FcellPtr* sPtr); 
        void Flist_move(Mat_DP& gradx_cytokine, Mat_DP& grady_cytokine, ECM& extraCellularMatrix); 
        void Fcell_move(DP gradx, DP grady, DP fdensity, DP cdensity, Mat_DP& collagen); 
        void Fcell_getinfo(Mat_DP& PDGF, Mat_DP& gradx_cytokine, Mat_DP& grady_cytokine);
        void output_fibroblast(FcellPtr* sPtr);

    protected:
        FcellPtr Froot;
        FcellPtr curPtr;
        int FN;
        DP speed;
    //ÕâÀï²»ÄÜÓÃFcellPtr *sPtr; *sPtr = NULL;ÕâÑùÐ´Ö®ºósystem("PAUSE");È«²¿Ê§Ð§
    //¶øÇÒbuildFlistÒ²Ê§Ð§ 
    //FcellPtr curPtr;     
    //curPtr = Flist;
    friend class ECM;
};


#endif
