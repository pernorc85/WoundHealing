#pragma once
#ifndef _ANGIOGENESIS_H_
#define _ANGIOGENESIS_H_

#include <list>
#include "nr.h"

struct Edge;
struct Node;
typedef Edge* EdgePtr;
typedef Node* NodePtr;
class ECM;
using namespace std;

typedef struct Endocyte Ecell;
typedef Ecell* EcellPtr;
struct Endocyte{
    Endocyte(DP xx_, DP yy_, DP theta_, int id_, int branching_):
        xx(xx_), yy(yy_), theta(theta_), id(id_), branching(branching_) {}
    DP xx, yy;
    DP theta;
    int id;
    int branching;
    EdgePtr curEdgePtr;
    EcellPtr next;
};

class Elist{
public:
    Elist(int xstep,int ystep):
            mXstep(xstep), mYstep(ystep), speed(15) {
        Eroot = new EcellPtr(NULL);
        cellmatrix_ptr = new Mat_DP(ystep, xstep);
        vas_density_ptr = new Mat_DP(ystep, xstep);
    }
    ~Elist() {
        delete cellmatrix_ptr;
        delete vas_density_ptr;
    }

    void buildElist(unsigned int ENinit, EdgePtr* edgeList);
    void Elist_move(Mat_DP &gradx_VEGF, Mat_DP &grady_VEGF, ECM &extraCellularMatrix, DP time_step); 
    void Ecell_move(EcellPtr cellPtr, DP gradx_VEGF, DP grady_VEGF, DP fdensity, DP cdensity, DP collagen, DP time_step);
    void Elist_branch(EdgePtr* edgeList, NodePtr* nodeList, ECM& extraCellularMatrix, Mat_DP& gradx_VEGF, Mat_DP& grady_VEGF, DP time_step);
    void Elist_looping(EdgePtr* edgeList, NodePtr* nodeList);
    void Ecell_newtip();

    void calculate_for_oxygen(Mat_DP& vas_density);
    
    void update_cellmatrix();
    void output_cellmatrix();
    void output_vasculature_density();
protected:
    int mXstep, mYstep;
    EcellPtr *Eroot;
    Mat_DP *cellmatrix_ptr;
    Mat_DP *vas_density_ptr;
    DP speed;//default speed of Endocyte 
    static int branch_id;
};



#endif
