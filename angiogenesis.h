#ifndef _ANGIOGENESIS_H_
#define _ANGIOGENESIS_H_

#include <list>
#include "nr.h"
#include "record_vasculature_network.h"

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
            mXstep(xstep), mYstep(ystep) {
        cellmatrix_ptr = new Mat_DP(ystep, xstep);
        vas_density_ptr = new Mat_DP(ystep, xstep);
    }
    ~Elist() {
        delete cellmatrix_ptr;
        delete vas_density_ptr;
    }

    void buildElist(EdgePtr* edgeList);
    void Ecell_move(EcellPtr cellPtr, DP gradx_VEGF, DP grady_VEGF, DP fdensity, DP cdensity, DP collagen, DP time_step);
    void Ecell_branch(EcellPtr* sPtr, EdgePtr* edgeList, NodePtr* nodeList, Mat_DP& collagen, Mat_DP& gradx_VEGF, Mat_DP& grady_VEGF, DP time_step);
    void Ecell_looping(EcellPtr* sPtr, EdgePtr* edgeList, NodePtr* nodeList);
    void Ecell_newtip();

    void calculate_for_oxygen(Mat_DP& vas_density);
    
    void update_cellmatrix(bool output);
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
