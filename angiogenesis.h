#ifndef _ANGIOGENESIS_H_
#define _ANGIOGENESIS_H_

#include <list>
#include "nr.h"
#include "record_vasculature_network.h"

using namespace std;

typedef struct Endocyte Ecell;
typedef Ecell* EcellPtr;
struct Endocyte{
    Endocyte(xx_, yy_, theta_, id_, branching_):
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
    Elist(int xstep,int ystep){
        cellmatrix_ptr = new Mat_DP(ystep, xstep);
        vas_density_ptr = new Mat_DP(ystep, xstep);
    }
    ~Elist() {
        delete cellmatrix_ptr;
        delete vas_density_ptr;
    }

    void buildElist(EcellPtr *sPtr, EdgePtr* edgeList);
    void Ecell_move(EcellPtr cellPtr, DP gradx_VEGF, DP grady_VEGF, DP fdensity, DP cdensity, DP collagen, Mat_DP& cellmatrix, DP speed);
    void Ecell_branch(EcellPtr* sPtr, EdgePtr* edgeList, NodePtr* nodeList, Mat_DP& collagen, Mat_DP& gradx_VEGF, Mat_DP& grady_VEGF);
    void Ecell_looping(EcellPtr* sPtr, EdgePtr* edgeList, NodePtr* nodeList, Mat_DP cellmatrix);

    void calculate_for_oxygen(Mat_DP& vas_density);
    
    void output_vasculuture(bool output);
    void output_vasculature_density();
protected:
    std::list<Endocyte> mECells;
    EcellPtr Eroot;
    Mat_DP *cellmatrix_ptr;
    Mat_DP *vas_density_ptr;
    static int branch_id;
};



#endif
