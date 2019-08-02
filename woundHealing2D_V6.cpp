//granulation-contraction and alignment_V6.cpp
//    VERSION: 6 (1 November 2015)
//    PURPOSE: C++ program for wound contraction model.
//             The main function initializes all parts of the model and carries out
//             simulation in a time interval defined by user.
//             Include functions for fibroblast migation and collagen alignment(align with both fibroblasts
//             and tension lines).
//    INPUT FILES:
//             None.
//    OUTPUT FILES:
//             None.
#define FALSE 0
#define TRUE 1 
//#define ANGIOGENESIS 1
#define MECHANICS 1
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include "BMP.h"
#include "ECM.h"
#include "Flist.h"
#include "chemokine.h"
#include "angiogenesis.h"
#ifdef MECHANICS
#include "step8_wound_mechanics.h"
#include "step8_wound_mechanics.cc"
#endif
#include "nr.h"

using namespace std; 

const int FDensity = 600;        
const int ENinit = 12;
const int xstep = 2000, ystep = 2000;//unit um
const int wound_radius = 350;
const int wound_long_axis = 400;
const int wound_short_axis = 300;
const DP time_total = 200;//unit hr
const DP time_step = 0.15;//unit hr
const DP D_P = 0.5;
const DP decay_P = 0.1;
const DP D_V = 1.0;
const DP decay_V = 0.1;
 
Mat_DP *fibroblast_density_ptr;
Mat_DP *speedfield_ptr;
Mat_DP *thetafield_ptr;

 
int main()
{   
    fibroblast_density_ptr = new Mat_DP(ystep, xstep);
    speedfield_ptr = new Mat_DP(ystep, xstep);
    thetafield_ptr = new Mat_DP(ystep, xstep);

    cout << "**********************************************************************" << endl;
    cout << "****************setup the continuum part of the model*****************" << endl;
    cout << "**********************************************************************" << endl;
    srand(time(NULL)); 
    int nSteps = (int)(time_total/time_step)+1;

    Chemokine PDGF(xstep,ystep,1.0,D_P,decay_P);
    //PDGF.initialize(wound_radius);
    PDGF.initialize_rectangle(800, 600);
  
#ifdef ANGIOGENESIS 
    Chemokine VEGF(xstep,ystep,1.0,D_V,decay_V);
    VEGF.initialize(wound_radius); 
    //VEGF.initialize(400, 300); 
#endif
//****************setup the discrete part of the model**************************
    cout << "**********************************************************************" << endl;
    cout << "****************setup the discrete part of the model*****************" << endl;
    cout << "**********************************************************************" << endl;
    Flist fibroblastList(xstep, ystep);
    //fibroblastList.initialize(FDensity, wound_radius);
    fibroblastList.initialize_rectangle(FDensity, 800, 600);
            
    ECM extraCellularMatrix(ystep,xstep);
    //extraCellularMatrix.initialize(wound_radius);
    extraCellularMatrix.initialize_rectangle(800, 600);
    
#ifdef ANGIOGENESIS 
    EdgePtr edgeList = NULL;
    NodePtr nodeList = NULL;
    Elist endoList(xstep, ystep);
    endoList.buildElist(ENinit, wound_radius, &edgeList);
#endif

#ifdef MECHANICS 
    //deallog.depth_console (0); 
    ElasticProblem<2> elastic_problem_2d(xstep, ystep);
#endif 
    //****************main loop*****************************************************

    for(size_t k=0; k < nSteps; k++){
        cout << "========================k=" << k << "=======================" << endl;
        //**********PDE for cytokine reaction-diffusion*********************************  
        PDGF.diffusion(time_step);
 
#ifdef MECHANICS         
        if(k % 10 == 0){
            fibroblastList.output_fibroblast();
            elastic_problem_2d.run_until_converge(extraCellularMatrix.collagen,
                                       extraCellularMatrix.collagen_density,
                                       *fibroblast_density_ptr);
        }
#endif        
 
        //**********collagen_alignment**************************************************       
#ifdef MECHANICS
        extraCellularMatrix.collagen_orientation(fibroblastList, elastic_problem_2d.tissue_displacement_x, 
                                                     elastic_problem_2d.tissue_displacement_y, time_step);
#else
        extraCellularMatrix.collagen_orientation_with_fibroblast(fibroblastList, time_step);
#endif
        cout << "collagen alignment done" << endl; 
        
        //**********interpolate for fibroblast******************************************                    
        PDGF.calculate_gradient();          
        cout << "PDGF gradient calculated" << endl;      
        fibroblastList.Flist_move(PDGF, extraCellularMatrix, time_step);         
        cout << "fibroblast move done" << endl;
    
        //******************************************************************************
#ifdef ANGIOGENESIS
        VEGF.diffusion(time_step);
        VEGF.calculate_gradient();
        endoList.Elist_move(*(VEGF.gradx), *(VEGF.grady), extraCellularMatrix, time_step);    
        cout << "endocyte move done" << endl;
        endoList.Elist_branch(&edgeList, &nodeList, extraCellularMatrix, *(VEGF.gradx), *(VEGF.grady), time_step);
        endoList.Elist_looping(&edgeList, &nodeList);
        endoList.update_cellmatrix();
#endif 
        //*************output the orientation of collagen and fibroblasts***************
        if(k % 10 == 0){
            extraCellularMatrix.output_collagen(fibroblastList);
            printf("output collagen...\n");
#ifdef ANGIOGENESIS
            endoList.output_cellmatrix();
            VEGF.output(k/10, "VEGF");
#endif
        }
        if(k % 50 == 0){
            PDGF.output(k/10, "PDGF");
        }
    }   
    return 0;    
}
 

