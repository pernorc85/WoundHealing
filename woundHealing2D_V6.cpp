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
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>
//#include "step8_wound_mechanics.h"
#include "ECM.h"
#include "Flist.h"
#include "chemokine.h"
#include "angiogenesis.h"
#include "nr.h"

#define FALSE 0
#define TRUE 1 
#define ANGIOGENESIS 1

using namespace std; 

const int FNinit = 1000;        
const int ENinit = 12;
const int xstep = 1000, ystep = 1000;//unit um
const int wound_radius = 350;
const DP time_total = 200;//unit hr
const DP time_step = 0.15;//unit hr
const DP D_P = 0.5;
const DP decay_P = 0.1;
const DP D_V = 1.0;
const DP decay_V = 0.1;
 
Mat_DP *fibroblast_density_ptr;
Mat_DP *speedfield_ptr;
Mat_DP *thetafield_ptr;

void output_woundcontour();

 
int main()
{   
    fibroblast_density_ptr = new Mat_DP(ystep, xstep);
    speedfield_ptr = new Mat_DP(ystep, xstep);
    thetafield_ptr = new Mat_DP(ystep, xstep);
//****************setup the continuum part of the model*************************
    cout << "initializing continuum part..." << endl;        
    srand(time(NULL)); 
    int nSteps = (int)(time_total/time_step)+1;

    Chemokine PDGF(ystep,xstep,1.0,D_P,decay_P);
    PDGF.initialize(wound_radius);
  
#ifdef ANGIOGENESIS 
    Chemokine VEGF(ystep,xstep,1.0,D_V,decay_V);
    VEGF.initialize(wound_radius); 
#endif
//****************setup the discrete part of the model**************************
    cout << "initializing discrete part..." << endl;
    Flist fibroblastList(xstep, ystep, FNinit);
            
    ECM extraCellularMatrix(ystep,xstep);
    extraCellularMatrix.initiate(wound_radius);

#ifdef ANGIOGENESIS 
    EdgePtr edgeList = NULL;
    NodePtr nodeList = NULL;
    Elist endoList(xstep, ystep);
    endoList.buildElist(ENinit, wound_radius, &edgeList);
#endif
 
    //deallog.depth_console (0); 
    //ElasticProblem<2> elastic_problem_2d;
    //
//****************main loop*****************************************************

    for(int k=0; k < nSteps; k++){
        cout << "========================k=" << k << "=======================" << endl;
        //**********PDE for cytokine reaction-diffusion*********************************  
        PDGF.diffusion(time_step);
        VEGF.diffusion(time_step);
 
        /* 
        if(k % 10 == 0){
            //starts arbitrarily large number
            first_iteration_tag = 1;
            residue_total = 100000;
            int iter_counter = 0;
            int continue_tag;
            do{ //condition of convergence???????????????????????????????????????????????????????
                cout <<endl << "iteration:"<<first_iteration_tag<<endl;
                elastic_problem_2d.run ();
                if(first_iteration_tag < 3){first_iteration_tag++; residue_total = 100000;}
                //first_iteration_tag = 1 : don't have solution
                //first_iteration_tag = 2 : have solution; don't have solution_pre
                else if(first_iteration_tag == 3)first_iteration_tag = 3;//have solution and solution_pre 
                iter_counter ++;
                    
                continue_tag = 1;
                if(iter_counter > 20 && iter_counter%10 == 0){
                    printf("Continues?");
                    scanf("%d", &continue_tag);
                }
                if(continue_tag == 0)break;     
            }while(residue_total > 500);
            output_woundcontour();
            printf("output wound contour...\n");
        }
        */   
 
        //**********collagen_alignment**************************************************       
        extraCellularMatrix.collagen_orientation(fibroblastList, time_step);
        cout << "collagen alignment done" << endl; 
        
        //**********interpolate for fibroblast******************************************                    
        PDGF.calculate_gradient();          
        cout << "PDGF gradient calculated" << endl;      
        fibroblastList.Flist_move(PDGF, extraCellularMatrix, time_step);         
        cout << "fibroblast move done" << endl;
    
        //******************************************************************************
#ifdef ANGIOGENESIS
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
 

/*
void output_woundcontour(){
    Mat_DP woundcontour(ystep,xstep);
    DP tdx,tdy;
    int i,j;
    for(i=0;i<ystep;i++){
        for(j=0;j<xstep;j++){
            woundcontour[i][j] = 0;
        }
    }
    
    for(i=0;i<ystep;i++){
        for(j=0;j<xstep;j++){
            //if(abs(abs(j-500.0)-500.0)<1.5 || abs(abs(i-500.0)-500.0)<1.5){ 
            if(abs(sqrt(pow(j-500.0,2)+pow((i-500.0)/0.8,2)) - 400.0)<1.5){
                tdx = tissue_displacement_x[i][j];
                tdy = tissue_displacement_y[i][j];
                if(i+(int)tdy>=0 && i+(int)tdy<=999 && j+(int)tdx>=0 && j+(int)tdy<=999){ 
                    woundcontour[i+(int)tdy][j+(int)tdx] = 1;
                }
            }
            //if( abs(abs(i-500.0)-300)<1.5 && abs(j-500.0) < 450){ 
//            if(abs(sqrt(pow(i-500.0,2)+pow(j-500.0,2)) - 400.0)<1.5){
//                tdx = tissue_displacement_x[i][j];
//                tdy = tissue_displacement_y[i][j];
//                woundcontour[i+(int)tdy][j+(int)tdx] = 1;
//            }

        }
    }
        
    char  file_name[14]="cont00000.BMP";
    static int out_i = 0;
    file_name[8] = 48+out_i%10;        file_name[7] = 48+out_i/10%10;
    file_name[6] = 48+out_i/10/10%10;  file_name[5] = 48+out_i/10/10/10%10;
    file_name[4] = 48+out_i/10/10/10/10%10;
    fp_BMP = fopen(file_name,"w");
    Prepare_BMP_Format(xstep,ystep);
    int jj;
    for(i=0;i<ystep;i++) {
        for(j=0;j<xstep;j++) {
            jj = woundcontour[i][j];
	        if(jj == 0 )  fwrite(&b_0,  1,3,fp_BMP);  
	        if(jj == 1)  fwrite(&b_16, 1,3,fp_BMP);             
        }
    }
    fclose(fp_BMP); 
    out_i++;
}
*/

