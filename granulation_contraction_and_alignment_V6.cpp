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
#include "ECM_Flist.h"
#include "chemokine.h"
#include "nr.h"

#define FALSE 0
#define TRUE 1 
#define FNinit 600*2 
#define MAXSTEP 100

using namespace std; 

int tstep;        
int xstep = 1000, ystep = 1000;//unit um

Mat_DP tensorK11(ystep,xstep);
Mat_DP tensorK12(ystep,xstep);
Mat_DP tensorK22(ystep,xstep);
Mat_DP speed(ystep,xstep);
Mat_DP theta(ystep,xstep);
 
 
int ADI(DP tlength, DP D, Mat_DP& PDGF, bool output); 
void Calculate_chemokine_gradient(Mat_DP& PDGF, Mat_DP& gradx_cytokine, Mat_DP& grady_cytokine);
void output_fibroblast(FcellPtr* sPtr);
void output_woundcontour();

 
int main()
{
    int i,j,k;
    
//****************setup the continuum part of the model*************************        
    DP xsize, ysize, tsize, theta, D_P;
    
    cout <<"Please enter the diffusion coefficient for PDGF."<<endl;
    cin >> D_P; 
    cout <<"Please enter the end of time."<<endl;
    cin >> tsize; 
    cout <<"xstep="<<xstep<<"ystep="<<ystep<<endl;
    
    srand(time(NULL)); 
    tstep = (int)(tsize/tlength)+1;

    Chemokine PDGF(ystep,xstep,1.0,tstep);
    PDGF.initialize(350.0);
    Mat_DP source(ystep,xstep);
    
//****************setup the discrete part of the model**************************
    Flist fibroblastList(FNinit);
    
    DP speed = 15;
    cout <<" ok for build fibroblasts list."<<endl;
    
    
    ECM extraCellularMatrix(ystep,xstep);
    extraCellularMatrix.initiate();
    
    //deallog.depth_console (0); 
    //ElasticProblem<2> elastic_problem_2d;

    for(k=0; k<tstep; k++){
//**********PDE for cytokine reaction-diffusion*********************************  
        if(k % 10 == 0)ADI(tlength, D_P, PDGF, TRUE);
        else ADI(tlength, D_P, PDGF, FALSE);  
        cout << "ADIout" << endl;
        
//         //starts arbitrarily large number
//        first_iteration_tag = 1;
//
//        residue_total = 100000;
//        int iter_counter = 0;
//        int continue_tag;
//        
//        if(k % 10 == 0){
//            do{ //condition of convergence???????????????????????????????????????????????????????
//                cout <<endl << "iteration:"<<first_iteration_tag<<endl;
//                elastic_problem_2d.run ();
//                if(first_iteration_tag < 3){first_iteration_tag++; residue_total = 100000;}
//                //first_iteration_tag = 1 : don't have solution
//                //first_iteration_tag = 2 : have solution; don't have solution_pre
//                else if(first_iteration_tag == 3)first_iteration_tag = 3;//have solution and solution_pre 
//                iter_counter ++;
//                    
//                continue_tag = 1;
//                if(iter_counter > 20 && iter_counter%10 == 0){
//                    printf("Continues?");
//                    scanf("%d", &continue_tag);
//                }
//                if(continue_tag == 0)break;     
//            }while(residue_total > 500);
//            output_woundcontour();
//            printf("output wound contour...\n");
//        }
        
 
//**********collagen_alignment**************************************************       
        extraCellularMatrix.collagen_orientation( *(PDGF.conc), fibroblastList);
        cout << "collagen alignment done" << endl; 
        
//**********interpolate for fibroblast******************************************                    
        PDGF.calculate_gradient();                 
        fibroblastList.Flist_move(*(PDGF.gradx), *(PDGF.grady), extraCellularMatrix);         
        cout << "fibroblast move done" << endl;
        
//*************output the orientation of collagen and fibroblasts***************
        if(k % 10 == 0){
            extraCellularMatrix.output_collagen(fibroblastList);
            printf("output collagen...\n");
        } 
//*************output the density of fibroblasts********************************               
//       if(k % 10 == 0){    
//            output_fibroblast(&Flist);
//            printf("output fibroblast...\n");
//        }                
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