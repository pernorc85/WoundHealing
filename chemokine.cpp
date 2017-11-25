#include "chemokine.h"
#include <cstdlib>
#include <cmath>
#include "nr.h"
#include "BMP.h"

using namespace std;
#define F(x)   0.5 - 1.0*x

Chemokine::Chemokine(int nx_, int ny_, double dx_, double dt_):
    xstep(nx_),ystep(ny_),dx(dx_),dt(dt_){
    conc = new Mat_DP(ystep,xstep);
    gradx = NULL;
    grady = NULL;
}

Chemokine::~Chemokine(){
    delete conc;
    if(gradx != NULL)delete gradx;
    if(grady != NULL)delete grady;
}

void Chemokine::initialize(double a/*major semi_axis*/, double b/*minor semi_axis*/){//if radius = 350
    int i,j;
    for(i=0;i<ystep;i++){//from bottom up
        for(j=0;j<xstep;j++){//from left to right
            if(pow(j-xstep/2,2)+pow((i-ystep/2)/0.8,2) >= (a + 100)*(a + 100) )(*conc)[i][j] = 0;//450^2 
            else if(pow(j-xstep/2,2)+pow((i-ystep/2)/0.8,2) < (a+100)*(a+100) && 
                    pow(j-xstep/2,2)+pow((i-ystep/2)/0.8,2) > (a-100)*(a-100) ){//450^2--250^2     
                (*conc)[i][j] = -sin( (sqrt(pow(j-xstep/2,2)+pow((i-ystep/2)/0.8,2))-a)/100*M_PI/2 )*0.5 + 0.5;
            } 
            else (*conc)[i][j] = 1;  
        }
    }
    isGradientCalculated = false;
}

void Chemokine::initialize(double radius){//if radius = 350
    int i,j;
    for(i=0;i<ystep;i++){//from bottom up
        for(j=0;j<xstep;j++){//from left to right
            double dist_from_center = sqrt(pow(j-xstep/2,2) + pow((i-ystep/2)/0.8, 2));
            if(dist_from_center >= radius + 100 ){
		(*conc)[i][j] = 0;//450^2 
            } else if(dist_from_center < radius+100 && dist_from_center > radius-100 ){     
                (*conc)[i][j] = -sin( (dist_from_center-radius)/100*M_PI/2 )*0.5 + 0.5;
            } else {
		(*conc)[i][j] = 1;  
            }
        }
    }
    isGradientCalculated = false;
}

void Chemokine::calculate_gradient(){
    if(gradx == NULL)gradx = new Mat_DP(ystep, xstep);
    if(grady == NULL)grady = new Mat_DP(ystep, xstep);
    for(int i=0;i<ystep;i++){
        for(int j=0;j<xstep;j++){
            if(j>=3 && j<xstep-3) (*gradx)[i][j] = (*conc)[i][j+3]-(*conc)[i][j-3];
            else (*gradx)[i][j] = 0.0;
            if(i>=3 && i<ystep-3) (*grady)[i][j] = (*conc)[i+3][j]-(*conc)[i-3][j];
            else (*grady)[i][j] = 0.0;
        }
    }
    isGradientCalculated = true;
}

void Chemokine::diffusion(double D, double tlength){
     //...call external function ADI...
}

void Chemokine::output(){
     printf("output chemokine profile...\n");
     static int out_i = 0;
     int jj;
     char file_name[15] = "chem00000.BMP";
     file_name[8] = 48+out_i%10;
     file_name[7] = 48+out_i/10%10;
     file_name[6] = 48+out_i/100%10;
     file_name[5] = 48+out_i/1000%10;
     BMP::Prepare_BMP_Format(ystep,xstep);
     BMP::output_BMP(file_name, 15, (*conc), xstep, ystep);
     out_i++;
}

void Chemokine::output_chemoFile(){
     FILE* dataFile;
     printf("output chemokine to text...\n");
     dataFile = fopen("chemokineBackup.txt","w");
     for(int i=0;i<ystep;i++){
         for(int j=0;j<xstep;j++){
             fprintf(dataFile, "%f\t", (*conc)[i][j]);
         }
         fprintf(dataFile,"\n");
     }
     fclose(dataFile);
}
     
            
    
