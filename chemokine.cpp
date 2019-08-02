#include "chemokine.h"
#include <cstdlib>
#include <cmath>
#include "nr.h"
#include "BMP.h"

using namespace std;
extern int ADI(DP, Mat_DP&, DP, DP);

Chemokine::Chemokine(int nx_, int ny_, double dx_, double D, double decay_rate_):
    xstep(nx_),ystep(ny_),dx(dx_),diffusion_coeff(D),decay_rate(decay_rate_){
    conc = new Mat_DP(ystep,xstep);
    gradx = NULL;
    grady = NULL;
}

Chemokine::~Chemokine(){
    delete conc;
    if(gradx != NULL)delete gradx;
    if(grady != NULL)delete grady;
}

void Chemokine::initialize(double a/*major semi_axis*/, double b/*minor semi_axis*/){
    int i,j;
    double c = sqrt(a*a - b*b);
    for(i=0;i<ystep;i++){//from bottom up
        for(j=0;j<xstep;j++){//from left to right
            double dist_from_focus = sqrt(pow(j-xstep/2,2) + pow(a/b*(i-ystep/2), 2));
            if(dist_from_focus > a + 50 ){
                (*conc)[i][j] = 0; 
            } else if(dist_from_focus < a + 50 && dist_from_focus > a - 50 ){     
                (*conc)[i][j] = -sin( (dist_from_focus - a)/50*M_PI/2 )*0.5 + 0.5;
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
            double dist_from_center = sqrt(pow(j-xstep/2,2) + pow(i-ystep/2, 2));
            if(dist_from_center >= radius + 25 ){
		(*conc)[i][j] = 0; 
            } else if(dist_from_center < radius+25 && dist_from_center > radius-75 ){     
                (*conc)[i][j] = -sin( (dist_from_center-radius+25)/50*M_PI/2 )*0.5 + 0.5;
            } else {
		(*conc)[i][j] = 1;  
            }
        }
    }
    isGradientCalculated = false;
}

void Chemokine::initialize_rectangle(double a, double b){
    int xcenter = xstep * 0.5;
    int ycenter = ystep * 0.5;
    for(int i=0; i<ystep; i++){//from bottom up
        for(int j=0; j<xstep; j++){
            if(abs(j-xcenter)>a*0.5+50 || abs(i-ycenter)>b*0.5+50) {
                (*conc)[i][j] = 0.0;
            } else if(abs(j-xcenter)>a*0.5-50 || abs(i-ycenter)>b*0.5-50) {
                (*conc)[i][j] = 0.0;
                if (abs(j-xcenter)<=a*0.5+50 && abs(j-xcenter)>a*0.5-50) {
                    double frame = abs(j-xcenter) - a*0.5;
                    (*conc)[i][j] = -sin( frame*0.02 * M_PI/2 )*0.5 + 0.5;
                }
                if(abs(i-ycenter)<=b*0.5+50 && abs(i-ycenter)>b*0.5-50){
                    double frame = abs(i-ycenter) - b*0.5;
                    double tmp = -sin( frame*0.02 * M_PI/2 )*0.5 + 0.5;
                    if ((*conc)[i][j] != 0.0 and tmp < (*conc)[i][j]) 
                        (*conc)[i][j] = tmp;
                    else if((*conc)[i][j] == 0.0)
                        (*conc)[i][j] = tmp;
                }
            } else {
                (*conc)[i][j] = 1.0;
            }
        }
    }
}

void Chemokine::initialize2(){
    int i, j;
    for(i=0;i<xstep;i++)
        for(j=0;j<ystep;j++)(*conc)[i][j]=0;
    for(i=xstep/2-490;i<xstep/2+490;i++)
        for(j=ystep/2-490;j<ystep/2+490;j++)(*conc)[i][j] = 1;
    for(i=xstep/2-480;i<xstep/2+480;i++)
        for(j=ystep/2-480;j<ystep/2+480;j++)(*conc)[i][j] = 1.3;
    for(i=xstep/2-460;i<xstep/2+460;i++)
        for(j=ystep/2-460;j<xstep/2+460;j++)(*conc)[i][j] = 1.5;
    for(i=xstep/2-430;i<xstep/2+430;i++)
        for(j=ystep/2-430;j<ystep/2+430;j++)(*conc)[i][j] = 1.7;
    for(i=xstep/2-400;i<xstep/2+400;i++)
        for(j=ystep/2-400;j<ystep/2+400;j++)(*conc)[i][j] = 2;
    for(i=xstep/2-350;i<xstep/2+350;i++)
        for(j=ystep/2-350;j<ystep/2+350;j++)(*conc)[i][j] = 2.3;
    for(i=xstep/2-300;i<xstep/2+300;i++)
        for(j=ystep/2-300;j<ystep/2+300;j++)(*conc)[i][j] = 2.5;
    for(i=xstep/2-200;i<xstep/2+200;i++)
        for(j=ystep/2-200;j<xstep/2+200;j++)(*conc)[i][j] = 3;
    for(i=xstep/2-100;i<xstep/2+100;i++)
        for(j=ystep/2-100;j<ystep/2+100;j++)(*conc)[i][j] = 4;
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

void Chemokine::diffusion(double time_step){
     ADI(time_step, *conc, diffusion_coeff, decay_rate); 
}

void Chemokine::output(int out_i, std::string GF_name){
     printf("output chemokine profile...\n");
     char file_name[21];
     sprintf(file_name, "output/%s%05d.BMP", GF_name.c_str(), out_i);
     BMP::output_BMP(file_name, 15, (*conc), xstep, ystep);
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
     
            
    
