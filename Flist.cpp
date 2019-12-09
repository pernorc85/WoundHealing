#include "Flist.h"
#include <string>
#include <cmath>
#include <iostream>
#include "chemokine.h"
#include "ECM.h"
#include "BMP.h"
#include "util.h"

using namespace std;   

extern Mat_DP *fibroblast_density_ptr;
extern Mat_DP *speedfield_ptr;
extern Mat_DP *thetafield_ptr;

Flist::Flist(int xstep, int ystep):
        mXstep(xstep), mYstep(ystep), speed(15) {}

void Flist::initialize(int FDensity, int wound_radius){
    DP x,y,theta;
    int FNumber = FDensity * (mXstep/1000.0) * (mYstep/1000.0);
    int fid = 0;
    for(int i=0;i<FNumber;i++){               
        x = rand()%mXstep;
        y = rand()%mYstep;
        theta = (rand()%360)*M_PI/180.0;
        

        if(pow(x-mXstep/2,2)+pow(y-mYstep/2,2) >= wound_radius*wound_radius){                
    	    mFCellMap.emplace(fid, Fibroblast(x,y, theta));
            fid++;    
        }       
    }         
    return;
}

void Flist::initialize(int FDensity, double a, double b){
    DP x,y,theta;
    int FNumber = FDensity * (mXstep/1000.0) * (mYstep/1000.0);
    int fid = 0;
    for(int i=0;i<FNumber;i++){               
        x = rand()%mXstep;
        y = rand()%mYstep;
        theta = (rand()%360)*M_PI/180.0;
        
        //if(pow(x-mXstep/2,2)+pow(y-mYstep,2) > 160000 || pow(x-mXstep,2)+pow(y-mYstep,2) < 90000 ){ 
        if(pow(x-mXstep/2,2)+pow(a/b*(y-mYstep/2),2) >= a*a){
    	    mFCellMap.emplace(fid, Fibroblast(x,y, theta));    
            fid++;
        }       
    }          
    return;
}

void Flist::initialize_rectangle(int FDensity, double a, double b){
    DP x,y,theta;
    int FNumber = FDensity * (mXstep/1000.0) * (mYstep/1000.0);
    int fid = 0;
    for(int i=0; i<FNumber; i++){
        x = rand()%mXstep;
        y = rand()%mYstep;
        theta = (rand()%360)*M_PI/180.0;

        if(abs(x-mXstep/2)>a*0.5 or abs(y-mYstep/2)>b*0.5){
            mFCellMap.emplace(fid, Fibroblast(x,y, theta));
            fid++;
        }
    }
    return;
}

void Flist::initialize_oval(int FDensity){
    DP x,y,theta;
    int FNumber = FDensity * (mXstep/1000.0) * (mYstep/1000.0);
    int fid = 0;
    for(int i=0; i<FNumber; i++){
        x = rand()%mXstep;
        y = rand()%mYstep;
        theta = (rand()%360)*M_PI/180.0;
        if(pow(x-mXstep/2,2)+pow(y-mYstep/2,2) >= 122500 || pow(x-mXstep/2,2)+pow((y-mYstep/2)/0.3,2) <=122500){  
            mFCellMap.emplace(fid, Fibroblast(x,y, theta));
            fid++;
        }
    }
    return;
}
                                         
void Flist::Flist_move(Chemokine& PDGF, ECM& extraCellularMatrix, 
                               const Mat_DP& tissue_displacement_x,
                               const Mat_DP& tissue_displacement_y, DP time_step) {
    DP gradx, grady;
    DP fdensity, cdensity;
    for(auto &item : mFCellMap) {
        gradx = (*PDGF.gradx)[(int)item.second.yy][(int)item.second.xx];
        grady = (*PDGF.grady)[(int)item.second.yy][(int)item.second.xx]; 
        fdensity = extraCellularMatrix.fibronectin_density[(int)(item.second.yy/5)][(int)(item.second.xx/5)];
        cdensity = extraCellularMatrix.collagen_density[(int)(item.second.yy/5)][(int)(item.second.xx/5)];
        DP conc = (*PDGF.conc)[(int)item.second.yy][(int)item.second.xx];
        if (conc != 0.0) cout << "conc = " << conc << endl;
        Fcell_activity_dynamics(&(item.second), conc, time_step);       

        Mat_DP F(2,2);
        F[0][0] = 1.0;
        F[0][1] = 0.0;
        F[1][0] = 0.0;
        F[1][1] = 1.0; 
        if (item.second.activity > 1e-3) {
            Fcell_move(&(item.second), gradx, grady, fdensity, cdensity, extraCellularMatrix.collagen, F, time_step); 
        }
    }
}

void Flist::Fcell_activity_dynamics(Fibroblast* curPtr, DP PDGF_conc, DP time_step) {
    //a fibroblast's activity will be stimulated by high level of PDGF, 
    //and will decay slowly over time
    //therefore when PDGF is diminished, fibroblast will remain active for a period of time
    double k1 = 0.5;
    double top_activity = 1.0;
    double k2 = 0.02;//half decay time = ln2 / k2 = 34.5hr;
    curPtr->activity += time_step * (
                            k1 * PDGF_conc * (1 - curPtr->activity / top_activity) 
                            - k2 * curPtr->activity );
}


void Flist::Fcell_move(Fibroblast* curPtr, DP gradx, DP grady, DP fdensity, DP cdensity, Mat_DP& collagen, Mat_DP& F, DP time_step) {
    DP rou1 = 0.3, rou2 = 0.4;
    //************************chemotaxis********************************************
    DP gradtheta = get_theta(gradx, grady);
    DP speed_taken = 0.0;
    
    if(cdensity < 0.7)speed_taken = speed * (1 + 2*fdensity)/3;
    else speed_taken = speed * (1 + 2*fdensity)/3 * (1.7 - cdensity);
    curPtr->speed = speed_taken * curPtr->activity;
    
    DP grad=min(1.0,1000*sqrt(pow(gradx,2)+pow(grady,2)));
    rou1=rou1*grad;
    rou2=rou2*cdensity;
    //cout<<"rou2="<<rou2<<endl;
    
    //***********************guidance of collagen***********************************      
    int gridy,gridx;
    DP l_y,l_x;
   
    int tempy = (int)(curPtr->yy*10);
    if( tempy%100 <= (100 - tempy%100) )gridy=(tempy-tempy%100)/100;
    else gridy=(tempy+100-tempy%100)/100;
    int tempx = (int)(curPtr->xx*10);
    if( tempx%100 <= (100 - tempx%100) )gridx=(tempx-tempx%100)/100;
    else gridx=(tempx+100-tempx%100)/100;
    
    DP collagen_temp,collagen_taken;
    DP collagen_sum = 0,l_sum = 0,collagen_avg=0;
    
    for(int j=max(gridy-2,0);j<=min(gridy+2,mYstep/10-1);j++){
        for(int l=max(gridx-2,0);l<min(gridx+2,mXstep/10-1);l++){
            collagen_taken = collagen[j][l];
            //be maticulous here-----------------------------------------------
            while(abs(collagen_taken - curPtr->theta) > M_PI/2) {
                if(collagen_taken > curPtr->theta) 
                    collagen_taken -= M_PI;
                else
                    collagen_taken += M_PI;
            }
            /*
            DP collagen_taken2; 
            collagen_temp = collagen[j][l];
            if(abs(curPtr->theta - collagen_temp) <= M_PI/2)//case 1&4; part of 2; part of 3; part of 6;
                collagen_taken2 = collagen_temp;
            else if(abs(curPtr->theta - collagen_temp) > M_PI/2 && abs(curPtr->theta - collagen_temp) <= M_PI){
                if(curPtr->theta <= M_PI/2 && collagen_temp > M_PI/2 && collagen_temp <= M_PI)//part of 2
                     collagen_taken2 = collagen_temp - M_PI;
                else collagen_taken2 = collagen_temp + M_PI;//part of 3; part of 6;
            }
            else if(curPtr->theta > M_PI*3/2 && collagen_temp < M_PI/2){//case 7
                if(collagen_temp+2*M_PI-curPtr->theta < M_PI/2)
                    collagen_taken2 = collagen_temp + 2*M_PI;
                else collagen_taken2 = collagen_temp + M_PI;  
            }
            else collagen_taken2 = collagen_temp + M_PI;//case 5&8 
            */
            //-----------------------------------------------------------------
            //calculate the weight---------------------------------------------
            l_y = (curPtr->yy-(j-2)*10) * (curPtr->yy-(j-1)*10) * (curPtr->yy-(j+1)*10) * (curPtr->yy-(j+2)*10)/20/10/10/20;
            l_x = (curPtr->xx-(l-2)*10) * (curPtr->xx-(l-1)*10) * (curPtr->xx-(l+1)*10) * (curPtr->xx-(l+2)*10)/20/10/10/20;
            l_sum += l_y*l_x;
            collagen_sum += collagen_taken*l_y*l_x;
        }
    }
    if(l_sum != 0)collagen_avg = collagen_sum/l_sum;

    if(gradtheta >= 0){
        if(gradtheta < M_PI/2 && curPtr->theta > 3/2*M_PI && curPtr->theta < 2*M_PI){
            curPtr->theta = (1 - rou1 - rou2)*curPtr->theta + rou2*collagen_avg + rou1*(gradtheta + 2*M_PI);
        } else if(gradtheta > 3/2*M_PI && gradtheta < 2*M_PI && curPtr->theta < M_PI/2){
            curPtr->theta = (1 - rou1 - rou2)*curPtr->theta + rou2*collagen_avg + rou1*(gradtheta - 2*M_PI);
        } else {
            curPtr->theta = (1 - rou1 - rou2)*curPtr->theta + rou2*collagen_avg + rou1*gradtheta;
        }
    }
    else if(gradtheta == -1) curPtr->theta = (1 - rou2)*curPtr->theta + rou2 * collagen_avg;

    while(curPtr->theta > 2*M_PI) {
        curPtr->theta -= 2*M_PI;
    }
    while(curPtr->theta < 0) {
        curPtr->theta += 2*M_PI;
    }

    //path length migrated need to be adjusted with F
    double path_elongation;
    {
        Mat_DP theta_vec(2,1);
        theta_vec[0][0] = cos(curPtr->theta);
        theta_vec[0][1] = sin(curPtr->theta);
        Mat_DP M = multiply(F, theta_vec);
        path_elongation = sqrt(M[0][0]*M[0][0] + M[0][1]*M[0][1]);
    }
    curPtr->xx += curPtr->speed * time_step * cos(curPtr->theta) / path_elongation;
    curPtr->yy += curPtr->speed * time_step * sin(curPtr->theta) / path_elongation;
     
    //(int)ÊÇÈ¡ÕûÊý²¿·Ö£¬²»ÊÇËÄÉáÎåÈë 
    if(curPtr->xx < 0)curPtr->xx = 0;
    if(curPtr->yy < 0)curPtr->yy = 0;
    if(curPtr->xx >= mXstep)curPtr->xx = mXstep-1;
    if(curPtr->yy >= mYstep)curPtr->yy = mYstep-1;
    
    return;          
}

void Flist::output_fibroblast(){
    Mat_DP &fibroblast_density = *fibroblast_density_ptr;
    Mat_DP &speedfield = *speedfield_ptr;
    Mat_DP &thetafield = *thetafield_ptr;
    
    for(size_t i = 0; i < mYstep; i++) {
        for(size_t j = 0; j < mXstep; j++) {
            fibroblast_density[i][j] = 0.0;
            speedfield[i][j] = 0.0;
            thetafield[i][j] = 0.0;
        }
    }

    cout << "output fibroblast........this may take a while.....\n" << endl;
    for(auto &item : mFCellMap) {
        DP cell_xx = item.second.xx;
        DP cell_yy = item.second.yy;
        DP cell_speed = item.second.speed;
        DP cell_theta = item.second.theta;
        int starty, endy, startx, endx, a = 45, b = 36;
        DP c = sqrt(a*a - b*b);
        starty = cell_yy<a ? 0 : (int)cell_yy - a;
        endy = (cell_yy > mYstep-a) ? mYstep-1 : (int)cell_yy + a;
        startx = cell_xx<a ? 0 : (int)cell_xx - a;
        endx = (cell_xx > mXstep-a) ? mXstep-1 : (int)cell_xx + a;

        DP focus1x = cell_xx - c*cos(cell_theta);
        DP focus1y = cell_yy - c*sin(cell_theta);
        DP focus2x = cell_xx + c*cos(cell_theta);
        DP focus2y = cell_yy + c*sin(cell_theta);

        DP distance;
        for(int i = starty; i <= endy; i++){
            for(int j = startx; j <= endx; j++){
                distance = sqrt(pow(i-focus1y,2) + pow(j-focus1x,2)) +
                           sqrt(pow(i-focus2y,2) + pow(j-focus2x,2));
                double current_density = 0.5 * cos(M_PI*pow(distance-2*c, 0.5)/pow(2*a-2*c, 0.5) ) + 0.5;
                if (distance <= 2*a && fibroblast_density[i][j] < current_density) {
                    fibroblast_density[i][j] = current_density;
                    speedfield[i][j] = cell_speed;
                    thetafield[i][j] = cell_theta;
                }
            }
        }
    }

    static int out_i = 0;
    char file_name[21];
    sprintf(file_name, "output/cell%05d.BMP", out_i);
    BMP::output_BMP(file_name, 14, fibroblast_density, mXstep, mYstep);
    out_i++;
}
