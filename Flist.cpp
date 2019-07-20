#include "Flist.h"
#include <string>
#include <cmath>
#include <iostream>
#include "chemokine.h"
#include "ECM.h"
#include "BMP.h"

using namespace std;   

extern Mat_DP *fibroblast_density_ptr;
extern Mat_DP *speedfield_ptr;
extern Mat_DP *thetafield_ptr;

Flist::Flist(int xstep, int ystep):
        mXstep(xstep), mYstep(ystep), speed(15) {}

void Flist::initialize(int FDensity, int wound_radius){
    DP x,y,theta;
    int FNumber = FDensity * (mXstep/1000.0) * (mYstep/1000.0);
    for(int i=0;i<FNumber;i++){               
        x = rand()%mXstep;
        y = rand()%mYstep;
        theta = (rand()%360)*M_PI/180.0;
        
        //if(pow(x-mXstep/2,2)+pow(y-mYstep,2) > 160000 || pow(x-mXstep,2)+pow(y-mYstep,2) < 90000 ){ 
        if(pow(x-mXstep/2,2)+pow(y-mYstep/2,2) >= wound_radius*wound_radius){// || pow((x-mXstep/2)/0.3,2)+pow(y-mYstep/2,2) >=202500){                      
    	    mFCells.push_back(Fibroblast(x,y, theta));    
        }       
    }          
    return;
}

void Flist::initialize(int FDensity, double a, double b){
    DP x,y,theta;
    int FNumber = FDensity * (mXstep/1000.0) * (mYstep/1000.0);
    for(int i=0;i<FNumber;i++){               
        x = rand()%mXstep;
        y = rand()%mYstep;
        theta = (rand()%360)*M_PI/180.0;
        
        //if(pow(x-mXstep/2,2)+pow(y-mYstep,2) > 160000 || pow(x-mXstep,2)+pow(y-mYstep,2) < 90000 ){ 
        if(pow(x-mXstep/2,2)+pow(a/b*(y-mYstep/2),2) >= a*a){
    	    mFCells.push_back(Fibroblast(x,y, theta));    
        }       
    }          
    return;
}

void Flist::Flist_move(Chemokine& PDGF, ECM& extraCellularMatrix, DP time_step)
{
    DP gradx, grady;
    DP fdensity, cdensity;
    std::list<Fibroblast>::iterator it = mFCells.begin();    
    for(;it != mFCells.end(); it++) {
        gradx = (*PDGF.gradx)[(int)it->yy][(int)it->xx];
        grady = (*PDGF.grady)[(int)it->yy][(int)it->xx]; 
        fdensity = extraCellularMatrix.fibronectin_density[(int)(it->yy/5)][(int)(it->xx/5)];
        cdensity = extraCellularMatrix.collagen_density[(int)(it->yy/5)][(int)(it->xx/5)];
        
        Fcell_move(it, gradx, grady, fdensity, cdensity, extraCellularMatrix.collagen, time_step);                         
    }
}

void Flist::Fcell_move(std::list<Fibroblast>::iterator curPtr, DP gradx, DP grady, DP fdensity, DP cdensity, Mat_DP& collagen, DP time_step)
{
    DP rou1 = 0.3, rou2 = 0.4;
    //************************chemotaxis********************************************
    DP gradtheta = -1;
    DP speed_taken = speed;
    
    if(cdensity < 0.7)speed_taken = speed * (1 + 2*fdensity)/3;
    else speed_taken = speed * (1 + 2*fdensity)/3 * (1.7 - cdensity);
    curPtr->speed = speed_taken;
    
    DP grad=min(1.0,1000*sqrt(pow(gradx,2)+pow(grady,2)));
    rou1=rou1*grad;
    rou2=rou2*cdensity;
    //cout<<"rou2="<<rou2<<endl;
    if(gradx != 0){
        if(gradx > 0 && grady >= 0)gradtheta = atan(grady/gradx);//1st phase 
        else if(gradx < 0 && grady >= 0)gradtheta = M_PI + atan(grady/gradx);//2nd phase
        else if(gradx < 0 && grady < 0)gradtheta = M_PI + atan(grady/gradx);//3nd phase
        else if(gradx > 0 && grady < 0)gradtheta = 2*M_PI + atan(grady/gradx);//4th phase
    } else{
         if(grady > 0)gradtheta =  M_PI/2;
         else if(grady < 0)gradtheta = M_PI*3/2;
         else gradtheta = -1;
    }    
    
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
        }else if(gradtheta > 3/2*M_PI && gradtheta < 2*M_PI && curPtr->theta < M_PI/2){
            curPtr->theta = (1 - rou1 - rou2)*curPtr->theta + rou2*collagen_avg + rou1*(gradtheta - 2*M_PI);
        }else {
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
    curPtr->xx += speed_taken * time_step * cos(curPtr->theta);
    curPtr->yy += speed_taken * time_step * sin(curPtr->theta);
     
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

    auto it = mFCells.begin();
    cout << "output fibroblast........this may take a while.....\n" << endl;
    for(; it != mFCells.end(); it++) {
        int starty, endy, startx, endx, a = 45, b = 36;
        DP c = sqrt(a*a - b*b);
        starty = it->yy<a ? 0 : (int)it->yy - a;
        endy = (it->yy > mYstep-a) ? mYstep-1 : (int)it->yy + a;
        startx = it->xx<a ? 0 : (int)it->xx - a;
        endx = (it->xx > mXstep-a) ? mXstep-1 : (int)it->xx + a;

        DP focus1x = it->xx - c*cos(it->theta);
        DP focus1y = it->yy - c*sin(it->theta);
        DP focus2x = it->xx + c*cos(it->theta);
        DP focus2y = it->yy + c*sin(it->theta);

        DP distance;
        for(int i = starty; i <= endy; i++){
            for(int j = startx; j <= endx; j++){
                distance = sqrt(pow(i-focus1y,2) + pow(j-focus1x,2)) +
                           sqrt(pow(i-focus2y,2) + pow(j-focus2x,2));
                double current_density = 0.5 * cos(M_PI*pow(distance-2*c, 0.5)/pow(2*a-2*c, 0.5) ) + 0.5;
                if (distance <= 2*a && fibroblast_density[i][j] < current_density) {
                    fibroblast_density[i][j] = current_density;
                    speedfield[i][j] = it->speed;
                    thetafield[i][j] = it->theta;
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
