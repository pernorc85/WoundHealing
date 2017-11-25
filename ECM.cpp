#include "ECM.h"
#include <string>
#include <cmath>
#include <iostream> 
#include "BMP.h" 
#include "Flist.h"
using namespace std;

void ECM::initiate(int wound_radius){
    int i,j;
    DP theta;
    for(i=0;i<ECM_ystep/5;i++){//from bottom up
        for(j=0;j<ECM_xstep/5;j++){//from left to right
            double dist_from_center = sqrt(pow(j*5-ECM_xstep/2, 2) + pow((i*5-ECM_ystep/2)/0.8, 2));
            if(dist_from_center >= 450){
                collagen_density[i][j] = 1; 
                fibronectin_density[i][j] = 0;
            }
            else if(dist_from_center < 450 && 
                    dist_from_center > 350 ){      
                collagen_density[i][j] = sin( (dist_from_center-400)/50*M_PI/2 )*0.5 + 0.5;
                fibronectin_density[i][j] = 1 - collagen_density[i][j];
            } 
            else{
                collagen_density[i][j] = 0;  
                fibronectin_density[i][j] = 1; 
            }    
        }
    }  
    
    for(i=0;i<ECM_ystep/10;i++){
        for(j=0;j<ECM_xstep/10;j++){
            theta = rand()%180*M_PI/180;//collagenºÍfibroblastµÄ·½Ïò¶¼ÓÃ»¡¶È±íÊ¾
            //ÓÃ360Ì«Âé·³ÁË£¬Ö±½ÓÓÃ180 
            collagen[i][j]=theta;
        }
    } 
}


void ECM::collagen_orientation(Flist& fibroblastList)
{
    int i,j;
    int L = 10;
    DP ftheta;
    DP omega, omegasum;
    DP kappa1 = 20, kappa2 = 5;
    DP pc = 0.44, dc = 0.44, df = 0.6;

    for(i=0;i<ECM_ystep/10;i++){//from bottom up
        for(j=0;j<ECM_xstep/10;j++){//from left to right
            ftheta = 0;
            omegasum = 0;
            std::list<Fibroblast>::iterator curPtr = fibroblastList.mFCells.begin();
            for(; curPtr != fibroblastList.mFCells.end(); curPtr++){
                if(abs(curPtr->yy-i*10)<=L && abs(curPtr->xx-j*10)<=L)
                    omega = (1-abs(curPtr->yy-i*10)/L) * (1-abs(curPtr->xx-j*10)/L);
                else omega = 0;
                omegasum += omega;
                //be maticulous here-----------------------------------------------
                if(omega != 0){
                    if(abs(curPtr->theta - collagen[i][j]) <= M_PI/2)//case 1&4; part of 2; part of 3; part of 6;
                        ftheta += omega * curPtr->theta;
                    else if(abs(curPtr->theta - collagen[i][j]) > M_PI/2 && abs(curPtr->theta - collagen [i][j]) <= M_PI){
                         if(curPtr->theta <= M_PI/2 && collagen[i][j] > M_PI/2 && collagen[i][j] <= M_PI)//part of 2
                             ftheta += omega * (curPtr->theta + M_PI);
                         else ftheta += omega * (curPtr->theta - M_PI);//part of 3; part of 6;
                    }
                    else if(curPtr->theta > M_PI*3/2 && collagen[i][j] < M_PI/2){//case 7
                         if(collagen[i][j]+ 2*M_PI -curPtr->theta < M_PI/2)
                             ftheta += omega * (curPtr->theta - 2*M_PI);
                         else ftheta += omega * (curPtr->theta - M_PI); 
                    }
                    else ftheta += omega * (curPtr->theta - M_PI);//case 5&8 
                }
                //----------------------------------------------------------------- 
            } 
               
            if(omegasum != 0){                
                ftheta /= omegasum;  
                DP increase = tlength*kappa1*abs(omegasum) * sin(ftheta-collagen[i][j]);
                collagen[i][j] += increase; 
                while(collagen[i][j] >= M_PI) {
                    collagen[i][j] -= M_PI;
                }
                while(collagen[i][j] < 0) {
                    collagen[i][j] += M_PI;
                }
            } 
            
//            td_x = tissue_displacement_x[i*10][j*10];
//            td_y = tissue_displacement_y[i*10][j*10];
//            if(td_x != 0){
//                if(td_x > 0 && td_y >= 0)tensiontheta = atan(td_y/td_x);//1st phase 
//                else if(td_x < 0 && td_y >= 0)tensiontheta = PI + atan(td_y/td_x);//2nd phase
//                else if(td_x < 0 && td_y < 0)tensiontheta = atan(td_y/td_x);//3nd phase
//                else if(td_x > 0 && td_y < 0)tensiontheta = PI + atan(td_y/td_x);//4th phase
//            }  
//            else{
//                if(td_y > 0)tensiontheta =  PI/2;
//                else if(td_y < 0)tensiontheta = M_PI/2;
//                else tensiontheta = -1;
//            }                
//   
//            if(tensiontheta != -1 && sqrt(pow(td_x,2)+pow(td_y,2))>200){
//                DP increase; 
//                if(abs(tensiontheta-collagen[i][j])<=PI/2)
//                    increase = tlength*kappa2*sqrt(pow(td_x/200,2)+pow(td_y/200,2))*sin(tensiontheta-collagen[i][j]);
//                else if(tensiontheta>PI/2)
//                    increase = tlength*kappa2*sqrt(pow(td_x/200,2)+pow(td_y/200,2))*sin(tensiontheta-PI-collagen[i][j]);
//                else increase = tlength*kappa2*sqrt(pow(td_x/200,2)+pow(td_y/200,2))*sin(tensiontheta+PI-collagen[i][j]);
//
//                collagen[i][j] += increase;
//                cout <<"increase="<<increase<<endl;
//                collagen[i][j] = (  (int)(collagen[i][j]*180/PI+180)%180  )*PI/180;
//            }                                                         
        }
    } 
          
    for(i=0;i<ECM_ystep/5;i++){//from bottom up
        for(j=0;j<ECM_xstep/5;j++){//from left to right
            omegasum = 0;
            std::list<Fibroblast>::iterator curPtr = fibroblastList.mFCells.begin();
            for(; curPtr != fibroblastList.mFCells.end(); curPtr++){
                if(abs(curPtr->yy-i*5)<=L && abs(curPtr->xx-j*5)<=L)
                    omega = (1-abs(curPtr->yy-i*5)/L) * (1-abs(curPtr->xx-j*5)/L);
                else omega = 0;
                omegasum += omega;
            } 
            collagen_density[i][j]    +=  tlength * (pc - dc * collagen_density[i][j]) * omegasum;   
            fibronectin_density[i][j] += -tlength * df * fibronectin_density[i][j] * omegasum;          
        }
    }
}




void ECM::output_collagen(Flist& fibroblastList)
{
    const DP dl=1;
    int k,r;
    DP x,y,xnew,ynew;
    cout << "xstep = " << ECM_xstep << "ystep = " << ECM_ystep <<endl;
    Mat_DP vf(ECM_ystep,ECM_xstep); 
    
    srand(time(NULL));
    for(int i=0;i<ECM_ystep;i++){
        for(int l=0;l<ECM_xstep;l++){            
            vf[i][l] = 0;                      
        }       
    }
                
    for(int j=0;j<ECM_ystep/10;j++){ 
        for(int l=0;l<ECM_xstep/10;l++){
            y = j*10;
            x = l*10;
            for(k=0;k<100;k++){
                //printf("x = %f, y = %f\n",x,y);                
                xnew = x + dl*cos(collagen[(int)(y/10)][(int)(x/10)]);
                if(xnew < 0)xnew = 0;
                if(xnew > ECM_xstep-1)xnew = ECM_xstep-1;
                ynew = y + dl*sin(collagen[(int)(y/10)][(int)(x/10)]);
                if(ynew < 0)ynew= 0;
                if(ynew > ECM_ystep-1)ynew = ECM_ystep-1;
                x = xnew;
                y = ynew;
                vf[(int)y%ECM_ystep][(int)x%ECM_xstep] = collagen_density[(int)(y/5)][(int)(x/5)];                
            }    
        }
    }
    
    
    Mat_DP cellmatrix(ECM_ystep,ECM_xstep); 
    
    for(int j=0; j<ECM_ystep; j++){
        for(int l=0; l<ECM_xstep; l++){
            cellmatrix[j][l] = 0;
        }
    }
    std::list<Fibroblast>::iterator curPtr = fibroblastList.mFCells.begin();
    for(; curPtr != fibroblastList.mFCells.end(); curPtr++){
        cellmatrix[(int)curPtr->yy][(int)curPtr->xx] = 1;
        if((int)curPtr->xx >= 2 && (int)curPtr->xx < ECM_xstep -2 &&
           (int)curPtr->yy >= 2 && (int)curPtr->yy < ECM_ystep -2){
            for(int j = (int)curPtr->yy-2;j <= (int)curPtr->yy+2;j++){
                for(int l = (int)curPtr->xx-2;l <= (int)curPtr->xx+2;l++){
                    cellmatrix[j][l] = 1;
                }
            } 
        }
        //cout << (int)curPtr->xx << (int)curPtr->yy << endl;
    }
    
    
    //************output to picture*****************************************
    //Ö»Êä³öÒ»·ùÍ¼ 
    static int out_i = 0;
    char file_name[14]="coll00000.BMP";
    sprintf(file_name, "coll%05d.BMP", out_i++);
    BMP::output_BMP2(file_name, 14, vf, cellmatrix, ECM_xstep, ECM_ystep);
   
    return;
}
    
