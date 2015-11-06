#include <string>
#include <cmath>
#include <iostream> 
#include "ECM_Flist_BMP.h"
using namespace std;   


Flist::Flist()
{
    Froot = NULL;
    curPtr = Froot;
    speed = 15;
}

Flist::Flist(int FNinit)
{
    
    int i; 
    DP x,y,theta;
    FcellPtr newPtr,prePtr,curPtr;
    
    Froot = NULL;
    speed = 15;
    for(i=0;i<FNinit;i++){               
        x = rand()%1000;
        y = rand()%1000;
        theta = (rand()%360)*PI/180;
        
        //if(pow(x-500,2)+pow(y-500,2) > 160000 || pow(x-500,2)+pow((y-500)/1,2) < 90000 ){ 
        if(pow(x-500.0,2)+pow((y-500.0)/0.8,2) >= 122500){// || pow((x-500.0)/0.3,2)+pow(y-500.0,2) >=202500){                      
       	    newPtr = (FcellPtr)malloc(sizeof(Fcell));
     	    newPtr->xx = x;
     	    newPtr->yy = y; 
    	    newPtr->theta = theta;
            newPtr->next = NULL;
    	    
     	    prePtr = NULL;
            curPtr = Froot;
            while(curPtr != NULL){
                prePtr = curPtr;
    	        curPtr = curPtr->next;  
            }      
            if(prePtr == NULL){
                newPtr->next =Froot;
                Froot = newPtr;
            }
            else{
                prePtr->next = newPtr;
                newPtr->next = curPtr;
            } 
        }       
    }          
    return;
}

void Flist::Flist_move(Mat_DP& gradx_cytokine, Mat_DP& grady_cytokine, ECM& extraCellularMatrix)
{
    DP gradx, grady;
    curPtr = Froot;         
    while(curPtr != NULL){ 
        gradx = gradx_cytokine[(int)curPtr->yy][(int)curPtr->xx];
        grady = grady_cytokine[(int)curPtr->yy][(int)curPtr->xx]; 
        DP fdensity = extraCellularMatrix.fibronectin_density[(int)curPtr->yy][(int)curPtr->xx];
        DP cdensity = extraCellularMatrix.collagen_density[(int)curPtr->yy][(int)curPtr->xx];
        
        Fcell_move(gradx, grady, fdensity, cdensity, extraCellularMatrix.collagen);                         
        curPtr = curPtr -> next; 
    }
}

void Flist::Fcell_move(DP gradx, DP grady, DP fdensity, DP cdensity, Mat_DP& collagen)
{
    DP rou1 = 0.5, rou2 = 0.3;
//************************chemotaxis********************************************
    DP gradtheta;
    DP speed_taken;
    
    if(cdensity < 0.7)speed_taken = speed * (1 + 2*fdensity)/3;
    else speed_taken = speed * (1 + 2*fdensity)/3 * (1.7 - cdensity);
    curPtr->speed = speed_taken;
    
    DP grad=min(1.0,1000*sqrt(pow(gradx,2)+pow(grady,2)));
    rou2=rou2*grad;
    rou1=rou1*cdensity;
//    cout<<"rou2="<<rou2<<endl;
    if(gradx != 0){
        if(gradx > 0 && grady >= 0)gradtheta = atan(grady/gradx);//1st phase 
        else if(gradx < 0 && grady >= 0)gradtheta = PI + atan(grady/gradx);//2nd phase
        else if(gradx < 0 && grady < 0)gradtheta = PI + atan(grady/gradx);//3nd phase
        else if(gradx > 0 && grady < 0)gradtheta = 2*PI + atan(grady/gradx);//4th phase
    }  
    else{
         if(grady > 0)gradtheta =  PI/2;
         else if(grady < 0)gradtheta = PI*3/2;
         else gradtheta = -1;
    }    
    
//***********************guidance of collagen***********************************      
    int i,j,k,x,y;
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
    
    for(j=max(gridy-2,0);j<=min(gridy+2,1000/10-1);j++){
        for(k=max(gridx-2,0);k<min(gridx+2,1000/10-1);k++){
            collagen_temp = collagen[j][k];
            //be maticulous here-----------------------------------------------
            if(abs(curPtr->theta - collagen_temp) <= PI/2)//case 1&4; part of 2; part of 3; part of 6;
                collagen_taken = collagen_temp;
            else if(abs(curPtr->theta - collagen_temp) > PI/2 && abs(curPtr->theta - collagen_temp) <= PI){
                if(curPtr->theta <= PI/2 && collagen_temp > PI/2 && collagen_temp <= PI)//part of 2
                     collagen_taken = collagen_temp - PI;
                else collagen_taken = collagen_temp + PI;//part of 3; part of 6;
            }
            else if(curPtr->theta > PI*3/2 && collagen_temp < PI/2){//case 7
                if(collagen_temp+2*PI-curPtr->theta < PI/2)
                    collagen_taken = collagen_temp + 2*PI;
                else collagen_taken = collagen_temp + PI;  
            }
            else collagen_taken = collagen_temp + PI;//case 5&8 
            //-----------------------------------------------------------------
            //calculate the weight---------------------------------------------
            l_y = (curPtr->yy-(j-2)*10) * (curPtr->yy-(j-1)*10) * (curPtr->yy-(j+1)*10) * (curPtr->yy-(j+2)*10)/20/10/10/20;
            l_x = (curPtr->xx-(k-2)*10) * (curPtr->xx-(k-1)*10) * (curPtr->xx-(k+1)*10) * (curPtr->xx-(k+2)*10)/20/10/10/20;
            l_sum += l_y*l_x;
            collagen_sum += collagen_taken*l_y*l_x;
        }
    }
    if(l_sum != 0)collagen_avg = collagen_sum/l_sum;

    if(gradtheta >= 0){
        if(gradtheta < PI/2 && curPtr->theta > 3/2*PI && curPtr->theta < 2*PI){
            curPtr->theta = (1 - rou1 - rou2)*curPtr->theta + rou1*collagen_avg + rou2*(gradtheta + 2*PI);
        }
        else if(gradtheta > 3/2*PI && gradtheta < 2*PI && curPtr->theta < PI/2){
            curPtr->theta = (1 - rou1 - rou2)*curPtr->theta + rou1*collagen_avg + rou2*(gradtheta - 2*PI);
        }
        else curPtr->theta = (1 - rou1 - rou2)*curPtr->theta + rou1*collagen_avg + rou2*gradtheta;
    }
    else if(gradtheta == -1) curPtr->theta = (1 - rou1)*curPtr->theta + rou1 * collagen_avg;

    curPtr->theta = (  (int)(curPtr->theta*180/PI+360)%360  )*PI/180; 
    curPtr->xx += speed_taken * tlength * cos(curPtr->theta);
    curPtr->yy += speed_taken * tlength * sin(curPtr->theta);
     
    //(int)ÊÇÈ¡ÕûÊý²¿·Ö£¬²»ÊÇËÄÉáÎåÈë 
    if(curPtr->xx < 0)curPtr->xx = 0;
    if(curPtr->yy < 0)curPtr->yy = 0;
    if(curPtr->xx >= 1000)curPtr->xx = 1000-1;
    if(curPtr->yy >= 1000)curPtr->yy = 1000-1;
    
    return;          
}
