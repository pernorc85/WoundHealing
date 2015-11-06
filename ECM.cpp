#include <string>
#include <cmath>
#include <iostream> 
#include "ECM_Flist_BMP.h"
#include "BMP.h" 
using namespace std;

void ECM::initiate(){
    int i,j;
    DP theta;
    for(i=0;i<ECM_ystep;i++){//from bottom up
        for(j=0;j<ECM_xstep;j++){//from left to right
            if(pow(j-500.0,2)+pow((i-500.0)/0.8,2) >= 202500){// || pow(j-500.0,2)+pow((i-500.0)/1,2) < 62500){
                collagen_density[i][j] = 1; 
                fibronectin_density[i][j] = 0;
            }
            else if(pow(j-500.0,2)+pow((i-500.0)/0.8,2) < 202500 && pow(j-500.0,2)+pow((i-500.0)/0.8,2) > 122500){      
                collagen_density[i][j] = sin( (sqrt(pow(j-500.0,2)+pow((i-500.0)/0.8,2))-400)/50*PI/2 )*0.5 + 0.5;
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
            theta = rand()%180*PI/180;//collagenºÍfibroblastµÄ·½Ïò¶¼ÓÃ»¡¶È±íÊ¾
            //ÓÃ360Ì«Âé·³ÁË£¬Ö±½ÓÓÃ180 
            collagen[i][j]=theta;
        }
    } 
}


void ECM::collagen_orientation(Mat_DP& PDGF, Flist& fibroblastList)
{
    int i,j;
    int L = 10;
    DP ftheta, f;
    DP omega, omegasum;
    DP kappa1 = 20, kappa2 = 5;
    DP td_x, td_y, tensiontheta;
    DP pc = 0.44, dc = 0.44, df = 0.6;
    FcellPtr curPtr; 

    for(i=0;i<ECM_ystep/10;i++){//from bottom up
        for(j=0;j<ECM_xstep/10;j++){//from left to right
            f = 0;
            omegasum = 0;
            curPtr = fibroblastList.Froot;
            while(curPtr != NULL){
                if(abs(curPtr->yy-i*10)<=L && abs(curPtr->xx-j*10)<=L)
                    omega = (1-abs(curPtr->yy-i*10)/L) * (1-abs(curPtr->xx-j*10)/L);
                else omega = 0;
                omegasum += omega;
                //be maticulous here-----------------------------------------------
                if(omega != 0){
                    if(abs(curPtr->theta - collagen[i][j]) <= PI/2)//case 1&4; part of 2; part of 3; part of 6;
                        f += omega * curPtr->theta;
                    else if(abs(curPtr->theta - collagen[i][j]) > PI/2 && abs(curPtr->theta - collagen [i][j]) <= PI){
                         if(curPtr->theta <= PI/2 && collagen[i][j] > PI/2 && collagen[i][j] <= PI)//part of 2
                             f += omega * (curPtr->theta + PI);
                         else f += omega * (curPtr->theta - PI);//part of 3; part of 6;
                    }
                    else if(curPtr->theta > PI*3/2 && collagen[i][j] < PI/2){//case 7
                         if(collagen[i][j]+ 2*PI -curPtr->theta < PI/2)
                             f += omega * (curPtr->theta - 2*PI);
                         else f += omega * (curPtr->theta - PI); 
                    }
                    else f += omega * (curPtr->theta - PI);//case 5&8 
                }
                //----------------------------------------------------------------- 
                curPtr = curPtr->next;
            } 
               
            if(omegasum != 0){                
                ftheta = f/omegasum;  
                DP increase = tlength*kappa1*abs(omegasum) * sin(ftheta-collagen[i][j]);
                collagen[i][j] += increase; 
                collagen[i][j] = (  (int)(collagen[i][j]*180/PI+180)%180  )*PI/180; 
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
//                else if(td_y < 0)tensiontheta = PI/2;
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
           
    for(i=0;i<ECM_ystep;i++){//from bottom up
        for(j=0;j<ECM_xstep;j++){//from left to right
            f = 0;
            omegasum = 0;
            curPtr = fibroblastList.Froot;
            while(curPtr != NULL){
                if(abs(curPtr->yy-i)<=L && abs(curPtr->xx-j)<=L)
                    omega = (1-abs(curPtr->yy-i)/L) * (1-abs(curPtr->xx-j)/L);
                else omega = 0;
                omegasum += omega;
                curPtr = curPtr->next;
            } 
            
            collagen_density[i][j] +=     tlength * (pc - dc * collagen_density[i][j]) * omegasum;   
            //PDGF[i][j]      -= PDGF[i][j]*tlength * (pc - dc * collagen_density[i][j]) * omegasum;
            fibronectin_density[i][j] += -tlength * df * fibronectin_density[i][j] * omegasum;          
        }
    }                               
}




void ECM::output_collagen(Flist& fibroblastList)
{
    const DP dl=1;
    int j,l,k,r;
    DP x,y,xnew,ynew;
    cout << "xstep = " << ECM_xstep << "ystep = " << ECM_ystep <<endl;
    Mat_DP vf(ECM_ystep,ECM_xstep); 
    
    srand(time(NULL));
    for(j=0;j<ECM_ystep;j++){
        for(l=0;l<ECM_xstep;l++){            
            vf[j][l] = 0;                      
        }       
    }
                
    for(j=0;j<ECM_ystep/10;j++){ 
        for(l=0;l<ECM_xstep/10;l++){
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
                vf[(int)y%ECM_ystep][(int)x%ECM_xstep] = collagen_density[(int)y][(int)x];                
            }    
        }
    }
    
    
    FcellPtr curPtr;
    Mat_DP cellmatrix(ECM_ystep,ECM_xstep); 
    
    for(j=0; j<ECM_ystep; j++){
        for(l=0; l<ECM_xstep; l++){
            cellmatrix[j][l] = 0;
        }
    }
    curPtr = fibroblastList.Froot;
    while(curPtr != NULL){
        cellmatrix[(int)curPtr->yy][(int)curPtr->xx] = 1;
        if((int)curPtr->xx >= 2 && (int)curPtr->xx < ECM_xstep -2 &&
           (int)curPtr->yy >= 2 && (int)curPtr->yy < ECM_ystep -2){
            for(j = (int)curPtr->yy-2;j <= (int)curPtr->yy+2;j++){
                for(l = (int)curPtr->xx-2;l <= (int)curPtr->xx+2;l++){
                    cellmatrix[j][l] = 1;
                }
            } 
        }
        //cout << (int)curPtr->xx << (int)curPtr->yy << endl;
        curPtr = curPtr->next;
    }
    
    
    //************output to picture*****************************************
    //Ö»Êä³öÒ»·ùÍ¼ 
    static int out_i = 0;
    char  file_name[14]="coll00000.BMP";
    file_name[8] = 48+out_i%10;        file_name[7] = 48+out_i/10%10;
    file_name[6] = 48+out_i/10/10%10;  file_name[5] = 48+out_i/10/10/10%10;
    file_name[4] = 48+out_i/10/10/10/10%10;             
    fp_BMP = fopen(file_name,"w");
    Prepare_BMP_Format(ECM_xstep,ECM_ystep);
    int jj;
    for(j=0;j<ECM_ystep;j++) {
        for(l=0;l<ECM_xstep;l++) { 
            jj = int((vf[j][l] - Minimal)/Interval); 
            if(cellmatrix[j][l] == 1) fwrite(&b_yellow, 1,3,fp_BMP);
            else{
    	        if(jj==0)   fwrite(&b0, 1,3,fp_BMP);
                if(jj==1)   fwrite(&b1, 1,3,fp_BMP);
                if(jj==2)   fwrite(&b2, 1,3,fp_BMP);
    	        if(jj==3)   fwrite(&b3, 1,3,fp_BMP);
    	        if(jj==4)   fwrite(&b4, 1,3,fp_BMP);
    	        if(jj==5)   fwrite(&b5, 1,3,fp_BMP);
    	        if(jj==6)   fwrite(&b6, 1,3,fp_BMP);
    	        if(jj==7)   fwrite(&b7, 1,3,fp_BMP);
    	        if(jj==8)   fwrite(&b8, 1,3,fp_BMP);
    	        if(jj==9)   fwrite(&b9, 1,3,fp_BMP);
    	        if(jj==10)  fwrite(&b10, 1,3,fp_BMP);
    	        if(jj==11)  fwrite(&b11, 1,3,fp_BMP);
    	        if(jj==12)  fwrite(&b12, 1,3,fp_BMP);
    	        if(jj==13)  fwrite(&b13, 1,3,fp_BMP);
    	        if(jj==14)  fwrite(&b14, 1,3,fp_BMP);
    	        if(jj==15)  fwrite(&b15, 1,3,fp_BMP);
    	        if(jj < 0 )  fwrite(&b_0,  1,3,fp_BMP);  
    	        if(jj > 15)  fwrite(&b_16, 1,3,fp_BMP); 
            }
        }
    }                   		    
    fclose(fp_BMP);    
    out_i++;
 
    return;
}
    
