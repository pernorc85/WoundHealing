#include "ECM.h"
#include <string>
#include <cmath>
#include <tuple>
#include <iostream> 
#include "BMP.h" 
#include "Flist.h"
#include "util.h"

using namespace std;

void ECM::initialize(double a/*major semi_axis*/, double b/*minor semi_axis*/){
    int i,j;
    double c =sqrt(a*a - b*b);
    for(i=0; i<ECM_ystep/5; i++){//from bottom up
        for(j=0; j<ECM_xstep/5; j++){
            double dist_from_focus = sqrt(pow(j*5-ECM_xstep/2, 2) + pow(a/b*(i*5-ECM_ystep/2), 2)); 
            if(dist_from_focus > a + 50) {
                collagen_density[i][j] = 1.0;
                fibronectin_density[i][j] = 0.0;
            } else if(dist_from_focus < a + 50 && dist_from_focus > a - 50) {
                collagen_density[i][j] = sin( (dist_from_focus-a)/50 * M_PI/2 )*0.5 + 0.5;
                fibronectin_density[i][j] = 1 - collagen_density[i][j];
            } else {
                collagen_density[i][j] = 0.0;
                fibronectin_density[i][j] = 1.0;
            }
        }
    }

    DP theta;
    for(i=0; i<ECM_ystep/10; i++){
        for(j=0; j<ECM_xstep/10; j++){
            theta = rand()%180 * M_PI/180;//collagenºÍfibroblastµÄ·½Ïò¶¼ÓÃ»¡¶È±íÊ¾
            //ÓÃ360Ì«Âé·³ÁË£¬Ö±½ÓÓÃ180 
            collagen[i][j] = theta;
            stretch_history[i][j] = 0.0;
        }
    } 
}

void ECM::initialize_rectangle(double a, double b){
    float xcenter = ECM_xstep * 0.5;
    float ycenter = ECM_ystep * 0.5;
    for(int i=0; i<ECM_ystep/5; i++){//from bottom up
        for(int j=0; j<ECM_xstep/5; j++){
            if(abs(j*5-xcenter)>a*0.5+50 || abs(i*5-ycenter)>b*0.5+50) {
                collagen_density[i][j] = 1.0;
                fibronectin_density[i][j] = 0.0;
            } else if(abs(j*5-xcenter)>a*0.5-50 || abs(i*5-ycenter)>b*0.5-50) {
                if(abs(j*5-xcenter)<=a*0.5+50 && abs(j*5-xcenter)>a*0.5-50){ 
                    double frame = abs(j*5-xcenter) - a*0.5;
                    collagen_density[i][j] = sin( frame*0.02 * M_PI/2 )*0.5 + 0.5;
                    fibronectin_density[i][j] = 1 - collagen_density[i][j];
                }
                if(abs(i*5-ycenter)<=b*0.5+50 && abs(i*5-ycenter)>b*0.5-50){
                    double frame = abs(i*5-ycenter) - b*0.5;
                    double tmp = sin( frame*0.02 * M_PI/2 )*0.5 + 0.5;
                    if(tmp > collagen_density[i][j]) {
                        collagen_density[i][j] = tmp; 
                        fibronectin_density[i][j] = 1 - collagen_density[i][j];
                    }
                }
            } else {
                collagen_density[i][j] = 0.0;
                fibronectin_density[i][j] = 1.0;
            }
        }
    }

    DP theta;
    for(int i=0; i<ECM_ystep/10; i++){
        for(int j=0; j<ECM_xstep/10; j++){
            theta = rand()%180 * M_PI/180;//collagenºÍfibroblastµÄ·½Ïò¶¼ÓÃ»¡¶È±íÊ¾
            //ÓÃ360Ì«Âé·³ÁË£¬Ö±½ÓÓÃ180 
            collagen[i][j] = theta;
            stretch_history[i][j] = 0.0;
        }
    } 
}

void ECM::initialize(int wound_radius){
    int i,j;
    DP theta;
    for(i=0;i<ECM_ystep/5;i++){//from bottom up
        for(j=0;j<ECM_xstep/5;j++){//from left to right
            double dist_from_center = sqrt(pow(j*5-ECM_xstep/2, 2) + pow(i*5-ECM_ystep/2, 2));
            if(dist_from_center >= wound_radius + 50){
                collagen_density[i][j] = 1; 
                fibronectin_density[i][j] = 0;
            } else if(dist_from_center < wound_radius+50 && dist_from_center > wound_radius-50 ){
                collagen_density[i][j] = sin( (dist_from_center-wound_radius)/50*M_PI/2 )*0.5 + 0.5;
                fibronectin_density[i][j] = 1 - collagen_density[i][j];
            } else{
                collagen_density[i][j] = 0;  
                fibronectin_density[i][j] = 1; 
            }
            /* 
            collagen_density[i][j] = 1;
            fibronectin_density[i][j] = 0;
            */ 
        }
    }  
    
    for(i=0; i<ECM_ystep/10; i++){
        for(j=0; j<ECM_xstep/10; j++){
            theta = rand()%180 * M_PI/180;//collagenºÍfibroblastµÄ·½Ïò¶¼ÓÃ»¡¶È±íÊ¾
            //ÓÃ360Ì«Âé·³ÁË£¬Ö±½ÓÓÃ180 
            collagen[i][j] = theta;
            stretch_history[i][j] = 0.0;
        }
    } 
}

void ECM::collagen_orientation(Flist& fList, const Mat_DP& tdx, const Mat_DP& tdy, double time_step) {
    collagen_orientation_with_fibroblast(fList, time_step);
    collagen_orientation_under_tension(tdx, tdy, time_step);
    return;
}

void ECM::collagen_orientation_with_fibroblast(Flist& fibroblastList, double time_step)
{
    int L = 30;
    DP ftheta;
    DP omega, omegasum;
    DP kappa1 = 5, kappa2 = 5;
    DP pc = 0.44, dc = 0.44, df = 0.6;

    for(size_t i=0; i<ECM_ystep/10; i++){//from bottom up
        for(size_t j=0; j<ECM_xstep/10; j++){//from left to right
            ftheta = 0;
            omegasum = 0;
            std::list<Fibroblast>::iterator curPtr = fibroblastList.mFCells.begin();
            //TODO: This needs acceleration
            //     Going through the entire list of fibroblasts is too slow for 5mm x 5mm modeling space
            for(; curPtr != fibroblastList.mFCells.end(); curPtr++){
                //if(abs(curPtr->yy-i*10)<=L && abs(curPtr->xx-j*10)<=L)
                //    omega = (1-abs(curPtr->yy-i*10)/L) * (1-abs(curPtr->xx-j*10)/L);
                double dist = sqrt( pow(curPtr->yy-i*10, 2) + pow(curPtr->xx-j*10, 2) );
                omega = (dist <=L)? 1-dist/L : 0;
                omegasum += omega;
                //be maticulous here-----------------------------------------------
                if(omega != 0){
                    DP f_taken = curPtr->theta;
                    while(abs(f_taken - collagen[i][j]) > M_PI/2) {
                        if(f_taken > collagen[i][j])
                            f_taken -= M_PI;
                        else
                            f_taken += M_PI;
                    }
                    ftheta += omega * f_taken;
                    /*   
                    DP f_taken2 = 0; 
                    if(abs(curPtr->theta - collagen[i][j]) <= M_PI/2)//case 1&4; part of 2; part of 3; part of 6;
                        f_taken2 = curPtr->theta;
                    else if(abs(curPtr->theta - collagen[i][j]) > M_PI/2 && abs(curPtr->theta - collagen [i][j]) <= M_PI){
                         if(curPtr->theta <= M_PI/2 && collagen[i][j] > M_PI/2 && collagen[i][j] <= M_PI)//part of 2
                             f_taken2 = curPtr->theta + M_PI;
                         else f_taken2 = curPtr->theta - M_PI;//part of 3; part of 6;
                    }
                    else if(curPtr->theta > M_PI*3/2 && collagen[i][j] < M_PI/2){//case 7
                         if(collagen[i][j]+ 2*M_PI -curPtr->theta < M_PI/2)
                             f_taken2 = curPtr->theta - 2*M_PI;
                         else f_taken2 = curPtr->theta - M_PI; 
                    }
                    else f_taken2 = curPtr->theta - M_PI;//case 5&8 
                    */ 
                }
                //----------------------------------------------------------------- 
            } 
               
            if(omegasum != 0){                
                ftheta /= omegasum;  
                DP increase = time_step * kappa1 * abs(omegasum) * sin(ftheta-collagen[i][j]);
                collagen[i][j] += increase; 
                while(collagen[i][j] >= M_PI) {
                    collagen[i][j] -= M_PI;
                }
                while(collagen[i][j] < 0) {
                    collagen[i][j] += M_PI;
                }
            } 
            
        }
    } 
          
    for(size_t i=0; i<ECM_ystep/5; i++){//from bottom up
        for(size_t j=0; j<ECM_xstep/5; j++){//from left to right
            omegasum = 0;
            std::list<Fibroblast>::iterator curPtr = fibroblastList.mFCells.begin();
            for(; curPtr != fibroblastList.mFCells.end(); curPtr++){
                if(abs(curPtr->yy-i*5)<=L && abs(curPtr->xx-j*5)<=L)
                    omega = (1-abs(curPtr->yy-i*5)/L) * (1-abs(curPtr->xx-j*5)/L);
                else omega = 0;
                omegasum += omega;
            } 
            collagen_density[i][j] +=  time_step * (pc - dc * collagen_density[i][j]) * omegasum;   
            fibronectin_density[i][j] += -time_step * df * fibronectin_density[i][j] * omegasum;          
        }
    }
}

/*
 * F = U*R;
 * Rotation R doesn't matter, since rotation doesn't cause stress
 */
std::tuple<double, double, double>
calculate_principal_strain(Mat_DP &F){
    Mat_DP C = multiply(transpose(F), F);//left Cauchy strain tensor, which is a positive definite matrix
    Mat_DP eigen_vecs(2,2);
    Vec_DP eig = eigen_decomposition(C, eigen_vecs);
    double stretch1 = sqrt(eig[0]);
    double stretch2 = sqrt(eig[1]);
    double principal_strain_direction;
    if(stretch1 > stretch2) {//actually this is always the case
        principal_strain_direction = atan(eigen_vecs[1][0] / eigen_vecs[0][0]);
        if (principal_strain_direction < 0) principal_strain_direction += M_PI;
        return std::make_tuple(stretch1, stretch2, principal_strain_direction);
    } else {
        principal_strain_direction = atan(eigen_vecs[1][1] / eigen_vecs[0][1]);
        if (principal_strain_direction < 0) principal_strain_direction += M_PI;
        return std::make_tuple(stretch2, stretch1, principal_strain_direction);
    }
}

/**
 * Collagen fibers re-orient in response to tissue tension, but not simply along tension line.
 * In unconstraint tissue, tissue will shrink in the direction perpendicular to tension, and collagen fibers will align along tension.
 * In constraint tissue, tissue cannot shrink, and collagen fiber will align perpendicular to tension to perform "strain avlidance". 
 */
void ECM::collagen_orientation_under_tension(const Mat_DP& tissue_displacement_x,
                                             const Mat_DP& tissue_displacement_y, double time_step){
    DP kappa1 = 5, kappa2 = 1;
    for(int i=0; i<ECM_ystep/10; i++){//from bottom up
        for(int j=0; j<ECM_xstep/10; j++){//from left to right
            Mat_DP F(2,2); //deformation gradient
            F[0][0] = tissue_displacement_x[i*10][(j+1)*10>=ECM_xstep ? j*10 : (j+1)*10] - tissue_displacement_x[i*10][(j-1)*10<0 ? 0 : (j-1)*10];
            F[0][1] = tissue_displacement_x[(i+1)*10>=ECM_ystep ? i*10 : (i+1)*10][j*10] - tissue_displacement_x[(i-1)*10<0 ? 0 : (i-1)*10][j*10];
            F[1][0] = tissue_displacement_y[i*10][(j+1)*10>=ECM_xstep ? j*10 : (j+1)*10] - tissue_displacement_y[i*10][(j-1)*10<0 ? 0 : (j-1)*10];
            F[1][1] = tissue_displacement_y[(i+1)*10>=ECM_ystep ? i*10 : (i+1)*10][j*10] - tissue_displacement_y[(i-1)*10<0 ? 0 : (i-1)*10][j*10];
            F[0][0] /= 20;
            F[0][1] /= 20;
            F[1][0] /= 20;
            F[1][1] /= 20;
            F[0][0] += 1;
            F[1][1] += 1;
            //tissue is stretched if stretch > 1.0, compressed if stretch is < 1.0
            //reorientation only happens when stretch > 1.2
            double stretch1, stretch2, tensiontheta;
            std::tie(stretch1, stretch2, tensiontheta) = calculate_principal_strain(F);
            stretch_history[i][j] *= 0.95;
            if(stretch1 > 1.2){
                stretch_history[i][j] += 1.0; 
            }
            cout << "stretch1 = " << stretch1 << endl;
            if(stretch_history[i][j] > 10.0){
                /*
                if (stretch2 > 1.0) {//constraint tissue, collagen align perpendicular to principal strain line
                    tensiontheta += M_PI * 0.5;
                    if (tensiontheta > M_PI) tensiontheta -= M_PI;
                }
                */
                cout << "tensiontheta = " << tensiontheta << endl;
                DP increase; 
                if(abs(tensiontheta-collagen[i][j])<=M_PI/2)
                    increase = time_step * kappa2 * sqrt(stretch1) * sin(tensiontheta-collagen[i][j]);
                else if(tensiontheta>M_PI/2)
                    increase = time_step * kappa2 * sqrt(stretch1) * sin(tensiontheta-M_PI-collagen[i][j]);
                else increase = time_step * kappa2 * sqrt(stretch1) * sin(tensiontheta+M_PI-collagen[i][j]);

                collagen[i][j] += increase;
                cout <<"increase="<<increase<<endl;
                while(collagen[i][j] >= M_PI) {
                    collagen[i][j] -= M_PI;
                }
                while(collagen[i][j] < 0) {
                    collagen[i][j] += M_PI;
                }
            } 
        }
    }
    return;                                                        
}

void ECM::output_collagen(Flist& fibroblastList){
    const DP dl=1.0;
    DP x,y,xnew,ynew;
    cout << "xstep = " << ECM_xstep << "ystep = " << ECM_ystep <<endl;
    Mat_DP vf(ECM_ystep,ECM_xstep); 
    
    srand(time(NULL));
    for(int i=0;i<ECM_ystep;i++){
        for(int j=0;j<ECM_xstep;j++){            
            vf[i][j] = 0;                      
        }       
    }
                
    for(int j=0; j<ECM_ystep/10; j++){ 
        for(int l=0; l<ECM_xstep/10; l++){
            y = j*10;
            x = l*10;
            for(int k=0; k<80; k++){
                //printf("x = %f, y = %f\n",x,y);                
                xnew = x + dl*cos(collagen[int(y/10)][int(x/10)]);
                if(xnew < 0)xnew = 0;
                if(xnew > ECM_xstep-1)xnew = ECM_xstep-1;
                ynew = y + dl*sin(collagen[int(y/10)][int(x/10)]);
                if(ynew < 0)ynew= 0;
                if(ynew > ECM_ystep-1)ynew = ECM_ystep-1;
                x = xnew;
                y = ynew;
                vf[(int)y][(int)x] = collagen_density[(int)(y/5)][(int)(x/5)];                
            }
            y = j*10;
            x = l*10;    
            for(int k=0; k<20; k++){
                xnew = x - dl*cos(collagen[int(y/10)][int(x/10)]); 
                if(xnew < 0)xnew = 0;
                if(xnew > ECM_xstep-1)xnew = ECM_xstep-1;
                ynew = y - dl*sin(collagen[int(y/10)][int(x/10)]);
                if(ynew < 0)ynew = 0;
                if(ynew > ECM_ystep-1)ynew = ECM_ystep-1;
                x = xnew;
                y = ynew;
                vf[(int)y][(int)x] = collagen_density[int(y/5)][int(x/5)];
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
            for(int j=(int)curPtr->yy-2; j<=(int)curPtr->yy+2; j++){
                for(int l=(int)curPtr->xx-2; l<=(int)curPtr->xx+2; l++){
                    cellmatrix[j][l] = 1;
                }
            } 
        }
        //cout << (int)curPtr->xx << (int)curPtr->yy << endl;
    }
    
    
    //************output to picture*****************************************
    //Ö»Êä³öÒ»·ùÍ¼ 
    static int out_i = 0;
    char file_name[21];
    sprintf(file_name, "output/coll%05d.BMP", out_i++);
    BMP::output_BMP2(file_name, 14, vf, cellmatrix, ECM_xstep, ECM_ystep);
   
    return;
}
    
