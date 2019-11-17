//Inplementation for reaction-diffution PDE with the form 
//du/dt-D*delta(u)=0 
//Neumann boundary condition == no flux boundary condition
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>//pow
#include "nr.h"
//#include "BMP.h"
using namespace std;  
 
 
const DP dx = 1;//Unit: um

void NR::tridag(Vec_I_DP &a, Vec_I_DP &b, Vec_I_DP &c, Vec_I_DP &r, Vec_O_DP &u)
{     
     int j;
     DP bet;
     //for(int k=0;k<10;k++)cout << r[k] << " ";
     //cout << endl;
     //cout << "tridag" << endl;
     
     int n=a.size();
     Vec_DP gam(n);
     if(b[0] == 0.0)nrerror("Error 1 in tridag");
     bet = b[0];
     u[0] = r[0]/bet; 
     for(j=1;j<n;j++){
         gam[j] = c[j-1]/bet;
         bet = b[j] - a[j]*gam[j];
         if(bet == 0.0) nrerror("Error 2 in tridag");
         u[j] = (r[j]-a[j]*u[j-1])/bet; 
     }
     for(j=(n-2);j>=0;j--)
         u[j] -= gam[j+1]*u[j+1];

}
 
int ADI(DP time_step, Mat_DP& GF, DP D){
    int j,l,k;
    
    DP alpha = D/2*time_step/pow(dx,2);
    
    int xstep, ystep, nSteps = 1;    
    ystep = GF.nrows();//?????????????
    xstep = GF.ncols();//?????????????
    //for (k=0;k<n;k++) cout << setw(12) << u[k];
    
    Mat_DP u(ystep,xstep),uhalf(ystep,xstep),unext(ystep,xstep);
    Vec_DP uhalftemp(xstep),unexttemp(xstep);
    Vec_DP diagx(xstep),superdx(xstep),subdx(xstep),rhsx(xstep);
    Vec_DP diagy(ystep),superdy(ystep),subdy(ystep),rhsy(ystep);
    
//*******************Neumann boundary condition*********************************   
    for(j=0;j<xstep;j++)diagx[j] = 1+2*alpha;
    superdx[0] = -2*alpha;
    for(j=1;j<xstep;j++)superdx[j] = -alpha;
    subdx[xstep-1] = -2*alpha;
    for(j=0;j<xstep-1;j++)subdx[j] = -alpha;
    
    for(l=0;l<ystep;l++)diagy[l] = 1+2*alpha;
    superdy[0] = -2*alpha;
    for(l=1;l<ystep;l++)superdy[l] = -alpha;
    subdy[ystep-1] = -2*alpha;
    for(l=0;l<ystep-1;l++)subdy[l] = -alpha;
//*******************Neumann boundary condition*********************************

//**********Initial condition**************************
    u = GF;        
    for(k=0;k<nSteps;k++){
        //***********µÚ¶þ¸ö°ë²½*************************************************
        //******
        //@@@@@@
        //******
        //******
        //******
        // - alpha*unext[j][l-1] + (1+2*alpha)*unext[j][l] - alpha*unext[j][l+1] 
        //     = alpha*uhalf[j-1][l] + (1-2*alpha)*uhalf[j][l] + alpha*uhalf[j+1][l]; 
        for(j=0;j<ystep;j++){                   
            for(l=0;l<xstep;l++){
                if(j == 0)rhsx[l] = (1-2*alpha)*u[j][l] + 2*alpha*u[j+1][l];
                else if(j == ystep-1)rhsx[l] = 2*alpha*u[j-1][l] + (1-2*alpha)*u[j][l];
                else rhsx[l] = alpha*u[j-1][l] + (1-2*alpha)*u[j][l] + alpha*u[j+1][l];
            }
            NR::tridag(subdx,diagx,superdx,rhsx,uhalftemp);
            for(l=0;l<xstep;l++)uhalf[j][l] = uhalftemp[l];
        }
        cout << "ok1" << endl;
       
        //***********µÚÒ»¸ö°ë²½*************************************************
        //*@****
        //*@****
        //*@****
        //*@****
        //*@****
        // - alpha*uhalf[j-1][l] + (1+2*alpha)*uhalf[j][l] - alpha*uhalf[j+1][l] 
        //     = alpha*u[j][l-1] + (1-2*alpha)*u[j][l] + alpha*u[j][l+1]; 
        for(l=0;l<xstep;l++){                                          
            for(j=0;j<ystep;j++){
                if(l == 0)rhsy[j] = (1-2*alpha)*uhalf[j][l] + 2*alpha*uhalf[j][l+1];
                else if(l == xstep-1)rhsy[j] = 2*alpha*uhalf[j][l-1] + (1-2*alpha)*uhalf[j][l]; 
                else rhsy[j] = alpha*uhalf[j][l-1] + (1-2*alpha)*uhalf[j][l] + alpha*uhalf[j][l+1];
            }
            NR::tridag(subdy,diagy,superdy,rhsy,unexttemp);
            for(j=0;j<ystep;j++)unext[j][l] = unexttemp[j];            
        }
        cout << "ok2" << endl;
        u = unext;
    }

            
    GF = u;                    
    return EXIT_SUCCESS;
}
