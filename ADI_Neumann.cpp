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
 
 
char  file_name[15]="self00000.BMP";
const DP dx = 1;//Unit: um
const DP dy = 1;//Unit: um
const DP dt = 0.01;//Unit: hr


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
 
int ADI(DP tlength, DP D, Mat_DP& PDGF, bool output)
{
    FILE* dataFile;
    int j,l,k;
    static int out_i = 0;
    
    DP alpha = D/2*dt/pow(dx,2);
    
    int xstep, ystep, tstep;    
    ystep = PDGF.nrows();//?????????????
    xstep = PDGF.ncols();//?????????????
    tstep = (int)(tlength/dt);
    //for (k=0;k<n;k++) cout << setw(12) << u[k];
    cout << "xstep" << xstep << "ystep" << ystep << "tstep" << tstep << "D" << D << endl;
    
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
    u = PDGF;        
    for(k=0;k<tstep;k++){
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
        u = unext;
    }

/*    
    if(output == 1){
        int jj;
        file_name[8] = 48+out_i%10;        file_name[7] = 48+out_i/10%10;
        file_name[6] = 48+out_i/10/10%10;  file_name[5] = 48+out_i/10/10/10%10;
        fp_BMP = fopen(file_name,"w") ;
        Prepare_BMP_Format(xstep,ystep) ;
	    for(j=0;j<ystep;j++) {
		    for(l=0;l<xstep;l++) {
                jj = int((10*u[j][l] - Minimal)/Interval); 
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
		        if(jj < 0)    fwrite(&b_0,  1,3,fp_BMP);
		        if(jj > 15)   fwrite(&b_16, 1,3,fp_BMP);
            }
	    }
        out_i++;
	    fclose(fp_BMP); 
    }
*/ 
            
    PDGF = u;                    
    return EXIT_SUCCESS;
}
