#include <string>
#include <cmath>
#include <iostream>
#include <map>
#include "BMP.h" 
using namespace std; 

namespace BMP{
          
unsigned char b0[3], b1[3], b2[3], b3[3];
unsigned char b4[3], b5[3], b6[3], b7[3];
unsigned char b8[3], b9[3], b10[3],b11[3];
unsigned char b12[3],b13[3],b14[3],b15[3];
unsigned char b_0[3],b_16[3];
unsigned char b_yellow[3], b_green[3], b_red[3];

unsigned char b_green1[3], b_green2[3], b_green3[3];	
unsigned char b_green4[3], b_green5[3];

unsigned char b_red1[3], b_red2[3], b_red3[3];	
unsigned char b_red4[3], b_red5[3];	

float Maximal=  2.0 ;
float Minimal=  0.0 ;
float Interval=(Maximal-Minimal)/16.0;

inline void Prepare_BMP_Format(std::ofstream &fp_BMP, int width, int height) {

    int i;
//The BMP format requires that N_x should be an integer times of 4

    struct BMP_54_Bytes {                        //  offset   size in unit of byte
        unsigned char        Bmp_id[2] ;         // 0        2 (bytes)
        unsigned long        Bmp_size ;          // 2        4 
        unsigned long        Researved;          // 6        4
        unsigned long        Dada_offset;        // 10       4
        unsigned long        Size_header_inf ;   // 14       4   (always = 40) ;
        unsigned long        Width;              // 18       4
        unsigned long        Height;             // 22       4
        unsigned short int   Color_planes   ;    // 26       2
        unsigned short int   Bits_pixel    ;     // 28       2
        unsigned long        Compress_type ;     // 30       4
        unsigned long        image_size ;        // 34       4
        unsigned long        X_pixels_meter ;    // 38       4
        unsigned long        Y_pixels_meter ;    // 42       4
        unsigned long        Color_used     ;    // 46       4
        unsigned long        Color_important;    // 50       4
    } pic ;

    pic.Bmp_id[0] =         'B';
    pic.Bmp_id[1] =         'M';         
    pic.Bmp_size  =         width * height*3 + 54;   // !!!!!!
    pic.Researved =         0;
    pic.Dada_offset =       54;
	  
    pic.Size_header_inf=    40;
    pic.Width=              width;             // !!!!!!  width/4 must be an integer
    pic.Height=             height;            // !!!!!!
    pic.Color_planes=       1;
    pic.Bits_pixel=         24;
    pic.Compress_type=      0;                 // 0 for no-compress
    pic.image_size=         0;
    pic.X_pixels_meter=     11808*2 ;          // 600*600 resolution 
    pic.Y_pixels_meter=     11808*2 ;    
    pic.Color_used =        0;
    pic.Color_important =   0;

    cout << "BMP ok0" << endl;
/* write the header information to BMP file  */
    fp_BMP.write((const char*)&pic.Bmp_id,2);
    fp_BMP.write((const char*)&pic.Bmp_size,4);
    fp_BMP.write((const char*)&pic.Researved,4);
    fp_BMP.write((const char*)&pic.Dada_offset,4);

    cout << "BMP ok1" << endl;
    fp_BMP.write((const char*)&pic.Size_header_inf,4);
    fp_BMP.write((const char*)&pic.Width,4);
    fp_BMP.write((const char*)&pic.Height,4);
    fp_BMP.write((const char*)&pic.Color_planes,2);
    fp_BMP.write((const char*)&pic.Bits_pixel,2); 
    fp_BMP.write((const char*)&pic.Compress_type,4);
    fp_BMP.write((const char*)&pic.image_size,4);
    fp_BMP.write((const char*)&pic.X_pixels_meter,4);
    fp_BMP.write((const char*)&pic.Y_pixels_meter,4);
    fp_BMP.write((const char*)&pic.Color_used,4);
    fp_BMP.write((const char*)&pic.Color_important,4);

/* prepare grayscale map: full range: 0 -- 255  
	                          used 16 grey scales         */
    for(i=0;i<3;i++) {
        b0[i]=15;       b1[i]=31;   	b2[i]=47;      b3[i]=63;
        b4[i]=79;       b5[i]=95;   	b6[i]=111;     b7[i]=127;
	b8[i]=143;      b9[i]=159;  	b10[i]=175;    b11[i]=191;
	b12[i]=207;     b13[i]=223;  	b14[i]=239;    b15[i]=255;
	b_0[i]  = 15 ;  // for scale levels below b0 
	b_16[i] = 255 ; // for scale levels above b15
    }

    for(i=0;i<3;i++) {
        b_green1[i] = i==1 ? 255 : 250;
        b_green2[i] = i==1 ? 255 : 200;
        b_green3[i] = i==1 ? 255 : 150;
        b_green4[i] = i==1 ? 255 : 100;
        b_green5[i] = i==1 ? 255 : 50;
    }

    for(i=0;i<3;i++) {
        b_red1[i] = i==2 ? 255 : 250;
        b_red2[i] = i==2 ? 255 : 200;
        b_red3[i] = i==2 ? 255 : 150;
        b_red4[i] = i==2 ? 255 : 100;
        b_red5[i] = i==2 ? 255 : 50;
    }
 
    b_yellow[0] = 0, b_yellow[1] = 255, b_yellow[2] = 255;
    b_green[0] = 0, b_green[1] = 255, b_green[2] = 0;
    b_red[0] = 255, b_red[1] = 0, b_red[2] = 0;
}

void output_BMP(char* file_name, int n, const Mat_DP& data, int xstep, int ystep){
    int jj;
    std::ofstream fp_BMP{ file_name, std::ios_base::binary };
    Prepare_BMP_Format(fp_BMP, xstep,ystep);
    for(int i=0;i<ystep;i++){
        for(int j=0;j<xstep;j++){
            jj = int( (data[i][j] - Minimal)/Interval );
            if(jj==0)   fp_BMP.write((const char*)&b0, 3);
            if(jj==1)   fp_BMP.write((const char*)&b1, 3);
            if(jj==2)   fp_BMP.write((const char*)&b2, 3);
            if(jj==3)   fp_BMP.write((const char*)&b3, 3);
            if(jj==4)   fp_BMP.write((const char*)&b4, 3);
            if(jj==5)   fp_BMP.write((const char*)&b5, 3);
            if(jj==6)   fp_BMP.write((const char*)&b6, 3);
            if(jj==7)   fp_BMP.write((const char*)&b7, 3);
            if(jj==8)   fp_BMP.write((const char*)&b8, 3);
            if(jj==9)   fp_BMP.write((const char*)&b9, 3);
            if(jj==10)   fp_BMP.write((const char*)&b10, 3);
            if(jj==11)   fp_BMP.write((const char*)&b11, 3);
            if(jj==12)   fp_BMP.write((const char*)&b12, 3);
            if(jj==13)   fp_BMP.write((const char*)&b13, 3);
            if(jj==14)   fp_BMP.write((const char*)&b14, 3);
            if(jj==15)   fp_BMP.write((const char*)&b15, 3);
            if(jj < 0 )  fp_BMP.write((const char*)&b_0, 3);  
            if(jj > 15)  fp_BMP.write((const char*)&b_16,3); 
        }
    }
    fp_BMP.close();
    return;
}

void output_BMP2(char* file_name,int n, const Mat_DP& data, const Mat_DP& cellmatrix, int xstep,int ystep){
    int jj;
    std::ofstream fp_BMP{ file_name, std::ios_base::binary };
    Prepare_BMP_Format(fp_BMP, xstep,ystep);
    for(int i=0;i<ystep;i++){
        for(int j=0;j<xstep;j++){
            jj = data[i][j] >= 0.0 ? int( (data[i][j] - Minimal)/Interval ) : -1;
            if (data[i][j] > 0.0 and data[i][j] < Interval) jj = 1;
            if(cellmatrix[i][j] == 1) fp_BMP.write((const char*)&b_yellow, 3);
            else{
                if(jj==0)   fp_BMP.write((const char*)&b0, 3);
                if(jj==1)   fp_BMP.write((const char*)&b1, 3);
                if(jj==2)   fp_BMP.write((const char*)&b2, 3);
                if(jj==3)   fp_BMP.write((const char*)&b3, 3);
                if(jj==4)   fp_BMP.write((const char*)&b4, 3);
                if(jj==5)   fp_BMP.write((const char*)&b5, 3);
                if(jj==6)   fp_BMP.write((const char*)&b6, 3);
                if(jj==7)   fp_BMP.write((const char*)&b7, 3);
                if(jj==8)   fp_BMP.write((const char*)&b8, 3);
                if(jj==9)   fp_BMP.write((const char*)&b9, 3);
                if(jj==10)   fp_BMP.write((const char*)&b10, 3);
                if(jj==11)   fp_BMP.write((const char*)&b11, 3);
                if(jj==12)   fp_BMP.write((const char*)&b12, 3);
                if(jj==13)   fp_BMP.write((const char*)&b13, 3);
                if(jj==14)   fp_BMP.write((const char*)&b14, 3);
                if(jj==15)   fp_BMP.write((const char*)&b15, 3);
                if(jj==-1)   fp_BMP.write((const char*)&b_green1, 3);
                if(jj < 0 )  fp_BMP.write((const char*)&b_0, 3);  
                if(jj > 15)  fp_BMP.write((const char*)&b_16, 3); 
            }
        }
    }
    fp_BMP.close();
    return;
}

/*
void output_BMP3(char* file_name,int n, const map<pair<int, int>, int>& data, int xstep,int ystep){
    int jj;
    FILE *fp_BMP = fopen(file_name, "w");
    Prepare_BMP_Format(fp_BMP, xstep,ystep);
    for(int i=0;i<ystep;i++){
        for(int j=0;j<xstep;j++){
            std::pair<int, int> key = std::make_pair(i, j);
            if(data.find(key) != data.end()) {
                jj = data.at(key);
                if(jj==0)   fwrite(&b_green1, 1,3,fp_BMP);
                if(jj==10)   fwrite(&b_green2, 1,3,fp_BMP);
                if(jj==20)   fwrite(&b_green3, 1,3,fp_BMP);
                if(jj==30)   fwrite(&b_green4, 1,3,fp_BMP);
                if(jj==40)   fwrite(&b_green5, 1,3,fp_BMP);
                if(jj==50)   fwrite(&b_red1, 1,3,fp_BMP);
                if(jj==60)   fwrite(&b_red2, 1,3,fp_BMP);
                if(jj==70)   fwrite(&b_red3, 1,3,fp_BMP);
                if(jj==80)   fwrite(&b_red4, 1,3,fp_BMP);
                if(jj==90)   fwrite(&b_red5, 1,3,fp_BMP);
                if(jj >= 100)  fwrite(&b_red5, 1,3,fp_BMP); 
            } else {
                fwrite(&b_0, 1,3,fp_BMP);
            } 
        }
    }
    fclose(fp_BMP);
    return;
} */
}//end of namespace BMP
