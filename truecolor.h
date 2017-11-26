#ifndef _TRUE_COLOR_H_
#define _TRUE_COLOR_H_

void Get_True_Color(float f, float max, float min,unsigned char& true_color1, unsigned char& true_color2, unsigned char& true_color3)
{
    int x;
    int r,g,b;
    int detal;
    x= (int) (  ((f-min)*1276)     /(max -min) );

    //if(1276 < num_color)
    int num_color = 1276;
    
    detal=1276/num_color;
    x= x- x%detal;
    
    cout<<"x="<<x<<endl;
    
    if(x<0)x=0;
    if(x>1275)x=1275;
    if(x<=255){
        r=255-x;
        g=0;
        b=255;
        //return RGB(r,g,b);
        true_color1 = b; true_color2 = g; true_color3 = r; 
    }
    
    if(x>255 && x<=510){
        r=0;
        g=x - 255;
        b=255;
        //return RGB(r,g,b);
        true_color1 = b; true_color2 = g; true_color3 = r; 
    }
    
    if(x>510 && x<=765){
        r=0;
        g=255;
        b=765-x;
        //return RGB(r,g,b);
        true_color1 = b; true_color2 = g; true_color3 = r; 
    }
    
    if(x>765 && x<=1020){
        r=x-765;
        g=255;
        b=0;
        //return RGB(r,g,b);
        true_color1 = b; true_color2 = g; true_color3 = r; 
    }
    
    if(x>1020){
        r=255;
        g=1275-x;
        b=0;
        //return RGB(r,g,b);
        true_color1 = b; true_color2 = g; true_color3 = r; 
    } 
    return;
}

#endif
