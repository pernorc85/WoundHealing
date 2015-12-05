#ifndef _FLIST_H_
#define _FLIST_H_

typedef struct Fibroblast Fcell;
typedef Fcell* FcellPtr;
struct Fibroblast{
    DP xx,yy,zz; 
    //DP theta;//polar coordinate: angle in x-y plane [0,360]
    //DP phi;  //polar coordinate: angle away from x-y plane [-90,90]
    DP angle_xx, angle_yy, angle_zz;
    //angle_x = cos(phi)*cos(theta)
    //angle_y = cos(phi)*sin(theta)
    //angle_z = sin(phi)
    DP speed;
    FcellPtr next;
};

class Flist
{
    public:
        Flist(); 
        Flist(int FNinit);
        
        void buildFlist(FcellPtr* sPtr); 
        void Flist_move(ReactionDiffusionEquation& PDGF_FEM,
                        ECM& extraCellularMatrix); 
        void Fcell_move(DP gradx, DP grady, DP gradz, DP fdensity, DP cdensity, 
                        Mat3D_DP& collagen_angle_xx,
                        Mat3D_DP& collagen_angle_yy,
                        Mat3D_DP& collagen_angle_zz); 
        void output_fibroblast();

    protected:
        FcellPtr Froot;
        FcellPtr curPtr;
        int FN;
        DP speed;
    //ÕâÀï²»ÄÜÓÃFcellPtr *sPtr; *sPtr = NULL;ÕâÑùÐ´Ö®ºósystem("PAUSE");È«²¿Ê§Ð§
    //¶øÇÒbuildFlistÒ²Ê§Ð§ 
    friend class ECM;
};

#endif
