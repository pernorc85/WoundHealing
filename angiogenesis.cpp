//angiogenesis.cpp
//    VERSION: (1 August 2010)
//    PURPOSE: Source file for the wound healing angiogenesis model.
//             Include all functions for the angiogenesis process:
//             Tip formation, tip migration, branching, anastamosis.
//             Also includes functions for numeric calculation related to angiogenesis
//             and output routine for vasculature morphology.
//    INPUT FILES:
//             None.
//    OUTPUT FILES:
//             "vasc*****.BMP" 

#include <cmath>
#include "angiogenesis.h"
#include "record_vasculature_network.h"
#include "BMP.h"
#include "ECM.h"


int Elist::branch_id=1;

DP speed_function(DP grad);
DP angle_substract(DP theta1,DP theta2);
void feedback_from_oxygen_to_VEGF(Mat_DP& vasculature_density,Mat_DP& VEGF);
void regression(Vec_DP xlist,int xcount,Vec_DP ylist,int ycount,double* b,double* xycor); 



//################################################################################
//## buildElist:                                                                 #
//##                                                                             #
//##         This function build a list of endothelial cells(vascular tips)      #
//##         placed at the wound boundary.                                       #
//##                                                                             #
//##         Program currently set to place 12 tips evenly spaced at the round   #
//##         shaped wound boundary with a direction inward.                      #
//################################################################################
void Elist::buildElist(unsigned int ENinit, double wound_radius, EdgePtr* edgeList)
{
    DP x,y,theta;
    EcellPtr newPtr=NULL,prePtr=NULL,curPtr=NULL;
    
    for(int i=0;i<ENinit;i++){
        x = mXstep/2 + (wound_radius + 100)*cos(M_PI/6*i);
        y = mYstep/2 + (wound_radius + 100)*sin(M_PI/6*i);
        if(i<=6)theta = M_PI/6*(i+6);
        else theta = M_PI/6*(i-6);
        cout << "x=" << x << " y=" << y << " theta=" << theta << endl; 
        {/*                
            EdgePtr new_edge = (EdgePtr)malloc(sizeof(Edge));
            new_edge->ID = running_edge_ID;
            cout<<"creating edge "<<running_edge_ID<<endl;            
            new_edge->left_node = NULL;
            //will be left as NULL forever.
            //In the final matrix for pressure and flow calculation, 
            //these edges without any left_node will be implemented specifically.
            new_edge->right_node = NULL;
            new_edge->xn = 0;
            new_edge->yn = 0;
            add_to_edgeList(edgeList,new_edge,running_edge_ID);
            running_edge_ID++;
        */}
        newPtr = new Endocyte(x, y, theta, branch_id++, 0);
        newPtr->curEdgePtr = NULL;
        newPtr->next = NULL;
        prePtr = NULL;
        curPtr = *Eroot;
        while(curPtr != NULL){
            prePtr = curPtr;
	    curPtr = curPtr->next;
        }      
        if(prePtr == NULL){
            newPtr->next =*Eroot;
            *Eroot = newPtr;
        }else{
            prePtr->next = newPtr;
            newPtr->next = curPtr;
        }          
    }      
    return;
}

void Elist::Elist_move(Mat_DP &gradx_VEGF, Mat_DP &grady_VEGF, ECM &extraCellularMatrix, DP time_step) {
    EcellPtr ecurPtr = *Eroot;
    Mat_DP &collagen = extraCellularMatrix.collagen;
    Mat_DP &fibronectin_density = extraCellularMatrix.fibronectin_density;
    Mat_DP &collagen_density = extraCellularMatrix.collagen_density;
    while(ecurPtr != NULL){
        DP gradx = gradx_VEGF[(int)ecurPtr->yy][(int)ecurPtr->xx];
        DP grady = grady_VEGF[(int)ecurPtr->yy][(int)ecurPtr->xx];
        DP collagen_theta = collagen[(int)(ecurPtr->yy/10)][(int)(ecurPtr->xx/10)];
        DP fdensity = fibronectin_density[(int)(ecurPtr->yy/5)][(int)(ecurPtr->xx/5)];
        DP cdensity = collagen_density[(int)(ecurPtr->yy/5)][(int)(ecurPtr->xx/5)];
        Ecell_move(ecurPtr, gradx, grady, cdensity, collagen_theta, time_step);
        ecurPtr = ecurPtr->next;
    } 
}

//################################################################################
//## Ecell_move:                                                                 #
//##                                                                             #
//##         This function implement the migration of endothelial cells          #
//##         following the guicance of chemokine and collagen.                   #
//##         Tips that are subject to branch will be marked.                     #
//################################################################################
void Elist::Ecell_move(EcellPtr cellPtr, DP gradx, DP grady, DP cdensity, DP collagen_temp, DP time_step)
{ 
    {/*
        if((cellPtr->curEdgePtr)->xn<200){
            (cellPtr->curEdgePtr)->trace_x[cellPtr->curEdgePtr->xn] = (int)cellPtr->xx;
            (cellPtr->curEdgePtr)->xn++;
        }
        if((cellPtr->curEdgePtr)->yn<200){
            (cellPtr->curEdgePtr)->trace_y[cellPtr->curEdgePtr->yn] = (int)cellPtr->yy;
            (cellPtr->curEdgePtr)->yn++;
        }
    */}
    Mat_DP& cellmatrix = *cellmatrix_ptr; 
    DP rou1 = 0.3, rou2 = 0.5;
    DP collagen_taken = collagen_temp;
    DP gradtheta = -1, gradtheta_taken; 
    
    DP grad=min(1.0,1000*sqrt(pow(gradx,2)+pow(grady,2)));
    //rou1 = rou1 * grad;
    //rou2 = rou2 * cdensity;
    DP speed_taken = 10;
    
    if(gradx != 0){
        if(gradx > 0 && grady >= 0)gradtheta = atan(grady/gradx);//1st phase
        else if(gradx < 0 && grady >= 0)gradtheta = M_PI + atan(grady/gradx);//2nd phase
        else if(gradx < 0 && grady <= 0)gradtheta = M_PI + atan(grady/gradx);//3nd phase
        else if(gradx > 0 && grady <= 0)gradtheta = 2*M_PI + atan(grady/gradx);//4th phase
        //if(gradtheta < 0 || gradtheta > 2*PI) cout<< "error" << endl;
    } else {
        if(grady > 0)gradtheta =  M_PI/2;
        else if(grady < 0)gradtheta = M_PI*3/2;
        else gradtheta = -1;
    }
 
    gradtheta_taken = gradtheta;
    while(abs(cellPtr->theta - gradtheta_taken) > M_PI) {
        if (gradtheta_taken < cellPtr->theta)
            gradtheta_taken += 2*M_PI;
        else if (gradtheta_taken > cellPtr->theta)
            gradtheta_taken -= 2*M_PI;
    }

    while(abs(cellPtr->theta - collagen_taken) > M_PI/2) {
        if(collagen_taken > cellPtr->theta){
            collagen_taken -= M_PI;
        }else{
            collagen_taken += M_PI;
        }
    }
    

    //make a decision if to branch 
    int branch = 0;
    if(//(collagen_taken - cellPtr->theta) * (gradtheta_taken - cellPtr->theta) < 0 and
        abs(collagen_taken - gradtheta_taken) > M_PI/2 and
        cellPtr->branching == 0)branch = 1;
    //make sure there is enough space to branch
    //==there is only one trail  
    int i,j;
    DP ydist, xdist, theta;
    if(branch == 1){
        for(i=max(0, (int)cellPtr->yy-8); i<min(mYstep-1, (int)cellPtr->yy+8); i++){
            for(j=max(0, (int)cellPtr->xx-8); j<min(mXstep-1, (int)cellPtr->xx+8); j++){
                if((int)(cellmatrix[i][j])%10 == 1 || (int)(cellmatrix[i][j])%10 == 3){
                    ydist = i - (int)cellPtr->yy; xdist = j - (int)cellPtr->xx;
                    if(xdist > 0 && ydist >= 0)theta = atan(ydist/xdist);//1st phase
                    else if(xdist < 0 && ydist >= 0)theta = M_PI + atan(ydist/xdist);//2nd phase
                    else if(xdist < 0 && ydist <= 0)theta = M_PI + atan(ydist/xdist);//3nd phase
                    else if(xdist > 0 && ydist <= 0)theta = 2*M_PI + atan(ydist/xdist);//4th phase
                    else{
                        if(ydist > 0)theta =  M_PI/2;
                        else if(ydist < 0)theta = M_PI*3/2;
                        else theta = -1;
                    }
                
                    if(theta != -1 && abs(angle_substract(theta, cellPtr->theta))< M_PI*2/3 )branch = 0;
                }     
            }
        }
    } 
    //--------------------------------------------------------------------
     
    if(branch == 1){
        cellPtr->branching = 25;
        cout << "branching marked " <<endl;
    } else {
        if(gradtheta >= 0){
            cellPtr->theta = (1 - rou1 - rou2)*cellPtr->theta + rou1*gradtheta_taken + rou2*collagen_taken;
        } else if(gradtheta == -1) {
            cellPtr->theta = (1 - rou2)*cellPtr->theta + rou2 * collagen_taken;
        }

        if(cellPtr->theta > 2*M_PI) {
            cellPtr->theta -= 2*M_PI;
        }
        if(cellPtr->theta < 0) {
            cellPtr->theta += 2*M_PI;
        }
        cellPtr->xx += speed_taken * time_step * cos(cellPtr->theta);
        cellPtr->yy += speed_taken * time_step * sin(cellPtr->theta);
        
        if(cellPtr->xx < 0)cellPtr->xx = 0;
        if(cellPtr->yy < 0)cellPtr->yy = 0;
        if(cellPtr->xx >= mXstep)cellPtr->xx = mXstep-1;
        if(cellPtr->yy >= mYstep)cellPtr->yy = mYstep-1;  
    }
    return;          
}

//################################################################################
//## Ecell_branch:                                                               #
//##                                                                             #
//##         This function implements the branching of a vascular tip            #
//##         After branching, one daughter tip will occupy the storage space of  #
//##         the original tip, the other daughter tip will be created as a new   #
//##         Ecell.                                                              #
//################################################################################
void Elist::Elist_branch(EdgePtr* edgeList, NodePtr* nodeList, ECM& extraCellularMatrix, 
                         Mat_DP& gradx_VEGF, Mat_DP& grady_VEGF, DP time_step)
{ 
    Mat_DP &collagen = extraCellularMatrix.collagen;
    DP collagen_taken, gradx, grady, gradtheta = -1;
    EcellPtr newPtr, prePtr, curPtr;
    curPtr = *Eroot; 
    int gridy,gridx;
    
    while(curPtr != NULL){
        if(curPtr->branching == 25){  
            printf("vascular tip branch.\n");
            collagen_taken = collagen[(int)(curPtr->yy/10)][(int)(curPtr->xx/10)]; 
            //be maticulous here-----------------------------------------------
            while(abs(curPtr->theta - collagen_taken) > M_PI/2) {
               if(collagen_taken > curPtr->theta){
                   collagen_taken -= M_PI;
               }else{
                   collagen_taken += M_PI;
               }
            }
            //-----------------------------------------------------------------
            gradx = gradx_VEGF[(int)curPtr->yy][(int)curPtr->xx];
            grady = grady_VEGF[(int)curPtr->yy][(int)curPtr->xx];
            if(gradx != 0){
                if(gradx > 0 && grady >= 0)gradtheta = atan(grady/gradx);//1st phase
                else if(gradx < 0 && grady >= 0)gradtheta = M_PI + atan(grady/gradx);//2nd phase
                else if(gradx < 0 && grady <= 0)gradtheta = M_PI + atan(grady/gradx);//3nd phase
                else if(gradx > 0 && grady <= 0)gradtheta = 2*M_PI + atan(grady/gradx);//4th phase
                //cout << "gradx=" << gradx << "grady=" << grady << "gradtheta=" << gradtheta << endl; 
            } else {
                if(grady > 0)gradtheta =  M_PI/2;
                else if(grady < 0)gradtheta = M_PI*3/2;
                else gradtheta = -1;
            }
           
            curPtr->xx += speed * time_step *cos(collagen_taken);
            curPtr->yy += speed * time_step *sin(collagen_taken); 
            if(curPtr->xx < 0)curPtr->xx = 0;
            if(curPtr->yy < 0)curPtr->yy = 0;
            if(curPtr->xx >= mXstep)curPtr->xx = mXstep-1;
            if(curPtr->yy >= mYstep)curPtr->yy = mYstep-1; 
            curPtr->theta = collagen_taken;
            curPtr->branching--;  
            
            DP xx = curPtr->xx + speed * time_step * cos(gradtheta);
            DP yy = curPtr->yy + speed * time_step * sin(gradtheta); 
            if(xx < 0)xx = 0;
            if(yy < 0)yy = 0;
            if(xx >= mXstep)xx = mXstep-1;
            if(yy >= mYstep)yy = mYstep-1; 
            newPtr = new Endocyte(xx, yy, gradtheta, branch_id, 24);
            branch_id++;
            newPtr->next = NULL; 
            
 
            {/*
                //***************create one new node and two new edges**************************
                printf("Size of node %d\n",sizeof(Node));
                NodePtr new_node = (NodePtr)malloc(sizeof(Node));
                new_node->ID = running_node_ID;
                cout<<"creating node "<<running_node_ID<<endl;
                new_node->xx = curPtr->xx;
                new_node->yy = curPtr->yy;
                if(curPtr->curEdgePtr == NULL)printf("curEdgePtr is NULL!!!\n");
                else{ new_node->parent_edge = curPtr->curEdgePtr; }
                new_node->child_edge1 = NULL;
                new_node->child_edge2 = NULL;
                new_node->next = NULL;
                add_to_nodeList(nodeList,new_node,running_node_ID);
                running_node_ID++;
                
                EdgePtr new_edge1 = (EdgePtr)malloc(sizeof(Edge));
                new_edge1->ID = running_edge_ID;
                cout<<"creating edge "<<running_edge_ID<<endl;
                new_edge1->left_node = new_node;
                new_edge1->right_node = NULL;
                new_edge1->xn = 0;
                new_edge1->yn = 0;
                add_to_edgeList(edgeList,new_edge1,running_edge_ID);
                running_edge_ID++;
                
                EdgePtr new_edge2 = (EdgePtr)malloc(sizeof(Edge));
                new_edge2->ID = running_edge_ID;
                cout<<"creating edge "<<running_edge_ID<<endl;
                new_edge2->left_node = new_node;
                new_edge2->right_node = NULL;
                new_edge2->xn = 0;
                new_edge2->yn = 0;
                add_to_edgeList(edgeList,new_edge2,running_edge_ID);
                running_edge_ID++;
                
                
                //***************associate two new edges with the new node**********************
                new_node->child_edge1 = new_edge1;
                new_node->child_edge2 = new_edge2;
                //***************associate the new node with the parent edge********************
                if(curPtr->curEdgePtr == NULL)printf("curEdgePtr is NULL!!!\n");
                else (curPtr->curEdgePtr)->right_node = new_node;
                //***************associate the two new edges with two tip epithelial cells******
                curPtr->curEdgePtr = new_edge1;
                newPtr->curEdgePtr = new_edge2;
            */}  
            prePtr = curPtr;
            curPtr = curPtr->next;    
            prePtr->next = newPtr;
            if(curPtr != NULL)newPtr->next = curPtr;
        }
                
        else{
            if(curPtr->branching > 0)curPtr->branching--;
            prePtr = curPtr;  
            curPtr = curPtr->next;
        }              
    }      
    return;
}  

//################################################################################
//## Ecell_looping:                                                              #
//##                                                                             #
//##         This function implements the anastamosis of two branches.           #
//##         Merging of two tips and merging of tip with sprout are treated      #
//##         seperately.                                                         #
//################################################################################
void Elist::Elist_looping(EdgePtr* edgeList, NodePtr* nodeList)
{ 
    Mat_DP& cellmatrix = *cellmatrix_ptr;
    EcellPtr tempPtr, prePtr, curPtr,tPtr;
//merging of two tips   
    curPtr = *Eroot; 
    prePtr = curPtr;
    curPtr = curPtr->next;
    while(curPtr != NULL){
        tPtr = *Eroot;
        while((int)curPtr->xx!=(int)tPtr->xx || (int)curPtr->yy!=(int)tPtr->yy)
            tPtr = tPtr->next;
        if(tPtr != curPtr && curPtr->id > 12){
            cout<<"tip-tip delete branch"<<curPtr->id<<endl;
            tempPtr = curPtr;
            prePtr->next = curPtr->next;
            free(tempPtr);
            curPtr = prePtr->next;
        } 
        else{
            prePtr = curPtr;
            curPtr = curPtr->next; 
        }    
    } 
//merging of tip with sprout
    curPtr = *Eroot;
    prePtr = curPtr;
    curPtr = curPtr->next;
    
    while(curPtr != NULL){                
        if((int)(cellmatrix[(int)curPtr->yy][(int)curPtr->xx])%10 == 1 && curPtr->id>12){
            cout<<"tip-sprout delete branch"<<curPtr->id<<endl;

            
            {/*
            //***************create one new node********************************************
                NodePtr new_node = (NodePtr)malloc(sizeof(Node));
                new_node->ID = running_node_ID;
                new_node->xx = curPtr->xx;
                new_node->yy = curPtr->yy;
                new_node->parent_edge = curPtr->curEdgePtr;
                new_node->child_edge1 = NULL;
                new_node->child_edge2 = NULL;
                add_to_nodeList(nodeList,new_node,running_node_ID);
                running_node_ID++;
            //***************identify the original vessel***********************************
                EdgePtr old_edge;
                int search_ID = edgeID_matrix[(int)curPtr->yy][(int)curPtr->xx];
                if(search_ID == -1)printf("Does not return a valid edge.\n");
                else search_edgeList(edgeList, old_edge, search_ID);
            //***************break the original vessel into two vessels****************                
                EdgePtr new_edge = (EdgePtr)malloc(sizeof(Edge));
                new_edge->ID = running_edge_ID;
                new_edge->left_node = new_node;
                new_edge->right_node = old_edge->right_node;
                old_edge->right_node = new_node;
                new_edge->xn = 0;
                new_edge->yn = 0;
                add_to_edgeList(edgeList,new_edge,running_edge_ID);
                running_edge_ID++;
            //***************associate the new node with the parent edge and the two new edges*************
                (curPtr->curEdgePtr)->right_node = new_node;
            */}
            tempPtr = curPtr;
            prePtr->next = curPtr->next;
            free(tempPtr);
            curPtr = prePtr->next;
            cout<<"deletion done."<<endl;

        }                
        else{
            prePtr = curPtr;
            curPtr = curPtr->next; 
        }  
    }  
    
    return;
}  

//################################################################################
//## update_cellmatrix:                                                          #
//##                                                                             #
//##         Update the matrix recording vasculature                             #
//##         Vasculature is recorded as the trace of migrating tip cells.        #
//##         Therefore updates only happen to the newly occupied grid by tip     #
//##         cells.                                                              #
//##                                                                             #
//################################################################################
void Elist::update_cellmatrix()
{ 
    Mat_DP &cellmatrix = *cellmatrix_ptr; 
    int i_range,j_range;
    
    for(int i=0;i<mYstep;i++) {
        for(int j=0;j<mXstep;j++) {
            if((int)(cellmatrix[i][j])%10 == 3)cellmatrix[i][j] -= 2;//11,21,31,41,... represents old vasculature
            if((int)(cellmatrix[i][j])%10 == 4)cellmatrix[i][j] -= 2;//12,22,32,42,... represents grids surrounding old vasculature
            if(cellmatrix[i][j] == 1){
                for(i_range = i-2;i_range <= i+2;i_range++){
                    for(j_range = j-2;j_range <= j+2;j_range++){
                        if(cellmatrix[i_range][j_range] == 0)cellmatrix[i_range][j_range] = 2;
                    }
                } 
            } 
            
            if(cellmatrix[i][j] == 5){//5 represents vessel sections with new tip formation potential
                for(i_range = i-2;i_range <= i+2;i_range++){
                    for(j_range = j-2;j_range <= j+2;j_range++){
                        if(cellmatrix[i_range][j_range] == 0)cellmatrix[i_range][j_range] = 6;
                        //6 represents grids surrounding the vessel section with new tip formation potential
                    }
                } 
            } 
            if(cellmatrix[i][j] == 7){//7 represents position of new tip cell
                for(i_range = i-2;i_range <= i+2;i_range++){
                    for(j_range = j-2;j_range <= j+2;j_range++){
                        if(cellmatrix[i_range][j_range] == 0)cellmatrix[i_range][j_range] = 8;
                        //8 represents grids surrounding new tip cell
                    }
                } 
            }  
                          
        }
    }
        
    EcellPtr ecurPtr = *Eroot; 
    
    while(ecurPtr != nullptr){
        int xx = (int)ecurPtr->xx;
        int yy = (int)ecurPtr->yy;
        cellmatrix[yy][xx] = ecurPtr->id*10 + 3;//3 represents newly formed vasculature
        for(int i=max(0, yy-2); i<=min(mYstep-1, yy+2); i++){
            for(int j=max(0, xx-2); j<=min(mXstep-1, xx+2); j++){
                if(cellmatrix[i][j] == 0)cellmatrix[i][j] = ecurPtr->id*10 + 4;//4 represents grids surrounding new vasculature
            }
        } 
        ecurPtr = ecurPtr->next;
    }
}


void Elist::output_cellmatrix(){ 
    Mat_DP &cellmatrix = *cellmatrix_ptr;
    unsigned char b_0[3];
    unsigned char b_16[3];
    unsigned char b_green[3];
    unsigned char b_red[3];
    for(int i=0; i<3; i++){
        b_0[i] = 15;
        b_16[i] = 255;
        b_green[i] = 0;
        b_red[i] = 0;
    }
    b_green[1] = 255;
    b_red[2] = 255;

    static int out_i = 0;
    char filename[21];
    sprintf(filename, "output/vasc%05d.BMP", out_i++);
    FILE* fp_BMP = fopen(filename,"w");
    BMP::Prepare_BMP_Format(fp_BMP, mXstep,mYstep);
    for(int i=0;i<mYstep;i++) {
        for(int j=0;j<mXstep;j++) {
            if(cellmatrix[i][j] == 0)fwrite(&b_0, 1,3,fp_BMP);
            if((int)(cellmatrix[i][j])%10 >= 1 && (int)(cellmatrix[i][j])%10 <= 4)fwrite(&b_16,1,3,fp_BMP); 
            if(cellmatrix[i][j] == 5 || cellmatrix[i][j] == 6)fwrite(&b_green,1,3,fp_BMP);
            if(cellmatrix[i][j] == 7 || cellmatrix[i][j] == 8)fwrite(&b_red,1,3,fp_BMP);
        }
    }

    fclose(fp_BMP); 
}



void Elist::calculate_for_oxygen(Mat_DP& vasculature_density)
{
     Mat_DP& cellmatrix = *cellmatrix_ptr;
     int i,j,i_range,j_range;
     DP omega;
     for(i=0;i<mYstep; i++){
         for(j=0; j<mXstep; j++){
             vasculature_density[i][j]=0;
             for(i_range=max(0,i-20);i_range<min(mYstep-1,i+20);i_range++){
                 for(j_range=max(0,j-int( sqrt( 20*20-(i-i_range)*(i-i_range) ) ) );j_range<min(mXstep-1,j+int( sqrt( 20*20-(i-i_range)*(i-i_range) ) ) );j_range++){
                     if(cellmatrix[i_range][j_range]==1 || cellmatrix[i_range][j_range]==2){
                         omega = (1-abs(i_range-i)/20) * (1-abs(j_range-j)/20);
                         vasculature_density[i][j] += omega;
                     }
                 }
             }
             vasculature_density[i][j] = vasculature_density[i][j]/20;             
        }
    }
    return;
}
  

void Elist::output_vasculature_density()
{
    static int out_i = 0;
    char  file_name[14];
    sprintf(file_name, "vden%05d.BMP", out_i++);
    BMP::output_BMP(file_name, 14, *vas_density_ptr, mYstep, mXstep);
    return;
}
     


void feedback_from_oxygen_to_VEGF(Mat_DP& vasculature_density,Mat_DP& VEGF)
{
     for(int j=0;j<1000;j++){
         for(int l=0;l<1000;l++){
             VEGF[j][l] = VEGF[j][l] - vasculature_density[j][l];
         }
     } 
     return;
}     
                 
//################################################################################
//## Ecell_newtip:                                                               #
//##                                                                             #
//##         This function implements the formation of new vascular tip.         #
//##                                                                             #
//##                                                                             #
//################################################################################
void Elist::Ecell_newtip(Mat_DP &activator, Mat_DP &inhibitor)
{
    DP Athreshold = 0.5;
    DP Ithreshold = 0.3;
    Mat_DP &cellmatrix = *cellmatrix_ptr;
    //find the section of vessel with a potential of new tips formation====================
    int j_range,l_range; 
    Mat_DP vesseltheta(mYstep,mXstep);
    for(int j=0;j<mYstep;j++){
        for(int l=0;l<mXstep;l++){
            vesseltheta[j][l] = -1;
            if(cellmatrix[j][l] == 1){                    
                //find the adjacent cellmatrix[j][l]==1
                //if only one line
                //use regression to find the direction of vessel
                Vec_DP xlist(800),ylist(800);
                int count=0;
                for(j_range=max(0,j-20);j_range<min(mYstep,j+20);j_range++){
                    for(l_range=max(0,l-20);l_range<min(mXstep,l+20);l_range++){
                        if((int)(cellmatrix[j_range][l_range])%10 ==1 || cellmatrix[j_range][l_range]==5){
                            ylist[count]=(double)j_range;
                            xlist[count]=(double)l_range;
                            count++;
                        }
                    }
                }
                cout<< "j=" << j << "l=" << l << "count = " << count << endl;
                //system("PAUSE");
                DP b,xycor=0;
                if(count>10)regression(xlist,count,ylist,count,&b,&xycor);
                if(abs(xycor)>0.8 and activator[j][l] > Athreshold and inhibitor[j][l] < Ithreshold){
                    cellmatrix[j][l]=5;
                    if(b>=0)vesseltheta[j][l]=atan(b);
                    if(b<0)vesseltheta[j][l]=M_PI+atan(b);    
                }
            }
            
        }
    }
    
    
    
    //==============find the exact position for new tip formation======================    
    for(int j=0;j<mYstep;j++){
        for(int l=0;l<mXstep;l++){
            if(cellmatrix[j][l] == 5){                    
            //find the adjacent cellmatrix[j][l]==5
            //for a series of 5, find the middle one=======================
            //use as the location of new tip
                int count_5=0, count_7=0;
                for(j_range=max(0,j-20);j_range<min(mYstep, j+20);j_range++){
                    for(l_range=max(0,l-20);l_range<min(mYstep, l+20);l_range++){
                        if(cellmatrix[j_range][l_range]==5)count_5++;
                    }
                }
                
                for(j_range=max(0,j-30);j_range<min(mYstep, j+30);j_range++){
                    for(l_range=max(0,l-30);l_range<min(mXstep, l+30);l_range++){
                        if(cellmatrix[j_range][l_range]==7)count_7++;
                    }
                }
                
                //cout<< "j=" << j << "l=" << l << "count = " << count << endl;
                //system("PAUSE");
                if(count_5>35 && count_7==0)cellmatrix[j][l]=7;   
            } 
        }
    }
    

    
    //========build Elist===========================================================
    EcellPtr newPtr,prePtr,curPtr;
       
    for(int j=0;j<mYstep;j++){
        for(int l=0;l<mXstep;l++){
            if(cellmatrix[j][l] == 7){                                                   
                double theta = vesseltheta[j][l] + M_PI/2;
                if (theta > 2 * M_PI) theta -= 2 * M_PI;
                if (theta < 0) theta += 2 * M_PI;
                newPtr = new Endocyte(l, j, theta, branch_id++, 24);       
                newPtr->next = NULL;

                prePtr = NULL;
                curPtr = *Eroot;
                while(curPtr != NULL){
                    prePtr = curPtr;
                    curPtr = curPtr->next;
                }
                if (prePtr == NULL){
                    newPtr->next = *Eroot;
                    *Eroot = newPtr;
                }else{
                    prePtr->next = newPtr;
                    newPtr->next = curPtr;
                }            
            }
        }
    }
    return; 
} 

DP speed_function(DP grad){
    DP Th1 = 1.5;
    DP Th2 = 6;
    if(grad < Th1) {
        return 1;
    } else if (grad < Th2) {
        return 1 + 3 * (grad - Th1) / (Th2 - Th1);
    } else {
        return 4;
    }
}
 
DP angle_substract(DP theta1,DP theta2){
    if(theta1<M_PI/2 && theta1>=0 && theta2>3/2*M_PI && theta2<2*M_PI)return theta1 + 2*M_PI -theta2;
    else if(theta2<M_PI/2 && theta2>=0 && theta1>3/2*M_PI && theta1<2*M_PI)return -(theta2 + 2*M_PI -theta1);
    else return theta1-theta2;
}

void regression(Vec_DP xlist,int xcount,Vec_DP ylist,int ycount,double* b,double* xycor){
     int i,j;
     double xmean=0,ymean=0;
     double xvar=0,yvar=0,xdev=0,ydev=0,xycov=0;
     for(i=0;i<xcount;i++){
         xmean+=xlist[i];
         ymean+=ylist[i];
     }
     xmean=xmean/xcount;
     ymean=ymean/ycount;
     
     for(i=0;i<xcount;i++){
        xvar+=pow(xlist[i]-xmean,2);
        yvar+=pow(ylist[i]-ymean,2);
        xycov+=(xlist[i]-xmean)*(ylist[i]-ymean);
     }
     xdev=sqrt(xvar);
     ydev=sqrt(yvar);
     printf("xvar=%f,yvar=%f,xdev=%f,ydev=%f\n",xvar,yvar,xdev,ydev);
     *xycor=xycov/(xdev*ydev);
     *b=xycov/xvar;
     printf("xycor=%f,b=%f\n",*xycor,*b);
     return;
}
