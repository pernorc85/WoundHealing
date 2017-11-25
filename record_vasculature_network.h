#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>//pow
#include "nr.h"
#include "truecolor.h"
#include "ECM.h"


extern int xstep, ystep;
static int running_node_ID=0;
static int running_edge_ID=0;


typedef struct Edge* EdgePtr;
typedef struct Node* NodePtr;
struct Node
{
        int ID;
        DP xx;
        DP yy;
        EdgePtr parent_edge;
        EdgePtr child_edge1;
        EdgePtr child_edge2;
        NodePtr next;       
};

struct Edge
{
        int ID;
        NodePtr left_node;
        NodePtr right_node; 
        int trace_x[200];
        int xn;
        int trace_y[200];       
        int yn;
        EdgePtr next;
};
        
NodePtr* NodeList;
EdgePtr* EdgeList;        
Mat_DP edgeID_matrix(1000,1000);
        

//in the main program
//
//if(branch)
//{
////***************create one new node and two new edges**************************
//    create_node(new_node, running_node_ID, position_x, position_y);
//    child_edge1 = new edge;
//    create_edge(child_edge1, running_edge_ID, new_node, NULL);
//    child_edge2 = new edge;
//    create_edge(child_edge2, running_edge_ID, new_node, NULL);
////***************associate two new edges with the new node**********************
//    associate_edge(new_node, parent_edge, child_edge1, child_edge_2);
////***************associate the new node with the parent edge********************
//    record_edge(parent_edge, new_node,2);
//}

//if(looping)
//{
////***************create one new node********************************************
//    create_node(new_node, running_node_ID, position_x, position_y);
////***************break the original vessel into two vessels****************
//    identify_edge(edgeID,edgePtr);
//    record_edge(edgePtr, newnode);
//    edge1 = new edge;
//    node2 = edgePtr->right_node;
//    create_edge(edge, newnode, node2);
////***************associate the new node with the parent edge and the two new edges*************
//    associate_edge(new_node, parent_edge, edgePtr, edge1);
//}

//void organize_into_matrix(NodePtr* nodelist, int node_count, EdgePtr* edgelist, int edge_count)
//{
//    matrix(node_count+edge_count, node_count+edge+count);
//    //from i=0 to i=edge_count, row i is associated with the ith vessel......
//    matrix(pressure of left node-presssure of right node = flow*resistance){
//        //we don't know the direction of flow...
//    //from i=edge_count+1 to i=edge_count+node_count
//    //row i is associated with the ith node......
//    matrix(summation of flows at this node is zero);
//    return;
//}
    
void add_to_edgeList(EdgePtr* edgeList,EdgePtr new_edge,int running_edge_ID);
void search_edgeList(EdgePtr* edgeList,EdgePtr old_edge,int edge_ID); 
void add_to_nodeList(NodePtr* nodeList,NodePtr new_node,int running_node_ID);  

void add_to_edgeList(EdgePtr* edgeList,EdgePtr new_edge,int running_edge_ID){
    cout<<"\tAdding edge to the edgeList...";
    EdgePtr prePtr, curPtr;
    int count=0;
    curPtr = *edgeList;
    while(count < running_edge_ID && curPtr != NULL){
        prePtr = curPtr;
        curPtr = curPtr->next;
        count++;
    }
    if(count == running_edge_ID && curPtr == NULL)cout<<"Correct order."<<endl;
    else cout<<"Incorrect order!!!"<<endl;
    if(prePtr == NULL){
        new_edge->next =*edgeList;
        *edgeList = new_edge;
        cout<<"\t...edge added to the start of the edgeList."<<endl;
    }
    else{
        prePtr->next = new_edge;
        new_edge->next = curPtr;
        cout<<"\t...edge added to the end of the edgeList."<<endl;
    } 
    return;
}

void search_edgeList(EdgePtr* edgeList,EdgePtr old_edge,int edge_ID){
     cout<<"\tSearching edge with ID "<<edge_ID<<endl;
     EdgePtr curPtr;
     curPtr = *edgeList;
     while(edge_ID != curPtr->ID){
         curPtr = curPtr->next;
     }
     old_edge = curPtr;
     return;
}

void add_to_nodeList(NodePtr* nodeList,NodePtr new_node,int running_node_ID){
    cout<<"\tAdding node to the nodeList...";
    NodePtr prePtr, curPtr;
    int count=0;
    curPtr = *nodeList;
    cout<<"ok...";
    while(count < running_node_ID && curPtr != NULL){
        prePtr = curPtr;
        curPtr = curPtr->next;
        count++;
    }
    if(count == running_node_ID && curPtr == NULL)cout<<"Correct order."<<endl;
    else cout<<"Incorrect order!!!"<<endl;
    if(prePtr == NULL){
        new_node->next =*nodeList;
        *nodeList = new_node;
        cout<<"\t...node added to the start of the nodeList."<<endl;
    }
    else{
        prePtr->next = new_node;
        new_node->next = curPtr;
        cout<<"\t...node added to the end of the nodeList."<<endl;
    } 
    return;   
}

void output_edgeID_matrix(EdgePtr* edgeList, int output){
     cout<<"output edgeID_matrix..."<<endl;
     int i,j,l;
     for(j=0;j<ystep;j++) {
	    for(l=0;l<xstep;l++) {
            edgeID_matrix[j][l] = -1;
        }
     }
     
     static int out_i=0;
     EdgePtr curPtr;
     curPtr = *edgeList;
     while(curPtr != NULL){
         for(i=0;i<curPtr->xn;i++){
             cout<<curPtr->trace_y[i]<<" "<<curPtr->trace_x[i]<<endl; 
             edgeID_matrix[(curPtr->trace_y)[i]][(curPtr->trace_x)[i]] = curPtr->ID;
         }
         curPtr = curPtr->next;
     }

    if(output == 1){
        char  file_name[15]="edge00000.BMP";
        int jj;
       	unsigned char true_color[3];
        file_name[8] = 48+out_i%10;        file_name[7] = 48+out_i/10%10;
        file_name[6] = 48+out_i/10/10%10;  file_name[5] = 48+out_i/10/10/10%10;
        fp_BMP = fopen(file_name,"w") ;
        BMP::Prepare_BMP_Format(1000,1000) ;
        for(j=0;j<ystep;j++) {
    	    for(l=0;l<xstep;l++) {
                jj = edgeID_matrix[j][l];
                if(jj>=0){
                    Get_True_Color(jj, 200, 0, true_color[0], true_color[1], true_color[2]);
                    cout<<"truecolor"<<(int)true_color[0]<<" "<<(int)true_color[1]<<" "<<(int)true_color[2]<<endl;
    	            fwrite(&true_color, 1,3,fp_BMP);
                }
                else{
                    fwrite(&b0, 1,3,fp_BMP);
                }
            }
        }
        out_i++;
        fclose(fp_BMP); 
    }

    return;                           
}    
     
