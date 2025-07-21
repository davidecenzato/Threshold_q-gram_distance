#include "Profile.h"
#include <iostream>
#include <stdlib.h>
using namespace std;

//--------------------------------------------------------------------------
//creator
Profile::Profile(int thr){
	threshold = thr;
	root = NULL;
	update = false;
	total_status = (int *)malloc(sizeof(int) * (threshold+2) * (threshold+2));
	int i=0;
	for(i=0;i<(threshold+2) * (threshold+2);i++){
		total_status[i]=0;
	}
	number_of_nodes = 0;
	
	nil = (node *)malloc(sizeof(node));
	nil->color = 0;
}

//--------------------------------------------------------------------------
//destructor
Profile::~Profile(){
	
}

//--------------------------------------------------------------------------
//insert
void Profile::insert(unsigned long long int qgram_ins, int n_str){
	update = false;//when I insert a node, the total_status is not update
	
	//if I have to create the root
	if(root == NULL){
		root = (node *)malloc(sizeof(node));
		root->qgram = qgram_ins;
		if (n_str == 1){
        	root->size1 = 1;
        	root->size2 = 0;
        }else{
        	root->size1 = 0;
        	root->size2 = 1;
        }
        root->sx = nil;
        root->dx = nil;
        root->father = nil;
        root->color = 0; // the root is black
        number_of_nodes = 1;
        return;
	}
	
	
	//if the root exists
	node *y=NULL;
	node *x=root;
	while (x != nil){
		if (x->qgram == qgram_ins){ //if we have find the q-gram and we have to change only size1 or size2
			if (n_str == 1){
        		if (x->size1 < threshold + 1){
        		 x->size1 = x->size1 + 1;
        		 }
			}else{
        		if (x->size2 < threshold + 1){
        		 x->size2 = x->size2 + 1;
        		 }
        	}
        	return;
		}
		
		y=x;
		if (qgram_ins < x->qgram){
			x = x->sx;
		}
		else{
			x = x->dx;
		}
	}
	//if we are here, we are in this situation:
	//x is nil and y is a q_gram and we have to create a new node and insert in x
	
	//Create and insert the new node z
	number_of_nodes = number_of_nodes+1;
	node *z = (node *)malloc(sizeof(node));
	z->qgram = qgram_ins;
	if (n_str == 1){
        z->size1 = 1;
        z->size2 = 0;
	}else{
        z->size1 = 0;
        z->size2 = 1;
	}
	z->sx = nil;
    z->dx = nil;
    
    z->father = y;
	z->color = 1;
	if (qgram_ins < y->qgram)
		y->sx=z;
	else
		y->dx=z;
		
	//fixing the rules of RB trees
	RB_insert_fixup(z);
	return;
}


//--------------------------------------------------------------------------
//left-rotate
void Profile::left_rotate(node *x){
	node *y = x->dx;
	x->dx=y->sx;
	
	if(y->sx != nil)
		y->sx->father = x;
	y->father = x->father;
	if (x->father == nil)
		root = y;
	else{
		if (x==x->father->sx)
			x->father->sx=y;
		else
			x->father->dx=y;
	}
	y->sx = x;
	x->father = y;
}
//--------------------------------------------------------------------------
//right-rotate
void Profile::right_rotate(node *x){
	node *y = x->sx;
	x->sx = y->dx;
	
	if(y->dx != nil)
		y->dx->father = x;
	y->father = x->father;
	if (x->father == nil)
		root = y;
	else{
		if (x==x->father->dx)
			x->father->dx=y;
		else
			x->father->sx=y;
	}
	y->dx = x;
	x->father = y;
}
//--------------------------------------------------------------------------
//RB-insert fixup
void Profile::RB_insert_fixup(node *z){
	node *y=NULL;
	while(z->father->color == 1){
		if (z->father == z->father->father->sx){
			y = z->father->father->dx;
			if (y->color == 1){
				z->father->color=0;
				y->color=0;
				z->father->father->color=1;
				z = z->father->father;
			}else{
				if(z==z->father->dx){
					z=z->father;
					left_rotate(z);
				}
				z->father->color=0;
				z->father->father->color=1;
				right_rotate(z->father->father);
			}
			
		}else{
			y = z->father->father->sx;
			if (y->color == 1){
				z->father->color=0;
				y->color=0;
				z->father->father->color=1;
				z = z->father->father;
			}else{
				if(z==z->father->sx){
					z=z->father;
					right_rotate(z);
				}
				z->father->color=0;
				z->father->father->color=1;
				left_rotate(z->father->father);
			}
			
			
		}
		
		
	}
	root->color=0;
}






//--------------------------------------------------------------------------
//get the total status
int* Profile::get_total_status(){
	if(!update){ //if you want it immediatly after calculateDistance, it requires O(1)
		calculateDistance();
	}
	return total_status;
}

//--------------------------------------------------------------------------
//calculate the distance
unsigned long long int Profile::calculateDistance(){
	
	//put zeros in total_status
	int i=0;
	for(i=0;i<(threshold+2) * (threshold+2);i++){
		total_status[i]=0;
	}
	
	//calculate the distance recursively, starting from the root
	long long int dist;
	dist = calculateDistanceRecursive(root);
	update = true;//when I calculate the distance, the total_status is update
	return dist;
}

//--------------------------------------------------------------------------
//calculate the distance recursively
unsigned long long int Profile::calculateDistanceRecursive(node *n){
	long long int res = 0;
	if (n == nil) //we arrive to the leaf
		return 0;
	
	if(n->size1 != n->size2) res = 1;
	res = res + calculateDistanceRecursive (n->sx) + calculateDistanceRecursive (n->dx);
	total_status[(n->size1*(threshold+2) )+n->size2] = total_status[(n->size1*(threshold+2) )+n->size2] + 1;
	return res;
}

//--------------------------------------------------------------------------
//get the number of nodes
unsigned long long int Profile::getNumber_of_nodes(){
	return number_of_nodes;
}
