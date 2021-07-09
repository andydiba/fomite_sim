/*
Copyright 2020 Andrew Di Battista andrew.di.battista@ultraleap.com

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/



#include <stdio.h>
#include <stdlib.h>
#include "queue.h"


struct Queue{

	void* data;
	unsigned long front;
	unsigned long CAP;
	unsigned long count;
	unsigned int size_b;

};

void* get_queue_data(Queue* q){
	return q->data;
}

unsigned int sizeOfQueueStruct(void){return sizeof(Queue);}

unsigned long queue_count(Queue* q){return q->count;} 

Queue* new_queue(unsigned int sizeOfItem){

	Queue* q;

	q = malloc(sizeof(Queue));

	q->size_b = sizeOfItem; 
	q->data = malloc(INIT_QUEUE_LENGTH*sizeOfItem); 
	q->front = 0;
	q->CAP = INIT_QUEUE_LENGTH;
	q->count = 0;

	return q;
}


void delete_queue(Queue* q){  
	//queue Object can still be re-initialised with new_queue
	//free data arrays....
	if(q->data != NULL)
	free(q->data);

	q->data = NULL;
	q->front = 0;
	q->count = 0;
	q->CAP = 0;
	q->size_b = -1;

}

bool isEmpty_queue(Queue* q){return (q->count==0);}

void stats_queue(Queue* q){

	printf("\n queue stats: count: %lu, front: %lu, capacity: %lu \n",
	q->count,q->front,q->CAP);	
}




void enqueue(Queue* q, void* item){

	//size of item already intialised size_b = number of bytes
	//ASSSUMES new_queue

	unsigned long pos,i,j;
	unsigned int k,nBytes;
	char* qptr;
	char* iptr;


	qptr = (char*)(q->data); //a byte pointer
	iptr = (char*)(item);

	nBytes = q->size_b;
	
	if ( (q->count) <  (q->CAP) ){
		//i.e. there's room for more, go to 'end'
		pos = (q->front + q->count)%(q->CAP);		//logical index
		pos = pos*(nBytes);				//index in bytes

		for(k=0;k<nBytes;k++)
		{
			qptr[pos+k] = iptr[k];		
		}
    	
    		q->count = q->count+1;

	}else{
	//%count == CAP  i.e. already full so double-up, copy and swap
   	    
	    unsigned long newMax = 2*(q->CAP);
	    char* temp = malloc(newMax*nBytes);

	    //copy existing
	    for(i=0;i<q->CAP;i++)
	    {
		pos = (q->front+i)%(q->CAP); //logical index
		pos = pos*nBytes;	//index in bytes

		j = i*nBytes; //byte index 

		for(k=0;k<nBytes;k++)
		{
		   temp[j+k] = qptr[pos+k];		
		}

	    }
		
	    //append new item
	    j = (q->CAP)*nBytes; 		   
	    for(k=0;k<nBytes;k++)
	    {
		temp[j+k] = iptr[k];		
	    }	  

	    
	    q->count = q->count+1;
	    q->CAP = newMax;
	
	    free(q->data);
	
	    q->data = temp;
	    q->front = 0; 		    //new queue starts at index 0 again
	}

}


void dequeue(Queue* q, void* dest){

	unsigned long front,pos;
	unsigned int k, nBytes;
	
	char* qptr;
	char* dptr;

	
	
	 if(q->count > 0)
    	 {	

		qptr = (char*)(q->data); //a byte pointer
		dptr = (char*)(dest);

		nBytes = q->size_b;

		front = q->front; //logical index

		pos = front*nBytes; //byte index
  		
		
		for(k=0;k<nBytes;k++)
		{
		    dptr[k] = qptr[pos+k];		
		}	  
  	
		//new front
    		q->front = (front + 1)%(q->CAP); 
   		q->count = q->count -1;
 	 }else{
		
            printf("\nWarning: attempted read from empty queue. Check first!");	
    	 }
  	
}  

void peek_queue(Queue* q, void* dest, unsigned long index){
	//similar to dequeue but does not affect front or count
	unsigned long offset,pos;
	unsigned int k, nBytes;
	
	char* qptr;
	char* dptr;

	 if(index < q->count)
    	 {	
		qptr = (char*)(q->data); 			//a byte pointer
		dptr = (char*)(dest);

		nBytes = q->size_b;

		offset = ( (q->front)+index )%(q->CAP); 	//logical index

		pos = offset*nBytes; 				//byte index
  			
		for(k=0;k<nBytes;k++)
		{
		    dptr[k] = qptr[pos+k];		
		}	  
  	
 	 }else{
		
            printf("\nWarning: Index out of bounds or empty queue!");	
    	 }
}  
    
 





	

