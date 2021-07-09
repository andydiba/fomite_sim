#ifndef QUEUE_H
#define QUEUE_H


#include <stdbool.h>


#define INIT_QUEUE_LENGTH 10

typedef struct Queue Queue; 

Queue* new_queue(unsigned int );

void* get_queue_data(Queue*);

unsigned int sizeOfQueueStruct(void);

unsigned long queue_count(Queue*);

void stats_queue(Queue* );

bool isEmpty_queue(Queue*);

void delete_queue(Queue* );

void* getQueue_data(void);

void enqueue(Queue* , void*);

void dequeue(Queue*, void*);

void peek_queue(Queue*, void*,unsigned long);


#endif