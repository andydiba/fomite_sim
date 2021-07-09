#ifndef PROBABILITY_H
#define PROBABILITY_H

#define PI 3.141592653589793


#define uniform() ( (double)rand()/RAND_MAX )


//Truncated Gaussian Parameters
typedef struct tgParam{

double md; //mode
double sd; //standard deviation
double mn; //min
double mx; //max

}tgParam;

double normTrunc(double ,tgParam*);

double sampleNormT(tgParam*);

double erf(double );

bool accept(double );

bool hasHappenedByNow(double  , double );

unsigned int numEvents(double );


#endif