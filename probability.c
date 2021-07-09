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
#include <math.h>
#include <stdbool.h>

#include "probability.h"


//returns TRUE with probability 'p'
bool accept(double p){

	double pb = uniform(); //[0,1]
  
  	if(pb<p)
  	 return true;//accept

	return false;//reject
}


//use an Exponential cummulative-distribution function (rate L) to decide if an event has happened by time 't'
bool hasHappenedByNow(double L , double t){
  
      double p = 1-exp(-L*t);
      return accept(p);
  
}

//Error function : needed for truncated gaussian computation
double erf(double x){

	double a1 =  0.254829592;
   	double a2 = -0.284496736;
   	double a3 =  1.421413741;
   	double a4 = -1.453152027;
   	double a5 =  1.061405429;
   	double p  =  0.3275911;
  
  	//Save the sign of x
    	int sgn = 1;
    	if (x < 0)
        	sgn = -1;
   
        x = fabs(x);
  
  
  	//A&S formula 7.1.26
   	double t = 1.0/(1.0 + p*x);
   	double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    	return sgn*y;

}

//Truncated (Gaussian) Normal Distribution
double normTrunc(double x,tgParam* P){


	   double md = P->md;
	   double sd = P->sd;
	   double mn = P->mn;
	   double mx = P->mx;


    
            if(x<mn || x>mx)
               return 0;
            else if (sd == 0){
             	
		return (x==md);
		
            }else{  
              
              double X1 = (x-md)/sd;
              double Xb = (mx-md)/sd;
              double Xa = (mn-md)/sd;
              
              double PHI1 = (1/sqrt(2*PI))*exp(-0.5*pow(X1,2));
              double PHIb = 0.5*(1+ erf(Xb/sqrt(2)) );
              double PHIa = 0.5*(1+ erf(Xa/sqrt(2)) );
              
              return (1/sd)*PHI1/(PHIb - PHIa);
            }  
            
}

//Sample from trunc norm distribution using 'rejection' sampling method
double sampleNormT(tgParam* P){
  
        bool done = false;
        unsigned int max_count = 1000;
        unsigned int count=0;


	double md = P->md;
	double sd = P->sd;
	double mn = P->mn;
	double mx = P->mx;
        
	double maxPDF,x,y,gx,r1,r2;


	if (sd == 0){
              if (md<=mx && md>=mn)
		{	return md;}

		return 0;	
        }  


        if (md>mx)
             maxPDF = normTrunc(mx,P);
        else if (md<mn)
             maxPDF = normTrunc(mn,P);
        else
             maxPDF = normTrunc(md,P);
        
        
        double xdom = mx-mn;
        
        while(!done && count<max_count)
        {
	
	  r1 = uniform();
	  r2 = uniform();

          x = mn+xdom*r1;
          y = maxPDF*r2;

          gx = normTrunc(x,P); 
          
          if(y <= gx)
              done = true;
                         
           count++;
        }
  
        if(count==max_count)
           printf("\nOops! max count reached!\n");

	return x;
        
}


//generates samples from Poisson distribution, in 1 unit of time
unsigned int numEvents(double L){
	
	unsigned int k;
	double t,r;	

	 if(L==0)
         { 
		k = 0;
    	 }else{
  
              t=0; //time of event
              k=0; //number of events (event index)

              //S=[]; //arrival times
	
 	         while(t<1)
		  {
 	            r = uniform();
 	            t = t - log(r)/L;
 	            k++;
 	           //S=[S; t];
 	         }  
                 
 	         //cut off last iteration
 	         //S=S(1:end-1);
 	         k=k-1;
	 }
	return k;

}


