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
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>



//import basic queue structure and special probability functions
#include "queue.h"
#include "probability.h"



//person status macros
#define INFECTIOUS 0
#define SUSCEPTIBLE 1
#define INFECTED 2



//define 'pathogen' structure
typedef struct Pathogen{

	double shed[2]; //pathgen shedding sampled from uniform distribution: U(shed[0],shed[1])
  	double ID50; //could store exponent i.e. 10^ID50?
  	//double muh;  //pathogen half-life (on skin)
  	double muf;  //pathogen half-life (on fomite)
}Pathogen;

//define 'people'
typedef struct Person{
	int status;	//Infectious?, Susceptible ? etc.
	double fdose;   //number of pathogens (dose) on finger
	unsigned int jumps; //controls the number of new location person can visit

	//simulation stats
	int gap; //for infected individuals, how many users after infectious contamination occured?

}Person;

//define 'fomite' or touchscreen interface
typedef struct Fomite{
      tgParam* use_stats;  //number of touches sampled from trunc. norm. dist.
      double L; 	   //rate of use/availability  (per minute)
      double dose;         //pathogen dose level on surface

      double CR;	   //cleaning rate (number of cleans (per fomite) per day

      //simulation stats
      int gap_count; 	   //resets when infectious user, counts up otherwise
}Fomite;

//Define a 'location', also used for initial population 'pool'
typedef struct Location{
	Queue* A;			//'Arrival' queue
	Queue* D;			//'Departure' queue

	Fomite* foms;
	unsigned int numFomites;

	double INTERACTION_RATE;      //default set to 1; has no effect.

	double L; 	//average rate (per minute) of people leaving location

	/*when mutiple locations are simualted, i_next stores the index of the
	other potential locations, p_next is the probability of going there. i.e.
	it's a Markov Chain. numNext is simply set to the total number of possible
	destinations from current location: for the inital pool = NUM_LOCATION else
	 NUM_LOCATIONS-1*/

	unsigned int* i_next;
	double* p_next;
	unsigned int numNext;

	//simulation stats
	unsigned int NI; 		//number of infections at this location


}Location;

//Markov Chains: given a person at Location[current], where to go next?
unsigned int getNewLocation(Location* loc, unsigned int current){

	double r,sum;
	unsigned int i,N;

	N = loc->numNext;

	if(N==0) //current location points to nowhere but maybe nLoc still !=0. Loop back to A queue of self
	return current;

	//else try a new location
	r = uniform();
	sum = 0;
	for(i=0; i<N;i++)
	{
		sum=sum+(loc->p_next[i]);

		if(r<sum)
		return (loc->i_next[i]);

	}

	printf("\n Warning: getNewLocation error");
	return current; //should never reach this

}



//Computes pathogen shedding (from finger) of infectious individual.
double DoseOnFinger(Pathogen* p){

	double dose;
	double r = uniform();

	dose = p->shed[0]+ (p->shed[1]-p->shed[0])*r;
	dose = dose/(p->ID50); //normalised
	dose = 0.014*dose; //finger tip proportion ~1.4or1.5%

	return dose;

}


//some macros used exclusively in selfInoc
#define INOC_PERIOD_MIN 20  //minutes time for innoculation must happen
#define FACE_TOUCH_RATE_HR 15  //face touch rate of typical person

//bool selfInoc(Person* person, Pathogen* pathogen,tgParam* face){
bool selfInoc(double fdose, Pathogen* pathogen,tgParam* face){


	//double fdose = person->fdose;
	//double muh = pathogen->muh;

	double inoc_period_hr = (double)INOC_PERIOD_MIN/60;


       // double avgSurvival_Fraction = (1/( log(muh)*inoc_period_hr) )*( pow(muh,inoc_period_hr)  -1);
        double avgSurvival_Fraction=1; //could use for simualting survival/decay on hands

	//Number of face touch event based on hourly rate normalised over 20 minutes
        unsigned int nt  = numEvents((double)FACE_TOUCH_RATE_HR*inoc_period_hr);

	double  dep  = sampleNormT(face);  //transfer (i,e, deposited) onto face

	double  mdose =  (fdose*avgSurvival_Fraction)*(1 - pow( (1-dep), nt)  ); //dose to mucosal membranes

	//probability of infection mdoses are in ID50 so each contributes a half chance
	double  p = 1 - exp(-0.7*mdose);

	//remove mdose from hands
	//person->fdose -= mdose;

	return accept(p);

}



//main simulation loop

void fomite_sim(Location* pool,
	Location* loc, unsigned int numLoc,
	tgParam* DEPOSIT,
	tgParam* PICKUP,
	tgParam* FACE,
	Pathogen* pathogen,
	double ARRIVAL_RATE_MIN
)
{
	unsigned int i,m,n,p,f, num_fom,nl, totalN;

	double tStep = 60; //seconds so Step of 60 s = 1 min

	double tStep_m =  tStep/60;  //number of minutes per step
  	double tStep_hr = tStep/3600;  //num hrs in a tStep  (fractional)

	double num_touch, transfer, tempDose,chooseF;

	//double current_time;

  	unsigned int SEC_PER_DAY = 24*60*60;

	unsigned int nSteps = SEC_PER_DAY/tStep;


	Fomite* cfom; //pointer to current fomite
	Person dummy; //buffer to store current person


	totalN = queue_count(pool[0].D); //total population...check later that they all made it out of the pool

	for (n=0;n<nSteps;n++)
	{
		// current_time = (double)n/nSteps;//% =  fraction of day


		if(n%(unsigned int)(1*tStep_m)==0)  //every  1 minute
		{
			//at each time step, get 'm' people from initial 'pool' and send them to first 'location'

			m = numEvents(ARRIVAL_RATE_MIN);
			//printf("\n%d",m);


			for(p=0;p<m;p++)
			{
				if(!isEmpty_queue(pool[0].D))
				{

				    dequeue(pool[0].D,&dummy);

				    nl = getNewLocation(&pool[0],0);

				    if(dummy.jumps!=0){
					dummy.jumps = dummy.jumps -1;

                       			enqueue(loc[nl].A,&dummy);
				    }

				}
			}
		}




		for(p=0;p<numLoc;p++)
		{

			/*At each location, get 'm' available touchscreens and grab as many people (if available)
			 out of 'arrival' queue to interact with them*/
			 num_fom = loc[p].numFomites;

			 double INTERACTION_RATE = loc[p].INTERACTION_RATE;
			 double QUEUE_RATE = (loc[p].foms[0].L)*num_fom;

			 if(!isEmpty_queue(loc[p].A))
			 {

				m = numEvents(QUEUE_RATE);
				//printf("\n\tm:%u",m);


				for(i=0;i<m;i++)
				{

					if(!isEmpty_queue(loc[p].A))
			 		{



					   dequeue(loc[p].A,&dummy);


						if(accept(INTERACTION_RATE))
						{
							//choose one of the available fomites to interact with at random
							chooseF = uniform();
							chooseF=ceil(chooseF*num_fom)-1;
							f=(unsigned int)chooseF;

							cfom = &(loc[p].foms[f]);//currently selected fomite


							//number of touch events for interaction with fomite
							num_touch = round(sampleNormT(cfom->use_stats));
							//printf("\n nt: %f",num_touch);

							if( dummy.status == INFECTIOUS)
                             				{

						          //% left after all touches
						          transfer  = pow( (1-sampleNormT(DEPOSIT) ) ,num_touch );
						          // % transferred
							  transfer= 1-transfer;

						          cfom->dose = dummy.fdose*transfer;
							//      dummy.fdose = dummy.fdose*(1-transfer); //removed from finger
							  dummy.fdose = DoseOnFinger(pathogen); //reload!?

							   cfom->gap_count=0;

							}else{//Susceptible or already infected

							  cfom->gap_count++;
							//% left after all touches
						          transfer  = pow( (1-sampleNormT(PICKUP) ) ,num_touch );
						          // % transferred
							  transfer= 1-transfer;

							   tempDose = dummy.fdose+ (cfom->dose)*transfer;//fdose===0
							// dummy.fdose = dummy.fdose+ (cfom->dose)*transfer;
                                  			  cfom->dose = (cfom->dose)*(1-transfer);


							 //infectious dose ?
							  if(dummy.status == SUSCEPTIBLE){

							  	//if(selfInoc(&dummy,pathogen,FACE)){
								if(selfInoc(tempDose,pathogen,FACE)){

							             dummy.status = INFECTED;
								     loc[p].NI++;
 								     dummy.gap = cfom->gap_count;

								}
							  }
							}
						}//interaction


					    enqueue(loc[p].D,&dummy); //arrivals go to departure

					 }//go to fomite- empty?
				}

			}//empty queue

			//Cleaing and pathogen half-life die-off

			for (f=0;f<num_fom;f++) //swap for queue to fom...
			{

				double CR;
				cfom = &(loc[p].foms[f]);//current fomite
				//decay rates
				  cfom->dose =  (cfom->dose)*pow( (pathogen->muf), tStep_hr);

				  CR = cfom->CR;
                       		//Cleaning ? -at 98% efficiency
                        	if( accept(CR/nSteps) )
				{
                                  cfom->dose  =  (cfom->dose)*(1-0.98);
				}

			}//for each fomite

		}//each location


		//Moving to new locations: Dpeartures from one location can go to another's arrival queue

		if(n%(unsigned int)(1*tStep_m)==0 )  //every  minutes
		{
			for(p=0;p<numLoc;p++)
			{
				unsigned long c = queue_count(loc[p].D);

				for(f=0;f<c;f++)
				{
				//if(!isEmpty_queue(loc[p].D)) //superfluous
				    dequeue(loc[p].D,&dummy);

				    if(dummy.jumps==0)//leave simualtion-> go to pool arrivals!
				    {
					enqueue(pool[0].A,&dummy);

				    }else if(hasHappenedByNow(loc[p].L, 1*tStep_m) )
				    {	//has left location to new arrival at another
					nl = getNewLocation(&loc[p],p);

					dummy.jumps = dummy.jumps -1;
					enqueue(loc[nl].A,&dummy);

				    }else{
					//stay at location: back of same queue
					enqueue(loc[p].D,&dummy);
				    }
				}
			}
		}



	}//time


	//Check that everyone made it out of the initial pool
	if(totalN != queue_count(pool[0].A))
	printf("\nWarning: sim did not get through entire pool");

	//round everybody up left to poolA
	for(p=0;p<numLoc;p++)
			{
				unsigned long c = queue_count(loc[p].D);

				for(f=0;f<c;f++)
				{
				    dequeue(loc[p].D,&dummy);
				    enqueue(pool[0].A,&dummy);
				}

				c = queue_count(loc[p].A);
				for(f=0;f<c;f++)
				{
				    dequeue(loc[p].A,&dummy);
				    enqueue(pool[0].A,&dummy);
				}

			}


}//eof function





/*MAIN SIMULATION Starting point*/


//locations
#define NUM_LOCATIONS 2
#define NUM_REALISATIONS 10000  //10000 or do 1000 x 10 for less powerful computers

#define NUM_PEOPLE 12000  //e.g. 12000 for airport simualtion, 100 for default

int main(int argc, char* argv[])
{


	/*File I/O: write to a .csv file*/
	FILE* fp;
	char filename[40];

	unsigned int i,j,k,X;


	if(argc>1)
	{
		sprintf(filename,"sim_%s.csv",argv[1]);
		fp = fopen(filename,"w");

		if (fp == NULL){
			printf("Could not open file for writing");
			return -1;
		}

		//fprintf(fp,"#x=disease_prevalence,nLocations:1,nTui/Loc:1;DEFAULT_TUI,nRealisation:10000,N:100,\n");
		//fprintf(fp,"#x=pathogen_survival_with_less Luse = 0.2,nLocations:1,nTui/Loc:1;DEFAULT_TUI,nRealisation:10000,N:100,\n");
		//fprintf(fp,"#x=pathogen_ID50,nLocations:1,nTui/Loc:1;DEFAULT_TUI,nRealisation:10000,N:100,\n");
		//fprintf(fp,"#x=numtouches,nLocations:1,nTui/Loc:1;DEFAULT_TUI with variable mode,nRealisation:10000,N:100,\n");
		//fprintf(fp,"#x=numFomitesperLocation,nLocations:1,nTui/Loc:variable;DEFAULT_TUI,nRealisation:10000,N:100,\n");
		//fprintf(fp,"#x=cleaningRate(daily),nLocations:1,nTui/Loc:1;DEFAULT_TUI,nRealisation:10000,N:100,\n");
	//	fprintf(fp,"#x=add locations,nLocations:1,nTui/Loc:1;DEFAULT_TUI,nRealisation:10000,N:100,\n");
		//fprintf(fp,"#AirportSim x=cleaning,nLocations:2,nTui/Loc:36,24;Check-in,Bagdrop,nRealisation:10000,N:12000,\n");
		fprintf(fp,"#AirportSim x=replace_touchless,nLocations:2,nTui/Loc:36,24;Check-in,Bagdrop,nRealisation:10000,N:12000,\n");

	}else{
		printf("\nRequired: output file name (example \"mydata.csv\") \n");
		return 0;
	}

	fprintf(fp,"x,R,CI,GAP");

	for(i=0;i<NUM_LOCATIONS;i++)
		fprintf(fp,",L%d",i);

	fprintf(fp,"\n");

	srand(time(0)); //seed random number generator





	 //Finger-to-fomite
	tgParam DEPOSIT;
  	DEPOSIT.md = 0.05;
  	DEPOSIT.sd = 0.1;
  	DEPOSIT.mn = 0;
  	DEPOSIT.mx = 0.6;


	 //Fomite-to-finger
	tgParam PICKUP;
	PICKUP.md = 0.20;
	PICKUP.sd = 0.2;
	PICKUP.mn = 0;
	PICKUP.mx = 0.6;

	//finger-to-face
	tgParam FACE;
	FACE.md = 0.35;
	FACE.sd = 0.1;
	FACE.mn = 0;
	FACE.mx = 1;

	//fomite interfaces

	tgParam BAG_DROP;
	BAG_DROP.md = 5;
	BAG_DROP.sd = 1.5;
	BAG_DROP.mn = 4;
	BAG_DROP.mx = 8;

	tgParam CHECK_IN;
	CHECK_IN.md = 10;
	CHECK_IN.sd = 3.5;
	CHECK_IN.mn = 8;
	CHECK_IN.mx = 30;

	tgParam TOUCH_LESS;
	TOUCH_LESS.md = 0;
	TOUCH_LESS.sd = 0;
	TOUCH_LESS.mn = 0;
	TOUCH_LESS.mx = 0;



	/*tgParam MCD_MENU;
	MCD_MENU.md = 12;
	MCD_MENU.sd = 3.5;
	MCD_MENU.mn = 10;
	MCD_MENU.mx = 40;
	*/

	//default TUI
	tgParam TOUCH_SCREEN_1;
	TOUCH_SCREEN_1.md = 5;
	TOUCH_SCREEN_1.sd = 1.5;
	TOUCH_SCREEN_1.mn = 4;
	TOUCH_SCREEN_1.mx = 30;


	//Default pathogen parameters
	Pathogen pathogen;
	pathogen.shed[0] = pow(10,4);
	pathogen.shed[1] = pow(10,6);
  pathogen.ID50 = pow(10,2);//2
  	//pathogen.muh = 0.50;  //1 hr half life
  pathogen.muf = 0.80;  //~3 hr half life  0.8


	//other parameters...
 	double  PREVALENCE = 0.02;

  unsigned int N = NUM_PEOPLE;//12000;//100;//1000*12;  //number of people per day

	//rate at which people leave the intial pool and enter the simulation
	double ARRIVAL_RATE_MIN = 20;//20;  0.2 = 1 every 5(default)   //0.5 1 every 2 min



//For each parameter setting use a dummy parameter 'x' and set it appropriately to parameters further down in th code


//Replace proportion with touchless interfaces
for(X=0;X<101;X=X+20){

	double x = (double)X/100;

//disease prevalence
/*for(X=0;X<11;X++){

	double x = (double)X;

	x = pow(1.5849,x);
	x=x/100;

	PREVALENCE = x;
*/

/*//Pathogen survival
for(X=0;X<100;X=X+5){

	double x = (double)X;

	x=x/100;
	x=x*x*x;

	pathogen.muf = x;  //half life( minutes) is -0.7/log(x) *60
	//have changes Luse 0.5(default), 0.1, 0.2
*/

/*//ID50
for(X=0;X<11;X++){

	double x = (double)X;

	x=pow(10,x);

	pathogen.ID50 = x;
*/

/*//number of touches per TUI
for(X=0;X<11;X++){

	double x = (double)X;
	x=5+2*x;
	//change mu for fomite
	TOUCH_SCREEN_1.md = x;
*/

/*//number of TUIs
for(X=0;X<11;X++){

	double x = (double)X;
	//change nFomites (see inner location loops...
	x=1+2*x;
*/

//Cleaning Rate
/*for(X=0;X<12;X++){

	double x = (double)X;
	//change CR rates in Fomite setup
	x = pow(2,x)-1; //CR
*/

/*//Additional locations...output will be a single value
for(X=0;X<1;X++){

	double x = (double)NUM_LOCATIONS;
		//note: adjust jumps number..and pool/loc next locations

*/

	//Simualtion statistics output

	//avg. number of infections at each location
	double NI[NUM_LOCATIONS];

	for(i=0;i<NUM_LOCATIONS;i++)
		NI[i]=0.0;


	double R_value,CI,gap; //averaged over j- realisations..
	R_value =0;
	CI = 0;
	gap=0;

	printf("\nsimulating parameter x: %f \n", x);


//number of realisations
for(j=0;j<NUM_REALISATIONS;j++){

	//for each new realisation, reset locations,fomites, people etc...

	Person dummy;

	Location loc[NUM_LOCATIONS];

	Location pool[1];

	pool[0].D = new_queue(sizeof(Person));
	pool[0].A = new_queue(sizeof(Person));

	//load pool with people
	for(i=0;i<N;i++)
	{
		if(accept(PREVALENCE))
		{

		 dummy.status = INFECTIOUS;
			 dummy.fdose = DoseOnFinger(&pathogen);

		}else{

		 dummy.status = SUSCEPTIBLE;
		 dummy.fdose = 0.0;
		 dummy.gap = 0;

		}


		//#############################################################
		dummy.jumps = NUM_LOCATIONS;  //jumps===num locations visited
		//#############################################################

		enqueue(pool[0].D,&dummy);

	}

	//set initial destinations
	pool[0].i_next = malloc(NUM_LOCATIONS*sizeof(unsigned int));
	pool[0].p_next = malloc(NUM_LOCATIONS*sizeof(double));
	pool[0].numNext = NUM_LOCATIONS;


	for(k=0;k<NUM_LOCATIONS;k++)
	{

		//####################################################################
		//set up Pool destinations
		pool[0].i_next[k] = k;

		if(k==0)
		{pool[0].p_next[k] = 1.0; } //flow into first location only e.g. check-in machines
		else
		{pool[0].p_next[k] = 0.0;}


		//pool[0].p_next[k] = (double)1.0/(NUM_LOCATIONS-1); //equal probability to any location
		//####################################################################


		Fomite* newfom = NULL;

		//Initialise a new location
		loc[k].A = new_queue(sizeof(Person));
		loc[k].D = new_queue(sizeof(Person));

		loc[k].INTERACTION_RATE = 1.0;
		loc[k].NI=0;

		//#####################################################################
		//	e.g. loc[k].L = 1000;//1000 people per minute move on,  or
		loc[k].L =0.05;	//0.05 -> avergae time of 20 minutes once in destination..

		//#####################################################################
		loc[k].i_next = NULL;
		loc[k].p_next = NULL;

		loc[k].numNext = NUM_LOCATIONS-1;

		loc[k].i_next = malloc((NUM_LOCATIONS-1)*sizeof(unsigned int));
		loc[k].p_next = malloc((NUM_LOCATIONS-1)*sizeof(double));

		//Set up markov chain between destinations
		unsigned int count=0;
		for(i=0;i<NUM_LOCATIONS;i++)
		{
			if(i!=k)
			{
				//#####################################################################
				loc[k].i_next[count] = i;
				loc[k].p_next[count] = 1.0;//(double)1/(NUM_LOCATIONS-1);
				//loc[k].p_next[count] = (double)1.0/(NUM_LOCATIONS-1);
				//##################################################################
				//printf("\n\t%f",loc[k].p_next[count]);
				count++;
			}

		}




		loc[k].foms = NULL;

		//####################################################################
		if(k==0)
		loc[k].numFomites = 36;  //location[0] = check-in machine
		else
		loc[k].numFomites = 24;  //lcoation[1] = bag drop

		//loc[k].numFomites = 1;// default scenario
		//###################################################################


		//create some fomites at that location
		newfom = malloc( loc[k].numFomites*sizeof(Fomite));

		//for all fomites at this location
		for(i=0;i<loc[k].numFomites;i++)
		{

			//###################################################################
			newfom[i].L = 0.5;//0.5/minute >2minutes of interaction
			//#####################################################################
			newfom[i].dose = 0;
			newfom[i].gap_count =0;

			double prop = x*loc[k].numFomites; //replacing proportion of screens with 'touch-less'


			//#################################################################

			//newfom[i].use_stats = &TOUCH_SCREEN_1;   //default


			if(k==0)
			{
				newfom[i].use_stats = &CHECK_IN;

				if(i<prop)
				newfom[i].use_stats = &TOUCH_LESS;


			}else{
				newfom[i].use_stats = &BAG_DROP;
				if(i<prop)
				newfom[i].use_stats = &TOUCH_LESS;
			}


			//Cleaning Rate!!###############################################
			newfom[i].CR = 0.0;
			//newfom[i].CR = x;
			//##################################################################

		}

		loc[k].foms = newfom;

	}//eof locations


//-----------------------------------------------------

	//DATA analysis...

	unsigned long E,I,G,count; //infectious I and newly infected E
	E=0;
	I=0;
	G=0;

	double tempR,tempCI,tempGap;

	//clock_t begin = clock();

	//main loop
	fomite_sim(
			pool,
			loc,
			NUM_LOCATIONS,
			&DEPOSIT,&PICKUP,&FACE,&pathogen,
			ARRIVAL_RATE_MIN
	);

	/*clock_t end = clock();
	double time_spent =0.0;
	time_spent +=(double)(end-begin)/CLOCKS_PER_SEC;
	printf("\n time elaspsed: is %f seconds.\n", time_spent);
	*/



	//should all be in A queue of pool!
	count = queue_count(pool[0].A);

	for(i=0;i<count;i++)
		{
			peek_queue(pool[0].A,&dummy,i);


			if(dummy.status == INFECTIOUS)
			{	I++;

			}else if(dummy.status==INFECTED){
				E++;
				G+=dummy.gap; //total gap over this realisation
			}
		}



	tempCI=0.0;
	if(I!=N)//if not all infectious
	tempCI = (double)E/(N-I);

	tempR =0.0;
	if(I>0)
	tempR = (double)E/I;

	tempGap=0.0;
	if(E>0)
	tempGap =(double)G/E; //average gap per infected

	//running averages
	R_value = (R_value*j + tempR)/(1+j);
	CI = (CI*j + tempCI)/(1+j);
	gap =(gap*j + tempGap)/(1+j);  //average gap over realisations

	/*printf("\n\t Number Infectious: %lu, Infected: %lu \n", I,E);
	printf("\n\t R_value: %f, CummulativeIncidence: %f %%\n", R_value,100*CI);
	printf("\n\t Avergae Gap: %f\n", gap);
	*/


	//clean-up heap
     	for(k=0;k<NUM_LOCATIONS;k++){
		NI[k]+=loc[k].NI;

	/*	printf("\nLocation: %d [A,D]\n",k);
		stats_queue(loc[k].A);
		stats_queue(loc[k].D);
 	*/	delete_queue(loc[k].A);
		delete_queue(loc[k].D);

		free(loc[k].foms);
		free(loc[k].i_next);
		free(loc[k].p_next);

	}
	/*printf("\nPool:[A,D]\n");
	stats_queue(pool[0].A);
	stats_queue(pool[0].D);
	*/
	delete_queue(pool[0].D);
	delete_queue(pool[0].A);
	free(pool[0].i_next);
	free(pool[0].p_next);



}//j -number of realisations

	//write to csv file

	fprintf(fp,"%f,%f,%f,%f",x,R_value,CI,gap);//append NI info...

	for(i=0;i<NUM_LOCATIONS;i++)
	fprintf(fp,",%f",NI[i]/NUM_REALISATIONS);//append NI info...

	fprintf(fp,"\n");

	printf("\n\t R_value: %f, CummulativeIncidence: %f%%, avg.Gap: %f\n", R_value,100*CI,gap);



}//x
	fclose(fp);
	printf("\nSimulation Complete: output stored in %s\n",filename);

	return 0;

}
