#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "r1279.h"

#define MAXREACT 2;// Maximum number of reactants
#define MAXPROD 2;// Maximum number of products
#define ELEMENTS 3;// Number of different elements ("element zero" also counts) 

/* Uncomment for debugging */
/*#define DEBUG */

typedef struct node_s {
	double tau_dat;
	struct node *next;
} NODE;
NODE *list_create(double);
NODE *list_insert_beginning(NODE *, double);
int count(NODE *);
NODE *in_middle(NODE *, int, double);
NODE *del_first_node(NODE *);
void display(NODE *);
double h0th(double);
double h1st(double, int);
double h2ndI(double, int, int);
int *rctn_order(int, int *, int **);
double compute_h(int, int *, double *, int **, int *, double, int);
double *propensity(int, int *, int **, double *, int *, int *, double, int);
double minimum(int, double *, int *, int, int *, NODE **, int *, double, int *, int *);
void indexx(int, double *, int *);


int main(void) {

/****************************************************************************************/
/* Exact stochastic simulation program of chemical systems with delays.		        */
/* This program uses Anderson Algorithm: 					        */
/* Anderson, J. Chem. Phys. 127 (2007) -https://doi.org/10.1063/1.2799998)		*/
/* To compile in Mac OS X (10.11.6), execute:						*/
/* $ gcc secs.c r1279.c nrutil.c indexx.c -o anderson_Model2_Hurdle anderson_Model2_Hurdle.c*/
/* 											*/
/* (To execute it: $ ./anderson_Model2_Hurdle < p53_data.dat )				*/
/****************************************************************************************/

 int maxrct = (int) MAXREACT;
 int maxprd = (int) MAXPROD;
 int maxelements = (int) ELEMENTS;

 int i, j, k, l, m;

 FILE *constants_file, *molec_file, *P0_file, *output_file;
 constants_file = fopen("constants.dat","r");
 molec_file = fopen("molecules.dat","r");
 P0_file = fopen("P0_Hill.dat", "r");
 output_file = fopen("out_p53.dat", "w");

 int RMAX = 20; /* Maxim number of reactions to be considered */
// int iseed = 1859794256;
 long int iseed;
 double MAXTIME = 43200.0;//Maximum simulation time, in seconds (=12h)

/* Reading p53_data.dat file (we fed it through the command line) */
 int *readdata;
 int **data;
 readdata = calloc(3+(maxrct+maxprd)*2, sizeof(int));
 data = calloc(RMAX, sizeof(int *));
 for (i=0; i<RMAX; i++) {
	data[i] = calloc(3+(maxrct+maxprd)*2, sizeof(int));
 }
 int counter=1;
 int R;
/****************************************************************************************/
/* WARNING!!! For model#2 proposed for the gene regulatory network of p53-Mdm2,	 	*/
/* propensity function for reaction P -> P + M (where P and M denote the number of 	*/
/* molecules for p53 and Mdm2 proteins, respectively) is computed taking into account	*/
/* Mdm2 activated production due to p53 proteins (probably dimers) binding to Mdm2	*/
/* promoter. In order to model this, we used a monotonically increasing Hill function.	*/
/* Parameters for this Hill function are cnt_react[3], Hill's coefficient (Hill) and    */
/* the half activation threshold (PHill0).						*/
/****************************************************************************************/
 double Hill;
 int nnum0, PHill0;
 double nnum1;
 while(fscanf(P0_file,"%i %lf\n", &nnum0, &nnum1) != EOF) {
	PHill0 = nnum0;
	Hill = nnum1;
 }
#ifdef DEBUG
printf("DEBUG: PHill0 = %i, Hill = %g\n", PHill0, Hill);
#endif

/****************************************************************************************/
/* WARNING! When changing the network configuration file, the most probable scenario is	*/
/* that the number of columns containing reactant/product labels also change (and the 	*/
/* same will happen with deltas). Hence, the line following this comment will need to 	*/
/* be modified accordingly.								*/
/* Keep in mind that the file containing the network configuration is parsed to the 	*/
/* executable through the command line:							*/
/* "./anderson < p53_data.dat"								*/
/****************************************************************************************/
 while(scanf("%i %i %i %i %i %i %i %i %i %i %i\n",&readdata[0],&readdata[1],&readdata[2],&readdata[3],&readdata[4],&readdata[5],&readdata[6],&readdata[7],&readdata[8],&readdata[9],&readdata[10])!=EOF) {
	for (i=0; i<(3+(maxrct+maxprd)*2); i++) {
		data[counter][i]=readdata[i];
#ifdef DEBUG
printf("DEBUG : readdata[%i] = %i\n", i, readdata[i]);
#endif
	}
	counter++;
 }
 /*20070905*/
 free(readdata);
 R=counter-1;
/* Note that readdata[0] = reaction number; 
 * readdata[1] = number of reactants; 
 * readdata[2] = number of products
 * Starting from readdata[3] and up to 3+readdata[1] we have reactant's labels,
 * and from readdata[3+readdata[1]+1] and up to 3+readdata[1]+1+readdata[2] we have product's labels	*/
 int *num_react; /* Vector containing the number of reactants in the i-th reaction */
 num_react = calloc(R+1, sizeof(int));

 int *num_prod;
 num_prod = calloc(R+1, sizeof(int));

 int **deltain;
 deltain = calloc(R+1, sizeof(int *));
 for (k=1; k<R+1; k++) {
	deltain[k] = calloc(num_react[k], sizeof(int));
 }
 int **deltaout;
 deltaout = calloc(R+1, sizeof(int *));
 for (k=1; k<R+1; k++) {
	deltaout[k] = calloc(num_prod[k], sizeof(int));
 }

 for (i=1; i<R+1; i++) {
	num_react[i]=data[i][1];
	num_prod[i]=data[i][2];
#ifdef DEBUG
printf("DEBUG basic: data[%i][1] = %i, data[%i][2] = %i\n", i, data[i][1], i, data[i][2]);
printf("DEBUG initial: num_react[%i] = %i, num_prod[%i] = %i\n", i, num_react[i], i, num_prod[i]);
#endif
 } 

 int **react; /* Double pointer containing reactant's labels */
 react = calloc(R+1, sizeof(int *));
 for (i=1; i<R+1; i++) {
	react[i] = calloc(num_react[i], sizeof(int));
 } /* Hence, react[i][j] = label reactant j in reaction i */

 int **prod;
 prod = calloc(R+1, sizeof(int *));
 for (j=1; j<R+1; j++) {
	prod[j] = calloc(num_prod[j], sizeof(int));
 }
 for (i=1; i<R+1; i++) {
	for (j=0; j<num_react[i]; j++) {
		react[i][j]=data[i][3+j];
	}
	for (k=0; k<num_prod[i]; k++) {
		prod[i][k]=data[i][3+maxrct+k];
	}
	for (l=0; l<num_react[i]; l++) {
		deltain[i][l]=data[i][3+maxrct+maxprd+l];
	}
	for (m=0; m<num_prod[i]; m++) {
		deltaout[i][m]=data[i][3+2*maxrct+maxprd+m];
	}
 }
 
 double *cnt_react;
 cnt_react = calloc(R+1, sizeof(double));
 cnt_react[0] = 0.0;
 int *type_delay; //Vector containing reaction type depending on delay
 type_delay = calloc(R+1, sizeof(int));
 double *delay; //Vector containing the delay (in seconds)
 delay = calloc(R+1, sizeof(double));
 int num0, num2;
 double num1, num3;
 int counter2=0;
 while(fscanf(constants_file,"%i %lf %i %lf\n", &num0, &num1, &num2, &num3) != EOF) {
	cnt_react[num0] = num1;
	type_delay[num0] = num2;
	delay[num0] = num3;
	counter2++;
 }
#ifdef DEBUG
printf("DEBUG: R = %i\n", R);
#endif
 if (counter2 != R) {
	printf("ERROR: data in files network.dat and constants.dat do not match\n");
	printf("ERROR: different number of lines\n");
	exit(-1);
 }

/* DEBUGGING: checking whether the program reads input data files properly */
#ifdef DEBUG
 printf("DEBUG: R = %i\n", R);
 for (i=1; i<R+1; i++) {
	printf("DEBUG: num_react[%i] = %i, num_prod[%i] = %i\n", i, num_react[i], i, num_prod[i]);
	for (j=0; j<num_react[i]; j++) {
		printf("DEBUG: react[%i][%i] = %i\n", i,j,react[i][j]);
	}
	for (k=0; k<num_prod[i]; k++) {
		printf("DEBUG: prod[%i][%i] = %i\n", i,k,prod[i][k]);
	}
	for (l=0; l<num_react[i]; l++) {
		printf("DEBUG: deltain[%i][%i] = %i\n", i,l,deltain[i][l]);
	}
	for (m=0; m<num_prod[i]; m++) {
		printf("DEBUG: deltaout[%i][%i] = %i\n", i,m,deltaout[i][m]);
	}
	printf("DEBUG: cnt_react[%i] = %g, type_delay[%i] = %i, delay[%i] = %g\n", i,cnt_react[i],i,type_delay[i],i,delay[i]);
 }
#endif

/* Reading the initial number of molecules for each specie */
 int *molecules;/* Initial number of molecules for each specie */
 molecules = calloc(maxelements, sizeof(int));
 int el1;
 int el2;
 int counter3 = 0;
 while(fscanf(molec_file,"%i %i\n", &el1, &el2) != EOF) {
	molecules[el1] = el2;
	counter3++;
 }
 if (counter3 != maxelements) {
	printf("ERROR: data in file molecules.dat, the number of elements doesn't match with the number of lines\n");
	exit(-1);
 }
#ifdef DEBUG
for (i=0; i< maxelements; i++) {
	printf("DEBUG molecules: molecules[%i] = %i\n", i, molecules[i]);
}
#endif

/********************************************************************************************************/
/* Generating vector structures s[i] that will allow us to follow which reaction fires next, taking 	*/
/* into account delayed reactions).									*/
/********************************************************************************************************/
 NODE *s[R+1];
 NODE *temp;
 int *delay_list;
 int howmany=0;
 for (i=1; i <= R; i++) {
	if (type_delay[i] != 0) {
		s[i] = list_create(1.0e19);
		howmany++;
#ifdef DEBUG
printf("DEBUG: s[%i]\n", i);
display(s[i]);
#endif
	}
 }

 if (howmany != 0) {
	delay_list = calloc(howmany+1, sizeof(int));
	j=1;
	for (i=1; i <= R; i++) {
		if (type_delay[i] != 0) {
			delay_list[j]=i;
			j++;
		}
	}
 }
 else {
	printf("\n\nNone of the reactions present delays.\n");
	delay_list = calloc(1, sizeof(int));
	delay_list[0] = 0;
/*	exit(-1);*/
 }

/* Checking everything is fine									*/
 int links; //Number of members in the chain
#ifdef DEBUG
 for (i=1; i <= howmany; i++) {
	links=count(s[delay_list[i]]);
	printf("DEBUG: s[%i] links = %i\n", delay_list[i], links);
 }
 for (i=1; i <= howmany; i++) printf("DEBUG: delay_list[%i] = %i\n", i, delay_list[i]);
#endif
/* Set P[i]=0 and T[i]=0										*/
 double *P, *T;
 P = calloc(R+1, sizeof(double));
 T = calloc(R+1, sizeof(double));
 for (i=0; i <= R; i++) {
	P[i] = 0.0;
	T[i] = 0.0;
 }

/* Defining time                                                                			*/
 double time = 0.0;
 setseed(&iseed); //Generating a random seed for the random number generator
 setr1279(iseed); //Initializing random number generator with a random seed 
#ifdef DEBUG
printf("DEBUG: iseed = %li\n", iseed);
#endif
 fprintf(output_file, "%10.10f %i %i\n", time, molecules[1], molecules[2]);//Data corresponding to time 0
 double thurdle=0.0; //Dummy time used to print data to output file given a fixed time interval (dhurdle) (in seconds)
 double dhurdle=30.0; //Time increment for printing data to output file amb el que printo les dades (in seconds)
 thurdle = thurdle + dhurdle;
 
/********************************************************************************************************/
/* Computing propensity functions (hazard) for each reaction						*/
/********************************************************************************************************/
 double *h; //Vector contaning current value of propensity function for each reaction
 int *type;
 type = rctn_order(R, num_react, react); //Vector contaning the order of each reaction 
 h = propensity(R, num_react, react, cnt_react, molecules, type, Hill, PHill0);
#ifdef DEBUG
for (i=1; i <= R; i++) printf("DEBUG propensity: type[%i] = %i, h[%i] = %g\n", i, type[i], i, h[i]);
#endif

/* (17/10/2008) If generated random number = 0 => log(0) -> -inf !!! Modifying r1279.c code to exclude zeros */
 double rndm;

/* Generating R independent, uniform (0,1) random numbers and setting P[i]=log(1.0/r[i])=-1.0*log(r[i]) */
 for (i=1; i <= R; i++) {
	rndm=r1279();
	while ( rndm == 0.0 ) rndm=r1279();
	P[i] = -1.0*log(rndm);
#ifdef DEBUG
printf("DEBUG: P[%i] = %g\n", i, P[i]);
#endif
 }
 double *deltat;
 deltat = calloc(R+1, sizeof(double));
 double mintime;
 int *rctnexec;
 int exec;//Executed reaction
 rctnexec = &exec;
 int *indxr, *indxr2;
 indxr = calloc(R+1, sizeof(int));
 indxr2 = calloc(howmany+1, sizeof(int));
 int *cmpltn;
 int completion; //If completion = 0 we are dealing with an initialization, if completion = 1 we have a completion
 cmpltn = &completion;
 NODE *dummy;

 while (thurdle <= MAXTIME) {
	for (i=1; i <= R; i++) {
		if (h[i] != 0.0) {
			deltat[i] = (P[i] - T[i])/h[i];
		}
		else {
			deltat[i] = 1.0e19;
		}
 	}
#ifdef DEBUG
for (i=1; i <= R; i++) printf("DEBUG: deltat[%i] = %g\n", i, deltat[i]);
#endif

	mintime = minimum(R, deltat, indxr, howmany, delay_list, s, indxr2, time, rctnexec, cmpltn);
#ifdef DEBUG
printf("DEBUG: mintime = %g, exec = %i, completion = %i\n", mintime, exec, completion);
#endif
	time = time + mintime;

/* Checking which kind of reaction has initiated/completed: a non-delayed reaction (ND), with Completion Delay	*/
/* (CD) or with Initiation and Completion Delay (ICD).								*/
	if (completion == 1) {//We have a finalization (completion)
#ifdef DEBUG
printf("DEBUG: We have a completion\n");
#endif
		if (type_delay[exec] == 1) {//Completion Delay (CD)
			for (i=0; i < num_react[exec]; i++) {
				molecules[react[exec][i]] = molecules[react[exec][i]]+deltain[exec][i];
			}
			for (i=0; i < num_prod[exec]; i++) {
	 			molecules[prod[exec][i]] = molecules[prod[exec][i]]+deltaout[exec][i];
			} 
		}
		if (type_delay[exec] == 2) {//Initiation and Completion Delay (ICD)
			for (i=0; i < num_prod[exec]; i++) {
	 			molecules[prod[exec][i]] = molecules[prod[exec][i]]+deltaout[exec][i];
			}
		}
		//DELETE THE FIRST ROW OF S[EXEC]
		dummy=del_first_node(s[exec]);
		s[exec]=dummy;
#ifdef DEBUG
printf("DEBUG: Esborrem el primer node\n");
display(s[exec]);
#endif
	}
	if (completion == 0) {//We have an initialization
#ifdef DEBUG
printf("DEBUG: We have an initialization\n");
#endif
		if (type_delay[exec] == 0) {//No Delay (ND)
			for (i=0; i < num_react[exec]; i++) {
				molecules[react[exec][i]] = molecules[react[exec][i]]+deltain[exec][i];
			}
			for (i=0; i < num_prod[exec]; i++) {
	 			molecules[prod[exec][i]] = molecules[prod[exec][i]]+deltaout[exec][i];
			}
		}
		else {
			links=count(s[exec]);
#ifdef DEBUG
printf("DEBUG: s[%i] links = %i\n", exec, links);
#endif
			if (type_delay[exec] == 1) {//CD
		//UPDATE S[EXEC] BY INSERTING time+tau_delay INTO S[EXEC] IN THE SECOND TO LAST POSITION 
			   if (links == 1) {
#ifdef DEBUG
printf("DEBUG: Entering if CD for links = 1\n");
#endif
				dummy=list_insert_beginning(s[exec], time+delay[exec]);
				s[exec]=dummy;
#ifdef DEBUG
printf("DEBUG: Adding %g at the beginning of the list , links = %i\n", time+delay[exec], links);
display(s[exec]);
#endif
			   } else {
#ifdef DEBUG
printf("DEBUG: Entering if CD for links != 1\n");
printf("DEBUG: links-1 = %i\n", links-1);
#endif
				dummy=in_middle(s[exec],links-1,time+delay[exec]);
#ifdef DEBUG
printf("DEBUG: Adding %g after position %i\n", time+delay[exec], links-1);
display(s[exec]);
#endif
			   }
			}
			if (type_delay[exec] == 2) {//ICD
				for (i=0; i < num_react[exec]; i++) {
					molecules[react[exec][i]] = molecules[react[exec][i]]+deltain[exec][i];
				}
		//UPDATE S[EXEC] BY INSERTING time+tau_delay INTO S[EXEC] IN THE SECOND TO LAST POSITION (PENULTIMA POS.)
			   	if (links == 1) {
#ifdef DEBUG
printf("DEBUG: Entering if ICD for links = 1\n");
#endif
					dummy=list_insert_beginning(s[exec], time+delay[exec]);
					s[exec]=dummy;
#ifdef DEBUG
printf("DEBUG: Adding %g at the beginning of the list, links = %i\n", time+delay[exec], links);
display(s[exec]);
#endif
			  	} else {
#ifdef DEBUG
printf("DEBUG: Entering if ICD for links != 1\n");
printf("DEBUG: links-1 = %i\n", links-1);
#endif
					dummy=in_middle(s[exec],links-1,time+delay[exec]);
#ifdef DEBUG
printf("DEBUG: Adding %g after position %i\n", time+delay[exec], links-1);
display(s[exec]);
#endif
			   	}
			}
		}
	}

	for (i=1; i <= R; i++) {
		T[i] = T[i] + h[i]*mintime;
	}
	if (completion == 0) {
		rndm=r1279(); 
		while ( rndm == 0.0 ) rndm=r1279();//Original r1279.c code was modified to avoid zeros when sorting random numbers; keeping this just in case r1279.c versions get messed up
		P[exec] = P[exec]-log(rndm);
	}
 
/* Recomputing propensity functions, h[i]				*/
	h = propensity(R, num_react, react, cnt_react, molecules, type, Hill, PHill0);
#ifdef DEBUG
for (i=1; i <= R; i++) printf("DEBUG recomp. propensity: type[%i] = %i, h[%i] = %g\n", i, type[i], i, h[i]);
#endif

/* Printing output (time, number of molecules for p53, number of molecules for Mdm2)	*/
 	if (time > thurdle) {
        	fprintf(output_file, "%10.2f %i %i\n", time, molecules[1], molecules[2]);
		thurdle = thurdle + dhurdle;
 	}

 }

 fclose(constants_file);
 fclose(molec_file);
 fclose(output_file);

 return 0;

}

NODE *list_create(double tau_dat) {
/* Subroutine to create the initial element of the list					*/
	NODE *node;
	if (!(node=malloc(sizeof(NODE)))) return NULL;
	node->tau_dat=tau_dat;
	node->next=NULL;
	return node;
}

NODE *list_insert_beginning(NODE *list, double tau_dat) {
/* Subroutine to insert a node with value "tau_dat" on the first position in the list	*/
	NODE *newnode;
	newnode=list_create(tau_dat);
	newnode->next = list;
	return newnode;
}

int count(NODE *list) {
/* Subroutine that counts the number of nodes in the list				*/
	int c=0;
	while (list!=NULL) {
		c++;
		list=list->next;
	}
	return c;
}

NODE *in_middle(NODE *node, int loc, double tau_dat) {
/* Subroutine that inserts a node with value "tau_dat" AFTER position "loc" in the list 	*/
/* Returns the first node in the list.								*/
/* INPUT:	*node: first node in the list							*/
/*		loc: position on which we want to insert the new node				*/
/*		tau_dat: value that the inserted node will contain				*/
/* OUTPUT:	*node: first node in the list							*/
	NODE *temp, *n;
	int c=1, flag=0;
	temp=node;
	if(node==NULL) {
		printf("\n\nLink List Is Empty. Can't Insert.\n");
		exit(-1);
	}
	else {
	  while(temp!=NULL) {
		if (c==loc) {
			n = (NODE *)malloc(sizeof(NODE));
			n->tau_dat=tau_dat;
			n->next=temp->next;
			temp->next=n;
			flag=1;
		}
		c++;
		temp=temp->next;
	  }
	  if (flag==0) {
		printf("\n\nNode Specified Doesn't Exist. Can't Enter The Data.\n");
		exit(-1);
	  }
	}
	return node;

}

NODE *del_first_node(NODE *node) {
/* Subroutine that removes the first node in the list and returns the "new" first node 		*/
	if (node == NULL) {
		printf("\n\nEmpty Linked List. Can't Delete The Data.\n");
		exit(-1);
	}
	else {
		NODE *shift;
		shift=node->next;
#ifdef DEBUG
printf("DEBUG: shift->tau_dat = %g\n", shift->tau_dat);
#endif
		free(node);
		return shift;
	}
}


void display(NODE *list) {
	if (list==NULL) printf("\n\nEmpty Link List. Can't Display The Data.\n\n");
	while(list!=NULL) {
		printf("%g\n",list->tau_dat);
		list=list->next;
	}
}

double h0th(double constant) {
/* Function that computes the propensity function for an order zero reaction 				*/
	double value;
	value = constant;
	return value;
}

double h1st(double constant, int molec1) {
/* Function computing the propensity function for a first order reaction				*/
	double value;
	value = constant*(double)molec1;
	return value;
}

double h2ndI(double constant, int molec1, int molec2) {
/* Function that computes the propensity function for a second order reaction (type I)			*/
	double value;
	value = constant*(double)molec1*(double)molec2;
	return value;
}

int *rctn_order(int R, int *num_react, int **react) {
/* Function that returns array type[1..R] that contains the order of each reaction. 		*/
/* If type[i] = 0 => zero order reaction; if type[i] = 1 => first order reaction;		*/
/* If type[i] = 2 => second order reaction (type I).						*/

	int i;
	int *type;
	type = calloc(R+1, sizeof(int));
	type[0] = 0;
	for (i=1; i<R+1; i++) {
		if ((num_react[i] == 1) && (react[i][0] !=0)) {
			type[i] = 1;
		} else if ((num_react[i] == 1) && (react[i][0] ==0)) {
			type[i] = 0;
		} else if (num_react[i] == 2) {
			type[i] = 2;
		} 
		else {
			printf("Reaction order is different than zero, first, and second order (type I).\n");
			exit(-1);
		}
	}
	return type;
#ifdef DEBUG
for (i=1; i<R+1; i++) printf("DEBUG reaction type: type[%i] = %i\n", i, type[i]);
#endif
}

double compute_h(int rctn, int *type, double *cnt_react, int **react, int *molecules, double Hill, int PHill0) {
/* Function that computes propensity function value for a given reaction			*/
/* Input:	rctn = num. that identifies the reaction for which we want to compute h (h[i])	*/
/*		type[1..R] = vector containing the order of each reaction (0, 1 or 2 type I)	*/
/*		cnt_react[1..R] = vector containing rate constants for each reaction		*/
/*		react[1..R][0,1] = array containing reactant labels for each reaction		*/
/*		molecules[1..MAXREACT] = vector containing the number of molecules present in 	*/
/*					 the system for each element				*/
/*		Hill = Hill function coefficient (for reaction 3)				*/
/*		PHill0 = p53 threshold  (half activation Hill function threshold, for reaction 3)*/

	double value, dummy;
	double (*hazard)();

	if (rctn == 3) {//Computing propensity function according to the monotonically increasing Hill function
		dummy = (double)PHill0/(double)(molecules[react[rctn][0]]);
#ifdef DEBUG
printf("DEBUG: dummy = %g\n", dummy);
#endif
		value = cnt_react[rctn]/(1.0+pow(dummy,Hill));
	} else {
		if (type[rctn] == 0) {
			hazard = h0th;
			value = hazard(cnt_react[rctn]);
		}
		else if (type[rctn] == 1) {
			hazard = h1st;
			value = hazard(cnt_react[rctn], molecules[react[rctn][0]]);
		}
		else if (type[rctn] == 2) {
			hazard = h2ndI;
			value = hazard(cnt_react[rctn], molecules[react[rctn][0]], molecules[react[rctn][1]]);
		}
		else {
			printf("ERROR: order of reaction different from 1st and 2ndI.\n");
			exit(-1);
		}
	}
	return value;
}

double *propensity(int R, int *num_react, int **react, double *cnt_react, int *molecules, int *type, double Hill, int PHill0) {
/* Function returning propensity function values for each reaction 				*/
/* Input:	R = num. of reactions in the network						*/
/*		num_react[1..R] = num. of reactants for each reaction				*/
/*		react[1..R][num_react[1..R]] = reactant labels for each reaction 		*/
/*		cnt_react[1..R] = rate constants for each reaction				*/
/*		molecules[react[1..R][j] = num. of molecules for each reactant			*/
/*		Hill = Hill function coefficient (for reaction 3)				*/
/*		PHill0 = p53 threshold  (half activation Hill function threshold, for reaction 3)*/
/* Output:	h[1..R] = vector containing propensity function values for each reaction	*/

/* Labelling each reaction to identify whether they are order zero, first order, 	*/
/* second order type I or second order type II.						*/

	int i;
	double aux;
	double *h;
	h = calloc(R+1, sizeof(double));
	h[0] = 0.0;
/* Computing propensity function h[i] (where i is the reaction number)			*/
	for (i=1; i <= R; i++) {
		aux = compute_h(i, type, cnt_react, react, molecules, Hill, PHill0);
		h[i] = aux;
	}
	return h;
}

double minimum(int R, double *delta, int *indxr, int count, int *del_list, NODE **node, int *indxr2, double t, int *rctnexec, int *cmpltn) {
/* Subroutine computing the minimum value of array delta[1,..,R], the minimum value of node[]->tau_dat, */
/* then computes the value of (node[]->tau_dat)_min - t and compares this valude with (delta[])_min	*/
/* returning the minimum value among these two, and returning also (using pointer rctnexec the number	*/
/* of the executed reaction).										*/
/* INPUT:	R = number of reactions in the network 							*/
/*		delta[1,..,R] = vector containing time increments for each reaction 			*/
/*		indxr[1,..,R] = vector containing the ordered list of reactions, from minor to major 	*/
/*				delta									*/
/*		count = tells us how many reactions are delayed						*/
/*		del_list[1,..,count] = list of delayed reactions in the network				*/
/*		node[1,..,R] = structure (singly linked list) containing delays	of reactions		*/
/*		indxr2[] = vector containing the ordered list of delayed reactions, from low to 	*/
/*			   high node (delay)								*/
/*		t = time										*/
/*		rctnexec = pointer returning the number of executed reaction 				*/
/*		cmpltn = pointer returning 0 or 1 depending on whether the fired reaction is an 	*/
/*			 initialization or a finalization						*/
/* OUTPUT:	value = min{delta[], (node[]->tau_dat)-t }						*/	
	int i, j, flag;
	double value, val1, preval2, val2;
	indexx(R, delta, indxr);
	val1 = delta[indxr[1]]; //Minimum value of delta vector
	if (count > 0) {
		double *vec2;
		vec2 = calloc(count+1, sizeof(double));
		for (i = 1; i <= count; i++) {
			vec2[i] = node[del_list[i]]->tau_dat;
		}
		indexx(count, vec2, indxr2);
		preval2 = vec2[indxr2[1]]; //Minimum value of delays
		val2 = preval2 - t;
		if (val1 < val2) {//Initiating a reaction
			*rctnexec = indxr[1];
			*cmpltn = 0;
			return val1;
		}
		else if (val2 < val1) {//Completing a reaction
	   	  flag = 0;
	   	  for(i=1; i <= count; i++) {
			if (node[del_list[i]]->tau_dat == preval2) {
				*rctnexec = del_list[i];
				flag = 1;
				*cmpltn = 1;
				return val2;
			} 
            	  }
	   	  if (flag == 0) {
			printf("\n\nError. Reaction Index Not Found.\n");
			exit(-1);
	    	  }

		}
		else {
			printf("val1 = %g, val2 = %g, which reaction should I choose!?\n", val1, val2);
			exit(-1);
		}
	} else {//None of the reactions present delay => there's no need for a vector structure s[j]
		*rctnexec = indxr[1];//Initiating a reaction 
		*cmpltn = 0;
		return val1;

	}


}
