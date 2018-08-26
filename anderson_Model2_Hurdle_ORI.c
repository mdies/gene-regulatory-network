#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/*#include <malloc/malloc.h>*/
#include <math.h>
#include "r1279.h"

#define MAXREACT 2;// Num maxim de reactants
#define MAXPROD 2;// Num maxim de productes
#define ELEMENTS 3;// Num d'elements diferents totals (comptant el zero)

/* Uncomment for debugging */
#define DEBUG 

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
/* Programa de simulacio estocastica de sistemes quimics amb retards.			*/
/* Algoritme d'Anderson.								*/
/* Per compilar el programa:								*/
/* $ gcc -lm secs.c r1279.c nrutil.c indexx.c -o anderson /usr/lib/libm.a anderson.c	*/
/* 											*/
/* (per executar-lo: $ ./anderson < p53_data.dat )					*/
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

 int RMAX = 20; /* Maxim nombre de reaccions a considerar */
// int iseed = 1859794256;
 long int iseed;
// double MAXTIME = 61200.0;//Temps maxim de la simulacio, en segons (17 hores)
 double MAXTIME = 10000.0;//Temps maxim de la simulacio, en segons 
// double MAXTIME = 7200000.0; // \Omega <=100
// double MAXTIME = 2160000.0; // \Omega=200, 500
// double MAXTIME = 1080000.0; // \Omega=1000

/* Llegim el fitxer p53_data.dat que li hem passat per la linia de comandes */
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
/* ATENCIO!!! Per al model 2 proposat per a la xarxa genetica de p53-Mdm2, la		*/
/* "propensity function" de la reaccio P -> P + M es calcula tenint en compte		*/
/* l'activacio de la produccio de Mdm2 donada per la unio de proteines p53 (dimers	*/
/* segurament) a la regio del promotor de Mdm2 del DNA. Aixo es fa mitjançant una 	*/
/* funcio de Hill monotonicament creixent. Els parametres d'aquesta funcio de Hill son	*/
/* cnt_react[3], el coeficient de Hill (Hill) i PHill0 (el "threshold").		*/
/****************************************************************************************/
// double Hill=4.0;
 double Hill;
// int PHill0 = 9;
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
/* Quan canvies les dades del fitxer de configuracio de la xarxa, el mes probable 	*/
/* es que tambe canviin les longituds de les columnes de les etiquetes de reactants	*/
/* i/o de productes (i per tant tambe les deltes); aixo implica que cal que refacis	*/
/* la linia de lectura que segueix a aquest comentari.				   	*/
/* Recorda que el fitxer de configuracio de la xarxa el passem per teclat:		*/
/* "./anderson < p53_data.dat"								*/
/****************************************************************************************/
 while(scanf("%i %i %i %i %i %i %i %i %i %i %i\n",&readdata[0],&readdata[1],&readdata[2],&readdata[3],&readdata[4],&readdata[5],&readdata[6],&readdata[7],&readdata[8],&readdata[9],&readdata[10])!=EOF) {
	for (i=0; i<(3+(maxrct+maxprd)*2); i++) {
		data[counter][i]=readdata[i];
#ifdef DEBUG
printf("DEBUG lectura dades: readdata[%i] = %i\n", i, readdata[i]);
#endif
	}
	counter++;
 }
 /*20070905*/
 free(readdata);
 R=counter-1;
/* Sabem que readdata[0] = numero de la reaccio; readdata[1] = numero de reactants; readdata[2] = numero de productes	*/
/* A partir de readdata[3] tenim les etiquetes dels reactants, i fins a 3+readdata[1], a readdata[3+readdata[1]+1] tenim*/
/* les etiquetes dels productes, i fins a 3+readdata[1]+1+readdata[2] 							*/
 int *num_react; /* Vector que contindra el num de reactants de la reaccio i-essima */
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
printf("DEBUG inicial: num_react[%i] = %i, num_prod[%i] = %i\n", i, num_react[i], i, num_prod[i]);
#endif
 } 

 int **react; /* Doble punter que contindra les etiquetes dels reactants */
 react = calloc(R+1, sizeof(int *));
 for (i=1; i<R+1; i++) {
	react[i] = calloc(num_react[i], sizeof(int));
 } /* Aixi, react[i][j] = label reactant j in reaction i */

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
 /*20070905*/
/*
 crec que sobra un dels 3 free
 int s = 0, st = 0;
 int longTmp = 0;
 for(s = 0; s < sizeof(data); s++) {
   longTmp = sizeof(data[s]);
//   for(st = 0; st < sizeof(longTmp); st++) {
//     free(data[s][st]);
   }
   free(data[s]);
 }
 free(data)*/

 double *cnt_react;
 cnt_react = calloc(R+1, sizeof(double));
 cnt_react[0] = 0.0;
 int *type_delay; //Vector que contindra el tipus de reaccio segons el delay
 type_delay = calloc(R+1, sizeof(int));
 double *delay; //Vector que contindra el delay (en segons)
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

/* DEBUGEM: mirem si el programa llegeix els fitxers de dades correctament */
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

/* Finalment ja nomes cal llegir el numero inicial de molecules per cada especie */
 int *molecules;/* Numero de molecules de cada especie */
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
/* Generem l"estructura de vectors" s[i]								*/
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
			delay_list[j]=i;//TODO: definir CD_rct[i] (contindra les etiquetes de les reaccions
			j++;		//amb Completion Delay; ICD_rct[i] (contindra les etiquetes de les
		}			//reaccions amb Initiation and Completion Delay)
	}
 }
 else {
	printf("\n\nNone of the reactions present delays.\n");
	delay_list = calloc(1, sizeof(int));
	delay_list[0] = 0;
/*	exit(-1);*/
 }

/* Comprovem que tot es correcte									*/
 int links; //Numero de membres de la cadena
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

/* Definim el temps                                                               			*/
 double time = 0.0;
 setseed(&iseed); //Generem una llavor a l'atzar
 setr1279(iseed); //Inicialitzem el generador de nombres aleatoris amb la llavor obtinguda
#ifdef DEBUG
printf("DEBUG: iseed = %li\n", iseed);
#endif
 fprintf(output_file, "%10.10f %i %i\n", time, molecules[1], molecules[2]);//Dades corresponents a temps 0
 double thurdle=0.0; //Temps que fixare per a printar les dades cada cert interval de temps (en segons)
 double dhurdle=30.0; //Increment de temps amb el que printo les dades (en segons)
 thurdle = thurdle + dhurdle;
 
/********************************************************************************************************/
/* Calculem les "propensity functions" (hazard) per a cadascuna de les reaccions			*/
/********************************************************************************************************/
 double *h; //Contindra el valor "actual" de la "propensity function" per a cada reaccio
 int *type;
 type = rctn_order(R, num_react, react); //Vector que conte l'ordre de cada reaccio
 h = propensity(R, num_react, react, cnt_react, molecules, type, Hill, PHill0);
#ifdef DEBUG
for (i=1; i <= R; i++) printf("DEBUG propensity: type[%i] = %i, h[%i] = %g\n", i, type[i], i, h[i]);
#endif

/* (17/10/08) Si el num. aleatori sortejat = 0 => log(0) -> -inf !!! Excloem aquest num. doncs	*/
 double rndm;

/* Generate R independent, uniform (0,1) random numbers and set P[i]=log(1.0/r[i])=-1.0*log(r[i])*/
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
 int exec;//Reaccio executada
 rctnexec = &exec;
 int *indxr, *indxr2;
 indxr = calloc(R+1, sizeof(int));
 indxr2 = calloc(howmany+1, sizeof(int));
 int *cmpltn;
 int completion; //Si completion = 0 es tracta d'una inicialitzacio, si completion = 1 es tracta d'una completacio
 cmpltn = &completion;
 NODE *dummy;

/* while (time <= MAXTIME) {*/
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

/* Necessito crear una subrutina minimum(int R, double *delta, NODE *node[k], double t, int *reaction) que	*/
/* em retorni el valor minim de delta[1,..,R], node[1,..,howmany]-t ; i que em retoni	*/
/* tambe la reaccio que s'ha executat: numero de reaccio (mitjançant un apuntador).	*/
	mintime = minimum(R, deltat, indxr, howmany, delay_list, s, indxr2, time, rctnexec, cmpltn);
#ifdef DEBUG
printf("DEBUG: mintime = %g, exec = %i, completion = %i\n", mintime, exec, completion);
#endif
	time = time + mintime;

/* Hem de veure ara quin tipus de reaccio hem iniciat/completat: una reaccio sense retard (ND), amb Completion	*/
/* Delay (CD) o be amb Initiation and Completion Delay (ICD).							*/
	if (completion == 1) {//Tenim una finalitzacio (completacio)
#ifdef DEBUG
printf("DEBUG: Tenim una completacio\n");
#endif
		if (type_delay[exec] == 1) {//Completion Delay (CD)
			for (i=0; i < num_react[exec]; i++) {//TODO: FER UNA FUNCIO UPDATE_REACT
				molecules[react[exec][i]] = molecules[react[exec][i]]+deltain[exec][i];
			}
			for (i=0; i < num_prod[exec]; i++) {//TODO: FER UNA FUNCIO UPDATE_PROD
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
	if (completion == 0) {//Tenim una inicialitzacio
#ifdef DEBUG
printf("DEBUG: Tenim una inicialitzacio\n");
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
		//UPDATE S[EXEC] BY INSERTING time+tau_delay INTO S[EXEC] IN THE SECOND TO LAST POSITION (PENULTIMA POS.)
			   if (links == 1) {
#ifdef DEBUG
printf("DEBUG: Entro a l'if CD de links = 1\n");
#endif
				dummy=list_insert_beginning(s[exec], time+delay[exec]);
				s[exec]=dummy;
#ifdef DEBUG
printf("DEBUG: Afegim %g al principi de la llista, links = %i\n", time+delay[exec], links);
display(s[exec]);
#endif
			   } else {
#ifdef DEBUG
printf("DEBUG: Entro a l'if CD de links != 1\n");
printf("DEBUG: links-1 = %i\n", links-1);
#endif
				dummy=in_middle(s[exec],links-1,time+delay[exec]);
#ifdef DEBUG
printf("DEBUG: Afegim %g despres de la posicio %i\n", time+delay[exec], links-1);
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
printf("DEBUG: Entro a l'if ICD de links = 1\n");
#endif
					dummy=list_insert_beginning(s[exec], time+delay[exec]);
					s[exec]=dummy;
#ifdef DEBUG
printf("DEBUG: Afegim %g al principi de la llista, links = %i\n", time+delay[exec], links);
display(s[exec]);
#endif
			  	} else {
#ifdef DEBUG
printf("DEBUG: Entro a l'if ICD de links != 1\n");
printf("DEBUG: links-1 = %i\n", links-1);
#endif
					dummy=in_middle(s[exec],links-1,time+delay[exec]);
#ifdef DEBUG
printf("DEBUG: Afegim %g despres de la posicio %i\n", time+delay[exec], links-1);
display(s[exec]);
#endif
			   	}
			}
		}
	}

	for (i=1; i <= R; i++) {
		T[i] = T[i] + h[i]*mintime;
	}
	if (completion == 0) {//PREGUNTA, POTS POSAR-HO EN EL PRIMER IF? ES PER NO HAVER DE FER-LO DUES VEGADES!!!
		rndm=r1279();
		while ( rndm == 0.0 ) rndm=r1279();
		P[exec] = P[exec]-log(rndm);
	}
 
/* Recalculem les "propensity functions", h[i]				*/
	h = propensity(R, num_react, react, cnt_react, molecules, type, Hill, PHill0);
#ifdef DEBUG
for (i=1; i <= R; i++) printf("DEBUG recalc. propensity: type[%i] = %i, h[%i] = %g\n", i, type[i], i, h[i]);
#endif

/* Printem la sortida de dades (temps, num. molec. p53, num. molec. Mdm2)	*/
 	if (time > thurdle) {
        	fprintf(output_file, "%10.2f %i %i\n", thurdle, molecules[1], molecules[2]);
		thurdle = thurdle + dhurdle;
 	}

 }

 fclose(constants_file);
 fclose(molec_file);
 fclose(output_file);

 return 0;

}

NODE *list_create(double tau_dat) {
/* Subrutina que crea l element inicial de la llista		*/
	NODE *node;
	if (!(node=malloc(sizeof(NODE)))) return NULL;
	node->tau_dat=tau_dat;
	node->next=NULL;
	return node;
}

NODE *list_insert_beginning(NODE *list, double tau_dat) {
/* Subrutina que inserta un node amb valor "tau_dat" en el primer lloc de la llista	*/
	NODE *newnode;
	newnode=list_create(tau_dat);
	newnode->next = list;
	return newnode;
}

int count(NODE *list) {
/* Subrutina que compta el numero de nodes de la llista	*/
	int c=0;
//	while (list->next && list->next!=NULL) {
	while (list!=NULL) {
		c++;
		list=list->next;
	}
	return c;
}

NODE *in_middle(NODE *node, int loc, double tau_dat) {//TODO: FER-LA VOID!!!
//void in_middle(NODE *q, int loc, double tau_dat){
/* Subrutina que inserta un node amb valor "tau_dat" DESPRES de la posicio 	*/
/* "loc" de la llista. Retorna el primer node de la llista.			*/
/* INPUT:	*node: primer node de la llista				*/
/*		loc: posicio on volem insertar el nou node			*/
/*		tau_dat: valor que contindra el node insertat			*/
/* OUTPUT:	*node: primer node de la llista				*/
	NODE *temp, *n;
	int c=1, flag=0;
	temp=node;
	if(node==NULL) {
		printf("\n\nLink List Is Empty. Can't Insert.\n");
		exit(-1);
	}
	else {
//	  while(temp!=NULL && loc != 1) {//Amb la segona condicio no podriem insertar nodes en segona posicio!!!
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
/* Subrutina que esborra el primer node de la llista, i retorna el "nou"	*/
/* primer node.								*/
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
/* Funcio que calcula la "propensity function" per a una reaccio d'ordre zero				*/
	double value;
	value = constant;
	return value;
}

/*double h1st(double constant, int molec1, int molec2) {*/
double h1st(double constant, int molec1) {
/* Funcio que calcula la "propensity function" per a una reaccio de primer ordre			*/
	double value;
	value = constant*(double)molec1;
	return value;
}

double h2ndI(double constant, int molec1, int molec2) {
/* Funcio que calcula la "propensity function" per a una reaccio de segon ordre	tipus I			*/
	double value;
	value = constant*(double)molec1*(double)molec2;
	return value;
}

int *rctn_order(int R, int *num_react, int **react) {
/* Funcio que retorna la matriu type[1..R] que conte l'ordre de cada reaccio. 		*/
/* Si type[i] = 0 => reaccio d'ordre zero; si type[i] = 1 => reaccio de primer ordre;	*/
/* si type[i] = 2 => reaccio de segon ordre tipus I.					*/

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
			printf("Reaccio d'ordre diferent a ordre zero, primer ordre i segon ordre tipus I.\n");
			exit(-1);
		}
	}
	return type;
#ifdef DEBUG
for (i=1; i<R+1; i++) printf("DEBUG reaction type: type[%i] = %i\n", i, type[i]);
#endif
}

double compute_h(int rctn, int *type, double *cnt_react, int **react, int *molecules, double Hill, int PHill0) {
/* Funcio que calcula el valor de la "propensity function" per a una reaccio donada		*/
/* Input:	rctn = num. identificador de la reaccio de la qual volem calcular h (h[i])	*/
/*		type[1..R] = vector que conte l'ordre de cada reaccio (ordre 0, 1 o 2 tipus I)	*/
/*		cnt_react[1..R] = vector que conte les "rate constants" de cada reaccio		*/
/*		react[1..R][0,1] = vector que conte les etiquetes dels reactants de cada reaccio*/
/*		molecules[1..MAXREACT] = vector que conte el num. de molecules existents de 	*/
/*					 cada element						*/
/*		Hill = coeficient de la funcio de Hill (per a la reaccio 3)			*/
/*		PHill0 = "threshold" de p53 (per a la reaccio 3)				*/

	double value, dummy;
	double (*hazard)();

	if (rctn == 3) {//Calcular la propensity function segons la funcio de Hill monotonicament creixent
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
/* Funcio que retorna el valor de la "propensity function" per a cadascuna de les reaccions	*/
/* Input:	R = num. de reaccions de la xarxa						*/
/*		num_react[1..R] = num. de reactants de cada reaccio				*/
/*		react[1..R][num_react[1..R]] = etiquetes dels reactants de cada reaccio	*/
/*		cnt_react[1..R] = "rate constants" de cada reaccio				*/
/*		molecules[react[1..R][j] = num. de molecules de cadascun dels reactants	*/
/*		Hill = coeficient de la funcio de Hill (per a la reaccio 3)			*/
/*		PHill0 = "threshold" de p53 (per a la reaccio 3)				*/
/* Output:	h[1..R] = vector que conte el valor de la "propensity function" per a cada reaccio	*/

/* Etiquetem cadascuna de les reaccions per identificar si son d'ordre zero, de	*/
/* primer ordre, de segon ordre tipus I o de segon ordre tipus II.			*/
	int i;
	double aux;
	double *h;
	h = calloc(R+1, sizeof(double));
	h[0] = 0.0;
/* Calculem la "propensity function" h[i] (on i es el num. de la reaccio) (la notacio de l'article	*/
/* de Gibson es a[i]).											*/
	for (i=1; i <= R; i++) {
		aux = compute_h(i, type, cnt_react, react, molecules, Hill, PHill0);
		h[i] = aux;
	}
	return h;
}

double minimum(int R, double *delta, int *indxr, int count, int *del_list, NODE **node, int *indxr2, double t, int *rctnexec, int *cmpltn) {
/* Subrutina que calcula el valor minim de l'array delta[1,..,R], el valor minim de node[]->tau_dat, 	*/
/* calcula despres el valor de (node[]->tau_dat)_min - t i compara aquest valor amb el de (delta[])_min */
/* retornant aixi el valor minim de tots dos, i retornant tambe (mitjançant l'apuntador rctnexec el 	*/
/* numero de la reaccio executada).									*/
/* INPUT:	R = numero de reaccions de la xarxa							*/
/*		delta[1,..,R] = vector que conte els increments de temps per a cadascuna de les reaccions*/
/*		indxr[1,..,R] = vector que conte la llista ordenada de les reaccions, de menor a major	*/
/*				delta									*/
/*		count = quantes reaccions hi ha amb retard 						*/
/*		del_list[1,..,count] = llistat de les diferents reaccions amb retard que hi ha a la xarxa*/
/*		node[1,..,R] = estructura "singly linked list" que conte els retards de les reaccions	*/
/*		indxr2[] = vector que conte la llista ordenada de les reaccions amb retards, de menor	*/
/*			   a major node (retard)							*/
/*		t = temps										*/
/*		rctnexec = apuntador que retornara el numero de la reaccio executada			*/
/*		cmpltn = apuntador que retorna 0 o 1 segons si la reaccio es una inicialitzacio o una	*/
/*			 finalitzacio									*/
/* OUTPUT:	value = min{delta[], (node[]->tau_dat)-t }						*/
	int i, j, flag;
	double value, val1, preval2, val2;
	indexx(R, delta, indxr);
	val1 = delta[indxr[1]]; //Valor minim del vector delta
	if (count > 0) {
		double *vec2;
		vec2 = calloc(count+1, sizeof(double));
		for (i = 1; i <= count; i++) {
			vec2[i] = node[del_list[i]]->tau_dat;
		}
		indexx(count, vec2, indxr2);
		preval2 = vec2[indxr2[1]]; //Valor minim dels retards
		val2 = preval2 - t;
		if (val1 < val2) {//Inicio una reaccio
			*rctnexec = indxr[1];
			*cmpltn = 0;
			return val1;
		}
		else if (val2 < val1) {//Completo una reaccio
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
			printf("val1 = %g, val2 = %g, quina reaccio triem doncs!?\n", val1, val2);
			exit(-1);
		}
	} else {//No tenim reaccions amb retard => no tenim estructura de vectors s[j]
		*rctnexec = indxr[1];//Inicio una reaccio
		*cmpltn = 0;
		return val1;

	}


}



















