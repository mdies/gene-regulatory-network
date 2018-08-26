#include "nrutil.h"
#define SWAP(a,b) itemp=(a); (a)=(b); (b)=itemp;
#define M 7;
#define NSTACK 50;

void indexx(int n, double *arr, int *indx) {
/* Indexes an array arr[1...n], i.e., outputs the array indx[1...n] such that arr[indx[j]] is	*/
/* in ascending order for j = 1, 2, ..., N. The input quantities n and arr are not changed.	*/
/* See Numerical Recipes in C (2nd. Ed.) p.338-339.						*/
	int nstack_max = (int) NSTACK;
	int m_max = (int) M;
	int i, indxt, ir=n, itemp, j, k, l=1;
	int *istack;
	int jstack=0;
	double a;
	
	istack=ivector(1,nstack_max);

	for (j=1; j<=n; j++) indx[j]=j;
	for(;;) {
		if (ir-l < m_max) {
			for (j=l+1; j<=ir; j++) {
				indxt=indx[j];
				a=arr[indxt];
				for (i=j-1; i>=l; i--) {
					if (arr[indx[i]] <= a) break;
					indx[i+1]=indx[i];
				}
				indx[i+1]=indxt;
			}
			if (jstack == 0) break;
			ir=istack[jstack--];
			l=istack[jstack--];
		} else {
			k=(l+ir) >> 1;
			SWAP(indx[k],indx[l+1]);
			if (arr[indx[l]] > arr[indx[ir]]) {
				SWAP(indx[l],indx[ir]);
			}
			if (arr[indx[l+1]] > arr[indx[ir]]) {
				SWAP(indx[l+1],indx[ir]);
			}
			if (arr[indx[l]] > arr[indx[l+1]]) {
				SWAP(indx[l],indx[l+1]);
			}
			i=l+1;
			j=ir;
			indxt=indx[l+1];
			a=arr[indxt];
			for(;;) {
				do i++; while (arr[indx[i]] < a);
				do j--; while (arr[indx[j]] > a);
				if (j < i) break;
				SWAP(indx[i],indx[j]);
			}
			indx[l+1]=indx[j];
			indx[j]=indxt;
			jstack += 2;
			if (jstack > nstack_max) nrerror("NSTACK too small in indexx.");
			if (ir-i+1 >= j-l) {
				istack[jstack]=ir;
				istack[jstack-1]=i;
				ir=j-1;
			} else {
				istack[jstack]=j-1;
				istack[jstack-1]=l;
				l=i;
			}
		}
	}
	free_ivector(istack,1,nstack_max);
}

/*void SWAP(int *a, int *b){*/
/* Swaps the tree nodes a and b in order to update the index structure appropiately.		*/
/*        int c=*a;
        *a=*b;
        *b=c;
}*/

