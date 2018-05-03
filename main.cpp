/*
 * BWSD 
 *
 * Authors: Felipe A. Louza, Simon Gog 
 * contact: louza@ic.unicamp.br
 * 02/05/2018
 *
 */

#include <sdsl/wavelet_trees.hpp>                                                                        
#include <sdsl/int_vector.hpp>                                                                           
                                                                                                         
#include <iostream>                                                                                      
                                                                                                         
using namespace sdsl;                                                                                    
using namespace std;  

#include <cstdio>
#include <cerrno>
#include <ctime>
#include <climits>

#include "lib/utils.h"
#include "lib/file.h"
#include "lib/bwt.h"
#include "external/malloc_count/malloc_count.h"
#include "external/gsacak.h"


#ifndef DEBUG
	#define DEBUG 0 
#endif

#ifndef TIME
  #define TIME 1
#endif

/******************************************************************************/

unsigned char* cat_char(unsigned char** R, int k, int_t *n);

int compute_all_bwsd_wt(unsigned char *s, uint_t k, uint_t n);//algorithm 1

int compute_all_bwsd_rmq(unsigned char *s, uint_t k, uint_t n);//Simon's algorithm 

/******************************************************************************/

int main(int argc, char** argv){

	time_t t_start=0;clock_t c_start=0;

	int CHECK=0, MODE=0;

	if(argc!=6){
		dies(__func__,NULL);
	}

	unsigned char **R;
	int_t i, n=0;
	int   k;

	char* c_dir = argv[1];
	char* c_file = argv[2];

	sscanf(argv[3], "%d", &k);
	sscanf(argv[4], "%u", &MODE);
	sscanf(argv[5], "%u", &CHECK);


	file_chdir(c_dir);

	//disk access
	R = (unsigned char**) file_load_multiple(c_file, k, &n);
	if(!R){
		fprintf(stderr, "Error: less than %d strings in %s\n", k, c_file);
		return 0;
	}

	//concatenate strings
	unsigned char *str = NULL;
	str = cat_char(R, k, &n);

	printf("K = %" PRId32 "\n", k);
	printf("N = %" PRIdN " bytes\n", n);
	printf("sizeof(int) = %zu bytes\n", sizeof(int_t));

	#if DEBUG
		printf("R:\n");
		for(i=0; i<k; i++)
			printf("%" PRIdN ") %s (%zu)\n", i, R[i], strlen((char*)R[i]));
	#endif

	//free memory
	for(i=0; i<k; i++)
		free(R[i]);
	free(R);

	switch(MODE){

		case 1: printf("## BWSD_WT ##\n"); 
			time_start(&t_start, &c_start);
			compute_all_bwsd_wt(str, k, n);
			printf("TOTAL:\n");
			fprintf(stderr,"%.6lf\n", time_stop(t_start, c_start));
			break;

		case 2:	printf("## BWSD_RMQ ##\n"); //Simon's algorithm 
			time_start(&t_start, &c_start);
//		compute_all_bwsd_rmq(str, k, n);
			printf("TOTAL:\n");
			fprintf(stderr,"%.6lf\n", time_stop(t_start, c_start));
			break;

		case 3:	printf("## BWSD ##\n");
			time_start(&t_start, &c_start);
			//brute force algorithm
			printf("TOTAL:\n");
			fprintf(stderr,"%.6lf\n", time_stop(t_start, c_start));
			break;
		
		default: break;
	}

	// validate	
	if(CHECK==1){
	}

	free(str);


return 0;
}

/******************************************************************************/

unsigned char* cat_char(unsigned char** R, int k, int_t *n){

	(*n)++; //add 0 at the end
	int_t i, j, l=0;
	unsigned char *str = (unsigned char*) malloc((*n)*sizeof(unsigned char));

	for(i=0; i<k; i++){
		int_t m = strlen((char*)R[i]);
		for(j=0; j<m; j++){
			if(R[i][j]<255) str[l++] = R[i][j]+1;
//			str[l++] = R[i][j];
		}
		str[l++] = 1; //add 1 as separator
	}

	str[l++]=0;
	if(*n>l){
	  str = (unsigned char*) realloc(str, (l)*sizeof(unsigned char)); 
	}
	*n = l;

return str;
}
/******************************************************************************/

int compute_all_bwsd_wt(unsigned char *s, uint_t k, uint_t n){

	int_t i;

	#if TIME
    time_t t_start=0;clock_t c_start=0;
		time_start(&t_start, &c_start); 
	#endif

  uint_t *SA = (uint_t*) malloc(n*sizeof(uint_t));
  int_t *DA = (int_t*) malloc(n*sizeof(int_t));

  for(i=0; i<n; i++) SA[i]=DA[i]=0;

  gsacak(s, SA, NULL, DA, (uint_t)n); //construct SA+DA

	#if TIME
		printf("1. BUILD SA+DA:\n");
		fprintf(stderr,"%.6lf\n", time_stop(t_start, c_start)); 
	#endif

  #if DEBUG
    s[n-1]='#'+1;
    printf("i) DA\tSA\tBWT\tsuffixes\n");
    for(i=0; i<min(n,20); i++){
      printf("%" PRIdN ") %" PRIdN "\t%" PRIdN "\t", i, DA[i], SA[i]);
      char bwt = SA[i]>0?s[SA[i]-1]:'$';
      if(bwt==1)bwt='$'+1;
      printf("%c\t", bwt-1);
      for(int_t j=SA[i]; j< min(n,SA[i]+10); j++)
        if(s[j]==1) printf("$");
        else printf("%c", s[j]-1);
      printf("\n");
    }
  #endif


	int_vector<> da(n);
	for(i=0;i<n;i++) da[i]=DA[i];
	wt_int<> wt;

	construct_im(wt, da);
	
	#if TIME
		printf("2. BUILD WT(DA):\n");
		fprintf(stderr,"%.6lf\n", time_stop(t_start, c_start)); 
	#endif



	cout << "wt.size()="<< wt.size() << endl;
	cout << "wt.sigma ="<< wt.sigma << endl;
	if (wt.size() > 0) {
	    // access an element
	    cout << "wt[0]=" << wt[0] << endl;
	    // rank an element (exclude)
	    uint64_t r = wt.rank(wt.size(), wt[0]);
	    cout << "wt.rank(wt.size(), wt[0])=" << r  << endl;
	    // select element ()
	    cout << "wt.select(r, wt[0]) = " << wt.select(r, wt[0]) << endl;
	}	




  free(SA);
  free(DA);

return 0;
}
/******************************************************************************/
