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
	#define DEBUG 1 
#endif

#ifndef TIME
  #define TIME 1
#endif

typedef map<uint32_t, uint32_t> tMII;
typedef vector<tMII> tVMII;

/******************************************************************************/

unsigned char* cat_char(unsigned char** R, int k, int_t *n);

int compute_all_bwsd_wt(unsigned char *s, uint_t k, uint_t n, char* c_file);//algorithm 1

int compute_all_bwsd_rmq(unsigned char *s, uint_t k, uint_t n, char* c_file);//Simon's algorithm 

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
			compute_all_bwsd_wt(str, k, n, c_file);
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

int compute_all_bwsd_wt(unsigned char *s, uint_t k, uint_t n, char* c_file){

	int_t i;

	string dir = "sdsl";
	mkdir(dir.c_str());
	string id = c_file;
	id += "."+to_string(k);

	cache_config m_config(true, dir, id);
	int_vector<> da(n);

	#if TIME
	  time_t t_start=0;clock_t c_start=0;
		time_start(&t_start, &c_start); 
	#endif

	//COMPUTE DA:
	if(!load_from_cache(da, "da", m_config)){

		int_t *SA = new int_t[n];
		int_t *DA = new int_t[n];
		
		for(i=0; i<n; i++) SA[i]=DA[i]=0;
		gsacak(s, (uint_t*)SA, NULL, DA, (uint_t)n); //construct SA+DA
	
		for(i=0;i<n;i++) da[i]=DA[i];
		store_to_cache(da, "da", m_config);

		delete[] SA;
		delete[] DA;
	}

	#if TIME
		printf("#1. DA:\n");
		fprintf(stderr,"%.6lf\n", time_stop(t_start, c_start)); 
	#endif

	wt_int<> wt;
	//COMPUTE WT(DA):
	if(!load_from_cache(wt, "wt", m_config)){
		construct_im(wt, da);
		store_to_cache(wt, "wt", m_config);
	}

	#if TIME
		printf("#2. WT(DA):\n");
		fprintf(stderr,"%.6lf\n", time_stop(t_start, c_start)); 
	#endif

	/**/
	cout<<"DA:"<<endl
	for(int_t i=0; i<n-1; i++) cout << wt[i] << ", ";
	cout<<wt[n-1]<<"\t("<<n-1<<")"<<endl;
	/**/

	//for(int_t i=0; i<k-1; i++){
	for(int_t i=0; i<1; i++){

		int_t qs = 0, qe=0;
		uint64_t len = wt.rank(wt.size(), i);
	  
		tVMII F(k);

		#if DEBUG
			cout<<"i = "<<i<<", "<<len<<endl;
		#endif

		int_t *ell = new int_t[n];
		for(int_t j=i+1; j<k; j++) ell[j]=0;

		//forloop S^i[1..n_i]
		for(uint64_t p=1; p<len+1; p++){
		
			qe = wt.select(p, i);

			#if DEBUG
				cout << "[" <<qs <<", "<<qe << "]   \t0^1"<<endl;
			#endif


			//foreach S^j in j+1..k
			for(int_t j=i+1; j<k; j++){

				int_t kj = wt.rank(qe,j) - wt.rank(qs,j);
				#if DEBUG
					cout << "###"<< j << ":\t1^"<< kj << endl;
				#endif

				if(kj>0){
					F[j][kj]++; //1^kj
					F[j][1]++;  //0^1
					ell[j]=1;
				}
				else{
			
					if(p==1) F[j][1]++;
					else{
						F[j][ell[j]]--; //1^kj
						F[j][ell[j]+1]++;  //0^1
					}
					ell[j]++;
				}
			}

			qs=qe+1;
		}

		//last iteration
		{
			#if DEBUG
				cout << "[" <<qs <<", "<<n<< "]"<<endl;
			#endif
			for(int_t j=i+1; j<k; j++){

				int_t kj = wt.rank(n,j) - wt.rank(qs,j);
				if(kj>0) F[j][kj]++; //1^kj

				#if DEBUG
					cout << "***"<< j << ":\t1^"<< kj << endl;
				#endif
			}
		}
		
		cout<<"\n####\n";
		//output (tmp)
		for(int_t j=i+1; j<k; j++){

			for(int_t p=0; p<n;p++) 
				if(da[p]==i || da[p]==j) cout<<da[p]<<"^1 "; 

			cout<<"("<<i<<", "<<j<<")\n";

			for(tMII::iterator it=F[j].begin(); it!=F[j].end(); ++it)
				cout << "#^" << it->first << ":\t" << it->second <<endl;

			cout<<"####\n";
		}		
	}


return 0;
}
/******************************************************************************/

int compute_all_bwsd_rmq(unsigned char *s, uint_t k, uint_t n, char* c_file);{//Simon's algorithm 



return 0;
}

/******************************************************************************/
