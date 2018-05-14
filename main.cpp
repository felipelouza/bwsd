/*
 * BWSD 
 *
 * Authors: Felipe A. Louza, Simon Gog 
 * contact: louza@usp.br
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

#ifndef OUTPUT 
	#define OUTPUT 0 
#endif

#ifndef DEBUG
	#define DEBUG 0 
#endif

#ifndef TIME
  #define TIME 1
#endif

#ifndef SAVE_SPACE 
  #define SAVE_SPACE 1
#endif

typedef map<uint32_t, uint32_t> tMII;
typedef map<uint32_t, double> tMID;

typedef vector<tMID> tVMID;

#if SAVE_SPACE
typedef vector<tMII> tVMII;
#else
typedef uint_t** tVMII;
#endif

/******************************************************************************/

int compute_all_bwsd_wt(unsigned char** R, uint_t k, uint_t n, char* c_file);//algorithm 1

int compute_all_bwsd_rmq(unsigned char **R, uint_t k, uint_t n, char* c_file);//Simon's algorithm 

int compute_all_bwsd(unsigned char** R, uint_t k, uint_t n, char* c_file);//brute force

/******************************************************************************/

int main(int argc, char** argv){

	time_t t_start=0;clock_t c_start=0;

	int CHECK=0, MODE=0;

	if(argc!=6){
		dies(__func__,NULL);
	}

	unsigned char **R;
	int_t n=0;
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

	switch(MODE){

		case 1: printf("## BWSD_WT ##\n"); 
			time_start(&t_start, &c_start);
			compute_all_bwsd_wt(R, k, n, c_file);
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
			compute_all_bwsd(R, k, n, c_file);
			printf("TOTAL:\n");
			fprintf(stderr,"%.6lf\n", time_stop(t_start, c_start));
			break;
		
		default: break;
	}

	// validate	
	if(CHECK==1){
	}


return 0;
}

/******************************************************************************/

unsigned char* cat(unsigned char* s1, unsigned char* s2, uint_t *n){

	(*n) = (strlen((char*)s1)+1)+(strlen((char*)s2)+1)+1; //add 0 at the end

//	printf("%s$%s$#\n", s1, s2);
//	cout<<"n = "<<*n<<endl;

	int_t j, l=0;
	unsigned char *str = (unsigned char*) malloc((*n)*sizeof(unsigned char));

	/**/
	{
		int_t m = strlen((char*)s1);
		for(j=0; j<m; j++) if(s1[j]<255) str[l++] = s1[j]+1;
		str[l++] = 1; //add 1 as separator
	}

	{
		int_t m = strlen((char*)s2);
		for(j=0; j<m; j++) if(s2[j]<255) str[l++] = s2[j]+1;
		str[l++] = 1; //add 1 as separator
	}
	/**/

	str[l++]=0;
	if(*n>l){
	  str = (unsigned char*) realloc(str, (l)*sizeof(unsigned char)); 
	}
	*n = l;

return str;
}

/******************************************************************************/

unsigned char* cat_all(unsigned char** R, int k, uint_t *n){

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

double compute_distance(tMII &t, int s){

	double dist = 0.0;
	
	for(tMII::iterator it=t.begin(); it!=t.end(); ++it){
		if(it->second){
			#if OUTPUT==1			//computes D_M
				dist++;
			#elif OUTPUT==2		//computes D_E
				double tmp = (double)it->second/(double)s;
				dist += tmp*log2(tmp);
			#endif
		}
	}

	#if OUTPUT==1     //computes D_M
		dist = 1.0/dist -1.0;
	#elif OUTPUT==2   //computes D_E
		if(dist) 	dist *= -1.0;
	#endif

return dist;
}
/******************************************************************************/

int compute_all_bwsd_wt(unsigned char** R, uint_t k, uint_t n, char* c_file){

	int_t i;

	//Concatenate strings
	/**/
	unsigned char *str = cat_all(R, k, &n);

	printf("K = %" PRId32 "\n", k);
	printf("N = %" PRIdN " bytes\n", n);
	printf("OUTPUT = %d\n", OUTPUT);
	printf("sizeof(int) = %zu bytes\n", sizeof(int_t));

	#if SAVE_SPACE
		cout<<"SAVE_SPACE"<<endl;
	#endif

	#if DEBUG
		printf("R:\n");
		for(i=0; i<k; i++)
			printf("%" PRIdN ") %s (%zu)\n", i, R[i], strlen((char*)R[i]));
	#endif

	//free memory
	for(i=0; i<k; i++)
		free(R[i]);
	free(R);

	/**/

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
	/**/
	if(!load_from_cache(da, "da", m_config)){

		int_t *SA = new int_t[n];
		int_t *DA = new int_t[n];
		
		for(i=0; i<n; i++) SA[i]=DA[i]=0;
		gsacak(str, (uint_t*)SA, NULL, DA, (uint_t)n); //construct SA+DA
	
		for(i=0;i<n;i++) da[i]=DA[i];
		store_to_cache(da, "da", m_config);

		delete[] SA;
		delete[] DA;
	}
	
	free(str);

	#if OUTPUT
		tVMID	result(k);
	#endif

	#if TIME
		printf("#1. DA:\n");
		fprintf(stderr,"%.6lf\n", time_stop(t_start, c_start)); 
	#endif

	//COMPUTE WT(DA):
	/**/
	//wt_int<> wt;
	wm_int<> wt;
	if(!load_from_cache(wt, "wt", m_config)){
		construct_im(wt, da);
		store_to_cache(wt, "wt", m_config);
	}

	#if TIME
		printf("#2. WT(DA):\n");
		fprintf(stderr,"%.6lf\n", time_stop(t_start, c_start)); 
	#endif

	#if DEBUG	
		cout<<"DA:"<<endl;
		for(int_t i=0; i<n-1; i++) cout << wt[i] << ", ";
		cout<<wt[n-1]<<"\t("<<n-1<<")"<<endl;
	#endif

	//avoid second wt.rank()
	int_t *tmp = new int_t[k+1];

	//avoid wt.rank()-queries when next > qe
	int_t **pos = new int_t*[k+1];
	
	for(int_t i=0; i<k; i++){
		tmp[i]=0;
		uint64_t len = wt.rank(wt.size(), i);   
		pos[i] = new int[len+1];
	}
	
	for(int_t i=1; i<n; i++){
		pos[da[i]][tmp[da[i]]]=i;
		tmp[da[i]]++;
	}
	
	for(int_t i=0; i<k; i++) pos[i][tmp[i]]=n;

	int skip=0;
	int total=0;

/*
	int_t *size = new int_t[k];
	for(int_t i=0; i<k; i++) size[i]=wt.rank(wt.size(), i);
*/

	int_t *s= new int_t[k];
	int_t *ell = new int_t[k];

	for(int_t i=0; i<k-1; i++){

		//int_t qs, qe=0;
		int_t qe=0;
		uint64_t len_i = wt.rank(wt.size(), i);
		//uint64_t len_i = size[i];

		//init
		#if SAVE_SPACE
			tVMII t(k);
		#else
/*
			uint_t** t = new uint_t*[k];
			for(int_t j=i+1; j<k; j++){

				uint64_t len_j = wt.rank(wt.size(), j);
				uint_t total = 2*(min(len_i, len_j)+1);

				cout<<"("<<i<<", "<<j<<")\n";
				cout<<"("<<len_i<<", "<<len_j<<")\n##\n";

				t[j] = new uint_t [total];
				for(int_t r=0; r<total; r++) 
					t[j][r] = 0;
			}
*/
		#endif

		#if DEBUG == 2
			cout<<"i = "<<i<<", "<<len_i<<endl;
		#endif

		for(int_t j=i+1; j<k; j++) s[j] = tmp[j] = ell[j]=0;
	
//		for(int_t j=i+1; j<k; j++) next[j]= pos[j][0]; //wt.select(1, i);

		//forloop S^i[1..n_i]
		for(uint64_t p=1; p<len_i+1; p++){

			#if DEBUG == 2
				//cout << "[" <<qs <<", "<<qe << "]   \t0^1"<<endl;
				cout << "[" <<qe <<", "<< wt.select(p, i) << "]   \t0^1"<<endl;
			#endif
		
			qe = wt.select(p, i);//TODO: replace qe
			//if(qe!=wt.select(p, i)) cout<<"ERROR"<<endl;

			for(int_t j=i+1; j<k; j++){

total++;

/*
if(tmp[j]==size[j]){//jump
ell[j]++;
continue;
}
*/

//if(next[j]!=pos[j][tmp[j]]) cout<<"OPA\n";
	
//if(next[j]>qe){//jump
if(pos[j][tmp[j]]>qe){//jump
ell[j]++;
skip++;
continue;
}
				int_t occ = wt.rank(qe,j);
				int_t kj = occ - tmp[j];

				//int_t kj = wt.rank(qe,j) - wt.rank(qs,j);//TODO: reduce 1 rank-query
				//if(tmp[j]!=wt.rank(qs,j))cout<<"ERROR"<<endl;
				tmp[j]=occ;


//if(occ<size[j])
//	next[j]=pos[j][occ];
////	next[j]=wt.select(occ+1, j);
//else
//	next[j]=n;

/*
if(occ<size[j])
else
	next[j]=n;
*/

				#if DEBUG == 2
					cout << "###"<< j << ":\t1^"<< kj << endl;
				#endif

				if(kj>0){
					t[j][kj]++; //1^kj
					t[j][ell[j]]++;  //0^lj
					ell[j]=1;

					s[j]+=2;
				}
/*
				else{
					ell[j]++;
				}
*/
			}
		}

		//last iteration
		{
			#if DEBUG == 2
				cout << "[" <<qs <<", "<<n<< "]"<<endl;
			#endif
			for(int_t j=i+1; j<k; j++){

				//int_t kj = wt.rank(n,j) - wt.rank(qs,j);
				int_t kj = wt.rank(n,j) - tmp[j];
				if(kj>0){
					t[j][kj]++; //1^kj
					s[j]++;
				}
				t[j][ell[j]]++;  //0^lj
				s[j]++;

				#if DEBUG == 2
					cout << "***"<< j << ":\t1^"<< kj << endl;
				#endif

	
				#if OUTPUT
					result[i][j] = compute_distance(t[j], s[j]);
					#if DEBUG
						cout<<"["<<i<<", "<<j<<"]\t\ts="<<s[j]<<"\n";
						cout<<"D = "<<result[i][j]<<endl;
					#endif
				#endif

			}
		}

		#if DEBUG
			cout<<"sum = "<<skip<<" / "<<total<<" = "<<(double)skip/(double)total<<endl;
			cout<<"\n####\t";

			//output (tmp)
			for(int_t j=i+1; j<k; j++){

				//debugging
				cout<<"("<<i<<", "<<j<<")\n";
				for(int_t p=0; p<n;p++){
					if(da[p]==i){
						int_t count=1;
						while(da[++p]!=j && p<n)	if(da[p]==i)count++;
						cout<<i<<"^"<<count<<" "; 
						if(p==n) break;
					}
					if(da[p]==j){
						int_t count=1;
						while(da[++p]!=i && p<n)	if(da[p]==j)count++;
						cout<<j<<"^"<<count<<" "; 
						if(p==n) break;
						p--;
					}
				}
				cout<<endl;

				#if SAVE_SPACE
					for(tMII::iterator it=t[j].begin(); it!=t[j].end(); ++it){
						if(it->second){
							cout << "t_" << it->first << ":\t" << it->second <<endl;
						}
					}
				#endif

				cout<<"####\t";
			}		
			cout<<endl;
		
		#endif
	}

	#if TIME
		printf("#3. ALL-BWSD:\n");
		fprintf(stderr,"%.6lf\n", time_stop(t_start, c_start)); 
	#endif

	delete[] ell;
	delete[] s;
	delete[] tmp;
//	delete[] next;

	for(int_t i=0; i<k; i++) delete[] pos[i];
	delete[] pos;


	#if DEBUG
		#if OUTPUT
			for(int_t i=0; i<k; i++){
				for(int_t j=0; j<i; j++){
					printf("%.2lf\t", result[j][i]);
				}
				printf("0.00\t");
				for(int_t j=i+1; j<k; j++){
					printf("%.2lf\t", result[i][j]);
				}
				cout<<endl;
			}
		#endif
	#endif


	//checksum: for the sake of sanity
	#if OUTPUT
		double sum=0.0;
		for(int_t i=0; i<k; i++)
			for(int_t j=i+1; j<k; j++)
				sum+=result[i][j];

		printf("checksum = %lf\n",sum);
	#endif


return 0;
}

/******************************************************************************/

int compute_all_bwsd(unsigned char** R, uint_t k, uint_t n, char* c_file){//brute force

	printf("K = %" PRId32 "\n", k);
	printf("N = %" PRIdN " bytes\n", n);
	printf("sizeof(int) = %zu bytes\n", sizeof(int_t));
	printf("OUTPUT = %d\n", OUTPUT);

	#if DEBUG
		printf("R:\n");
		for(int_t i=0; i<k; i++)
			printf("%" PRIdN ") %s (%zu)\n", i, R[i], strlen((char*)R[i]));
	#endif

	tVMID	result(k);

	/**/
	for(int_t i=0; i<k-1; i++){

		#if DEBUG 
			cout<<endl;
		#endif
		tVMII t(k);


		//foreach S^j in j+1..k
		for(int_t j=i+1; j<k; j++){
	
			int_t s= 0;
			//concatenates
			unsigned char *str = cat(R[i], R[j], &n);
			//printf("%s\t(%d)\n\n", str,n);

			//build DA
			
			int_t *SA = new int_t[n];
			int_t *DA = new int_t[n+1];
			for(int_t p=0; p<n; p++) SA[p]=DA[p]=0;
			gsacak(str, (uint_t*)SA, NULL, DA, (uint_t)n); //construct SA+DA

			#if DEBUG 
				cout<<"####\t";
				cout<<"("<<i<<", "<<j<<")\n";
			#endif

			free(str);

			//output
			for(int_t p=1; p<n;){

					{
						int_t count=1;
						while(DA[++p]==0 && p<n) count++;
						#if DEBUG 
							cout<<i<<"^"<<count<<" "; 
						#endif

						t[j][count]++; //0^count
						s++;
					}
					if(p==n) break;

					{
						int_t count=1;
						while(DA[++p]==1 && p<n) count++;
						#if DEBUG
							cout<<j<<"^"<<count<<" "; 
						#endif

						t[j][count]++; //1^count
						s++;
					}
			}
			#if DEBUG
				cout<<endl;
			#endif

			#if OUTPUT
				result[i][j] = compute_distance(t[j], s);
			#endif

			#if DEBUG
			for(tMII::iterator it=t[j].begin(); it!=t[j].end(); ++it)
				if(it->second)
					cout << "#^" << it->first << ":\t" << it->second <<endl;
			#endif

			delete[] DA;
			delete[] SA;
		}//for j=i+1..n

		#if DEBUG
			cout<<"####\n";
		#endif
	}

	//free memory
	for(int_t i=0; i<k; i++)
		free(R[i]);
	free(R);

	#if DEBUG
		#if OUTPUT
			for(int_t i=0; i<k; i++){
				for(int_t j=0; j<i; j++){
					printf("%.2lf\t", result[j][i]);
				}
				printf("0.00\t");
				for(int_t j=i+1; j<k; j++){
					printf("%.2lf\t", result[i][j]);
				}
				cout<<endl;
			}
		#endif
	#endif

	//checksum: for the sake of sanity
	#if OUTPUT
		double sum=0.0;
		for(int_t i=0; i<k; i++)
			for(int_t j=i+1; j<k; j++)
				sum+=result[i][j];

		printf("checksum = %lf\n",sum);
	#endif

return 0;
}

/******************************************************************************/

int compute_all_bwsd_rmq(unsigned char *s, uint_t k, uint_t n, char* c_file){//Simon's algorithm 



return 0;
}
/******************************************************************************/
