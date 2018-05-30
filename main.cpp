/*
 * BWSD 
 *
 * Authors: Felipe A. Louza, Simon Gog 
 * contact: louza@usp.br
 * 22/05/2018
 *
 */

#include <sdsl/wavelet_trees.hpp>                                                                        
#include <sdsl/int_vector.hpp>                                                                           
#include <sdsl/bit_vectors.hpp>                                                                           
#include <sdsl/rank_support.hpp>
#include <sdsl/rmq_support.hpp>

#include <iostream>

using namespace sdsl;
using namespace std;  

#include <cstdio>
#include <cerrno>
#include <ctime>
#include <climits>
#include <map>
#include <string>

#include "lib/utils.h"
#include "lib/file.h"
#include "lib/bwt.h"
#include "external/gsacak.h"

#ifndef OUTPUT 
	#define OUTPUT 1 
#endif

#ifndef DEBUG
	#define DEBUG 0 
#endif

#ifndef TIME
  #define TIME 1
#endif

#ifndef WT
  #define WT 0
#endif

#ifndef OPT_VERSION
  #define OPT_VERSION 1
#endif

#ifndef SD_VECTOR 
  #define SD_VECTOR 1
#endif

//debug purposes
#ifndef WORST_CASE 
  #define WORST_CASE 1
#endif


//#define Result(i,j) (result[i][j])
#define Result(i,j) (Md[(k*(k-1)/2) - (k-i)*((k-i)-1)/2 + j - i - 1])

typedef map<uint32_t, uint32_t> tMII;
typedef map<uint32_t, double> tMID;

typedef vector<tMID> tVMID;
typedef vector<tMII> tVMII;

typedef unordered_map<size_t,size_t> tUMII;
typedef array<size_t,2> tAII; 


/******************************************************************************/

int compute_all_bwsd_rank(unsigned char** R, uint_t k, uint_t n, char* c_file);//algorithm 1

int compute_all_bwsd_rmq(unsigned char **R, uint_t k, uint_t n, char* c_file);//Simon's algorithm 

int compute_all_bwsd_nk(unsigned char **R, uint_t k, uint_t n, char* c_file);//Simon's algorithm 

int compute_all_bwsd(unsigned char** R, uint_t k, uint_t n, char* c_file);//brute force

/******************************************************************************/

int main(int argc, char** argv){

	time_t t_start=0;clock_t c_start=0;

	int MODE=0;

	if(argc!=5){
		dies(__func__,NULL);
	}

	unsigned char **R;
	int_t n=0;
	int   k;

	char* c_dir = argv[1];
	char* c_file = argv[2];

	sscanf(argv[3], "%d", &k);
	sscanf(argv[4], "%u", &MODE);

	file_chdir(c_dir);

	//disk access
	R = (unsigned char**) file_load_multiple(c_file, k, &n);
	if(!R){
		fprintf(stderr, "Error: less than %d strings in %s\n", k, c_file);
		return 0;
	}

	switch(MODE){

		case 1: printf("## BWSD_RANK ##\n"); 
			time_start(&t_start, &c_start);
			compute_all_bwsd_rank(R, k, n, c_file);
			printf("TOTAL:\n");
			fprintf(stderr,"%.6lf\n", time_stop(t_start, c_start));
			break;

		case 2:	printf("## BWSD_RMQ ##\n"); //Simon's algorithm 
			time_start(&t_start, &c_start);
      compute_all_bwsd_rmq(R, k, n, c_file);
			printf("TOTAL:\n");
			fprintf(stderr,"%.6lf\n", time_stop(t_start, c_start));
			break;

		case 3:	printf("## BWSD ##\n");
			time_start(&t_start, &c_start);
			//straightforward algorithm
			compute_all_bwsd(R, k, n, c_file);
			printf("TOTAL:\n");
			fprintf(stderr,"%.6lf\n", time_stop(t_start, c_start));
			break;
		
		default: break;
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

	int_t max=0;

	for(i=0; i<k; i++){
		int_t m = strlen((char*)R[i]);
		if(m>max) max=m;
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

	cout<<"max = "<<max<<endl;

return str;
}
/******************************************************************************/

double compute_distance(tMII &t, int_t s){

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

void push_if_not_empty(stack<tAII>& s, tAII range){
    if ( range[0] < range[1]+1 ) {
        s.push(range);
    }
}

/******************************************************************************/


int compute_all_bwsd_rank(unsigned char** R, uint_t k, uint_t n, char* c_file){

	int_t i;

	//Concatenate strings
	/**/
	unsigned char *str = cat_all(R, k, &n);

	printf("K = %" PRId32 "\n", k);
	printf("N = %" PRIdN " bytes\n", n);
	printf("OUTPUT = %d\n", OUTPUT);
	printf("sizeof(int) = %zu bytes\n", sizeof(int_t));

	#if OPT_VERSION 
		cout<<"OPT_VERSION"<<endl;
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

	#if WORST_CASE
		int_t j1, j2=0;
		for(int_t i=0; i<k; i++){
			j1 = i*(n/k)+1;
			j2 = (i+1)*(n/k)+1;
			for(int_t j=j1; j<j2 && j<n; j++)
				da[j]=i;
		}
		for(int_t j=j2; j<n; j++)
				da[j]=k;
		//cout<<"DA: "; for(int_t i=0; i<n; i++) cout<<da[i]<<" ";cout<<endl;
	#endif

	#if OUTPUT
		//tVMID	result(k);
		size_t m = (k*k-k)/2.0;
		double *Md = new double[m];
	#endif
 
	#if TIME
		printf("#1. DA:\n");
		fprintf(stderr,"%.6lf\n", time_stop(t_start, c_start)); 
	#endif

	#if WT
		//COMPUTE WT(DA):
		wm_int<> wt;
		if(!load_from_cache(wt, "wt", m_config)){
			construct_im(wt, da);
			store_to_cache(wt, "wt", m_config);
		}
	#else

		#if SD_VECTOR
			sd_vector_builder *B = new sd_vector_builder[k];
		
			int_t *size = new int_t[k];
			for(int_t i=0; i<k; i++) size[i]=0;
			for(int_t i=1; i<n; i++) size[da[i]]++;

			for(int_t i=0; i<k; i++)
				B[i] = sd_vector_builder(n,size[i]);

			delete[] size;
		#else
			bit_vector *B = new bit_vector[k];
			for(int_t i=0; i<k; i++)
				B[i] = bit_vector(n,0);
		#endif

		for(int_t i=1; i<n; i++)
			#if SD_VECTOR
				B[da[i]].set(i);
			#else
				B[da[i]][i]=1;
			#endif

		#if SD_VECTOR
			sd_vector<> *B_sd = new sd_vector<>[k];
			for(int_t i=0; i<k; i++)
				B_sd[i] = sd_vector<>(B[i]);
			delete[] B;
		#endif

		//rank-support
		#if SD_VECTOR
			rank_support_sd<1> *B_rank = new rank_support_sd<1>[k];
		#else
			rank_support_v<1> *B_rank = new rank_support_v<1>[k];
		#endif

		for(int_t i=0; i<k; i++){
			#if SD_VECTOR
				B_rank[i] = rank_support_sd<1>(&B_sd[i]);
			#else
				B_rank[i].set_vector(&B[i]);
				util::init_support(B_rank[i], &B[i]);
			#endif
		}

		//select-support
		#if OPT_VERSION == 0
			#if SD_VECTOR
				select_support_sd<1> *B_select= new select_support_sd<1>[k];	
			#else
				select_support_mcl<1,1> *B_select= new select_support_mcl<1,1>[k];
			#endif
			for(int_t i=0; i<k; i++){
				#if SD_VECTOR
					B_select[i].set_vector(&B_sd[i]);	
				#else
					B_select[i].set_vector(&B[i]); 
					util::init_support(B_select[i], &B[i]);
				#endif		
			}
		#endif

	#endif //#if WT==0

	#if TIME
		#if WT
			printf("#2. WT(DA):\n");
		#else
			#if SD_VECTOR
				printf("#2. BITMAP_SD(DA):\n");
			#else
				printf("#2. BITMAP(DA):\n");
			#endif
		#endif
		fprintf(stderr,"%.6lf\n", time_stop(t_start, c_start)); 
	#endif

	#if DEBUG	
		cout<<"DA:"<<endl;
		for(int_t i=0; i<n-1; i++) cout << da[i] << ", ";
		cout<<da[n-1]<<"\t("<<n-1<<")"<<endl;
	#endif


	//avoid second wt.rank()
	int_t *rank = new int_t[k+1];

	#if OPT_VERSION 
		//avoid wt.rank()-queries when next > qe
		int_t **pos = new int_t*[k+1];
		
		for(int_t i=0; i<k; i++){
			rank[i]=0;
			#if WT
				uint64_t len = wt.rank(wt.size(), i);   
			#else
				uint64_t len = B_rank[i](n);// wt.rank(wt.size(), i);   
			#endif
			pos[i] = new int[len+1];
		}
		
		for(int_t i=1; i<n; i++){
			pos[da[i]][rank[da[i]]]=i;
			rank[da[i]]++;
		}
		
		for(int_t i=0; i<k; i++) pos[i][rank[i]]=n+1;

		uint64_t skip=0, total=0;
	#endif

	int_t *s= new int_t[k];
	int_t *ell = new int_t[k];

	for(int_t i=0; i<k-1; i++){

		int_t qe=0;
		#if WT
			uint64_t len_i = wt.rank(wt.size(), i);
		#else
			uint64_t len_i = B_rank[i](n);
		#endif

		//init
		tVMII t(k);

		#if DEBUG == 2
			cout<<"i = "<<i<<", "<<len_i<<endl;
		#endif

		for(int_t j=i+1; j<k; j++) s[j] = rank[j] = ell[j]=0;
	
		//forloop S^i[1..n_i]
		for(uint64_t p=1; p<len_i+1; p++){

			#if DEBUG == 2
				cout << "[" <<qe <<", "<< wt.select(p, i) << "]   \t0^1"<<endl;
			#endif
		
			#if OPT_VERSION
				qe = pos[i][p-1];
			#else
				#if WT
					qe = wt.select(p, i);
				#else
					qe = B_select[i](p);
				#endif
			#endif

			for(int_t j=i+1; j<k; j++){

				#if OPT_VERSION
					total++;
					if(pos[j][rank[j]]>qe){//skip if next_j(rank[j]+1) > qe
						ell[j]++;
						skip++;
						continue;
					}
				#endif
	
				//int_t kj = wt.rank(qe,j) - wt.rank(qs,j);//reduce 1 rank-query
				#if WT
					int_t occ = wt.rank(qe,j);
				#else
					int_t occ = B_rank[j](qe);
				#endif

				int_t kj = occ - rank[j];
				rank[j]=occ;

				#if DEBUG == 2
					cout << "###"<< j << ":\t1^"<< kj << endl;
				#endif

				#if OPT_VERSION == 0
					if(kj>0){
				#endif
						t[j][kj]++; //1^kj
						t[j][ell[j]]++;  //0^lj
						ell[j]=1;

						s[j]+=2;
				#if OPT_VERSION == 0
					}
					else{
						ell[j]++;
					}
				#endif
			}
		}

		//last iteration
		{
			#if DEBUG == 2
				cout << "[" <<qs <<", "<<n<< "]"<<endl;
			#endif
			for(int_t j=i+1; j<k; j++){

				#if OPT_VERSION
					total++;

					if(pos[j][rank[j]]>n){//jump if next_j(rank[j]+1) > qe
						t[j][ell[j]]++;  //0^lj
						s[j]++;
						skip++;
					}
					else{
				#endif

				//int_t kj = wt.rank(n,j) - wt.rank(qs,j);
				#if WT
					int_t kj = wt.rank(n,j) - rank[j];
				#else
					int_t kj = B_rank[j](n) - rank[j];
				#endif

				#if OPT_VERSION == 0
					if(kj>0){
				#endif
						t[j][kj]++; //1^kj
						s[j]++;
				#if OPT_VERSION == 0
					}
				#endif
				t[j][ell[j]]++;  //0^lj
				s[j]++;

				#if DEBUG == 2
					cout << "***"<< j << ":\t1^"<< kj << endl;
				#endif
	
				#if OPT_VERSION
					}
				#endif

				#if OUTPUT
					Result(i,j) = compute_distance(t[j], s[j]);
					#if DEBUG
						cout<<"["<<i<<", "<<j<<"]\t\ts="<<s[j]<<"\n";
						cout<<"D = "<<Result(i,j)<<endl;
					#endif
				#endif

			}
		}

		#if DEBUG
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

				/*
				for(tMII::iterator it=t[j].begin(); it!=t[j].end(); ++it){
					if(it->second){
						cout << "t_" << it->first << ":\t" << it->second <<endl;
					}
				}
				*/

				cout<<"####\t";
			}		
			cout<<endl;
		
		#endif
	}

	#if TIME
		printf("#3. BWSD-RANK:\n");
		fprintf(stderr,"%.6lf\n", time_stop(t_start, c_start)); 
	#endif

	delete[] ell;
	delete[] s;
	delete[] rank;

	#if OPT_VERSION
		for(int_t i=0; i<k; i++) delete[] pos[i];
		delete[] pos;
	#endif

	#if WT==0
		delete[] B_rank;
		#if OPT_VERSION == 0
			delete[] B_select;
		#endif
		#if SD_VECTOR
			delete[] B_sd;
		#else
			delete[] B;
		#endif
	#endif

	#if DEBUG
		#if OUTPUT
			for(int_t i=0; i<k; i++){
				for(int_t j=0; j<i; j++){
					printf("%.2lf\t", Result(j,i));
				}
				printf("0.00\t");
				for(int_t j=i+1; j<k; j++){
					printf("%.2lf\t", Result(i,j));
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
				sum+=Result(i,j);

		printf("checksum = %lf\n",sum);
	#endif

	#if OUTPUT
		delete[] Md;
	#endif

	#if OPT_VERSION
		cout<<"skip = "<<skip<<" / "<<total<<" = "<<(double)skip/(double)total<<endl;
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

	//tVMID	result(k);
	size_t m = (k*k-k)/2.0;
	double *Md = new double[m];

	/**/
	for(int_t i=0; i<k-1; i++){

		#if DEBUG 
			cout<<endl;
		#endif


		//foreach S^j in j+1..k
		for(int_t j=i+1; j<k; j++){
	
			tMII t;

			int_t s= 0;
			//concatenates
			unsigned char *str = cat(R[i], R[j], &n);
			//printf("%s\t(%d)\n\n", str,n);

			//build DA
			
			int_t *SA = new int_t[n];
			int_t *DA = new int_t[n+1];
			for(int_t p=0; p<n; p++) SA[p]=DA[p]=0;
			gsacak(str, (uint_t*)SA, NULL, DA, (uint_t)n); //construct SA+DA

			#if WORST_CASE
				int_t j1=1;
				int_t j2=n/2+1;
				for(int_t i=j1; i<j2; i++) DA[i]=0;
				for(int_t i=j2; i<n;	i++) DA[i]=1;
				//cout<<"DA: "; for(int_t i=0; i<n; i++) cout<<DA[i]<<" ";cout<<endl;
			#endif

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

						t[count]++; //0^count
						s++;
					}
					if(p==n) break;

					{
						int_t count=1;
						while(DA[++p]==1 && p<n) count++;
						#if DEBUG
							cout<<j<<"^"<<count<<" "; 
						#endif

						t[count]++; //1^count
						s++;
					}
			}
			#if DEBUG
				cout<<endl;
			#endif

			#if OUTPUT
				Result(i,j) = compute_distance(t, s);
			#endif

			//#if DEBUG
			//for(tMII::iterator it=t[j].begin(); it!=t[j].end(); ++it)
			//	if(it->second)
			//		cout << "#^" << it->first << ":\t" << it->second <<endl;
			//#endif

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
					printf("%.2lf\t", Result(j,i));
				}
				printf("0.00\t");
				for(int_t j=i+1; j<k; j++){
					printf("%.2lf\t", Result(i,j));
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
				sum+=Result(i,j);

		printf("checksum = %lf\n",sum);
	#endif

	#if OUTPUT
		delete[] Md;
	#endif

return 0;
}

/******************************************************************************/

int compute_all_bwsd_rmq(unsigned char** S, uint_t k, uint_t n, char* c_file){//Simon's algorithm 

	int_t i;
	
	//Concatenate strings
	/**/
	unsigned char *str = cat_all(S, k, &n);
	
	printf("K = %" PRId32 "\n", k);
	printf("N = %" PRIdN " bytes\n", n);
	printf("OUTPUT = %d\n", OUTPUT);
	printf("sizeof(int) = %zu bytes\n", sizeof(int_t));
	
	#if DEBUG
		printf("R:\n");
		for(i=0; i<k; i++)
			printf("%" PRIdN ") %s (%zu)\n", i, S[i], strlen((char*)S[i]));
	#endif
	
	//free memory
	for(i=0; i<k; i++)
		free(S[i]);
	free(S);
	
	/**/
	
	string dir = "sdsl";
	mkdir(dir.c_str());
	string id = c_file;
	id += "."+to_string(k);
	
	cache_config m_config(true, dir, id);
	int_vector<> da(n);
	//vector<uint64_t> da(n);
	
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

	#if WORST_CASE
		int_t j1, j2=0;
		for(int_t i=0; i<k; i++){
			j1 = i*(n/k)+1;
			j2 = (i+1)*(n/k)+1;
			for(int_t j=j1; j<j2 && j<n; j++)
				da[j]=i;
		}
		for(int_t j=j2; j<n; j++)
				da[j]=k;
		//cout<<"DA: "; for(int_t i=0; i<n; i++) cout<<da[i]<<" ";cout<<endl;
	#endif

	#if OUTPUT
		//tVMID	result(k);
		size_t m = (k*k-k)/2.0;
		double *Md = new double[m];
	#endif
	
	#if TIME
		printf("#1. DA:\n");
		fprintf(stderr,"%.6lf\n", time_stop(t_start, c_start)); 
	#endif
	
	auto max_da = *std::max_element(da.begin(), da.end());
	
	#if DEBUG  
		cout << max_da << endl;
		auto print_array = [](int_vector<>& vec, string label) {
			cout << label << ":";
			for(const auto& x : vec) {
				cout << " " << setw(2) << x;
			}
		  cout << endl;
		}; 
		print_array(da, "da  ");
	#endif  
	
	int_vector<> P(da.size(), da.size());
	int_vector<> N(da.size(), 0);
	int_vector<> R(da.size(), 0);
	
	{
		int_vector<> last_occ(max_da+1, 0);
		for (size_t i=0; i < da.size(); ++i) {
			P[i] = last_occ[da[i]];
			last_occ[da[i]] = i+1;
		}
	}
	#if DEBUG
	  print_array(P, "P+1");
	#endif
	{
		vector<uint64_t> next_occ(max_da+1, da.size());
		for(size_t i=da.size(); i > 0; --i) {
			N[i-1] = next_occ[da[i-1]];
			next_occ[da[i-1]] = i-1;
		}
	}
	#if DEBUG
		print_array(N, "N  ");
	#endif
	{
		for(size_t i=0; i<da.size(); ++i){
			if ( P[i] > 0 ) {
				R[i] = R[P[i]-1]+1;
			}
		}
	}
	#if DEBUG
	  print_array(R, "R  ");
	#endif
	
	rmq_succinct_sct<true>  rmq_P(&P);
	rmq_succinct_sct<false> RMQ_N(&N);
	
	vector<size_t> seen(max_da+1,0); 
	stack<size_t>  seen_stack;
	vector<size_t> freq(max_da+1,0); 
	vector<size_t> last_occ(max_da+1, 0);
	
	vector<vector<tMII>> counts(max_da+1, vector<tMII>(max_da+1));
	vector<vector<int_t>> runs(max_da+1, vector<int_t>(max_da+1));  
	
	#if TIME
	  printf("#2. RMQ:\n");
		fprintf(stderr,"%.6lf\n", time_stop(t_start, c_start)); 
	#endif
	
	for(size_t i=1; i < da.size() + max_da; ++i){
		size_t d  = i < da.size() ? da[i] : i-da.size();
		size_t lb = last_occ[d];
		size_t rb = i < da.size() ? i-1 : da.size()-1;
		last_occ[d] = i+1;
		
		#if DEBUG
		  cout << "i="<<setw(2)<<i<<" d="<<setw(2)<<d<<" ["<<lb<<","<<rb<<"] ";
		#endif
		
		stack<tAII> ranges;
		push_if_not_empty(ranges, {lb,rb});
		while ( !ranges.empty() ) {
			size_t _lb = ranges.top()[0];
			size_t _rb = ranges.top()[1];
			ranges.pop();
			size_t _m = rmq_P(_lb, _rb);
			if ( P[_m] < lb+1 ) { // equiv P[_m]-1 < lb 
				seen_stack.push(da[_m]);
				freq[da[_m]] = R[_m];
				push_if_not_empty(ranges, {_lb, _m-1});
				push_if_not_empty(ranges, {_m+1, _rb});
			}
		}
		#if DEBUG
		  cout << " unique docs: " << seen_stack.size() << " ";
		#endif
		
		push_if_not_empty(ranges, {lb,rb});
		while ( !ranges.empty() ) {
			size_t _lb = ranges.top()[0];
			size_t _rb = ranges.top()[1];
			ranges.pop();
			size_t _m = RMQ_N(_lb, _rb);
			if ( N[_m] > rb ) {
				freq[da[_m]] = R[_m] - freq[da[_m]] + 1;
				push_if_not_empty(ranges, {_m+1, _rb});
				push_if_not_empty(ranges, {_lb, _m-1});
			}
		}
		
		while( !seen_stack.empty() ) {
			auto x = seen_stack.top();
			seen_stack.pop();
			#if DEBUG
			  cout << " (d="<<x<<", f="<<freq[x]<<")";
			#endif
			auto i1 = std::min(d,x);
			auto i2 = std::max(d,x);
			++counts[i1][i2][freq[x]];
			++runs[i1][i2];
		}
		#if DEBUG
		  cout << endl;
		#endif
	}
	
	for(size_t i=0; i<=max_da; ++i){
		for(size_t j=i+1; j<=max_da; ++j) {
			Result(i,j) = compute_distance(counts[i][j], runs[i][j]);
		}
	}
	
	#if TIME
		printf("#3. BWSD-RMQ:\n");
		fprintf(stderr,"%.6lf\n", time_stop(t_start, c_start)); 
	#endif
	
	#if DEBUG
		#if OUTPUT
			for(int_t i=0; i<k; i++){
				for(int_t j=0; j<i; j++){
					printf("%.2lf\t", Result(j,i));
				}
				printf("0.00\t");
				for(int_t j=i+1; j<k; j++){
					printf("%.2lf\t", Result(i,j));
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
				sum+=Result(i,j);
	
		printf("checksum = %lf\n",sum);
	#endif
    
	#if OUTPUT
		delete[] Md;
	#endif

return 0;
}
/******************************************************************************/

int compute_all_bwsd_nk(unsigned char** S, uint_t k, uint_t n, char* c_file){//Simon's algorithm 

	int_t i;
	
	//Concatenate strings
	/**/
	unsigned char *str = cat_all(S, k, &n);
	
	printf("K = %" PRId32 "\n", k);
	printf("N = %" PRIdN " bytes\n", n);
	printf("OUTPUT = %d\n", OUTPUT);
	printf("sizeof(int) = %zu bytes\n", sizeof(int_t));
	
	#if DEBUG
		printf("R:\n");
		for(i=0; i<k; i++)
			printf("%" PRIdN ") %s (%zu)\n", i, S[i], strlen((char*)S[i]));
	#endif
	
	//free memory
	for(i=0; i<k; i++)
		free(S[i]);
	free(S);
	
	/**/
	
	string dir = "sdsl";
	mkdir(dir.c_str());
	string id = c_file;
	id += "."+to_string(k);
	
	cache_config m_config(true, dir, id);
	int_vector<> da(n);
	//vector<uint64_t> da(n);
	
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
		//tVMID	result(k);
		size_t m = (k*k-k)/2.0;
		double *Md = new double[m];
	#endif
	
	#if TIME
		printf("#1. DA:\n");
		fprintf(stderr,"%.6lf\n", time_stop(t_start, c_start)); 
	#endif
	
	auto max_da = *std::max_element(da.begin(), da.end());
	
	#if DEBUG  
		cout << max_da << endl;
		auto print_array = [](int_vector<>& vec, string label) {
			cout << label << ":";
			for(const auto& x : vec) {
				cout << " " << setw(2) << x;
			}
			cout << endl;
		}; 
		print_array(da, "da  ");
	#endif  
	
	int_vector<> P(da.size(), da.size());
	int_vector<> R(da.size(), 0);
	
	{
		int_vector<> last_occ(max_da+1, 0);
		for (size_t i=0; i < da.size(); ++i) {
			P[i] = last_occ[da[i]];
			last_occ[da[i]] = i+1;
		}
	}
	#if DEBUG
	  print_array(P, "P+1");
	#endif
	{
		for(size_t i=0; i<da.size(); ++i){
			if ( P[i] > 0 ) {
				R[i] = R[P[i]-1]+1;
			}
		}
	}
	#if DEBUG
		print_array(R, "R  ");
	#endif
	
	vector<size_t> seen(max_da+1,0); 
	stack<size_t>  seen_stack;
	vector<size_t> freq(max_da+1,0); 
	vector<size_t> last_occ(max_da+1, 0);
	
	vector<vector<tMII>> counts(max_da+1, vector<tMII>(max_da+1));
	vector<vector<int_t>> runs(max_da+1, vector<int_t>(max_da+1));  
	
	#if TIME
		printf("#2. RMQ:\n");
		fprintf(stderr,"%.6lf\n", time_stop(t_start, c_start)); 
	#endif
	
	for(size_t i=1; i < da.size() + max_da; ++i){
		size_t d  = i < da.size() ? da[i] : i-da.size();
		size_t lb = last_occ[d];
		size_t rb = i < da.size() ? i-1 : da.size()-1;
		last_occ[d] = i+1;
		#if DEBUG
			cout << "i="<<setw(2)<<i<<" d="<<setw(2)<<d<<" ["<<lb<<","<<rb<<"] ";
		#endif
		
		for(size_t j=lb; j<=rb; ++j){
			if ( seen[da[j]] != i ) {
				seen_stack.push(da[j]);
				seen[da[j]] = i;
				freq[da[j]] = 1;
			} else {
				++freq[da[j]];
			}
		}
		#if DEBUG
		  cout << " unique docs: " << seen_stack.size() << " ";
		#endif
		
		while( !seen_stack.empty() ) {
			auto x = seen_stack.top();
			seen_stack.pop();
			#if DEBUG
				cout << " (d="<<x<<", f="<<freq[x]<<")";
			#endif
			auto i1 = std::min(d,x);
			auto i2 = std::max(d,x);
			++counts[i1][i2][freq[x]];
			++runs[i1][i2];
		}
		#if DEBUG
		  cout << endl;
		#endif
	}
	
	for(size_t i=0; i<=max_da; ++i){
		for(size_t j=i+1; j<=max_da; ++j) {
			Result(i,j) = compute_distance(counts[i][j], runs[i][j]);
		}
	}
	
	#if TIME
		printf("#3. BWSD-RMQ:\n");
		fprintf(stderr,"%.6lf\n", time_stop(t_start, c_start)); 
	#endif
	
	#if DEBUG
		#if OUTPUT
			for(int_t i=0; i<k; i++){
				for(int_t j=0; j<i; j++){
					printf("%.2lf\t", Result(j,i));
				}
				printf("0.00\t");
				for(int_t j=i+1; j<k; j++){
					printf("%.2lf\t", Result(i,j));
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
				sum+=Result(i,j);
	
		printf("checksum = %lf\n",sum);
	#endif
    
	#if OUTPUT
		delete[] Md;
	#endif

return 0;
}
/******************************************************************************/
