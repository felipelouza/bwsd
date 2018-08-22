// vim: noai:ts=2:sw=2

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
#include <fstream>

#include <omp.h>

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

#ifndef OMP
	#define OMP 1
#endif

#ifndef DEBUG
	#define DEBUG 0 
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

int compute_all_bwsd_rank(unsigned char** R, uint_t k, uint_t n, char* c_file, int dist, int output, int check, int print, int verbose);//algorithm 1

int compute_all_bwsd(unsigned char** R, uint_t k, uint_t n, char* c_file, int dist, int output, int check, int print, int verbose);//straightforward

int compute_all_bwsd_rmq_Nz(unsigned char **R, uint_t k, uint_t n, char* c_file, int dist, int output, int check, int print, int verbose);//algorithm 2

int compute_all_bwsd_rmq_Nk(unsigned char **R, uint_t k, uint_t n, char* c_file, int dist, int output, int check, int print, int verbose);//algorithm 2

/******************************************************************************/

void usage(char *name){
  printf("\n\tUsage: %s [options] FILE d\n\n",name);
  puts("Computes BWSD-based distances for the first d sequences of a collection.");
//  puts("Sequences from FILE are extracted according to FILE's");
  puts("Extensions currently supported extensions: .txt .fasta .fastq\n");
  puts("Available options:");
  puts("\t-h    this help message");
  puts("\t-A a  preferred algorithm to use (default is alg. 1 BIT_sd)");
  puts("\t-B b  BWSD-based distance to compute, options: 1. expectation (default), 2. Shannon entropy");
	puts("\t-T t  use t parallel threads (default 1)");
  puts("\t-o    write output matrix to FILE.output.bin");
  puts("\t-p    print the output matrix (for debug)");
  puts("\t-v    verbose output\n");
  exit(EXIT_FAILURE);
}

/******************************************************************************/


int main(int argc, char** argv){

	#if OMP
		double omp_start=0.0;
	#else
		time_t t_start=0;clock_t c_start=0;
	#endif

	extern char *optarg;
	extern int optind;
	int c;
	char *c_dir=NULL, *c_file=NULL, *c_input=NULL;

	int verbose=0, check=0, print=0, dist=1;
	int MODE=1;//preferred algorithm
	int k=0;
	int output=0; //outputs the matrix to FILE.output.bin
	#if OMP
		int n_threads=1;
	#endif
	
	while ((c=getopt(argc, argv, "vcpA:B:T:ho")) != -1) {
		switch (c) {
			case 'v':
				verbose++; break;
			case 'c':
				check++; break;
			case 'p':
				print++; break;
			case 'A':
				MODE = (int) atoi(optarg); break;
			case 'B':
				dist = (int) atoi(optarg); break;
			#if OMP
			case 'T':
				n_threads = (int) atoi(optarg); break;
			#endif
			case 'h':
				usage(argv[0]); break;      // show usage and stop
			case 'o':
				output++; break;
			case '?':
				exit(EXIT_FAILURE);
		}
	}

	if(optind+2==argc) {
		c_input=argv[optind++];
		k = (int) atoi(argv[optind++]);
	}
	else  usage(argv[0]);

	c_file= strrchr(c_input, '/')+1;
	c_dir = strndup(c_input, strlen(c_input)-strlen(c_file));

	unsigned char **R;
	int_t n=0;

	file_chdir(c_dir);

	//disk access
	R = (unsigned char**) file_load_multiple(c_file, &k, &n);
	if(!R){
		fprintf(stderr, "Error: less than %d strings in %s\n", k, c_file);
		return 0;
	}

	if(verbose){
		printf("########\n");
		printf("DIR = %s\n", c_dir);
		printf("FILE = %s\n", c_file);
		if(output) printf("OUTPUT = %s.output.bin\n", c_file);
		printf("Algorithm = %d\n", MODE);
		printf("Strings = %d\n", k);
		printf("N = %d\n", n);
		printf("BWSD-dist: ");
		if(dist==1) printf("Expectation (default)\n");
		else if (dist==2) printf("Shannon entropy\n");
		printf("sizeof(int) = %zu bytes\n", sizeof(int_t));
	}

	#if OMP
		omp_set_num_threads(n_threads);
	
		#pragma omp parallel
		{
		if(omp_get_thread_num()==0)
			if(verbose) printf("N_THREADS: %d\n", omp_get_num_threads());
		}
		if(verbose) printf("N_PROCS: %d\n", omp_get_num_procs());
	#endif
	
	if(verbose){
		printf("########\n");
	}

	#if OMP
		omp_start = omp_get_wtime();
	#else 
		time_start(&t_start, &c_start);
	#endif

	switch(MODE){

		case 1:  //Algorithm 1, O(Nk) time
			#if WT
				printf("## BWSD_WT ##\n");
			#elif SD_VECTOR
				printf("## BWSD_BIT_sd ##\n"); 
			#else
				printf("## BWSD_BIT ##\n");
			#endif
			compute_all_bwsd_rank(R, k, n, c_file, dist, output, check, print, verbose);
			break;

		case 2:	printf("## BWSD_RMQ_Nk ##\n"); //Algorithm 2, O(Nk) time
			compute_all_bwsd_rmq_Nk(R, k, n, c_file, dist, output, check, print, verbose);
			break;

		case 3:	printf("## BWSD_SF ##\n"); //Straightforward algorithm, O(Nk) time
			compute_all_bwsd(R, k, n, c_file, dist, output, check, print, verbose);
			break;
		
		case 4:	printf("## BWSD_RMQ_Nz ##\n"); //Algorithm 2, O(N+z) time
			compute_all_bwsd_rmq_Nz(R, k, n, c_file, dist, output, check, print, verbose);
			break;


		default: printf("ERROR: please, choose a valid -A options 1, 2, 3 and 4.\n");
			break;
	}
	
	printf("TOTAL:\n");
	#if OMP
		printf("TIME = %.6lf (in seconds)\n", omp_get_wtime()-omp_start);
		fprintf(stderr, "%.6lf\n", omp_get_wtime()-omp_start);
	#else
		fprintf(stderr,"%.6lf\n", time_stop(t_start, c_start));
	#endif

return 0;
}

/******************************************************************************/

unsigned char* cat(unsigned char* s1, unsigned char* s2, uint_t *n){

	(*n) = (strlen((char*)s1)+1)+(strlen((char*)s2)+1)+1; //add 0 at the end

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

unsigned char* cat_all(unsigned char** R, int k, uint_t *n, int verbose){

	(*n)++; //add 0 at the end
	int_t i, j, l=0;
	unsigned char *str = (unsigned char*) malloc((*n)*sizeof(unsigned char));

	int_t max=0;

	for(i=0; i<k; i++){
		int_t m = strlen((char*)R[i]);
		if(m>max) max=m;
		for(j=0; j<m; j++){
			if(R[i][j]<255) str[l++] = R[i][j]+1;
		}
		str[l++] = 1; //add 1 as separator
	}

	str[l++]=0;
	if(*n>l){
	  str = (unsigned char*) realloc(str, (l)*sizeof(unsigned char)); 
	}
	*n = l;

	if(verbose)	cout<<"longest string = "<<max<<endl;

return str;
}


/******************************************************************************/

double bwsd_expectation(tMII &t, int_t s){

	double value = 0.0;
	
	for(tMII::iterator it=t.begin(); it!=t.end(); ++it){
		if(it->second){
			value += (double)it->first*(double)it->second/(double)s;
		}
	}

	value = value-1.0;

return value;
}

/******************************************************************************/

double bwsd_shannon_entropy(tMII &t, int_t s){

	double value = 0.0;
	
	for(tMII::iterator it=t.begin(); it!=t.end(); ++it){
		if(it->second){
			double tmp = (double)it->second/(double)s;
			value += tmp*log2(tmp);
		}
	}

	if(value) value *= -1.0;

return value;
}
/******************************************************************************/

double compute_distance(tMII &t, int_t s, int dist){

	//implemented distances
	switch (dist) {
		case 1:
			return bwsd_expectation(t, s); break;
		case 2:
			return bwsd_shannon_entropy(t, s); break;
		default: break;
	}

return 0.0;
}
/******************************************************************************/

void push_if_not_empty(stack<tAII>& s, tAII range){
    if ( range[0] < range[1]+1 ) {
        s.push(range);
    }
}

/******************************************************************************/

int	print_output(double *Md, uint_t k){

	double zero=0.0;

	for(int_t i=0; i<k; i++){
		for(int_t j=0; j<i; j++){
			printf("%.3lf\t", Result(j,i));
		}
		printf("%.3lf\t", zero);
		for(int_t j=i+1; j<k; j++){
			printf("%.3lf\t", Result(i,j));
		}
		cout<<endl;
	}

return 0;
}

/******************************************************************************/

int	write_output(char* c_file, double* Md, size_t m){
	
	char c_out[256];
	sprintf(c_out, "%s.output.bin", c_file);
	ofstream f_out(c_out, ios::out | ios::binary);

	f_out.seekp(0);
	for(int i=0; i<m; i++) f_out.write(reinterpret_cast<char*>(&Md[i]), sizeof(double));
	f_out.close();

	cout << "writing " << m*sizeof(double) << " bytes to: "<<c_out<<"\n";

return 1;
}

/******************************************************************************/

int	read_output(char* c_file, double* Md, size_t m, uint_t k){

	char c_out[256];
	sprintf(c_out, "%s.output.bin", c_file);
	ifstream f_in(c_out, ios::in | ios::binary);

	cout << f_in.gcount() << " bytes read\t("<<m<<")\n";
	f_in.seekg(0);
	for(int i=0; i<m; i++) f_in.read(reinterpret_cast<char*>(&Md[i]), sizeof(double));

return 0;
}
			
/******************************************************************************/

int compute_all_bwsd_rank(unsigned char** R, uint_t k, uint_t n, char* c_file, int dist, int output, int check, int print, int verbose){

	int_t i;

	//Concatenate strings
	/**/
	unsigned char *str = cat_all(R, k, &n, verbose);

	#if OPT_VERSION 
		if(verbose) cout<<"OPT_VERSION"<<endl;
	#endif

	#if WORST_CASE 
		cout<<"WORST_CASE"<<endl;
	#endif

	#if DEBUG
		printf("R:\n");
		for(i=0; i<k; i++)
			printf("%" PRIdN ") %s (%zu)\n", i, R[i], strlen((char*)R[i]));
	#endif

	//free memory
	#if OMP
		#pragma omp parallel for 
	#endif
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

	#if OMP
		double omp_start=0.0;
		if(verbose) omp_start = omp_get_wtime();
	#else
		time_t t_start=0;clock_t c_start=0;
		if(verbose) time_start(&t_start, &c_start); 
	#endif

	//COMPUTE DA:
	/**/
	if(!load_from_cache(da, "da", m_config)){

		int_t *SA = new int_t[n];
		int_t *DA = new int_t[n];
		
		#if OMP
			#pragma omp parallel for 
		#endif
		for(i=0; i<n; i++) SA[i]=DA[i]=0;
		gsacak(str, (uint_t*)SA, NULL, DA, (uint_t)n); //construct SA+DA
	
		#if OMP
			#pragma omp parallel for 
		#endif
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
		da[j2-1]=k-1;
		n=j2;
		da.resize(n);
		store_to_cache(da, "da", m_config);
		#if DEBUG
			cout<<"DA: "; 
			for(int_t i=1; i<n; i++){
					if(da[i]!=da[i-1]) cout<<endl; 
					cout<<da[i]<<" ";
			}
			cout<<endl;		
		#endif
	#endif


	size_t m = (k*k-k)/2.0;
	double *Md = new double[m];
 
	if(verbose){
		printf("#1. DA:\n");
		#if OMP
			printf("TIME = %.6lf (in seconds)\n", omp_get_wtime()-omp_start);
			fprintf(stderr, "%.6lf\n", omp_get_wtime()-omp_start);
		#else
			fprintf(stderr,"%.6lf\n", time_stop(t_start, c_start));
		#endif
	}

	#if WT
		//COMPUTE WT(DA):
		wm_int<> wt;
		#if WORST_CASE
			construct_im(wt, da);
			store_to_cache(wt, "wt", m_config);
		#else
			if(!load_from_cache(wt, "wt", m_config)){
				construct_im(wt, da);
				store_to_cache(wt, "wt", m_config);
			}
		#endif
	#else

		#if SD_VECTOR
			sd_vector_builder *B = new sd_vector_builder[k+1];
		
			int_t *size = new int_t[k+1];
			#if OMP
				#pragma omp parallel for 
			#endif
			for(int_t i=0; i<k; i++) size[i]=0;
			for(int_t i=1; i<n; i++) size[da[i]]++;

			#if OMP
				#pragma omp parallel for 
			#endif
			for(int_t i=0; i<k; i++)
				B[i] = sd_vector_builder(n,size[i]);

			delete[] size;
		#else
			bit_vector *B = new bit_vector[k];
			#if OMP
				#pragma omp parallel for 
			#endif
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
			#if OMP
				#pragma omp parallel for 
			#endif
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

		#if OMP
			#pragma omp parallel for 
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
			#if OMP
				#pragma omp parallel for 
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

	#if WT
		if(verbose) printf("#2. WT(DA):\n");
	#else
		#if SD_VECTOR
			if(verbose) printf("#2. BITMAP_SD(DA):\n");
		#else
			if(verbose) printf("#2. BITMAP(DA):\n");
		#endif
	#endif

	if(verbose){
		#if OMP
			printf("TIME = %.6lf (in seconds)\n", omp_get_wtime()-omp_start);
			fprintf(stderr, "%.6lf\n", omp_get_wtime()-omp_start);
		#else
			fprintf(stderr,"%.6lf\n", time_stop(t_start, c_start));
		#endif
	}

	#if DEBUG	
		cout<<"DA:"<<endl;
		for(int_t i=0; i<n-1; i++) cout << da[i] << ", ";
		cout<<da[n-1]<<"\t("<<n-1<<")"<<endl;
	#endif



	#if OPT_VERSION 
	
		int_t *rank = new int_t[k+1];
		//avoid wt.rank()-queries when next > qe
		int_t **pos = new int_t*[k+1];
		
		#if OMP
			#pragma omp parallel for 
		#endif
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
		
		#if OMP
			#pragma omp parallel for 
		#endif
		for(int_t i=0; i<k; i++) pos[i][rank[i]]=n+1;

		uint64_t skip=0, total=1;

		delete[] rank;
	#endif


	#if OMP
		#if OPT_VERSION 
				#pragma omp parallel for reduction(+:skip) reduction(+:total) 
		#else
				#pragma omp parallel for 
		#endif
	#endif
	for(int_t i=0; i<k-1; i++){

		int_t *s= new int_t[k];
		int_t *ell = new int_t[k];
	
		//avoid second wt.rank()
		int_t *rank = new int_t[k+1];

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

		#if OMP
			#pragma omp parallel for 
		#endif
		for(int_t j=i+1; j<k; j++) s[j] = rank[j] = ell[j]=0;
	
		//forloop S^i[1..n_i]
		for(uint64_t p=1; p<len_i+1; p++){

			#if DEBUG == 2
				#if WT
					cout << "[" <<qe <<", "<< wt.select(p, i) << "]   \t0^1"<<endl;
				#endif
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
				cout << "[.." <<", "<<n<< "]"<<endl;
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

				Result(i,j) = compute_distance(t[j], s[j], dist);
				#if DEBUG
					cout<<"["<<i<<", "<<j<<"]\t\ts="<<s[j]<<"\n";
					cout<<"D = "<<Result(i,j)<<endl;
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

				cout<<"####\t";
			}		
			cout<<endl;
		
		#endif
	
		delete[] ell;
		delete[] s;
		delete[] rank;
	}

	if(verbose){
		printf("#3. BWSD-RANK:\n");
		#if OMP
			printf("TIME = %.6lf (in seconds)\n", omp_get_wtime()-omp_start);
			fprintf(stderr, "%.6lf\n", omp_get_wtime()-omp_start);
		#else
			fprintf(stderr,"%.6lf\n", time_stop(t_start, c_start));
		#endif
	}

	#if OPT_VERSION
		#if OMP
			#pragma omp parallel for 
		#endif
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

	if(output){
		write_output(c_file, Md, m);
		#if DEBUG
			for(int i=0; i<m; i++) Md[i]=0.0;
			read_output(c_file, Md, m, k);
		#endif
	}
	
	if(print){
		print_output(Md, k);
	}

	//checksum: for the sake of sanity
	if(check){
		double sum=0.0;
		for(int i=0; i<m; i++) sum+=Md[i];

		printf("checksum = %lf\n",sum);
	}

	delete[] Md;

	#if OPT_VERSION
		if(verbose) cout<<"skip = "<<skip<<" / "<<total<<" = "<<(double)skip/(double)total<<endl;
	#endif

return 0;
}

/******************************************************************************/

int compute_all_bwsd(unsigned char** R, uint_t k, uint_t n, char* c_file, int dist, int output, int check, int print, int verbose){ //straightforward

	#if WORST_CASE 
		if(verbose) cout<<"WORST_CASE"<<endl;
	#endif

	#if DEBUG
		printf("R:\n");
		for(int_t i=0; i<k; i++)
			printf("%" PRIdN ") %s (%zu)\n", i, R[i], strlen((char*)R[i]));
	#endif

	//tVMID	result(k);
	size_t m = (k*k-k)/2.0;
	double *Md = new double[m];

	/**/
	#if OMP
		#pragma omp parallel for 
	#endif
	for(int_t i=0; i<k-1; i++){

		#if DEBUG 
			cout<<endl;
		#endif


		//foreach S^j in j+1..k
		#if OMP
			#pragma omp parallel for 
		#endif
		for(int_t j=i+1; j<k; j++){
	
			tMII t;

			int_t s=0;
			uint_t n=0;
			//concatenates
			unsigned char *str = cat(R[i], R[j], &n);

			//build DA
			int_t *SA = new int_t[n];
			int_t *DA = new int_t[n+1];
			for(int_t p=0; p<n; p++) SA[p]=DA[p]=0;
			gsacak(str, (uint_t*)SA, NULL, DA, (uint_t)n); //construct SA+DA

			#if WORST_CASE
				if(n%2==0) n=n-1;
				int_t j1=1;
				int_t j2=n/2+1;
				for(int_t i=j1; i<j2; i++) DA[i]=0;
				for(int_t i=j2; i<n;	i++) DA[i]=1;
				#if DEBUG
					cout<<"DA: "; 
					for(int_t i=1; i<n; i++){
						if(DA[i]!=DA[i-1]) cout<<endl; 
						cout<<DA[i]<<" ";
					}
					cout<<endl;
				#endif
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

			Result(i,j) = compute_distance(t, s, dist);

			delete[] DA;
			delete[] SA;
		}//for j=i+1..n

		#if DEBUG
			cout<<"####\n";
		#endif
	}

	//free memory
	#if OMP
		#pragma omp parallel for 
	#endif
	for(int_t i=0; i<k; i++)
		free(R[i]);
	free(R);

	if(output){
		write_output(c_file, Md, m);
		#if DEBUG
			for(int i=0; i<m; i++) Md[i]=0.0;
			read_output(c_file, Md, m, k);
		#endif
	}

	if(print){
		print_output(Md, k);
	}

	//checksum: for the sake of sanity
	if(check){
		double sum=0.0;
		for(int i=0; i<m; i++) sum+=Md[i];

		printf("checksum = %lf\n",sum);
	}
	delete[] Md;

return 0;
}

/******************************************************************************/

int compute_all_bwsd_rmq_Nk(unsigned char** S, uint_t k, uint_t n, char* c_file, int dist, int output, int check, int print, int verbose){

	int_t i;
	
	//Concatenate strings
	/**/
	unsigned char *str = cat_all(S, k, &n, verbose);
	
	#if DEBUG
		printf("R:\n");
		for(i=0; i<k; i++)
			printf("%" PRIdN ") %s (%zu)\n", i, S[i], strlen((char*)S[i]));
	#endif
	
	#if OPT_VERSION 
		if(verbose) cout<<"OPT_VERSION"<<endl;
	#endif

	#if WORST_CASE 
		cout<<"WORST_CASE"<<endl;
	#endif

	//free memory
	#if OMP
		#pragma omp parallel for 
	#endif
	for(i=0; i<k; i++)
		free(S[i]);
	free(S);
	
	/**/
	
	string dir = "sdsl";
	mkdir(dir.c_str());
	string id = c_file;
	id += "."+to_string(k);
	
	cache_config m_config(true, dir, id);
	int_vector<64> da(n);
	//vector<uint64_t> da(n);
	
	#if OMP
		double omp_start=0.0;
		if(verbose) omp_start = omp_get_wtime();
	#else
		time_t t_start=0;clock_t c_start=0;
		if(verbose) time_start(&t_start, &c_start); 
	#endif
	
	//COMPUTE DA:
	/**/
	if(!load_from_cache(da, "da_rmq", m_config)){
	
		int_t *SA = new int_t[n];
		int_t *DA = new int_t[n];
		
		#if OMP
			#pragma omp parallel for 
		#endif
		for(i=0; i<n; i++) SA[i]=DA[i]=0;
		gsacak(str, (uint_t*)SA, NULL, DA, (uint_t)n); //construct SA+DA
	
		#if OMP
			#pragma omp parallel for 
		#endif
		for(i=0;i<n;i++) da[i]=DA[i];
		store_to_cache(da, "da_rmq", m_config);
	
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
		da[j2-1]=k-1;
		n=j2;
		da.resize(n);
		#if DEBUG
			cout<<"DA: "; 
			for(int_t i=1; i<n; i++){
				if(da[i]!=da[i-1]) cout<<endl; 
				cout<<da[i]<<" ";
			}
			cout<<endl;		
		#endif
	#endif

	size_t m = (k*k-k)/2.0;
	double *Md = new double[m];
	
	if(verbose){
		printf("#1. DA:\n");
		#if OMP
			printf("TIME = %.6lf (in seconds)\n", omp_get_wtime()-omp_start);
			fprintf(stderr, "%.6lf\n", omp_get_wtime()-omp_start);
		#else
			fprintf(stderr,"%.6lf\n", time_stop(t_start, c_start));
		#endif
	}
	
	auto max_da = *std::max_element(da.begin(), da.end());

	#if DEBUG  
		cout << max_da << endl;
		auto print_array = [](int_vector<64>& vec, string label) {
			cout << label << ":";
			for(const auto& x : vec) {
				cout << " " << setw(2) << x;
			}
		  cout << endl;
		}; 
		print_array(da, "da  ");
	#endif  
	
	int_vector<64> P(da.size(), da.size());
	int_vector<64> N(da.size()+1, 0);
	int_vector<64> R(da.size(), 0);
	
	{
		int_vector<64> last_occ(max_da+1, 0);
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
		N[da.size()]=da.size()+max_da;
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
	
	//vector<size_t> seen(max_da+1,0); 
	

//	vector<vector<tMII>> counts(max_da+1, vector<tMII>(max_da+1));
//	vector<vector<int_t>> runs(max_da+1, vector<int_t>(max_da+1));  

	if(verbose){
		printf("#2. RMQ:\n");
		#if OMP
			printf("TIME = %.6lf (in seconds)\n", omp_get_wtime()-omp_start);
			fprintf(stderr, "%.6lf\n", omp_get_wtime()-omp_start);
		#else
			fprintf(stderr,"%.6lf\n", time_stop(t_start, c_start));
		#endif
	}
	
	#if OMP
		#pragma omp parallel for 
	#endif
	for(size_t d=0; d < max_da; d++){
	
		int_t *ell = new int_t[max_da+1];
		stack<size_t>  seen_stack;
		vector<size_t> freq(max_da+1,0); 
		vector<size_t> last_occ(max_da+1, 0);

		//init
		tVMII t(max_da+1);
		vector<int_t> runs(max_da+1);  

		size_t p=d+1;
		#if WORST_CASE
			p = d*(n/k)+1;
		#endif

		#if OPT_VERSION
			#if OMP
				#pragma omp parallel for 
			#endif
			for(int_t j=d+1; j<k; j++) ell[j]=0;
		#endif

		for(; p < da.size()+max_da; p=N[p]){

			i=p;

			//size_t d  = i < da.size() ? da[i] : i-da.size();

			size_t lb = last_occ[d];
			size_t rb = i < da.size() ? i-1 : da.size()-1;
			last_occ[d] = i+1;
			
			#if DEBUG == 2
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
					if(da[_m]>d){
						seen_stack.push(da[_m]);
						freq[da[_m]] = R[_m];
					}
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
					if(da[_m]>d){
						freq[da[_m]] = R[_m] - freq[da[_m]] + 1;
					}
					push_if_not_empty(ranges, {_m+1, _rb});
					push_if_not_empty(ranges, {_lb, _m-1});
				}
			}
			
			while( !seen_stack.empty() ) {
				auto j = seen_stack.top();
				seen_stack.pop();
				#if DEBUG 
				  cout << d << " (d="<<j<<", f="<<freq[j]<<")\n";
				#endif

				#if OPT_VERSION
					if(ell[j]){
						t[j][ell[j]]++;  //0^lj
						++runs[j];
						ell[j]=0;
					}
				#endif

				t[j][freq[j]]++;
				++runs[j];
			}

			#if DEBUG
			  cout << endl;
			#endif

			#if OPT_VERSION
				if(p<da.size())
					for(int_t j=d+1; j<k; j++) ell[j]++;//considers that S_j does not occur in [lb, rb]
			#endif
		}

		//count runs for S_d
		#if OMP
			#pragma omp parallel for 
		#endif
		for(size_t j=d+1; j<max_da; j++){

			#if OPT_VERSION
				if(ell[j]){
					t[j][ell[j]]++;  //0^lj
					++runs[j];
				}
			#else
				size_t next_d = d+1;
				size_t next_j = j+1;
	
				#if WORST_CASE
					next_d = d*(n/k)+1;
					next_j = j*(n/k)+1;
				#endif
	
				while(next_d<da.size()){
					size_t sum=0;
	
					while(next_d < next_j){
						sum++;
						next_d = N[next_d];
					}
					while(next_j < next_d){
						next_j = N[next_j];
					}
					t[j][sum]++;
					++runs[j];
				}
			#endif
			Result(d,j) = compute_distance(t[j], runs[j], dist);
		}

		delete[] ell;
	}
	
	if(verbose){
		printf("#3. BWSD-RMQ:\n");
		#if OMP
			printf("TIME = %.6lf (in seconds)\n", omp_get_wtime()-omp_start);
			fprintf(stderr, "%.6lf\n", omp_get_wtime()-omp_start);
		#else
			fprintf(stderr,"%.6lf\n", time_stop(t_start, c_start));
		#endif
	}
	
	if(output){
		write_output(c_file, Md, m);
		#if DEBUG
			for(int i=0; i<m; i++) Md[i]=0.0;
			read_output(c_file, Md, m, k);
		#endif
	}

	if(print){
		print_output(Md, k);
	}

	//checksum: for the sake of sanity
	if(check){
		double sum=0.0;
		for(int i=0; i<m; i++) sum+=Md[i];
	
		printf("checksum = %lf\n",sum);
	}
    
	delete[] Md;

return 0;
}

/******************************************************************************/

int compute_all_bwsd_rmq_Nz(unsigned char** S, uint_t k, uint_t n, char* c_file, int dist, int output, int check, int print, int verbose){

	int_t i;
	
	//Concatenate strings
	/**/
	unsigned char *str = cat_all(S, k, &n, verbose);

	#if WORST_CASE 
		cout<<"WORST_CASE"<<endl;
	#endif
	
	#if DEBUG
		printf("R:\n");
		for(i=0; i<k; i++)
			printf("%" PRIdN ") %s (%zu)\n", i, S[i], strlen((char*)S[i]));
	#endif
	
	//free memory
	#if OMP
		#pragma omp parallel for 
	#endif
	for(i=0; i<k; i++)
		free(S[i]);
	free(S);
	
	/**/
	
	string dir = "sdsl";
	mkdir(dir.c_str());
	string id = c_file;
	id += "."+to_string(k);
	
	cache_config m_config(true, dir, id);
	int_vector<64> da(n);
	//vector<uint64_t> da(n);
	
	#if OMP
		double omp_start=0.0;
		if(verbose) omp_start = omp_get_wtime();
	#else
		time_t t_start=0;clock_t c_start=0;
		if(verbose) time_start(&t_start, &c_start); 
	#endif
	
	//COMPUTE DA:
	/**/
	if(!load_from_cache(da, "da_rmq", m_config)){
	
		int_t *SA = new int_t[n];
		int_t *DA = new int_t[n];
		
		#if OMP
			#pragma omp parallel for 
		#endif
		for(i=0; i<n; i++) SA[i]=DA[i]=0;
		gsacak(str, (uint_t*)SA, NULL, DA, (uint_t)n); //construct SA+DA
	
		#if OMP
			#pragma omp parallel for 
		#endif
		for(i=0;i<n;i++) da[i]=DA[i];
		store_to_cache(da, "da_rmq", m_config);
	
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
		da[j2-1]=k-1;
		n=j2;
		da.resize(n);
		#if DEBUG
			cout<<"DA: "; 
			for(int_t i=1; i<n; i++){
					if(da[i]!=da[i-1]) cout<<endl; 
					cout<<da[i]<<" ";
			}
			cout<<endl;		
		#endif
	#endif

	size_t m = (k*k-k)/2.0;
	double *Md = new double[m];
	
	if(verbose){
		printf("#1. DA:\n");
		#if OMP
			printf("TIME = %.6lf (in seconds)\n", omp_get_wtime()-omp_start);
			fprintf(stderr, "%.6lf\n", omp_get_wtime()-omp_start);
		#else
			fprintf(stderr,"%.6lf\n", time_stop(t_start, c_start));
		#endif
	}
	
	auto max_da = *std::max_element(da.begin(), da.end());
	
	#if DEBUG  
		cout << max_da << endl;
		auto print_array = [](int_vector<64>& vec, string label) {
			cout << label << ":";
			for(const auto& x : vec) {
				cout << " " << setw(2) << x;
			}
		  cout << endl;
		}; 
		print_array(da, "da  ");
	#endif  
	
	int_vector<64> P(da.size(), da.size());
	int_vector<64> N(da.size()+1, 0);
	int_vector<64> R(da.size(), 0);
	
	{
		int_vector<64> last_occ(max_da+1, 0);
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
		N[da.size()]=da.size()+max_da;
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
	
	//vector<size_t> seen(max_da+1,0); 
	
	vector<vector<tMII>> counts(max_da+1, vector<tMII>(max_da+1));
	vector<vector<int_t>> runs(max_da+1, vector<int_t>(max_da+1));  
	
	if(verbose){
		printf("#2. RMQ:\n");
		#if OMP
			printf("TIME = %.6lf (in seconds)\n", omp_get_wtime()-omp_start);
			fprintf(stderr, "%.6lf\n", omp_get_wtime()-omp_start);
		#else
			fprintf(stderr,"%.6lf\n", time_stop(t_start, c_start));
		#endif
	}
	
//replaced by the two following forloops
//for(size_t i=1; i < da.size() + max_da; ++i){

	#if OMP
		#pragma omp parallel for 
	#endif
	for(size_t d=0; d < max_da; d++){
	
		stack<size_t>  seen_stack;
		vector<size_t> freq(max_da+1,0); 
		vector<size_t> last_occ(max_da+1, 0);

		size_t p=d+1;
		#if WORST_CASE
			p = d*(n/k)+1;
		#endif

		for(; p < da.size()+max_da; p=N[p]){

			i=p;
//		size_t d  = i < da.size() ? da[i] : i-da.size();
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
	}
	
	#if OMP
		#pragma omp parallel for 
	#endif
	for(size_t i=0; i<max_da; ++i){
		#if OMP
			#pragma omp parallel for 
		#endif
		for(size_t j=i+1; j<max_da; ++j) {
			Result(i,j) = compute_distance(counts[i][j], runs[i][j], dist);
		}
	}

	
	if(verbose){
		printf("#3. BWSD-RMQ:\n");
		#if OMP
			printf("TIME = %.6lf (in seconds)\n", omp_get_wtime()-omp_start);
			fprintf(stderr, "%.6lf\n", omp_get_wtime()-omp_start);
		#else
			fprintf(stderr,"%.6lf\n", time_stop(t_start, c_start));
		#endif
	}
	
	if(output){
		write_output(c_file, Md, m);
		#if DEBUG
			for(int i=0; i<m; i++) Md[i]=0.0;
			read_output(c_file, Md, m, k);
		#endif
	}

	if(print){
		print_output(Md, k);
	}
	
	//checksum: for the sake of sanity
	if(check){
		double sum=0.0;
		for(int i=0; i<m; i++) sum+=Md[i];
	
		printf("checksum = %lf\n",sum);
	}
    
	delete[] Md;

return 0;
}


/******************************************************************************/

