#include <iostream>
#include <string>
#include <map>
#include <sdsl/rmq_support.hpp>

using namespace std;
using namespace sdsl;

typedef unordered_map<uint64_t,uint64_t> tUMII;

double compute_distance(const tUMII& m, size_t runs) {
    double dist = 0;
    for(const auto& it : m) {
        double tmp = (double)it.second / (double) runs;
        dist += tmp*log2(tmp);
    }
    dist *= -1.0;
    return dist;
}

int main(){
    vector<uint64_t> D = {0,1,2,3,2,2,1,2,2,3,4,2,1,2,3,4,2,1,2,1,0,5,4,3,2,3,2,1,3,2,3,4,5,3,2,3,3,1,1,1,2,0};
    // assume that D contains all values in 0..maxD-1
    auto maxD = *std::max_element(D.begin(), D.end());

    cout << maxD << endl;

    auto print_array = [](vector<uint64_t>& vec, string label) {
        cout << label << ":";
        for(const auto& x : vec) {
            cout << " " << setw(2) << x;
        }
        cout << endl;
    }; 

    print_array(D, "D  ");


    vector<uint64_t> P(D.size(), D.size());
    vector<uint64_t> N(D.size(), 0);
    vector<uint64_t> R(D.size(), 0);
    {
        vector<uint64_t> last_occ(maxD+1, 0);
        for (size_t i=0; i < D.size(); ++i) {
            P[i] = last_occ[D[i]];
            last_occ[D[i]] = i+1;
        }
    }
    print_array(P, "P+1");
    {
        vector<uint64_t> next_occ(maxD+1, D.size());
        for(size_t i=D.size(); i > 0; --i) {
            N[i-1] = next_occ[D[i-1]];
            next_occ[D[i-1]] = i;
        }
    }
    print_array(N, "N  ");
    {
        for(size_t i=0; i<D.size(); ++i){
            if ( P[i] > 0 ) {
                R[i] = R[P[i]-1]+1;
            }
        }
    }
    print_array(R, "R  ");

    rmq_succinct_sct<true>  rmq_P(&P);
    rmq_succinct_sct<false> RMQ_N(&P);

    vector<size_t> seen(maxD+1,0); 
    stack<size_t>  seen_stack;
    vector<size_t> freq(maxD+1,0); 
    vector<size_t> last_occ(maxD+1, 0);

    vector<vector<tUMII>> counts(maxD+1, vector<tUMII>(maxD+1));
    vector<vector<uint64_t>> runs(maxD+1, vector<uint64_t>(maxD+1));

    for(size_t i=1; i < D.size() + maxD; ++i){
        size_t d  = i < D.size() ? D[i] : i-D.size();
        size_t lb = last_occ[d];
        size_t rb = i-1;
        last_occ[d] = i+1;
        cout << "i="<<setw(2)<<i<<" d="<<setw(2)<<d<<" ["<<lb<<","<<rb<<"] ";
        for(size_t j=lb; j<=rb; ++j){
            if ( seen[D[j]] != i ) {
                seen_stack.push(D[j]);
                seen[D[j]] = i;
                freq[D[j]] = 1;
            } else {
                ++freq[D[j]];
            }
        }
        cout << " unique docs: " << seen_stack.size() << " ";
        while( !seen_stack.empty() ) {
            auto x = seen_stack.top();
            seen_stack.pop();
            cout << " (d="<<x<<", f="<<freq[x]<<")";
            auto i1 = std::min(d,x);
            auto i2 = std::max(d,x);
            ++counts[i1][i2][freq[x]];
            ++runs[i1][i2];
        }
        cout << endl;
    }

    for(size_t i=0; i<=maxD; ++i){
        for(size_t j=i+1; j<=maxD; ++j) {
            double dist = compute_distance(counts[i][j], runs[i][j]);
            cout << dist << " "; 
        }
        cout << endl;
    }
   


    return 0;
}
