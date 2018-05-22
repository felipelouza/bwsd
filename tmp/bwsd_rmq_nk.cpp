#include <iostream>
#include <string>
#include <map>
#include <sdsl/rmq_support.hpp>

using namespace std;
using namespace sdsl;

typedef unordered_map<uint64_t,uint64_t> tUMII;

double compute_distance(const tUMII& m, size_t runs, bool show=false) {
    double dist = 0;
    for(const auto& it : m) {
        if (show ) {
            cout << " ("<<it.second<<"/"<<runs<<")";
        }
        double tmp = (double)it.second / (double) runs;
        dist += tmp*log2(tmp);
    }
    dist *= -1.0;
    return dist;
}

int main(){
    vector<uint64_t> da = {0,1,2,3,2,2,1,2,2,3,4,2,1,2,3,4,2,1,2,1,0,5,4,3,2,3,2,1,3,2,3,4,5,3,2,3,3,1,1,1,2,0};
    // assume that D contains all values in 0..maxD-1
    auto max_da = *std::max_element(da.begin(), da.end());

    cout << max_da << endl;

    auto print_array = [](vector<uint64_t>& vec, string label) {
        cout << label << ":";
        for(const auto& x : vec) {
            cout << " " << setw(2) << x;
        }
        cout << endl;
    }; 

    print_array(da, "DA  ");


    vector<uint64_t> P(da.size(), da.size());
    vector<uint64_t> N(da.size(), 0);
    vector<uint64_t> R(da.size(), 0);
    {
        vector<uint64_t> last_occ(max_da+1, 0);
        for (size_t i=0; i < da.size(); ++i) {
            P[i] = last_occ[da[i]];
            last_occ[da[i]] = i+1;
        }
    }
    print_array(P, "P+1");
    {
        for(size_t i=0; i<da.size(); ++i){
            if ( P[i] > 0 ) {
                R[i] = R[P[i]-1]+1;
            }
        }
    }
    print_array(R, "R  ");

    vector<size_t> seen(max_da+1,0); 
    stack<size_t>  seen_stack;
    vector<size_t> freq(max_da+1,0); 
    vector<size_t> last_occ(max_da+1, 0);

    vector<vector<tUMII>> counts(max_da+1, vector<tUMII>(max_da+1));
    vector<vector<uint64_t>> runs(max_da+1, vector<uint64_t>(max_da+1));

    for(size_t i=1; i < da.size() + max_da; ++i){
        size_t d  = i < da.size() ? da[i] : i-da.size();
        size_t lb = last_occ[d];
        size_t rb = i < da.size() ? i-1 : da.size()-1;
        last_occ[d] = i+1;
        cout << "i="<<setw(2)<<i<<" d="<<setw(2)<<d<<" ["<<lb<<","<<rb<<"] ";
        for(size_t j=lb; j<=rb; ++j){
            if ( seen[da[j]] != i ) {
                seen_stack.push(da[j]);
                seen[da[j]] = i;
                freq[da[j]] = 1;
            } else {
                ++freq[da[j]];
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

    cout << "       ";
    for(size_t i=0; i<=max_da; ++i) {
        cout << setw(6) << i << " ";
    }
    cout << endl;
    for(size_t i=0; i<=max_da; ++i){
        cout << setw(6) << i << " "; 
        size_t w = 6+7*(i+1);
        for(size_t j=i+1; j<=max_da; ++j) {
            double dist = compute_distance(counts[i][j], runs[i][j]);
            cout << setw(w) << setprecision(3) << dist << " "; 
            w = 6;
        }
        cout << endl;
    }
    cout << endl;

    for(size_t i=0; i<=max_da; ++i){
        for(size_t j=i+1; j<=max_da; ++j) {
            cout << i << " " << j <<":";
            double dist = compute_distance(counts[i][j], runs[i][j], true);
            cout << endl;
        }
    }
   


    return 0;
}
