#include <vector>
#include <numeric>
#include <stack>
#include <eigen34/Eigen/Eigen>

#include "helpfunctions.hpp"

#include "solve_modulo.hpp"



// computes the gcd of the entries of a vector M
int gcdVec(Eigen::Matrix<int, -1, 1> &M) {
    int d = M(0,0);
    int r = M.rows();
    for(int i = 1; i < r; i++) {
        d = gcd(d, M(i, 0));
    }
    return d;
}

// finds one solution to ax=b mod n by bruteforce
int axbmodnOneSol(int a, int b, int n) {
    for(int i = 0; i < n; i++) {
        if( (a*i) % n == b) return i;
    }
    return -1;      //shouldn't happen
}

// finds all solutions to ax=b mod n
// assumes there is at least one solution
std::vector<int> axbmodn(int a, int b, int n) {
    while(b<0) b +=n;
    int x0 = axbmodnOneSol(a, b, n);
    std::vector<int> v;
    int d = gcd(a, n);
    v.reserve(d);
    v.push_back(x0);
    for(int t = 1; t < d; t++) {
        v.push_back(x0 + n*t/d);
    }
    return v;
}

bool checkIfRedundant(std::vector<Eigen::Matrix<int, -1, 1>> &L, Eigen::Matrix<int, -1, 1> &M, int Delta, bool generic) {
    if(!generic) {
        // in the non-generic case we only have to check three things:
        // M is not Delta*unit vector
        int r = M.rows();
        std::vector<Eigen::Matrix<int, -1, 1>> units;
        for(int i = 0; i < r; i++) {
            Eigen::Matrix<int, -1, 1> U = Eigen::Matrix<int, -1, 1>::Zero(r);
            U(i, 0) = Delta;
            units.push_back(U);
            units.push_back(-U);
        }
        if( std::find(units.begin(), units.end(), M) != units.end() ) return true;

        // M is not zero
        if(M.isZero()) return true;

        // M or -M are included already
        return (std::find(L.begin(), L.end(), M) != L.end()) || (std::find(L.begin(), L.end(), -M) != L.end());
    }
    // in the generic case we check whether M is a multiple of a column already included
    int d = gcdVec(M);
    for(int i = 1; i <= d; i++) {
        if(d%i == 0) {
            if(std::find(L.begin(), L.end(), M/i) != L.end()) {
                return true;
            }
        }
    }
    return false;
}

// finds all representants in Z^n of the vector v which has elements in Z_Delta
// so that the entries are nonzero, last entry is positive in generic case
// or just last entry non-negative in the non-generic case
void allRepresentants(std::vector<int> v, std::vector<Eigen::Matrix<int, -1, 1>> &L, int Delta, bool generic) {
    int r = v.size();    
    std::vector<std::vector<int>> w;
    // because we solve the system by backwards substitution, the first entry of v will be the last entry of the column
    // restrict that entry to be non-negative
    // in the non-generic case that does not yet guarantee that only one of M or -M will be included, as they both can
    // have the last entry be equal to 0
    if(v[0] == 0) {
        if(generic) w.push_back({Delta});
        else w.push_back({Delta, 0});
    }
    else {
        w.push_back({v[0]});
    }
    // all options for the other entries
    for(unsigned int i = 1; i < v.size(); i++) {
        if(v[i] == 0) {
            if(generic) w.push_back({Delta, -Delta});
            else w.push_back({Delta, 0, -Delta});
        }
        else {
            w.push_back({v[i], v[i]-Delta});
        }
    }
    // combine
    std::vector<std::vector<int>> K = cartesian(w);
    for(auto I: K) {
        Eigen::Matrix<int, -1, 1> M(r, 1);
        // M = (I_r, ..., I_1) where I_k = v_k mod Delta
        for(int i = 0; i < r; i++) {
            M(i, 0) = I[r-i-1];
        } 
        if(!checkIfRedundant(L, M, Delta, generic)) {
            L.push_back(M);
        }
    }
}

// Solves Ax=0 mod n with matrix A and returns all relevant integer solutions
std::vector<Eigen::Matrix<int,-1,1>> Ax0modn(Eigen::MatrixXi &A, int n, bool generic) {
    int r = A.cols();
    std::vector<Eigen::Matrix<int,-1,1>> L;
    L.reserve(std::pow(2, r-1)*n);
    
    // backwards substitution
    std::vector<int> w = axbmodn(A(r-1, r-1), 0, n);
    std::stack<std::vector<int>> buildingSols;
    for(int u: w) {
        buildingSols.push({u});
    }

    std::vector<int> v;
    while(!buildingSols.empty()) {
        // v = [v_1, ..., v_k] corresponds to solution (*, .., *, v_k, ..., v_1) where values * have not been computed yet
        // Compute one more value, then put back on stack or add to list L of all solutions
        v = buildingSols.top();
        buildingSols.pop();
        int k = v.size();
        int s = 0;
        for(int i = 0; i < k; i++) {
            s += -A(r-k-1, r-1-i)*v[i];
        }
        w = axbmodn(A(r-k-1, r-k-1), s, n);
        for(int u: w) {
            v.push_back(u);
            if(k+1 == r) {
                allRepresentants(v, L, n, generic);
            }
            else {
                buildingSols.push(v);
            }
            v.pop_back();
        }
    }
    return L;
}