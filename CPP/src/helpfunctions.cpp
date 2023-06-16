#include <vector>
#include <unordered_set>
#include <eigen34/Eigen/Eigen>
#include <iostream>
#include <numeric>



#include "helpfunctions.hpp"



const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, 0, ", ", "\n");

void printVector(std::vector<int> v, std::ostream& stream) {
    stream << "[ ";
    for(auto element : v)
        stream << element << ' ';
    stream << "]\n";
}

void printVecVector(std::vector<std::vector<int>> w, std::ostream& stream) {
    stream << "[" << std::endl;
    for(auto element: w) {
        printVector(element, stream);
    }
    stream << "]" << std::endl;
}

void printMatrix(Eigen::MatrixXi &M, std::ostream& stream) {
    if(M.cols()>1) stream << M.format(CSVFormat) << std::endl;
    else stream << M << std::endl;
}

void printMatVec(std::vector<Eigen::MatrixXi> &m, std::ostream& stream) {
    for(Eigen::MatrixXi element: m) {
        printMatrix(element, stream);
        stream << std::endl;
    }
}

void printMySet(std::unordered_set<std::vector<int>, VectorHash> &unorderedsetOfVectors, std::ostream& stream) {
    for (auto it = unorderedsetOfVectors.begin(); it != unorderedsetOfVectors.end(); it++) {
        printVector(*it, stream);
    }
}




// From https://stackoverflow.com/a/31169617
std::vector<std::vector<int>> cartesian( std::vector<std::vector<int> >& v ) {
    std::vector<std::vector<int>> w;
    auto product = []( long long a, std::vector<int>& b ) { return a*b.size(); };
    const long long N = accumulate( v.begin(), v.end(), 1LL, product );
    std::vector<int> u(v.size());
    for( long long n=0 ; n<N ; ++n ) {
        lldiv_t q { n, 0 };
        // div_t q {n, 0};
        for( long long i=v.size()-1 ; 0<=i ; --i ) {
            q = lldiv( q.quot, v[i].size() );
            u[i] = v[i][q.rem];
        }
        w.push_back(u);
    }
    return w;
}



std::vector<int> range(int n) {
    std::vector<int> v(n);
    std::iota(std::begin(v), std::end(v), 0); //0 is the starting number
    return v;
}




// Puts all "subsets" of given size `left` of the vector arr into the vector v
// from stackoverflow at some point
void subsetsHelp(std::vector<int> arr, int left, int index, std::vector<int> &l, std::vector<std::vector<int>> &v)  {
    if(left==0){
        v.push_back(l);
        return;
    }
    for(unsigned int i=index; i< arr.size();i++){
        l.push_back(arr[i]);
        subsetsHelp(arr,left-1,i+1,l, v);
        l.pop_back();
    }
}    


std::vector<std::vector<int>> subsets(std::vector<int> v, int k) {
    std::vector<int> lt;
    std::vector<std::vector<int>> w;
    subsetsHelp(v, k, 0, lt, w);
    return w;
}




size_t VectorHash::operator()(const std::vector<int>& v) const {
    std::hash<int> hasher;
    size_t seed = 0;
    for (int i : v) {
        seed ^= hasher(i) + 0x9e3779b9 + (seed<<6) + (seed>>2);
    }
    return seed;
}


void sortMatrix(Eigen::MatrixXi& A) {
    int m = A.cols();
    std::vector<int> v = range(m);
    for(int i = 0; i < m; i++) {
        int j = 0;
        while(j < A.rows() && A(j, i) == 0) j++; 
        if(j == A.rows()) continue; //if the column of A is zero
        if(A(j, i) > 0) {
            continue;
        }
        else if (A(j, i) < 0) {
            A.col(i) = -A.col(i);
            continue;
        }
    }
    std::sort(v.begin(), v.end(),
        [&](int x, int y) -> bool {
            int j = 0;
            while(j < A.rows() && A(j, x) == A(j, y)) j++;
            if(j==A.rows()) return true; 
            return A(j, x) < A(j, y);
        });
    Eigen::MatrixXi tmp = A(Eigen::indexing::all, v);
    A = tmp;
}



int gcd(int a, int b) {
    while (b != 0)  {
        int t = b;
        b = a % b;
        a = t;
    }
    return a;
}



Eigen::MatrixXi matrixFromColumns(std::vector<Eigen::Matrix<int,-1,1>> &L) {
    Eigen::MatrixXi D = L[0];
    D.conservativeResize(D.rows(), L.size());
    for(unsigned int j = 1; j < L.size(); j++) {
        D.col(j) = L[j].col(0);
    }
    return D;
}



// v[w]
std::vector<int> vectorIndexedByAnother(std::vector<int> v, std::vector<int> w) {
    std::vector<int> u;
    for(int i: w) {
        u.push_back(v[i]);
    }
    return u;
}




Eigen::MatrixXi trivialConstruction(int Delta, int r) {
    Eigen::MatrixXi A = Eigen::MatrixXi::Zero(r, r+1);
    for(int i=0; i<r; i++) {
        A(i, i) = 1;
        A(i, r) = Delta;
    }
    return A;
}