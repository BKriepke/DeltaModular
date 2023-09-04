#include <stack>
#include <vector>
#include <algorithm>
#include <eigen34/Eigen/Eigen>
#include <string>
#include <fstream>

#include "helpfunctions.hpp"

#include "hermite_forms.hpp"



std::vector<Eigen::MatrixXi> hnfsFromFile(int Delta, int r) {
    std::string filename = "../data/hnfs/"+std::to_string(r)+"_"+std::to_string(Delta)+".txt";
    std::ifstream fin (filename);

    std::vector<Eigen::MatrixXi> H;
    Eigen::MatrixXi tmp = Eigen::MatrixXi::Zero(r,r);
    // https://stackoverflow.com/a/36342394
    if (fin.is_open()) {
        while(!fin.eof()) {
            for (int row = 0; row < r; row++) {
                for (int col = 0; col < r; col++) {
                    int item;
                    fin >> item;
                    tmp(row, col) = item;
                }
            }
            if(!fin.eof()) H.push_back(tmp); //check that you didnt try to read end of file
        }
        fin.close();
    }
    return H;
}


std::vector<Eigen::MatrixXi> findHermiteForms(int Delta, int r) {
    std::vector<Eigen::MatrixXi> H = hnfsFromFile(Delta, r);
    if(H.size() > 0) return H;
    // else it didnt find the file, so we manually compute the relevant hermite forms
    // will be much less efficient as we do not check for as much equivalence
    std::cout << "Warning: File with inequivalent HNFs was not found." << std::endl;
    // Compute sorted diagonals that multiply to Delta
    std::vector<std::vector<int>> diags = sortedDiags(Delta, r);
    // Successively augment with columns of the form
    // [entries < d[n], d[n], 0, ..., 0]
    for(auto d: diags) {
        Eigen::MatrixXi v = Eigen::MatrixXi::Zero(r, 1);
        v(0,0) = d[0];
        std::stack<Eigen::Matrix<int, -1, -1>> matrixStack;
        matrixStack.push(v);

        while(!matrixStack.empty()) {
            Eigen::MatrixXi M = matrixStack.top();
            matrixStack.pop();
            int n = M.cols();
            std::vector<int> entries;
            if(n+1==r && d[n] == Delta) entries = range(d[n]/2+1); //last column and diagonal is (1,1,..,Delta), applies operation 2 described in paper
            else entries = range(d[n]);
            std::vector<std::vector<int>> L(n, entries);
            std::vector<std::vector<int>> K = cartesian(L);
            for(auto I: K) {
                if(n+1==r && d[n] == Delta && !std::is_sorted(I.begin(), I.end())) continue;  // operation 2
                // J is columns [entries < d[n], d[n], 0, ..., 0]
                Eigen::MatrixXi J = Eigen::MatrixXi::Zero(r, 1);
                for(int i = 0; i < n; i++) {
                    J(i, 0) = I[i];
                }
                J(n, 0) = d[n];
                Eigen::MatrixXi res(M.rows(), n+1);
                res << M, J;
                if(n+1 == r) {
                    H.push_back(res);
                }
                else {
                    matrixStack.push(res);
                }
            }
        }
        
    }
    removeEquivalent(H);
    return H;
}


// Finds all diagonals that multiply to Delta with minimum entry min
std::vector<std::vector<int>> findDiags(int Delta, int r, int min) {
    std::vector<std::vector<int>> v;
    std::vector<int> w;
    if(r == 1) {
        v.push_back({Delta});
        return v;
    }
    for(int i = min; (i*i) <= Delta; i++) {
        if( Delta % i == 0) {
            std::vector<std::vector<int>> u = findDiags(Delta/i, r-1, i);
            for(auto element: u) {
                w.clear();
                w.push_back(i);
                w.insert(w.end(), element.begin(), element.end());
                v.push_back(w);
            }
        }
    }
    return v;
}

// Finds all diagonals that multiply to Delta with minimum entry 1
std::vector<std::vector<int>> sortedDiags(int Delta, int r) {
    return findDiags(Delta, r, 1);
}


void removeEquivalent(std::vector<Eigen::MatrixXi> & v) {

    // check gcd condition from sorting diags
    v.erase(std::remove_if(v.begin(), v.end(),
                            [](Eigen::MatrixXi const & A) {
                                for(int i=0; i < A.cols()-1; i++) {
                                    if( gcd(A(i, i+1), A(i+1, i+1)) < A(i, i)) return true;
                                }
                                return false;
                            }), v.end());

    return;
}