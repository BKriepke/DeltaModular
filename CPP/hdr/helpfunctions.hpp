#ifndef HELPFUNCTIONS_H
#define HELPFUNCTIONS_H

#include <vector>
#include <unordered_set>
#include <eigen34/Eigen/Eigen>
#include <iostream>



struct VectorHash {
    size_t operator()(const std::vector<int>& v) const;
};



// prints

void printVector(std::vector<int> v, std::ostream& stream=std::cout);

void printVecVector(std::vector<std::vector<int>> w, std::ostream& stream=std::cout);

void printMatrix(Eigen::MatrixXi &M, std::ostream& stream=std::cout);

void printMatVec(std::vector<Eigen::MatrixXi> &m, std::ostream& stream=std::cout);

void printMySet(std::unordered_set<std::vector<int>, VectorHash> &unorderedsetOfVectors, std::ostream& stream=std::cout);





std::vector<std::vector<int>> cartesian( std::vector<std::vector<int> >& v);



// returns {0, 1, 2, ..., n-1} 
std::vector<int> range(int n);



// Puts all "subsets" of given size `left` of the vector arr into the vector v
// Helper function for subsets(...)
void subsetsHelp(std::vector<int> arr, int left, int index, std::vector<int> &l, std::vector<std::vector<int>> &v);

// gives all subsets of size k of vector v
std::vector<std::vector<int>> subsets(std::vector<int> v, int k);



void sortMatrix(Eigen::MatrixXi& A);



int gcd(int a, int b);




Eigen::MatrixXi matrixFromColumns(std::vector<Eigen::Matrix<int,-1,1>> &L);




std::vector<int> vectorIndexedByAnother(std::vector<int> v, std::vector<int> w);




Eigen::MatrixXi trivialConstruction(int Delta, int r);


#endif
