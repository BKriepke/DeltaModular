#ifndef MODULO_H
#define MODULO_H

#include <vector>
#include <numeric>
#include <eigen34/Eigen/Eigen>





int gcd(int a, int b);

int gcdVec(Eigen::Matrix<int, -1, 1> &M);

int axbmodnOneSol(int a, int b, int n);

std::vector<int> axbmodn(int a, int b, int n);

bool checkIfRedundant(std::vector<Eigen::Matrix<int, -1, 1>> &L, Eigen::Matrix<int, -1, 1> &M, int Delta, bool generic);

void allRepresentants(std::vector<int> v, std::vector<Eigen::Matrix<int, -1, 1>> &L, int Delta, bool generic);

std::vector<Eigen::Matrix<int,-1,1>> Ax0modn(Eigen::MatrixXi &A, int n, bool generic);

#endif