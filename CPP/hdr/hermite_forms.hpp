#ifndef HERMITEFORMS_H
#define HERMITEFORMS_H

#include <vector>
#include <eigen34/Eigen/Eigen>


std::vector<Eigen::MatrixXi> hnfsFromFile(int Delta, int r);



// finds rxr hermite forms with determinant Delta and sorted diagonal
std::vector<Eigen::MatrixXi> findHermiteForms(int Delta, int r);



// diags

std::vector<std::vector<int>> findDiags(int Delta, int r, int min);

std::vector<std::vector<int>> sortedDiags(int Delta, int r);



void removeEquivalent(std::vector<Eigen::MatrixXi> & v);

#endif