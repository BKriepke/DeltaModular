#ifndef HYPERGRAPHALG_H
#define HYPERGRAPHALG_H

#include <eigen34/Eigen/Eigen>
#include <vector>

#include "hclique.hpp"


float det(Eigen::MatrixXf &B);

bool checkIndRowsCols(Eigen::MatrixXf &M, int Delta, std::vector<int> &rows, std::vector<int> &cols, bool generic);

bool checkCols(Eigen::MatrixXf& M, int Delta, std::vector<int> &cols, bool generic);

HyperGraph buildHyperGraph(Eigen::MatrixXf& M, int Delta, bool generic);

#endif