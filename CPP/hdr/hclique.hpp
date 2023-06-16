#ifndef HCLIQUE_H
#define HCLIQUE_H

#include <vector>
#include <unordered_set>
#include <iostream>

#include "coloring.hpp"
#include "helpfunctions.hpp"

using MySet = std::unordered_set<std::vector<int>, VectorHash>;

class HyperGraph {
public:
    int n;
    int r; //largest hyperedges have size r
    std::vector<MySet> edges;  //stores the hyperedges by size, edges[k] stores edges of size k
    // the idea is that during the building of the hypergraph this allows to loop over edges of size s-1 to find all edges of size s
    // using that the hypergraph will be downward-closed

    // Constructor and destructor
    HyperGraph(int n, int r);

    void addEdge(std::vector<int> v);

    bool isEdge(std::vector<int> v);

    std::vector<int> InitialCandidates(std::vector<int> Q);

    std::vector<int> hCliqueMain(int maxAlreadyFound);

    std::vector<std::vector<int>> hCliqueMainAll(int targetSize);

    std::vector<int> NewCandidates(std::vector<int> C, Graph G, int i);

    void hCliqueExpand(std::vector<int> Q, std::vector<int> C, std::vector<int>& Qmax, int maxAlreadyFound);

    void hCliqueExpandAll(std::vector<int> Q, std::vector<int> C, std::vector<std::vector<int>>& Qmax, int targetSize);

    Graph buildGraph(std::vector<int> Q, std::vector<int> C);

    void print(std::ostream& stream=std::cout);
};

#endif