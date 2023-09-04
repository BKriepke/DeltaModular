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
    int r; //largest hyperedges have size <=r
    std::vector<MySet> edges;  //stores the hyperedges by size, edges[k] stores edges of size k
    // the idea is that during the building of the hypergraph this allows to loop over edges of size s-1 to find all edges of size s
    // using that the hypergraph will be downward-closed

    // Constructor
    HyperGraph(int n, int r);

    void addEdge(std::vector<int> v);

    bool isEdge(std::vector<int> v);

    // for hyperedge Q, compute candidates C, so that Q union {c} is clique for each c in C
    std::vector<int> InitialCandidates(std::vector<int> Q);

    // hclique main call, computes a single maximum clique in hypergraph
    std::vector<int> hCliqueMain(int maxAlreadyFound);

    // computes all maximum cliques in hypergraph. Assumes that size of maximum cliques is already known
    std::vector<std::vector<int>> hCliqueMainAll(int targetSize);

    // finds candidates c, so that Q union {a} union {c} is still clique
    // those are exactly the nodes which are connected to a in graph G, where a = C[i]
    std::vector<int> NewCandidates(std::vector<int> C, Graph G, int i);


    // given clique Q, checks whether it is bigger than the currently largest found clique Qmax
    // Then recursively calls function again with all possible extensions of Q
    void hCliqueExpand(std::vector<int> Q, std::vector<int> C, std::vector<int>& Qmax, int maxAlreadyFound);

    // targetSize is the size of maximum cliques in hypergraph
    // checks if clique Q is of that size and if yes, adds it to vector of all currently found maximum cliques
    // otherwise recursively calls function again with all possible extensions of Q
    void hCliqueExpandAll(std::vector<int> Q, std::vector<int> C, std::vector<std::vector<int>>& Qmax, int targetSize);

    // if a,b are in C, then for Q union {a,b} to be a clique, you need that {q_1, ..., q_r-2, a, b} is a hyperedge
    // for all choices q_1, ..., q_r-2 in Q. Compute a graph G with edges ab iff this condition is true
    // if we want to add candidates c_1, ... c_m to Q, then they need to form a clique in this graph G
    // therefore the clique number of G gives an upper bound on the size of the maximum clique in Q union C
    // this allows earlier pruning of branches of the search tree
    // we later approximate the clique number by a greedy coloring approach
    Graph buildGraph(std::vector<int> Q, std::vector<int> C);

    // prints number of nodes and number of edges by size
    void print(std::ostream& stream=std::cout);
};

#endif