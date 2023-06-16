#ifndef COLORING_H
#define COLORING_H

#include <vector>
 
// A class that represents an undirected graph
class Graph
{
    int n;    // No. of vertices
    std::vector<std::vector<int>> adj;
public:
    // Constructor and destructor
    Graph(int n);

    // function to add an edge to graph
    void addEdge(int v, int w);

    // checks if there is an edge between two nodes
    bool isEdge(int v, int w);
 
    // Returns greedy coloring number of the vertices
    int greedyColoringNumber();
};

#endif