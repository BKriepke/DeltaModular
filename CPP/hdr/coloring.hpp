#ifndef COLORING_H
#define COLORING_H

#include <vector>
 
// A class that represents an undirected graph
class Graph
{
    int n;    // No. of vertices
    std::vector<std::vector<int>> adj;      // adjacency matrix
public:
    // Constructor and destructor
    Graph(int n);

    // function to add an edge to graph
    void addEdge(int v, int w);

    // checks if there is an edge between two nodes
    bool isEdge(int v, int w);
 
    // Assigns colors (starting from 0) to all vertices and returns
    // the number of colors used
    int greedyColoringNumber();
};

#endif