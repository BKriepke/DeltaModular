// A C++ program to implement greedy algorithm for graph coloring
// From https://www.geeksforgeeks.org/graph-coloring-set-2-greedy-algorithm/
#include <vector>
#include <algorithm>

#include "helpfunctions.hpp"

#include "coloring.hpp"
 

Graph::Graph(int n) {
    this->n = n;
    adj = std::vector<std::vector<int>>(n, std::vector<int>(0));
}

void Graph::addEdge(int v, int w) {
    adj[v].push_back(w);
    adj[w].push_back(v);  // Note: the graph is undirected
}

bool Graph::isEdge(int v, int w) {
    return std::find(adj[v].begin(), adj[v].end(), w) != adj[v].end();
}
 
// Assigns colors (starting from 0) to all vertices and returns
// the number of colors used
int Graph::greedyColoringNumber()
{
    std::vector<int> result(n);

    std::vector<int> reorder = range(n);

    std::sort(reorder.begin(), reorder.end(), 
                [&](int a, int b) -> bool {
                    return adj[a].size() > adj[b].size();
                });
 
    // Assign the first color to first vertex
    result[0]  = 0;
 
    // Initialize remaining V-1 vertices as unassigned
    for (int u = 1; u < n; u++)
        result[u] = -1;  // no color is assigned to u
 
    // A temporary array to store the available colors. True
    // value of available[cr] would mean that the color cr is
    // assigned to one of its adjacent vertices
    std::vector<bool> available(n);
    for (int cr = 0; cr < n; cr++)
        available[cr] = false;
 
    // Assign colors to remaining V-1 vertices
    // for (int u = 1; u < n; u++)
    for(int j = 0; j < n; j++)
    {
        int u = reorder[j];
        // Process all adjacent vertices and flag their colors
        // as unavailable
        std::vector<int>::iterator i;
        for (i = adj[u].begin(); i != adj[u].end(); ++i)
            if (result[*i] != -1)
                available[result[*i]] = true;
 
        // Find the first available color
        int cr;
        for (cr = 0; cr < n; cr++)
            if (available[cr] == false)
                break;
 
        result[u] = cr; // Assign the found color
 
        // Reset the values back to false for the next iteration
        for (i = adj[u].begin(); i != adj[u].end(); ++i)
            if (result[*i] != -1)
                available[result[*i]] = false;
    }
 
    return (*std::max_element(result.begin() , result.end()))+1;
}