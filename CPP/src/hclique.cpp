#include <vector>
#include <algorithm>
#include <unordered_set>

#include "helpfunctions.hpp"
#include "coloring.hpp"

#include "hclique.hpp"

using MySet = std::unordered_set<std::vector<int>, VectorHash>;




HyperGraph::HyperGraph(int n, int r) {
    this->n = n;
    this->r = r;
    MySet v;
    for(int i = 0; i <= r; i++) {
        edges.push_back(v);
    }
    for(int i = 0; i < n; i++) {
        edges[1].insert({i});   //just the nodes itself
    }
}

void HyperGraph::addEdge(std::vector<int> v) {
    int s = v.size();
    edges[s].insert(v);
    return;
}

bool HyperGraph::isEdge(std::vector<int> v) {
    int s = v.size();
    return edges[s].count(v) > 0;
}

std::vector<int> HyperGraph::InitialCandidates(std::vector<int> Q) {
    std::vector<int> C;
    int m = Q.back();  //max index in Q

    std::vector<std::vector<int>> v;
    if (r==3) v = {{Q[0], Q[1]}, {Q[0], Q[2]}, {Q[1], Q[2]}};
    else {
        v = subsets(Q, r-1);
    }
    for(auto & subedge: v) subedge.push_back(-1);
    
    for(int i = m+1; i < n; i++) {
        bool isValid = true;
        for(std::vector<int> & subedge: v) {
            subedge[r-1] = i;
            if(!isEdge(subedge)) {
                isValid = false;
                break;
            }
        }
        if(isValid) {
            C.push_back(i);
        }
    }
    return C;
}

std::vector<int> HyperGraph::hCliqueMain(int maxAlreadyFound) {
    std::vector<int> Qmax;
    int k = r;
    // hypergraph might not contain edges of size r, find biggest edge size in H
    while(edges[k].size() == 0) k--;  // edges[1].size() > 0 so will stop
    if(k!=r) {
        // just return one edge of size k (if k==1 then returns a node)
        return *edges[k].begin();
    }
    for(auto it = edges[r].begin(); it != edges[r].end(); it++) {
        std::vector<int> Q = *it;
        std::vector<int> C = InitialCandidates(Q);
        hCliqueExpand(Q, C, Qmax, maxAlreadyFound);
    }
    return Qmax;
}

std::vector<int> HyperGraph::NewCandidates(std::vector<int> C, Graph G, int i) {
    int n = C.size();
    std::vector<int> Cnew;
    for(int j = i+1; j < n; j++) {
        int b = C[j];
        if(G.isEdge(i, j)) {
            Cnew.push_back(b);
        }
    }
    return Cnew;
}

void HyperGraph::hCliqueExpand(std::vector<int> Q, std::vector<int> C, std::vector<int>& Qmax, int maxAlreadyFound) {
    if(Q.size() > Qmax.size()) {
        Qmax = Q;
    }
    int n = C.size();
    int qs = Qmax.size();
    int m = std::max(qs, maxAlreadyFound);
    int q = Q.size();
    if(q + n <= m) {
        return;
    }
    Graph G = buildGraph(Q, C);
    int chi = G.greedyColoringNumber();
    if(q + chi <= m) {
        return;
    }
    for(int i = 0; i < n; i++) {
        std::vector<int> Cnew = NewCandidates(C, G, i);
        int a = C[i];
        Q.push_back(a);
        hCliqueExpand(Q, Cnew, Qmax, m);
        Q.pop_back();
    }
}

std::vector<std::vector<int>> HyperGraph::hCliqueMainAll(int targetSize) {
    std::vector<std::vector<int>> Qmax;
    int k = r;
    // hypergraph might not contain edges of size r, find biggest edge size in H
    while(edges[k].size() == 0) k--;  // edges[1].size() > 0 so will stop
    if(k!=r) {
        if(k < targetSize) return {};
        // else k == targetSize, k>targetSize not possible as targetSize was previously computed to be max possible
        // return all edges of size k
        Qmax.insert(Qmax.end(), edges[k].begin(), edges[k].end());
        return Qmax;
    }
    for(auto it = edges[r].begin(); it != edges[r].end(); it++) {
        std::vector<int> Q = *it;
        std::vector<int> C = InitialCandidates(Q);
        hCliqueExpandAll(Q, C, Qmax, targetSize);
    }
    return Qmax;
}

void HyperGraph::hCliqueExpandAll(std::vector<int> Q, std::vector<int> C, std::vector<std::vector<int>>& Qmax, int targetSize) {
    int q = Q.size();
    if(q == targetSize) {
        Qmax.push_back(Q);
        return;
    }
    int n = C.size();
    if(q + n < targetSize) {
        return;
    }
    Graph G = buildGraph(Q, C);
    int chi = G.greedyColoringNumber();
    if(q + chi < targetSize) {
        return;
    }
    for(int i = 0; i < n; i++) {
        std::vector<int> Cnew = NewCandidates(C, G, i);
        int a = C[i];
        Q.push_back(a);
        hCliqueExpandAll(Q, Cnew, Qmax, targetSize);
        Q.pop_back();
    }
}

Graph HyperGraph::buildGraph(std::vector<int> Q, std::vector<int> C) {
    int n = C.size();
    std::vector<int> v;
    Graph G(n);
    // could also use downward-closed property here and go through all subedges of Q of size <= r-2
    // not sure how much more efficient that would be
    std::vector<std::vector<int>> subsetsQ = subsets(Q, r-2);
    for(auto & subedge: subsetsQ) {
        subedge.push_back(-1);
        subedge.push_back(-1);
    }

    v.push_back(-1);
    v.push_back(-1);
    for(int i = 0; i < n; i++) {
        int a = C[i];
        v[0] = a;
        for(auto & subedge: subsetsQ) subedge[r-2] = a;
        for(int j = i+1; j < n; j++) {
            int b = C[j];
            v[1] = b;

            if(isEdge(v)) {
                bool isValid = true;

                for(auto & subedge: subsetsQ) {
                    subedge[r-1] = b;
                    if(!isEdge(subedge)) {
                        isValid = false;
                        break;
                    }
                }

                if(isValid) {
                    G.addEdge(i, j);
                }
            }
        }
    }

    return G;
}



void HyperGraph::print(std::ostream& stream) {
    stream << "Number of nodes: " << n << std::endl;
    for(int i = 2; i <= r; i++) {
        stream << "Number of edges of size " << i << ": " << edges[i].size() << std::endl;
        // printMySet(edges[i], stream);
        stream << std::endl;
    }
}