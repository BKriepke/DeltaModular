#include <eigen34/Eigen/Eigen>
#include <vector>

#include "hclique.hpp"
#include "helpfunctions.hpp"

#include "hypergraphBuild.hpp"



#define contains(x) count(x)



float det(Eigen::MatrixXf &B) {
	int r = B.cols();
	if(r==2) {
		return B(0,0)*B(1,1)-B(0,1)*B(1,0);
	}
	int s = 0;
	if(r==3) {
		// Sarrus
		s += B(0,0)*B(1,1)*B(2,2);
		s += B(0,1)*B(1,2)*B(2,0);
		s += B(0,2)*B(1,0)*B(2,1);
		s -= B(0,0)*B(1,2)*B(2,1);
		s -= B(0,1)*B(1,0)*B(2,2);
		s -= B(0,2)*B(1,1)*B(2,0);
		return s;
	}
	// if r>3 then use determinant from Eigen
	return B.determinant();
}

bool checkIndRowsCols(Eigen::MatrixXf &M, int Delta, std::vector<int> &rows, std::vector<int> &cols, bool generic) {
	int k = rows.size();
	int n = cols.size();
	assert(k == n);
	Eigen::MatrixXf B = M(rows, cols);
	int d = round(abs(det(B)));
	if(generic) return d != 0 && d <= std::pow(Delta, k);
	else return d <= std::pow(Delta, k);
}

bool checkCols(Eigen::MatrixXf& M, int Delta, std::vector<int> &cols, bool generic) {
	//Assumes that all smaller combinations of the columns have been checked
	
	// Now checks the options for the rows
	int r = M.rows();
	int k = cols.size();
	assert(k <= r);

	std::vector<std::vector<int>> v;

	if(r==3) {
		if(k==2) v = {{0, 1}, {0, 2}, {1, 2}};
		else v = {{0,1,2}};
	}
	else {
		std::vector<int> rows = range(r);
		v = subsets(rows, k);
	}
	for(auto it = v.begin(); it != v.end(); it++) {
		if (!checkIndRowsCols(M, Delta, *it, cols, generic)) {
			return false;
		}
	}
	return true;
}

HyperGraph buildHyperGraph(Eigen::MatrixXf& M, int Delta, bool generic) {
    int r = M.rows();   
    int n = M.cols();

    HyperGraph H(n, r);

    // add all edges, starting from size 1 (aka nodes)
	// uses that the graph is downward-closed
	for(int i = 1; i < r; i++) {
		//current edges have size i
		std::vector<std::vector<int>> subedgeIndices = subsets(range(i), i-1);
		for(std::vector<int> e: H.edges[i]) {

			int maxInEdge = e.back();   // e is sorted

			for(int c = maxInEdge+1; c < n; c++) {
				bool isValid = true;

				for(auto indices: subedgeIndices) {
					std::vector<int> subedge = vectorIndexedByAnother(e, indices);
					subedge.push_back(c);
					if(!H.isEdge(subedge)) {
						isValid = false;
						break;
					}
				}
				if(isValid) {
					e.push_back(c);
					if(checkCols(M, Delta, e, generic)) {
						H.addEdge(e);
					}
					e.pop_back();
				}
			}
		}
	}
    
    return H;
}