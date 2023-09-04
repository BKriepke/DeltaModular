#include <eigen34/Eigen/Eigen>
#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include <string>
#include <stdexcept>

#include <omp.h>

#include "hermite_forms.hpp"
#include "hclique.hpp"
#include "solve_modulo.hpp"
#include "hypergraphBuild.hpp"





// increases counter and prints progress update
void increaseCounter(int& counter, unsigned int s) {
    counter++;
    if(counter % 10 == 0) {
        std::cout << "HNFs checked: " << counter << " / " << s << std::endl;
    }
}

// measures time since last time this function was called. Useful for some testing
int timeSinceLast(std::chrono::_V2::system_clock::time_point &last) {
    auto now = std::chrono::system_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(now-last);
    last = now;
    return duration.count();
}

// reads g(Delta, r) or h(Delta, r) from file. Throws error if the value does not exist
int fetchMaxSize(int Delta, int r, bool generic) {
    std::string filename = "../data/";
    if(generic) filename += "Generic/";
    else filename += "NonGeneric/";
    filename += "r="+std::to_string(r)+"/values.txt";
    std::ifstream myfile (filename);
    int item;
    if(myfile.is_open()) {
        while (myfile >> item) {
            if(item == Delta) {
                myfile >> item;
                return item;
            }
            myfile >> item;
        }
        std::string errMsg;
        if(generic) errMsg += "g(";
        else errMsg += "h(";
        errMsg += std::to_string(Delta) + ","+std::to_string(r) + ") has not been computed yet";
        throw std::runtime_error(errMsg);
    }
    throw std::runtime_error(std::string("File '" + filename + "' with pre-computed values was not found"));
}

// prints the std::vector of matrices D obtained from oneDelta(...) to a file
// also prints g(Delta, r) or h(Delta, r) to file
void printToFile(int Delta, int r, bool generic, bool singleExample, std::vector<Eigen::MatrixXi> result) {
    std::fstream myfile;
    std::string filename = "../data/";
    if(generic) filename += "Generic/";
    else filename += "NonGeneric/";
    filename += "r="+std::to_string(r)+"/";

    if(singleExample) {
        myfile.open(filename+"values.txt", std::fstream::out | std::fstream::app);
        if(myfile.is_open()) {
            myfile << Delta << '\t' << result[0].cols() << std::endl;
        }
        myfile.close();
        filename += "OneExample/";
    }
    else filename += "AllExamples/";

    myfile.open(filename+"Delta="+std::to_string(Delta)+".txt", std::fstream::out);
    if(myfile.is_open()) {
        printMatVec(result, myfile);
    }
    myfile.close();
}

// either computes g(Delta, r) or h(Delta, r) and returns a (std::vector containing only one) matrix D achieving this size
// or reads the corresponding value from file and returns a list (std::vector) of matrices D achieving this size
// this list is guaranteed to contain at least one representative of every equivalence class
// the list will in general contain more than one representative per class
std::vector<Eigen::MatrixXi> oneDelta(int Delta, int r, bool generic, bool singleExample) {

    int targetSize = 0;
    if(!singleExample) targetSize = fetchMaxSize(Delta, r, generic) - r;

    std::cout << "Current Delta: " << Delta << std::endl;

    std::vector<Eigen::MatrixXi> HNFs = findHermiteForms(Delta, r);
    int s = HNFs.size();
    std::cout << "Number of HNFs: " << s << std::endl;

    unsigned int bestSoFar = 1;  // trivial construction
    Eigen::MatrixXi tmp1, tmp2;
    if(singleExample) tmp2 = trivialConstruction(Delta, r); // in case nothing better is found, so we at least return something
    else tmp2.resize(r, r+targetSize);
    int counter = 0;

    std::vector<Eigen::MatrixXi> result;

    // parallelize the for-loop. Works well in the generic setting, as most hypergraphs have roughly the same size and 
    // therefore most threads take the same time
    // does not work as well in the non-generic setting, most of the time will be spent waiting on a single thread
    // corresponding to one very big hypergraph with a lot of maximum cliques
    // might be fruitful to instead parallelize hclique itself, which is a DFS algorithm
    omp_set_num_threads(56);
    #pragma omp parallel for schedule(dynamic) 
    for(unsigned int i = 0; i < HNFs.size(); i++) {
        Eigen::MatrixXi A = HNFs[i];
        std::vector<Eigen::Matrix<int,-1,1>> L = Ax0modn(A, Delta, generic);
        Eigen::MatrixXi D = matrixFromColumns(L);
        Eigen::MatrixXf Df = D.cast<float>();           // Eigen has no determinant function for integer matrices
        HyperGraph H = buildHyperGraph(Df, Delta, generic);
        std::vector<int> v;
        std::vector<std::vector<int>> w;
        if(singleExample) {
            v = H.hCliqueMain(bestSoFar);
        }
        else {
            w = H.hCliqueMainAll(targetSize);
        }
        #pragma omp critical 
        {
            increaseCounter(counter, s);
            if(singleExample) {
                if(v.size() > bestSoFar) {
                    bestSoFar = v.size();
                    tmp1 = A*D(Eigen::indexing::all, v)/Delta;
                    tmp2.resize(r, r+bestSoFar);
                    tmp2 << A, tmp1;
                    std::cout << "Current max length found: " << bestSoFar + r << std::endl;
                    std::cout << "Current max length matrix: " << std::endl;
                    std::cout << tmp2 << std::endl;
                }
            }
            else {
                for(std::vector<int> v: w) {
                    tmp1 = A*D(Eigen::indexing::all, v)/Delta;
                    tmp2 << A, tmp1;
                    result.push_back(tmp2);
                }
            }
        }
    }
    if(singleExample) result.push_back(tmp2);

    return result;
}


int main() {
    Eigen::initParallel();

    // read in parameters from command line
    int startDelta, endDelta, r;
    bool generic, singleExample, print;
    std::cout << "Enter start Delta: ";
    std::cin >> startDelta;
    std::cout << "Enter end Delta: ";
    std::cin >> endDelta;
    std::cout << "Enter dimension r: ";
    std::cin >> r;
    std::cout << "Enter 0 for non-generic and 1 for generic: ";
    std::cin >> generic;
    std::cout << "Enter 0 for all examples and 1 for single example: ";
    std::cin >> singleExample;
    std::cout << "Do you want to print the results to a file? (1/0 for yes/no. Make sure to have to correct directory setup. This might overwrite existing files): ";
    std::cin >> print;

    for(int Delta = startDelta; Delta <= endDelta; Delta++) {
        auto start = std::chrono::system_clock::now();

        std::vector<Eigen::MatrixXi> result = oneDelta(Delta, r, generic, singleExample);
        if(print) printToFile(Delta, r, generic, singleExample, result);

        auto stop = std::chrono::system_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop-start);

        // prints results from computations to command line
        std::cout << "Delta = " << Delta << std::endl;
        std::cout << "Time taken in seconds: " << duration.count() << std::endl;
        if(singleExample) {
            std::cout << "Max length found: " << result[0].cols() << std::endl;
            std::cout << "One matrix with that length: " << std::endl;
            printMatVec(result);
        }
        else {
            std::cout << "Number of matrices found: " << result.size() << std::endl;
            printMatVec(result);
        }
        std::cout << std::endl;
    }


    return 0;
}