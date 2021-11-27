#include <iostream>

#include "include/Graph.hpp"

using namespace std;

int main() {
    cout << "Hello World!!!" << endl;
    Graph *g = new Graph(0);
    g->insertVertex(1);
    g->insertVertex(2);
    g->insertVertex(3);
    g->insertVertex(4);
    g->insertEdge(1, 2);
    g->insertEdge(1, 3);
    g->insertEdge(1, 4);
    g->insertEdge(2, 3);

    g->printAdjList();
    delete g;
    return 0;
}