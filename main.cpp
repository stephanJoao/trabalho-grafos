#include <iostream>

#include "include/Graph.hpp"

using namespace std;

int main() {
    cout << "Hello World!!!" << endl;
    // Graph *g = new Graph(0);
    // g->insertEdge(1, 2);
    // g->insertEdge(1, 3);
    // g->insertEdge(1, 4);
    // g->insertEdge(2, 3);
    
    Graph *g = new Graph(0, false, true, false);
    g->insertEdge(1, 2, 5);
    g->insertEdge(1, 3, 8);
    g->insertEdge(1, 6, 7);
    g->insertEdge(1, 8, 6);
    g->insertEdge(2, 3, 4);
    g->insertEdge(2, 4, 8);
    g->insertEdge(2, 5, 7);
    g->insertEdge(2, 8, 9);
    g->insertEdge(3, 5, 4);
    g->insertEdge(3, 6, 5);
    g->insertEdge(4, 5, 9);
    g->insertEdge(4, 7, 4);
    g->insertEdge(5, 6, 3);
    g->insertEdge(5, 7, 10);
    g->insertEdge(6, 7, 6);

    g->printAdjList();
    g->BFS(1);
    g->saveToDot("graph1.dot");
    delete g;
    return 0;
}