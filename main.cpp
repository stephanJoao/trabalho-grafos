#include <iostream>
#include <vector>

#include "include/Graph.hpp"

using namespace std;

int main() {
    cout << "Hello World!!!" << endl;
    
    Graph *g = new Graph(0, true, true, false);
    g->insertEdge(1, 2, 7);
    g->insertEdge(1, 3, 1);
    g->insertEdge(2, 4, 4);
    g->insertEdge(2, 6, 1);
    g->insertEdge(3, 2, 5);
    g->insertEdge(3, 5, 2);
    g->insertEdge(3, 6, 7);
    g->insertEdge(5, 4, 5);
    g->insertEdge(5, 2, 2);
    g->insertEdge(6, 5, 3);

    g->saveToDot("graph1.dot");
    g->Dijkstra(1, 2);
    g->printAdjList();
    g->BFS(1);

    delete g;

    return 0;
}