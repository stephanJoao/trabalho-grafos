#include <iostream>
#include <vector>

#include "include/Graph.hpp"

using namespace std;

int main() {
    cout << "Hello World!!!" << endl;
    
    Graph *g = new Graph(0, true, true, false);
    g->insertEdge(0, 1, 7);
    g->insertEdge(0, 2, 1);
    g->insertEdge(1, 3, 4);
    g->insertEdge(1, 5, 1);
    g->insertEdge(2, 1, 5);
    g->insertEdge(2, 4, 2);
    g->insertEdge(2, 5, 7);
    g->insertEdge(4, 3, 5);
    g->insertEdge(4, 1, 2);
    g->insertEdge(5, 4, 3);

    g->saveToDot("graph1.dot");
    g->Dijkstra(0, 1);
    //g->printAdjList();
    //g->BFS(1);

    delete g;

    return 0;
}