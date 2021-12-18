#include <iostream>
#include <vector>

#include "include/Graph.hpp"

using namespace std;

int main() {
    cout << "Hello World!!!" << endl;
    
    Graph *g = new Graph(0, true, true, false);
    g->insertEdge(0, 1, 4);
    g->insertEdge(0, 2, 2);
    g->insertEdge(1, 4, 3);
    g->insertEdge(1, 3, 1);
    g->insertEdge(1, 2, -3);
    g->insertEdge(2, 4, 2);
    g->insertEdge(2, 3, 3);
    g->insertEdge(3, 4, -2);
    g->insertEdge(3, 6, 4);
    g->insertEdge(4, 5, 3);
    g->insertEdge(4, 6, 3);
    g->insertEdge(5, 6, 1);


    g->saveToDot("graph1.dot");
    g->Dijkstra(0, 1);
    //g->printAdjList();
    //g->BFS(1);

    delete g;

    return 0;
}