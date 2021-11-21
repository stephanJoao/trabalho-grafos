#include <iostream>

#include "include/Vertex.hpp"
#include "include/Edge.hpp"

using namespace std;

int main() {
    Vertex *v1 = new Vertex(1);
    Vertex *v2 = new Vertex(2);
    v1->insertEdge(v2);
    delete v2;
    delete v1;
    return 0;
}