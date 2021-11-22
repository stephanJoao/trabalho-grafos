#include <iostream>
// #include <fstream>
// #include <stack>
// #include <queue>
// #include <list>
// #include <math.h>
// #include <cstdlib>
// #include <ctime>
// #include <float.h>
// #include <iomanip>

#include "../include/Graph.hpp"
#include "../include/Vertex.hpp"
#include "../include/Edge.hpp"

//* Constructors and destructors implementations

Graph::Graph(int order, bool directed, bool weighted_edge, bool weighted_vertex)
{
    this->order = order;
    this->directed = directed;
    this->weighted_edge = weighted_edge;
    this->weighted_vertex = weighted_vertex;
    this->first_vertex = nullptr;
    this->last_vertex = nullptr;
    this->number_edges = 0;
}

Graph::~Graph()
{
    Vertex *next_vertex = this->first_vertex;

    while (next_vertex != nullptr)
    {
        next_vertex->removeAllEdges();
        Vertex *aux_vertex = next_vertex->getNextVertex();
        delete next_vertex;
        next_vertex = aux_vertex;
    }
}

//* Getters and setters implementations

int Graph::getOrder()
{
    return this->order;
}

int Graph::getNumberEdges()
{
    return this->number_edges;
}

// Function that verifies if the graph is directed
bool Graph::getDirected()
{
    return this->directed;
}

// Function that verifies if the graph is weighted at the edges
bool Graph::getWeightedEdge()
{
    return this->weighted_edge;
}

// Function that verifies if the graph is weighted at the vertexs
bool Graph::getWeightedVertex()
{
    return this->weighted_vertex;
}

Vertex* Graph::getFirstVertex()
{
    return this->first_vertex;
}

Vertex* Graph::getLastVertex()
{
    return this->last_vertex;
}

Vertex* Graph::getVertex(int id)
{
    Vertex* v = this->first_vertex;  
    while(v != nullptr && v->getId() != id) {
        v = v->getNextVertex();
    }
    if(v == nullptr) {
        std::cerr << "Vertex not found!" << std::endl;
        return nullptr;
    }
    return v;
}

//* Other methods

/*
    The outdegree attribute of vertices is used as a counter for the number of edges in the graph.
    This allows the correct updating of the numbers of edges in the graph being directed or not.
 */

/**
 * Insert vertex in graph.
 *
 * @param id Edge id
 * @param directed If the graph is directed
 * @param target_vertex Target vertex
 */
void Graph::insertVertex(int id, float weight)
{
    Vertex *v = new Vertex(id, weight);
    if (this->first_vertex == nullptr)
        this->first_vertex = v;
    else
        this->last_vertex->setNextVertex(v);
    
    this->last_vertex = v;
}

/**
 * @brief Insert edge on the vertex edges list
 * 
 * @param id First vertex's ID
 * @param target_id Second Vertex's ID
 * @param weight Weight of the edge
 */
void Graph::insertEdge(int id, int target_id, float weight)
{
    Vertex *v = this->getVertex(id);
    if (this->directed){
        v->insertEdge(target_id, weight);
    } else {
        v->insertEdge(target_id, weight);
        v = this->getVertex(target_id);
        v->insertEdge(id, weight);
    }
}

// void Graph::removeVertex(int id)
// { 
//     Vertex *v = this->getVertex(id);
//     if (v != nullptr)
//         delete v;
//     else 
//         return;     
// }

/**
 * @brief Return if the vertex is in the graph
 * 
 * @param id Vertex's ID
 * @return true if
 * @return false 
 */
bool Graph::searchVertex(int id)
{
    if (this->first_vertex != nullptr)
    {
        for (Vertex *v = this->first_vertex; v != nullptr; v = v->getNextVertex())
            if (v->getId() == id)
                return true;
    }
    
    return false;
}

void Graph::printAdjList() 
{
    std::cout << "Imprimindo Grafo como uma lista de adjacência:\n" << std::endl;
    Vertex* v = first_vertex;
    while(v != nullptr) {
        std::cout << v->getId();
        Edge* e = v->getFirstEdge();
        while(e != nullptr) {
            std::cout << " -> " << e->getTargetId();
            e = e->getNextEdge();
        }
        std::cout << std::endl;
        v = v->getNextVertex();
    }
}

// //Function that prints a set of edges belongs breadth tree

// void Graph::breadthFirstSearch(std::ofstream &output_file){
    
// }



// float Graph::floydMarshall(int idSource, int idTarget){
    
// }

   

// float Graph::dijkstra(int idSource, int idTarget){
    
// }

// //function that prints a topological sorting
// void topologicalSorting(){

// }

// void breadthFirstSearch(std::ofstream& output_file){

// }

// Graph* getVertexInduced(int* listIdVertexs){

// }

// Graph* agmKuskal(){

// }
// Graph* agmPrim(){

// }
