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

//Function that verifies if the graph is directed
bool Graph::getDirected()
{
    return this->directed;
}

//Function that verifies if the graph is weighted at the edges
bool Graph::getWeightedEdge()
{
    return this->weighted_edge;
}

//Function that verifies if the graph is weighted at the vertexs
bool Graph::getWeightedVertex()
{
    return this->weighted_vertex;
}

Vertex *Graph::getFirstVertex()
{
    return this->first_vertex;
}

Vertex *Graph::getLastVertex()
{
    return this->last_vertex;
}

// // Other methods
// /*
//     The outdegree attribute of vertexs is used as a counter for the number of edges in the graph.
//     This allows the correct updating of the numbers of edges in the graph being directed or not.
// */
// void Graph::insertVertex(int id)
// {
    
// }

// void Graph::insertEdge(int id, int target_id, float weight)
// {

    
// }

// void Graph::removeVertex(int id){ 
    
// }

// bool Graph::searchVertex(int id)
// {
    
// }

// Vertex *Graph::getVertex(int id)
// {

    
// }


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
