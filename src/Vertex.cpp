#include <unordered_map>
#include <iostream>
#include <fstream>

#include "../include/Vertex.hpp"

//* Constructors and destructors implementations

Vertex::Vertex(int id, float weight)
{
    this->id = id;
    this->out_degree = 0;
    this->in_degree = 0;
    this->weight = weight;
}

Vertex::~Vertex()
{
    for (std::unordered_map<int, Edge*>::iterator it = edges.begin(); it != edges.end(); ++it)
        delete it->second;
}

//* Getters and setters implementations

int Vertex::getId()
{
    return this->id;
}

int Vertex::getInDegree()
{
    return this->in_degree;
}

int Vertex::getOutDegree()
{
    return this->out_degree;
}

float Vertex::getWeight()
{
    return this->weight;
}

void Vertex::setWeight(float weight)
{
    this->weight = weight;
}

std::unordered_map<int, Edge*> Vertex::getEdges() 
{
    return this->edges;
}

//* Other methods implementations

/**
 * Search edge in vertex.
 *
 * @param target Target vertex id
 */
bool Vertex::searchEdge(int target_id)
{
    return edges[target_id] != nullptr;
}

/**
 * Insert edge in vertex.
 *
 * @param target Target vertex id
 * @param weight Weight of the edge (defaults to 1 if not given)
 */
void Vertex::insertEdge(int target_id, float weight) 
{
    if(edges.count(target_id) == 0) {
        Edge* e = new Edge(target_id, weight);
        edges.insert({target_id, e});
    }
    else
        std::cout << "WARNING: " << id << " -> " << target_id << " already exists" << std::endl;
}

/**
 * Remove all edges from vertex.
 */
void Vertex::removeAllEdges()
{
    this->edges.clear();
}

/**
 * Remove edge from vertex.
 *
 * @param id Edge id
 * @param directed If the graph is directed
 * @param target_vertex Target vertex
 */
// int Vertex::removeEdge(int id, bool directed, Vertex* target_vertex)
// {
//     this->edges.erase(id);
// }

void Vertex::incrementInDegree()
{
    this->in_degree++;
}

void Vertex::incrementOutDegree()
{
    this->out_degree++;
}

void Vertex::decrementInDegree()
{
    this->in_degree--;
}

void Vertex::decrementOutDegree()
{
    this->out_degree--;
}

// Edge* Vertex::hasEdgeBetween(int target_id)
// {

//     for(Edge *auxEdge = this->first_edge; auxEdge != nullptr; auxEdge = auxEdge->getNextEdge())
//     {
//         if(auxEdge->getTargetVertex()->getId() == target_id)
//             return auxEdge;
//     }
//     return nullptr;
// }
