#include "../include/Vertex.hpp"

// Constructors and destructors implementations
Vertex::Vertex(int id)
{
    this->id = id;
    this->first_edge = nullptr;
    this->out_degree = 0;
    this->in_degree = 0; 
}

Vertex::~Vertex()
{
    Edge* next_edge = this->first_edge;

    while(next_edge != nullptr)
    {
        Edge* aux_edge = next_edge->getNextEdge();
        delete next_edge;
        next_edge = aux_edge;
    }
}

// Getters and setters implementations
Edge* Vertex::getFirstEdge()
{
    return this->first_edge;
}

int Vertex::getId()
{
    return id;
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
    return weight;
}

void Vertex::setWeight(float weight)
{
    this->weight = weight;
}

// Functions implementations

/**
 * Insert edge in vertex.
 *
 * @param target Pointer to target vertex 
 * @param weight Weight of the edge (if not received ...)
 */
void Vertex::insertEdge(Vertex* target, float weight = 1) 
{
    // checks if there is at least one edge in vertex
    if (this->first_edge != nullptr){
        //
        Edge* new_edge = new Edge(target, nullptr, weight);
        this->first_edge = new_edge;
        this->last_edge = new_edge;
    } else {
        //
        
    }
}