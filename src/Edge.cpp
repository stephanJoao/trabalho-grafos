#include "../include/Edge.hpp"
#include "../include/Vertex.hpp"

//* Constructors and destructors implementations

Edge::Edge(Vertex* target_vertex, Edge* next_edge, float weight)
{
    this->target_vertex = target_vertex;
    this->next_edge = next_edge;
    this->weight = weight;
}

Edge::~Edge()
{
    if (this->next_edge != nullptr){
        this->next_edge = nullptr;
    }
}

//* Getters and setters implementations

Vertex* Edge::getTargetVertex()
{
    return this->target_vertex;
}

Edge* Edge::getNextEdge()
{
    return this->next_edge;
}

void Edge::setNextEdge(Edge* edge)
{
    this->next_edge = edge;
}

float Edge::getWeight()
{
    return this->weight;
}

void Edge::setWeight(float weight)
{
    this->weight = weight;
}

//* Other methods implementations