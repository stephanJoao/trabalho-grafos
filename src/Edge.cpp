#include "../include/Edge.hpp"

// Constructors and destructors implementations
Edge::Edge(Vertex* target_node, Edge* next_edge, float weight = 1)
{
    this->target_node = target_node;
    this->next_edge = next_edge;
    this->weight = weight;
}

Edge::~Edge()
{
    if (this->next_edge != nullptr){
        delete this->next_edge;
        this->next_edge = nullptr;
    }
}

// Getters and setters implementations
Vertex* Edge::getTargetNode()
{
    return target_node;
}

Edge* Edge::getNextEdge()
{
    return next_edge;
}

void Edge::setNextEdge(Edge* edge)
{
    this->next_edge = edge;
}

float Edge::getWeight()
{
    return weight;
}

void Edge::setWeight(float weight)
{
    this->weight = weight;
}