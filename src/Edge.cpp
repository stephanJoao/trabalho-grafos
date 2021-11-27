#include "../include/Edge.hpp"

//* Constructors and destructors implementations

Edge::Edge(int target_id, Edge* next_edge, float weight)
{
    this->target_id = target_id;
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

int Edge::getTargetId()
{
    return this->target_id;
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