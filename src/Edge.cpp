#include "../include/Edge.hpp"

//* Constructors and destructors implementations

Edge::Edge(int target_id, float weight)
{
    this->target_id = target_id;
    this->weight = weight;
}

Edge::~Edge()
{
    
}

//* Getters and setters implementations

int Edge::getTargetId()
{
    return this->target_id;
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