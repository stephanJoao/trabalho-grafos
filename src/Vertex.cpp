#include "../include/Vertex.hpp"

Vertex::Vertex()
{
    //ctor
}

Vertex::~Vertex()
{
    //dtor
}

Vertex::Vertex(int value)
{
    this->value = value;
}

int Vertex::getValue()
{
    return value;
}