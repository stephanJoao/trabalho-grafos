#ifndef EDGE_HPP
#define EDGE_HPP

#include "Vertex.hpp"

class Edge
{
    // Attributes
    private:
        Vertex* target_node;
        float weight;
        Edge* next_edge;

    // Methods
    private:

    public:
        // Constructors and destructors
        Edge(Vertex* target_node, Edge* next_edge, float weight);
        ~Edge();

        // Getters and setters
        Vertex* getTargetNode();
        Edge* getNextEdge();
        void setNextEdge(Edge* edge);
        float getWeight();
        void setWeight(float weight);

        // Functions
};

#endif 