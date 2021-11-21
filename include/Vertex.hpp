#ifndef VERTEX_HPP
#define VERTEX_HPP

#include "Edge.hpp"

class Edge;
class Vertex
{
    // Attributes
    private:
        Edge* first_edge;
        Edge* last_edge;
        int id;
        unsigned int in_degree;
        unsigned int out_degree;
        float weight;

    // Methods
    private:

    public:
        // Constructors and destructors
        Vertex(int id);
        ~Vertex();

        // Getters and setters
        Edge* getFirstEdge();
        Edge* getLastEdge();
        int getId();
        int getInDegree();
        int getOutDegree();
        float getWeight();
        void setWeight(float weight);

        // Functions
        void insertEdge(Vertex* target, float weight = 1);
};

#endif 