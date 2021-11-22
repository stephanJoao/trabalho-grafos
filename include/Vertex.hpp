#ifndef VERTEX_HPP
#define VERTEX_HPP

#include "Edge.hpp"

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
        Vertex* next_vertex;

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
        Vertex* getNextVertex();
        void setNextVertex(Vertex* vertex);

        // Functions

        bool searchEdge(int target_id);
        void insertEdge(int target_id, float weight = 1);
        void removeAllEdges();
        // int removeEdge(int id, bool directed, Vertex* target_vertex);
        void incrementOutDegree();
        void decrementOutDegree();
        void incrementInDegree();
        void decrementInDegree();
        // Edge* hasEdgeBetween(int target_id);
        // Auxiliar methods
};

#endif 