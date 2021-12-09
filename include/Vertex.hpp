#ifndef VERTEX_HPP
#define VERTEX_HPP

#include <unordered_map>

#include "Edge.hpp"

class Vertex
{
    // Attributes
    private:
        std::unordered_map<int, Edge*> edges;
        int id;
        unsigned int in_degree;
        unsigned int out_degree;
        float weight;

    // Methods
    private:

    public:
        // Constructors and destructors

        Vertex(int id, float weight = 1);
        ~Vertex();

        // Getters and setters

        int getId();
        int getInDegree();
        int getOutDegree();
        float getWeight();
        void setWeight(float weight);
        std::unordered_map<int, Edge*> getEdges();

        // Other methods

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