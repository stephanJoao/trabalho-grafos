#ifndef EDGE_HPP
#define EDGE_HPP

class Vertex;
class Edge
{
    // Attributes
    private:
        Vertex* target_vertex;
        float weight;
        Edge* next_edge;

    // Methods
    private:

    public:
        // Constructors and destructors

        Edge(Vertex* target_vertex, Edge* next_edge, float weight = 1);
        ~Edge();

        // Getters and setters

        Vertex* getTargetVertex();
        Edge* getNextEdge();
        void setNextEdge(Edge* edge);
        float getWeight();
        void setWeight(float weight);

        // Other methods
};

#endif 