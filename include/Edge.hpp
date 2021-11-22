#ifndef EDGE_HPP
#define EDGE_HPP

class Edge
{
    // Attributes
    private:
        int target_id;
        float weight;
        Edge* next_edge;

    // Methods
    private:

    public:
        // Constructors and destructors

        Edge(int target_id, Edge* next_edge, float weight = 1);
        ~Edge();

        // Getters and setters

        int getTargetId();
        Edge* getNextEdge();
        void setNextEdge(Edge* edge);
        float getWeight();
        void setWeight(float weight);

        // Functions
};

#endif 