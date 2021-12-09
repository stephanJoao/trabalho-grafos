#ifndef EDGE_HPP
#define EDGE_HPP

class Edge
{
    // Attributes
    private:
        int target_id;
        float weight;

    // Methods
    private:

    public:
        // Constructors and destructors

        Edge(int target_id, float weight = 1);
        ~Edge();

        // Getters and setters

        int getTargetId();
        float getWeight();
        void setWeight(float weight);

        // Other methods
};

#endif 