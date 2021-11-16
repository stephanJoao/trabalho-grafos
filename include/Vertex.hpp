#ifndef VERTEX_HPP
#define VERTEX_HPP

class Vertex
{
    private:
        int value;
        /* data */
    public:
        Vertex();
        Vertex(int value);
        ~Vertex();

        int getValue();
};

#endif 