#ifndef GRAPH_H_INCLUDED
#define GRAPH_H_INCLUDED

// #include <fstream>
// #include <stack>
// #include <list>

#include <unordered_map>

#include "Vertex.hpp"

class Graph
{
    // Attributes
    private:
        int order;
        int number_edges;
        bool directed;
        bool weighted_edge;
        bool weighted_vertex;
        std::unordered_map<int, Vertex*> vertices;

    // Methods
    private:

    public:
        // Constructors and destructors

        Graph(int order, bool directed = false, bool weighted_edge = false, bool weighted_vertex = false);
        ~Graph();

        // Getters and setters

        int getOrder();
        int getNumberEdges();
        bool isDirected();
        bool isWeightedEdge();
        bool isWeightedVertex();
        Vertex* getVertex(int id);
    
        //Other methods

        void insertVertex(int id, float weight = 1);
        void insertEdge(int id, int target_id, float weight = 1, float vertex_weight = 1);
        // void removeVertex(int id);
        bool searchVertex(int id);

        void printAdjList();
        void BFS(int id);
        void saveToDot(std::string outfile_name);
    //     //methods phase1

    //     void topologicalSorting();
    //     void breadthFirstSearch(std::ofstream& output_file);
    //     Graph* getVertexInduced(int* listIdVertexs);
    //     Graph* agmKuskal();
    //     Graph* agmPrim();
    //     float floydMarshall(int idSource, int idTarget);
    //     float dijkstra(int idSource, int idTarget);

    //     //methods phase1

    //     float greed();
    //     float greedRandom();
    //     float greedRactiveRandom();
    // private:
    //     //Auxiliar methods

};

#endif // GRAPH_H_INCLUDED