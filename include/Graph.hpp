#ifndef GRAPH_H_INCLUDED
#define GRAPH_H_INCLUDED

// #include <fstream>
// #include <stack>
#include <list>

#include <map>
#include <set>
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
    
        // Other methods

        void insertVertex(int id, float weight = 1);
        void insertEdge(int id, int target_id, float edge_weight = 1, float source_vertex_weight = 1, float target_vertex_weight = 1);
        void insertMissingVertices();
        bool searchVertex(int id);
        void getInfo();
        void printAdjList();
        void saveToDot(std::string outfile_name, std::set<std::pair<int,int>> *red_edges = nullptr, std::set<std::pair<int,int>> *gray_edges = nullptr);

        // Assignment especific methods
        
        // A
        // B
        
        bool Dijkstra(int source_id, int target_id);        
        bool Floyd(int source_id, int target_id);
        
        // Prim
        
        bool MST_Kruskal();
        bool BFS(int id);
        void topologicalSorting();
        void auxTopologicalSorting(int id, std::map<int, int>& colors, std::list<int>& order);

};

#endif // GRAPH_H_INCLUDED
