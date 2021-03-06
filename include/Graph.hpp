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
        int clusters;
        
        std::string outfile_name;

    // Methods
    private:
    
    public:
        // Constructors and destructors

        Graph(int order, int clusters, bool directed = false, bool weighted_edge = false, bool weighted_vertex = false);
        ~Graph();

        // Getters and setters

        int getOrder();
        int getNumberEdges();
        bool isDirected();
        bool isWeightedEdge();
        bool isWeightedVertex();
        Vertex *getVertex(int id);
        void setOutfileName(std::string outfile_name);

        // Other methods

        void insertVertex(int id, float weight = 1);
        void insertEdge(int id, int target_id, float edge_weight = 1, float source_vertex_weight = 1, float target_vertex_weight = 1);
        void insertMissingVertices();
        bool searchVertex(int id);
        void getInfo();
        void printAdjList();
        void saveToDot(std::string function, std::set<std::pair<int, int>> *red_edges = nullptr, std::set<std::pair<int, int>> *gray_edges = nullptr);

        // First assignment especific methods

        void DirectTransitiveClosure(int id);
        void AuxDirectTransitiveClosure(int id, std::set<std::pair<int, int>> *);

        bool Dijkstra(int source_id, int target_id);
        bool Floyd(int source_id, int target_id);

        // Prim

        bool MST_Kruskal();
        bool BFS(int id);
        void topologicalSorting();
        void auxTopologicalSorting(int id, std::map<int, int> &colors, std::list<int> &order, bool *dag);

        // Second assignment especific methods
        
        int OldGreedy(int clusters, float alfa = 0);
        
        int Greedy(float alfa = 0);
        void printGreedyTxt(int cost);
        
        int GreedyRandomizedAdaptative(float alfa, int iterations);
        void printGreedyRATxt(float best_cost, float alfa, int best_it, int seed);
        
        int GreedyRandomizedAdaptativeReactive(float alfas[], int tam_alfa, int iterations, int stack);
        void printGreedyRARTxt(float best_cost, float best_alfa, int best_it, int seed);
    
};

#endif // GRAPH_H_INCLUDED
