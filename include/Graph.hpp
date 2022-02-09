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
        std::unordered_map<int, Vertex *> vertices;
        
        std::string outfile_name;

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
        int Greedy(int clusters, float alfa = 0);
        void printGreedyTxt(std::string file_name, std::string instance_name, int cost, double CPU_time, double wall_time);
        int GreedyRandomizedAdaptative(int clusters, float alfa, int* seed, int* best_it, int iterations);
        void printGreedyRandomizedAdaptativeTxt(std::string file_name, std::string instance_name, int iterations, float alfa, int seed, int best_cost, int best_it, double CPU_time, double wall_time);
    
};

#endif // GRAPH_H_INCLUDED
