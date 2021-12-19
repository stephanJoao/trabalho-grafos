#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <fstream>
#include <set>
// #include <stack>
#include <list>
// #include <math.h>
// #include <cstdlib>
// #include <ctime>
// #include <float.h>
// #include <iomanip>

#include "../include/Graph.hpp"
#include "../include/Vertex.hpp"
#include "../include/Edge.hpp"

//* Constructors and destructors implementations

Graph::Graph(int order, bool directed, bool weighted_edge, bool weighted_vertex)
{
    this->order = order;
    this->directed = directed;
    this->weighted_edge = weighted_edge;
    this->weighted_vertex = weighted_vertex;
    this->number_edges = 0;
}

Graph::~Graph()
{
    for (std::unordered_map<int, Vertex*>::iterator it = vertices.begin(); it != vertices.end(); ++it)
        delete it->second;
}

//* Getters and setters implementations

int Graph::getOrder()
{
    return this->order;
}

int Graph::getNumberEdges()
{
    return this->number_edges;
}

// Function that verifies if the graph is directed
bool Graph::isDirected()
{
    return this->directed;
}

// Function that verifies if the graph is weighted at the edges
bool Graph::isWeightedEdge()
{
    return this->weighted_edge;
}

// Function that verifies if the graph is weighted at the vertexs
bool Graph::isWeightedVertex()
{
    return this->weighted_vertex;
}

Vertex* Graph::getVertex(int id)
{
    return vertices[id];
}

//* Other methods

/*
    The outdegree attribute of vertices is used as a counter for the number of edges in the graph.
    This allows the correct updating of the numbers of edges in the graph being directed or not.
 */

/**
 * Insert vertex in graph.
 *
 * @param id Edge id
 * @param directed If the graph is directed
 * @param target_vertex Target vertex
 */
void Graph::insertVertex(int id, float weight)
{
    if (vertices.count(id) == 0) {
        Vertex *v = new Vertex(id, weight);
        vertices.insert({id, v});
    }
}

/**
 * @brief Insert edge on the vertex edges list
 * 
 * @param id First vertex's ID
 * @param target_id Second Vertex's ID
 * @param weight Weight of the edge
 */
void Graph::insertEdge(int id, int target_id, float edge_weight, float source_vertex_weight, float target_vertex_weight)
{
    if (id == target_id) {
        std::cerr << "Vertices cannot have equal ID!" << std::endl;
        return;
    }

    insertVertex(id, source_vertex_weight);
    insertVertex(target_id, target_vertex_weight);
    
    if (this->directed){
        vertices[id]->insertEdge(target_id, edge_weight);
    } else {
        vertices[id]->insertEdge(target_id, edge_weight);
        vertices[target_id]->insertEdge(id, edge_weight);
    }
}

// void Graph::removeVertex(int id)
// {

// }

/**
 * @brief Return if the vertex is in the graph
 * 
 * @param id Vertex's ID
 * @return true if
 * @return false 
 */
bool Graph::searchVertex(int id)
{
    return vertices[id] != nullptr;
}

// Print adjacency list
void Graph::printAdjList()
{
    std::cout << "Imprimindo Grafo como uma lista de adjacência:" << std::endl;
    for(std::unordered_map<int, Vertex*>::iterator itV = vertices.begin(); itV != vertices.end(); ++itV) {
        Vertex *v = itV->second;
        std::cout << v->getId();
        std::unordered_map<int, Edge*> edges = v->getEdges();
        for(std::unordered_map<int, Edge*>::iterator itE = edges.begin(); itE != edges.end(); ++itE) {
            Edge *e = itE->second;
            std::cout << " -> " << e->getTargetId();
        }
        std::cout << std::endl;
    }
}

struct SimpleEdge
{
    int a;
    int b;
    int weight;

    SimpleEdge(int a, int b, int weight, bool directed)
    {
        if(!directed && b < a)
            std::swap(a, b);
        
        this->a = a;
        this->b = b;
        this->weight = weight;
    }

    bool operator <(const SimpleEdge & e) const
    {
        return weight < e.weight;
    }

    bool operator==(const SimpleEdge& e) const
    {
        if (a != e.a)
            return false;
        if (b != e.b)
            return false;
        return true;
    }
};

inline bool operator<(const SimpleEdge& e1, const SimpleEdge& e2)
{
    return e1.a < e2.a;
}

/**
 * @brief Kruskal's Minimum Spanning Tree
 * 
 * @return mst_edges Set of tree edges
 */
std::set<std::pair<int, int>>* Graph::MST_Kruskal()
{
    std::cout << "Arvore Geradora Minima de Kruskal" << std::endl;
    
    if (vertices.empty()) {
        std::cerr << "No vertices in graph!" << std::endl;
        return nullptr;
    }

    std::unordered_map<int, int> subsets;
    std::list<SimpleEdge> sorted_edges;
    std::set<std::pair<int, int>> *mst_edges = new std::set<std::pair<int, int>>;
    int cost = 0;

    int i = 0;
    for(std::unordered_map<int, Vertex*>::iterator itV = vertices.begin(); itV != vertices.end(); ++itV) {
        subsets.insert({itV->second->getId(), i});
        i++;

        std::unordered_map<int, Edge*> edges = itV->second->getEdges();
        for(std::unordered_map<int, Edge*>::iterator itE = edges.begin(); itE != edges.end(); ++itE) {
            sorted_edges.push_back(SimpleEdge(itV->second->getId(), itE->second->getTargetId(), itE->second->getWeight(), directed));
        }
    }

    sorted_edges.sort();
    sorted_edges.unique();
    
    i = 0;
    while(i < order-1 && !sorted_edges.empty())
    {
        SimpleEdge e = sorted_edges.front();
        sorted_edges.pop_front();

        int smallest = subsets.at(e.a);
        int greatest = subsets.at(e.b);
        if(smallest != greatest)
        {
            mst_edges->insert(std::pair<int, int>(e.a, e.b));
            cost += e.weight;

            if(greatest < smallest)
                std::swap(smallest, greatest);
            
            for(std::unordered_map<int, int>::iterator it = subsets.begin(); it != subsets.end(); ++it) {
                if(it->second == greatest)
                    it->second = smallest;
            }
            i++;
        }
    }

    std::cout << "Custo: " << cost << std::endl;
    return mst_edges;
}

/**
 * @brief Print the vertices in order of visit according to the Breadth First Search algorithm. 
 * Color legend:
 * - White: unvisited
 * - Gray: to visit
 * - Black: visited
 * 
 * @param id ID of the starting vertex
 * @param back_edges Set of back edges (will be overwritten)
 * @return tree_edges Set of tree edges
 */
std::set<std::pair<int, int>>* Graph::BFS(int id, std::set<std::pair<int, int>>* back_edges)
{
    std::cout << "Caminhamento em largura:" << std::endl;

    if (!searchVertex(id)) {
        std::cerr << "Vertex not found!" << std::endl;
        return nullptr;
    }

    // Set of tree edges (which will be printed in a different color in the .dot)
    std::set<std::pair<int, int>>* tree_edges = new std::set<std::pair<int, int>>;

    std::unordered_map<int, char> colored_vertices;
    for(std::unordered_map<int, Vertex*>::iterator it = vertices.begin(); it != vertices.end(); ++it) {
        colored_vertices.insert({it->second->getId(), 'w'}); // Insert as not visited vertex (white)
    }

    colored_vertices.at(id) = 'g'; // Vertex v is the first to be visited (gray)

    // Queue of vertices to visit
    std::queue<int> toVisit;

    // Enqueue current vertex and mark it as "to visit"
    toVisit.push(id);

    while(!toVisit.empty())
    {
        // Dequeue top vertex in queue
        id = toVisit.front();
        std::cout << id << " ";
        toVisit.pop();

        // Check for unvisited vertices in adjacency list and enqueue them
        std::unordered_map<int, Edge*> edges = vertices[id]->getEdges();
        for(std::unordered_map<int, Edge*>::iterator itE = edges.begin(); itE != edges.end(); ++itE)
        {
            int target_id = itE->second->getTargetId();
            if(colored_vertices.at(target_id) == 'w')
            {
                std::pair<int, int> pair_vertices(id, target_id);
                if(pair_vertices.second < pair_vertices.first)
                    std::swap(pair_vertices.first, pair_vertices.second);
                tree_edges->insert(pair_vertices);
                colored_vertices.at(target_id) = 'g';
                toVisit.push(target_id);
            }
            else if(colored_vertices.at(target_id) == 'g')
            {
                std::pair<int, int> pair_vertices(id, target_id);
                if(pair_vertices.second < pair_vertices.first)
                    std::swap(pair_vertices.first, pair_vertices.second);
                back_edges->insert(pair_vertices);
            }
        }

        colored_vertices.at(id) = 'b'; // Mark current vertex as visited (black)
    }
    std::cout << std::endl;

    return tree_edges;
}

/**
 * @brief Save the graph in .dot
 * 
 * @param outfile_name Name of the file to which the graph will be saved
 */
void Graph::saveToDot(std::string outfile_name, std::set<std::pair<int,int>> *red_edges, 
std::set<std::pair<int,int>> *dashed_edges)
{
    std::ofstream outfile(outfile_name);
    std::string edgeop; // "->" or "--"

    if (this->directed)
    {
        outfile << "digraph Grafo {" << std::endl;
        edgeop = " -> ";
    }
    else
    {
        outfile << "graph Grafo {" << std::endl;
        edgeop = " -- ";
    }

    if(isWeightedVertex()) {
        for(std::unordered_map<int, Vertex*>::iterator itV = vertices.begin(); itV != vertices.end(); ++itV) {
            Vertex *v = itV->second;
            outfile << v->getId() << " [weight=" << v->getWeight() << "];" << std::endl;
        }
        std::cout << std::endl;
    }

    // Set of "edges" that have already been written to the file
    std::set<std::pair<int,int>> included;

    // Iterate through vertices and edges
    for(std::unordered_map<int, Vertex*>::iterator itV = vertices.begin(); itV != vertices.end(); ++itV) {
        Vertex *v = itV->second;
        std::unordered_map<int, Edge*> edges = v->getEdges();
        for(std::unordered_map<int, Edge*>::iterator itE = edges.begin(); itE != edges.end(); ++itE) {
            Edge *e = itE->second;
            
            // Create "edge" in which the smallest vertex comes first
            std::pair<int,int> pair_vertices(v->getId(), e->getTargetId());
            if(pair_vertices.second < pair_vertices.first)
                std::swap(pair_vertices.first, pair_vertices.second);

            // If this "edge" hasn't been included yet
            // Then mark as included and write it to the file
            if(!included.count(pair_vertices)) {
                included.insert(pair_vertices);
                outfile << v->getId() << edgeop << e->getTargetId();
                if(isWeightedEdge()) {
                    outfile << " [weight=" << e->getWeight();
                    if(red_edges != nullptr && red_edges->count(pair_vertices))
                        outfile << " color=\"red\"";
                    else if(dashed_edges != nullptr && dashed_edges->count(pair_vertices))
                        outfile << " style=\"dashed\"";
                    outfile << "]";
                } else { 
                    if(red_edges != nullptr && red_edges->count(pair_vertices))
                        outfile << " [color=\"red\"]";
                    else if(dashed_edges != nullptr && dashed_edges->count(pair_vertices))
                        outfile << " [style=\"dashed\"]";
                }
                outfile << ";" << std::endl;
            }
        }
    }
    
    outfile << "}" << std::endl;
}

// //Function that prints a set of edges belongs breadth tree

// void Graph::breadthFirstSearch(std::ofstream &output_file){
    
// }



// float Graph::floydMarshall(int idSource, int idTarget){
    
// }

   

// float Graph::dijkstra(int idSource, int idTarget){
    
// }

// //function that prints a topological sorting
// void topologicalSorting(){

// }

// void breadthFirstSearch(std::ofstream& output_file){

// }

// Graph* getVertexInduced(int* listIdVertexs){

// }

// Graph* agmKuskal(){

// }
// Graph* agmPrim(){

// }
