#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <fstream>
#include <set>
// #include <stack>
// #include <list>
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
        std::cerr << "Vertices cannot have equal ID! (Selfloop not allowed)" << std::endl;
        return;
    }

    if (id >= order) {
        order = id;
    }
    if (target_id >= order) {
        order = target_id;
    }

    insertVertex(id, source_vertex_weight);
    insertVertex(target_id, target_vertex_weight);
    
    if (this->directed){
        // if(vertices.find(id) != vertices.end()){
        //     std::cout << "Vertex already present" << std::endl;
        //     if(vertices[id]->getEdges().find(target_id) != vertices[id]->getEdges().end())
        //         std::cout << "Edge already present" << std::endl;
        // }
        vertices[id]->insertEdge(target_id, edge_weight);
    } else {
        // if(vertices.find(id) != vertices.end()){
        //     std::cout << "Vertex already present" << std::endl;
        //     if(vertices[id]->getEdges().find(target_id) != vertices[id]->getEdges().end())
        //         std::cout << "Edge already present" << std::endl;
        // }
        vertices[id]->insertEdge(target_id, edge_weight);
        // if(vertices.find(target_id) != vertices.end()){
        //     std::cout << "Vertex already present" << std::endl;
        //     if(vertices[target_id]->getEdges().find(id) != vertices[target_id]->getEdges().end())
        //         std::cout << "Edge already present" << std::endl;
        // }
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

void Graph::getInfo() {
    std::cout << "\nVERIFYING PROVIDED GRAPH" << std::endl;
    std::cout << "  Provided order = " << order << std::endl;
    std::cout << "  Real number of vertices = " << vertices.size() << std::endl;
    std::cout << "  Missing vertices (disconnected):\n        ";
    for(int i = 0; i < order+1; i++) {
        if (vertices.count(i) == 0) {
            std::cout << i << " ";
        }
    }
    std::cout << "\n";
}

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

/**
 * @brief Save the graph in .dot
 * 
 * @param outfile_name Name of the file to which the graph will be saved
 */
void Graph::saveToDot(std::string outfile_name)
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
            int greatest = v->getId();
            int smallest = e->getTargetId();
            if(greatest < smallest)
                std::swap(greatest, smallest);

            std::pair<int,int> pair_vertices(smallest, greatest);

            // If this "edge" hasn't been included yet
            // Then mark as included and write it to the file
            if(!included.count(pair_vertices)) {
                included.insert(pair_vertices);
                outfile << v->getId() << edgeop << e->getTargetId();
                if(isWeightedEdge()) {
                    outfile << " [weight=" << e->getWeight() << "]";
                }
                outfile << ";" << std::endl;
            }
        }
    }
    
    outfile << "}" << std::endl;
}

typedef struct
{
    int length;
    std::vector<int> path;
}Length;


/**
 * @brief Calculates the smallest path between two vertices using Dijkstra's algorithm
 //TODO Salvar grafo em .dot seguindo alguma regra pro caminho mínimo (cor diferente da aresta talvez)
 //TODO Grafos de exemplo não ponderados com segmentation fault
 * 
 * @param source_id ID of the starting vertex
 * @param target_id ID of the target vertex
 */
void Graph::Dijkstra(int source_id, int target_id)
{
    //TODO verificação da existência dos nós
    
    int n = vertices.size();

    // Initializes non iterated vertices
    //std::cout << "Initializes non iterated vertices" << std::endl;
    std::vector<int> non_iterated;
    for(std::unordered_map<int, Vertex*>::iterator it = vertices.begin(); it != vertices.end(); ++it) {
        Vertex *v = it->second;
        if(v->getId() != source_id)
            non_iterated.push_back(v->getId());
    }
    
    // Initializes iterated vertices
    //std::cout << "Initializes iterated vertices" << std::endl;
    std::vector<int> iterated;
    iterated.push_back(source_id);

    // Initializes vector of lengths and paths (called pi)
    //std::cout << "Initializes pi" << std::endl;    
    Length *pi = new Length[n];
    //std::cout << "n = " << order+1 << std::endl;
    for(int i = 0; i < n; i++) {
        pi[i].length = std::numeric_limits<int>::max();
    }
    pi[source_id].length = 0;
    std::unordered_map<int, Edge*> edges = vertices[source_id]->getEdges();
    for(std::unordered_map<int, Edge*>::iterator it = edges.begin(); it != edges.end(); ++it) {
        Edge *e = it->second;
        pi[e->getTargetId()].length = e->getWeight();
        pi[e->getTargetId()].path.push_back(source_id);
        pi[e->getTargetId()].path.push_back(e->getTargetId());
    }

    std::cout << "Inicia algoritmo" << std::endl;
    // Algorithm
    while(non_iterated.size() > 0) {
        std::cout << "Chegou aqui?" << std::endl;
        for(int i = 0; i < non_iterated.size(); i++) 
            std::cout << non_iterated[i] << " ";
        std::cout << "\nPI: ";
        for(int i = 0; i < vertices.size(); i++) 
            std::cout << pi[i].length << " ";
        std::cout << "\n";
        
        // Find j with minimal path cost in non_iterated
        int j = 0;
        int j_length = std::numeric_limits<int>::max();
        int j_id = 0;
        for(int i = 0 ; i < non_iterated.size(); i++) {
            if(pi[non_iterated[i]].length < j_length) {
                j = i;
                j_id = non_iterated[j];
                j_length = pi[j_id].length;
            }
        }
        // Removes this j from the non_iterated
        std::cout << "Erases " << non_iterated[j] << " at " << j << std::endl;
        non_iterated.erase(non_iterated.begin() + j);
        std::cout << "Erased, new size: " << non_iterated.size() << std::endl;

        // Iterate through j adjacencies
        std::unordered_map<int, Edge*> edges = vertices[j_id]->getEdges();
        for(std::unordered_map<int, Edge*>::iterator it = edges.begin(); it != edges.end(); ++it) {
            Edge *e = it->second;
            int pi_aux;
            if(pi[j_id].length != std::numeric_limits<int>::max())
                pi_aux = pi[j_id].length + e->getWeight();
            else
                pi_aux = std::numeric_limits<int>::max();
            if(pi_aux < pi[e->getTargetId()].length) {
                pi[e->getTargetId()].length = pi_aux;
                pi[e->getTargetId()].path = pi[j_id].path;
                pi[e->getTargetId()].path.push_back(e->getTargetId());
                bool non = false;
                for(int i = 0 ; i < non_iterated.size(); i++){
                    if(non_iterated[i] == e->getTargetId()) {
                        non = true;
                        break;
                    }
                }
                if(!non)
                    non_iterated.push_back(e->getTargetId());
            }
            std::cout << "Iteratin' through j vertices" << std::endl;
        }
        std::cout << "Iterated through j vertices" << std::endl;
    }

    for(int i = 0; i < n; i++) {
        std::cout << pi[i].length << " ";
    }
    std::cout << "\n";
    
    for(int j = 0; j < n; j++)
    {
        if(pi[j].path.size() == 0)
            std::cout << "Caminho vazio";
        for(int i = 0; i < pi[j].path.size(); i++) {
            std::cout << pi[j].path[i] << " ";
        }
        std::cout << "\n";
    }

    std::set<std::pair<int, int>> *path = new std::set<std::pair<int, int>>();
    for(int i = 0; i < pi[source_id].path.size() - 1; i++){
        int x = pi[source_id].path[i], y = pi[source_id].path[i];
        if(pi[source_id].path[i] < pi[source_id].path[i + 1])
            y = pi[source_id].path[i + 1];
        else
            x = pi[source_id].path[i + 1];
        path->insert({x, y});
    }
    //saveToDot("dijkstra.dot");


    delete [] pi;
}

/**
 * @brief Calculates the smallest path between two vertices using Floyd's algorithm
 //TODO Salvar grafo em .dot seguindo alguma regra pro caminho mínimo (cor diferente da aresta talvez)
 * 
 * @param source_id ID of the starting vertex
 * @param target_id ID of the target vertex
 */
void Graph::Floyd(int source_id, int target_id)
{
    //TODO verificação da existência dos nós
    
    int n = vertices.size();

    Length A[n][n];
    Length A_previous[n][n];

    // Creation of matrix A_0
    for(int i = 0; i < n; i++){
        for(int j = 0; j < n; j++) {
            if(i != j) {
                A[i][j].length = 1;
                if(vertices.count(i) == 0){
                    A[i][j].length = std::numeric_limits<int>::max();
                }
                else {
                    if(vertices[i]->getEdges().count(j) == 0)
                        A[i][j].length = std::numeric_limits<int>::max();
                    else                
                        A[i][j].length = vertices[i]->getEdges()[j]->getWeight(); 
                }
            }
            else {
                A[i][j].length = 0;
            }
        }
    }

    for(int x = 0; x < n; x++) {
        for(int y = 0; y < n; y++) {
            printf("%12d ", A[x][y].length);
        }
        printf("\n");
    }
    printf("\n");

    // Iterations through A matrix according to vertices allowed
    for(int k = 0; k < n; k++) {
        for(int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) {
                A_previous[i][j].length = A[i][j].length;
            }
        }
        for(int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) {
                if(i != j) {
                    if(A_previous[i][j].length > A_previous[i][k].length + A_previous[k][j].length)
                        A[i][j].length = A_previous[i][k].length + A_previous[k][j].length;
                }
            }
        }
        for(int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) {
                printf("%2d ", A[i][j].length);
            }
            printf("\n");
        }
        printf("\n");
    }

    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            printf("%2d ", A[i][j].length);
        }
        printf("\n");
    }
}

/**
 * @brief Prints the vertices in order according to the Breadth First Search algorithm
 * 
 * @param id ID of the starting vertex
 */
void Graph::BFS(int id)
{
    std::cout << "Caminhamento em largura:" << std::endl;
    if (!searchVertex(id)) {
        std::cerr << "Vertex not found!" << std::endl;
        return;
    }
    
    // Hash Set of visited vertices
    std::unordered_set<int> visited;
 
    // Queue of vertices to visit
    std::queue<int> toVisit;
 
    // Enqueue current vertex and add it to the visited vertices set
    toVisit.push(id);
    visited.insert(id);
 
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
            if(!visited.count(target_id))
            {
                visited.insert(target_id);
                toVisit.push(target_id);
            }
        }
    }
    std::cout << std::endl;
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
