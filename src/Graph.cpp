#include <iostream>
#include <fstream>
#include <unordered_map>
#include <queue>
#include <set>
#include <list>
#include <stack>
#include <list>
#include <map>

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

    bool operator ==(const SimpleEdge& e) const
    {
        if (a != e.a)
            return false;
        if (b != e.b)
            return false;
        return true;
    }
};

/**
 * @brief Kruskal's Minimum Spanning Tree
 * 
 * @return mst_edges Set of tree edges
 */
bool Graph::MST_Kruskal()
{
    std::cout << "Arvore Geradora Minima de Kruskal" << std::endl;
    
    if (vertices.empty()) {
        std::cerr << "No vertices in graph!" << std::endl;
        return false;
    }

    std::unordered_map<int, int> subsets;
    std::list<SimpleEdge> sorted_edges;
    std::set<std::pair<int, int>> *mst_edges = new std::set<std::pair<int, int>>;
    int cost = 0;

    // Create subsets ranging from 0 to n-1
    int i = 0;
    for(std::unordered_map<int, Vertex*>::iterator itV = vertices.begin(); itV != vertices.end(); ++itV) {
        subsets.insert({itV->second->getId(), i}); // Put each vertex in its subset
        i++;

        // Traverse the vertex's edges and add them to the list of edges
        std::unordered_map<int, Edge*> edges = itV->second->getEdges();
        for(std::unordered_map<int, Edge*>::iterator itE = edges.begin(); itE != edges.end(); ++itE) {
            // Create a "simple edge" containing the origin vertex's id, the target vertex's id, and the edge weight
            sorted_edges.push_back(SimpleEdge(itV->second->getId(), itE->second->getTargetId(), itE->second->getWeight(), this->directed));
        }
    }

    sorted_edges.sort(); // Sort the edges in asceding order according to their weights
    sorted_edges.unique(); // Remove duplicates (in undirected graph)
    
    i = 0;
    while(i < order-1 && !sorted_edges.empty())
    {
        SimpleEdge e = sorted_edges.front(); // Get first edge in list
        sorted_edges.pop_front(); // Pop edge from the list of edges

        // Get the IDs of the subsets the vertices a and b are in
        int smallest = subsets.at(e.a);
        int greatest = subsets.at(e.b);
        if(smallest != greatest)
        {
            // std::cout << "(" << e.a << ", " << e.b << ") ";
            mst_edges->insert(std::pair<int, int>(e.a, e.b));
            cost += e.weight; // sum up the edge weight to the MST cost

            if(greatest < smallest)
                std::swap(smallest, greatest);
            
            // Merge the subsets (changing the greatest id to the smallest id)
            for(std::unordered_map<int, int>::iterator it = subsets.begin(); it != subsets.end(); ++it) {
                if(it->second == greatest)
                    it->second = smallest;
            }
            i++;
        }
    }

    std::cout << std::endl << "Custo: " << cost << std::endl;
    saveToDot("mst_kruskal.dot", mst_edges);

    delete mst_edges;
    return true;
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
bool Graph::BFS(int id)
{
    std::cout << "Caminhamento em largura:" << std::endl;

    if (!searchVertex(id)) {
        std::cerr << "Vertex not found!" << std::endl;
        return false;
    }

    // Sets of tree and back edges (which will be printed in different colors in the .dot)
    std::set<std::pair<int, int>>* tree_edges = new std::set<std::pair<int, int>>;
    std::set<std::pair<int, int>>* back_edges = new std::set<std::pair<int, int>>;

    std::unordered_map<int, char> colored_vertices;
    for(std::unordered_map<int, Vertex*>::iterator it = vertices.begin(); it != vertices.end(); ++it) {
        colored_vertices.insert({it->second->getId(), 'w'}); // Insert as not visited (white)
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
        // std::cout << id << " ";
        toVisit.pop();

        // Check for unvisited vertices in adjacency list and enqueue them
        std::unordered_map<int, Edge*> edges = vertices[id]->getEdges();
        for(std::unordered_map<int, Edge*>::iterator itE = edges.begin(); itE != edges.end(); ++itE)
        {
            int target_id = itE->second->getTargetId();

            // If the target vertex is white, then it is unvisited
            if(colored_vertices.at(target_id) == 'w')
            {
                std::pair<int, int> pair_vertices(id, target_id);
                if(!this->directed && pair_vertices.second < pair_vertices.first)
                    std::swap(pair_vertices.first, pair_vertices.second);
                tree_edges->insert(pair_vertices);
                colored_vertices.at(target_id) = 'g'; // mark it as "to visit"
                toVisit.push(target_id);
            }
            else if(colored_vertices.at(target_id) == 'g')
            {
                std::pair<int, int> pair_vertices(id, target_id);
                if(!this->directed && pair_vertices.second < pair_vertices.first)
                    std::swap(pair_vertices.first, pair_vertices.second);
                back_edges->insert(pair_vertices);
            }
        }

        colored_vertices.at(id) = 'b'; // Mark current vertex as visited (black)
    }
    std::cout << std::endl;

    saveToDot("bfs_tree.dot", tree_edges, back_edges);

    delete tree_edges;
    delete back_edges;
    return true;
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
            if(!this->directed && pair_vertices.second < pair_vertices.first)
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

typedef struct
{
    int length;
    std::vector<int> path;
}Length;


/**
 * @brief Prints the smallest path between two vertices 
 //TODO Salvar grafo em .dot seguindo alguma regra pro caminho mínimo (cor diferente da aresta talvez)
 //TODO Grafos de exemplo não ponderados com segmentation fault
 * 
 * @param source_id ID of the starting vertex
 * @param target_id ID of the target vertex
 */
void Graph::Dijkstra(int source_id, int target_id)
{
    //TODO verificação da existência dos nós
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
    
    Length *pi = new Length[vertices.size() + 1];
    for(int i = 0; i < vertices.size(); i++) {
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

    //std::cout << "Inicia algoritmo" << std::endl;
    // Algorithm
    while(non_iterated.size() > 0) {
        // std::cout << "Chegou aqui?" << std::endl;
        // for(int i = 0; i < non_iterated.size(); i++) 
        //     std::cout << non_iterated[i] << " ";
        // std::cout << "\nPI: ";
        // for(int i = 0; i < vertices.size(); i++) 
        //     std::cout << pi[i].length << " ";
        // std::cout << "\n";
        
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
        //std::cout << "Erases " << non_iterated[j] << " at " << j << std::endl;
        non_iterated.erase(non_iterated.begin() + j);
        //std::cout << "Erased, new size: " << non_iterated.size() << std::endl;

        // Iterate through j adjacencies
        std::unordered_map<int, Edge*> edges = vertices[j_id]->getEdges();
        for(std::unordered_map<int, Edge*>::iterator it = edges.begin(); it != edges.end(); ++it) {
            Edge *e = it->second;
            int pi_aux;
            pi_aux = pi[j_id].length + e->getWeight();
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
            //std::cout << "Iteratin' through j vertices" << std::endl;
        }
        //std::cout << "Iterated through j vertices" << std::endl;
    }

    for(int i = 0; i < vertices.size(); i++) {
        std::cout << pi[i].length << " ";
    }
    std::cout << "\n";
    // for(int j = 0; j < vertices.size(); j++)
    // {
    //     for(int i = 0; i < pi[j].path.size(); i++) {
    //         std::cout << pi[j].path[i] << " ";
    //     }
    //     std::cout << "\n";
    // }

    delete [] pi;
}

void Graph::auxTopologicalSorting(int id, std::map<int, int>& colors, std::list<int>& order)
{

    // Iterator for adjacent list
    std::unordered_map<int, Edge*>::iterator itE;

    // Root vertex becomes gray
    colors.at(id) = 1;

    std::unordered_map<int, Edge*> edges = vertices[id]->getEdges();

    for (itE = edges.begin(); itE != edges.end(); ++itE){
      
        // gets ID of adjascent vertex
        int target_id = itE->second->getTargetId();

        if (colors.at(target_id) == 0){
            // adjascent vertex unvisited
            auxTopologicalSorting(target_id, colors, order);
        } else if (colors.at(target_id) == 1) {
            // graph has a return edge so can't have a topological sorting
            std::cout << "Grafo não é um DAG! (grafo ciclico em (" << id << ", " << target_id << ")" << std::endl;
            return;
        }

    }

    // Vertex is terminated -> color becomes black and is added to list
    colors.at(id) = 2;
    order.push_front(id);
}

/**
 * @brief Prints the topological sorting from a DAG (based on the Depth First Search algorithm)
 * 
 */
void Graph::topologicalSorting(){

    std::cout << "Ordenação topológica do grafo:" << std::endl;
    if (!this->directed) {
        std::cout << "Grafo não é um DAG! (não direcionado)" << std::endl;
        return;
    }

    // Iterator for list of vertices
    std::unordered_map<int, Vertex*>::iterator it;
  
    // Map of the vertices Id and a color that indicates wether the vertice is unvisited (0 - white), being visited (1 - gray) or terminated (2 - black)
    std::map<int, int> colors;

    // List of vertices sorted in topological order
    std::list<int> topologicalOrder;

    // Initially all vertices receive color white
    for (it = vertices.begin(); it != vertices.end(); ++it){
        colors.insert({it->second->getId(), 0});
        /* std::cout << colors.at(it->second->getId()); */
    }

    // Loop through all unvisited vertices
    for (it = vertices.begin(); it != vertices.end(); ++it){
        if (colors.at(it->second->getId()) == 0){
            /* std::cout << it->second->getId() << " eh branco" << std::endl; */
            auxTopologicalSorting(it->second->getId(), colors, topologicalOrder);
        }
    }

    while (!topologicalOrder.empty()){
      std::cout << topologicalOrder.front() << " ";
      topologicalOrder.pop_front();
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
