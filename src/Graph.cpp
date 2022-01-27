#include <iostream>
#include <fstream>
#include <unordered_map>
#include <queue>
#include <set>
#include <list>
#include <map>
#include <string>

#include "../include/Graph.hpp"
#include "../include/Vertex.hpp"
#include "../include/Edge.hpp"

#define INFINITY std::numeric_limits<int>::max()

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
    for (std::unordered_map<int, Vertex *>::iterator it = vertices.begin(); it != vertices.end(); ++it)
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

Vertex *Graph::getVertex(int id)
{
    return vertices[id];
}

void Graph::setOutfileName(std::string outfile_name)
{
    this->outfile_name = "./output/" + outfile_name;

    // If filename doesn't have an extension then add ".dot" in the end
    if (outfile_name.rfind('.') == std::string::npos)
        this->outfile_name += ".dot";
}

//* Other methods

/**
 * Insert vertex in graph.
 *
 * @param id Edge id
 * @param directed If the graph is directed
 * @param target_vertex Target vertex
 */
void Graph::insertVertex(int id, float weight)
{
    if (vertices.count(id) == 0)
    {
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
    if (id == target_id)
    {
        std::cerr << "Vertices cannot have equal ID! (Selfloop not allowed)" << id << std::endl;
        return;
    }

    if (id >= order)
    {
        order = id;
    }
    if (target_id >= order)
    {
        order = target_id;
    }

    insertVertex(id, source_vertex_weight);
    insertVertex(target_id, target_vertex_weight);

    if (this->directed)
    {
        // if(vertices.find(id) != vertices.end()){
        //     std::cout << "Vertex already present" << std::endl;
        //     if(vertices[id]->getEdges().find(target_id) != vertices[id]->getEdges().end())
        //         std::cout << "Edge already present" << std::endl;
        // }
        vertices[id]->insertEdge(target_id, edge_weight);
    }
    else
    {
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

/**
 * @brief Insert missing vertices considering the graph having vertices from 0 to n - 1.
 */
void Graph::insertMissingVertices()
{
    for (int i = 0; i < order; i++)
    {
        if (vertices.count(i) == 0)
            insertVertex(i);
    }
}

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

void Graph::getInfo()
{
    std::cout << "\nVERIFYING PROVIDED GRAPH" << std::endl;
    std::cout << "  Provided order = " << order << std::endl;
    std::cout << "  Real number of vertices = " << vertices.size() << std::endl;
    std::cout << "  Missing vertices (disconnected):\n        ";
    for (int i = 0; i < order + 1; i++)
    {
        if (vertices.count(i) == 0)
        {
            std::cout << i << " ";
        }
    }
    std::cout << "\n";
}

void Graph::printAdjList()
{
    std::cout << "Imprimindo Grafo como uma lista de adjacência:" << std::endl;
    for (std::unordered_map<int, Vertex *>::iterator itV = vertices.begin(); itV != vertices.end(); ++itV)
    {
        Vertex *v = itV->second;
        std::cout << v->getId();
        std::unordered_map<int, Edge *> edges = v->getEdges();
        for (std::unordered_map<int, Edge *>::iterator itE = edges.begin(); itE != edges.end(); ++itE)
        {
            Edge *e = itE->second;
            std::cout << " -> " << e->getTargetId();
        }
        std::cout << std::endl;
    }
}

/**
 * @brief Save the graph in .dot
 * 
 */
void Graph::saveToDot(std::string function, std::set<std::pair<int, int>> *red_edges, std::set<std::pair<int, int>> *dashed_edges)
{
    std::ofstream outfile;
    outfile.open(outfile_name, std::ios::app);
    std::string edgeop; // "->" or "--"

    if (this->directed)
    {
        outfile << "digraph " << function << " {" << std::endl;
        edgeop = " -> ";
    }
    else
    {
        outfile << "graph " << function << " {" << std::endl;
        edgeop = " -- ";
    }

    if (isWeightedVertex())
    {
        for (std::unordered_map<int, Vertex *>::iterator itV = vertices.begin(); itV != vertices.end(); ++itV)
        {
            Vertex *v = itV->second;
            outfile << "\t" << v->getId() << " [weight=" << v->getWeight() << "];" << std::endl;
        }
        std::cout << std::endl;
    }

    // Set of "edges" that have already been written to the file
    std::set<std::pair<int, int>> included;

    // Iterate through vertices and edges
    for (std::unordered_map<int, Vertex *>::iterator itV = vertices.begin(); itV != vertices.end(); ++itV)
    {
        Vertex *v = itV->second;
        std::unordered_map<int, Edge *> edges = v->getEdges();
        for (std::unordered_map<int, Edge *>::iterator itE = edges.begin(); itE != edges.end(); ++itE)
        {
            Edge *e = itE->second;

            // Create "edge" in which the smallest vertex comes first
            std::pair<int, int> pair_vertices(v->getId(), e->getTargetId());
            if (!this->directed && pair_vertices.second < pair_vertices.first)
                std::swap(pair_vertices.first, pair_vertices.second);

            // If this "edge" hasn't been included yet
            // Then mark as included and write it to the file
            if (!included.count(pair_vertices))
            {
                included.insert(pair_vertices);
                outfile << "\t" << v->getId() << edgeop << e->getTargetId();
                if (isWeightedEdge())
                {
                    outfile << " [weight=" << e->getWeight();
                    if (red_edges != nullptr && red_edges->count(pair_vertices))
                        outfile << " color=\"red\"";
                    else if (dashed_edges != nullptr && dashed_edges->count(pair_vertices))
                        outfile << " style=\"dashed\"";
                    outfile << "]";
                }
                else
                {
                    if (red_edges != nullptr && red_edges->count(pair_vertices))
                        outfile << " [color=\"red\"]";
                    else if (dashed_edges != nullptr && dashed_edges->count(pair_vertices))
                        outfile << " [style=\"dashed\"]";
                }
                outfile << ";" << std::endl;
            }
        }
    }

    outfile << "}" << std::endl;
}

//* Assignment specific methods

typedef struct
{
    int length;
    std::vector<int> path;
} Length;

/**
 * @brief Generates a vertex induced graph from the direct transitive closure
 * 
 * @param id ID of the vertex to find transitive closure
 */
void Graph::DirectTransitiveClosure(int id)
{
    if (this->directed)
    {
        std::set<std::pair<int, int>> *visited = new std::set<std::pair<int, int>>;

        // Check for unvisited vertices in adjacency list and enqueue them
        std::unordered_map<int, Edge *> edges = vertices[id]->getEdges();
        for (std::unordered_map<int, Edge *>::iterator itE = edges.begin(); itE != edges.end(); ++itE)
        {
            int target_id = itE->second->getTargetId();
            std::pair<int, int> pair_vertices(id, target_id);
            if (visited->count(pair_vertices) == 0)
            {
                visited->insert(pair_vertices);
                // std::cout << "(" << pair_vertices.first << ", " << pair_vertices.second << "), " << std::endl;
            }
            this->AuxDirectTransitiveClosure(target_id, visited);
        }

        std::cout << "Fecho transitivo direto: " << std::endl;
        std::cout << "  { ";
        for (std::set<std::pair<int, int>>::iterator it = visited->begin(); it != visited->end(); ++it)
            std::cout << "(" << it->first << ", " << it->second << "), ";
        std::cout << " }\n\n";

        saveToDot("DirectTransitiveClosure", visited);
    }
    else
    {
        std::cout << "Grafo não direcionado" << std::endl;
        return;
    }
}
void Graph::AuxDirectTransitiveClosure(int id, std::set<std::pair<int, int>> *visited)
{
    // Check for unvisited vertices in adjacency list 
    std::unordered_map<int, Edge*> edges = vertices[id]->getEdges();
    for (std::unordered_map<int, Edge *>::iterator itE = edges.begin(); itE != edges.end(); ++itE)
    {
        int target_id = itE->second->getTargetId();
        std::pair<int, int> pair_vertices(id, target_id);
        if (visited->count(pair_vertices) == 0)
        {
            // std::cout << "(" << pair_vertices.first << ", " << pair_vertices.second << "), " << std::endl;
            visited->insert(pair_vertices);
        }
        
        this->AuxDirectTransitiveClosure(target_id, visited);
    }
}

/**
 * @brief Generates a vertex induced graph from the indirect transitive closure
 * 
 * @param id ID of the vertex to find transitive closure
 */


/**
 * @brief Calculates the smallest path between two vertices using Dijkstra's algorithm, then appends it in a .dot file
 * 
 * @param source_id ID of the starting vertex
 * @param target_id ID of the target vertex
 */
bool Graph::Dijkstra(int source_id, int target_id)
{
    // Verifies if the vertices exist in the graph (returns false -> not executed)
    if (vertices.count(source_id) == 0 || vertices.count(target_id) == 0)
    {
        std::cerr << "ERROR: Dijkstra not executed, one or both of the vertices are not present in the graph" << std::endl;
        return false;
    }

    int n = vertices.size();

    // Initializes non iterated vertices
    std::vector<int> non_iterated;
    for (std::unordered_map<int, Vertex *>::iterator it = vertices.begin(); it != vertices.end(); ++it)
    {
        Vertex *v = it->second;
        if (v->getId() != source_id)
            non_iterated.push_back(v->getId());
    }

    // Initializes iterated vertices
    std::vector<int> iterated;
    iterated.push_back(source_id);

    // Initializes vector of lengths and paths (called pi)
    std::unordered_map<int, Length> pi;
    for (std::unordered_map<int, Vertex *>::iterator it = vertices.begin(); it != vertices.end(); ++it)
    {
        Vertex *v = it->second;
        if (v->getId() == source_id)
        {
            Length l;
            l.length = 0;
            pi.insert({v->getId(), l});
        }
        else
        {
            Length l;
            l.length = INFINITY;
            pi.insert({v->getId(), l});
        }
    }
    std::unordered_map<int, Edge *> edges = vertices[source_id]->getEdges();
    for (std::unordered_map<int, Edge *>::iterator it = edges.begin(); it != edges.end(); ++it)
    {
        Edge *e = it->second;
        pi[e->getTargetId()].length = e->getWeight();
        pi[e->getTargetId()].path.push_back(source_id);
        pi[e->getTargetId()].path.push_back(e->getTargetId());
    }

    // Algorithm
    while (non_iterated.size() > 0)
    {
        // Find j with minimal path cost in non_iterated
        int j = 0;
        int j_length = INFINITY;
        int j_id = 0;
        for (int i = 0; i < non_iterated.size(); i++)
        {
            if (pi[non_iterated[i]].length < j_length)
            {
                j = i;
                j_id = non_iterated[j];
                j_length = pi[j_id].length;
            }
        }
        // Removes this j from the non_iterated
        non_iterated.erase(non_iterated.begin() + j);

        // Iterate through j adjacencies
        std::unordered_map<int, Edge *> edges = vertices[j_id]->getEdges();
        for (std::unordered_map<int, Edge *>::iterator it = edges.begin(); it != edges.end(); ++it)
        {
            Edge *e = it->second;
            int pi_aux;
            if (pi[j_id].length != INFINITY)
                pi_aux = pi[j_id].length + e->getWeight();
            else
                pi_aux = INFINITY;
            if (pi_aux < pi[e->getTargetId()].length)
            {
                pi[e->getTargetId()].length = pi_aux;
                pi[e->getTargetId()].path = pi[j_id].path;
                pi[e->getTargetId()].path.push_back(e->getTargetId());
                bool non = false;
                for (int i = 0; i < non_iterated.size(); i++)
                {
                    if (non_iterated[i] == e->getTargetId())
                    {
                        non = true;
                        break;
                    }
                }
                if (!non)
                    non_iterated.push_back(e->getTargetId());
            }
        }
    }

    std::set<std::pair<int, int>> *path = new std::set<std::pair<int, int>>;
    if (pi[target_id].path.size() > 0)
    {
        std::cout << "  Path found by Dijkstra:" << std::endl;
        std::cout << "      {";
        for (int i = 0; i < pi[target_id].path.size() - 1; i++)
        {
            int x = pi[target_id].path[i], y = pi[target_id].path[i];
            if (pi[target_id].path[i] < pi[target_id].path[i + 1])
            {
                y = pi[target_id].path[i + 1];
                std::cout << "(" << x << ", " << y << ")";
            }
            else
            {
                x = pi[target_id].path[i + 1];
                std::cout << "(" << y << ", " << x << ")";
            }
            if (i < pi[target_id].path.size() - 2)
                std::cout << ", ";
            path->insert({x, y});
        }
        std::cout << "}" << std::endl;
    }
    saveToDot("Dijkstra", path);
    delete path;

    return true;
}

/**
 * @brief Calculates the smallest path between two vertices using Floyd's algorithm, then appends it in a .dot file
 * 
 * @param source_id ID of the starting vertex
 * @param target_id ID of the target vertex
 */
bool Graph::Floyd(int source_id, int target_id)
{
    // Verifies if the vertices exist in the graph (returns false -> not executed)
    if (vertices.count(source_id) == 0 || vertices.count(target_id) == 0)
    {
        std::cerr << "ERROR: Floyd not executed, one or both of the vertices are not present in the graph" << std::endl;
        return false;
    }

    int n = vertices.size();

    Length **A = new Length *[n];
    for (int i = 0; i < n; i++)
    {
        A[i] = new Length[n];
    }

    // Creation of matrix A_0
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i != j)
            {
                if (vertices.count(i) == 0)
                {
                    A[i][j].length = INFINITY;
                }
                else
                {
                    if (vertices[i]->getEdges().count(j) == 0)
                    {
                        A[i][j].length = INFINITY;
                    }
                    else
                    {
                        A[i][j].length = vertices[i]->getEdges()[j]->getWeight();
                        A[i][j].path.push_back(i);
                        A[i][j].path.push_back(j);
                    }
                }
            }
            else
            {
                A[i][j].length = 0;
            }
        }
    }

    // Iterations through A matrix according to K index
    for (int k = 0; k < n; k++)
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                int A_calc;
                if (i != j)
                {
                    if (A[i][k].length != INFINITY && A[k][j].length != INFINITY)
                        A_calc = A[i][k].length + A[k][j].length;
                    else
                        A_calc = INFINITY;
                    if (A[i][j].length > A_calc)
                    {
                        A[i][j].length = A_calc;
                        A[i][j].path = A[i][k].path;
                        for (int l = 1; l < A[k][j].path.size(); l++)
                            A[i][j].path.push_back(A[k][j].path[l]);
                    }
                }
            }
        }
    }

    std::set<std::pair<int, int>> *path = new std::set<std::pair<int, int>>;
    if (A[source_id][target_id].path.size() > 0)
    {
        std::cout << "  Path found by Floyd:" << std::endl;
        std::cout << "      {";
        for (int i = 0; i < A[source_id][target_id].path.size() - 1; i++)
        {
            int x = A[source_id][target_id].path[i], y = A[source_id][target_id].path[i];
            if (A[source_id][target_id].path[i] < A[source_id][target_id].path[i + 1])
            {
                y = A[source_id][target_id].path[i + 1];
                std::cout << "(" << x << ", " << y << ")";
            }
            else
            {
                x = A[source_id][target_id].path[i + 1];
                std::cout << "(" << y << ", " << x << ")";
            }
            if (i < A[source_id][target_id].path.size() - 2)
                std::cout << ", ";
            path->insert({x, y});
        }
        std::cout << "}" << std::endl;
    }

    saveToDot("Floyd", path);
    delete path;

    // Deletion of used dynamic matrixes
    for (int i = 0; i < n; i++)
    {
        delete[] A[i];
    }
    delete[] A;

    return true;
}

struct SimpleEdge
{
    int a;
    int b;
    int weight;

    SimpleEdge(int a, int b, int weight, bool directed)
    {
        if (!directed && b < a)
            std::swap(a, b);

        this->a = a;
        this->b = b;
        this->weight = weight;
    }

    bool operator<(const SimpleEdge &e) const
    {
        return weight < e.weight;
    }

    bool operator==(const SimpleEdge &e) const
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
 */
bool Graph::MST_Kruskal()
{
    if (vertices.empty())
    {
        std::cerr << "ERROR: Kruskal not executed, vertex not found" << std::endl;
        return false;
    }

    std::cout << "  Kruskal's Minimum Spanning Tree:" << std::endl;

    std::unordered_map<int, int> subsets;
    std::list<SimpleEdge> sorted_edges;
    std::set<std::pair<int, int>> *mst_edges = new std::set<std::pair<int, int>>;
    int cost = 0;

    // Create subsets ranging from 0 to n-1
    int i = 0;
    for (std::unordered_map<int, Vertex *>::iterator itV = vertices.begin(); itV != vertices.end(); ++itV)
    {
        subsets.insert({itV->second->getId(), i}); // Put each vertex in its subset
        i++;

        // Traverse the vertex's edges and add them to the list of edges
        std::unordered_map<int, Edge *> edges = itV->second->getEdges();
        for (std::unordered_map<int, Edge *>::iterator itE = edges.begin(); itE != edges.end(); ++itE)
        {
            // Create a "simple edge" containing the origin vertex's id, the target vertex's id, and the edge weight
            sorted_edges.push_back(SimpleEdge(itV->second->getId(), itE->second->getTargetId(), itE->second->getWeight(), this->directed));
        }
    }

    sorted_edges.sort();   // Sort the edges in asceding order according to their weights
    sorted_edges.unique(); // Remove duplicates (in undirected graph)

    std::cout << "    {";
    i = 0;
    while (i < order - 1 && !sorted_edges.empty())
    {
        SimpleEdge e = sorted_edges.front(); // Get first edge in list
        sorted_edges.pop_front();            // Pop edge from the list of edges

        // Get the IDs of the subsets the vertices a and b are in
        int smallest = subsets.at(e.a);
        int greatest = subsets.at(e.b);
        if (smallest != greatest)
        {
            std::cout << "(" << e.a << ", " << e.b << ") "; //print to screen
            mst_edges->insert(std::pair<int, int>(e.a, e.b));
            cost += e.weight; // sum up the edge weight to the MST cost

            if (greatest < smallest)
                std::swap(smallest, greatest);

            // Merge the subsets (changing the greatest id to the smallest id)
            for (std::unordered_map<int, int>::iterator it = subsets.begin(); it != subsets.end(); ++it)
            {
                if (it->second == greatest)
                    it->second = smallest;
            }
            i++;
        }
    }

    std::cout << "}" << std::endl
              << "Cost: " << cost << std::endl;
    saveToDot("Kruskal", mst_edges);

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
 */
bool Graph::BFS(int id)
{
    if (!searchVertex(id))
    {
        std::cerr << "ERROR: BFS not executed, vertex not found!" << std::endl;
        return false;
    }

    std::cout << "  Breadth First Search:" << std::endl;

    // Sets of tree and back edges (which will be printed in different colors in the .dot)
    std::set<std::pair<int, int>> *back_edges = new std::set<std::pair<int, int>>;

    std::unordered_map<int, char> colored_vertices;
    for (std::unordered_map<int, Vertex *>::iterator it = vertices.begin(); it != vertices.end(); ++it)
    {
        colored_vertices.insert({it->second->getId(), 'w'}); // Insert as not visited (white)
    }

    colored_vertices.at(id) = 'g'; // Vertex v is the first to be visited (gray)

    // Queue of vertices to visit
    std::queue<int> toVisit;

    // Enqueue current vertex and mark it as "to visit"
    toVisit.push(id);

    std::cout << "    ";
    while (!toVisit.empty())
    {
        // Dequeue top vertex in queue
        id = toVisit.front();
        std::cout << id << " "; //print to screen
        toVisit.pop();

        // Check for unvisited vertices in adjacency list and enqueue them
        std::unordered_map<int, Edge *> edges = vertices[id]->getEdges();
        for (std::unordered_map<int, Edge *>::iterator itE = edges.begin(); itE != edges.end(); ++itE)
        {
            int target_id = itE->second->getTargetId();

            // If the target vertex is white, then it is unvisited
            if (colored_vertices.at(target_id) == 'w')
            {
                std::pair<int, int> pair_vertices(id, target_id);
                if (!this->directed && pair_vertices.second < pair_vertices.first)
                    std::swap(pair_vertices.first, pair_vertices.second);
                colored_vertices.at(target_id) = 'g'; // mark it as "to visit"
                toVisit.push(target_id);
            }
            else if (colored_vertices.at(target_id) == 'g')
            {
                std::pair<int, int> pair_vertices(id, target_id);
                if (!this->directed && pair_vertices.second < pair_vertices.first)
                    std::swap(pair_vertices.first, pair_vertices.second);
                back_edges->insert(pair_vertices);
            }
        }

        colored_vertices.at(id) = 'b'; // Mark current vertex as visited (black)
    }
    std::cout << std::endl;

    saveToDot("BFS", back_edges);

    delete back_edges;
    return true;
}

void Graph::auxTopologicalSorting(int id, std::map<int, int> &colors, std::list<int> &order, bool *dag)
{

    // Iterator for adjacent list
    std::unordered_map<int, Edge *>::iterator itE;

    // Root vertex becomes gray
    colors.at(id) = 1;

    std::unordered_map<int, Edge *> edges = vertices[id]->getEdges();

    for (itE = edges.begin(); itE != edges.end(); ++itE)
    {

        // gets ID of adjascent vertex
        int target_id = itE->second->getTargetId();

        if (colors.at(target_id) == 0)
        {
            // adjascent vertex unvisited
            auxTopologicalSorting(target_id, colors, order, dag);
        }
        else if (colors.at(target_id) == 1)
        {
            // graph has a return edge so can't have a topological sorting
            // std::cout << "Grafo não é um DAG! (grafo ciclico em (" << id << ", " << target_id << ")" << std::endl;
            *dag = false;
            return;
        }
    }

    // Vertex is terminated -> color becomes black and is added to list
    colors.at(id) = 2;
    order.push_front(id);
}

/**
 * @brief Prints out a topological sorting of the graph
 * 
 
 */
void Graph::topologicalSorting()
{

    bool dag = true;

    std::cout << "Ordenação topológica do grafo:" << std::endl;
    if (!this->directed)
    {
        std::cout << "Grafo não é um DAG! (não direcionado)" << std::endl;
        return;
    }

    // Iterator for list of vertices
    std::unordered_map<int, Vertex *>::iterator it;

    // Map of the vertices Id and a color that indicates wether the vertice is unvisited (0 - white), being visited (1 - gray) or terminated (2 - black)
    std::map<int, int> colors;

    // List of vertices sorted in topological order
    std::list<int> topologicalOrder;

    // Initially all vertices receive color white
    for (it = vertices.begin(); it != vertices.end(); ++it)
    {
        colors.insert({it->second->getId(), 0});
        /* std::cout << colors.at(it->second->getId()); */
    }

    // Loop through all unvisited vertices
    for (it = vertices.begin(); it != vertices.end(); ++it)
    {
        if (colors.at(it->second->getId()) == 0)
        {
            /* std::cout << it->second->getId() << " eh branco" << std::endl; */
            auxTopologicalSorting(it->second->getId(), colors, topologicalOrder, &dag);
        }
        if (!dag)
        {
            std::cout << "Grafo não é um DAG!" << std::endl;
            return;
        }
    }

    while (!topologicalOrder.empty())
    {
        std::cout << topologicalOrder.front() << " ";
        topologicalOrder.pop_front();
    }

    std::cout << std::endl;
}