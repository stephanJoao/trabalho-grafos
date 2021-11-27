#include "../include/Vertex.hpp"

//* Constructors and destructors implementations

Vertex::Vertex(int id, float weight)
{
    this->id = id;
    this->first_edge = nullptr;
    this->out_degree = 0;
    this->in_degree = 0;
    this->weight = weight;
    this->next_vertex = nullptr;
}

Vertex::~Vertex()
{
    Edge* next_edge = this->first_edge;

    while(next_edge != nullptr)
    {
        Edge* aux_edge = next_edge->getNextEdge();
        delete next_edge;
        next_edge = aux_edge;
    }
}

//* Getters and setters implementations

Edge* Vertex::getFirstEdge()
{
    return this->first_edge;
}

int Vertex::getId()
{
    return this->id;
}

int Vertex::getInDegree()
{
    return this->in_degree;
}

int Vertex::getOutDegree()
{
    return this->out_degree;
}

float Vertex::getWeight()
{
    return this->weight;
}

void Vertex::setWeight(float weight)
{
    this->weight = weight;
}

Vertex* Vertex::getNextVertex()
{
    return this->next_vertex;
}
    
void Vertex::setNextVertex(Vertex* vertex)
{
    this->next_vertex = vertex;
}

//* Other methods implementations

/**
 * Search edge in vertex.
 *
 * @param target Target vertex id
 */
bool Vertex::searchEdge(int target_id){
     // If there is at least one edge in the vertex
    if(this->first_edge != nullptr){
        // Searching for a specific edge of target id equal to target_id
        for(Edge* aux = this->first_edge; aux != nullptr; aux = aux->getNextEdge())
            if(aux->getTargetId() == target_id)
                return true;
    }

    return false;
}

/**
 * Insert edge in vertex.
 *
 * @param target Target vertex id
 * @param weight Weight of the edge (defaults to 1 if not given)
 */
void Vertex::insertEdge(int target_id, float weight) 
{
    // If there is at least one edge in vertex
    if (this->first_edge != nullptr){
        // Chain the new edge after the last edge
        Edge* new_edge = new Edge(target_id, nullptr, weight);
        this->last_edge->setNextEdge(new_edge);
        this->last_edge = new_edge;
    } else {
        // If there is no edge, just add the new edge
        this->first_edge = new Edge(target_id, nullptr, weight);
        this->last_edge = this->first_edge;
    }
}

/**
 * Remove all edges from vertex.
 */
void Vertex::removeAllEdges()
{
    // If there is at least one edge in vertex
    if(this->first_edge != nullptr){
        Edge* next = nullptr;
        Edge* aux = this->first_edge;
        // Removing all edges
        while(aux != nullptr){
            next = aux->getNextEdge();
            delete aux;
            aux = next;
        }
    }

    this->first_edge = nullptr;
    this->last_edge = nullptr;
}

/**
 * Remove edge from vertex.
 *
 * @param id Edge id
 * @param directed If the graph is directed
 * @param target_vertex Target vertex
 */
int Vertex::removeEdge(int id, bool directed, Vertex* target_vertex){
    // Check if edge exists in vertex
    if(this->searchEdge(id)){
        Edge* aux = this->first_edge;
        Edge* previous = nullptr;
        // Searching for the edge to be removed
        while(aux->getTargetId() != id){
            previous = aux;
            aux = aux->getNextEdge();
        }
        // Keeping the integrity of the edge list
        if(previous != nullptr)
            previous->setNextEdge(aux->getNextEdge());
        else
            this->first_edge = aux->getNextEdge();

        if(aux == this->last_edge)
            this->last_edge = previous;

        if(aux->getNextEdge() == this->last_edge)
            this->last_edge = aux->getNextEdge();

        delete aux;
        // Decrement out degree if the graph is directed
        if(directed){
            this->decrementOutDegree();
        } else {
            this->decrementInDegree();
            target_vertex->decrementInDegree();
        }

        return 1;
    }

    return 0;
}

void Vertex::incrementInDegree()
{
    this->in_degree++;
}

void Vertex::incrementOutDegree()
{
    this->out_degree++;
}

void Vertex::decrementInDegree()
{
    this->in_degree--;
}

void Vertex::decrementOutDegree()
{
    this->out_degree--;
}

// Edge* Vertex::hasEdgeBetween(int target_id)
// {

//     for(Edge *auxEdge = this->first_edge; auxEdge != nullptr; auxEdge = auxEdge->getNextEdge())
//     {
//         if(auxEdge->getTargetId() == target_id)
//             return auxEdge;
//     }
//     return nullptr;
// }
