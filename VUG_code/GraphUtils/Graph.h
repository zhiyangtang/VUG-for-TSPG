#ifndef GRAPH_H
#define GRAPH_H
#include "Utils.h"

struct PerEdge {
    EdgeID edgeId;
    VertexID fromId, toId;
    Time time;
};


struct PerNeighbor {
    EdgeID edgeId;
    VertexID neighbor;
    Time time;
    PerNeighbor(): edgeId(0), neighbor(0), time(0){}
    PerNeighbor(EdgeID a, VertexID b, Time c): edgeId(a), neighbor(b), time(c){}
    bool operator<(const PerNeighbor& a)const {
        return neighbor<a.neighbor;
    }
    bool operator<<(const PerNeighbor& a)const {
        return time<a.time;
    }
};


struct Essential {
    Time time;
    int size;
    VertexID* vertices = nullptr;
};


class Graph {

    public:

        Graph(const char* inputGraphFilename);
        VertexID VN;                                                // |V| of graph
        EdgeID EN;                                                  // |E| of graph
        Time min_time, max_time;                                                  
        vector<PerEdge> edges;                                      // store all edges in graph
        PerNeighbor *inNeighbors_time, *outNeighbors_time;          // neighbors of each vertex, length=EN
        EdgeID *inNeighborsLocator, *outNeighborsLocator;           // locate where to find the neighbors of a vertex, length=VN
        /*
        Compressed Sparse Row (CSR)
        If V={0,1,2,3} and E={0->1, 0->2, 0->3, 1->2, 1->3, 2->3} then:
        outNeighbors = [1,2,3,2,3,3] and outNeighborsLocator = [0,3,5],
        which means the out-neighbors of vertex v are: outNeighbors[outNeighborsLocator[v]:outNeighborsLocator[v+1]].
        For example, the out-neighbors of vertex 0 are: outNeighbor[0:3], i.e., [1,2,3]
        */

    private:
        const char* graphFilename;
        void loadGraphFile();

};


bool sortByFromIdTime(PerEdge& a, PerEdge& b) {
    if (a.fromId==b.fromId) {
        if (a.time==b.time)
            return a.toId<b.toId;
        return a.time<b.time;
    }
    return a.fromId<b.fromId;
}


bool sortFromIdTime(PerEdge& a, PerEdge& b) {
    if (a.fromId==b.fromId) {
        if (a.time==b.time)
            return a.toId<b.toId;
        return a.time>b.time;
    }
    return a.fromId<b.fromId;
}


bool sortByToIdTime(PerEdge& a, PerEdge& b) {
    if (a.toId==b.toId) {
        if (a.time==b.time)
            return a.fromId<b.fromId;
        return a.time<b.time;
    }
    return a.toId<b.toId;
}


bool sortToIdTime(PerEdge& a, PerEdge& b) {
    if (a.toId==b.toId) {
        if (a.time==b.time)
            return a.fromId<b.fromId;
        return a.time>b.time;
    }
    return a.toId<b.toId;
}


bool sortByEdgeTime(PerEdge& a, PerEdge& b) {
    if (a.time == b.time) {
        if (a.fromId == b.fromId)
            return a.toId<b.toId;
        return a.fromId<b.fromId;
    }
    return a.time<b.time;
}


bool sortByEdgeID(PerEdge& a, PerEdge& b) {
    return a.edgeId < b.edgeId;
}

#endif