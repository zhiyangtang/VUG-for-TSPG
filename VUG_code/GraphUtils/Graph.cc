#ifndef GRAPH_CC
#define GRAPH_CC
#include "Graph.h"
using namespace std;

Graph::Graph(const char* inputGraphFilename) {
    graphFilename = inputGraphFilename;
    loadGraphFile();
}


void Graph::loadGraphFile() {

    double startTime = getCurrentTimeInMs();
    printf("Loading graph file: %s ...\n", graphFilename);

    VertexID n = 0;
    EdgeID m = 0;

    ifstream ifs(graphFilename);
    for (string line; getline(ifs, line);) {
        istringstream iss(line);
        unsigned type;
        PerEdge e;
        e.edgeId = m;

        iss >> e.fromId >> e.toId >> e.time; 

        if (m == 0) {
            min_time = e.time;
            max_time = e.time;
        }
        else {
            if (e.time < min_time)
                min_time = e.time;
            if (e.time > max_time)
                max_time = e.time;
        }

        if (e.fromId >= n || e.toId >= n)
        {
            n = max(e.fromId, e.toId) + 1;
        }

        edges.push_back(e);
        m++;
    }

    VN = n;
    EN = m;

    max_time = (max_time - min_time) / 86400 + 2;   // from UNIX timestamp to day

    inNeighbors_time = new PerNeighbor[EN];
    inNeighborsLocator = new EdgeID[VN+1];
    outNeighbors_time = new PerNeighbor[EN];
    outNeighborsLocator = new EdgeID[VN+1];

    sort(edges.begin(), edges.end(), sortByToIdTime);
    EdgeID curLocator = 0;
    VertexID curId = 0, id;
    EdgeID i = 0;
    
    while (i<EN) {
        id = edges[i].toId;
        while (curId<=id) {
            inNeighborsLocator[curId] = curLocator;
            curId++;
        }
        while (i<EN && edges[i].toId==id) {
            edges[i].time = (edges[i].time - min_time) / 86400 + 1;
            inNeighbors_time[curLocator] = {edges[i].edgeId, edges[i].fromId, edges[i].time};
            curLocator++;
            i++;
        }
    }
    while (curId<=VN) {
        inNeighborsLocator[curId] = curLocator;
        curId++;
    }

    sort(edges.begin(), edges.end(), sortByFromIdTime);
    curLocator = 0;
    curId = 0;
    i = 0;
    while (i<EN) {
        id = edges[i].fromId;
        while (curId<=id) {
            outNeighborsLocator[curId] = curLocator;
            curId++;
        }
        while (i<EN && edges[i].fromId==id) {
            outNeighbors_time[curLocator] = {edges[i].edgeId, edges[i].toId, edges[i].time};
            curLocator++;
            i++;
        }
    }
    while (curId<=VN) {
        outNeighborsLocator[curId] = curLocator;
        curId++;
    }
    sort(edges.begin(), edges.end(), sortByEdgeID); 
    
    double timeCost = getCurrentTimeInMs() - startTime;

    printf("- Finish. |V|=%d and |E|=%d, time cost: %.2f ms\n", VN, EN, timeCost);
    logFile << VN << "," << EN << "," << str(timeCost) << ",";
}

#endif