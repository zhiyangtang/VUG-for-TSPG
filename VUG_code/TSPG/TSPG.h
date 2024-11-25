#ifndef TSPG_H
#define TSPG_H
#include "../GraphUtils/Graph.cc"

class TSPG {
    public:
    TSPG(Graph* inputGraph);
    EdgeID executeQuery(VertexID source, VertexID target, short begin, short end);      // execute for each query
    void answerAllQueries(vector<PerQuery>& queries);                                   // experiments for answering all queries
    void cleanUp();                                                                     // free memory after running all queries

    private:
    // basic graph information
    Graph* graph;
    VertexID VN;                                              // |V| of graph
    EdgeID EN;                                                // |E| of graph
    vector<PerEdge> edges;                                    // store all edges in graph
    PerNeighbor *inNeighbors_time, *outNeighbors_time;        // neighbors of each vertex, length=EN
    EdgeID *inNeighborsLocator, *outNeighborsLocator;         // locate where to find the neighbors of a vertex, length=VN

    double T1 = 0, T2 = 0, T3 = 0;
    EdgeID max = 0;

    // initialize and refresh memory for queries
    void init(); 
    inline void refreshMemory();

    // adaptive time-bi-directional BFS
    VertexID s, t;
    Time t_begin, t_end;
    VertexID *dropflag;
    EdgeID *isvisitedEdge, tpgEnd;
    PerEdge *tpg_graph;
    void earliestTimeBFS();
    void latestTimeBFS();
    void prunedVertices();
    void prunedEdges();

    // essential vertices
    Essential *forwardEV, *backwardEV;
    VertexID forwardEVEnd, backwardEVEnd;
    EdgeID *forwardLocator, *new_forwardLocator, *backwardLocator, *new_backwardLocator, upperEnd;
    int *forwardFlag, *forwardEVLen, *backwardFlag, *backwardEVLen;
    PerEdge *upper_bound_graph, *new_sortOutEdges, *new_sortInEdges;
    EdgeID *edgesForVerification, edgesForVerificationEnd;
    VertexID forwardLocatorEnd, new_forwardLocatorEnd, backwardLocatorEnd, new_backwardLocatorEnd;
    void forward();
    void backward();
    void UpperBound();

    // exact graph
    VertexID *departures, *arrivals;
    VertexID departure, arrival, departuresEnd, arrivalsEnd;
    Time departuretime, arrivaltime;
    VertexID maxStack;
    EdgeID *hasVerified;
    VertexID *isinStack, *stack, stackEnd;
    EdgeID *stacke, stackeEnd;
    void verifyUndeterminedEdge();
    void DFS(VertexID &u, Time &time);
    void addEdge();
    bool tryAddEdges(VertexID& d, VertexID& a, Time& dt, Time& at);
    bool forwardSearch(VertexID& from, VertexID& to, EdgeID& e, Time& time, Time& fromtime);
    bool backwardSearch(VertexID& from, EdgeID& e, Time& time);
    bool backward2Search(VertexID& from, VertexID& to, EdgeID& e, Time& time, Time& totime);
    bool forward2Search(VertexID& to, EdgeID& e, Time& time);
    void addToResult();

    EdgeID *results, resultEnd;

    // write statistics of queries to file    
    #ifdef WRITE_STATISTICS
        void initStatisticStorage(int queryNumber);
        double getCurrentSpaceCost();
        void cleanUpStatisticStorage();  
        int queryId = 0;
        EdgeID *numOfUpperbound, *numOfAnswers, *numOfTPAnswers;
        VertexID maxFrontierSize;        
        double *spaceCosts;  
    #endif   
};


#endif