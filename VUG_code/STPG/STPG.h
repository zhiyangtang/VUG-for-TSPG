#ifndef STPG_H
#define STPG_H
#include "../GraphUtils/Graph.cc"

class STPG {
    public:
    STPG(Graph* inputGraph);
    EdgeID executeQuery(VertexID source, VertexID target, Time begin, Time end);      // execute for each query
    void answerAllQueries(vector<PerQuery>& queries);                                 // experiments for answering all queries
    void cleanUp();                                                                   // free memory after running all queries

    private:
    Graph* graph;
    VertexID VN;                                              // |V| of graph
    EdgeID EN;                                                // |E| of graph
    vector<PerEdge> edges;                                    // store all edges in graph
    PerNeighbor *inNeighbors_time, *outNeighbors_time;        // neighbors of each vertex, length=EN
    EdgeID *inNeighborsLocator, *outNeighborsLocator;         // locate where to find the neighbors of a vertex, length=VN

    void init(); 
    inline void refreshMemory();

    VertexID s, t;
    Time t_begin, t_end;
    VertexID *dropflag;
    EdgeID *isvisitedEdge, tpgEnd;
    PerEdge *tpg_graph;
    void earliestTimeBFS();
    void latestTimeBFS();
    void prunedVertices();
    void prunedEdges();

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

    VertexID *departures, *arrivals;
    VertexID maxStack;
    EdgeID *hasVerified;
    VertexID *isinStack, *stack, stackEnd;
    EdgeID *stacke, stackeEnd;
    void verifyUndeterminedEdge();
    bool forwardSearch(VertexID& from, VertexID& to, EdgeID& e, Time& time, Time& fromtime);
    bool backwardSearch(VertexID& from, EdgeID& e, Time& time);
    bool backward2Search(VertexID& from, VertexID& to, EdgeID& e, Time& time, Time& totime);
    bool forward2Search(VertexID& to, EdgeID& e, Time& time);
    void addToResult(VertexID& from, VertexID& to, Time& time);

    EdgeID *results, resultEnd;
 
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