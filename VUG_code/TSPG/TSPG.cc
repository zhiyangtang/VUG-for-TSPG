#ifndef TSPG_CC
#define TSPG_CC
#include "TSPG.h"
using namespace std;

TSPG::TSPG(Graph* inputGraph) {
    // basic graph information
    graph = inputGraph;
    VN = graph->VN;
    EN = graph->EN;
    edges = graph->edges;
    inNeighborsLocator = graph->inNeighborsLocator;
    inNeighbors_time = graph->inNeighbors_time;
    outNeighborsLocator = graph->outNeighborsLocator;
    outNeighbors_time = graph->outNeighbors_time;
    
    // initialization
    init();
}

// initialization
void TSPG::init() {
    
    nextFrontier = new VertexID[VN];
    forwardFrontier = new VertexID[VN];
    backwardFrontier = new VertexID[VN];
    
    earliestTime = new Time[VN]();
    latestTime = new Time[VN]();
    isvisitedEdge = new EdgeID[EN];
    
    dropflag = new VertexID[VN]();
    tpg_graph = new PerEdge[EN];
    upper_bound_graph = new PerEdge[EN];

    forwardLocator = new EdgeID[VN + 1];
    forwardFlag = new int[VN]();
    forwardEVLen = new int[VN]();

    backwardLocator = new EdgeID[VN + 1];
    backwardFlag = new int[VN]();
    backwardEVLen = new int[VN]();

    new_forwardLocator = new EdgeID[VN + 1];
    new_backwardLocator = new EdgeID[VN + 1];

    edgesForVerification = new EdgeID[EN];
    departures = new VertexID[VN]();
    arrivals = new VertexID[VN]();

    hasVerified = new EdgeID[EN]();
    stack = new VertexID[VN];
    stacke = new EdgeID[VN];
    isinStack = new VertexID[VN]();

    results = new EdgeID[EN];
}

// experiments for answering all queries
void TSPG::answerAllQueries(vector<PerQuery>& queries) {
    if (queries.size()==0) {
        printf("!No query\n");
        return;
    }
    
    // answers and statistics file
    #ifdef WRITE_ANSWERS
        resultFile.open("../"+answerPath+queryFilename+".answer",ios::app);
        resultFile << "number of edges,edge ids" << endl;
    #endif
    #ifdef WRITE_STATISTICS
        queryId = 0;
        initStatisticStorage(queries.size());
    #endif

    // initialization
    printf("Running STPG ...\n");
    double startTime = getCurrentTimeInMs();

    // double ratio=0, ratio1=0;

    for (PerQuery& query : queries) {

        resultEnd = executeQuery(query.source, query.target, query.star, query.end);
        
        // sort and write results to file
        #ifdef WRITE_ANSWERS
            std::sort(results, results + resultEnd);
            resultFile << resultEnd;
            for (EdgeID j = 0; j < resultEnd; j++)
                resultFile << "," << results[j];
            resultFile << endl;
        #endif
        #ifdef WRITE_STATISTICS
            // if (numOfTPAnswers[queryId] == 0)
            //     ratio += 0;
            // else
            //     ratio+=(numOfAnswers[queryId]/numOfTPAnswers[queryId]);
            // if (numOfUpperbound[queryId] == 0)
            //     ratio1 += 0;
            // else
            //     ratio1+=(numOfAnswers[queryId]/numOfUpperbound[queryId]);
            statisticsFile << spaceCosts[queryId] << "," << numOfTPAnswers[queryId] << "," << numOfUpperbound[queryId] << "," << numOfAnswers[queryId] << endl;
            queryId++;
        #endif
    }

    // output logs
    double timeCost = getCurrentTimeInMs() - startTime;
    printf("- Finish. Time cost: %.2f ms\n", timeCost);
    logFile << str(timeCost) << endl;

    // cout << ratio << endl;
    // cout << ratio1 << endl;
    //cout << T1 << endl;
    // cout << T1+T2 << endl;
    //cout << T3 << endl;
    // cout << max << endl;

    // output answers and statistics file
    #ifdef WRITE_STATISTICS
        cleanUpStatisticStorage();
    #endif
    #ifdef WRITE_ANSWERS
        resultFile.close();
    #endif
}

// execute for each query
EdgeID TSPG::executeQuery(VertexID source, VertexID target, short begin, short end) {
    
    // initialization
    s = source;
    t = target;
    t_begin = begin;
    t_end = end;
    resultEnd = 0;

    double startTime1 = getCurrentTimeInMs();
    // earliest and latest time for TPG
    earliestTimeBFS();
    latestTimeBFS();
    prunedVertices();
    prunedEdges();
    T1 += getCurrentTimeInMs() - startTime1;

    double startTime2 = getCurrentTimeInMs();
    // essential vertices for upper_bound graph of STPG
    forward();
    backward();
    UpperBound();
    T2 += getCurrentTimeInMs() - startTime2;

    //double startTime3 = getCurrentTimeInMs();
    // exact graph of STPG
    if (edgesForVerificationEnd) {
        verifyUndeterminedEdge();
    }
    // stackEnd = 0;
    // stackeEnd = 0;
    // maxStack = 0;
    // resultEnd = 0;
    // DFS(s,t_begin);
    //T3 += getCurrentTimeInMs() - startTime3;
    
    refreshMemory();

    // statistics
    #ifdef WRITE_STATISTICS
        numOfAnswers[queryId] = resultEnd;
        spaceCosts[queryId] = getCurrentSpaceCost();
    #endif

    return resultEnd;
}

void TSPG::earliestTimeBFS() {
    earliestTime[s] = offset;
    earliestTime[t] = offset;
    forwardFrontier[0] = s;
    forwardFrontierEnd = 1;
    memset(isvisitedEdge, 0, sizeof(EdgeID)*EN);

    while (forwardFrontierEnd > 0) {
        nextFrontierEnd = 0;

        for (VertexID i = 0; i < forwardFrontierEnd; i++) {
            VertexID& u = forwardFrontier[i];

            // iterate each edge u->v
            for (EdgeID edgeLocator = outNeighborsLocator[u]; edgeLocator < outNeighborsLocator[u + 1]; edgeLocator++) {
                Time& time = outNeighbors_time[edgeLocator].time;
                VertexID& v = outNeighbors_time[edgeLocator].neighbor;
                EdgeID& edgeID = outNeighbors_time[edgeLocator].edgeId;

                if (isvisitedEdge[edgeID] == 1 || time < t_begin || v == t || v == s) continue;   
                else if (time > t_end) break;
                else {
                    if (time + offset > earliestTime[u]) {
                        isvisitedEdge[edgeID] = 1;
                        // earliestTime[v] does not exist and update earliestTime[v]
                        if (earliestTime[v] <= offset) {
                            earliestTime[v] = time + offset;
                            // push to next frontier
                            if (time != t_end) {
                                nextFrontier[nextFrontierEnd] = v;
                                nextFrontierEnd++;
                            }
                        }
                        // earliestTime[v] exists
                        else {
                            // update to an earlier time and push to next frontier
                            if (time + offset < earliestTime[v]) {
                                earliestTime[v] = time + offset;
                                VertexID* re = find(nextFrontier, nextFrontier + nextFrontierEnd, v);
                                if (re == nextFrontier + nextFrontierEnd) {
                                    nextFrontier[nextFrontierEnd] = v;
                                    nextFrontierEnd++;
                                }
                            }
                        }
                    }
                }
            }
        }

        // swap frontier
        VertexID* tmp = forwardFrontier;
        forwardFrontier = nextFrontier;
        nextFrontier = tmp;
        forwardFrontierEnd = nextFrontierEnd;

        // statistics
        #ifdef WRITE_STATISTICS
            maxFrontierSize = max(forwardFrontierEnd, maxFrontierSize);
        #endif
    }
}

void TSPG::latestTimeBFS() {
    latestTime[s] = offset + t_end + 1;
    latestTime[t] = offset + t_end + 1;
    backwardFrontier[0] = t;
    backwardFrontierEnd = 1;
    memset(isvisitedEdge, 0, sizeof(EdgeID)*EN);

    while(backwardFrontierEnd > 0) {
        nextFrontierEnd = 0;

        for (VertexID i = 0; i < backwardFrontierEnd; i++) {
            VertexID& u = backwardFrontier[i];

            // iterate each edge u->v
            for (EdgeID edgeLocator = inNeighborsLocator[u]; edgeLocator < inNeighborsLocator[u + 1]; edgeLocator++) {
                Time& time = inNeighbors_time[edgeLocator].time;
                VertexID& v = inNeighbors_time[edgeLocator].neighbor;
                EdgeID& edgeID = outNeighbors_time[edgeLocator].edgeId;

                if (isvisitedEdge[edgeID] == 1 || time < t_begin || v == t || v == s) continue;   
                else if (time > t_end) break;
                else {
                    if (time + offset < latestTime[u]) {
                        isvisitedEdge[edgeID] = 1;
                        // earliestTime[v] does not exist and update earliestTime[v]
                        if (latestTime[v] <= offset) {
                            latestTime[v] = time + offset;
                            // push to next frontier
                            if (time != t_begin) {
                                nextFrontier[nextFrontierEnd] = v;
                                nextFrontierEnd++;
                            }
                        }
                        // earliestTime[v] exists
                        else {
                            // update to an earlier time and push to next frontier
                            if (time + offset > latestTime[v]) {
                                latestTime[v] = time + offset;
                                VertexID* re = find(nextFrontier, nextFrontier + nextFrontierEnd, v);
                                if (re == nextFrontier + nextFrontierEnd) {
                                    nextFrontier[nextFrontierEnd] = v;
                                    nextFrontierEnd++;
                                }
                            }
                        }
                    }
                }
            }
        }

        // swap frontier
        VertexID* tmp = backwardFrontier;
        backwardFrontier = nextFrontier;
        nextFrontier = tmp;
        backwardFrontierEnd = nextFrontierEnd;

        // statistics
        #ifdef WRITE_STATISTICS
            maxFrontierSize = max(backwardFrontierEnd, maxFrontierSize);
        #endif
    }
}

void TSPG::prunedVertices() {
    for (int i = 0; i < VN; ++i) {
        if (i != s && i != t) {
            if (earliestTime[i] <= offset || latestTime[i] <= offset || earliestTime[i] >= latestTime[i]) {
                dropflag[i] = offset + 1;
            }
        }
    }
}

void TSPG::prunedEdges() {
    tpgEnd = 0;
    for (int i = 0; i < VN; ++i) {
        if (dropflag[i] != offset + 1) {
            for (EdgeID edgeLocator = outNeighborsLocator[i]; edgeLocator < outNeighborsLocator[i + 1]; edgeLocator++) {
                EdgeID& id = outNeighbors_time[edgeLocator].edgeId;
                VertexID& v = outNeighbors_time[edgeLocator].neighbor;
                Time& time = outNeighbors_time[edgeLocator].time;

                if (dropflag[v] != offset + 1 && v != s && i != t && time >= t_begin && time <= t_end && time + offset > earliestTime[i] && time + offset < latestTime[v]) {
                    tpg_graph[tpgEnd++] = edges[id];
                }
            }
        }
    }

    std::sort(tpg_graph, tpg_graph + tpgEnd, sortByEdgeTime);

    // statistics
    #ifdef WRITE_STATISTICS
        numOfTPAnswers[queryId] = tpgEnd;
    #endif
}

void TSPG::forward() {
    
    PerEdge* sortInEdges = new PerEdge[tpgEnd];
    std::copy(tpg_graph, tpg_graph + tpgEnd, sortInEdges);
    std::sort(sortInEdges, sortInEdges + tpgEnd, sortByToIdTime);

    forwardEVEnd = 0;
    forwardLocatorEnd = 0;
    forwardEV = new Essential[tpgEnd];
    forwardEVEnd += tpgEnd;

    EdgeID curLocator = 0;
    VertexID curId = 0, id;
    int i = 0;

    while (i < tpgEnd) {
        id = sortInEdges[i].toId;
        while (curId<=id) {
            forwardLocator[curId] = curLocator;
            curId++;
        }
        forwardLocatorEnd++;
        while (i < tpgEnd && sortInEdges[i].toId==id) {
            forwardEV[curLocator].time = sortInEdges[i].time;
            curLocator++;
            i++;
        }
    }
    while (curId <= VN) {
        forwardLocator[curId] = curLocator;
        curId++;
    }

    delete[] sortInEdges;

    for (int i = 0; i < tpgEnd; ++i) {
        Time& time = tpg_graph[i].time;
        VertexID& from = tpg_graph[i].fromId;
        VertexID& to = tpg_graph[i].toId;

        if (forwardFlag[to] != offset + 2 && to != t) {
            // one hop
            if (from == s) {
                forwardFlag[to] = offset + 2;
                if (forwardEV[forwardLocator[to]].vertices == nullptr) {
                    forwardEV[forwardLocator[to]].vertices = new VertexID(to);
                    forwardEV[forwardLocator[to]].size = 1;
                    forwardEVLen[to] = forwardLocator[to] + offset;
                    forwardEVEnd += 1;
                }

                else {
                    forwardEVLen[to] += 1;
                    forwardEV[forwardEVLen[to] - offset].vertices = new VertexID(to);
                    forwardEV[forwardEVLen[to] - offset].size = 1;
                    forwardEVEnd += 1;
                }
            }

            else {
                if (forwardFlag[to] != offset + 1) {
                    int locator = forwardEVLen[from] - offset;
                    while (time == forwardEV[locator].time)
                        locator--;
                    int size = forwardEV[locator].size;
                    int newsize = size + 1;
                    VertexID* ev = new VertexID[newsize];
                    std::copy(forwardEV[locator].vertices, forwardEV[locator].vertices + size, ev);
                    ev[size] = to;
                    std::sort(ev, ev + newsize);
                    
                    if (forwardEV[forwardLocator[to]].vertices == nullptr) {
                        forwardEV[forwardLocator[to]].vertices = ev;
                        forwardEV[forwardLocator[to]].size = newsize;
                        forwardEVLen[to] = forwardLocator[to] + offset;
                        forwardEVEnd += newsize;
                    }

                    else {
                        int forsize = forwardEV[forwardEVLen[to] - offset].size;
                        VertexID* tmp = new VertexID[min(newsize, forsize)];
                        VertexID* end = std::set_intersection(ev, ev + newsize, forwardEV[forwardEVLen[to] - offset].vertices, forwardEV[forwardEVLen[to] - offset].vertices + forsize, tmp);
                        VertexID* result = new VertexID[end - tmp];
                        std::copy(tmp, end, result);

                        forwardEVLen[to] += 1;
                        forwardEV[forwardEVLen[to] - offset].vertices = result;
                        forwardEV[forwardEVLen[to] - offset].size = end - tmp;
                        forwardEVEnd = forwardEVEnd + end - tmp;

                        if (end - tmp == 1)
                            forwardFlag[to] = offset + 1;
                        
                        delete[] tmp;
                        delete[] ev;
                    }
                }
            }
        }
    }
}

void TSPG::backward() {

    PerEdge* sortOutEdges = new PerEdge[tpgEnd];
    std::copy(tpg_graph, tpg_graph + tpgEnd, sortOutEdges);
    std::sort(sortOutEdges, sortOutEdges + tpgEnd, sortFromIdTime);

    backwardEVEnd = 0;
    backwardLocatorEnd = 0;
    backwardEV = new Essential[tpgEnd];
    backwardEVEnd += tpgEnd;

    EdgeID curLocator = 0;
    VertexID curId = 0, id;
    int i = 0;

    while (i < tpgEnd) {
        id = sortOutEdges[i].fromId;
        while (curId <= id) {
            backwardLocator[curId] = curLocator;
            curId++;
        }
        backwardLocatorEnd++;
        while (i < tpgEnd && sortOutEdges[i].fromId == id) {
            backwardEV[curLocator].time = sortOutEdges[i].time;
            curLocator++;
            i++;
        }
    }
    while (curId <= VN) {
        backwardLocator[curId] = curLocator;
        curId++;
    }

    delete[] sortOutEdges;

    for (int i = tpgEnd - 1; i >= 0; --i) {

        Time& time = tpg_graph[i].time;
        VertexID& from = tpg_graph[i].fromId;
        VertexID& to = tpg_graph[i].toId;

        if (backwardFlag[from] != offset + 2 && from != s) {
            // one hop
            if (to == t) {
                backwardFlag[from] = offset + 2;

                if (backwardEV[backwardLocator[from]].vertices == nullptr) {
                    backwardEV[backwardLocator[from]].vertices = new VertexID(from);
                    backwardEVLen[from] = backwardLocator[from] + offset;
                    backwardEV[backwardLocator[from]].size = 1;
                    backwardEVEnd++;
                }

                else {
                    backwardEVLen[from] += 1;
                    backwardEV[backwardEVLen[from] - offset].vertices = new VertexID(from);
                    backwardEV[backwardEVLen[from] - offset].size = 1;
                    backwardEVEnd++;
                }
            }

            else {
                if (backwardFlag[from] != offset + 1) {
                    int locator = backwardEVLen[to] - offset;
                    while (time == backwardEV[locator].time)
                        locator--;

                    int size = backwardEV[locator].size;
                    int newsize = size + 1;
                    VertexID* ev = new VertexID[newsize];
                    std::copy(backwardEV[locator].vertices, backwardEV[locator].vertices + size, ev);
                    ev[size] = from;
                    std::sort(ev, ev + newsize);

                    if (backwardEV[backwardLocator[from]].vertices == nullptr) {
                        backwardEV[backwardLocator[from]].vertices = ev;
                        backwardEVLen[from] = backwardLocator[from] + offset;
                        backwardEV[backwardLocator[from]].size = newsize;
                        backwardEVEnd += newsize;
                    }

                    else {
                        int forsize = backwardEV[backwardEVLen[from] - offset].size;
                        VertexID* tmp = new VertexID[min(newsize, forsize)];
                        VertexID* end = std::set_intersection(ev, ev + newsize, backwardEV[backwardEVLen[from] - offset].vertices, backwardEV[backwardEVLen[from] - offset].vertices + forsize, tmp);
                        VertexID* result = new VertexID[end - tmp];
                        std::copy(tmp, end, result);

                        backwardEVLen[from] += 1;
                        backwardEV[backwardEVLen[from] - offset].vertices = result;
                        backwardEV[backwardEVLen[from] - offset].size = end - tmp;
                        backwardEVEnd = backwardEVEnd + end - tmp;

                        if (end - tmp == 1) 
                            backwardFlag[from] = offset + 1;
                        
                        delete[] tmp;
                        delete[] ev;
                    }
                }
            }
        }

    }
}

void TSPG::UpperBound() {
    
    upperEnd = 0;
    edgesForVerificationEnd = 0;
    departuresEnd = 0, arrivalsEnd = 0;
    new_forwardLocatorEnd = 0, new_backwardLocatorEnd = 0;
    int forsize, backsize;

    for (int i = 0; i < tpgEnd; ++i) {
        
        if (tpg_graph[i].fromId == s || tpg_graph[i].toId == t) {
            results[resultEnd++] = tpg_graph[i].edgeId;
            upper_bound_graph[upperEnd++] = edges[tpg_graph[i].edgeId];
            hasVerified[tpg_graph[i].edgeId] = offset + 1; 
        }
        else {
            VertexID *forward, *backward;
            for (EdgeID edgeLocator = forwardLocator[tpg_graph[i].fromId + 1] - 1; edgeLocator >= forwardLocator[tpg_graph[i].fromId]; --edgeLocator) {
                if (forwardEV[edgeLocator].time < tpg_graph[i].time) {
                    if (forwardEV[edgeLocator].vertices != nullptr) {
                        forward = forwardEV[edgeLocator].vertices;
                        forsize = forwardEV[edgeLocator].size;
                        break;
                    }
                }
            }
            for (EdgeID edgeLocator = backwardLocator[tpg_graph[i].toId + 1] - 1; edgeLocator >= backwardLocator[tpg_graph[i].toId]; --edgeLocator) {
                if (backwardEV[edgeLocator].time > tpg_graph[i].time) {
                    if (backwardEV[edgeLocator].vertices != nullptr) {
                        backward = backwardEV[edgeLocator].vertices;
                        backsize = backwardEV[edgeLocator].size;
                        break;
                    }
                }
            }

            VertexID *tmp = new VertexID[min(forsize, backsize)];
            VertexID *end = std::set_intersection(forward, forward + forsize, backward, backward + backsize, tmp);

            if (end - tmp == 0) {
                upper_bound_graph[upperEnd++] = tpg_graph[i];
                
                if (forwardFlag[tpg_graph[i].fromId] == offset + 2 && backwardFlag[tpg_graph[i].toId] == offset + 2 && forsize == 1 && backsize == 1 && forwardEV[forwardLocator[tpg_graph[i].fromId]].time < tpg_graph[i].time && backwardEV[backwardLocator[tpg_graph[i].toId]].time > tpg_graph[i].time) {
                    results[resultEnd++] = tpg_graph[i].edgeId;
                    hasVerified[tpg_graph[i].edgeId] = offset + 1;
                    if (departures[tpg_graph[i].toId] <= offset) {
                        departures[tpg_graph[i].toId] = offset + 1;
                        departuresEnd++;
                    }
                    if (arrivals[tpg_graph[i].fromId] <= offset) {
                        arrivals[tpg_graph[i].fromId] = offset + 1;
                        arrivalsEnd++;
                    }
                }

                else if (forwardFlag[tpg_graph[i].fromId] == offset + 2 && forsize == 1 && forwardEV[forwardLocator[tpg_graph[i].fromId]].time < tpg_graph[i].time) {
                    results[resultEnd++] = tpg_graph[i].edgeId;
                    hasVerified[tpg_graph[i].edgeId] = offset + 1;
                    if (departures[tpg_graph[i].toId] <= offset) {
                        departures[tpg_graph[i].toId] = offset + 1;
                        departuresEnd++;
                    }
                }
                
                else if (backwardFlag[tpg_graph[i].toId] == offset + 2 && backsize == 1 && backwardEV[backwardLocator[tpg_graph[i].toId]].time > tpg_graph[i].time) {
                    results[resultEnd++] = tpg_graph[i].edgeId;
                    hasVerified[tpg_graph[i].edgeId] = offset + 1;
                    if (arrivals[tpg_graph[i].fromId] <= offset) {
                        arrivals[tpg_graph[i].fromId] = offset + 1;
                        arrivalsEnd++;
                    }
                }

                else {
                    edgesForVerification[edgesForVerificationEnd++] = tpg_graph[i].edgeId;
                }
            }

            delete[] tmp;
        } 
    }

    new_sortInEdges = new PerEdge[upperEnd];
    std::copy(upper_bound_graph, upper_bound_graph + upperEnd, new_sortInEdges);
    std::sort(new_sortInEdges, new_sortInEdges + upperEnd, sortByToIdTime);
    
    EdgeID curLocator = 0;
    VertexID curId = 0, id;
    int i = 0;

    while (i < upperEnd) {
        id = new_sortInEdges[i].toId;
        while (curId <= id) {
            new_forwardLocator[curId] = curLocator;
            curId++;
        }
        new_forwardLocatorEnd++;
        while (i < upperEnd && new_sortInEdges[i].toId == id) {
            curLocator++;
            i++;
        }
    }
    while (curId <= VN) {
        new_forwardLocator[curId] = curLocator;
        curId++;
    }

    new_sortOutEdges = new PerEdge[upperEnd];
    std::copy(upper_bound_graph, upper_bound_graph + upperEnd, new_sortOutEdges);
    std::sort(new_sortOutEdges, new_sortOutEdges + upperEnd, sortFromIdTime);

    curLocator = 0;
    curId = 0;
    i = 0;

    while (i < upperEnd) {
        id = new_sortOutEdges[i].fromId;
        while (curId <= id) {
            new_backwardLocator[curId] = curLocator;
            curId++;
        }
        new_backwardLocatorEnd++;
        while (i < upperEnd && new_sortOutEdges[i].fromId == id) {
            curLocator++;
            i++;
        }
    }
    while (curId <= VN) {
        new_backwardLocator[curId] = curLocator;
        curId++;
    }

    // statistics
    #ifdef WRITE_STATISTICS
        numOfUpperbound[queryId] = upperEnd;
    #endif

    // if (max < upperEnd)
    //     max = upperEnd;

}

void TSPG::verifyUndeterminedEdge() {                              
    VertexID from, to;
    Time time;
    
    for (int i = 0; i < edgesForVerificationEnd; ++i) {
        
        if (hasVerified[edgesForVerification[i]] <= offset) {
            
            stackEnd = 0, stackeEnd = 0;
            maxStack = 0;
            from = edges[edgesForVerification[i]].fromId;
            to = edges[edgesForVerification[i]].toId;
            time = edges[edgesForVerification[i]].time;

            if (t_end - time < time - t_begin) {
                if (forwardSearch(from, to, edgesForVerification[i], time, time)) {
                    addToResult(); 
                }

                else {
                    for (EdgeID edgeLocators = new_backwardLocator[from]; edgeLocators < new_backwardLocator[from + 1]; ++edgeLocators) {
                        if (new_sortOutEdges[edgeLocators].time < time)
                            break;
                        else {
                            if (new_sortOutEdges[edgeLocators].time == time) {
                                if (new_sortOutEdges[edgeLocators].toId == to) {
                                    hasVerified[new_sortOutEdges[edgeLocators].edgeId] = offset + 1;
                                }
                            }
                        }
                    }
                } 
            }
            
            else {
                if (backward2Search(from, to, edgesForVerification[i], time, time)) {
                    addToResult();
                }
                else {
                    for (EdgeID edgeLocators = new_backwardLocator[from]; edgeLocators < new_backwardLocator[from + 1]; ++edgeLocators) {
                        if (new_sortOutEdges[edgeLocators].time < time)
                            break;
                        else {
                            if (new_sortOutEdges[edgeLocators].time == time) {
                                if (new_sortOutEdges[edgeLocators].toId == to) {
                                    hasVerified[new_sortOutEdges[edgeLocators].edgeId] = offset + 1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

bool TSPG::forwardSearch(VertexID& from, VertexID& to, EdgeID& e, Time& time, Time& fromtime) {
    bool ans = false;
    stack[stackEnd++] = to;
    isinStack[to] = offset + 1;
    stacke[stackeEnd++] = e;

    if (arrivals[to] == offset + 1) {
        arrival = to;
        arrivaltime = time;
        
        if (backwardSearch(from, e, fromtime)) {
            return true;
        }

        else {
            for (EdgeID edgeLocators = new_backwardLocator[to]; edgeLocators < new_backwardLocator[to + 1]; ++edgeLocators) {
                VertexID& u = new_sortOutEdges[edgeLocators].toId;
                if (isinStack[u] <= offset && u != from && new_sortOutEdges[edgeLocators].time > time) {
                    ans = forwardSearch(from, u, new_sortOutEdges[edgeLocators].edgeId, new_sortOutEdges[edgeLocators].time, fromtime);
                    if (ans) {
                        return true;
                    }                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
                }
                else if (new_sortOutEdges[edgeLocators].time <= time)
                    break;
            }

            // EdgeID *search = new EdgeID[new_backwardLocator[to + 1]-new_backwardLocator[to]+1]();
            // for (EdgeID edgeLocators = new_backwardLocator[to]; edgeLocators < new_backwardLocator[to + 1]; ++edgeLocators) {
            //     VertexID& u = new_sortOutEdges[edgeLocators].toId;
            //     if (isinStack[u] <= offset && u != from && new_sortOutEdges[edgeLocators].time > time && !hasVerified[new_sortOutEdges[edgeLocators].edgeId]) {
            //         search[edgeLocators-new_backwardLocator[to]] = 1;
            //         ans = forwardSearch(from, u, new_sortOutEdges[edgeLocators].edgeId, new_sortOutEdges[edgeLocators].time, fromtime);
            //         if (ans) {
            //             delete[] search;
            //             return true;
            //         }
            //     }
            //     else if (new_sortOutEdges[edgeLocators].time <= time)
            //         break;
            // }
            // for (EdgeID edgeLocators = new_backwardLocator[to]; edgeLocators < new_backwardLocator[to + 1]; ++edgeLocators) {
            //     VertexID& u = new_sortOutEdges[edgeLocators].toId;
            //     if (isinStack[u] <= offset && u != from && new_sortOutEdges[edgeLocators].time > time && !search[edgeLocators-new_backwardLocator[to]]) {
            //         search[edgeLocators-new_backwardLocator[to]] = 1;
            //         ans = forwardSearch(from, u, new_sortOutEdges[edgeLocators].edgeId, new_sortOutEdges[edgeLocators].time, fromtime);
            //         if (ans) {
            //             delete[] search;
            //             return true;
            //         }
            //     }
            //     else if (new_sortOutEdges[edgeLocators].time <= time)
            //         break;
            // }
            // delete[] search;
        }
    }

    else {

        for (EdgeID edgeLocators = new_backwardLocator[to]; edgeLocators < new_backwardLocator[to + 1]; ++edgeLocators) {
            VertexID& u = new_sortOutEdges[edgeLocators].toId;
            if (isinStack[u] <= offset && u != from && new_sortOutEdges[edgeLocators].time > time) {
                ans = forwardSearch(from, u, new_sortOutEdges[edgeLocators].edgeId, new_sortOutEdges[edgeLocators].time, fromtime);
                if (ans) {
                    return true;
                }
            }
            else if (new_sortOutEdges[edgeLocators].time <= time)
                break;
        }

        // EdgeID *search = new EdgeID[new_backwardLocator[to + 1]-new_backwardLocator[to]+1]();
        // for (EdgeID edgeLocators = new_backwardLocator[to]; edgeLocators < new_backwardLocator[to + 1]; ++edgeLocators) {
        //     VertexID& u = new_sortOutEdges[edgeLocators].toId;
        //     if (isinStack[u] <= offset && u != from && new_sortOutEdges[edgeLocators].time > time && !hasVerified[new_sortOutEdges[edgeLocators].edgeId]) {
        //         search[edgeLocators-new_backwardLocator[to]] = 1;
        //         ans = forwardSearch(from, u, new_sortOutEdges[edgeLocators].edgeId, new_sortOutEdges[edgeLocators].time, fromtime);
        //         if (ans) {
        //             delete[] search;
        //             return true;
        //         }
        //     }
        //     else if (new_sortOutEdges[edgeLocators].time <= time)
        //         break;
        // }
        // for (EdgeID edgeLocators = new_backwardLocator[to]; edgeLocators < new_backwardLocator[to + 1]; ++edgeLocators) {
        //     VertexID& u = new_sortOutEdges[edgeLocators].toId;
        //     if (isinStack[u] <= offset && u != from && new_sortOutEdges[edgeLocators].time > time && !search[edgeLocators-new_backwardLocator[to]]) {
        //         search[edgeLocators-new_backwardLocator[to]] = 1;
        //         ans = forwardSearch(from, u, new_sortOutEdges[edgeLocators].edgeId, new_sortOutEdges[edgeLocators].time, fromtime);
        //         if (ans) {
        //             delete[] search;
        //             return true;
        //         }
        //     }
        //     else if (new_sortOutEdges[edgeLocators].time <= time)
        //         break;
        // }
        // delete[] search;
    }

    stackEnd--;
    stackeEnd--;
    isinStack[to] = offset;

    return false;
}

bool TSPG::backwardSearch(VertexID& from, EdgeID& e, Time& time) {
    bool ans = false;
    stack[stackEnd++] = from;
    isinStack[from] = offset + 1;
    stacke[stackeEnd++] = e;

    maxStack = max(stackEnd, maxStack);

    if (departures[from] == offset + 1) {
        departure = from;
        departuretime = time;
        if (tryAddEdges(departure, arrival, departuretime, arrivaltime)) {
            return true;
        }
        else {
            
            for (EdgeID edgeLocators = new_forwardLocator[from]; edgeLocators < new_forwardLocator[from + 1]; ++edgeLocators) {
                VertexID& u = new_sortInEdges[edgeLocators].fromId;
                if (isinStack[u] <= offset && new_sortInEdges[edgeLocators].time < time) {
                    ans = backwardSearch(u, new_sortInEdges[edgeLocators].edgeId, new_sortInEdges[edgeLocators].time);
                    if (ans) {
                        return true;
                    }
                }
                else if (new_sortInEdges[edgeLocators].time >= time)
                    break;
            }
            
            // EdgeID *search = new EdgeID[new_forwardLocator[from + 1]-new_forwardLocator[from]+1]();
            // for (EdgeID edgeLocators = new_forwardLocator[from]; edgeLocators < new_forwardLocator[from + 1]; ++edgeLocators) {
            //     VertexID& u = new_sortInEdges[edgeLocators].fromId;
            //     if (isinStack[u] <= offset && new_sortInEdges[edgeLocators].time < time && !hasVerified[new_sortInEdges[edgeLocators].edgeId]) {
            //         search[edgeLocators-new_forwardLocator[from]] = 1;
            //         ans = backwardSearch(u, new_sortInEdges[edgeLocators].edgeId, new_sortInEdges[edgeLocators].time);
            //         if (ans) {
            //             delete[] search;
            //             return true;
            //         }
            //     }
            //     else if (new_sortInEdges[edgeLocators].time >= time)
            //         break;
            // }
            // for (EdgeID edgeLocators = new_forwardLocator[from]; edgeLocators < new_forwardLocator[from + 1]; ++edgeLocators) {
            //     VertexID& u = new_sortInEdges[edgeLocators].fromId;
            //     if (isinStack[u] <= offset && new_sortInEdges[edgeLocators].time < time && !search[edgeLocators-new_forwardLocator[from]]) {
            //         search[edgeLocators-new_forwardLocator[from]] = 1;
            //         ans = backwardSearch(u, new_sortInEdges[edgeLocators].edgeId, new_sortInEdges[edgeLocators].time);
            //         if (ans) {
            //             delete[] search;
            //             return true;
            //         }
            //     }
            //     else if (new_sortInEdges[edgeLocators].time >= time)
            //         break;
            // }
            // delete[] search;
        }
    }

    else {
        
        for (EdgeID edgeLocators = new_forwardLocator[from]; edgeLocators < new_forwardLocator[from + 1]; ++edgeLocators) {
            VertexID& u = new_sortInEdges[edgeLocators].fromId;
            if (isinStack[u] <= offset && new_sortInEdges[edgeLocators].time < time) {
                ans = backwardSearch(u, new_sortInEdges[edgeLocators].edgeId, new_sortInEdges[edgeLocators].time);
                if (ans) {
                    return true;
                }
            }
            else if (new_sortInEdges[edgeLocators].time >= time)
                break;
        }

        // EdgeID *search = new EdgeID[new_forwardLocator[from + 1]-new_forwardLocator[from]+1]();
        // for (EdgeID edgeLocators = new_forwardLocator[from]; edgeLocators < new_forwardLocator[from + 1]; ++edgeLocators) {
        //     VertexID& u = new_sortInEdges[edgeLocators].fromId;
        //     if (isinStack[u] <= offset && new_sortInEdges[edgeLocators].time < time && !hasVerified[new_sortInEdges[edgeLocators].edgeId]) {
        //         search[edgeLocators-new_forwardLocator[from]] = 1;
        //         ans = backwardSearch(u, new_sortInEdges[edgeLocators].edgeId, new_sortInEdges[edgeLocators].time);
        //         if (ans) {
        //             delete[] search;
        //             return true;
        //         }
        //     }
        //     else if (new_sortInEdges[edgeLocators].time >= time)
        //         break;
        // }
        // for (EdgeID edgeLocators = new_forwardLocator[from]; edgeLocators < new_forwardLocator[from + 1]; ++edgeLocators) {
        //     VertexID& u = new_sortInEdges[edgeLocators].fromId;
        //     if (isinStack[u] <= offset && new_sortInEdges[edgeLocators].time < time && !search[edgeLocators-new_forwardLocator[from]]) {
        //         search[edgeLocators-new_forwardLocator[from]] = 1;
        //         ans = backwardSearch(u, new_sortInEdges[edgeLocators].edgeId, new_sortInEdges[edgeLocators].time);
        //         if (ans) {
        //             delete[] search;
        //             return true;
        //         }
        //     }
        //     else if (new_sortInEdges[edgeLocators].time >= time)
        //         break;
        // }
        // delete[] search;
    }

    stackEnd--;
    stackeEnd--;
    isinStack[from] = offset;

    return false;
}

bool TSPG::backward2Search(VertexID& from, VertexID& to, EdgeID& e, Time& time, Time& totime) {
    bool ans = false;
    stack[stackEnd++] = from;
    isinStack[from] = offset + 1;
    stacke[stackeEnd++] = e;

    if (departures[from] == offset + 1) {
        departure = from;
        departuretime = time;
        if (forward2Search(to, e, totime)) {
            return true;
        }
        else {
            
            for (EdgeID edgeLocators = new_forwardLocator[from]; edgeLocators < new_forwardLocator[from + 1]; ++edgeLocators) {
                VertexID& u = new_sortInEdges[edgeLocators].fromId;
                if (isinStack[u] <= offset && u != to && new_sortInEdges[edgeLocators].time < time) {
                    ans = backward2Search(u, to, new_sortInEdges[edgeLocators].edgeId, new_sortInEdges[edgeLocators].time, totime);
                    if (ans) {
                        return true;
                    }
                }
                else if (new_sortInEdges[edgeLocators].time >= time)
                    break;
            }

            // EdgeID *search = new EdgeID[new_forwardLocator[from + 1]-new_forwardLocator[from]+1]();
            // for (EdgeID edgeLocators = new_forwardLocator[from]; edgeLocators < new_forwardLocator[from + 1]; ++edgeLocators) {
            //     VertexID& u = new_sortInEdges[edgeLocators].fromId;
            //     if (isinStack[u] <= offset && u != to && new_sortInEdges[edgeLocators].time < time && !hasVerified[new_sortInEdges[edgeLocators].edgeId]) {
            //         search[edgeLocators-new_forwardLocator[from]] = 1;
            //         ans = backward2Search(u, to, new_sortInEdges[edgeLocators].edgeId, new_sortInEdges[edgeLocators].time, totime);
            //         if (ans) {
            //             delete[] search;
            //             return true;
            //         }
            //     }
            //     else if (new_sortInEdges[edgeLocators].time >= time)
            //         break;
            // }
            // for (EdgeID edgeLocators = new_forwardLocator[from]; edgeLocators < new_forwardLocator[from + 1]; ++edgeLocators) {
            //     VertexID& u = new_sortInEdges[edgeLocators].fromId;
            //     if (isinStack[u] <= offset && u != to && new_sortInEdges[edgeLocators].time < time && !search[edgeLocators-new_forwardLocator[from]]) {
            //         search[edgeLocators-new_forwardLocator[from]] = 1;
            //         ans = backward2Search(u, to, new_sortInEdges[edgeLocators].edgeId, new_sortInEdges[edgeLocators].time, totime);
            //         if (ans) {
            //             delete[] search;
            //             return true;
            //         }
            //     }
            //     else if (new_sortInEdges[edgeLocators].time >= time)
            //         break;
            // }
            // delete[] search;
        }
    }

    else {
        
        for (EdgeID edgeLocators = new_forwardLocator[from]; edgeLocators < new_forwardLocator[from + 1]; ++edgeLocators) {
            VertexID& u = new_sortInEdges[edgeLocators].fromId;
            if (isinStack[u] <= offset && u != to && new_sortInEdges[edgeLocators].time < time) {
                ans = backward2Search(u, to, new_sortInEdges[edgeLocators].edgeId, new_sortInEdges[edgeLocators].time, totime);
                if (ans) {
                    return true;
                }
            }
            else if (new_sortInEdges[edgeLocators].time >= time)
                break;
        }

        // EdgeID *search = new EdgeID[new_forwardLocator[from + 1]-new_forwardLocator[from]+1]();
        // for (EdgeID edgeLocators = new_forwardLocator[from]; edgeLocators < new_forwardLocator[from + 1]; ++edgeLocators) {
        //     VertexID& u = new_sortInEdges[edgeLocators].fromId;
        //     if (isinStack[u] <= offset && u != to && new_sortInEdges[edgeLocators].time < time && !hasVerified[new_sortInEdges[edgeLocators].edgeId]) {
        //         search[edgeLocators-new_forwardLocator[from]] = 1;
        //         ans = backward2Search(u, to, new_sortInEdges[edgeLocators].edgeId, new_sortInEdges[edgeLocators].time, totime);
        //         if (ans) {
        //             delete[] search;
        //             return true;
        //         }
        //     }
        //     else if (new_sortInEdges[edgeLocators].time >= time)
        //         break;
        // }
        // for (EdgeID edgeLocators = new_forwardLocator[from]; edgeLocators < new_forwardLocator[from + 1]; ++edgeLocators) {
        //     VertexID& u = new_sortInEdges[edgeLocators].fromId;
        //     if (isinStack[u] <= offset && u != to && new_sortInEdges[edgeLocators].time < time && !search[edgeLocators-new_forwardLocator[from]]) {
        //         search[edgeLocators-new_forwardLocator[from]] = 1;
        //         ans = backward2Search(u, to, new_sortInEdges[edgeLocators].edgeId, new_sortInEdges[edgeLocators].time, totime);
        //         if (ans) {
        //             delete[] search;
        //             return true;
        //         }
        //     }
        //     else if (new_sortInEdges[edgeLocators].time >= time)
        //         break;
        // }
        // delete[] search;
    }

    stackEnd--;
    stackeEnd--;
    isinStack[from] = offset;

    return false;
}

bool TSPG::forward2Search(VertexID& to, EdgeID& e, Time& time) {
    bool ans = false;
    stack[stackEnd++] = to;
    isinStack[to] = offset + 1;
    stacke[stackeEnd++] = e;

    maxStack = max(stackEnd, maxStack);

    if (arrivals[to] == offset + 1) {
        arrival = to;
        arrivaltime = time;
        if (tryAddEdges(departure, arrival, departuretime, arrivaltime)) {
            return true;
        }
        else {

            for (EdgeID edgeLocators = new_backwardLocator[to]; edgeLocators < new_backwardLocator[to + 1]; ++edgeLocators) {
                VertexID& u = new_sortOutEdges[edgeLocators].toId;
                if (isinStack[u] <= offset && new_sortOutEdges[edgeLocators].time > time) {
                    ans = forward2Search(u, new_sortOutEdges[edgeLocators].edgeId, new_sortOutEdges[edgeLocators].time);
                    if (ans) {
                        return true;
                    }
                }
                else if (new_sortOutEdges[edgeLocators].time <= time)
                    break;
            }

            // EdgeID *search = new EdgeID[new_backwardLocator[to + 1]-new_backwardLocator[to]+1]();
            // for (EdgeID edgeLocators = new_backwardLocator[to]; edgeLocators < new_backwardLocator[to + 1]; ++edgeLocators) {
            //     VertexID& u = new_sortOutEdges[edgeLocators].toId;
            //     if (isinStack[u] <= offset && new_sortOutEdges[edgeLocators].time > time && !hasVerified[new_sortOutEdges[edgeLocators].edgeId]) {
            //         search[edgeLocators-new_backwardLocator[to]] = 1;
            //         ans = forward2Search(u, new_sortOutEdges[edgeLocators].edgeId, new_sortOutEdges[edgeLocators].time);
            //         if (ans) {
            //             delete[] search;
            //             return true;
            //         }
            //     }
            //     else if (new_sortOutEdges[edgeLocators].time <= time)
            //         break;
            // }
            // for (EdgeID edgeLocators = new_backwardLocator[to]; edgeLocators < new_backwardLocator[to + 1]; ++edgeLocators) {
            //     VertexID& u = new_sortOutEdges[edgeLocators].toId;
            //     if (isinStack[u] <= offset && new_sortOutEdges[edgeLocators].time > time && !search[edgeLocators-new_backwardLocator[to]]) {
            //         search[edgeLocators-new_backwardLocator[to]] = 1;
            //         ans = forward2Search(u, new_sortOutEdges[edgeLocators].edgeId, new_sortOutEdges[edgeLocators].time);
            //         if (ans) {
            //             delete[] search;
            //             return true;
            //         }
            //     }
            //     else if (new_sortOutEdges[edgeLocators].time <= time)
            //         break;
            // }
            // delete[] search;
        }
    }

    else {
        
        for (EdgeID edgeLocators = new_backwardLocator[to]; edgeLocators < new_backwardLocator[to + 1]; ++edgeLocators) {
            VertexID& u = new_sortOutEdges[edgeLocators].toId;
            if (isinStack[u] <= offset && new_sortOutEdges[edgeLocators].time > time) {
                ans = forward2Search(u, new_sortOutEdges[edgeLocators].edgeId, new_sortOutEdges[edgeLocators].time);
                if (ans) {
                    return true;
                }
            }
            else if (new_sortOutEdges[edgeLocators].time <= time)
                break;
        }

        // EdgeID *search = new EdgeID[new_backwardLocator[to + 1]-new_backwardLocator[to]+1]();
        // for (EdgeID edgeLocators = new_backwardLocator[to]; edgeLocators < new_backwardLocator[to + 1]; ++edgeLocators) {
        //     VertexID& u = new_sortOutEdges[edgeLocators].toId;
        //     if (isinStack[u] <= offset && new_sortOutEdges[edgeLocators].time > time && !hasVerified[new_sortOutEdges[edgeLocators].edgeId]) {
        //         search[edgeLocators-new_backwardLocator[to]] = 1;
        //         ans = forward2Search(u, new_sortOutEdges[edgeLocators].edgeId, new_sortOutEdges[edgeLocators].time);
        //         if (ans) {
        //             delete[] search;
        //             return true;
        //         }
        //     }
        //     else if (new_sortOutEdges[edgeLocators].time <= time)
        //         break;
        // }
        // for (EdgeID edgeLocators = new_backwardLocator[to]; edgeLocators < new_backwardLocator[to + 1]; ++edgeLocators) {
        //     VertexID& u = new_sortOutEdges[edgeLocators].toId;

        //     if (isinStack[u] <= offset && new_sortOutEdges[edgeLocators].time > time && !search[edgeLocators-new_backwardLocator[to]]) {
        //         search[edgeLocators-new_backwardLocator[to]] = 1;
        //         ans = forward2Search(u, new_sortOutEdges[edgeLocators].edgeId, new_sortOutEdges[edgeLocators].time);
        //         if (ans) {
        //             delete[] search;
        //             return true;
        //         }
        //     }
        //     else if (new_sortOutEdges[edgeLocators].time <= time)
        //         break;
        // }
        // delete[] search;
    }

    stackEnd--;
    stackeEnd--;
    isinStack[to] = offset;

    return false;
}

bool TSPG::tryAddEdges(VertexID& d, VertexID& a, Time& dt, Time& at) {
    
    for (EdgeID edgeLocator = new_forwardLocator[d]; edgeLocator < new_forwardLocator[d + 1]; edgeLocator++) {
        VertexID& InD = new_sortInEdges[edgeLocator].fromId;
        if (isinStack[InD] <= offset && new_sortInEdges[edgeLocator].time < dt && forwardFlag[InD] == offset + 2) {
            stack[stackEnd++] = InD;
            isinStack[InD] = offset + 1;
            for (EdgeID edgeLocators = new_backwardLocator[a]; edgeLocators < new_backwardLocator[a + 1]; edgeLocators++) {
                VertexID& OutA = new_sortOutEdges[edgeLocators].toId;
                if (isinStack[OutA] <= offset && new_sortOutEdges[edgeLocators].time > at && backwardFlag[OutA] == offset + 2) {
                    stackEnd--;
                    isinStack[InD] = offset;
                    return true;
                }
            }
            stackEnd--;
            isinStack[InD] = offset;
        }
    }

    return false;
}

void TSPG::addToResult() {
    
    vector<PerEdge> curEdges;

    for (int j = 0; j < stackeEnd; ++j) {
        isinStack[stack[j]] = offset;
        curEdges.push_back(edges[stacke[j]]);
    }

    std::sort(curEdges.begin(), curEdges.end(), sortByEdgeTime);

    Time pTime, bTime;

    for (int j = 0; j < curEdges.size(); ++j) {
        if (hasVerified[curEdges[j].edgeId] <= offset) {
            if (j == 0) {
                pTime = curEdges[j].time;
                if (curEdges[j].time == curEdges[j + 1].time) {
                    if (j + 1 == curEdges.size() - 1)
                        bTime = curEdges[j].time;
                    else
                        bTime = curEdges[j + 2].time - 1;
                }
                else
                    bTime = curEdges[j + 1].time - 1;
            } 
            
            else if (j == curEdges.size() - 1){
                pTime = curEdges[j - 1].time + 1;
                bTime = curEdges[j].time;
            }

            else {
                pTime = curEdges[j - 1].time + 1;
                if (curEdges[j].time == curEdges[j + 1].time) {
                    if (j + 1 == curEdges.size() - 1)
                        bTime = curEdges[j].time;
                    else
                        bTime = curEdges[j + 2].time - 1;
                }
                else
                    bTime = curEdges[j + 1].time - 1;
            }

            for (EdgeID edgeLocators = new_backwardLocator[curEdges[j].fromId]; edgeLocators < new_backwardLocator[curEdges[j].fromId + 1]; edgeLocators++) {
                if (new_sortOutEdges[edgeLocators].time < pTime)
                    break;
                else {
                    if (new_sortOutEdges[edgeLocators].time <= bTime && new_sortOutEdges[edgeLocators].toId == curEdges[j].toId) {
                        if (hasVerified[curEdges[j].edgeId] <= offset) {
                            results[resultEnd++] = curEdges[j].edgeId;
                            hasVerified[curEdges[j].edgeId] = offset + 1;
                        }
                    }
                }
            }
        }
    }

    // for (int j = 0; j < stackeEnd; ++j) {
    //     isinStack[stack[j]] = offset; 
    //     if (hasVerified[stacke[j]] <= offset) {
    //         results[resultEnd++] = stacke[j];
    //         hasVerified[stacke[j]] = offset + 1;
    //     }
    // }
}

void TSPG::DFS(VertexID &u, Time &time) {
    stack[stackEnd++] = u;
    isinStack[u] = offset + 1;

    if (u == s) time = time - 1;

    if (u == t) {
        addEdge();
    }
    else {
        for (EdgeID edgeLocators = new_backwardLocator[u]; edgeLocators < new_backwardLocator[u + 1]; ++edgeLocators) {
            VertexID& v = new_sortOutEdges[edgeLocators].toId;
            if (isinStack[v] <= offset && new_sortOutEdges[edgeLocators].time > time) {
                stacke[stackeEnd++] = new_sortOutEdges[edgeLocators].edgeId;
                DFS(v, new_sortOutEdges[edgeLocators].time);
            }           
        }
    }
    stackEnd--;
    stackeEnd--;
    isinStack[u] = 0;
}

void TSPG::addEdge() {
    maxStack = max(stackEnd, maxStack);
    for (int k = 0; k < stackeEnd; ++k) {
        if (!hasVerified[stacke[k]]) {
            results[resultEnd++] = stacke[k];
            hasVerified[stacke[k]] = 1;
        }
    }
}

// refresh memory for new query
inline void TSPG::refreshMemory() {
    // refresh offset
    offset += t_end + 1;
    
    if (offset >= 4294964295) {
        offset = 0;
        
        memset(earliestTime, 0, sizeof(Time)*VN);
        memset(latestTime, 0, sizeof(Time)*VN);
        memset(dropflag, 0, sizeof(VertexID)*VN);

        memset(forwardFlag, 0, sizeof(int)*VN);
        memset(forwardEVLen, 0, sizeof(int)*VN);
        memset(backwardFlag, 0, sizeof(int)*VN);
        memset(backwardEVLen, 0, sizeof(int)*VN);

        memset(departures, 0, sizeof(VertexID)*VN);
        memset(arrivals, 0, sizeof(VertexID)*VN);

        memset(hasVerified, 0, sizeof(EdgeID)*EN);
        memset(isinStack, 0, sizeof(VertexID)*VN);
    }

    // refresh storages
    for (int i = 0; i < tpgEnd; ++i) {
        delete[] forwardEV[i].vertices;
        forwardEV[i].vertices = nullptr;
    }
    delete[] forwardEV;

    for (int i = 0; i < tpgEnd; ++i) {
        delete[] backwardEV[i].vertices;
        backwardEV[i].vertices = nullptr;
    }
    delete[] backwardEV;

    delete[] new_sortInEdges;
    delete[] new_sortOutEdges;
    
    // statistics
    #ifdef WRITE_STATISTICS
        maxFrontierSize = 1;
    #endif
}

// free up memories
void STPG::cleanUp() {
    
    delete[] nextFrontier;
    delete[] forwardFrontier;
    delete[] backwardFrontier;

    delete[] earliestTime;
    delete[] latestTime;
    delete[] isvisitedEdge;

    delete[] dropflag;
    delete[] tpg_graph;
    delete[] upper_bound_graph;

    delete[] forwardLocator;
    delete[] forwardFlag;
    delete[] forwardEVLen;
    delete[] backwardLocator;
    delete[] backwardFlag;
    delete[] backwardEVLen;
    delete[] new_forwardLocator;
    delete[] new_backwardLocator;

    delete[] edgesForVerification;
    delete[] departures;
    delete[] arrivals;

    delete[] hasVerified;
    delete[] stack;
    delete[] stacke;
    delete[] isinStack;

    delete[] results;
}

// for statistics
#ifdef WRITE_STATISTICS

    // return space cost of current query
    double TSPG::getCurrentSpaceCost() {
        double spaceCost = 0;

        // earliest and latest time
        spaceCost += sizeof(VertexID)*maxFrontierSize;          // nextFrontier
        spaceCost += sizeof(VertexID)*maxFrontierSize;          // forwardFrontier
        spaceCost += sizeof(VertexID)*maxFrontierSize;          // backwardFrontier
        
        spaceCost += sizeof(Time)*VN;                           // earliestTime
        spaceCost += sizeof(Time)*VN;                           // latestTime
        spaceCost += sizeof(EdgeID)*EN;                         // isvisitedEdge

        // TPG
        spaceCost += sizeof(VertexID)*VN;                       // dropflag
        spaceCost += sizeof(PerEdge)*tpgEnd;                    // tpg_graph
        spaceCost += sizeof(PerEdge)*tpgEnd;                    // sortInEdges
        spaceCost += sizeof(PerEdge)*tpgEnd;                    // sortOutEdges

        // essential vertices
        spaceCost += sizeof(VertexID)*forwardEVEnd;             // forwardEV
        spaceCost += sizeof(VertexID)*backwardEVEnd;            // backwardEV
        spaceCost += sizeof(int)*forwardLocatorEnd;             // forwardFlag
        spaceCost += sizeof(EdgeID)*forwardLocatorEnd;          // forwardLocator
        spaceCost += sizeof(int)*backwardLocatorEnd;            // backwardFlag
        spaceCost += sizeof(EdgeID)*backwardLocatorEnd;         // backwardLocator
        spaceCost += sizeof(int)*forwardLocatorEnd;             // forwardEVLen
        spaceCost += sizeof(int)*backwardLocatorEnd;            // backwardEVLen
        
        // UpperBound Graph
        spaceCost += sizeof(VertexID)*departuresEnd;            // departures
        spaceCost += sizeof(VertexID)*arrivalsEnd;              // arrivals
        spaceCost += sizeof(EdgeID)*edgesForVerificationEnd;    // Undetermined edges
        spaceCost += sizeof(PerEdge)*upperEnd;                  // upper_bound_graph
        spaceCost += sizeof(PerEdge)*upperEnd;                  // new_sortInEdges
        spaceCost += sizeof(PerEdge)*upperEnd;                  // new_sortOutEdges
        spaceCost += sizeof(EdgeID)*new_forwardLocatorEnd;      // new_forwardLocator
        spaceCost += sizeof(EdgeID)*new_backwardLocatorEnd;     // new_backwardLocator

        // verifyUndeterminedEdge
        spaceCost += sizeof(EdgeID)*edgesForVerificationEnd;    // hasVerified
        spaceCost += sizeof(VertexID)*new_forwardLocatorEnd;    // isinStack
        spaceCost += sizeof(VertexID)*maxStack;                 // stack
        spaceCost += sizeof(EdgeID)*maxStack;                   // stacke

        // results
        spaceCost += sizeof(EdgeID)*resultEnd;                  // results

        return spaceCost;
    }
    
    void TSPG::initStatisticStorage(int queryNumber){
        numOfAnswers = new EdgeID[queryNumber]();
        numOfUpperbound = new EdgeID[queryNumber]();
        numOfTPAnswers = new VertexID[queryNumber]();
        spaceCosts = new double[queryNumber]();
        // open statistics file
        statisticsFile.open("../" +statisticsPath+extractFilename(queryFilename)+".csv",ios::app);
        statisticsFile << "Space cost (bytes),# TPG edges,# Upper-bound edges,# TSPG edges" << endl;
    }

    void TSPG::cleanUpStatisticStorage(){
        delete[] numOfAnswers;
        delete[] numOfUpperbound;
        delete[] numOfTPAnswers;
        delete[] spaceCosts;

        statisticsFile.close();
    }

#endif

#endif