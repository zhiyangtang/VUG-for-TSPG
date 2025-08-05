#ifndef STPG_CC
#define STPG_CC
#include "STPG.h"
using namespace std;

STPG::STPG(Graph* inputGraph) {

    graph = inputGraph;
    VN = graph->VN;
    EN = graph->EN;
    edges = graph->edges;
    inNeighborsLocator = graph->inNeighborsLocator;
    inNeighbors_time = graph->inNeighbors_time;
    outNeighborsLocator = graph->outNeighborsLocator;
    outNeighbors_time = graph->outNeighbors_time;

    init();
}


void STPG::init() {
    
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


void STPG::answerAllQueries(vector<PerQuery>& queries) {
    if (queries.size()==0) {
        printf("!No query\n");
        return;
    }
    
    #ifdef WRITE_ANSWERS
        resultFile.open("../"+answerPath+queryFilename+".answer",ios::app);
        resultFile << "number of edges,edge ids" << endl;
    #endif
    #ifdef WRITE_STATISTICS
        queryId = 0;
        initStatisticStorage(queries.size());
    #endif

    printf("Running STPG ...\n");
    double startTime = getCurrentTimeInMs();

    for (PerQuery& query : queries) {

        resultEnd = executeQuery(query.source, query.target, query.star, query.end);

        #ifdef WRITE_ANSWERS
            std::sort(results, results + resultEnd);
            resultFile << resultEnd;
            for (EdgeID j = 0; j < resultEnd; j++) {
                resultFile << "," << results[j];
            }
            resultFile << endl;
        #endif
        #ifdef WRITE_STATISTICS
            statisticsFile << spaceCosts[queryId] << "," << numOfTPAnswers[queryId] << "," << numOfUpperbound[queryId] << "," << numOfAnswers[queryId] << endl;
            queryId++;
        #endif    
    }

    double timeCost = getCurrentTimeInMs() - startTime;
    printf("- Finish. Time cost: %.2f ms\n", timeCost);
    logFile << str(timeCost) << endl;

    #ifdef WRITE_STATISTICS
        cleanUpStatisticStorage();
    #endif
    #ifdef WRITE_ANSWERS
        resultFile.close();
    #endif
}


EdgeID STPG::executeQuery(VertexID source, VertexID target, Time begin, Time end) {
    
    s = source;
    t = target;
    t_begin = begin;
    t_end = end;
    resultEnd = 0;

    earliestTimeBFS();
    latestTimeBFS();
    prunedVertices();
    prunedEdges();

    forward();
    backward();
    UpperBound();

    if (edgesForVerificationEnd) {
        verifyUndeterminedEdge();
    }

    refreshMemory();

    #ifdef WRITE_STATISTICS
        numOfAnswers[queryId] = resultEnd;
        spaceCosts[queryId] = getCurrentSpaceCost();
    #endif

    return resultEnd;
}


void STPG::earliestTimeBFS() {
    earliestTime[s] = offset;
    forwardFrontier[0] = s;
    forwardFrontierEnd = 1;
    memset(isvisitedEdge, 0, sizeof(EdgeID)*EN);

    while (forwardFrontierEnd > 0) {
        nextFrontierEnd = 0;

        for (VertexID i = 0; i < forwardFrontierEnd; i++) {
            VertexID& u = forwardFrontier[i];

            for (EdgeID edgeLocator = outNeighborsLocator[u]; edgeLocator < outNeighborsLocator[u + 1]; edgeLocator++) {
                Time& time = outNeighbors_time[edgeLocator].time;
                VertexID& v = outNeighbors_time[edgeLocator].neighbor;
                EdgeID& edgeID = outNeighbors_time[edgeLocator].edgeId;
                if (isvisitedEdge[edgeID] == 1 || time < t_begin || v == s) continue;   
                else if (time > t_end) break;
                else {
                    if (time + offset > earliestTime[u]) {
                        isvisitedEdge[edgeID] = 1;
                        if (earliestTime[v] <= offset) {
                            earliestTime[v] = time + offset;
                            if (time != t_end && v != t) {
                                nextFrontier[nextFrontierEnd++] = v;
                            }
                        }
                        else {
                            if (time + offset < earliestTime[v]) {
                                earliestTime[v] = time + offset;
                                VertexID* re = find(nextFrontier, nextFrontier + nextFrontierEnd, v);
                                if (re == nextFrontier + nextFrontierEnd && v != t) {
                                    nextFrontier[nextFrontierEnd++] = v;
                                }
                            }
                        }
                    }
                }
            }
        }

        VertexID* tmp = forwardFrontier;
        forwardFrontier = nextFrontier;
        nextFrontier = tmp;
        forwardFrontierEnd = nextFrontierEnd;

        #ifdef WRITE_STATISTICS
            maxFrontierSize = max(forwardFrontierEnd, maxFrontierSize);
        #endif
    }
}


void STPG::latestTimeBFS() {
    latestTime[t] = offset + t_end + 1;
    backwardFrontier[0] = t;
    backwardFrontierEnd = 1;
    memset(isvisitedEdge, 0, sizeof(EdgeID)*EN);

    while(backwardFrontierEnd > 0) {
        nextFrontierEnd = 0;

        for (VertexID i = 0; i < backwardFrontierEnd; i++) {
            VertexID& u = backwardFrontier[i];

            for (EdgeID edgeLocator = inNeighborsLocator[u]; edgeLocator < inNeighborsLocator[u + 1]; edgeLocator++) {
                Time& time = inNeighbors_time[edgeLocator].time;
                VertexID& v = inNeighbors_time[edgeLocator].neighbor;
                EdgeID& edgeID = outNeighbors_time[edgeLocator].edgeId;

                if (isvisitedEdge[edgeID] == 1 || time < t_begin || v == t || v == s) continue;   
                else if (time > t_end) break;
                else {
                    if (time + offset < latestTime[u]) {
                        isvisitedEdge[edgeID] = 1;
                        if (latestTime[v] <= offset) {
                            latestTime[v] = time + offset;
                            if (time != t_begin) {
                                nextFrontier[nextFrontierEnd++] = v;
                            }
                        }
                        else {
                            if (time + offset > latestTime[v]) {
                                latestTime[v] = time + offset;
                                VertexID* re = find(nextFrontier, nextFrontier + nextFrontierEnd, v);
                                if (re == nextFrontier + nextFrontierEnd) {
                                    nextFrontier[nextFrontierEnd++] = v;
                                }
                            }
                        }
                    }
                }
            }
        }

        VertexID* tmp = backwardFrontier;
        backwardFrontier = nextFrontier;
        nextFrontier = tmp;
        backwardFrontierEnd = nextFrontierEnd;

        #ifdef WRITE_STATISTICS
            maxFrontierSize = max(backwardFrontierEnd, maxFrontierSize);
        #endif
    }
}


void STPG::prunedVertices() {
    for (int i = 0; i < VN; ++i) {
        if (i != s && i != t) {
            if (earliestTime[i] <= offset || latestTime[i] <= offset || earliestTime[i] >= latestTime[i]) {
                dropflag[i] = offset + 1;
            }
        }
    }
}


void STPG::prunedEdges() {
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

    #ifdef WRITE_STATISTICS
        numOfTPAnswers[queryId] = tpgEnd;
    #endif
}


void STPG::forward() {
    
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
            if (from == s) {
                forwardFlag[to] = offset + 2;
                if (!departures[to])
                    departures[to] = time;
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


void STPG::backward() {

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
            if (to == t) {
                backwardFlag[from] = offset + 2;
                if (!arrivals[from])
                    arrivals[from] = time;
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


void STPG::UpperBound() {
    
    upperEnd = 0;
    edgesForVerificationEnd = 0;
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
                
                if (forwardFlag[tpg_graph[i].fromId] == offset + 2 && backwardFlag[tpg_graph[i].toId] == offset + 2 && departures[tpg_graph[i].fromId] < tpg_graph[i].time && arrivals[tpg_graph[i].toId] > tpg_graph[i].time) {
                    results[resultEnd++] = tpg_graph[i].edgeId;
                    hasVerified[tpg_graph[i].edgeId] = offset + 1;
                }

                else if (forwardFlag[tpg_graph[i].fromId] == offset + 2 && departures[tpg_graph[i].fromId] < tpg_graph[i].time) {
                    results[resultEnd++] = tpg_graph[i].edgeId;
                    hasVerified[tpg_graph[i].edgeId] = offset + 1;
                }
                
                else if (backwardFlag[tpg_graph[i].toId] == offset + 2 && arrivals[tpg_graph[i].toId] > tpg_graph[i].time) {
                    results[resultEnd++] = tpg_graph[i].edgeId;
                    hasVerified[tpg_graph[i].edgeId] = offset + 1;
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

    #ifdef WRITE_STATISTICS
        numOfUpperbound[queryId] = upperEnd;
    #endif
}

void STPG::verifyUndeterminedEdge() {                              
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
                    addToResult(from, to, time); 
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
                    addToResult(from, to, time);
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


bool STPG::forwardSearch(VertexID& from, VertexID& to, EdgeID& e, Time& time, Time& fromtime) {
    
    bool ans = false;
    stack[stackEnd++] = to;
    isinStack[to] = offset + 1;
    stacke[stackeEnd++] = e;

    if (to == t) {
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
    }

    stackEnd--;
    stackeEnd--;
    isinStack[to] = offset;

    return false;
}


bool STPG::backwardSearch(VertexID& from, EdgeID& e, Time& time) {
    
    bool ans = false;
    stack[stackEnd++] = from;
    isinStack[from] = offset + 1;
    stacke[stackeEnd++] = e;

    maxStack = max(stackEnd, maxStack);

    if (from == s) {
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
    }

    stackEnd--;
    stackeEnd--;
    isinStack[from] = offset;

    return false;
}


bool STPG::backward2Search(VertexID& from, VertexID& to, EdgeID& e, Time& time, Time& totime) {
    
    bool ans = false;
    stack[stackEnd++] = from;
    isinStack[from] = offset + 1;
    stacke[stackeEnd++] = e;

    if (from == s) {
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
    }

    stackEnd--;
    stackeEnd--;
    isinStack[from] = offset;

    return false;
}


bool STPG::forward2Search(VertexID& to, EdgeID& e, Time& time) {
    
    bool ans = false;
    stack[stackEnd++] = to;
    isinStack[to] = offset + 1;
    stacke[stackeEnd++] = e;

    maxStack = max(stackEnd, maxStack);

    if (to == t) {
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
    }

    stackEnd--;
    stackeEnd--;
    isinStack[to] = offset;

    return false;
}


void STPG::addToResult(VertexID& from, VertexID& to, Time& time) {
    
    int i, j; 
    if (t_end - time < time - t_begin) {
        isinStack[stack[stackEnd - 1]] = offset; 
        isinStack[stack[stackEnd - 2]] = offset;
        for (i = stackeEnd - 3; i >= 3; --i) {
            isinStack[stack[i]] = offset;
            if (stack[i] == from) {
                for (EdgeID edgeLocators = new_backwardLocator[stack[i]]; edgeLocators < new_backwardLocator[stack[i] + 1]; ++edgeLocators) {
                    VertexID& u = new_sortOutEdges[edgeLocators].toId;
                    if (u == to && new_sortOutEdges[edgeLocators].time > edges[stacke[i + 1]].time && new_sortOutEdges[edgeLocators].time < edges[stacke[1]].time) {
                        if (hasVerified[new_sortOutEdges[edgeLocators].edgeId] <= offset) {
                            results[resultEnd++] = new_sortOutEdges[edgeLocators].edgeId;
                            hasVerified[new_sortOutEdges[edgeLocators].edgeId] = offset + 1;
                        }
                    }
                }
                break;
            }
            else {
                for (EdgeID edgeLocators = new_backwardLocator[stack[i]]; edgeLocators < new_backwardLocator[stack[i] + 1]; ++edgeLocators) {
                    VertexID& u = new_sortOutEdges[edgeLocators].toId;
                    if (u == from) {
                        if (u == stack[i - 1] && new_sortOutEdges[edgeLocators].time > edges[stacke[i + 1]].time && new_sortOutEdges[edgeLocators].time < edges[stacke[0]].time) {
                            if (hasVerified[new_sortOutEdges[edgeLocators].edgeId] <= offset) {
                                results[resultEnd++] = new_sortOutEdges[edgeLocators].edgeId;
                                hasVerified[new_sortOutEdges[edgeLocators].edgeId] = offset + 1;
                            }
                        }
                    }
                    else {
                        if (u == stack[i - 1] && new_sortOutEdges[edgeLocators].time > edges[stacke[i + 1]].time && new_sortOutEdges[edgeLocators].time < edges[stacke[i - 1]].time) {
                            if (hasVerified[new_sortOutEdges[edgeLocators].edgeId] <= offset) {
                                results[resultEnd++] = new_sortOutEdges[edgeLocators].edgeId;
                                hasVerified[new_sortOutEdges[edgeLocators].edgeId] = offset + 1;
                            }
                        }
                    }
                }
            }
        }
        
        isinStack[stack[i - 1]] = offset;
        isinStack[stack[i - 2]] = offset;
        
        for (i = 0; i < stackeEnd - 5; ++i) {
            isinStack[stack[i]] = offset;
            if (stack[i + 2] == t) {
                break;
            }
            else {
                for (EdgeID edgeLocators = new_backwardLocator[stack[i]]; edgeLocators < new_backwardLocator[stack[i] + 1]; ++edgeLocators) {
                    VertexID& u = new_sortOutEdges[edgeLocators].toId;
                    if (u == stack[i + 1] && new_sortOutEdges[edgeLocators].time > edges[stacke[i]].time && new_sortOutEdges[edgeLocators].time < edges[stacke[i + 2]].time) {
                        if (hasVerified[new_sortOutEdges[edgeLocators].edgeId] <= offset) {
                            results[resultEnd++] = new_sortOutEdges[edgeLocators].edgeId;
                            hasVerified[new_sortOutEdges[edgeLocators].edgeId] = offset + 1;
                        }
                    }
                }
            }
        }
    }
    else {
        isinStack[stack[0]] = offset;
        for (i = 1; i <= stackEnd - 5; ++i) {
            isinStack[stack[i]] = offset;
            if (stack[i + 1] == s) {
                break;
            }
            else {
                for (EdgeID edgeLocators = new_backwardLocator[stack[i]]; edgeLocators < new_backwardLocator[stack[i] + 1]; ++edgeLocators) {
                    VertexID& u = new_sortOutEdges[edgeLocators].toId;
                    if (u == stack[i - 1] && new_sortOutEdges[edgeLocators].time > edges[stacke[i + 1]].time && new_sortOutEdges[edgeLocators].time < edges[stacke[i - 1]].time) {
                        if (hasVerified[new_sortOutEdges[edgeLocators].edgeId] <= offset) {
                            results[resultEnd++] = new_sortOutEdges[edgeLocators].edgeId;
                            hasVerified[new_sortOutEdges[edgeLocators].edgeId] = offset + 1;
                        }
                    }
                }
            }
        }
        
        isinStack[stack[++i]] = offset;
        
        for (j = i + 1; j <= stackEnd - 3; ++j) {
            isinStack[stack[j]] = offset;
            if (stack[j] == to) {
                for (EdgeID edgeLocators = new_backwardLocator[stack[0]]; edgeLocators < new_backwardLocator[stack[0] + 1]; ++edgeLocators) {
                    VertexID& u = new_sortOutEdges[edgeLocators].toId;
                    if (u == stack[j] && new_sortOutEdges[edgeLocators].time > edges[stacke[1]].time && new_sortOutEdges[edgeLocators].time < edges[stacke[j + 1]].time) {
                        if (hasVerified[new_sortOutEdges[edgeLocators].edgeId] <= offset) {
                            results[resultEnd++] = new_sortOutEdges[edgeLocators].edgeId;
                            hasVerified[new_sortOutEdges[edgeLocators].edgeId] = offset + 1;
                        }
                    }
                }
            }
            
            if (stack[j + 2] == t)
                break;
            
            for (EdgeID edgeLocators = new_backwardLocator[stack[i]]; edgeLocators < new_backwardLocator[stack[i] + 1]; ++edgeLocators) {
                VertexID& u = new_sortOutEdges[edgeLocators].toId;
                if (stack[j] == to) {
                    if (u == stack[j + 1] && new_sortOutEdges[edgeLocators].time > edges[stacke[j]].time && new_sortOutEdges[edgeLocators].time < edges[stacke[j + 2]].time) {
                        if (hasVerified[new_sortOutEdges[edgeLocators].edgeId] <= offset) {
                            results[resultEnd++] = new_sortOutEdges[edgeLocators].edgeId;
                            hasVerified[new_sortOutEdges[edgeLocators].edgeId] = offset + 1;
                        }
                    }
                }
                else {
                    if (u == stack[j + 1] && new_sortOutEdges[edgeLocators].time > edges[stacke[0]].time && new_sortOutEdges[edgeLocators].time < edges[stacke[j + 2]].time) {
                        if (hasVerified[new_sortOutEdges[edgeLocators].edgeId] <= offset) {
                            results[resultEnd++] = new_sortOutEdges[edgeLocators].edgeId;
                            hasVerified[new_sortOutEdges[edgeLocators].edgeId] = offset + 1;
                        }
                    }
                }
            }
        }
        
        isinStack[stack[j + 1]] = offset;
        isinStack[stack[j + 2]] = offset;
    }
}


inline void STPG::refreshMemory() {

    offset += t_end + 1;

    memset(departures, 0, sizeof(VertexID)*VN);
    memset(arrivals, 0, sizeof(VertexID)*VN);
    
    if (offset >= 4294964295) {
        offset = 0;
        
        memset(earliestTime, 0, sizeof(Time)*VN);
        memset(latestTime, 0, sizeof(Time)*VN);
        memset(dropflag, 0, sizeof(VertexID)*VN);

        memset(forwardFlag, 0, sizeof(int)*VN);
        memset(forwardEVLen, 0, sizeof(int)*VN);
        memset(backwardFlag, 0, sizeof(int)*VN);
        memset(backwardEVLen, 0, sizeof(int)*VN);

        memset(hasVerified, 0, sizeof(EdgeID)*EN);
        memset(isinStack, 0, sizeof(VertexID)*VN);
    }

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
    
    #ifdef WRITE_STATISTICS
        maxFrontierSize = 1;
    #endif
}


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


#ifdef WRITE_STATISTICS

    double STPG::getCurrentSpaceCost() {
        double spaceCost = 0;

        spaceCost += sizeof(VertexID)*maxFrontierSize;          // nextFrontier
        spaceCost += sizeof(VertexID)*maxFrontierSize;          // forwardFrontier
        spaceCost += sizeof(VertexID)*maxFrontierSize;          // backwardFrontier
        
        spaceCost += sizeof(Time)*VN;                           // earliestTime
        spaceCost += sizeof(Time)*VN;                           // latestTime
        spaceCost += sizeof(EdgeID)*EN;                         // isvisitedEdge

        spaceCost += sizeof(VertexID)*VN;                       // dropflag
        spaceCost += sizeof(PerEdge)*tpgEnd;                    // tpg_graph
        spaceCost += sizeof(PerEdge)*tpgEnd;                    // sortInEdges
        spaceCost += sizeof(PerEdge)*tpgEnd;                    // sortOutEdges

        spaceCost += sizeof(VertexID)*forwardEVEnd;             // forwardEV
        spaceCost += sizeof(VertexID)*backwardEVEnd;            // backwardEV
        spaceCost += sizeof(int)*forwardLocatorEnd;             // forwardFlag
        spaceCost += sizeof(EdgeID)*forwardLocatorEnd;          // forwardLocator
        spaceCost += sizeof(int)*backwardLocatorEnd;            // backwardFlag
        spaceCost += sizeof(EdgeID)*backwardLocatorEnd;         // backwardLocator
        spaceCost += sizeof(int)*forwardLocatorEnd;             // forwardEVLen
        spaceCost += sizeof(int)*backwardLocatorEnd;            // backwardEVLen
        
        spaceCost += sizeof(VertexID)*VN;                       // departures
        spaceCost += sizeof(VertexID)*VN;                       // arrivals
        spaceCost += sizeof(EdgeID)*edgesForVerificationEnd;    // Undetermined edges
        spaceCost += sizeof(PerEdge)*upperEnd;                  // upper_bound_graph
        spaceCost += sizeof(PerEdge)*upperEnd;                  // new_sortInEdges
        spaceCost += sizeof(PerEdge)*upperEnd;                  // new_sortOutEdges
        spaceCost += sizeof(EdgeID)*new_forwardLocatorEnd;      // new_forwardLocator
        spaceCost += sizeof(EdgeID)*new_backwardLocatorEnd;     // new_backwardLocator

        spaceCost += sizeof(EdgeID)*edgesForVerificationEnd;    // hasVerified
        spaceCost += sizeof(VertexID)*new_forwardLocatorEnd;    // isinStack
        spaceCost += sizeof(VertexID)*maxStack;                 // stack
        spaceCost += sizeof(EdgeID)*maxStack;                   // stacke

        spaceCost += sizeof(EdgeID)*resultEnd;                  // results

        return spaceCost;
    }
    
    void STPG::initStatisticStorage(int queryNumber){
        numOfAnswers = new EdgeID[queryNumber]();
        numOfUpperbound = new EdgeID[queryNumber]();
        numOfTPAnswers = new VertexID[queryNumber]();
        spaceCosts = new double[queryNumber]();

        statisticsFile.open("../" +statisticsPath+extractFilename(queryFilename)+".csv",ios::app);
        statisticsFile << "Space cost (bytes),# TPG edges,# Upper-bound edges,# STPG edges" << endl;
    }

    void STPG::cleanUpStatisticStorage(){
        delete[] numOfAnswers;
        delete[] numOfUpperbound;
        delete[] numOfTPAnswers;
        delete[] spaceCosts;

        statisticsFile.close();
    }

#endif

#endif