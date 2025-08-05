#include "../../GraphUtils/Graph.cc"
#include<cstdlib>
using namespace std;

int main(int argc, char *argv[]) {
    
    if(argc < 3) {
        cout << "Usage: ./GenerateQueries <Graph File> <number of queries>" << endl;
        exit(1);
    }
    graphFilename = extractFilename(argv[1]); 
    int numOfQueries = stoi(argv[2]);

    Graph* graph = new Graph(("../"+graphFilename).c_str());
    VertexID VN = graph->VN;
    PerNeighbor* outNeighbors = graph->outNeighbors_time;
    EdgeID* outNeighborsLocator = graph->outNeighborsLocator;

    vector<PerQuery> queries(numOfQueries);
    int currentCount = 0;
    VertexID visitedVerticesEnd, forwardFrontierEnd, nextFrontierEnd;

	srand(2025);

    cout << "Generating queries ..." << endl;
    vector<Time> attend_time(VN, 0);
    vector<VertexID> visitedVertices(VN);
    vector<VertexID> forwardFrontier(VN);
    vector<VertexID> nextFrontier(VN);
    Time begin_time, end_time;

    while (currentCount < numOfQueries) {
       
        VertexID s = rand()%VN;
        if (outNeighborsLocator[s] == outNeighborsLocator[s + 1])
            continue;
        
        visitedVerticesEnd = 0;
        std::fill(attend_time.begin(), attend_time.end(), 0);

        begin_time = outNeighbors[outNeighborsLocator[s] + rand()%(outNeighborsLocator[s + 1] - outNeighborsLocator[s])].time;
        end_time = begin_time + 10; // end_time = begin_time + theta

        if (end_time > 86399)   // end_time > |T|
            continue;

        VertexID t;

        attend_time[s] = 0;
        forwardFrontier[0] = s;
        forwardFrontierEnd = 1;

        while (forwardFrontierEnd > 0) {
            nextFrontierEnd = 0;
            for (VertexID i = 0; i < forwardFrontierEnd; i++) {
                VertexID& u = forwardFrontier[i];
                for (EdgeID edgeLocator = outNeighborsLocator[u]; edgeLocator < outNeighborsLocator[u + 1]; edgeLocator++) {
                    Time& time = outNeighbors[edgeLocator].time;
                    VertexID& v = outNeighbors[edgeLocator].neighbor;
                    if (time < begin_time || v == s || u == v) continue;
                    else if (time > end_time) break;
                    else {
                        if (time > attend_time[u]) {
                            if (attend_time[v] == 0) {
                                attend_time[v] = time;
                                visitedVertices[visitedVerticesEnd] = v;
                                visitedVerticesEnd++;
                                if (time != end_time) {
                                    auto end = nextFrontier.begin() + nextFrontierEnd;
                                    if (std::find(nextFrontier.begin(), end, v) == end) {
                                        nextFrontier[nextFrontierEnd] = v;
                                        nextFrontierEnd++;
                                    }
                                }
                            }
                            else {
                                if (time < attend_time[v]) {
                                    attend_time[v] = time;
                                    auto end = nextFrontier.begin() + nextFrontierEnd;
                                    if (std::find(nextFrontier.begin(), end, v) == end) {
                                        nextFrontier[nextFrontierEnd] = v;
                                        nextFrontierEnd++;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            forwardFrontier.swap(nextFrontier);
            forwardFrontierEnd = nextFrontierEnd;
        }
        
        if (!visitedVerticesEnd)    // ensure s can reach t
            continue;

        t = visitedVertices[rand()%visitedVerticesEnd];

        PerQuery quer;
        quer.source = s, quer.target = t, quer.star = begin_time, quer.end = end_time;
        queries[currentCount] = quer;

        currentCount++;    
    }

    ofstream queryFile;
    queryFile.open("../"+graphFilename+".query");
    for (auto& query : queries)
        queryFile << query.source << " " << query.target << " " << query.star << " " << query.end << endl;
    queryFile.close();

    return 0;
}
