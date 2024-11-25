#include "TSPG.cc"
using namespace std;

int main(int argc, char *argv[]) {    

    // program input parameters
    if(argc < 3) {
        cout << "Usage: ./Run <Graph File> <Query File>" << endl;
        exit(1);
    }
    graphFilename = extractFilename(argv[1]); 
    queryFilename = extractFilename(argv[2]); 

    // basic logs
    logFile.open("../"+logPath, ios::app);
    outputBasicLogs("TSPG");

    // initialize the graph
    Graph* graph = new Graph(("../"+datasetPath+graphFilename).c_str());

    // initialize the queries
    vector<PerQuery> queries;
    loadQueries(("../" + datasetPath + queryFilename).c_str(), queries);

    // STPG
    STPG* method = new STPG(graph);
    method->answerAllQueries(queries);
    method->cleanUp();


    return 0;
}

