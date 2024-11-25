#ifndef UTILS_H
#define UTILS_H

#include "../Config.h"

#include <bits/stdc++.h>
#include <sys/time.h>
using namespace std;

// default graph
string graphFilename = "TestGraph.txt";

// default query file
string queryFilename = "TestGraph.query";

// default time interval
short timestart = 2;
short timeend = 7;


#if (defined _WIN32) || (defined _WIN64)
    #define WINDOWS 1
#else
    #define MAX 4294964295
#endif


#define min(a, b) ((a)<(b)?(a):(b))
#define max(a, b) ((a)>(b)?(a):(b))
#define abs(a) ((a)>0?(a):(-(a)))
int getCeil(double value) {
    return int(ceil(value)+0.1);
}
int getFloor(double value) {
    return int(floor(value)+0.1);
}

// get current time in millisecond
inline double getCurrentTimeInMs() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec*1000 + tv.tv_usec * 1e-3;
}

// get current time for writing logs
string getCurrentLogTime( )
{
    time_t now = time(0);
    tm *ltm = localtime(&now);
    return  to_string(1900 + ltm->tm_year) + "-" + to_string(1 + ltm->tm_mon) + "-" + to_string(ltm->tm_mday) + " " +
            to_string(ltm->tm_hour) + ":" + to_string(ltm->tm_min) + ":" + to_string(ltm->tm_sec);
}

// number to string
string str(float value, short fixed=2) {
    string result = to_string(value);
    if (abs(value-int(value))<1e-6)
        result = to_string(int(value));
    if (result.find(".")!=string::npos) 
        return result.substr(0, result.find(".") + fixed + 1);
    return result;
}

// extract filename
string extractFilename(string s) {
    short startPos=0, endPos=s.length();
    if (s.rfind("/")!=string::npos) 
        startPos = s.rfind("/")+1;
    return s.substr(startPos, endPos);
}

// output file stream
ofstream logFile, resultFile, statisticsFile, File;

// basic logs output
void outputBasicLogs(string methodName) {

    // to console
    printf("################ %s ################\n", methodName.c_str());
    printf("- Dataset: %s\n", graphFilename.c_str());
    printf("- Query File: %s\n", queryFilename.c_str());
    printf("- Start Time: %s\n", getCurrentLogTime().c_str());

    // to log file
    logFile<<methodName<<","<<getCurrentLogTime()<<","<<graphFilename<<","<<queryFilename<<",";
}


// get string of paramters
string getParaString() {
    return "-["+str(timestart)+","+str(timeend)+"]";
}

// storing each query
struct PerQuery {
    VertexID source, target;
    Time star, end;
};

// load query file
void loadQueries(const char* queryFilename, vector<PerQuery>& queries) {
    printf("Loading query file ...\n");
    VertexID fromId, toId;
    Time star, end;
    FILE* f = fopen(queryFilename, "r");
    while (fscanf(f, "%d %d %d %d", &fromId, &toId, &star, &end) != EOF) 
        queries.push_back({fromId, toId, star, end});
    printf("- Finish. %d queries loaded\n", queries.size());
    logFile<<queries.size()<<",";
    fclose(f);
}

// common variables
unsigned int offset = 0;
VertexID *forwardFrontier, forwardFrontierEnd, *backwardFrontier, backwardFrontierEnd, *nextFrontier, nextFrontierEnd;
Time *earliestTime, *latestTime;

#endif

