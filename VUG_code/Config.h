#ifndef CONFIG_H
#define CONFIG_H

#include <bits/stdc++.h>
using namespace std;

/*  
 * Input and Output  
*/

// whether to write answers to file
#define WRITE_ANSWERS

// answer file relative path
string answerPath = "Results/Answers/";

// whether to write statistics logs to file
#define WRITE_STATISTICS

// statistics file relative path
string statisticsPath = "Results/Statistics/";

// dataset relative path
string datasetPath = "Datasets/";

// log file relative path
string logPath = "Results/Logs.csv";

/*
Variable Type 
*/

// integer type for Vertex IDs (related to graph max size)
typedef unsigned int VertexID;

// integer type for Edge IDs (related to graph max size)
typedef unsigned int EdgeID;

// integer type for Time IDs
typedef unsigned int Time;

#endif