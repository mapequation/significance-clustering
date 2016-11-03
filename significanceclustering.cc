#include "significanceclustering.h"

using namespace std;
using std::cout;
using std::cin;
using std::endl;

unsigned stou(char *s){
  return strtoul(s,(char **)NULL,10);
}


// Call: trade <seed> <Ntries>
int main(int argc,char *argv[]){

  if( argc == 1){
    cout  << CALL_SYNTAX;
    exit(-1);
  }
 
  unsigned int seed = 1234;
  double conf = 0.95;
  string partitionsFileName = "noname";
  string weightsFileName = "noname";
  string nodeOutFileName = "noname";
  string moduleOutFileName = "noname";

  int argNr = 1;
  while(argNr < argc){
    if(to_string(argv[argNr]) == "-h"){
      cout << CALL_SYNTAX;
      exit(-1);
    }
    else if(to_string(argv[argNr]) == "-s"){
      argNr++;
      seed = atoi(argv[argNr]);
      argNr++;
    }
    else if(to_string(argv[argNr]) == "-c"){
      argNr++;
      conf = atof(argv[argNr]);
      argNr++;
    }
    else if(to_string(argv[argNr]) == "-w"){
      argNr++;
      weightsFileName = to_string(argv[argNr]);
      argNr++;
    }
    else{

      if(argv[argNr][0] == '-'){
        cout << "Unknown command: " << to_string(argv[argNr]) << endl << CALL_SYNTAX;
        exit(-1);
      }
      else if(partitionsFileName == "noname"){
        partitionsFileName = to_string(argv[argNr]);
        argNr++;
      }
      else if(nodeOutFileName == "noname"){
        nodeOutFileName = to_string(argv[argNr]);
        argNr++;
      }
      else if(moduleOutFileName == "noname"){
        moduleOutFileName = to_string(argv[argNr]);
        argNr++;
      }
      else{
        cout << "Too many parameters without flags." << endl << CALL_SYNTAX;
        exit(-1);
      }
    }
  }

  if(partitionsFileName == "noname"){
    cout << "No partitions file provided." << endl << CALL_SYNTAX;
    exit(-1);
  }
  else if(nodeOutFileName == "noname"){
     cout << "No node out file provided." << endl << CALL_SYNTAX;
     exit(-1);
  }
  else if(moduleOutFileName == "noname"){
     cout << "No module out file provided." << endl << CALL_SYNTAX;
     exit(-1);
  }

  cout << "Setup:" << endl;
  cout << "-->Using seed: " << seed << endl;
  cout << "-->Confidence level: " << conf << endl;
  cout << "-->Reading from partitions file: " << partitionsFileName << endl;
  if(weightsFileName == "noname")
    cout << "-->No weights file provided. Using uniform weights of nodes in signifificance clustering." << endl;
  else
    cout << "-->Using weights of nodes in significance clustering from file: " << weightsFileName << endl;
  cout << "-->Writing results to files: " << nodeOutFileName << " and " << moduleOutFileName << endl;
  
  // int Nsamples = 0;
  vector<int > rawPartition;
  vector<vector<int > > bootPartitions;

  int Nnodes = 0;
  int NbootSamples = 0;

  std::mt19937 mtrand(seed);
  ifstream partitionsFile(partitionsFileName);

  // Read partitions file
  readPartitionsFile(rawPartition,bootPartitions,partitionsFile,Nnodes,NbootSamples);

  vector<double> weights(Nnodes,1.0/Nnodes);
  if(weightsFileName != "noname"){
    // Read wights from file if provided
    ifstream weightsFile(weightsFileName);
    readWeightsFile(weights,weightsFile);
  }
  
  // Store modules by weight and nodes in modules by weight
  multimap<double,treeNode,greater<double> > treeMap;
  generateTreeMap(rawPartition,weights,treeMap);

  // Calculate significance clusters
  vector<bool> significantVec = vector<bool>(Nnodes);
  findConfCore(treeMap,bootPartitions,significantVec,conf,mtrand);
  vector<pair<int,int> > mergers;
  findConfModules(treeMap,bootPartitions,significantVec,mergers,conf);

  printSignificanceClustering(significantVec,mergers,treeMap,nodeOutFileName,moduleOutFileName);

}
