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
  
  int Nsamples = 0;
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

  
  // // Print significance map in .smap format for the Map Generator at www.mapequation.org
  // oss.str("");
  // oss << networkName << ".smap";
  // outfile.open(oss.str().c_str());
  // outfile << "# modules: " << Nmod << endl;
  // outfile << "# modulelinks: " << sortedLinks.size() << endl;
  // outfile << "# nodes: " << Nnode << endl;
  // outfile << "# links: " << network.Nlinks << endl;
  // outfile << "# codelength: " << greedy->codeLength/log(2.0) << endl;
  // outfile << "*Directed" << endl;
  // outfile << "*Modules " << Nmod << endl;
  // k = 0;
  // for(multimap<double,treeNode,greater<double> >::iterator it = treeMap.begin(); it != treeMap.end(); it++){
  //   outfile << k+1 << " \"" << it->second.members.begin()->second.second << ",...\" " << it->first << " " << it->second.exit << endl;
  //   k++;
  // }
  // outfile << "*Insignificants " << mergers.size() << endl;
  // for(vector<pair<int,int> >::iterator it = mergers.begin(); it != mergers.end(); it++)
  //   outfile << it->first+1 << ">" << it->second+1 << endl;
  // outfile << "*Nodes " << Nnode << endl;
  // k = 1;
  // for(multimap<double,treeNode,greater<double> >::iterator it = treeMap.begin(); it != treeMap.end(); it++){
  //   string s;
  //   s.append(to_string(k));
  //   printSignificantTree(s,it,&outfile,significantVec);
  //   k++;
  // }
  // outfile << "*Links " << sortedLinks.size() << endl;
  // for(multimap<double,pair<int,int>,greater<double> >::iterator it = sortedLinks.begin();it != sortedLinks.end();it++)   
  //   outfile << it->second.first << " " << it->second.second << " " << 1.0*it->first << endl;
  // outfile.close();






  // // Sort network
  // multimap<double,int,greater<int> > sortedNodes;
  // for(map<int,double>::iterator it = degree.begin(); it != degree.end(); it++)
  //   sortedNodes.insert(make_pair(it->second,it->first));

  // vector<int> toNewIndex(Nnodes);
  // vector<int> toOldIndex(Nnodes);
  // vector<double> cumDegree(Nnodes);
  // vector<double> normDegree(Nnodes);
  // double cum = 0.0;
  // int i = 0;
  // for(multimap<double,int,greater<int> >::iterator it = sortedNodes.begin(); it != sortedNodes.end(); it++){
  //   cum += it->first/totalDegree;
  //   cumDegree[i] = cum;
  //   normDegree[i] = it->first/totalDegree;
  //   toNewIndex[it->second] = i;
  //   toOldIndex[i] = it->second;
  //   i++;
  // }

  // vector<vector<double> > cumLinkWeights(Nnodes);
  // vector<vector<int> > cumLinks(Nnodes);
  // for(int i=0;i<Nnodes;i++){
  //   int oldIndex = toOldIndex[i];
  //   int Nlinks = net[oldIndex].size();
  //   cumLinkWeights[i] = vector<double>(Nlinks);
  //   cumLinks[i] = vector<int>(Nlinks);
  //   cum = 0.0;
  //   int j = 0;
  //   for(multimap<double,int,greater<double> >::iterator it = net[oldIndex].begin(); it != net[oldIndex].end(); it++){
  //     cum += it->first/degree[oldIndex];
  //     cumLinkWeights[i][j] = cum;
  //     cumLinks[i][j] = toNewIndex[it->second];
  //     j++;
  //   }
  // }


  // string oufilename = linklist;
  // size_t lastPeriod = oufilename.find_last_of(".");
  // if(lastPeriod != string::npos)
  //   oufilename = string(oufilename.begin(),oufilename.begin() + lastPeriod);
  // oufilename += "_sampledMarkovEntropy.txt";

  // ofstream outfile;
  // outfile.open(oufilename);
  // outfile << "MarkovTime MarkovEntropy" << endl;
  // for(int i=0;i<Nsamples;i++){
  //   double time = pow(10.0,log10(mintime) + 1.0*i*(log10(maxtime)-log10(mintime))/(Nsamples-1));
  //   double h = markovEntropy(mtrand,time,Nnodes,cumDegree,normDegree,cumLinkWeights,cumLinks);
  //   cout << time << " " << h << endl;
  //   outfile << time << " " << h << endl;
  // } 

  // outfile.close();

}
