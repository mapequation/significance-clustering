#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <random>
#include <functional>
#include <map>
#include <set>
using namespace std;


unsigned stou(char *s);

template <class T>
inline std::string to_string (const T& t){
  std::stringstream ss;
  ss << t;
  return ss.str();
}

std::uniform_real_distribution<double> realdist(0,1);
std::uniform_int_distribution<double> bitdist(0,1);

const string CALL_SYNTAX = "Call: ./sigclu [-s <seed>] [-c <confidencelevel>] [-w <weightsfile>] partitionsfile nodeoutfile moduleoutfile\n"
  "seed: Any positive integer.\n"
  "confidencelevel: The confidence as a fraction, default is 0.95.\n"
  "partitionsfile: Each column represents a partition, the first for the raw partition and the remaining for bootstrap partitions. Row number corresponds to node id.\n"
  "nodeoutfile: 1 or 0 if a node does or does not belong to the significant core of its module.\n"
  "moduleoutfile: moduleId1 moduleId2 means that the significant core of moduleId1 cooccurs with moduleId2 more than a fraction 1 - conf of the samples.\n"
  "weightsfile: One column for weights of each node.  Row number corresponds to node id.\n";

vector<string> tokenize(const string& str,string& delimiters)
{

  vector<string> tokens;

  // skip delimiters at beginning.
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);

  // find first "non-delimiter".
  string::size_type pos = str.find_first_of(delimiters, lastPos);

  while (string::npos != pos || string::npos != lastPos){

    // found a token, add it to the vector.
    tokens.push_back(str.substr(lastPos, pos - lastPos));

    // skip delimiters.  Note the "not_of"
    lastPos = str.find_first_not_of(delimiters, pos);

    // find next "non-delimiter"
    pos = str.find_first_of(delimiters, lastPos);

  }

  return tokens;
}

void find_and_replace(string &str,string searchString,string replaceString){

  string::size_type pos = 0;
  while ( (pos = str.find(searchString, pos)) != string::npos ) {
    str.replace( pos, searchString.size(), replaceString );
    pos++;
  }
  
}

class treeNode{
 public:
  // double exit;
  int moduleId;
  multimap<double,pair<int,string>,greater<double> > members;
  multimap<double,treeNode,greater<double> > nextLevel;
};

void readPartitionsFile(vector<int> &rawPartition,vector<vector<int > > &bootPartitions,ifstream &partitionsFile,int &Nnodes,int &NbootSamples){

  cout << "Reading partitions file " << flush;  

  string line;  
  string buf;

  // Count number of nodes and boot partitions
  getline(partitionsFile,line);
    Nnodes++; // First line corresponds to first node
  istringstream read(line);
  while(read >> buf)
      NbootSamples++;
  NbootSamples--; // The first column is the raw partition
  cout << "with 1 + " << NbootSamples << " partitions " << flush;

  // Count remaining nodes
  while(getline(partitionsFile,line) != 0)
    Nnodes++;
  cout << "of " << Nnodes << " nodes..." << flush;

  rawPartition = vector<int>(Nnodes);
  bootPartitions = vector<vector<int> >(NbootSamples,vector<int>(Nnodes));

  // Restart from beginning of file
  partitionsFile.clear();
  partitionsFile.seekg(0, ios::beg);

  // Read partitions data    
  int nodeNr = 0;
  while(getline(partitionsFile,line) != 0){
    istringstream read(line);
    read >> buf;
    rawPartition[nodeNr] = atoi(buf.c_str());
    int i = 0;
    while(read >> buf){
      bootPartitions[i][nodeNr] = atoi(buf.c_str()); 
      i++;
    }
    nodeNr++;
  }
  partitionsFile.close();
  cout << "done!" << endl;

  // for(int i=0;i<NbootSamples;i++){
  //   for(int j=0;j<Nnodes;j++)
  //     cout << bootPartitions[i][j] << " ";
  //   cout << endl;
  // }

}

void readWeightsFile(vector<double> &weights,ifstream &weightsFile){

  cout << "Reading weight file " << flush;  

  string line;  
  string buf;

  double totWeight = 0.0;
  int nodeNr = 0;
  while(getline(weightsFile,line) != 0){
    istringstream read(line);
    read >> buf;
    weights[nodeNr] = atof(buf.c_str());
    totWeight += weights[nodeNr];
    nodeNr++;
  }

  cout << "with total weight " << totWeight << ", normalizing to 1..." << flush;
  int Nnodes = weights.size();
  if(Nnodes != nodeNr){
    cout << nodeNr << "nodes in weights file, expecting " << Nnodes << " from partitions file, exiting..." << endl;
    abort();
  }

  // Normalize to 1
  for(int i=0;i<nodeNr;i++)
    weights[i] /= totWeight;

  cout << "done!" << endl;
}

void printSignificanceClustering(vector<bool> &significantVec,vector<pair<int,int> > &mergers,multimap<double,treeNode,greater<double> > &treeMap,string nodeOutFileName,string moduleOutFileName){

  int M = treeMap.size();
  vector<int> moduleId = vector<int>(M);
  int i = 0;
  for(multimap<double,treeNode>::iterator it = treeMap.begin(); it != treeMap.end(); it++,i++)
    moduleId[i] = it->second.moduleId;


  ofstream outfile(nodeOutFileName);
  int Nnodes = significantVec.size();
  outfile << "# 1 or 0 if a node does or does not belong to the significant core of its module." << endl;
  outfile << "Significant" << endl;
  for(int i=0;i<Nnodes;i++)
    outfile << significantVec[i] << endl;
  outfile.close();

  outfile.open(moduleOutFileName);
  outfile << "# moduleId1 moduleId2 means that the significant core of moduleId1 cooccurs with moduleId2 more than a fraction 1 - conf of the samples." << endl;
  outfile << "moduleId1 moduleId2" << endl;
  for(vector<pair<int,int> >::iterator it = mergers.begin(); it != mergers.end(); it++)
    outfile << moduleId[it->first] << " " << moduleId[it->second] << endl;
  outfile.close();


}

void generateTreeMap(vector<int> &rawPartition, vector<double> &weights, multimap<double,treeNode,greater<double> > &treeMap){

  int Nnodes = rawPartition.size();
  map<int,double> moduleWeights;
  map<int,vector<int> > moduleMembers;
  for(int i=0;i<Nnodes;i++){
    moduleWeights[rawPartition[i]] += weights[i];
    moduleMembers[rawPartition[i]].push_back(i);
  }

  multimap<double,treeNode,greater<double> >::iterator it_tM;
  map<int,double>::iterator modW_it = moduleWeights.begin();
  map<int,vector<int> >::iterator modM_it = moduleMembers.begin();
  for(; modM_it != moduleMembers.end(); modW_it++,modM_it++){
    int Nmembers = modM_it->second.size();
    treeNode tmp_tN;
    it_tM = treeMap.insert(make_pair(modW_it->second,tmp_tN));
    it_tM->second.moduleId = modM_it->first;
    for(int j=0;j<Nmembers;j++){
      it_tM->second.members.insert(make_pair(weights[modM_it->second[j]],make_pair(modM_it->second[j],to_string(modM_it->second[j]))));
    }
  }


}


// void printSignificantTree(string s,multimap<double,treeNode,greater<double> >::iterator it_tM,ofstream *outfile,vector<bool> &significantVec){
  
//   multimap<double,treeNode,greater<double> >::iterator it;
//   if(it_tM->second.nextLevel.size() > 0){
//     int i=1;
//     for(it = it_tM->second.nextLevel.begin(); it != it_tM->second.nextLevel.end(); it++){
//       string cpy_s(s + to_string(i));
//       printSignificantTree(cpy_s,it,outfile,significantVec);
//       i++;
//     }
//   }
//   else{
//     int i = 1;
//     for(multimap<double,pair<int,string>,greater<double> >::iterator mem = it_tM->second.members.begin(); mem != it_tM->second.members.end(); mem++){
//       string cpy_s(s + (significantVec[mem->second.first] == true ? (":") : (";") ) + to_string(i) + " \"" + mem->second.second + "\" " + to_string(mem->first));
//       (*outfile) << cpy_s << endl;
//       i++;
//     }  
//   }
// }

void findConfCore(multimap<double,treeNode,greater<double> > &treeMap,vector<vector<int > > &bootPartitions,vector<bool> &significantVec,double conf,std::mt19937 &mtrand){
  
  double epsilon = 1.0e-5;

  int Nboots = bootPartitions.size();
  int Nremove = static_cast<int>((1.0-conf)*Nboots+0.5);
  
  int Nnode = bootPartitions[0].size();
  vector<double> size = vector<double>(Nnode);
  
  int M = treeMap.size();
  vector<double> moduleSize = vector<double>(M);
  vector<int> moduleId = vector<int>(M);

  vector<vector<int> > modSortMembers = vector<vector<int> >(M);
  int i = 0;
  
  for(multimap<double,treeNode>::iterator it = treeMap.begin(); it != treeMap.end(); it++){
    moduleId[i] = it->second.moduleId;
    int Nmem = it->second.members.size();
    modSortMembers[i] = vector<int>(Nmem);
    int j=0;
    for(multimap<double,pair<int,string>,greater<double> >::iterator it2 = it->second.members.begin(); it2 != it->second.members.end(); it2++){
      int mem = it2->second.first;
      modSortMembers[i][j] = mem;
      // size[mem] = 1.0*it2->first;
      size[mem] = (1.0-epsilon)*it2->first + epsilon/Nnode;
      moduleSize[i] += size[mem];
      j++;
    }
    i++;
  }
  
  
  cout << endl << "MCMC to maximize significant core of modules, enumerated by size with module id in parenthesis:" << endl;
  for(int i=0;i<M;i++){
    
    cout << "module " << i+1 << " (" << moduleId[i] << "): " << flush;
    
    int N = modSortMembers[i].size();
    std::uniform_int_distribution<double> intNdist(0,N-1);
    vector<bool> confState = vector<bool>(N);
    vector<bool> maxConfState = vector<bool>(N);
    double maxScore = -1.0;
    
    double confSize = 0.0;
    int confN = 0;
    double maxConfSize = 0.0;
    int maxConfN = 0;
    double score = 0.0;
    int penalty = 0;
    
    double pW = 10.0*moduleSize[i];
    
    if(N != 1){
      
      int maxModNr = 0;
      for(int j=0;j<Nboots;j++)
        for(int k=0;k<N;k++)
          if(bootPartitions[j][modSortMembers[i][k]] > maxModNr )
            maxModNr = bootPartitions[j][modSortMembers[i][k]];
      maxModNr++;
      
      // Initiate weights of module assignments
      
      // Keep track of the size order of modules
      vector<multimap<double,int,greater<double> > > sortModSizes = vector<multimap<double,int,greater<double> > >(Nboots);
      // Keep track of which modules that are included
      vector<vector<pair<int,multimap<double,int,greater<double> >::iterator > > > mapModSizes = vector<vector<pair<int,multimap<double,int,greater<double> >::iterator > > >(Nboots);
      for(int j=0;j<Nboots;j++){
        mapModSizes[j] = vector<pair<int,multimap<double,int,greater<double> >::iterator > >(maxModNr);
        for(int k=0;k<maxModNr;k++){
          mapModSizes[j][k].first = 0;
        }
      }
      
      // Randomized start
      for(int j=0;j<N;j++){
        if(bitdist(mtrand)){
          confState[j] = true;
          confSize += size[modSortMembers[i][j]];
          confN++;
          
          for(int k=0;k<Nboots;k++){
            
            int modNr = bootPartitions[k][modSortMembers[i][j]];
            if(mapModSizes[k][modNr].first == 0){
              double newSize = size[modSortMembers[i][j]];
              multimap<double,int,greater<double> >::iterator it = sortModSizes[k].insert(make_pair(newSize,modNr));
              mapModSizes[k][modNr].second = it;
              mapModSizes[k][modNr].first++;
            }
            else{
              double newSize = size[modSortMembers[i][j]] + mapModSizes[k][modNr].second->first;
              multimap<double,int,greater<double> >::iterator it = sortModSizes[k].insert(mapModSizes[k][modNr].second,make_pair(newSize,modNr));
              sortModSizes[k].erase(mapModSizes[k][modNr].second);
              mapModSizes[k][modNr].second = it;
              mapModSizes[k][modNr].first++;
            }
            
          }
        }
        else
          confState[j] = false;
      }
      
      multimap<double,pair<int,int> > scoreRank;
      score = 0.0;
      // Calculate penalty
      for(int j=0;j<Nboots;j++){
        double tmpScore = 0.0;
        int tmpPenalty = 0;
        if(!sortModSizes[j].empty()){
            int modNr = sortModSizes[j].begin()->second; 
            tmpScore = sortModSizes[j].begin()->first;
            tmpPenalty = confN - mapModSizes[j][modNr].first; //penalty is the number of nodes not in biggest field
        }
        scoreRank.insert(make_pair(tmpScore-pW*tmpPenalty,make_pair(tmpPenalty,j)));
        score += tmpScore;
        penalty += tmpPenalty;
      }
      
      // Remove worst results
      multimap<double,pair<int,int> >::iterator it = scoreRank.begin();
      for(int j=0;j<Nremove;j++){
        int bootNr = it->second.second;
        double tmpScore = 0.0;
        int tmpPenalty = 0;
        if(!sortModSizes[bootNr].empty()){
          tmpScore = sortModSizes[bootNr].begin()->first;
          tmpPenalty = it->second.first;
        }
        score -= tmpScore;
        penalty -= tmpPenalty;
        it++;
      }
      
      
      //Monte Carlo to maximize confident size
      int Niter = static_cast<int>(pow(1.0*N,1.0));
      if(Niter < 100)
        Niter = 100;
      int attempts = 0;
      int switches = 0;
      bool search = true;
      while(search){
        
        double T = 1.0;
        
        do{
          
          attempts = 0;
          switches = 0;
          for(int j=0;j<Niter;j++){
            
            multimap<double,pair<int,int> >().swap(scoreRank);
            double newConfSize = confSize;
            int newConfN = confN;
            double newScore = 0.0;
            int newPenalty = 0;
            
            int flip = intNdist(mtrand);
            int nodeNr = modSortMembers[i][flip];
            
            if(confState[flip]){ // Remove one node from confident subset
              newConfSize -= size[nodeNr];
              newConfN--;
              for(int k=0;k<Nboots;k++){
                int modNr = bootPartitions[k][nodeNr];
                double tmpScore = sortModSizes[k].begin()->first;
                int tmpPenalty = newConfN - mapModSizes[k][sortModSizes[k].begin()->second].first;
                if(mapModSizes[k][modNr].second == sortModSizes[k].begin()){
                  tmpScore -= size[nodeNr];
                  tmpPenalty = newConfN - (mapModSizes[k][modNr].first-1);
                  if(sortModSizes[k].size() > 1){
                    multimap<double,int,greater<double> >::iterator it = sortModSizes[k].begin();
                    it++;
                    if(it->first > tmpScore){ // Check if second in ranking is larger
                      tmpScore = it->first;
                      tmpPenalty = newConfN - mapModSizes[k][it->second].first;
                    }
                  }
                }
                scoreRank.insert(make_pair(tmpScore-pW*tmpPenalty,make_pair(tmpPenalty,k)));
                newScore += tmpScore;
                newPenalty += tmpPenalty;
              }
            }
            else{ // Add one node to confident subset
              newConfSize += size[nodeNr];
              newConfN++;
              for(int k=0;k<Nboots;k++){
                
                int modNr = bootPartitions[k][nodeNr];
                double tmpScore = 0.0;
                int tmpPenalty = 0;
                if(sortModSizes[k].empty()){ // No nodes in confident subset
                  tmpScore = size[nodeNr];
                  tmpPenalty = 0;
                }
                else{
                  
                  if(mapModSizes[k][modNr].first == 0){ // First confident node in module
                    if(size[nodeNr] > sortModSizes[k].begin()->first){
                      tmpScore = size[nodeNr];
                      tmpPenalty = newConfN - 1;
                    }
                    else{
                      tmpScore = sortModSizes[k].begin()->first;
                      tmpPenalty = newConfN - mapModSizes[k][sortModSizes[k].begin()->second].first;
                    }
                  }
                  else{ // Not first confident node in module
                    
                    if(mapModSizes[k][modNr].second == sortModSizes[k].begin()){
                      tmpScore = sortModSizes[k].begin()->first + size[nodeNr];
                      tmpPenalty = newConfN - (mapModSizes[k][sortModSizes[k].begin()->second].first+1);
                    }
                    else{
                      
                      if(mapModSizes[k][modNr].second->first + size[nodeNr] > sortModSizes[k].begin()->first){
                        tmpScore = mapModSizes[k][modNr].second->first + size[nodeNr];
                        tmpPenalty = newConfN - (mapModSizes[k][modNr].first+1);
                      }
                      else{
                        tmpScore = sortModSizes[k].begin()->first;
                        tmpPenalty = newConfN - mapModSizes[k][sortModSizes[k].begin()->second].first;
                      }
                      
                    }
                    
                  }
                }
                
                scoreRank.insert(make_pair(tmpScore-pW*tmpPenalty,make_pair(tmpPenalty,k)));
                newScore += tmpScore;
                newPenalty += tmpPenalty;
              }
            }
            
           // Remove worst results
           multimap<double,pair<int,int> >::iterator it = scoreRank.begin();
           for(int j=0;j<Nremove;j++){
             int bootNr = it->second.second;
             if(!sortModSizes[bootNr].empty()){
               double tmpScore = sortModSizes[bootNr].begin()->first;
               int tmpPenalty = it->second.first;
               newScore -= tmpScore;
               newPenalty -= tmpPenalty;
             }
             
             it++;
           }
      
            if(exp(((newScore-pW*newPenalty)-(score-pW*penalty))/T) > realdist(mtrand)){
              
              // Update data structures
              if(confState[flip]){ // Remove one node from confident subset
                
                for(int k=0;k<Nboots;k++){
                  int modNr = bootPartitions[k][nodeNr];
                  mapModSizes[k][modNr].first--;
                  if(mapModSizes[k][modNr].first == 0){ // Remove last confident node from module
                    sortModSizes[k].erase(mapModSizes[k][modNr].second);
                  }
                  else{
                    double newSize = mapModSizes[k][modNr].second->first - size[nodeNr];
                    multimap<double,int,greater<double> >::iterator it = sortModSizes[k].insert(mapModSizes[k][modNr].second,make_pair(newSize,modNr));
                    sortModSizes[k].erase(mapModSizes[k][modNr].second);
                    mapModSizes[k][modNr].second = it;
                  }
                  
                }
              }
              else{ // Add one node from confident subset
                
                for(int k=0;k<Nboots;k++){
                  
                  int modNr = bootPartitions[k][nodeNr];
                  if(mapModSizes[k][modNr].first == 0){
                    double newSize = size[nodeNr];
                    multimap<double,int,greater<double> >::iterator it = sortModSizes[k].insert(make_pair(newSize,modNr));
                    mapModSizes[k][modNr].second = it;
                    mapModSizes[k][modNr].first++;
                  }
                  else{
                    double newSize = size[nodeNr] + mapModSizes[k][modNr].second->first;
                    multimap<double,int,greater<double> >::iterator it = sortModSizes[k].insert(mapModSizes[k][modNr].second,make_pair(newSize,modNr));
                    sortModSizes[k].erase(mapModSizes[k][modNr].second);
                    mapModSizes[k][modNr].second = it;
                    mapModSizes[k][modNr].first++;
                  }
                  
                }
                
              }
              
              confSize = newConfSize;
              confN = newConfN;
              confState[flip] = !confState[flip];
              penalty = newPenalty;
              score = newScore;
              switches++;
            }
            attempts++;
            
            if(penalty == 0 && score > maxScore){
              for(int k=0;k<N;k++)
                maxConfState[k] = confState[k];
              maxScore = score;
              maxConfSize = confSize;
              maxConfN = confN;
            }
            
          }
          T *= 0.99;
          //cout << T << " " <<  1.0*switches/attempts << " " << score << " " << penalty << "    " << confSize << " " << confN << endl;
        } while(switches > 0);
        
        if(maxScore > 0.0)
          search = false;
        
      }
      
    }
    else{
      maxConfState[0] = true;
      maxConfSize = moduleSize[i];
      maxConfN = 1;
    }
    
    
    cout << maxConfN << "/" << N << " confident nodes and " << maxConfSize << "/" << moduleSize[i] << " (" << 100*maxConfSize/moduleSize[i] << " percent) of the flow." << endl;
    
    for(int j=0;j<N;j++)
      significantVec[modSortMembers[i][j]] = maxConfState[j];
    
  }
  
  //  i = 0;
  //  for(multimap<double,treeNode>::reverse_iterator it = treeMap.rbegin(); it != treeMap.rend(); it++){
  //    int j=0;
  //    for(multimap<double,treeNode>::reverse_iterator it2 = it->second.nextLevel.rbegin(); it2 != it->second.nextLevel.rend(); it2++){
  //      it2->second.significant = significantVec[modSortMembers[i][j]];
  //      j++;
  //    }
  //    i++;
  //  }
  
}

void findConfModules(multimap<double,treeNode,greater<double> > &treeMap,vector<vector<int > > &bootPartitions,vector<bool> &significantVec,vector<pair<int,int> > &mergers,double conf){
  
  
  int M = treeMap.size();
  int N = bootPartitions[0].size();
  int Nboots = bootPartitions.size();
  vector<int> moduleId = vector<int>(M);


  // Calculate total size of confident journals in field
  vector<vector<int> > significantNodes = vector<vector<int> >(N);
  vector<double> confFieldSize = vector<double>(M,0.0);
  vector<multimap<int,pair<int,vector<int> >,greater<int> > > coExist = vector<multimap<int,pair<int,vector<int> >,greater<int> > >(M);
  
  

  int clusterNr = 0;
  for(multimap<double,treeNode,greater<double> >::iterator it = treeMap.begin();  it != treeMap.end(); it++){
    moduleId[clusterNr] = it->second.moduleId;
    for(multimap<double,pair<int,string>,greater<double> >::iterator it2 = it->second.members.begin(); it2 != it->second.members.end(); it2++){
      if(significantVec[it2->second.first])
        significantNodes[clusterNr].push_back(it2->second.first);
    }
    clusterNr++;
  }
  
  vector<vector<int> > coexistCount = vector<vector<int> >(M);
  for(int i=0;i<M;i++)
    coexistCount[i] = vector<int>(Nboots,0);
  
  cout << endl << "Now calculate number of times two modules are clustered together, enumerated by size with module id in parenthesis:" << endl;
  int i=0;
  for(multimap<double,treeNode,greater<double> >::iterator it1 = treeMap.begin();  it1 != treeMap.end(); it1++){ // i
    int j=i;
    for(multimap<double,treeNode,greater<double> >::iterator it2 = it1;  it2 != treeMap.end(); it2++){ // j
      
      if(it1 != it2){
        
        int coEx = 0;
        for(int k=0;k<Nboots;k++){
          int modNr = bootPartitions[k][significantNodes[j][0]];
          bool joined = true;
          int iNnode = significantNodes[i].size();
          int jNnode = significantNodes[j].size();
          
          for(int l=0;l<jNnode;l++){
            if(bootPartitions[k][significantNodes[j][l]] != modNr){
              joined = false;
              break;
            }
          }
          
          if(joined){
            for(int l=0;l<iNnode;l++){
              if(bootPartitions[k][significantNodes[i][l]] != modNr){
                joined = false;
                break;
              }
            }
          }
          
          if(joined){
            coEx++;
            coexistCount[i][k]++;
            coexistCount[j][k]++;
            
            multimap<int,pair<int,vector<int> >,greater<int> >::iterator it = coExist[i].find(j);
            if(it != coExist[i].end()){
              it->second.first++;
              it->second.second.push_back(k);
            }
            else{
              vector<int> tmp;
              tmp.push_back(k);
              coExist[i].insert(make_pair(j,make_pair(1,tmp)));
            }
            it = coExist[j].find(i);
            if(it != coExist[j].end()){
              it->second.first++;
              it->second.second.push_back(k);
            }
            else{
              vector<int> tmp;
              tmp.push_back(k);
              it = coExist[j].insert(make_pair(i,make_pair(1,tmp)));
            }
          }
          
        }
        
        
      }
      
      j++;
    }
    
    // Re-sort co-exist strcuture
    multimap<int,pair<int,vector<int> >,greater<int> > tmp = coExist[i];
    multimap<int,pair<int,vector<int> >,greater<int> >().swap(coExist[i]);
    for(multimap<int,pair<int,vector<int> >,greater<int> >::iterator it = tmp.begin(); it != tmp.end(); it++){
      coExist[i].insert(make_pair(it->second.first,make_pair(it->first,it->second.second)));
    }
    
    i++;
    
  }
  
  vector<int> mergeVec = vector<int>(M);
  i=0;
  for(multimap<double,treeNode,greater<double> >::iterator it1 = treeMap.begin();  it1 != treeMap.end(); it1++){
    int singleN = 0;
    for(int k=0;k<Nboots;k++)
      if(coexistCount[i][k] == 0)
        singleN++;
    cout << "Module " << i+1 << " (" << moduleId[i] << ") is standalone " << singleN << "/" << Nboots << " times";
    if(singleN == Nboots){
      cout << "." << endl;
    }
    else{
      cout << " and clustered together with: ";
      for(multimap<int,pair<int,vector<int> >,greater<int> >::iterator it2 = coExist[i].begin(); it2 != coExist[i].end(); it2++)
        cout << it2->second.first+1 << " (" << moduleId[it2->second.first] << ") " << it2->first << " times, ";
      cout << endl;
    }
    i++;
  }
  
  // Find merges
  for(int i=M-1;i>=0;i--){
    
    bool searchMerge = true;
    mergeVec[i] = i;
    while(searchMerge){
      int singleN = 0;
      for(int k=0;k<Nboots;k++)
        if(coexistCount[i][k] == 0)
          singleN++;
      if(1.0*singleN/Nboots < conf){
        int mergeWith = coExist[i].begin()->second.first;
        
        for(int j=0;j<coExist[i].begin()->first;j++){
          coexistCount[i][coExist[i].begin()->second.second[j]]--;
          coexistCount[mergeWith][coExist[i].begin()->second.second[j]]--;
        }
        coExist[i].erase(coExist[i].begin());
        
        for(multimap<int,pair<int,vector<int> >,greater<int> >::iterator it = coExist[mergeWith].begin(); it != coExist[mergeWith].end(); it++){
          if(it->second.first == i){
            coExist[mergeWith].erase(it);
            break;
          }
          
        }
        
        if(mergeWith < i){
          mergeVec[i] = mergeWith;
          mergers.push_back(make_pair(i,mergeWith));
          searchMerge = false;
        }
        
      }
      else{
        searchMerge = false;
      }
    }
  }
  
  if(mergers.size() > 0){
    cout << endl << "Module associations of modules that are not significantly standalone:" << endl;
    
    for(vector<pair<int,int> >::iterator it = mergers.begin(); it != mergers.end(); it++)
      cout << it->first+1 << " (" << moduleId[it->first] << ") -> " << it->second+1 << " (" << moduleId[it->second] << ")" << endl;
  }
  else{
    cout << "All modules are significantly standalone." << endl;
  }

}

