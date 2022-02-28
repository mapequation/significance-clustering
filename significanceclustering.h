#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <vector>

using namespace std;

template <class T>
inline std::string to_string(const T &t) {
  std::stringstream ss;
  ss << t;
  return ss.str();
}

std::uniform_real_distribution<double> realdist(0, 1);
std::uniform_int_distribution<int> bitdist(0, 1);

class TreeNode {
public:
  int moduleId{};
  multimap<double, pair<int, string>, greater<double>> members;
};

void readPartitionsFile(vector<int> &rawPartition,
                        vector<vector<int>> &bootPartitions,
                        ifstream &partitionsFile,
                        int &Nnodes,
                        int &NbootSamples) {
  cout << "Reading partitions file " << flush;

  string line;
  string buf;

  // Count number of nodes and boot partitions
  getline(partitionsFile, line);
  Nnodes++; // First line corresponds to first node
  {         // braces to avoid shadowing "read"
    istringstream read(line);
    while (read >> buf)
      NbootSamples++;
  }
  NbootSamples--; // The first column is the raw partition
  cout << "with 1 + " << NbootSamples << " partitions " << flush;

  // Count remaining nodes
  while (getline(partitionsFile, line)) {
    Nnodes++;
  }
  cout << "of " << Nnodes << " nodes..." << flush;

  rawPartition = vector<int>(Nnodes);
  bootPartitions = vector<vector<int>>(NbootSamples, vector<int>(Nnodes));

  // Restart from beginning of file
  partitionsFile.clear();
  partitionsFile.seekg(0, ios::beg);

  // Read partitions data
  int nodeNr = 0;
  while (getline(partitionsFile, line)) {
    istringstream read(line);
    read >> buf;
    rawPartition[nodeNr] = atoi(buf.c_str());
    int i = 0;
    while (read >> buf) {
      bootPartitions[i][nodeNr] = atoi(buf.c_str());
      i++;
    }
    nodeNr++;
  }
  partitionsFile.close();
  cout << "done!\n";
}

void readWeightsFile(vector<double> &weights, ifstream &weightsFile) {
  cout << "Reading weight file " << flush;

  string line;
  string buf;

  double totWeight = 0.0;
  int nodeNr = 0;
  while (getline(weightsFile, line)) {
    istringstream read(line);
    read >> buf;
    weights[nodeNr] = atof(buf.c_str());
    totWeight += weights[nodeNr];
    nodeNr++;
  }

  cout << "with total weight " << totWeight << ", normalizing to 1..." << flush;
  int Nnodes = weights.size();
  if (Nnodes != nodeNr) {
    cout << nodeNr << "nodes in weights file, expecting " << Nnodes
         << " from partitions file, exiting...\n";
    abort();
  }

  // Normalize to 1
  for (int i = 0; i < nodeNr; i++)
    weights[i] /= totWeight;

  cout << "done!\n";
}

void printSignificanceClustering(vector<pair<bool, double>> &significantVec,
                                 vector<pair<int, int>> &mergers,
                                 multimap<double, TreeNode, greater<double>> &treeMap,
                                 const string &nodeOutFileName,
                                 const string &moduleOutFileName) {
  int M = treeMap.size();
  auto moduleId = vector<int>(M);
  int i = 0;
  for (auto it = treeMap.begin(); it != treeMap.end(); it++, i++)
    moduleId[i] = it->second.moduleId;

  ofstream outfile(nodeOutFileName);
  int Nnodes = significantVec.size();

  outfile << "# First column: 1 in significant core, 0 otherwise. Second "
             "column: the relative co-clustering with significant core nodes "
             "over all bootstrap clusterings\n";
  outfile << "Significant SignificantScore\n";

  for (int i = 0; i < Nnodes; i++)
    outfile << significantVec[i].first << " " << significantVec[i].second << '\n';

  outfile.close();

  outfile.open(moduleOutFileName);
  outfile << "# moduleId1 moduleId2 means that the significant core of moduleId1 "
             "cooccurs with moduleId2 more than a fraction 1 - conf of the "
             "samples.\n";
  outfile << "moduleId1 moduleId2\n";

  for (auto &merger : mergers)
    outfile << moduleId[merger.first] << " " << moduleId[merger.second] << '\n';

  outfile.close();
}

void generateTreeMap(vector<int> &rawPartition,
                     vector<double> &weights,
                     multimap<double, TreeNode, greater<double>> &treeMap) {
  int Nnodes = rawPartition.size();
  map<int, double> moduleWeights;
  map<int, vector<int>> moduleMembers;

  for (int i = 0; i < Nnodes; i++) {
    moduleWeights[rawPartition[i]] += weights[i];
    moduleMembers[rawPartition[i]].push_back(i);
  }

  auto modW_it = moduleWeights.begin();
  auto modM_it = moduleMembers.begin();

  for (; modM_it != moduleMembers.end(); modW_it++, modM_it++) {
    int Nmembers = modM_it->second.size();
    TreeNode tmp_tN;
    auto it_tM = treeMap.insert(make_pair(modW_it->second, tmp_tN));
    it_tM->second.moduleId = modM_it->first;

    for (int j = 0; j < Nmembers; j++) {
      it_tM->second.members.insert(make_pair(
          weights[modM_it->second[j]],
          make_pair(modM_it->second[j], to_string(modM_it->second[j]))));
    }
  }
}

void findConfCore(multimap<double, TreeNode, greater<double>> &treeMap,
                  vector<vector<int>> &bootPartitions,
                  vector<pair<bool, double>> &significantVec,
                  double conf,
                  std::mt19937 &mtrand) {
  double epsilon = 1.0e-5;

  int Nboots = bootPartitions.size();
  int Nremove = static_cast<int>((1.0 - conf) * Nboots + 0.5);

  int Nnode = bootPartitions[0].size();
  auto size = vector<double>(Nnode);

  int M = treeMap.size();
  auto moduleSize = vector<double>(M);
  auto moduleId = vector<int>(M);

  auto modSortMembers = vector<vector<int>>(M);
  int i = 0;

  for (auto &it : treeMap) {
    moduleId[i] = it.second.moduleId;
    int Nmem = it.second.members.size();
    modSortMembers[i] = vector<int>(Nmem);
    int j = 0;
    for (auto &member : it.second.members) {
      int mem = member.second.first;
      modSortMembers[i][j] = mem;
      size[mem] = (1.0 - epsilon) * member.first + epsilon / Nnode;
      moduleSize[i] += size[mem];
      j++;
    }
    i++;
  }

  cout << "\nMCMC to maximize significant core of modules, enumerated by size with module id in parenthesis:\n";
  for (int i = 0; i < M; i++) {
    cout << "module " << i + 1 << " (" << moduleId[i] << "): " << flush;

    int N = modSortMembers[i].size();
    std::uniform_int_distribution<int> intNdist(0, N - 1);
    auto confState = vector<bool>(N);
    auto maxConfState = vector<bool>(N);
    double maxScore = -1.0;

    double confSize = 0.0;
    int confN = 0;
    double maxConfSize = 0.0;
    int maxConfN = 0;
    double score = 0.0;
    int penalty = 0;

    double pW = 10.0 * moduleSize[i];

    if (N != 1) {
      int maxModNr = 0;

      for (int j = 0; j < Nboots; j++) {
        for (int k = 0; k < N; k++) {
          if (bootPartitions[j][modSortMembers[i][k]] > maxModNr) {
            maxModNr = bootPartitions[j][modSortMembers[i][k]];
          }
        }
      }

      maxModNr++;

      // Initiate weights of module assignments

      // Keep track of the size order of modules
      auto sortModSizes = vector<multimap<double, int, greater<double>>>(Nboots);

      // Keep track of which modules that are included
      auto mapModSizes = vector<vector<pair<int, multimap<double, int, greater<double>>::iterator>>>(Nboots);
      vector<vector<pair<int, double>>> bestMapModSizes(Nboots);

      for (int j = 0; j < Nboots; j++) {
        mapModSizes[j] = vector<pair<int, multimap<double, int, greater<double>>::iterator>>(maxModNr);
        bestMapModSizes[j] = vector<pair<int, double>>(maxModNr);
        for (int k = 0; k < maxModNr; k++) {
          mapModSizes[j][k].first = 0;
        }
      }

      // Randomized start
      for (int j = 0; j < N; j++) {
        if (bitdist(mtrand)) {
          confState[j] = true;
          confSize += size[modSortMembers[i][j]];
          confN++;

          for (int k = 0; k < Nboots; k++) {
            int modNr = bootPartitions[k][modSortMembers[i][j]];

            if (mapModSizes[k][modNr].first == 0) {
              double newSize = size[modSortMembers[i][j]];
              auto it = sortModSizes[k].insert(make_pair(newSize, modNr));
              mapModSizes[k][modNr].second = it;
              mapModSizes[k][modNr].first++;
            } else {
              double newSize = size[modSortMembers[i][j]] + mapModSizes[k][modNr].second->first;
              auto it = sortModSizes[k].insert(mapModSizes[k][modNr].second,
                                               make_pair(newSize, modNr));
              sortModSizes[k].erase(mapModSizes[k][modNr].second);
              mapModSizes[k][modNr].second = it;
              mapModSizes[k][modNr].first++;
            }
          }
        } else {
          confState[j] = false;
        }
      }

      multimap<double, pair<int, int>> scoreRank;
      score = 0.0;

      // Calculate penalty
      for (int j = 0; j < Nboots; j++) {
        double tmpScore = 0.0;
        int tmpPenalty = 0;
        if (!sortModSizes[j].empty()) {
          int modNr = sortModSizes[j].begin()->second;
          tmpScore = sortModSizes[j].begin()->first;
          // penalty is the number of nodes not in biggest field
          tmpPenalty = confN - mapModSizes[j][modNr].first;
        }
        scoreRank.insert(make_pair(tmpScore - pW * tmpPenalty, make_pair(tmpPenalty, j)));
        score += tmpScore;
        penalty += tmpPenalty;
      }

      // Remove worst results
      auto it = scoreRank.begin();
      for (int j = 0; j < Nremove; j++) {
        int bootNr = it->second.second;
        double tmpScore = 0.0;
        int tmpPenalty = 0;
        if (!sortModSizes[bootNr].empty()) {
          tmpScore = sortModSizes[bootNr].begin()->first;
          tmpPenalty = it->second.first;
        }
        score -= tmpScore;
        penalty -= tmpPenalty;
        it++;
      }

      // Monte Carlo to maximize confident size
      int Niter = static_cast<int>(pow(1.0 * N, 1.0));
      if (Niter < 100)
        Niter = 100;
      int attempts;
      int switches;

      bool search = true;
      while (search) {
        double T = 1.0;

        do {
          attempts = 0;
          switches = 0;

          for (int j = 0; j < Niter; j++) {
            multimap<double, pair<int, int>>().swap(scoreRank);
            double newConfSize = confSize;
            int newConfN = confN;
            double newScore = 0.0;
            int newPenalty = 0;

            int flip = intNdist(mtrand);
            int nodeNr = modSortMembers[i][flip];

            if (confState[flip]) {
              // Remove one node from confident subset
              newConfSize -= size[nodeNr];
              newConfN--;

              for (int k = 0; k < Nboots; k++) {
                int modNr = bootPartitions[k][nodeNr];
                double tmpScore = sortModSizes[k].begin()->first;
                int tmpPenalty = newConfN - mapModSizes[k][sortModSizes[k].begin()->second].first;

                if (mapModSizes[k][modNr].second == sortModSizes[k].begin()) {
                  tmpScore -= size[nodeNr];
                  tmpPenalty = newConfN - (mapModSizes[k][modNr].first - 1);

                  if (sortModSizes[k].size() > 1) {
                    auto it = sortModSizes[k].begin();
                    it++;

                    // Check if second in ranking is larger
                    if (it->first > tmpScore) {
                      tmpScore = it->first;
                      tmpPenalty = newConfN - mapModSizes[k][it->second].first;
                    }
                  }
                }

                scoreRank.insert(make_pair(tmpScore - pW * tmpPenalty, make_pair(tmpPenalty, k)));
                newScore += tmpScore;
                newPenalty += tmpPenalty;
              }
            } else {
              // Add one node to confident subset
              newConfSize += size[nodeNr];
              newConfN++;

              for (int k = 0; k < Nboots; k++) {
                int modNr = bootPartitions[k][nodeNr];
                double tmpScore = 0.0;
                int tmpPenalty = 0;

                if (sortModSizes[k].empty()) {
                  // No nodes in confident subset
                  tmpScore = size[nodeNr];
                  tmpPenalty = 0;
                } else {
                  if (mapModSizes[k][modNr].first == 0) {
                    // First confident node in module
                    if (size[nodeNr] > sortModSizes[k].begin()->first) {
                      tmpScore = size[nodeNr];
                      tmpPenalty = newConfN - 1;
                    } else {
                      tmpScore = sortModSizes[k].begin()->first;
                      tmpPenalty = newConfN - mapModSizes[k][sortModSizes[k].begin()->second].first;
                    }
                  } else {
                    // Not first confident node in module
                    if (mapModSizes[k][modNr].second == sortModSizes[k].begin()) {
                      tmpScore = sortModSizes[k].begin()->first + size[nodeNr];
                      tmpPenalty = newConfN - (mapModSizes[k][sortModSizes[k].begin()->second].first + 1);
                    } else {
                      if (mapModSizes[k][modNr].second->first + size[nodeNr] > sortModSizes[k].begin()->first) {
                        tmpScore = mapModSizes[k][modNr].second->first + size[nodeNr];
                        tmpPenalty = newConfN - (mapModSizes[k][modNr].first + 1);
                      } else {
                        tmpScore = sortModSizes[k].begin()->first;
                        tmpPenalty = newConfN - mapModSizes[k][sortModSizes[k].begin()->second].first;
                      }
                    }
                  }
                }

                scoreRank.insert(make_pair(tmpScore - pW * tmpPenalty, make_pair(tmpPenalty, k)));
                newScore += tmpScore;
                newPenalty += tmpPenalty;
              }
            }

            // Remove worst results
            auto it = scoreRank.begin();
            for (int j = 0; j < Nremove; j++) {
              int bootNr = it->second.second;
              if (!sortModSizes[bootNr].empty()) {
                double tmpScore = sortModSizes[bootNr].begin()->first;
                int tmpPenalty = it->second.first;
                newScore -= tmpScore;
                newPenalty -= tmpPenalty;
              }

              it++;
            }

            if (exp(((newScore - pW * newPenalty) - (score - pW * penalty)) / T) > realdist(mtrand)) {
              // Update data structures
              if (confState[flip]) {
                // Remove one node from confident subset
                for (int k = 0; k < Nboots; k++) {
                  int modNr = bootPartitions[k][nodeNr];
                  mapModSizes[k][modNr].first--;

                  // Remove last confident node from module
                  if (mapModSizes[k][modNr].first == 0) {
                    sortModSizes[k].erase(mapModSizes[k][modNr].second);
                  } else {
                    double newSize = mapModSizes[k][modNr].second->first - size[nodeNr];
                    auto it = sortModSizes[k].insert(mapModSizes[k][modNr].second, make_pair(newSize, modNr));
                    sortModSizes[k].erase(mapModSizes[k][modNr].second);
                    mapModSizes[k][modNr].second = it;
                  }
                }
              } else {
                // Add one node from confident subset
                for (int k = 0; k < Nboots; k++) {
                  int modNr = bootPartitions[k][nodeNr];

                  if (mapModSizes[k][modNr].first == 0) {
                    double newSize = size[nodeNr];
                    auto it = sortModSizes[k].insert(make_pair(newSize, modNr));
                    mapModSizes[k][modNr].second = it;
                    mapModSizes[k][modNr].first++;
                  } else {
                    double newSize = size[nodeNr] + mapModSizes[k][modNr].second->first;
                    auto it = sortModSizes[k].insert(mapModSizes[k][modNr].second, make_pair(newSize, modNr));
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

            if (penalty == 0 && score > maxScore) {
              for (int k = 0; k < N; k++) {
                maxConfState[k] = confState[k];
              }

              for (int k = 0; k < Nboots; k++) {
                for (int l = 0; l < maxModNr; l++) {
                  bestMapModSizes[k][l].first = mapModSizes[k][l].first;

                  if (mapModSizes[k][l].first > 0)
                    bestMapModSizes[k][l].second = mapModSizes[k][l].second->first;
                  else
                    bestMapModSizes[k][l].second = 0.0;
                }
              }

              maxScore = score;
              maxConfSize = confSize;
              maxConfN = confN;
            }
          }

          T *= 0.99;

        } while (switches > 0);

        if (maxScore > 0.0)
          search = false;
      }

      for (int j = 0; j < N; j++) {
        double p = 0.0;

        for (int k = 0; k < Nboots; k++) {
          // Node j in cluster i belongs to
          // cluster modNr in resampled network k
          int modNr = bootPartitions[k][modSortMembers[i][j]];
          p += bestMapModSizes[k][modNr].second / maxConfSize;
        }

        p /= 1.0 * Nboots;
        significantVec[modSortMembers[i][j]].second = p;
      }
    } else {
      maxConfState[0] = true;
      maxConfSize = moduleSize[i];
      maxConfN = 1;
      significantVec[modSortMembers[i][0]].second = 1.0;
    }

    cout << maxConfN << "/" << N << " confident nodes and " << maxConfSize
         << "/" << moduleSize[i] << " (" << 100 * maxConfSize / moduleSize[i]
         << " percent) of the flow.\n";

    for (int j = 0; j < N; j++) {
      significantVec[modSortMembers[i][j]].first = maxConfState[j];
    }
  }
}

void findConfModules(multimap<double, TreeNode, greater<double>> &treeMap,
                     vector<vector<int>> &bootPartitions,
                     vector<pair<bool, double>> &significantVec,
                     vector<pair<int, int>> &mergers,
                     double conf) {
  int M = treeMap.size();
  int N = bootPartitions[0].size();
  int Nboots = bootPartitions.size();
  auto moduleId = vector<int>(M);

  // Calculate total size of confident journals in field
  auto significantNodes = vector<vector<int>>(N);
  auto confFieldSize = vector<double>(M, 0.0);
  auto coExist = vector<multimap<int, pair<int, vector<int>>, greater<int>>>(M);

  int clusterNr = 0;
  for (auto &it : treeMap) {
    moduleId[clusterNr] = it.second.moduleId;
    for (auto &member : it.second.members) {
      if (significantVec[member.second.first].first)
        significantNodes[clusterNr].push_back(member.second.first);
    }
    clusterNr++;
  }

  auto coexistCount = vector<vector<int>>(M);
  for (int i = 0; i < M; i++)
    coexistCount[i] = vector<int>(Nboots, 0);

  cout << "\nNow calculate number of times two modules are clustered together, "
          "enumerated by size with module id in parenthesis:\n";

  int i = 0;
  for (auto it1 = treeMap.begin(); it1 != treeMap.end(); it1++) { // i
    int j = i;
    for (auto it2 = it1; it2 != treeMap.end(); it2++) { // j
      if (it1 != it2) {
        int coEx = 0;
        for (int k = 0; k < Nboots; k++) {
          int modNr = bootPartitions[k][significantNodes[j][0]];
          bool joined = true;
          int iNnode = significantNodes[i].size();
          int jNnode = significantNodes[j].size();

          for (int l = 0; l < jNnode; l++) {
            if (bootPartitions[k][significantNodes[j][l]] != modNr) {
              joined = false;
              break;
            }
          }

          if (joined) {
            for (int l = 0; l < iNnode; l++) {
              if (bootPartitions[k][significantNodes[i][l]] != modNr) {
                joined = false;
                break;
              }
            }
          }

          if (joined) {
            coEx++;
            coexistCount[i][k]++;
            coexistCount[j][k]++;

            auto it = coExist[i].find(j);
            if (it != coExist[i].end()) {
              it->second.first++;
              it->second.second.push_back(k);
            } else {
              vector<int> tmp;
              tmp.push_back(k);
              coExist[i].insert(make_pair(j, make_pair(1, tmp)));
            }
            it = coExist[j].find(i);
            if (it != coExist[j].end()) {
              it->second.first++;
              it->second.second.push_back(k);
            } else {
              vector<int> tmp;
              tmp.push_back(k);
              // FIXME never used
              it = coExist[j].insert(make_pair(i, make_pair(1, tmp)));
            }
          }
        }
      }

      j++;
    }

    // Re-sort co-exist strcuture
    auto tmp = coExist[i];
    multimap<int, pair<int, vector<int>>, greater<int>>().swap(coExist[i]);
    for (auto &it : tmp) {
      coExist[i].insert(make_pair(it.second.first, make_pair(it.first, it.second.second)));
    }

    i++;
  }

  auto mergeVec = vector<int>(M);
  i = 0;

  for (auto it1 = treeMap.begin(); it1 != treeMap.end(); it1++) {
    int singleN = 0;
    for (int k = 0; k < Nboots; k++) {
      if (coexistCount[i][k] == 0)
        singleN++;
    }

    cout << "Module " << i + 1 << " (" << moduleId[i] << ") is standalone "
         << singleN << "/" << Nboots << " times";

    if (singleN == Nboots) {
      cout << ".\n";
    } else {
      cout << " and clustered together with: ";
      for (auto &it2 : coExist[i]) {
        cout << it2.second.first + 1 << " (" << moduleId[it2.second.first]
             << ") " << it2.first << " times, ";
      }
      cout << '\n';
    }
    i++;
  }

  // Find merges
  for (int i = M - 1; i >= 0; i--) {
    bool searchMerge = true;
    mergeVec[i] = i;

    while (searchMerge) {
      int singleN = 0;
      for (int k = 0; k < Nboots; k++) {
        if (coexistCount[i][k] == 0)
          singleN++;
      }

      if (1.0 * singleN / Nboots < conf) {
        int mergeWith = coExist[i].begin()->second.first;

        for (int j = 0; j < coExist[i].begin()->first; j++) {
          coexistCount[i][coExist[i].begin()->second.second[j]]--;
          coexistCount[mergeWith][coExist[i].begin()->second.second[j]]--;
        }

        coExist[i].erase(coExist[i].begin());

        for (auto it = coExist[mergeWith].begin(); it != coExist[mergeWith].end(); it++) {
          if (it->second.first == i) {
            coExist[mergeWith].erase(it);
            break;
          }
        }

        if (mergeWith < i) {
          mergeVec[i] = mergeWith;
          mergers.emplace_back(i, mergeWith);
          searchMerge = false;
        }

      } else {
        searchMerge = false;
      }
    }
  }

  if (!mergers.empty()) {
    cout << "\nModule associations of modules that are not significantly standalone:\n";

    for (auto &merger : mergers) {
      cout << merger.first + 1 << " (" << moduleId[merger.first] << ") -> "
           << merger.second + 1 << " (" << moduleId[merger.second] << ")\n";
    }
  } else {
    cout << "All modules are significantly standalone.\n";
  }
}
