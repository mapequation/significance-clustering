#include "significanceclustering.h"

using namespace std;

const string CALL_SYNTAX =
    "Call: ./sigclu [-s <seed>] [-c <confidencelevel>] [-w <weightsfile>] "
    "partitionsfile nodeoutfile moduleoutfile\n"
    "seed: Any positive integer.\n"
    "confidencelevel: The confidence as a fraction, default is 0.95.\n"
    "partitionsfile: Each column represents a partition, the first for the raw "
    "partition and the remaining for bootstrap partitions. Row number "
    "corresponds to node id.\n"
    "nodeoutfile: 1 or 0 if a node does or does not belong to the significant "
    "core of its module.\n"
    "moduleoutfile: moduleId1 moduleId2 means that the significant core of "
    "moduleId1 cooccurs with moduleId2 more than a fraction 1 - conf of the "
    "samples.\n"
    "weightsfile: One column for weights of each node.  Row number corresponds "
    "to node id.\n";

int main(int argc, char *argv[]) {
  if (argc == 1) {
    cout << CALL_SYNTAX;
    exit(-1);
  }

  unsigned int seed = 1234;
  double conf = 0.95;
  string partitionsFileName = "noname";
  string weightsFileName = "noname";
  string nodeOutFileName = "noname";
  string moduleOutFileName = "noname";

  int argNr = 1;
  while (argNr < argc) {
    if (to_string(argv[argNr]) == "-h") {
      cout << CALL_SYNTAX;
      exit(-1);
    } else if (to_string(argv[argNr]) == "-s") {
      argNr++;
      seed = atoi(argv[argNr]);
      argNr++;
    } else if (to_string(argv[argNr]) == "-c") {
      argNr++;
      conf = atof(argv[argNr]);
      argNr++;
    } else if (to_string(argv[argNr]) == "-w") {
      argNr++;
      weightsFileName = to_string(argv[argNr]);
      argNr++;
    } else {

      if (argv[argNr][0] == '-') {
        cout << "Unknown command: " << to_string(argv[argNr]) << '\n'
             << CALL_SYNTAX;
        exit(-1);
      } else if (partitionsFileName == "noname") {
        partitionsFileName = to_string(argv[argNr]);
        argNr++;
      } else if (nodeOutFileName == "noname") {
        nodeOutFileName = to_string(argv[argNr]);
        argNr++;
      } else if (moduleOutFileName == "noname") {
        moduleOutFileName = to_string(argv[argNr]);
        argNr++;
      } else {
        cout << "Too many parameters without flags.\n"
             << CALL_SYNTAX;
        exit(-1);
      }
    }
  }

  if (partitionsFileName == "noname") {
    cout << "No partitions file provided.\n"
         << CALL_SYNTAX;
    exit(-1);
  } else if (nodeOutFileName == "noname") {
    cout << "No node out file provided.\n"
         << CALL_SYNTAX;
    exit(-1);
  } else if (moduleOutFileName == "noname") {
    cout << "No module out file provided.\n"
         << CALL_SYNTAX;
    exit(-1);
  }

  cout << "Setup:\n";
  cout << "-->Using seed: " << seed << '\n';
  cout << "-->Confidence level: " << conf << '\n';
  cout << "-->Reading from partitions file: " << partitionsFileName << '\n';

  if (weightsFileName == "noname")
    cout << "-->No weights file provided. Using uniform weights of nodes in signifificance clustering.\n";
  else
    cout << "-->Using weights of nodes in significance clustering from file: " << weightsFileName << '\n';

  cout << "-->Writing results to files: " << nodeOutFileName << " and " << moduleOutFileName << '\n';

  vector<int> rawPartition;
  vector<vector<int>> bootPartitions;

  int Nnodes = 0;
  int NbootSamples = 0;

  std::mt19937 mtrand(seed);
  ifstream partitionsFile;

  try {
    partitionsFile.open(partitionsFileName);
  } catch (const std::exception &e) {
    cout << "Could not open " << partitionsFileName << ", no such file.\n";
    exit(-1);
  }

  // Read partitions file
  readPartitionsFile(rawPartition, bootPartitions, partitionsFile, Nnodes, NbootSamples);

  vector<double> weights(Nnodes, 1.0 / Nnodes);
  if (weightsFileName != "noname") {
    // Read wights from file if provided
    ifstream weightsFile(weightsFileName);
    readWeightsFile(weights, weightsFile);
  }

  // Store modules by weight and nodes in modules by weight
  multimap<double, TreeNode, greater<double>> treeMap;
  generateTreeMap(rawPartition, weights, treeMap);

  // Calculate significance clusters
  vector<pair<bool, double>> significantVec(Nnodes);
  // find core in each module
  findConfCore(treeMap, bootPartitions, significantVec, conf, mtrand);
  vector<pair<int, int>> mergers;
  // for cores to be significant,
  // modules should significantly stand alone
  findConfModules(treeMap, bootPartitions, significantVec, mergers, conf);

  printSignificanceClustering(significantVec, mergers, treeMap, nodeOutFileName, moduleOutFileName);
}
