#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <utility>
#include <algorithm>
#include <string>
#include <random>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <omp.h>

using namespace std;

static size_t signatureWidth; // Signature size (in bits)
static size_t signatureSize;  // Signature size (in uint64_t)
static size_t kmerLength;     // Kmer length
static float density;         // % of sequence set as bits
static bool fastaOutput;      // Output fasta or csv

vector<pair<string, string>> loadFasta(const char *path)
{
  vector<pair<string, string>> sequences;

  FILE *fp = fopen(path, "r");
  if (!fp) {
    fprintf(stderr, "Failed to load %s\n", path);
    exit(1);
  }
  for (;;) {
    char seqNameBuf[8192];
    if (fscanf(fp, " >%[^\n]\n", seqNameBuf) < 1) break;
    string sequenceBuf;

    for (;;) {
      int c = fgetc(fp);
      if (c == EOF || c == '>') {
        ungetc(c, fp);
        break;
      }
      if (isalpha(c)) {
        sequenceBuf.push_back(c);
      }
    }
    sequences.push_back(make_pair(string(seqNameBuf), sequenceBuf));
  }
  fclose(fp);

  return sequences;
}

void generateSignature(uint64_t *output, const pair<string, string> &fasta)
{
  // Generate a signature from the kmers contained within
  
  string fastaSequence = fasta.second;
  // If the sequence is shorter than the kmer length, pad it with Xs
  while (fastaSequence.size() < kmerLength) {
    fastaSequence.push_back('X');
  }
  
  ranlux24_base rng;
  uniform_int_distribution<int> dist(-64 * signatureSize, signatureSize * 64 - 1);
  vector<int> unflattenedSignature(signatureSize * 64);
  int setBits = density * signatureSize * 64;
  //fprintf(stderr, "%s\n", fastaSequence.c_str());
  
  for (size_t i = 0; i < fastaSequence.size() - kmerLength + 1; i++) {
    seed_seq rngSeed(begin(fastaSequence) + i, begin(fastaSequence) + i + kmerLength);
    rng.seed(rngSeed);
    string kmer(begin(fastaSequence) + i, begin(fastaSequence) + i + kmerLength);
    //fprintf(stderr, "- %s\n", kmer.c_str());
    
    for (int j = 0; j < setBits; j++) {
      int bitPos = dist(rng);
      if (bitPos >= 0) {
        unflattenedSignature[bitPos] += 1;
      } else {
        unflattenedSignature[bitPos + 64 * signatureSize] -= 1;
      }
    }
  }
  fill(output, output + signatureSize, 0);
  for (size_t i = 0; i < signatureSize * 64; i++) {
    if (unflattenedSignature[i] > 0) {
      output[i / 64] |= (uint64_t)1 << (i % 64);
    }
  }
}

vector<uint64_t> convertFastaToSignatures(const vector<pair<string, string>> &fasta)
{
  vector<uint64_t> output;
  // Allocate space for the strings
  
  output.resize(fasta.size() * signatureSize);
  
  #pragma omp parallel for schedule(dynamic)
  for (size_t i = 0; i < fasta.size(); i++) {
    generateSignature(&output[signatureSize * i], fasta[i]);
  }
  
  return output;
}

void outputClusters(const vector<size_t> &clusters)
{
  for (size_t sig = 0; sig < clusters.size(); sig++)
  {
    printf("%llu,%llu\n", static_cast<unsigned long long>(sig), static_cast<unsigned long long>(clusters[sig]));
  }
}

void outputFastaClusters(const vector<size_t> &clusters, const vector<pair<string, string>> &fasta)
{
  fprintf(stderr, "Writing out %zu records\n", clusters.size());
  for (size_t sig = 0; sig < clusters.size(); sig++)
  {
    printf(">%llu\n%s\n", static_cast<unsigned long long>(clusters[sig]), fasta[sig].second.c_str());
  }
}
/*
vector<size_t> clusterSignatures(const vector<uint64_t> &sigs)
{
  auto rng = ranlux24_base();
  
  auto dist = uniform_int_distribution<size_t>(0, clusterCount - 1);
  size_t sigCount = sigs.size() / signatureSize;
  vector<size_t> clusters(sigCount);
  
  for (size_t i = 0; i < sigCount; i++) {
    clusters[i] = dist(rng);
  }
  
  return clusters;
}
*/

// Parameters
size_t ktree_order = 10;
size_t ktree_capacity = 1000000;

void dbgPrintSignature(const uint64_t *sig)
{
  fprintf(stderr, "%p: ", sig);
  for (size_t i = 0; i < signatureSize * 64; i++) {
    if (sig[i / 64] & (1ull << (i % 64))) {
      fprintf(stderr, "1");
    } else {
      fprintf(stderr, "0");
    }
  }
  fprintf(stderr, "\n");
}

void dbgPrintMatrix(const uint64_t *matrix)
{
  size_t ktree_csig_height = (ktree_order + 63) / 64;
  for (size_t i = 0; i < signatureSize * 64; i++) {
    fprintf(stderr, "%03zu:", i);
    for (size_t j = 0; j < ktree_csig_height * 64; j++) {
      auto val = matrix[i * ktree_csig_height + (j / 64)];
      if (val & (1ull << (j % 64))) {
        fprintf(stderr, "1");
      } else {
        fprintf(stderr, "0");
      }
    }
    fprintf(stderr, "\n");
    if (i >= 5) {
      fprintf(stderr, "...............\n");
      break;
    }
  }
}

template<class RNG>
vector<uint64_t> createRandomSigs(RNG &&rng, const vector<uint64_t> &sigs)
{
  constexpr size_t clusterCount = 2;
  vector<uint64_t> clusterSigs(signatureSize * clusterCount);
  size_t signatureCount = sigs.size() / signatureSize;
  uniform_int_distribution<size_t> dist(0, signatureCount - 1);
  bool finished = false;
  
  unordered_set<string> uniqueSigs;
  for (size_t i = 0; i < signatureCount; i++) {
    size_t sig = dist(rng);
    string sigData(signatureSize * sizeof(uint64_t), ' ');
    memcpy(&sigData[0], &sigs[sig * signatureSize], signatureSize * sizeof(uint64_t));
    uniqueSigs.insert(sigData);
    if (uniqueSigs.size() >= clusterCount) {
      finished = true;
      break;
    }
  }
  
  size_t i = 0;
  for (const auto &sig : uniqueSigs) {
    memcpy(&clusterSigs[i * signatureSize], sig.data(), signatureSize * sizeof(uint64_t));
    i++;
  }
  
  if (!finished) {
    if (uniqueSigs.size() != 1) {
      fprintf(stderr, "This should not happen\n");
      exit(1);
    }
    for (size_t i = 0; i < signatureSize; i++) {
      clusterSigs.push_back(clusterSigs[i]);
    }
  }
  
  return clusterSigs;
}

vector<vector<size_t>> createClusterLists(const vector<size_t> &clusters)
{
  constexpr size_t clusterCount = 2;
  vector<vector<size_t>> clusterLists(clusterCount);
  for (size_t i = 0; i < clusters.size(); i++) {
    clusterLists[clusters[i]].push_back(i);
  }
  return clusterLists;
}

vector<uint64_t> createClusterSigs(const vector<vector<size_t>> &clusterLists, const vector<uint64_t> &sigs)
{
  constexpr size_t clusterCount = 2;
  vector<uint64_t> clusterSigs(signatureSize * clusterCount);
  //#pragma omp parallel
  {
    vector<int> unflattenedSignature(signatureWidth);
    //#pragma omp for
    for (size_t cluster = 0; cluster < clusterLists.size(); cluster++) {
      fill(begin(unflattenedSignature), end(unflattenedSignature), 0);
      
      for (size_t signature : clusterLists[cluster]) {
        const uint64_t *signatureData = &sigs[signatureSize * signature];
        for (size_t i = 0; i < signatureWidth; i++) {
          uint64_t signatureMask = (uint64_t)1 << (i % 64);
          if (signatureMask & signatureData[i / 64]) {
            unflattenedSignature[i] += 1;
          } else {
            unflattenedSignature[i] -= 1;
          }
        }
      }
      
      uint64_t *flattenedSignature = &clusterSigs[cluster * signatureSize];
      for (size_t i = 0; i < signatureWidth; i++) {
        if (unflattenedSignature[i] > 0) {
          flattenedSignature[i / 64] |= (uint64_t)1 << (i % 64);
        }
      }
    }
  }
  return clusterSigs;
}

void reclusterSignatures(vector<size_t> &clusters, const vector<uint64_t> &meanSigs, const vector<uint64_t> &sigs)
{
  set<size_t> allClusters;
  for (size_t sig = 0; sig < clusters.size(); sig++) {
    const uint64_t *sourceSignature = &sigs[sig * signatureSize];
    size_t minHdCluster = 0;
    size_t minHd = numeric_limits<size_t>::max();

    for (size_t cluster = 0; cluster < 2; cluster++) {
      const uint64_t *clusterSignature = &meanSigs[cluster * signatureSize];
      size_t hd = 0;
      for (size_t i = 0; i < signatureSize; i++) {
        hd += __builtin_popcountll(sourceSignature[i] ^ clusterSignature[i]);
      }
      if (hd < minHd) {
        minHd = hd;
        minHdCluster = cluster;
      }
    }
    clusters[sig] = minHdCluster;
    allClusters.insert(minHdCluster);
  }
  
  if (allClusters.size() == 1) {
    // We can't have everything in the same cluster.
    // If this did happen, just split them evenly
    for (size_t sig = 0; sig < clusters.size(); sig++) {
      clusters[sig] = sig % 2;
    }
  }
}

// There are two kinds of ktree nodes- branch nodes and leaf nodes
// Both contain a signature matrix, plus their own signature
// (the root node signature does not matter and can be blank)
// Branch nodes then contain 'order' links to other nodes
// Leaf nodes do not.
// However, as leaf nodes may become branch nodes, we allocate
// the space anyway.
// As the space to be used is determined at runtime, we use
// parallel arrays, not structs

struct KTree {
  size_t root = numeric_limits<size_t>::max(); // # of root node
  vector<size_t> childCounts; // n entries, number of children
  vector<int> isBranchNode; // n entries, is this a branch node
  vector<size_t> childLinks; // n * o entries, links to children
  vector<size_t> parentLinks; // n entries, links to parents
  vector<uint64_t> means; // n * signatureSize entries, node signatures
  vector<uint64_t> matrices; // n * (o / 64) * signatureSize * 64 entries
  vector<omp_lock_t> locks; // n locks
  size_t order;
  size_t capacity = 0; // Set during construction, currently can't change
  size_t matrixHeight;
  size_t matrixSize;
  
  void reserve(size_t capacity) {
    // For safety, only call this at startup currently
    if (this->capacity != 0) {
      fprintf(stderr, "Reserve can only be called from 0 capacity\n");
      exit(1);
    }
    this->capacity = capacity;
    matrixHeight = (order + 63) / 64;
    matrixSize = matrixHeight * signatureSize * 64;
    
    #pragma omp parallel
    {
      #pragma omp single
      {
        childCounts.resize(capacity);
      }
      #pragma omp single
      {
        isBranchNode.resize(capacity);
      }
      #pragma omp single
      {
        childLinks.resize(capacity * order);
      }
      #pragma omp single
      {
        parentLinks.resize(capacity);
      }
      #pragma omp single
      {
        locks.resize(capacity);
      }
      #pragma omp single
      {
        matrices.resize(capacity * matrixSize);
      }
      #pragma omp single
      {
        means.resize(capacity * signatureSize);
      }
    }
  }
  
  KTree(size_t order_, size_t capacity) : order{order_} {    
    reserve(capacity);
  }
  
  size_t calcHD(const uint64_t *a, const uint64_t *b) const
  {
    size_t c = 0;
    for (size_t i = 0; i < signatureSize; i++) {
    c += __builtin_popcountll(a[i] ^ b[i]);
    }
    return c;
  }
  
  size_t traverse(const uint64_t *signature) const
  {
    size_t node = root;
    while (isBranchNode[node]) {
      size_t lowestHD = numeric_limits<size_t>::max();
      size_t lowestHDchild = 0;
      
      for (size_t i = 0; i < childCounts[node]; i++) {
        size_t child = childLinks[node * order + i];
        size_t hd = calcHD(&means[child * signatureSize], signature);
        if (hd < lowestHD) {
          lowestHD = hd;
          lowestHDchild = child;
        }
      }
      node = lowestHDchild;
    }
    return node;
  }
  
  void addSigToMatrix(uint64_t *matrix, size_t child, const uint64_t *sig) const
  {
    size_t childPos = child / 64;
    size_t childOff = child % 64;
    
    //fprintf(stderr, "Adding this signature:\n");
    //dbgPrintSignature(sig);
    //fprintf(stderr, "To this matrix:\n");
    //dbgPrintMatrix(matrix);
    
    for (size_t i = 0; i < signatureSize * 64; i++) {
      matrix[i * matrixHeight + childPos] |= ((sig[i / 64] >> (i % 64)) & 0x01) << childOff;
    }
    //fprintf(stderr, "Resulting in:\n");
    //dbgPrintMatrix(matrix);
  }
  void removeSigFromMatrix(uint64_t *matrix, size_t child) const
  {
    size_t childPos = child / 64;
    size_t childOff = child % 64;
    
    uint64_t mask = ~(1ull << childOff);
    
    //fprintf(stderr, "Removing the %zuth child from matrix\n", child);    
    for (size_t i = 0; i < signatureSize * 64; i++) {
      matrix[i * matrixHeight + childPos] &= mask;
    }
    //fprintf(stderr, "Resulting in:\n");
    //dbgPrintMatrix(matrix);
  }
  
  void recalculateSig(size_t node)
  {
    size_t children = childCounts[node];
    uint64_t *matrix = &matrices[node * matrixSize];
    uint64_t *sig = &means[node * signatureSize];
    fill(sig, sig + signatureSize, 0ull);
    
    auto threshold = (children / 2) + 1;
    
    for (size_t i = 0; i < signatureSize * 64; i++) {
      size_t c = 0;
      for (size_t j = 0; j < matrixHeight; j++) {
        auto val = matrix[i * matrixHeight + j];
        c += __builtin_popcountll(val);
      }
      if (c >= threshold) {
        sig[i / 64] |= 1ull << (i % 64);
      }
    }
    //fprintf(stderr, "Mean sig:\n");
    //dbgPrintSignature(sig);
  }
  
  void recalculateUp(size_t node)
  {
    size_t limit = 10;
    //fprintf(stderr, "RecalculateUp %zu\n", node);
    while (node != root) {
      recalculateSig(node);
      node = parentLinks[node];
      if (omp_test_lock(&locks[node])) {
        omp_unset_lock(&locks[node]);
      } else {
        break;
      }
      
      // Put a limit on how far we go up
      // At some point it stops mattering, plus this helps avoid inf loops
      // caused by cycles getting into the tree structure
      limit--;
      if (limit == 0) return;
      //fprintf(stderr, "-> %zu\n", node);
    }
  }
  
  size_t getNewNodeIdx(vector<size_t> &insertionList)
  {
    if (insertionList.empty()) {
      fprintf(stderr, "ERROR: ran out of insertion points\n");
      exit(1);
    }
    size_t idx = insertionList.back();
    insertionList.pop_back();
    
    // Initialise lock
    omp_init_lock(&locks[idx]);
    return idx;
  }
  
  template<class RNG>
  void splitNode(RNG &&rng, size_t node, const uint64_t *sig, vector<size_t> &insertionList, size_t link)
  {
    //fprintf(stderr, "Splitting node %zu\n", node);
    // Add 'sig' to the current node, splitting it in the process
    //fprintf(stderr, "Adding signature:\n");
    //dbgPrintSignature(sig);
    size_t nodeSigs = childCounts[node] + 1;
    vector<uint64_t> sigs(nodeSigs * signatureSize);
    memcpy(&sigs[ childCounts[node] * signatureSize], sig, signatureSize * sizeof(uint64_t));
    
    for (int i = 0; i < childCounts[node]; i++) {
      uint64_t *currentSig = &sigs[i * signatureSize];
      uint64_t *matrix = &matrices[node * matrixSize];
      for (size_t j = 0; j < signatureSize * 64; j++) {
        currentSig[j / 64] |= ((matrix[j * matrixHeight + i / 64] >> (i % 64)) & 1) << (j % 64);
      }
    }
    
    /*
    fprintf(stderr, "Signatures converted for clustering:\n");
    for (size_t i = 0; i < nodeSigs; i++) {
      uint64_t *currentSig = &sigs[i * signatureSize];
      dbgPrintSignature(currentSig);
    }
    */
    
    vector<uint64_t> meanSigs = createRandomSigs(rng, sigs);
    vector<size_t> clusters(nodeSigs);
    vector<vector<size_t>> clusterLists;
    for (int iteration = 0; iteration < 4; iteration++) {
      //fprintf(stderr, "Iteration %d\n", iteration);
      reclusterSignatures(clusters, meanSigs, sigs);
      clusterLists = createClusterLists(clusters);
      meanSigs = createClusterSigs(clusterLists, sigs);
    }
    
    /*
    // Display clusters (debugging purposes)
    for (const auto &clusterList : clusterLists) {
      fprintf(stderr, "Cluster:\n");
      for (size_t seqIdx : clusterList) {
        uint64_t *currentSig = &sigs[seqIdx * signatureSize];
        dbgPrintSignature(currentSig);
      }
    }
    */
    
    // Create the sibling node
    size_t sibling = getNewNodeIdx(insertionList);
    
    size_t newlyAddedIdx = childCounts[node];
    
    childCounts[sibling] = clusterLists[1].size();
    isBranchNode[sibling] = isBranchNode[node];
    {
      size_t siblingIdx = 0;
      for (size_t seqIdx : clusterLists[1]) {
        if (seqIdx < newlyAddedIdx) {
          childLinks[sibling * order + siblingIdx] = childLinks[node * order + seqIdx];
        } else {
          childLinks[sibling * order + siblingIdx] = link;
        }
        // If this is a branch node, relink the child to the new parent
        if (isBranchNode[sibling]) {
          parentLinks[childLinks[sibling * order + siblingIdx]] = sibling;
        }
        addSigToMatrix(&matrices[sibling * matrixSize], siblingIdx, &sigs[seqIdx * signatureSize]);
        siblingIdx++;
      }
    }
    memcpy(&means[sibling * signatureSize], &meanSigs[1 * signatureSize], signatureSize * sizeof(uint64_t));
    
    // Fill the current node with the other cluster of signatures
    {
      fill(&matrices[node * matrixSize], &matrices[node * matrixSize] + matrixSize, 0ull);
      size_t nodeIdx = 0;
      for (size_t seqIdx : clusterLists[0]) {
        if (seqIdx < newlyAddedIdx) {
          childLinks[node * order + nodeIdx] = childLinks[node * order + seqIdx];
        } else {
          childLinks[node * order + nodeIdx] = link;
        }
        // If this is a branch node, relink the child to the new parent
        if (isBranchNode[node]) {
          parentLinks[childLinks[node * order + nodeIdx]] = node;
        }
        addSigToMatrix(&matrices[node * matrixSize], nodeIdx, &sigs[seqIdx * signatureSize]);
        nodeIdx++;
      }
    }
    childCounts[node] = clusterLists[0].size();
    
    // Is this the root level?
    if (node == root) {
      //fprintf(stderr, "Node being split is root node\n");
      
      // Create a new root node
      size_t newRoot;
      newRoot = getNewNodeIdx(insertionList);
      
      // Link this node and the sibling to it
      parentLinks[node] = newRoot;
      parentLinks[sibling] = newRoot;

      childCounts[newRoot] = 2;
      isBranchNode[newRoot] = 1;
      childLinks[newRoot * order + 0] = node;
      childLinks[newRoot * order + 1] = node;
      addSigToMatrix(&matrices[newRoot * matrixSize], 0, &meanSigs[0 * signatureSize]);
      addSigToMatrix(&matrices[newRoot * matrixSize], 1, &meanSigs[1 * signatureSize]);
      
      root = newRoot;
    } else {
      
      // First, update the reference to this node in the parent with the new mean
      size_t parent = parentLinks[node];
      
      // Lock the parent
      omp_set_lock(&locks[parent]);
      
      size_t idx = numeric_limits<size_t>::max();
      for (size_t i = 0; i < childCounts[parent]; i++) {
        if (childLinks[parent * order + i] == node) {
          idx = i;
          break;
        }
      }
      if (idx == numeric_limits<size_t>::max()) {
        //fprintf(stderr, "Error: node %zu is not its parent's (%zu) child\n", node, parent);
        
        // Abort. Unlock the parent and get out of here
        omp_unset_lock(&locks[parent]);
        return;
        
        //exit(1);
      }
      
      removeSigFromMatrix(&matrices[parent * matrixSize], idx);
      addSigToMatrix(&matrices[parent * matrixSize], idx, &meanSigs[0 * signatureSize]);
      
      // Connect sibling node to parent
      parentLinks[sibling] = parent;
      
      // Now add a link in the parent node to the sibling node
      if (childCounts[parent] + 1 < order) {
        addSigToMatrix(&matrices[parent * matrixSize], childCounts[parent], &meanSigs[1 * signatureSize]);
        childLinks[parent * order + childCounts[parent]] = sibling;
        childCounts[parent]++;
        
        // Update signatures (may change?)
        recalculateUp(parent);
      } else {
        splitNode(rng, parent, &meanSigs[1 * signatureSize], insertionList, sibling);
      }
      // Unlock the parent
      omp_unset_lock(&locks[parent]);
    }
    
    //fprintf(stderr, "Split finished\n");
  }
  
  template<class RNG>
  void insert(RNG &&rng, const uint64_t *signature, vector<size_t> &insertionList)
  {
    // Warning: ALWAYS INSERT THE FIRST NODE SINGLE-THREADED
    // We don't have any protection from this because it would slow everything down to do so
    if (root == numeric_limits<size_t>::max()) {
      root = getNewNodeIdx(insertionList);
      childCounts[root] = 0;
      isBranchNode[root] = 0;
    }
    
    size_t insertionPoint = traverse(signature);
    
    //fprintf(stderr, "Inserting at %zu\n", insertionPoint);
    omp_set_lock(&locks[insertionPoint]);
    if (childCounts[insertionPoint] < order) {
      addSigToMatrix(&matrices[insertionPoint * matrixSize], childCounts[insertionPoint], signature);
      childCounts[insertionPoint]++;
    } else {
      splitNode(rng, insertionPoint, signature, insertionList, 0);
    }
    omp_unset_lock(&locks[insertionPoint]);
    
    //fprintf(stderr, "Node %zu now has %zu leaves\n", insertionPoint, childCounts[insertionPoint]);
  }
  
  void destroyLocks(size_t node)
  {
    omp_destroy_lock(&locks[node]);
    if (isBranchNode[node]) {
      for (size_t i = 0; i < childCounts[node]; i++) {
        destroyLocks(childLinks[node * order + i]);
      }
    }
  }
  void destroyLocks()
  {
    destroyLocks(root);
  }
};

void compressClusterList(vector<size_t> &clusters)
{
  unordered_map<size_t, size_t> remap;
  for (size_t &clus : clusters) {
    if (remap.count(clus)) {
      clus = remap[clus];
    } else {
      size_t newClus = remap.size();
      remap[clus] = newClus;
      clus = newClus;
    }
  }
  fprintf(stderr, "Output %zu clusters\n", remap.size());
}

vector<size_t> clusterSignatures(const vector<uint64_t> &sigs)
{
  size_t sigCount = sigs.size() / signatureSize;
  vector<size_t> clusters(sigCount);
  KTree tree(ktree_order, ktree_capacity);
  
  size_t firstNodes = 1;
  if (firstNodes > sigCount) firstNodes = sigCount;
  
  vector<size_t> insertionList;
  for (size_t i = 0; i < firstNodes; i++) {
    insertionList.push_back(firstNodes - i - 1);
  }
  
  default_random_engine rng;
  // Insert first 1 nodes single-threaded
  for (size_t i = 0; i < firstNodes; i++) {
    tree.insert(rng, &sigs[i * signatureSize], insertionList);
  }
  
  // What's the next free insertion point?
  size_t nextFree = insertionList.back();
  
  #pragma omp parallel
  {
    default_random_engine rng;
    vector<size_t> insertionList;
    
    #pragma omp for
    for (size_t i = nextFree; i < ktree_capacity; i++) {
      insertionList.push_back(ktree_capacity - i - 1);
    }
    
    #pragma omp for
    for (size_t i = firstNodes; i < sigCount; i++) {
      tree.insert(rng, &sigs[i * signatureSize], insertionList);
    }
  }
  
  // We've created the tree. Now reinsert everything
  #pragma omp parallel for
  for (size_t i = 0; i < sigCount; i++) {
    size_t clus = tree.traverse(&sigs[i * signatureSize]);
    clusters[i] = clus;
  }
  
  // We want to compress the cluster list down
  compressClusterList(clusters);
  
  // Recursively destroy all locks
  tree.destroyLocks();
  
  return clusters;
}

int main(int argc, char **argv)
{
  if (argc < 2) {
    fprintf(stderr, "Usage: %s (options) [fasta input]\n", argv[0]);
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -sw [signature width]\n");
    fprintf(stderr, "  -k [kmer length]\n");
    fprintf(stderr, "  -d [signature density]\n");
    fprintf(stderr, "  -o [tree order]\n");
    fprintf(stderr, "  -c [starting capacity]\n");
    fprintf(stderr, "  --fasta-output\n");
    return 1;
  }
  signatureWidth = 256;
  kmerLength = 5;
  density = 1.0f / 21.0f;
  fastaOutput = false;
  
  string fastaFile = "";
  
  for (int a = 1; a < argc; a++) {
    string arg(argv[a]);
    if (arg == "-sw") signatureWidth = atoi(argv[++a]);
    else if (arg == "-k") kmerLength = atoi(argv[++a]);
    else if (arg == "-d") density = atof(argv[++a]);
    else if (arg == "-o") ktree_order = atoi(argv[++a]);
    else if (arg == "-c") ktree_capacity = atoi(argv[++a]);
    else if (arg == "--fasta-output") fastaOutput = true;
    else if (fastaFile.empty()) fastaFile = arg;
    else {
      fprintf(stderr, "Invalid or extra argument: %s\n", arg.c_str());
      exit(1);
    }
  }
    
  if (signatureWidth <= 0 || signatureWidth % 64 != 0) {
    fprintf(stderr, "Error: signature width is not a multiple of 64\n");
    return 1;
  }
  if (kmerLength <= 0) {
    fprintf(stderr, "Error: kmer length must be a positive nonzero integer\n");
    return 1;
  }
  if (density < 0.0f || density > 1.0f) {
    fprintf(stderr, "Error: density must be a positive value between 0 and 1\n");
    return 1;
  }

  signatureSize = signatureWidth / 64;
  
  fprintf(stderr, "Loading fasta...");
  auto fasta = loadFasta(fastaFile.c_str());
  fprintf(stderr, " loaded %llu sequences\n", static_cast<unsigned long long>(fasta.size()));
  fprintf(stderr, "Converting fasta to signatures...");
  auto sigs = convertFastaToSignatures(fasta);
  fprintf(stderr, " done\n");
  fprintf(stderr, "Clustering signatures...\n");
  auto clusters = clusterSignatures(sigs);
  fprintf(stderr, "Writing output\n");
  if (!fastaOutput) {
    outputClusters(clusters);
  } else {
    outputFastaClusters(clusters, fasta);
  }
  
  return 0;
}
