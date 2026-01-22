#include "FusedLasso_graph.h"

void  Graph::distance(const bool fromSource) 
{
  list<int>::iterator listIt;
  // set the distances to maximum
  for(listIt = maxFlowNodeList.begin(); listIt != maxFlowNodeList.end(); ++listIt) {
    dist[*listIt] = maxFlowNodeList.size() + 2;
  }
  queue<int> next;
  int u; // number of node that is currently being looked at
  GraphEdge* e;
  Node::iterator edgeIt;
  
  // initialize the nodes that are connected to the source
  for(listIt = maxFlowNodeList.begin(); listIt != maxFlowNodeList.end(); ++listIt) {
    if(fromSource) {
      if(sourceCapacity[*listIt] - tolerance > sourceFlow[*listIt]) {
        dist[*listIt] = 1;
        next.push(*listIt);
      }
    }
    else {
      if(sinkCapacity[*listIt] - tolerance > sinkFlow[*listIt]) {
        dist[*listIt] = 1;
        next.push(*listIt);
      }
    }
  }   
  
  while(!next.empty())
  {
    u=next.front();
    next.pop();
    // checking node u
    // go through all edges in node u;
    for(edgeIt= nodes[u].begin(); edgeIt != nodes[u].end(); ++edgeIt)
    {
      if(!subNodes[(*edgeIt)->to]) {
        continue; // nothing to do
      }
      if(fromSource) { // if fromSource, then do it forward, otherwise backwards
        e=*edgeIt;
      }
      else {
        e = (*edgeIt)->backwards;
      }
      if(e->flow < e->capacity-tolerance) // is reachable
      {
        if(dist[(*edgeIt) -> to]>dist[u]+1) // distance can be updated
        {
          dist[(*edgeIt)->to] = dist[u]+1;
          next.push((*edgeIt)->to);
        }
      }
    }
  }
};

bool Graph::push(const int from, GraphEdge *e)
{
  bool isToActive;
  // push as much flow as possible (also check the backwards flow first)
  double addFlowBack = min(exFlow[from], e->backwards->flow);
  double addFlow = min(exFlow[from]-addFlowBack, e->capacity - e->flow);
  e->flow +=addFlow;
  e->backwards->flow -= addFlowBack;
  
  exFlow[from] -= addFlow + addFlowBack;
  // save if to is not already an active node
  isToActive = (exFlow[e->to] > tolerance);
  exFlow[e->to] += addFlow + addFlowBack;
  if(!isToActive) // add it to the active nodes if not active already and not source or sink
  {
    insertActiveNode(e->to);
  }
  return(exFlow[from] > tolerance); // return true if the node is still active
}


bool Graph::pushToSource(const int from) {
  double addFlowBack = min(exFlow[from], sourceFlow[from]);
  exFlow[from] -= addFlowBack;
  sourceFlow[from] -= addFlowBack;
  return(exFlow[from] > tolerance);
}

bool Graph::pushToSink(const int from) {
  double addFlow = min(exFlow[from], sinkCapacity[from] - sinkFlow[from]);
  exFlow[from] -= addFlow;
  sinkFlow[from] += addFlow;
  return(exFlow[from] > tolerance);
}


int Graph::findDist(int nodeNum)
{
  Node::iterator edgeIt;
  
  // check if nodeNum is connected to the source
  int newDist;
  if(sourceFlow[nodeNum] > tolerance) {
    newDist = maxFlowNodeList.size() + 2;
  }
  else {
    newDist = 2 * maxFlowNodeList.size() + 2;
  }
  
  for(edgeIt = nodes[nodeNum].begin(); edgeIt != nodes[nodeNum].end(); ++edgeIt)
  {
    if(!subNodes[(*edgeIt)->to]) {
      continue;
    }
    // check that there is residual capacity left
    if((*edgeIt)->backwards->flow + (*edgeIt)->capacity - (*edgeIt)->flow > tolerance)
    {
      newDist = min(newDist, dist[(*edgeIt)->to]+1);
    }
  }
  return(newDist);
}


void Graph::preprocess()
{
  // clean the active nodes;
  activeByDist.assign(2* maxFlowNodeList.size() + 3, list<int>(0));
  level=-1;
  
  // compute the distance label from the sink using the current flow
  distance();
  
  list<int>::iterator listIt;
  // go through all nodes connected to the source and push the maximum flow
  for(listIt = maxFlowNodeList.begin(); listIt != maxFlowNodeList.end(); ++listIt)
  {
    // what is the excess in the connected node
    exFlow[*listIt] = sourceCapacity[*listIt] - sourceFlow[*listIt];
    // set flow to maximum
    sourceFlow[*listIt] = sourceCapacity[*listIt];
    // if positive excess flow, insert it as an active node
    if(exFlow[*listIt] > tolerance)
    {
      insertActiveNode(*listIt);
    }
  }
}

bool Graph::pushRelabel(const int i)
{
  // go through all edges from node i and look for admissible ones
  // that is d(i) > d(j)+1 and flow(i,j) < cap(i,j)
  bool didPush = false; // saves if a successfull push has occured
  Node::iterator edgeIt;
  // check if we can push to the sink
  if(sinkCapacity[i] - sinkFlow[i] > tolerance) {
    if(!pushToSink(i)) {
      return(false);               
    }
    didPush = true;
  }
  
  // check if we can push back to the source
  if(dist[i] > maxFlowNodeList.size() + 1 && sourceFlow[i] > 0) {
    if(!pushToSource(i)) {
      return(false);
    }
    didPush = true;
  }
  
  for(edgeIt = nodes[i].begin(); edgeIt != nodes[i].end(); ++edgeIt)
  {
    if(!subNodes[(*edgeIt)->to]) {
      continue; // nothing to do
    }
    
    // check if the edge is admissible
    if((dist[i] > dist[(*edgeIt)->to]) && ((*edgeIt)->capacity > (*edgeIt)->flow + tolerance))
    {
      didPush=true;
      if(!push(i, *edgeIt)) // node became inactive
      {
        return(false);
      }
    }
  }
  if(!didPush) // no successfull push was possible; relabel
  {
    dist[i] = findDist(i);
  }
  
  return(true); // node still active as there is still an excess
}


bool Graph::getLargestActiveNode(int& nodeNum)
{
  // check that there is an active node with distance level
  if(level<0) // nothing in here
  {
    return(false);
  }
  if(activeByDist[level].empty()) // no there isn't, step down until there is
  {
    do
    {
      --level;
    }
    while((level>=0) && (activeByDist[level].empty()) ); // check first if level still >=0, then if list is empty
    if(level<0) // no active nodes left
    {
      return(false);
    }
  }
  // save the largest element and delete it;
  nodeNum = activeByDist[level].front();
  activeByDist[level].pop_front();
  return(true); // saved active node
}

void Graph::insertActiveNode(const int nodeNum)
{
  // check if new node is the largest
  if(dist[nodeNum]>level)
  {
    level = dist[nodeNum];
  }
  // insert the node
  activeByDist[dist[nodeNum]].push_back(nodeNum);
}


bool Graph::findMaxFlowPrelabelPush()
{
  int activeNodeNum;
  preprocess();
  while(getLargestActiveNode(activeNodeNum))
  {
    //        cout << level << " ";
    if(pushRelabel(activeNodeNum)) // node remains active; if inactive, does not need to be replaced
    {
      insertActiveNode(activeNodeNum);
    }
  }
  
  return(true);
}



bool Graph::findMaxFlowEdmondsKarp() {
  vector<int> parentTable(nodes.size(), -1);
  vector<int> parentEdgeIndex(nodes.size(), -1);
  vector<int> nodeWithParent;
  vector<double> maxFlow(nodes.size(), 0);
  int endNode;
  
  list<int>::iterator listIt;
  double lastFlow;
  
  for(listIt = maxFlowNodeList.begin(); listIt != maxFlowNodeList.end(); ++listIt) {
    lastFlow = 1;
    while(lastFlow > 0) {
      breadthFirstSearch(*listIt, parentTable, parentEdgeIndex, nodeWithParent, endNode, maxFlow);
      lastFlow = addFlow(*listIt, endNode, parentTable, parentEdgeIndex, maxFlow); 
      //printGraph(cout);
    }
  }
  return(true);
}






void Graph::addEdge(const int from, const int to, const double capacityFromTo, const double capacityToFrom, const double flowFromTo, const double flowToFrom)
{
  GraphEdge* e1 = new(GraphEdge);
  GraphEdge* e2 = new(GraphEdge);
  
  // set where they are going to
  e1->to = to;
  e2->to = from;
  
  // set the backward edge
  e1->backwards = e2;
  e2->backwards = e1;
  
  // set the flow (also works as derivative)
  e1->flow = flowFromTo;
  e2->flow = flowToFrom;
  
  // set the capacity
  e1->capacity = capacityFromTo;
  e2->capacity = capacityToFrom;
  
  nodes[from].push_back(e1);
  nodes[to].push_back(e2);
};



list<int> Graph::connectedTo(const list<int>& nodeList)
{
  list<int> conn; // nodes it is connected to;
  list<int>::const_iterator listIt;
  Node::iterator edgeIt;
  
  // first, set all nodes that are currently of interest
  for(listIt = nodeList.begin(); listIt != nodeList.end(); ++listIt) {
    subNodesCT[*listIt] = true;
  }
  
  // step through all nodes in the set
  for(listIt = nodeList.begin(); listIt != nodeList.end(); ++listIt)
  {
    if(*listIt < nodes.size()) {
      for(edgeIt = nodes[*listIt].begin(); edgeIt != nodes[*listIt].end(); ++edgeIt) {
        if(!subNodesCT[(*edgeIt)->to]) {
          conn.push_back((*edgeIt)->to);
        }
      }
    }
  }
  
  // clean up all nodes that are currently of interest
  clearSubNodesCT(nodeList);
  
  conn.sort();
  conn.unique();
  return(conn);
}


void Graph::breadthFirstSearch(const int startNode, vector<int>&parentTable, vector<int>& parentEdgeIndex, vector<int>& nodeWithParent, int& endNode, vector<double>& maxFlow) {
  vector<int> nodesToSearch;
  
  // clear the parent table
  for(unsigned int i = 0; i < nodeWithParent.size(); ++i) {
    parentTable[nodeWithParent[i]] = -1;
  }
  nodeWithParent.clear();
  
  list<int> nodeToSearch;
  nodesToSearch.push_back(startNode);
  maxFlow[startNode] = sourceCapacity[startNode] - sourceFlow[startNode]; 
  
  if(maxFlow[startNode] == 0) {
    endNode = startNode;
    return;
  }
  
  int u;
  GraphEdge* e;
  while(nodesToSearch.size() > 0) {
    u = nodesToSearch.back();
    if(sinkCapacity[u] - sinkFlow[u] > 0) {
      endNode = u;
      return;
    }
    
    nodesToSearch.pop_back();
    for(int edgeIndex = 0; edgeIndex <  nodes[u].size(); ++edgeIndex) {
      e = nodes[u][edgeIndex];
      if(!subNodes[e->to]) {
        continue; // not part of the maxFlowGraph
      }
      
      // check if u is connected to sink with spare capacity
      // go through the edges for things to search
      double addFlow = e->capacity - e->flow + e->backwards->flow;
      if(addFlow > 0 && parentTable[e->to] == -1) {
        parentTable[e->to] = u;
        parentEdgeIndex[e->to] = edgeIndex;
        nodeWithParent.push_back(e->to);
        maxFlow[e->to] = min(maxFlow[u], addFlow);
        nodesToSearch.push_back(e->to);
      }
    }
    
  }
  
  // did not find anything
  endNode = startNode;
  maxFlow[startNode] = 0;
  return;
}


double Graph::addFlow(const int startNode, const int endNode, const vector<int>& parentTable, const vector<int>& parentEdgeIndex, const vector<double>& maxFlow) {
  
  double totalFlow = min(maxFlow[endNode], sinkCapacity[endNode] - sinkFlow[endNode]);
  
  // add the flow
  sourceFlow[startNode] += totalFlow;
  sinkFlow[endNode] += totalFlow;
  if(startNode == endNode) {
    return totalFlow; // nothing to do anymore
  }
  int curNode = endNode;
  int parentEdge;
  int parentNode;
  GraphEdge* e;
  while(curNode != startNode) {
    parentNode = parentTable[curNode];
    parentEdge = parentEdgeIndex[curNode];
    e = nodes[parentNode][parentEdge];
    
    double restFlow = totalFlow;
    if(e->backwards->flow > 0) {
      if(e->backwards->flow < totalFlow) {
        restFlow -= e->backwards->flow;
        e->backwards->flow = 0;
        e->flow += restFlow;
      }
      else {
        e->backwards->flow -= totalFlow;
      }
    }
    else {
      e->flow += totalFlow;
    }
    
    curNode = parentNode;
  }
  
  return totalFlow;
}



void Graph::printGraph(ostream& outStream) const
{
  Node::const_iterator edgeIt; // iterator over the edges
  list<int>::const_iterator listIt;
  
  // print the numbers for the source and sink
  for(listIt = maxFlowNodeList.begin(); listIt!= maxFlowNodeList.end(); ++listIt) {
    if(sourceCapacity[*listIt] > 0) {
      outStream << *listIt << ": Source Cap: " << sourceCapacity[*listIt] << " Flow: " <<
        sourceFlow[*listIt] << endl;
    }    
    if(sinkCapacity[*listIt] > 0) {
      outStream << *listIt << ": Sink Cap: " << sinkCapacity[*listIt] << " Flow: " <<
        sinkFlow[*listIt] << endl;
    }    }
  
  for(unsigned int node = 0; node < nodes.size(); ++node)
  {
    outStream << "Node Number: " << node << endl;
    outStream << "Dist: " << dist[node] << endl;
    outStream << "ExFlow: " << exFlow[node] << endl;
    outStream << "Edges:" << endl;
    for(edgeIt = nodes[node].begin(); edgeIt != nodes[node].end(); ++edgeIt)
    {
      outStream << "To: " << (*edgeIt)->to << " Cap: " << (*edgeIt)->capacity <<
        " Flow: " << (*edgeIt)->flow  << endl;
    }
    outStream << endl;
  }
  outStream << endl;
}


void Graph::constructorWorker(const vector<list<int> >& conn, const vector<list<double> >& weights) {
  unsigned int numberOfNodes = conn.size();
  list<int> connOneNode;
  list<double> weightsOneNode;
  unsigned int numOfConn;
  unsigned int node1, node2;
  double weight;
  
  nodes.clear(); nodes.resize(numberOfNodes);
  sourceCapacity.clear(); sourceCapacity.resize(numberOfNodes, 0);
  sourceFlow.clear(); sourceFlow.resize(numberOfNodes, 0);
  sinkCapacity.clear(); sinkCapacity.resize(numberOfNodes, 0);
  sinkFlow.clear(); sinkFlow.resize(numberOfNodes, 0);
  subNodesICG.clear(); subNodesICG.resize(numberOfNodes, false);
  subNodesCWSM.clear(); subNodesCWSM.resize(numberOfNodes, false);
  subNodesCT.clear(); subNodesCT.resize(numberOfNodes, false);
  subNodes.clear(); subNodes.resize(numberOfNodes, false);
  exFlow.clear(); exFlow.resize(numberOfNodes, 0);
  dist.clear(); dist.resize(numberOfNodes, 0);
  activeByDist.clear();
  level = 0;
  maxFlowNodeList.clear();
  
  // go through all the nodes
  for(unsigned int i=0; i<numberOfNodes; ++i)
  {
    connOneNode = conn[i]; // get the connection of the current node
    weightsOneNode = weights[i]; // and with the weights
    numOfConn = connOneNode.size(); // how many are there
    node1 = i; // number of node currently working on
    for(unsigned int j=0; j<numOfConn; ++j) // go through all connections of this node
    {
      node2 = connOneNode.front();
      connOneNode.pop_front();
      weight = weightsOneNode.front();
      weightsOneNode.pop_front();
      
      // only check nodes that have a larger number (to avoid doubles and nodes to itself)
      if(node2>node1)
      {
        addEdge(node1, node2, weight, weight);
      }
    }
  }
}


Graph::Graph(Rcpp::List connList, Rcpp::List connWeights)
{
  Rcpp::IntegerVector connOneNode;
  Rcpp::NumericVector weightOneNode;    
  int numberOfNodes = connList.length();
  int numOfConn;
  int node1,node2;
  
  vector<list<int> > conn(numberOfNodes);
  vector<list<double> > weights(numberOfNodes);
  list<int> foo;
  list<double> bar;
  
  // go through all the nodes
  for(int i=0; i < numberOfNodes; ++i)
  {
    connOneNode = connList[i]; // get the connection of the current node
    weightOneNode = connWeights[i]; // and for the weights
    numOfConn = connOneNode.length(); // how many are there
    foo.clear();
    bar.clear();
    node1 = i; // number of node currently working on
    for(int j = 0; j < numOfConn; ++j) // go through all connections of this node
    {
      foo.push_back(connOneNode[j] - 1);
      bar.push_back(weightOneNode[j]);
      
    }
    conn[i] = foo;
    weights[i] = bar;
  }
  
  constructorWorker(conn, weights);
}

// Graph::Graph(const vector<list<int> >& conn, const vector<list<double> >& weights)
// {
//   constructorWorker(conn, weights);
// }

Graph::~Graph()
{
  // go through all nodes and delete all the edges 
  Nodes::iterator nodeIt;
  Node::iterator edgeIt;
  
  for(nodeIt = nodes.begin(); nodeIt != nodes.end(); ++nodeIt) {
    for(edgeIt = nodeIt->begin(); edgeIt != nodeIt->end(); ++edgeIt) {
      delete(*edgeIt);
      *edgeIt = NULL;
    }
  }
}

void Graph::addNode() {
  nodes.push_back(Node(0));
  subNodesICG.push_back(false);
  subNodesCWSM.push_back(false);
  subNodesCT.push_back(false);
  subNodes.push_back(false);
  sourceCapacity.push_back(0);
  sinkCapacity.push_back(0);
  sourceFlow.push_back(0);
  sinkFlow.push_back(0);
  exFlow.push_back(0);
  dist.push_back(0);
}



void Graph::makeCopy(const Graph& pg) {
  // make a deep copy of the penaltygraph
  Node::const_iterator edgeIt;
  GraphEdge* edgePtr;
  GraphEdge* edgePtrBack;
  
  // copy all the elements that can just be copied
  subNodesICG = pg.subNodesICG;
  subNodesCWSM = pg.subNodesCWSM;
  subNodes = pg.subNodes;
  subNodesCT = pg.subNodesCT;
  maxFlowNodeList = pg.maxFlowNodeList;
  sourceCapacity = pg.sourceCapacity;
  sourceFlow = pg.sourceFlow;
  sinkCapacity = pg.sinkCapacity;
  sinkFlow = pg.sinkFlow;
  exFlow = pg.exFlow;
  dist = pg.dist;
  activeByDist = pg.activeByDist;
  level = pg.level;
  
  // prepare nodes vector for copy
  nodes.clear(); nodes.resize(pg.nodes.size(), vector<GraphEdge*>());
  
  // go through all the edges 
  for(int node = 0; node < pg.nodes.size(); ++node) {
    for(edgeIt = pg.nodes[node].begin(); edgeIt != pg.nodes[node].end(); ++edgeIt) {
      if((*edgeIt)->to > node) {
        edgePtr = *edgeIt;
        edgePtrBack = edgePtr->backwards;
        addEdge(node, edgePtr->to, edgePtr->capacity, edgePtrBack->capacity, edgePtr->flow, edgePtrBack->flow); 
      }
    }
  }
}


Graph::Graph(const Graph& pg) {
  makeCopy(pg);
}


Graph& Graph::operator= (const Graph& pg) {
  makeCopy(pg);
  return(*this);
}

vector<list<int> > Graph::identifyConnectedGroups(const list<int>& x) {
  list<int> currentGroup;
  list<int> neighbours;
  vector<list<int> > allGroups(0);
  list<int> stack;
  list<int>::const_iterator listIt;
  list<int>::const_iterator neighIt;
  
  // set the subNodes that appear in the list
  for(listIt = x.begin(); listIt != x.end(); ++listIt) {
    subNodesICG[*listIt] = true;
  }
  
  for(listIt = x.begin(); listIt != x.end(); ++listIt) {
    if(subNodesICG[*listIt]) {
      subNodesICG[*listIt] = false;
      // first add an element of x to the new set
      currentGroup.clear();
      currentGroup.push_back(*listIt);
      stack.push_back(*listIt);
      // now for all new elements
      // add their neighbours that are also in x 
      // then remove them from x
      while(stack.size() > 0) {
        neighbours = connectedTo(stack);
        stack.clear();
        for(neighIt = neighbours.begin(); neighIt != neighbours.end(); ++neighIt) {
          if(subNodesICG[*neighIt]) {
            subNodesICG[*neighIt] = false;
            currentGroup.push_back(*neighIt); // in currentGroup, everything can only be inserted once
            stack.push_back(*neighIt);
          }
        }
      }
      allGroups.push_back(currentGroup);
    }
  }
  // Clean up: ensure all nodes in x have subNodesICG set to false
  clearSubNodesICG(x);
  return(allGroups);
}






list<int> Graph::connectedWithSameValue(int node, const vector<double>& values, double accuracy) {
  list<int> stack;
  list<int> stackConnect;
  list<int>::iterator listIt;
  list<int> connected;
  list<int> neighbours;
  double curVal = values[node];
  
  stack.push_back(node);
  connected.push_back(node);
  subNodesCWSM[node] = true;
  
  while(stack.size() > 0) {
    stackConnect = connectedTo(stack);
    stack.clear();
    for(listIt = stackConnect.begin(); listIt != stackConnect.end(); ++listIt) {
      // Check if values are equal (special handling for zero)
      bool valuesEqual = (curVal == 0 && values[*listIt] == 0) || 
                         (curVal != 0 && fabs(values[*listIt] - curVal) < accuracy/10);
      
      if(valuesEqual && !subNodesCWSM[*listIt]) {
        connected.push_back(*listIt);
        subNodesCWSM[*listIt] = true;
        stack.push_back(*listIt);
      }
    }
  }
  
  // Clean up the subnodes
  clearSubNodesCWSM(connected);
  return connected;
}


vector<list<int> > Graph::splitGroup(const list<int>& groupNodes, vector<double>& nodePull, double adjustment, double lambda2, bool invert) {
  vector<list<int> > foo(0);
  list<double> saveNodePulls;
  list<int> mfgNodes;
  list<int>::const_iterator listIt;
  list<double>::iterator saveIt;
  
  // adjust everything
  for(listIt = groupNodes.begin(); listIt != groupNodes.end(); ++listIt) {
    saveNodePulls.push_back(nodePull[*listIt]);
    nodePull[*listIt] += adjustment;
    if(invert) {
      nodePull[*listIt] *= -1;
    }      
    nodePull[*listIt] /= lambda2; 
  }
  
  initializeMaxFlow(groupNodes, nodePull); 
  findMaxFlowEdmondsKarp();
  mfgNodes = reachableFromSource();
  foo = identifyConnectedGroups(mfgNodes);
  
  // now reset the nodePulls to their old values
  for(listIt = groupNodes.begin(), saveIt = saveNodePulls.begin(); listIt != groupNodes.end(); ++listIt, ++saveIt) {
    nodePull[*listIt] = *saveIt;
  }
  // now the nodes need to be sorted
  return foo;
}

list<int> Graph::reachableFromSource()
{
  list<int> reachable;
  distance(true); // get the distance from the source; nodes.size() means that it can't be
  // reached
  // reachable nodes are all nodes in parent except source and sink
  list<int>::iterator listIt;
  for(listIt = maxFlowNodeList.begin(); listIt != maxFlowNodeList.end(); ++listIt) {
    if((int) dist[*listIt] <= maxFlowNodeList.size())
    {
      // node is reachable; insert node in external notation
      reachable.push_back(*listIt); 
    }
  }
  return(reachable);
}


list<int> Graph::getComplement(list<int> &x, list<int>& y)
{
  list<int> complement;
  x.sort();
  y.sort();
  list<int>::const_iterator xIt;
  list<int>::const_iterator yIt;
  
  xIt = x.begin();
  for(yIt = y.begin(); yIt != y.end(); ++yIt) {
    while(xIt != x.end() && *xIt < *yIt) {
      ++xIt;
    }
    if(xIt == x.end() || *yIt != *xIt) {
      complement.push_back(*yIt);
    }
  }
  
  return(complement);
}


void Graph::initializeMaxFlow(const list<int>& groupNodes, const vector<double>& nodePull) {
  list<int>::const_iterator listIt;
  Node::iterator edgeIt;
  
  for(listIt = maxFlowNodeList.begin(); listIt != maxFlowNodeList.end(); ++listIt) {
    subNodes[*listIt] = false;
  }
  
  maxFlowNodeList = groupNodes;
  // set the appropriate subNodes to true
  
  
  for(listIt = groupNodes.begin(); listIt != groupNodes.end(); ++listIt) {
    for(edgeIt = nodes[*listIt].begin(); edgeIt != nodes[*listIt].end(); ++edgeIt) {
      (*edgeIt)->flow = 0;
    }
    // now clean up everything w.r.t. source and sink
    sourceFlow[*listIt] = 0;
    sinkFlow[*listIt] = 0;
    exFlow[*listIt] = 0;
    dist[*listIt] = 0;
    subNodes[*listIt] = true;
    if(nodePull[*listIt] > 0) {
      sourceCapacity[*listIt] = nodePull[*listIt];
      sinkCapacity[*listIt] = 0;
    }
    else {
      sourceCapacity[*listIt] = 0;
      sinkCapacity[*listIt] = -nodePull[*listIt];
    }
    
  }
  
  activeByDist.clear();
  
  level = -1;
}


double Graph::getFlowFromSource() {
  double res = 0;
  list<int>::const_iterator listIt;
  for(listIt = maxFlowNodeList.begin(); listIt != maxFlowNodeList.end(); ++listIt) {
    res += sourceFlow[*listIt];
  }
  return res;
}

double Graph::getFlowIntoSink() {
  double res = 0;
  list<int>::const_iterator listIt;
  for(listIt = maxFlowNodeList.begin(); listIt != maxFlowNodeList.end(); ++listIt) {
    res += sinkFlow[*listIt];
  }
  return res;
}



list<int> Graph::allNodes()
{
  list<int> all;
  
  for(int i = 0; i < nodes.size(); ++i) {
    all.push_back(i);
  }
  
  return(all);
}


void Graph::getFusedConnectionsWeights(const vector<int>& fusions, int numFusions, vector<vector<int> > &connections, vector<vector<double> > &weights) {
  list<int>::iterator listIt;
  list<int>::iterator mapListIt;
  Node::iterator nodeIt;
  
  // prepare the parameters
  connections.clear(); connections.resize(numFusions);
  weights.clear(); weights.resize(numFusions);
  
  // create a mapping of fused groups to original groups
  vector<list<int> > mapFusedToOriginal(numFusions);
  for(int i = 0; i < fusions.size(); ++i) {
    mapFusedToOriginal[fusions[i]].push_back(i);
  }
  
  list<int> usedFusedGroups;
  vector<int> fusedGroupPosMap(numFusions, -1);
  int toGroup;
  // now build up the weights and connections of the fused groups
  for(int curGrp = 0; curGrp < numFusions; ++curGrp) {
    for(mapListIt = mapFusedToOriginal[curGrp].begin(); mapListIt != mapFusedToOriginal[curGrp].end(); ++mapListIt) {
      for(nodeIt = nodes[*mapListIt].begin(); nodeIt != nodes[*mapListIt].end(); ++nodeIt) {
        toGroup = fusions[(*nodeIt)->to];
        // check if it is to a different group
        if(toGroup != curGrp) {
          // check if it was found before
          if(fusedGroupPosMap[toGroup] == -1) { // add it
            fusedGroupPosMap[toGroup] = usedFusedGroups.size();
            usedFusedGroups.push_back(toGroup);
            connections[curGrp].push_back(toGroup);
            weights[curGrp].push_back(0);
          }
          // now add the weight in
          weights[curGrp][fusedGroupPosMap[toGroup]] += (*nodeIt)->capacity; 
        }
      }
    }
    
    // clean up
    for(listIt = usedFusedGroups.begin(); listIt != usedFusedGroups.end(); ++listIt) {
      fusedGroupPosMap[*listIt] = -1;
    }
    usedFusedGroups.clear();
  }
}


void Graph::makePullAdjustment(const vector<double>& beta, vector<double>& nodePull, double lambda2, double accuracy) {
  Node::iterator nodeIt;
  for(int pos = 0; pos < nodePull.size(); ++pos) {
    for(nodeIt = nodes[pos].begin(); nodeIt != nodes[pos].end(); ++nodeIt) {
      if(beta[pos] - beta[(*nodeIt)->to] > accuracy) {
        nodePull[pos] += lambda2 * (*nodeIt)->capacity;
      }
      else if(beta[pos] - beta[(*nodeIt)->to] < -accuracy) {
        nodePull[pos] -= lambda2 * (*nodeIt)->capacity;
      }
    }
  }
  
}




void Graph::printIntList(ostream& outStream, const list<int> x) {
  list<int>::const_iterator listIt;
  for(listIt = x.begin(); listIt != x.end(); ++listIt) {
    outStream << *listIt << ", ";
  }
  outStream << endl;
}

// ============================================================================
// HELPER FUNCTIONS FOR SUBNODE MANAGEMENT
// 

void Graph::initializeSubNodes() {
  int size = nodes.size();
  subNodesICG.assign(size, false);
  subNodesCWSM.assign(size, false);
  subNodesCT.assign(size, false);
  subNodes.assign(size, false);
}

void Graph::clearSubNodesICG(const list<int>& nodeList) {
  list<int>::const_iterator listIt;
  for(listIt = nodeList.begin(); listIt != nodeList.end(); ++listIt) {
    subNodesICG[*listIt] = false;
  }
}

void Graph::clearSubNodesCWSM(const list<int>& nodeList) {
  list<int>::const_iterator listIt;
  for(listIt = nodeList.begin(); listIt != nodeList.end(); ++listIt) {
    subNodesCWSM[*listIt] = false;
  }
}

void Graph::clearSubNodesCT(const list<int>& nodeList) {
  list<int>::const_iterator listIt;
  for(listIt = nodeList.begin(); listIt != nodeList.end(); ++listIt) {
    subNodesCT[*listIt] = false;
  }
}


