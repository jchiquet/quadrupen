#include "FusedLasso_class.h"

using namespace std;

FusedLasso::FusedLasso(const SparseMatrix &X, const vector<double> &y, const vector<double> &wObs, const vector<double> &beta, const vector<double>& wLambda1, const Graph& graph, int maxIterInner, int maxIterOuter, double accuracy, int maxActivateVars, double lambda1, double lambda2, regEnum regType) {
  
  // Store input data
  this->X = X;
  this->y = y;
  this->wObs = wObs;
  this->beta = beta;
  this->n = X.n;
  this->p = X.p;
  this->wLambda1 = wLambda1;
  this->regType = regType;
  this->pg = graph;
  this->lambda1 = lambda1;
  this->lambda2 = lambda2;
  
  // Store algorithm parameters
  this->maxIterInner = maxIterInner;
  this->maxIterOuter = maxIterOuter;
  this->accuracy = accuracy;
  this->maxActivateVars = maxActivateVars;
  
  // Initialize fusion information: all variables start unfused
  this->fusedGroupSize.resize(p, 1);
  this->fusions.resize(p);
  for(int i = 0; i < p; ++i) {
    this->fusions[i] = i;
  }
  
  // Initialize quadratic derivative object based on regression type
  switch(regType) {
  case GAUSSIAN: 
    quadratic = counted_ptr<QuadraticDerivative>(
        new QuadraticDerivativeDiagonal(this->X, this->y, this->wObs, this->beta));
    break;
  case BINOMIAL:
    quadratic = counted_ptr<QuadraticDerivative>(
        new QuadraticDerivativeLogistic(this->X, this->y, this->wObs, this->beta));
  }
  
  lastDetailedFusions = fusions;
  fusionLevel = 0;
}

FusedLasso::~FusedLasso() {}

// assumes that beta has been updated from the fusedBeta before -> translateFusedToOriginal()
vector<double> FusedLasso::getPulls() {
  vector<double> pulls(beta.size());
  for(unsigned int i=0; i < pulls.size(); ++i) {
    pulls[i] = quadratic->getDerivative(i);
  }
  return pulls;
}

// gives the current fusions based on the current betas (no splits will be performed)
void FusedLasso::getEqualFusions(vector<int>& newFusions, vector<int>& newFusedGroupSize, bool zeroSingle, double accFactor) {
  newFusions.clear(); newFusedGroupSize.clear();
  newFusions.resize(beta.size(), -1);
  list<int> nodes;
  list<int>::iterator listIt;
  int groupNum = 0;
  
  double accHere = accuracy * accFactor;
  
  for(unsigned int i = 0; i < beta.size(); ++i) {
    if(newFusions[i] == -1) {
      nodes = pg.connectedWithSameValue(i, beta, accHere);
      if(zeroSingle && beta[i]==0) {
        for(listIt = nodes.begin(); listIt != nodes.end(); ++listIt) {
          newFusions[*listIt] = groupNum;
          groupNum++;
          newFusedGroupSize.push_back(1);
        }
      }
      else {
        newFusedGroupSize.push_back(nodes.size());
        for(listIt = nodes.begin(); listIt != nodes.end(); ++listIt) {
          newFusions[*listIt] = groupNum;
        }
        ++groupNum;
      }
    }
  }
}

bool FusedLasso::areFusionsEqual(vector<int> &fusion1, vector<int> &fusion2) {
  if(fusion1.size() != fusion2.size()) {
    return false;
  }
  
  for(unsigned int i = 0; i < fusion1.size(); ++i) {
    if(fusion1[i] != fusion2[i]) {
      return(false);
    }
  }
  
  return true;
}

void FusedLasso::getSplitFusionsActive(vector<int>& newFusions, vector<int>& newFusedGroupSize) {
  
  newFusions.clear(); newFusedGroupSize.clear();
  newFusions.resize(p, -1);
  
  // first, we need to get the pulls correctly
  vector<double> nodePull = getPulls();
  // adjust the nodePulls appropriately
  makePullAdjustment(beta, nodePull, lambda2);
  vector<list<int> > groups;
  list<int> nodes;
  list<int>::iterator listIt;
  int newNumFused = 0;
  
  for(unsigned int i = 0; i < p; ++i) {
    if(newFusions[i] == -1) { // not yet grouped
      nodes = pg.connectedWithSameValue(i, beta, accuracy);
      if(beta[i] == 0) { // all variables are singles
        for(listIt = nodes.begin(); listIt != nodes.end(); ++listIt) {
          beta[*listIt] = 0; // sets some betas not more than accuracy different
          // from 0 back to zero
          newFusedGroupSize.push_back(1);
          newFusions[*listIt] = newNumFused;
          newNumFused++;
        }
        continue;
      }
      else if(beta[i] > 0) {
        groups = pg.splitGroup(nodes, nodePull, lambda1, lambda2, false);
      }
      else {
        groups = pg.splitGroup(nodes, nodePull, -lambda1, lambda2, false);
      }
      
      try {
        addComplementaryGroups(groups, nodes);
      }
      catch(std::exception &ex) {
        forward_exception_to_r(ex);
      }
      
      // now mark all the sets in here
      for(int j = 0; j < groups.size(); ++j) {
        
        // setting this here is not completely correct 
        newFusedGroupSize.push_back(groups[j].size());
        for(listIt = groups[j].begin(); listIt != groups[j].end(); ++listIt) {
          newFusions[*listIt] = newNumFused;
        }
        newNumFused++;
      }
      
    }
  }
}

void FusedLasso::getSplitFusionsInactive(vector<int>& newFusions, vector<int>& newFusedGroupSize) {
  
  newFusions.clear(); newFusedGroupSize.clear();
  newFusions.resize(p, -1);
  
  // first, we need to get the pulls correctly
  vector<double> nodePull = getPulls();
  // adjust the nodePulls appropriately
  makePullAdjustment(beta, nodePull, lambda2);
  vector<list<int> > groups;
  list<int> nodes;
  list<int>::iterator listIt;
  int newNumFused = 0;
  
  for(unsigned int i = 0; i < p; ++i) {
    if(newFusions[i] == -1) { // not yet grouped
      nodes = pg.connectedWithSameValue(i, beta, accuracy);
      if(beta[i] == 0) {
        vector<list<int> > foo;
        groups = pg.splitGroup(nodes, nodePull, -lambda1, lambda2, false);
        foo = pg.splitGroup(nodes, nodePull, lambda1, lambda2, true);
        groups.insert(groups.end(), foo.begin(), foo.end());
      }
      else { // leave everything the same as before
        groups.resize(1);
        groups[0] = nodes;
      }
      
      try {
        addComplementaryGroups(groups, nodes);
      }
      catch(std::exception &ex) {
        forward_exception_to_r(ex);
      }
      
      // now mark all the sets in here
      for(int j = 0; j < groups.size(); ++j) {
        // setting this here is not completely correct 
        newFusedGroupSize.push_back(groups[j].size());
        for(listIt = groups[j].begin(); listIt != groups[j].end(); ++listIt) {
          newFusions[*listIt] = newNumFused;
        }
        newNumFused++;
      }
    }
  }
}


void FusedLasso::sortAllLists(vector<list<int> >& x) {
  for(unsigned int i = 0; i < x.size(); ++i) {
    x[i].sort();
  }
}


bool listComp(const list<int>& list1, const list<int>& list2) {
  int val1, val2;
  val1 = list1.front();
  val2 = list2.front();
  return val1 < val2;
}

void FusedLasso::addComplementaryGroups(vector<list<int> >& groups, list<int>& nodes) {
  list<int> complNodes = nodes;
  vector<list<int> > foo;
  vector<list<int> > groupsCopy = groups;
  vector<list<int> > groupsSave = groups;
  list<int> groupsMerged;
  list<int>::iterator listIt;
  
  // first find the complementary nodes
  for(unsigned int i = 0; i < groups.size(); ++i) {
    groupsMerged.splice(groupsMerged.end(), groupsCopy[i]);
  }
  
  groupsMerged.sort(); 
  int groupsMergedSize = groupsMerged.size();
  groupsMerged.unique();
  if(groupsMergedSize > groupsMerged.size()) {
    throw std::out_of_range("A node was grouped twice..."); 
  }
  
  complNodes = pg.getComplement(groupsMerged, nodes);
  
  // now make the rest into groups and add it to the groups
  foo = pg.identifyConnectedGroups(complNodes);
  groups.insert(groups.end(), foo.begin(), foo.end()); 
  sortAllLists(groups);
  sort(groups.begin(), groups.end(), listComp); 
}

vector<double> FusedLasso::translateOriginalToFused() {
  vector<double> fusedBeta(fusedGroupSize.size());
  for(unsigned int i = 0; i < fusions.size(); ++i) {
    fusedBeta[fusions[i]] = beta[i];
  }
  return fusedBeta;
}

bool FusedLasso::identifyNewFusionsHuber() {
  vector<int> newFusions;
  vector<int> newFusedGroupSize;
  
  double huberParam = 1000;
  
  // create the object to run the coordinate descent algorithm
  vector<int> singleFusions(p);
  for(unsigned int i = 0; i < p; ++i) {
    singleFusions[i] = i;
  }
  
  vector<vector<int> > connectionsSingle;
  vector<vector<double> > wLambda2Single;
  
  pg.getFusedConnectionsWeights(singleFusions, p, connectionsSingle, wLambda2Single);
  
  // now first run Huber, then normal
  FusedLassoCoordinate flcHuber(quadratic, wLambda1, connectionsSingle, wLambda2Single, 100, accuracy, 100000, lambda1, lambda2, Huber, huberParam);
  bool resultHuber = flcHuber.runAlgorithm();
  
  FusedLassoCoordinate flc(quadratic, wLambda1, connectionsSingle, wLambda2Single, maxIterInner, accuracy, maxActivateVars, lambda1, lambda2, L1);
  bool result = flc.runAlgorithm();
  beta = flc.getBetaOriginal(singleFusions);
  
  getEqualFusions(newFusions, newFusedGroupSize);
  
  if(areFusionsEqual(fusions, newFusions)) {
    return(false);
  }
  else {
    fusions = newFusions;
    fusedGroupSize = newFusedGroupSize;
    return(true);
  }
  
}


bool FusedLasso::identifyNewFusions(double lastIterChange, const int maxFusionLevel) {
  vector<int> newFusions;
  vector<int> newFusedGroupSize;
  
  // If no fusion levels allowed, do nothing
  if(maxFusionLevel == -1) {
    return(false);
  }
  
  // Update fusion level based on convergence
  if(lastIterChange > accuracy) {
    fusionLevel = 0; // Reset if not converged
  }
  else {
    ++fusionLevel; // Increment if converged
  }
  
  // Try fusions at different levels
  if(fusionLevel == 0 && fusionLevel <= maxFusionLevel) {
    // Level 0: Identify variables with equal values
    getEqualFusions(newFusions, newFusedGroupSize);
    if(areFusionsEqual(fusions, newFusions)) {
      ++fusionLevel; // Move to next level if no change
    }
  }
  
  if(fusionLevel == 1 && fusionLevel <= maxFusionLevel) {
    // Level 1: Split inactive variables
    getSplitFusionsInactive(newFusions, newFusedGroupSize);
    if(areFusionsEqual(fusions, newFusions)) {
      ++fusionLevel; // Move to next level if no change
    }
  }
  
  if(fusionLevel == 2 && fusionLevel <= maxFusionLevel) {
    // Level 2: Split active variables
    getSplitFusionsActive(newFusions, newFusedGroupSize);
    if(areFusionsEqual(fusions, newFusions)) {
      ++fusionLevel; // Move to next level if no change
    }
  } 
  
  // Check if we should stop
  if(fusionLevel > maxFusionLevel || fusionLevel == 3) {
    return(false); // Convergence reached
  }
  else {
    fusions = newFusions;
    fusedGroupSize = newFusedGroupSize;
    return(true); // New fusions found
  }
}


bool FusedLasso::runFused(penEnum penType) {
  
  vector<double> betaFused = translateOriginalToFused();
  vector<vector<int> > fusedGroups(fusedGroupSize.size());
  for(unsigned int i = 0; i < fusions.size(); ++i) {
    fusedGroups[fusions[i]].push_back(i);
  }
  
  SparseMatrix fusedX  = X.createFusedX(fusedGroups);
  vector<double> wLambda1Coordinate(fusedGroupSize.size());
  for(unsigned int i = 0; i < wLambda1Coordinate.size(); ++i) {
    wLambda1Coordinate[i] = 0;
    for(int j = 0; j < fusedGroups[i].size(); ++j) {
      wLambda1Coordinate[i] += wLambda1[fusedGroups[i][j]];
    }
  }
  
  vector<vector<int> > connectionsFused;
  vector<vector<double> > wLambda2Fused;
  
  pg.getFusedConnectionsWeights(fusions, fusedGroupSize.size(), connectionsFused, wLambda2Fused);
  
  counted_ptr<QuadraticDerivative> quadDer(new QuadraticDerivativeDiagonal(fusedX, y, wObs, betaFused));
  FusedLassoCoordinate flc(quadDer, wLambda1Coordinate, connectionsFused, wLambda2Fused, maxIterInner, accuracy/10, maxActivateVars, lambda1, lambda2, penType);
  
  bool result = flc.runAlgorithm();
  innerIterNum += flc.getIterNum();
  
  // now translate the result back
  beta = flc.getBetaOriginal(fusions);
  for(unsigned int i = 0; i< beta.size(); ++i) {
    quadratic->updateBeta(i, beta[i]);
  }
  outerIterNum++;
  return result;
  
}


bool FusedLasso::runFusedGeneral(penEnum penType) {
  
  if(outerIterNum > maxIterOuter) {
    return(false);
  }
  
  vector<double> betaFused = translateOriginalToFused();
  vector<vector<int> > fusedGroups(fusedGroupSize.size());
  for(unsigned int i = 0; i < fusions.size(); ++i) {
    fusedGroups[fusions[i]].push_back(i);
  }

  SparseMatrix fusedX = X.createFusedX(fusedGroups);
  vector<double> wLambda1Fused(fusedGroupSize.size());
  for(unsigned int i = 0; i < fusedGroups.size(); ++i) {
    wLambda1Fused[i] = 0;
    for(int j = 0; j < fusedGroups[i].size(); ++j) {
      wLambda1Fused[i] += wLambda1[fusedGroups[i][j]];
    }
  }
  
  vector<vector<int> > connectionsFused;
  vector<vector<double> > wLambda2Fused;
  
  pg.getFusedConnectionsWeights(fusions, fusedGroupSize.size(), connectionsFused, wLambda2Fused);
  
  vector<double> oldBeta;
  double maxBetaChange = accuracy + 1;
  counted_ptr<QuadraticDerivative> quadDer;
  counted_ptr<FusedLassoCoordinate> flc;
  bool result = true;
  
  while(maxBetaChange > accuracy && result && outerIterNum <= maxIterOuter) {
    oldBeta = betaFused;
    quadDer = counted_ptr<QuadraticDerivative>(new QuadraticDerivativeLogistic(fusedX, y, wObs, betaFused));
    if(quadDer -> isExtreme()) {
      return(false);
    }
    
    flc = counted_ptr<FusedLassoCoordinate> (new FusedLassoCoordinate(quadDer, wLambda1Fused, connectionsFused, wLambda2Fused, maxIterInner, accuracy/10 , maxActivateVars, lambda1, lambda2, penType));
    result =  flc->runAlgorithm();
    betaFused = flc->getBeta();
    innerIterNum += flc->getIterNum();

    maxBetaChange = maxDiffDoubleVec(oldBeta, betaFused);
    outerIterNum++;

  }
  
  beta = flc->getBetaOriginal(fusions);
  
  quadratic = counted_ptr<QuadraticDerivative>(new QuadraticDerivativeLogistic(X, y, wObs, beta)); 
  
  if(quadDer -> isExtreme()) {
    return(false);
  }
  
  return (result && outerIterNum < maxIterOuter);
}

bool FusedLasso::runAlgorithmL2() {
  outerIterNum = 0;
  innerIterNum = 0;
  bool lastRunOK = true;
  // set single fusions
  fusedGroupSize.resize(p, 1);
  for(unsigned int i = 0; i < p; ++i) {
    fusions[i] = i;
  }
  
  if(regType != GAUSSIAN) {
    lastRunOK = runFusedGeneral(L2);
  } else {
    lastRunOK = runFused(L2);
  }
  
  if(lastRunOK) {
    return true;
  } else {
    return false;
  }
}


bool FusedLasso::runAlgorithmHuber() {
  
  outerIterNum = 0;
  innerIterNum = 0;
  bool lastRunOK = true;
  double lastIterChange;
  
  vector<double> oldBeta;
  
  if(regType != GAUSSIAN) {
    lastRunOK = runFusedGeneral(L1);
  } else {
    lastRunOK = runFused(L1);
  }
  
  // check if we have an extreme result now
  if(quadratic->isExtreme()) {
    return false;
  }
  
  //    while(outerIterNum < maxIterOuter && lastRunOK) {
  oldBeta = beta;
  
  bool newFusionDifferent = identifyNewFusionsHuber();
  
  vector<int> newFusions = fusions;
  vector<int> newFusedGroupSize = fusedGroupSize;
  
  do { 
    oldBeta = beta;
    fusions = newFusions;
    fusedGroupSize = newFusedGroupSize;
    if(regType != GAUSSIAN) {
      lastRunOK = runFusedGeneral(L1);
    }
    else {
      lastRunOK = runFused(L1);
    }
    
    // check if we have an extreme result now
    if(quadratic->isExtreme()) {
      return false;
    }
    getEqualFusions(newFusions, newFusedGroupSize);
    
    fusions = newFusions;
    fusedGroupSize = newFusedGroupSize;
    lastIterChange = maxDiffDoubleVec(oldBeta, beta);
  } while(!areFusionsEqual(fusions, newFusions) && lastIterChange > accuracy);
  
  if(lastRunOK && outerIterNum < maxIterOuter) {
    return true;
  }
  else {
    return false;
  }
}

bool FusedLasso::runAlgorithm(const int maxFusionLevel) {
  outerIterNum = 0;
  innerIterNum = 0;
  bool lastRunOK = true;
  bool newFusionDifferent;
  vector<double> oldBeta; 
  fusionLevel = 0;
  
  while(outerIterNum < maxIterOuter && lastRunOK) {
    oldBeta = beta;
    
    // Run optimization on current fusion structure
    if(regType != GAUSSIAN) {
      lastRunOK = runFusedGeneral(L1);
    } else {
      lastRunOK = runFused(L1);
    }
    
    // Check for extreme results
    if(quadratic->isExtreme()) {
      lastRunOK = false;
      break;
    }
    
    // Try to identify new fusions
    double lastIterChange = maxDiffDoubleVec(oldBeta, beta);
    newFusionDifferent = identifyNewFusions(lastIterChange, maxFusionLevel);
    
    if(!newFusionDifferent) {
      // No new fusions found: algorithm has converged
      break;
    }
  }
  
  return (lastRunOK && outerIterNum < maxIterOuter);
}

SparseMatrix FusedLasso::runAlgorithm(
    penEnum penType,
    const int maxFusionLevel, 
    const vector<double>& lambda1Vec, 
    const vector<double>& lambda2Vec, 
    const int maxNonZero, 
    vector<bool>& success, 
    vector<int>& outerIterNumVec, 
    vector<int>& innerIterNumVec, 
    bool verbose) {
  
  SparseMatrix betaSols(beta.size());
  success.clear(); success.resize(lambda1Vec.size());
  outerIterNumVec.clear(); outerIterNumVec.resize(lambda1Vec.size());
  innerIterNumVec.clear(); innerIterNumVec.resize(lambda1Vec.size());
  
  for(unsigned int i = 0; i < lambda1Vec.size(); ++i) {
    setNewLambdas(lambda1Vec[i], lambda2Vec[i]);
    switch(penType) {
    case L1: if(maxFusionLevel >= -1) {
      success[i] = runAlgorithm(maxFusionLevel);        
    }
    else {
      success[i] = runAlgorithmHuber();
    }
    break;
    case L2: success[i] = runAlgorithmL2(); break;
    case Huber: success[i] = runAlgorithmHuber(); break;
    }
    outerIterNumVec[i] = getOuterIterNum();
    innerIterNumVec[i] = getInnerIterNum();
    betaSols.addColumn(getBeta());
    if(numNonZero(getBeta()) > maxNonZero) {
      break; // number really saved can be read from the returned matrix
    }
    if(!success[i]) {
      break; // algorithm did not finish successfully
    }
  }
  return betaSols;
}

void FusedLasso::setNewLambdas(double lambda1, double lambda2) {
  // when setting the new lambdas, only the lambdas
  // have to be reset; no other changes are necessary
  this->lambda1 = lambda1;
  this->lambda2 = lambda2;
}

double FusedLasso::findMaxLambda1(const list<int>& exemptVars) {
  vector<vector<int> > connectionsEmpty(p,vector<int>(0));
  vector<vector<double> > wLambda2Empty(p, vector<double>(0));
  vector<double> wLambda1Extreme(p, 1e6);
  list<int>::const_iterator listIt;
  for(listIt = exemptVars.begin(); listIt != exemptVars.end(); ++listIt) {
    wLambda1Extreme[*listIt] = 1e-4;
  }
  
  FusedLassoCoordinate flc(quadratic, wLambda1Extreme, connectionsEmpty, wLambda2Empty, maxIterInner, accuracy, maxActivateVars, 1, 1, L1);
  
  bool result = flc.runAlgorithm();
  
  this->beta = flc.getBeta();
  
  if(regType != GAUSSIAN) {
    quadratic = counted_ptr<QuadraticDerivative>(new QuadraticDerivativeLogistic(this->X, this->y, this->wObs, this->beta));
  }
  
  double highLambda1 = 0;
  double foo;
  // now search through the 0 variables which will start next
  for(unsigned int i = 0; i < p; ++i) {
    if(quadratic->getBeta(i) == 0) {
      foo = fabs(quadratic->getDerivative(i) / wLambda1[i]);
      if(highLambda1 < foo) {
        highLambda1 = foo;
      }
    }
  }
  return highLambda1;
}


// identify lambda2
double FusedLasso::findMaxLambda2(double lambda1_here) {
  if(p > 100000) {
    REprintf("Too big to estimate lambda2; returning 1");
    return n;
  }
  
  outerIterNum=0; // reset as this is run several times
  // first, set beta equal to 0.1
  for(unsigned int i = 0; i < p; i++) {
    beta[i] = 0.1;
    quadratic->updateBeta(i, 0.1);
  }
  
  vector<int> newFusions;
  vector<int> newFusedGroupSize;
  
  lambda1 = lambda1_here;
  lambda2 = n * 2e6;
  setNewLambdas(lambda1, lambda2);
  getEqualFusions(newFusions, newFusedGroupSize, false, 0.001/accuracy);
  runAlgorithmHuber();
  
  vector<double> oldBeta(p,0);
  for(unsigned int i=0; i < p; ++i) {
    oldBeta[i] = beta[i];
  }
  bool hasChanged;
  hasChanged=false;
  
  lambda2/=2;
  for(unsigned int i = 0; i < 30; ++i) {
    setNewLambdas(lambda1, lambda2);
    runAlgorithmHuber();
    for(unsigned int i = 0; i < p; ++i) {
      if(fabs(oldBeta[i] - beta[i]) > accuracy * 10) {
        hasChanged = true;
      }
    }
    if(hasChanged) {
      return(lambda2*2);
    } 
    else {
      lambda2 /= 2;
    }
  } 
  return(lambda2);
}


void FusedLasso::getBeta(SparseMatrix& betaMat) {
  betaMat.addColumn(beta);
}

void FusedLasso::printBetaAndGroup(ostream& outStream) {
  outStream << "========== Fusions =================" << endl;
  
  int curGroup = -1;
  for(unsigned int i = 0; i < p; ++i) {
    if(fusions[i] > curGroup) { // a new group has started
      curGroup++;
      outStream << " Group: " << i << " Beta: " << beta[i] ;
    }
  }
  
  outStream << endl << "========== End Fusions =============" << endl;
}

void FusedLasso::makePullAdjustment(const vector<double>& beta, vector<double>& nodePull, double lambda2) {
  pg.makePullAdjustment(beta, nodePull, lambda2, accuracy);
}

double FusedLasso::calcGroupAverage(int pos, vector<double>& nodePull, vector<double>& beta) {
  double curBeta = beta[pos];
  double average = 0;
  double groupSize = 0;
  while(fabs(beta[pos] - curBeta) < accuracy) {
    groupSize++;
    average += nodePull[pos];
    ++pos;
  }
  return average/groupSize;
} 

void FusedLasso::printGroups(vector<list<int> > x, ostream& outStream) {
  list<int>::iterator listIt;
  for(unsigned int i = 0; i < x.size(); ++i) {
    for(listIt = x[i].begin(); listIt != x[i].end(); ++listIt) {
      outStream << "i: " << i << " : " << *listIt << endl;
    }
  }
}

void FusedLasso::checkSolution(ostream& outStream) {
  vector<double> nodePull = getPulls();
  vector<double> nodePullOriginal = nodePull;
  vector<double> nodePullAdj = nodePull;
  makePullAdjustment( beta, nodePullAdj, lambda2);
  vector<int> newFusions;
  vector<int> newFusedGroupSize;
  
  // make all groups that are equal
  getEqualFusions(newFusions, newFusedGroupSize);
  
  // first adjust those for lambda1
  int curPos = 0;
  int curGroup = newFusions[curPos];
  double adjustment;
  
  outStream << "Lambda1: " << lambda1 << endl;
  outStream << "Lambda2: " << lambda2 << endl;
  
  while(curPos < nodePull.size()) {
    
    if(beta[curPos] > 0) {
      nodePull[curPos] += lambda1;
      curPos++;
    }
    else if(beta[curPos] < 0) {
      nodePull[curPos] -= lambda1;
      curPos++;
    }
    else { // is 0, so just subtract the average of the group
      curGroup = newFusions[curPos];
      adjustment = calcGroupAverage(curPos, nodePullAdj, beta);
      if(fabs(adjustment) > lambda1) {
        outStream << "=============== Problem =================" << endl;
        outStream << "Group " << curGroup << " at position " << curPos <<
          " has adjustment " << adjustment << endl;
      }
      while(newFusions[curPos] == curGroup && curPos < nodePull.size()) {
        nodePull[curPos] -= adjustment;
        curPos++;
      }
    }
  }
  
  // now run the maximum flow graph and check that 
  list<int> allNodes = pg.allNodes();
  // now adjust nodePull by dividing by lambda2
  for(unsigned int i = 0; i < nodePull.size(); ++i) {
    nodePull[i] /= lambda2;
    outStream << i << ": " << nodePull[i] << endl;
  }
  
  pg.initializeMaxFlow(allNodes, nodePull);
  pg.findMaxFlowEdmondsKarp();
  GraphEdge* e, *eBack;
  
  // now print out a compact version of the solution
  for(unsigned int i = 0; i < beta.size(); ++i) {
    e = pg.nodes[i].back();
    eBack = e->backwards;
    outStream << "Node " << i << " Beta: " << beta[i] << " Deriv: " << nodePullOriginal[i] <<
      " DerivAdj: " << nodePull[i];
    if(i < beta.size()-1) {
      outStream << "t[i,i+1]: " << e->flow << "Back: " << eBack->flow;
    }
    outStream << endl;
  }    
  
  outStream << "================= Whole MaxflowGraph ==================" << endl;
  //    pg.printGraph(cout);
  outStream << "==================== END WHOLE =======================" << endl;
}

