#include "FusedLasso_optim.h"

using namespace std;

FusedLassoCoordinate::FusedLassoCoordinate(counted_ptr<QuadraticDerivative> quadDer, const vector<double>& wLambda1, const vector<vector<int> > &connections, const vector<vector<double> > &wLambda2, int maxIterInner, double accuracy, int maxActivateVars, double lambda1, double lambda2, penEnum penType, double huberParam) : quadratic(quadDer) {
 
    // Set parameters
    this->p = quadDer->getp();
    this->n = quadDer->getn();
    this->lambda1 = lambda1;
    this->lambda2 = lambda2;
    this->accuracy = accuracy;
    this->huberParam = huberParam;
    this->penType = penType;
    this->maxIterInner = maxIterInner;
    this->maxActivateVars = maxActivateVars;
    this->beta = quadDer->getBetaVec();
    
    // Pre-multiply wLambda1 by lambda1
    this->wLambda1Mult.resize(wLambda1.size());
    for(size_t i = 0; i < wLambda1.size(); ++i) {
        wLambda1Mult[i] = wLambda1[i] * lambda1;
    }
    
    // Initialize active variables (all non-zero beta)
    active.clear();
    active.reserve(min(n, p));
    for(int i = 0; i < p; ++i) {
        if(beta[i] != 0) {
            active.push_back(i);
            quadDer->activate(i);
        }
    }
    
    // Store connections and pre-multiply wLambda2 by lambda2
    this->connections = connections;
    this->wLambda2Mult = wLambda2;
    for(size_t i = 0; i < wLambda2Mult.size(); ++i) {
        for(size_t j = 0; j < wLambda2Mult[i].size(); ++j) {
            wLambda2Mult[i][j] *= lambda2;
        }
    }
    
    quadratic = quadDer;
}



bool FusedLassoCoordinate::activateBeta(int pos) {
    // first check if the element is already included
    vector<int>::iterator vecIt = find(active.begin(), active.end(), pos);
    if(vecIt == active.end()) {// could not find the element
        active.push_back(pos);
        quadratic->activate(pos);
        return true;
    }
    return false;
}


bool FusedLassoCoordinate::deactivateBeta(int pos) {
    vector<int>::iterator vecIt = find(active.begin(), active.end(), pos);
    if(vecIt != active.end()) {
        active.erase(vecIt); // delete the element
        return true;
    }
    return false;
}

void FusedLassoCoordinate::singleStep(int pos) {
    
    // here, view the problem as only having variable pos
    // find the optimal position and move beta there
    double deriv = quadratic->getDerivative(pos);
    double slope = quadratic->getHessian(pos); 

    switch(penType) {
        case L1: Steps::findRootL1(beta, deriv, slope, pos, wLambda1Mult[pos], connections[pos], wLambda2Mult[pos]); break;
        case Huber: Steps::findRootHuber(beta, deriv, slope, pos, wLambda1Mult[pos], connections[pos], wLambda2Mult[pos], huberParam); break;
        case L2: Steps::findRootL2(beta, deriv, slope, pos, wLambda1Mult[pos], connections[pos], wLambda2Mult[pos]); break;
    }
    quadratic->updateBeta(pos, beta[pos]);
}

double FusedLassoCoordinate::singleIteration(const int iterNum) {
    // Return early if no active variables
    if(active.size() == 0) {
        return 0;
    }

    // Perform single steps for all active variables
    // and track the maximum change
    vector<int>::iterator vecIt;
    double maxNorm = 0;
    
    for(vecIt = active.begin(); vecIt != active.end(); ++vecIt) {
        double oldBeta = beta[*vecIt];
        singleStep(*vecIt);
        
        double change = fabs(beta[*vecIt] - oldBeta);
        if(change > maxNorm) { 
            maxNorm = change; 
        }
    }

    return maxNorm;
}


int FusedLassoCoordinate::activateVariables() {
    // go through all variables that are equal to 0
    // calculate their penalty adjusted derivative and sort them
    // by absolute value of derivative
    // then active the top variables, but at most maxVars
    multimap<double, int> activationCandidates;
    vector<double> derivVec = quadratic->getDerivativeVec(); 
    vector<bool> isActive(beta.size(), false);
    for(unsigned int i = 0; i < active.size(); ++i) {
        isActive[active[i]] = true;
    }

    double adjDeriv;
    for(unsigned int pos = 0; pos < beta.size(); ++pos) {
        if(!isActive[pos]) {
            double zeroPenalty = wLambda1Mult[pos];

            adjDeriv = derivVec[pos];
            derivAdjustment(adjDeriv, zeroPenalty, pos);
            // now if in absolute value > lambda1, add to list of
            // possible activations
            if(fabs(adjDeriv) > zeroPenalty) {
                activationCandidates.insert(make_pair(fabs(adjDeriv), pos));
            }
        }
    }

    // now go through the activation candidates
    multimap<double, int>::reverse_iterator mapIt;
    int activateCount;
    for(mapIt = activationCandidates.rbegin(), activateCount = 0; mapIt !=activationCandidates.rend() && activateCount < maxActivateVars; ++mapIt, ++ activateCount) {
        active.push_back(mapIt->second);       
        quadratic->activate(mapIt->second); 
    }

    // and sort the vector of active variables
    sort(active.begin(), active.end());

    return activateCount;
}

void FusedLassoCoordinate::derivAdjustment(double& adjDeriv, double& zeroPenalty, int pos) {
    // Adjust derivative based on penalty type and connected variables
    
    switch(penType) {
        case L1:
            // L1 penalty adjustment
            for(size_t i = 0; i < connections[pos].size(); ++i) {
                double connBeta = beta[connections[pos][i]];
                if(connBeta > 0) {
                    adjDeriv -= wLambda2Mult[pos][i];
                }
                else if(connBeta < 0) {
                    adjDeriv += wLambda2Mult[pos][i];
                }
                else {
                    zeroPenalty += wLambda2Mult[pos][i];
                }
            }
            break;
            
        case Huber:
            // Huber penalty adjustment
            for(size_t i = 0; i < connections[pos].size(); ++i) {
                double connBeta = beta[connections[pos][i]];
                double threshold = 1.0 / huberParam;
                if(connBeta > threshold) {
                    adjDeriv -= wLambda2Mult[pos][i];
                }
                else if(connBeta < -threshold) {
                    adjDeriv += wLambda2Mult[pos][i];
                }
                else {
                    // Smooth region: quadratic penalty
                    adjDeriv -= huberParam * connBeta * wLambda2Mult[pos][i];
                }
            }
            break;
            
        case L2:
            // L2 penalty adjustment
            for(size_t i = 0; i < connections[pos].size(); ++i) {
                adjDeriv -= beta[connections[pos][i]] * 2 * wLambda2Mult[pos][i];
            }
            break;
    }
}

int FusedLassoCoordinate::deactivateVariables() {
    // go through all active variables and remove the
    // ones from being active that are equal to 0

    int deleteCount = 0;
    vector<int>::iterator vecIt;
    for(vecIt = active.begin(); vecIt != active.end(); ) {
       if(beta[*vecIt] == 0) { // yes, deactivate
          vecIt = active.erase(vecIt); // due to deletion don't need to move one element ahead
          deleteCount++;
       }
       else {
           vecIt++; 
       }
    }
    return deleteCount;
}

bool FusedLassoCoordinate::runAlgorithm() {
    // Main optimization loop with variable activation
    bool finished = false;
    int newActivations;
    iterNum = 0;
    double curError;

    while(!finished && iterNum < maxIterInner) {
        curError = 1;
        
        // Inner loop: optimize over active variables until convergence
        while(curError > accuracy && iterNum < maxIterInner) {
            curError = singleIteration(iterNum);
            R_CheckUserInterrupt(); // Allow user to cancel computation

            iterNum++;
            if(curError > 1e5) {
                return(false); // Divergence detected
            }
        }

        // Try to activate new variables
        newActivations = activateVariables();

        if(newActivations == 0) {
            // No new activations: algorithm has converged
            finished = true;
        }
        if(iterNum >= maxIterInner) {
            // Maximum iterations reached
            finished = false;
        }
    }

    return(finished);
}

int FusedLassoCoordinate::getIterNum() const {
    return iterNum;
}

vector<double> FusedLassoCoordinate::getBetaOriginal(vector<int>& fusions) {
    vector<double> betaOrig(fusions.size());

    for(unsigned int i = 0; i < fusions.size(); ++i) {
        if(fusions[i] > beta.size() - 1) {
            Rcpp::stop("The node has a group that is too large. Node: %i", i, "\n") ;
        }
        betaOrig[i] = beta[fusions[i]];
    }
    return betaOrig;
}

void FusedLassoCoordinate::printBetaActive(ostream& outStream) {
    outStream << "===================== Beta Active =====================" << endl;
    for(unsigned int i = 0; i < beta.size(); ++i) {
        if(beta[i] != 0) {
            outStream << i << ":" << beta[i] << " | ";
        }
    }
    outStream << endl << "======================== End Beta Active ================" << endl;
}

void FusedLassoCoordinate::printDerivActive(ostream& outStream) {
    outStream << "===================== Deriv Active =====================" << endl;
    for(unsigned int i = 0; i < beta.size(); ++i) {
        if(beta[i] != 0) {
            outStream << i << ":" << quadratic->getDerivative(i) << " | ";
        }
    }
    outStream << endl << "======================== End Deriv Active ================" << endl;
}

void FusedLassoCoordinate::printBeta(ostream& outStream) {
    outStream << "======================= Beta ===============" << endl;
    for(unsigned int i = 0; i < beta.size(); ++i) {
        outStream << beta[i] << ", ";
    }
    outStream << endl << "================ End Beta ================" << endl;
}

void FusedLassoCoordinate::printDerivs(ostream& outStream) {
    outStream << "============= Derivs =============" << endl;
    for(unsigned int i = 0; i < beta.size(); ++i) {
        outStream << quadratic->getDerivative(i) << " " ;
    }
    outStream << endl << "============= End Derivs =========" << endl;
}

void FusedLassoCoordinate::printConnectionsWeights(ostream& outStream) {
    outStream << "============ Conn Weights ============" << endl;
    for(unsigned int i = 0; i < connections.size(); ++i) {
        outStream << "Node:" << i << endl; 
        for(unsigned int j = 0; j < connections[i].size(); ++j) {
            outStream << connections[i][j] << " " << wLambda2Mult[i][j] << ";";
        }
        outStream << endl;
    }

}

void FusedLassoCoordinate::printPosSingleStepInfo(int pos, ostream& outStream) {
    outStream << "Info at Step: " << pos << endl;
    outStream << "Deriv: " << quadratic->getDerivative(pos) << endl;    
    outStream << "Hessian: " << quadratic->getHessian(pos) << endl;
    outStream << "OldBeta: " << beta[pos] << endl;
    outStream << "wLambda1: " << wLambda1Mult[pos] << endl;
    outStream << "Connections: " << endl;
    printVector(connections[pos], outStream);
    outStream << "wLambda2Mult: " << endl;
    printVector(wLambda2Mult[pos], outStream);
}

