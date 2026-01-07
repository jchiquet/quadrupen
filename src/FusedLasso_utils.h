#ifndef _FUSEDLASSO_UTILS_H_
#define _FUSEDLASSO_UTILS_H_

#include <Rcpp.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <utility>
#include <vector>
#include <limits>
#include <string>
#include <list>
#include <set>
#include <map>
#include <queue>
#include <math.h>

using namespace std;

#define Abs(x)    ((x) < 0 ? -(x) : (x))
#define Max(a, b) ((a) > (b) ? (a) : (b))
#define Min(a, b) ((a) < (b) ? (a) : (b))

const double tolerance =1.0e-8;
const double infinite = std::numeric_limits<double>::max();
const int infiniteInt = std::numeric_limits<int>::max();

double RelDif(double a, double b);
double RelDifNoAbs(double a, double b);

double maxDiffDoubleVec(const vector<double>&x, const vector<double>& y);

inline int signum(double x) {return((x>0)-(x<0));};

int numNonZero(const vector<double>& x);

// read in a matrix from file filename
void readDoubleMatrix(vector<double>& X, int& n, int& p, string filename);

// helper function for readDoubleMatrix -> transposes the input
// needed as in R matrices are in column major order
vector<double> transpose(vector<double>& X, int n, int p);

// checks if all positions in x belongs to the same group in fusions
bool checkIfSubset(const set<int>& x, const vector<int>& fusions);

void printVector(const vector<double>& x, ostream& outStream); 

void printVector(const vector<int>& x, ostream& outStream); 

void printList(const list<double>& x, ostream& outStream); 

void printList(const list<int>& x, ostream& outStream); 

void printMatrix(vector<double>& X, int n, int p, ostream& outStream);

/* CLASS STEPS
 * Given the current betas, it finds the root of the function
 */

class Steps {
  inline static bool updateBeta(vector<double>& beta, const int pos, double& deriv, const double& slope, const double& zeroStepSize, const vector<int>& stepIdx, const vector<double>& stepSize);
  
public:
  // most important application -> find the root
  
  static void findRoot(vector<double> &beta, double derivQuad, const double slope, const int pos, const double zeroStepSize, const vector<int>& stepIdx, const vector<double>& stepSize);
  
  static void findRootL1(vector<double> &beta, double derivQuad, const double hessian, const int pos, const double zeroStepSize, const vector<int>& stepIdx, const vector<double>& stepSize);
  
  static void findRootL2(vector<double> &beta, double derivQuad, const double hessian, const int pos, const double zeroStepSize, const vector<int>& l2Idx, const vector<double>& l2Mult);
  
  static void findRootHuber(vector<double> &beta, double derivQuad, const double hessian, const int pos, const double zeroStepSize, const vector<int>& huberIdx, const vector<double>& huberMult, double huberParam);
  
};

/*
 * counted_ptr - simple reference counted pointer.
 *
 * The is a non-intrusive implementation that allocates an additional
 * int and pointer for every counted object.
 */

template <class X> class counted_ptr
{
public:
  typedef X element_type;
  
  explicit counted_ptr(X* p = 0) // allocate a new counter
    : itsCounter(0) {if (p) itsCounter = new counter(p);}
  ~counted_ptr()
  {release();}
  counted_ptr(const counted_ptr& r) throw()
  {acquire(r.itsCounter);}
  counted_ptr& operator=(const counted_ptr& r)
  {
    if (this != &r) {
      release();
      acquire(r.itsCounter);
    }
    return *this;
  }
  
  X& operator*()  const throw()   {return *itsCounter->ptr;}
  X* operator->() const throw()   {return itsCounter->ptr;}
  X* get()        const throw()   {return itsCounter ? itsCounter->ptr : 0;}
  bool unique()   const throw()
  {return (itsCounter ? itsCounter->count == 1 : true);}
  
private:
  
  struct counter {
    counter(X* p = 0, unsigned c = 1) : ptr(p), count(c) {}
    X*          ptr;
    unsigned    count;
  }* itsCounter;
  
  void acquire(counter* c) throw()
  { // increment the count
    itsCounter = c;
    if (c) ++c->count;
  }
  
  void release()
  { // decrement the count, delete if it is 0
    if (itsCounter) {
      if (--itsCounter->count == 0) {
        delete itsCounter->ptr;
        delete itsCounter;
      }
      itsCounter = 0;
    }
  }
};


#endif // _FUSEDLASSO_UTILS_H_
