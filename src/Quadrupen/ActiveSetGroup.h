/*
 * Author: Julien CHIQUET
 */
#ifndef _ActiveSetGroup_H
#define _ActiveSetGroup_H

using namespace Rcpp;
using namespace arma;
using namespace std;

#include "RegressionData.h"
#include "ActiveSet.h"

template <typename matrix>
class ActiveSetGroup: public ActiveSet<matrix> {

public:
  
  using ActiveSet<matrix>::A_ ;
  
  // VARIABLES FOR HANDLING THE ACTIVE SET
  uvec G_             ; // set of currently activated groups
  uvec is_grp_in_     ; // indicator of active groups (0/1)
  uvec grp_sizes_     ; // vector of group sizes
  vector<uvec> group_ ; // vector of current parameters
  
public:
  
  ActiveSetGroup() {} ;
  ActiveSetGroup(const RegressionData<matrix> &data, const uvec& group_ind, const bool use_chol) ;

  // ACTIVE SET HANDLING
  void add_group(uword, const RegressionData<matrix> &) ; // add a group of variables in the active set
  void del_group(uword, vec&) ; // remove the group of variable activated in position ind_grp_class
  void del_groups(uvec igrps_out, vec& beta) ;
  const uword size_grp() const { return G_.n_elem ; }

};

template <typename matrix>
ActiveSetGroup<matrix>::ActiveSetGroup(const RegressionData<matrix>& data, const uvec& group_ind, const bool use_chol) :
  ActiveSet<matrix>(data, use_chol) {

    // Vector of group and group sizes
    uvec grp = unique(group_ind) ;
    uword nb_grp =  grp.n_elem ;
    grp_sizes_.zeros(nb_grp) ;
    for (auto it = group_ind.begin(); it != group_ind.end(); it++) {
      grp_sizes_[*it - 1]++;
    }
    uword current_elt = 0 ;
    for (auto it = grp_sizes_.begin(); it != grp_sizes_.end(); it++) {
      group_.push_back(regspace<uvec>(current_elt, current_elt + *it - 1)) ;
      current_elt = current_elt + *it ;
    }

    is_grp_in_.zeros(grp_sizes_.n_elem) ;
}

template <typename matrix>
void ActiveSetGroup<matrix>::add_group(uword grp_in, const RegressionData<matrix>& data) {
  
  G_.resize(size_grp() + 1) ;
  G_.tail(1) = grp_in       ;
  is_grp_in_[grp_in] = 1    ;

  ActiveSet<matrix>::add_vars(group_[grp_in], data) ;
}

template <typename matrix>
void ActiveSetGroup<matrix>::del_group(uword igrp_out, vec& beta) {
  
  uword ifirst_var_in = 0 ;
  for (uword i=0; i < igrp_out; i++) ifirst_var_in += grp_sizes_(G_[i]) ;
  uword n_vars_to_del = grp_sizes_(G_[igrp_out]);
  uvec ivars_out = regspace<uvec>(ifirst_var_in + n_vars_to_del - 1, ifirst_var_in);
  
  this->del_vars(ivars_out, beta) ;
  
  is_grp_in_[G_[igrp_out]] = 0 ;
  G_.shed_row(igrp_out)        ;
}

template <typename matrix>
void ActiveSetGroup<matrix>::del_groups(uvec igrps_out, vec& beta) {
  if (igrps_out.is_empty()) return;
  
  uvec sorted_indices = sort(igrps_out, "descend");
  for (uword i = 0; i < sorted_indices.n_elem; ++i) {
    uword idx = sorted_indices(i);
    if (idx < G_.n_elem) {
      this->del_group(idx, beta);
    }
  }
}

#endif

