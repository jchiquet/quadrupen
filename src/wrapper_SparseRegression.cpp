/*
 * Author: Julien CHIQUET
 *         MIA Paris-Saclay
 */

#include "Quadrupen/SparseRegularizer.h"

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List elastic_net_dense_cpp(
    const Environment &dataModel,
    const bool        &intercept,
    const List        &regParam,
    const List        &control
    ) {
  RegressionData<mat> data(dataModel, intercept, as<bool>(control["normalize"])) ;
  SparseRegularizer<mat,SparseNorm::L1> enet(std::move(data), regParam, control) ;
  return enet.to_list(enet.solution_path(control)) ;
}

// [[Rcpp::export]]
List elastic_net_sparse_cpp(
    const Environment &dataModel,
    const bool        &intercept,
    const List        &regParam,
    const List        &control
) {
  RegressionData<sp_mat> data(dataModel, intercept, as<bool>(control["normalize"])) ;
  SparseRegularizer<sp_mat,SparseNorm::L1> enet(std::move(data), regParam, control) ;
  return enet.to_list(enet.solution_path(control)) ;
}

// [[Rcpp::export]]
List mcp_dense_cpp(
    const Environment &dataModel,
    const bool        &intercept,
    const List        &regParam,
    const List        &control
) {
  RegressionData<mat> data(dataModel, intercept, as<bool>(control["normalize"])) ;
  SparseRegularizer<mat,SparseNorm::MCP> enet(std::move(data), regParam, control) ;
  return enet.to_list(enet.solution_path(control)) ;
}

// [[Rcpp::export]]
List mcp_sparse_cpp(
    const Environment &dataModel,
    const bool        &intercept,
    const List        &regParam,
    const List        &control
) {
  RegressionData<sp_mat> data(dataModel, intercept, as<bool>(control["normalize"])) ;
  SparseRegularizer<sp_mat,SparseNorm::MCP> enet(std::move(data), regParam, control) ;
  return enet.to_list(enet.solution_path(control)) ;
}

// [[Rcpp::export]]
List scad_dense_cpp(
    const Environment &dataModel,
    const bool        &intercept,
    const List        &regParam,
    const List        &control
) {
  RegressionData<mat> data(dataModel, intercept, as<bool>(control["normalize"])) ;
  SparseRegularizer<mat,SparseNorm::SCAD> enet(std::move(data), regParam, control) ;
  return enet.to_list(enet.solution_path(control)) ;
}

// [[Rcpp::export]]
List scad_sparse_cpp(
    const Environment &dataModel,
    const bool        &intercept,
    const List        &regParam,
    const List        &control
) {
  RegressionData<sp_mat> data(dataModel, intercept, as<bool>(control["normalize"])) ;
  SparseRegularizer<sp_mat,SparseNorm::SCAD> enet(std::move(data), regParam, control) ;
  return enet.to_list(enet.solution_path(control)) ;
}
