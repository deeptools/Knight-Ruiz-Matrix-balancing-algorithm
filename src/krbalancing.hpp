/*Conversion of the Knight-Ruiz matrix balancing algorithm from matlab to c++*/
#ifndef KRBALANCING_HPP
#define KRBALANCING_HPP

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>
#include <pybind11/numpy.h>
#include <pybind11/eigen.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include <set>
#include <map>
#include <math.h>
#include <algorithm>
#include <omp.h>
#include <limits>
namespace py = pybind11;

typedef Eigen::Matrix<int64_t, Eigen::Dynamic, 1> VectorXint64;

using SparseMatrixCol = Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t>;
size_t num_threads = 10; //TODO add it as an argument to be set by user

class kr_balancing
{
public:
  kr_balancing(const int64_t &input_rows, const int64_t &input_cols,
               const int64_t &input_nnz,
               const Eigen::Ref<VectorXint64> input_outer,
               const Eigen::Ref<VectorXint64> input_inner,
               const Eigen::Ref<Eigen::VectorXd> input_values);
  ~kr_balancing() {}
  void computeKR();
  void outer_loop();
  void inner_loop();
  void compute_normalised_matrix(bool &);
  void rescale_norm_vector();
  const SparseMatrixCol *get_normalised_matrix(bool &rescale);
  const SparseMatrixCol *get_normalisation_vector(bool &rescale);

private:
  std::vector<double> res;
  int fl = 0; //0 = on , 1 = off
  int Delta = 3;
  double delta = 0.1;
  double tol = 1e-6;
  double g = 0.9;
  double etamax = 0.1;
  Eigen::MatrixXd e;
  SparseMatrixCol A;
  SparseMatrixCol rho_km1;
  SparseMatrixCol rho_km2;
  int k;
  Eigen::VectorXd y;
  SparseMatrixCol p;
  SparseMatrixCol Z;
  double innertol;
  int i = 0; //Outer itteration count
  int MVP = 0;
  SparseMatrixCol v;
  SparseMatrixCol x;
  Eigen::SparseVector<double> rk;
  SparseMatrixCol output;
  bool rescaled;

};

#endif
