/*Conversion of the Knight-Ruiz matrix balancing algorithm from matlab to c++*/
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
//#include <Eigen/Sparse>
//#include <Eigen/Dense>
#include <limits>
namespace py = pybind11;
using SparseMatrixR = Eigen::SparseMatrix<double,Eigen::ColMajor>;


class kr_balancing{
     public:

      kr_balancing(const SparseMatrixR & input){
        std::cout<< "read input"<<std::endl;
        A = input;
        e.resize(A.rows(),1);
        e.setOnes();
        /*Replace zeros with 0.00001 on the main diagonal*/
        SparseMatrixR I;
        I.resize(A.rows(),A.cols());
        I.setIdentity();
        I = I*0.00001;
        A = A + I;
        x0 = e.sparseView();
      }
      size_t outter_loop_count = 0;
      void outter_loop(){
        double stop_tol = tol*0.5;
        double eta = etamax;
        x = x0;
        double rt = std::pow(tol,2);
        v = x.cwiseProduct(A*x);
        rk = Eigen::VectorXd::Constant(v.rows(),1) - v; //Vecotr of 1 - v
        rho_km1 = rk.conjugate().transpose()*rk;
        rho_km2 = rho_km1;
        double rout = rho_km1.coeff(0,0);
        double rold = rout;
        if(fl == 1) std::cout<< 'intermediate convergence statistics is off'<<std::endl;
        while(rout > rt){//outer itteration
          outter_loop_count ++;
          i=i+1; k=0; y=e.sparseView();
          innertol = std::max(std::pow(eta,2)*rout,rt);
          inner_loop();

          x=x.cwiseProduct(y);

          v=x.cwiseProduct(A*x);
          rk = Eigen::VectorXd::Constant(v.rows(),1) - v;
          rho_km1 = rk.conjugate().transpose()*rk;
          rout = rho_km1.coeff(0,0);
          MVP = MVP + k + 1;
          //Update inner iteration stopping criterion.
          double rat = rout/rold; rold = rout; double res_norm = std::sqrt(rout);
          double eta_o = eta; eta = g*rat;
          if(g*std::pow(eta_o,2) > 0.1){
            eta = std::max(eta,g*std::pow(eta_o,2));
          }
          eta = std::max(std::min(eta,etamax),stop_tol/res_norm);
          if(fl == 1){
            res.push_back(res_norm);
          }
          if(outter_loop_count %50 ==0) std::cout << "outer loop number "<< outter_loop_count <<std::endl;
         }
      }
      void inner_loop(){
        size_t inner_loop_count = 0;
        while (rho_km1.coeff(0,0) > innertol){ // Inner itteration (conjugate gradient method)
          inner_loop_count ++;
          k++;
          if(k == 1){
            Z = rk.cwiseQuotient(v);
            p=Z;
            rho_km1 = rk.conjugate().transpose()*Z;
          }else{
            Eigen::VectorXd beta=rho_km1.cwiseQuotient(rho_km2);
            p = Z + (beta(0)*p);

          }
          //Update search direction efficiently.
          SparseMatrixR w=x.cwiseProduct(A*(x.cwiseProduct(p))) + v.cwiseProduct(p);
          SparseMatrixR alpha_mat = rho_km1.cwiseQuotient(p.conjugate().transpose()*w);
          double alpha = alpha_mat.coeff(0,0);
          Eigen::VectorXd ap = alpha * p;
          //Test distance to boundary of cone.
          Eigen::VectorXd ynew = y + ap;

          if(ynew.minCoeff() <= delta){
            if(delta == 0) break;
            Eigen::Matrix<bool, Eigen::Dynamic , Eigen::Dynamic> ind_helper = (ap.array()<0);
            Eigen::VectorXi ind = Eigen::VectorXi::LinSpaced(ind_helper.size(),0,ind_helper.size()-1);
            ind.conservativeResize(std::stable_partition(
              ind.data(), ind.data()+ind.size(), [&ind_helper](int i){return ind_helper(i);})-ind.data());
            Eigen::VectorXd y_copy = y;
            for(i = 0; i < ind.size(); i++){ //TODO proper masking? //ind_vec.unaryExpr(x);
               y_copy(ind(i)) = (delta - y_copy(ind(i)))/ap(i);
            }
            double gamma = y_copy.minCoeff();
            y = y + gamma*ap;
            break;
          }
          if(ynew.minCoeff() >= Delta){
            Eigen::Matrix<bool, Eigen::Dynamic , Eigen::Dynamic> ind_helper = (ynew.array() > Delta);
            Eigen::VectorXi ind = Eigen::VectorXi::LinSpaced(ind_helper.size(),0,ind_helper.size()-1);
            ind.conservativeResize(std::stable_partition(
              ind.data(), ind.data()+ind.size(), [&ind_helper](int i){return ind_helper(i);})-ind.data());
            Eigen::VectorXd y_copy = y;
            for(i = 0; i < ind.size(); i++){ //TODO proper masking? //ind_vec.unaryExpr(x);
                 y_copy(ind(i)) = (Delta - y_copy(ind(i)))/ap(i);
            }
            double gamma = y_copy.minCoeff();
            y = y + gamma*ap;
            break;
          }

          y = ynew;
          rk = rk - (alpha*w); rho_km2 = rho_km1;
          Z = rk.cwiseQuotient(v); rho_km1 = rk.conjugate().transpose()*Z;
          if(inner_loop_count %100 ==0) std::cout << "inner loop number "<< inner_loop_count <<std::endl;


        }//End of the inner 'while'
      }
      SparseMatrixR  get_output(){
        std::cout << "export output of shape "<< A.cwiseProduct(x*x.transpose()).rows() <<std::endl;
        return A.cwiseProduct(x*x.transpose());
      }
     private:
       std::vector<double> res;
       unsigned int fl = 0; //0 = on , 1 = off
       unsigned int Delta = 3;
       double delta = 0.1;
       double tol = 1e-6;
       double g = 0.9;
       double etamax = 0.1;
       SparseMatrixR x0;
       Eigen::MatrixXd e;
       SparseMatrixR A;
       SparseMatrixR rho_km1;
       SparseMatrixR rho_km2;
       unsigned int k;
       Eigen::VectorXd y;
       SparseMatrixR p;
       SparseMatrixR Z;
       double innertol;
       unsigned int i = 0; //Outer itteration count
       unsigned int MVP = 0;
       SparseMatrixR v;
       SparseMatrixR x;
       Eigen::SparseVector<double> rk;
};


PYBIND11_MODULE(KRBalancing, m) {
  py::class_<kr_balancing>(m, "kr_balancing")
    .def(py::init< const SparseMatrixR & >())
    .def("outter_loop", &kr_balancing::outter_loop)
    .def("get_output",&kr_balancing::get_output, py::return_value_policy::reference_internal);

}

//c++ -O3 -Wall -I /path/to/eigen -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes) KRBalancing.cpp -o KRBalancing$(python3-config --extension-suffix)
