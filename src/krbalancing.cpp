#include "krbalancing.hpp"

kr_balancing::kr_balancing(const  SparseMatrixCol & input){
    A = input;
    e.resize(A.rows(),1);
    e.setOnes();
    /*Replace zeros with 0.00001 on the main diagonal*/
    SparseMatrixCol I;
    I.resize(A.rows(),A.cols());
    I.setIdentity();
    I = I*0.00001;
    A = A + I;
    rescaled = false;
}


void kr_balancing::computeKR(){
  x = e.sparseView();
  assert(x.isVector());
  v = x.cwiseProduct(A*x);
  rk = Eigen::VectorXd::Constant(v.rows(),1) - v; //Vecotr of 1 - v
  assert(rk.isVector());
  rho_km1 = rk.conjugate().transpose()*rk;
  assert(rho_km1.coeff(0,0)>=0);
  rho_km2 = rho_km1;
  outer_loop();
}


void kr_balancing::outer_loop(){
  size_t outer_loop_count = 0;
  double stop_tol = tol*0.5;
  double eta = etamax;
  double rt = std::pow(tol,2);
  double rout = rho_km1.coeff(0,0);
  double rold = rout;
  if(fl == 1) std::cout<< 'intermediate convergence statistics is off'<<std::endl;
  size_t nochange = 0;
  while(rout > rt){//outer itteration
    outer_loop_count ++;
    i=i+1; k=0; y=e.sparseView();
    innertol = std::max(std::pow(eta,2)*rout,rt);
    inner_loop();
    //assert(rho_km1.coeff(0,0) > 0);
    assert(rho_km1.coeff(0,0) != std::numeric_limits<double>::infinity());
    x=x.cwiseProduct(y);
    assert(x.isVector());
    v=x.cwiseProduct((A*x));
    rk = Eigen::VectorXd::Constant(v.rows(),1) - v;
    rho_km1 = rk.conjugate().transpose()*rk;
    //assert(rho_km1.coeff(0,0) > 0);
    assert(rho_km1.coeff(0,0) != std::numeric_limits<double>::infinity());
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
    if(outer_loop_count %50 ==0) std::cout << "outer loop number "<< outer_loop_count <<std::endl;
  }//End of outer loop
  //    if (nochange >= 100) {
  //       std::cout << "nochange >=100" <<std::endl;
  //       return 0;
  //    }else{
  //    return &x;
  //    }
}


void kr_balancing::inner_loop(){
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
    SparseMatrixCol w=x.cwiseProduct(A*(x.cwiseProduct(p))) + v.cwiseProduct(p);
    SparseMatrixCol alpha_mat = rho_km1.cwiseQuotient(p.conjugate().transpose()*w);
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
      Eigen::VectorXd y_copy;
      y_copy.resize(ind.size());
      for(size_t i = 0; i < ind.size(); i++){
        y_copy(i) = (delta - y(ind(i)))/ap(ind(i));
      }
      double gamma = y_copy.minCoeff();
      y = y + gamma*ap;
      break;
    }
    if(ynew.maxCoeff() >= Delta){
      Eigen::Matrix<bool, Eigen::Dynamic , Eigen::Dynamic> ind_helper = (ynew.array() > Delta);
      Eigen::VectorXi ind = Eigen::VectorXi::LinSpaced(ind_helper.size(),0,ind_helper.size()-1);
      ind.conservativeResize(std::stable_partition(
      ind.data(), ind.data()+ind.size(), [&ind_helper](int i){return ind_helper(i);})-ind.data());
      Eigen::VectorXd y_copy;
      y_copy.resize(ind.size());
      for(size_t i = 0; i < ind.size(); i++){
        y_copy(i) = (Delta - y(ind(i)))/ap(ind(i));
      }
      double gamma = y_copy.minCoeff();
      y = y + (gamma*ap);
      break;
    }
    y = ynew;
    rk = rk - (alpha*w); rho_km2 = rho_km1;
    Z = rk.cwiseQuotient(v); rho_km1 = rk.conjugate().transpose()*Z;
  }//End of the inner 'while'
}


void kr_balancing::compute_normalised_matrix(bool & rescale){
  assert(A.rows() == A.cols());
  if(rescale ==true && rescaled == false){
    rescale_norm_vector();
    rescaled = true;
  }else{
    A = SparseMatrixCol(A.triangularView<Eigen::Upper>());
    std::cout << "here!"<<std::endl;
  }
  #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
  for (int k=0; k<A.outerSize(); ++k){
    #pragma omp critical
    for(SparseMatrixCol::InnerIterator it(A,k); it; ++it){
      it.valueRef() = it.value()*x.coeff(it.row(),0)*x.coeff(it.col(),0);
    }
  }
}


//Rescaling the normalisation factor (x)
void kr_balancing::rescale_norm_vector(){
  double original_sum = 0.0;
  double norm_vector_sum = 0.0;
  assert(A.rows() == A.cols());
  A = SparseMatrixCol(A.triangularView<Eigen::Upper>());

  #pragma omp parallel for num_threads(num_threads)
  for (int k=0; k<A.outerSize(); ++k){
    #pragma omp critical
    for(SparseMatrixCol::InnerIterator it(A,k); it; ++it){
      if(it.row()==it.col()){
        original_sum += it.value();
        norm_vector_sum += it.value()*x.coeff(it.row(),0)*x.coeff(it.col(),0);
      }else{
        original_sum += it.value()*2;
        norm_vector_sum += it.value()*x.coeff(it.row(),0)*x.coeff(it.col(),0)*2;
      }
    }
  }

  std::cout << "normalisation factor is "<< std::sqrt(norm_vector_sum/original_sum) <<std::endl;
  x /= std::sqrt(norm_vector_sum/original_sum);

}


const SparseMatrixCol* kr_balancing::get_normalised_matrix(bool & rescale){
  compute_normalised_matrix(rescale);
  return &A;
}


const SparseMatrixCol* kr_balancing::get_normalisation_vector(bool & rescale){
  if(rescale ==true && rescaled == false){
    rescale_norm_vector();
    rescaled = true;
  }
  return &x;
}


PYBIND11_MODULE(krbalancing, m) {
  py::class_<kr_balancing>(m, "kr_balancing")
    .def(py::init< const SparseMatrixCol & >())
    .def("computeKR", &kr_balancing::computeKR)
    .def("get_normalisation_vector",&kr_balancing::get_normalisation_vector, py::return_value_policy::reference_internal)
    .def("get_normalised_matrix",&kr_balancing::get_normalised_matrix, py::return_value_policy::reference_internal);

}

//c++ -O3 -Wall -I /path/to/eigen -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes) KRBalancing.cpp -o KRBalancing$(python3-config --extension-suffix)
