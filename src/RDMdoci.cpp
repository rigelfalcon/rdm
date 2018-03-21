#include "RDMdoci.hpp"
namespace rdm{
void RDMdoci::initialize(size_t K, size_t npairs) {

    if (K < npairs) {
        throw std::overflow_error("Invalid argument: too many electrons to place into the given number of spatial orbitals");
    }

    this->npairs = npairs;
    this->K = K;
    // Set the number of basis functions
    auto nbf_ = boost::math::binomial_coefficient<double>(this->K, this->npairs);
    if (nbf_ > 4294967295.0) {
        // before casting into unsigned long, we have to make sure that it can fit
        throw std::overflow_error("The number of basis functions for the sector is too high to be cast into unsigned long.");
    }
    this->nbf = static_cast<unsigned long>(nbf_);
    this->ad_mat = bmqc::AddressingScheme(this->K, this->npairs);

}

RDMdoci::RDMdoci(std::vector<double> coefs, size_t K, size_t npairs):RDM_class(coefs) {
    initialize(K, npairs);
}

RDMdoci::RDMdoci(double* coefs, size_t length, size_t K, size_t npairs):RDM_class(coefs,length) {
    initialize(K, npairs);
}

RDMdoci::RDMdoci(Eigen::VectorXd coefs, size_t K, size_t npairs):RDM_class(coefs) {
    initialize(K, npairs);
}

void RDMdoci::compute1RDM() {
    this->oneRDMaa = Eigen::MatrixXd::Zero(this->K,this->K);
    size_t bf_base = this->ad_mat.generateBitVector_long(0); //first basis function
    for (size_t i = 0; i < nbf; i++) {
        if(i>0){
            bf_base = bmqc::next_long_permutation(bf_base);
        }

        for (size_t j = 0; j < this->K; j++) {
            if (bf_base & (1<<j)) {
                oneRDMaa(j,j) += coefs(i)*coefs(i);
            }

        }
    }
    oneRDMbb = oneRDMaa;


}

void RDMdoci::compute2RDM() {
    this->twoRDMaaaa = Eigen::Tensor<double, 4>(this->K,this->K,this->K,this->K);
    this->twoRDMaaaa.setZero();
    this->twoRDMabba = Eigen::Tensor<double, 4>(this->K,this->K,this->K,this->K);
    this->twoRDMabba.setZero();


    size_t bf_base = this->ad_mat.generateBitVector_long(0); //first basis function
    for (size_t i = 0; i < nbf; i++) {
        if(i>0){
            bf_base = bmqc::next_long_permutation(bf_base);
        }
        for (size_t j = 0; j < this->K; j++) {
            if (bmqc::annihilation(bf_base,j)){
                twoRDMabba(j,j,j,j) += coefs(i)*coefs(i);

                for(size_t l = 0; l < j; l++){
                    if (bmqc::creation(bf_base,l)){
                        size_t address = this->ad_mat.fetchAddress(bf_base);
                        double coef = coefs(i)*coefs(address);
                        twoRDMabba(l,l,j,j) += coef;
                        twoRDMabba(j,j,l,l) += coef;
                        bf_base -= 1<<l;//flip back (so we don't need to copy the set)
                    }else { //if we can't create we can annihilated (but only on the diagonal, no transform required
                        double coef = coefs(i) * coefs(i);
                        twoRDMaaaa(j, l, j, l) += coef;
                        twoRDMaaaa(j, l, l, j) -= coef;
                        twoRDMaaaa(l, j, j, l) -= coef;
                        twoRDMaaaa(l, j, l, j) += coef;

                        twoRDMabba(j, l, j, l) += coef;
                        twoRDMabba(l, j, l, j) += coef;

                    }

                }
                bf_base += 1<<j;
            }


        }
    }
    this->twoRDMbbbb=twoRDMaaaa;
    this->twoRDMbaab=twoRDMbaab;



}

void RDMdoci::compute2RDMchemical() {
    this->twoRDMaaaa = Eigen::Tensor<double, 4>(this->K,this->K,this->K,this->K);
    this->twoRDMaaaa.setZero();
    this->twoRDMabba = Eigen::Tensor<double, 4>(this->K,this->K,this->K,this->K);
    this->twoRDMabba.setZero();


    size_t bf_base = this->ad_mat.generateBitVector_long(0); //first basis function
    for (size_t i = 0; i < nbf; i++) {
        if(i>0){
            bf_base = bmqc::next_long_permutation(bf_base);
        }
        for (size_t j = 0; j < this->K; j++) {
            if (bmqc::annihilation(bf_base,j)){
                twoRDMabba(j,j,j,j) += coefs(i)*coefs(i);

                for(size_t l = 0; l < j; l++){
                    if (bmqc::creation(bf_base,l)){
                        size_t address = this->ad_mat.fetchAddress(bf_base);
                        double coef = coefs(i)*coefs(address);
                        twoRDMabba(l,j,l,j) += coef;
                        twoRDMabba(j,l,j,l) += coef;
                        bf_base -= 1<<l;//flip back (so we don't need to copy the set)
                    }else { //if we can't create we can annihilated (but only on the diagonal, no transform required
                        double coef = coefs(i) * coefs(i);
                        twoRDMaaaa(l, l, j, j) += coef;
                        twoRDMaaaa(j, l, l, j) -= coef;
                        twoRDMaaaa(l, j, j, l) -= coef;
                        twoRDMaaaa(j, j, l, l) += coef;

                        twoRDMabba(l, l, j, j) += coef;
                        twoRDMabba(j, j, l, l) += coef;

                    }

                }
                bf_base += 1<<j;
            }


        }
    }
    this->twoRDMbbbb=twoRDMaaaa;
    this->twoRDMbaab=twoRDMbaab;



}
}