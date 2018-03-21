#include "RDM_class.hpp"
namespace rdm {

const Eigen::MatrixXd &RDM_class::getOneRDMaa() const {
    return oneRDMaa;
}

const Eigen::MatrixXd &RDM_class::getOneRDMbb() const {
    return oneRDMbb;
}

const Eigen::Tensor<double, 4> &RDM_class::getTwoRDMaaaa() const {
    return twoRDMaaaa;
}

const Eigen::Tensor<double, 4> &RDM_class::getTwoRDMabba() const {
    return twoRDMabba;
}

const Eigen::Tensor<double, 4> &RDM_class::getTwoRDMbbbb() const {
    return twoRDMbbbb;
}

RDM_class::RDM_class(std::vector<double> coefs) {
    this->coefs = Eigen::Map<Eigen::VectorXd>(coefs.data(), coefs.size());
}

RDM_class::RDM_class(Eigen::VectorXd coefs) {
    this->coefs = std::move(coefs);
}

RDM_class::RDM_class(double *coefs, size_t length) {
    this->coefs = Eigen::Map<Eigen::VectorXd>(coefs, length);

}

const Eigen::VectorXd &RDM_class::getCoefs() const {
    return coefs;
}

RDM_class::RDM_class() {}

const Eigen::Tensor<double, 4> &RDM_class::getTwoRDMbaab() const {
    return twoRDMbaab;
}
}