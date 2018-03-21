#ifndef RDM_RDM_HPP
#define RDM_RDM_HPP
#include <bitlong.hpp>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>
#include <boost/math/special_functions.hpp>
#include "bmqc.hpp"
namespace rdm{
class RDM_class {
protected:
    bool orc = false; //onerdms computed?
    bool trc = false; //twordms computed?
    double precision = 16; //double precision
    Eigen::VectorXd coefs; //CI vector
    Eigen::MatrixXd oneRDMaa; //reduced density matrix of one electron operators (alpha alpha)
    Eigen::MatrixXd oneRDMbb; //reduced density matrix of one electron operators (beta beta)
    Eigen::Tensor<double, 4> twoRDMaaaa; //reduced density matrix two electron operators
    Eigen::Tensor<double, 4> twoRDMabba; //reduced density matrix two electron operators
    Eigen::Tensor<double, 4> twoRDMbaab; //reduced density matrix two electron operators
    Eigen::Tensor<double, 4> twoRDMbbbb; //reduced density matrix two electron operators
//constructorss
    RDM_class(std::vector<double> coefs);
    RDM_class(Eigen::VectorXd coefs);
    RDM_class(double* coefs, size_t length);
public:
    RDM_class();
    virtual void compute1RDM()=0;
    virtual void compute2RDM()=0;
public:
    //Getters
    const Eigen::VectorXd &getCoefs() const;

    const Eigen::MatrixXd &getOneRDMaa() const;

    const Eigen::MatrixXd &getOneRDMbb() const;

    const Eigen::Tensor<double, 4> &getTwoRDMaaaa() const;

    const Eigen::Tensor<double, 4> &getTwoRDMabba() const;

    const Eigen::Tensor<double, 4> &getTwoRDMbaab() const;

    const Eigen::Tensor<double, 4> &getTwoRDMbbbb() const;
};
}

#endif //RDM_RDM_HPP
