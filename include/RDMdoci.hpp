#ifndef RDM_RDMDOCI_HPP
#define RDM_RDMDOCI_HPP

#include "RDM_class.hpp"
namespace rdm{
class RDMdoci: public RDM_class {
private:
    size_t K;
    size_t npairs;
    size_t nbf;
    bmqc::AddressingScheme ad_mat;
    void initialize(size_t K, size_t npairs); //helper constructor;


public:
    /**
     * Constructors
     * for different type of representations of coefficients vectors.
     */
    RDMdoci(std::vector<double> coefs,size_t K, size_t npairs); //std::vector<double>
    RDMdoci(double* coefs, size_t length, size_t K, size_t npairs); //double[]
    RDMdoci(Eigen::VectorXd coefs,size_t K, size_t npairs); //Eigen::VectorXd

    void compute1RDM() override;
    void compute2RDM() override;
    void compute2RDMchemical();
};
}

#endif //RDM_RDMDOCI_HPP
