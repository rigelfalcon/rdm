#define BOOST_TEST_MODULE "RDMdoci"

#include "RDMdoci.hpp"

#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>
#define threshold 1.0e-06

BOOST_AUTO_TEST_CASE ( RDMdoci_test ) {

    std::vector<double> test_set = {0.9852125942454278, 0.005366388840493838, 0.014549652019718093, 0.09666112425601708, 0.07827561044365904, 0.11681271122761872};
    rdm::RDMdoci test_rdmdoci(test_set,4,2);
    test_rdmdoci.compute1RDM();
    test_rdmdoci.compute2RDM();

    Eigen::MatrixXd rdm_check = test_rdmdoci.getOneRDMaa();
    BOOST_CHECK(std::abs(rdm_check(0,0) - 0.9800160269314305)< threshold);
    BOOST_CHECK(std::abs(rdm_check(1,1) - 0.9769826194240284)< threshold);
    BOOST_CHECK(std::abs(rdm_check(2,2) - 0.0138857000074293)< threshold);
    BOOST_CHECK(std::abs(rdm_check(3,3) - 0.0291156536371116)< threshold);

    Eigen::Tensor<double, 4> abba = test_rdmdoci.getTwoRDMabba();
    Eigen::Tensor<double, 4> aaaa = test_rdmdoci.getTwoRDMaaaa();

    //FIXME test if 2RDM are correct non-visually


}

