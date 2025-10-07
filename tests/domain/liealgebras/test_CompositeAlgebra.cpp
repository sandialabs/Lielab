#include <Lielab.hpp>

#include <catch2/catch_all.hpp>

#include <iostream>

const Lielab::domain::cn ycn1 = Lielab::domain::cn::from_vector({1.0, 2.0, 3.0, 4.0});
const Lielab::domain::glc yglc1 = Lielab::domain::glc::from_vector({5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0});
const Lielab::domain::glr yglr1 = Lielab::domain::glr::from_vector({13.0, 14.0, 15.0, 16.0});
const Lielab::domain::rn yrn1 = Lielab::domain::rn::from_vector({17.0, 18.0});
const Lielab::domain::se yse1 = Lielab::domain::se::from_vector({19.0, 20.0, 21.0});
const Lielab::domain::so yso1 = Lielab::domain::so::from_vector({22.0, 23.0, 24.0});
const Lielab::domain::sp ysp1 = Lielab::domain::sp::from_vector({25.0, 26.0, 27.0});
const Lielab::domain::su ysu1 = Lielab::domain::su::from_vector({28.0, 29.0, 30.0});

TEST_CASE("CompositeAlgebra to_string", "[domain]")
{
    using namespace Lielab::domain;

    const CompositeAlgebra y1 = CompositeAlgebra({ycn1, yglc1, yglr1, yrn1, yse1, yso1, ysp1, ysu1});
    std::vector<CompositeAlgebra::TYPES> elements;
    elements.push_back(ycn1);
    elements.push_back(yglc1);
    elements.push_back(yglr1);
    elements.push_back(yrn1);
    elements.push_back(yse1);
    elements.push_back(yso1);
    elements.push_back(ysp1);
    elements.push_back(ysu1);
    const CompositeAlgebra y2 = CompositeAlgebra(elements);

    CHECK(y1.to_string() == "c^2 ⊕ gl(2, C) ⊕ gl(2, R) ⊕ r^2 ⊕ se(2) ⊕ so(3) ⊕ sp(2, R) ⊕ su(2)");
    CHECK(y2.to_string() == "c^2 ⊕ gl(2, C) ⊕ gl(2, R) ⊕ r^2 ⊕ se(2) ⊕ so(3) ⊕ sp(2, R) ⊕ su(2)");
}

TEST_CASE("CompositeAlgebra main_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const CompositeAlgebra xblank = CompositeAlgebra();
    CHECK(xblank.get_dimension() == 0);

    const CompositeAlgebra x0 = CompositeAlgebra(0);
    CHECK(x0.get_dimension() == 0);
    const CompositeAlgebra x1 = CompositeAlgebra(1);
    CHECK(x1.get_dimension() == 2);
    const CompositeAlgebra x10 = CompositeAlgebra(10);
    CHECK(x10.get_dimension() == 200);
}

TEST_CASE("CompositeAlgebra list_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const CompositeAlgebra y1 = CompositeAlgebra({ycn1, yglc1, yglr1, yrn1, yse1, yso1, ysp1, ysu1});
    std::vector<CompositeAlgebra::TYPES> elements;
    elements.push_back(ycn1);
    elements.push_back(yglc1);
    elements.push_back(yglr1);
    elements.push_back(yrn1);
    elements.push_back(yse1);
    elements.push_back(yso1);
    elements.push_back(ysp1);
    elements.push_back(ysu1);
    const CompositeAlgebra y2 = CompositeAlgebra(elements);
}

TEST_CASE("CompositeAlgebra basis_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const CompositeAlgebra xm10 = CompositeAlgebra::basis(-1, 0);
    CHECK(xm10.get_dimension() == 0);
    const Eigen::VectorXd xm10bar = xm10.get_vector();
    CHECK(xm10bar.size() == 0);

    const CompositeAlgebra x00 = CompositeAlgebra::basis(0, 0);
    CHECK(x00.get_dimension() == 0);
    const Eigen::VectorXd x00bar = x00.get_vector();
    CHECK(x00bar.size() == 0);

    const CompositeAlgebra x10 = CompositeAlgebra::basis(1, 0);
    CHECK(x10.get_dimension() == 0);
    const Eigen::VectorXd x10bar = x10.get_vector();
    CHECK(x10bar.size() == 0);

    const CompositeAlgebra x01 = CompositeAlgebra::basis(0, 1);
    CHECK(x01.get_dimension() == 2);
    const Eigen::VectorXd x01bar = x01.get_vector();
    REQUIRE(x01bar.size() == 2);
    CHECK(x01bar(0) == 1.0);
    CHECK(x01bar(1) == 0.0);

    const CompositeAlgebra x11 = CompositeAlgebra::basis(1, 1);
    CHECK(x11.get_dimension() == 2);
    const Eigen::VectorXd x11bar = x11.get_vector();
    REQUIRE(x11bar.size() == 2);
    CHECK(x11bar(0) == 0.0);
    CHECK(x11bar(1) == 1.0);

    const CompositeAlgebra x21 = CompositeAlgebra::basis(2, 1);
    CHECK(x21.get_dimension() == 2);
    const Eigen::VectorXd x21bar = x21.get_vector();
    REQUIRE(x21bar.size() == 2);
    CHECK(x21bar(0) == 0.0);
    CHECK(x21bar(1) == 0.0);

    const CompositeAlgebra x02 = CompositeAlgebra::basis(0, 2);
    CHECK(x02.get_dimension() == 8);
    const Eigen::VectorXd x02bar = x02.get_vector();
    REQUIRE(x02bar.size() == 8);
    CHECK(x02bar(0) == 1.0);
    CHECK(x02bar(1) == 0.0);
    CHECK(x02bar(2) == 0.0);
    CHECK(x02bar(3) == 0.0);
    CHECK(x02bar(4) == 0.0);
    CHECK(x02bar(5) == 0.0);
    CHECK(x02bar(6) == 0.0);
    CHECK(x02bar(7) == 0.0);
}

TEST_CASE("CompositeAlgebra from_shape_initializer", "[domain]")
{
    using namespace Lielab::domain;

    const CompositeAlgebra x0 = CompositeAlgebra::from_shape(0);
    CHECK(x0.get_dimension() == 0);
    const Eigen::MatrixXcd x0hat = x0.get_matrix();
    CHECK(x0hat.rows() == 0);
    CHECK(x0hat.cols() == 0);

    const CompositeAlgebra x1 = CompositeAlgebra::from_shape(1);
    CHECK(x1.get_dimension() == 2);
    const Eigen::MatrixXcd x1hat = x1.get_matrix();
    CHECK(x1hat.rows() == 1);
    CHECK(x1hat.cols() == 1);

    const CompositeAlgebra x2 = CompositeAlgebra::from_shape(2);
    CHECK(x2.get_dimension() == 8);
    const Eigen::MatrixXcd x2hat = x2.get_matrix();
    CHECK(x2hat.rows() == 2);
    CHECK(x2hat.cols() == 2);
}

TEST_CASE("CompositeAlgebra get_dimension", "[domain]")
{
    using namespace Lielab::domain;

    CompositeAlgebra zero(0), one(1), two(2), three(3), four(4), five(5), six(6), seven(7), eight(8);

    // Dimensions
    CHECK(zero.get_dimension() == 0);
    CHECK(one.get_dimension() == 2);
    CHECK(two.get_dimension() == 8);
    CHECK(three.get_dimension() == 18);
    CHECK(four.get_dimension() == 32);
    CHECK(five.get_dimension() == 50);
    CHECK(six.get_dimension() == 72);
    CHECK(seven.get_dimension() == 98);
    CHECK(eight.get_dimension() == 128);
}

TEST_CASE("CompositeAlgebra set/get_vector", "[domain]")
{
    /*!
    * Tests the set/get_vector operation.
    */

    using namespace Lielab::domain;

    CompositeAlgebra y1 = CompositeAlgebra({ycn1, yglc1, yglr1, yrn1, yse1, yso1, ysp1, ysu1});
    y1.set_vector({30.0, 29.0, 28.0, 27.0, 26.0, 25.0, 24.0, 23.0, 22.0, 21.0,
                   20.0, 19.0, 18.0, 17.0, 16.0, 15.0, 14.0, 13.0, 12.0, 11.0,
                   10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0});
    
    const Eigen::VectorXd y1bar = y1.get_vector();
    REQUIRE(y1bar.size() == 30);
    CHECK(y1bar(0) == 30.0);
    CHECK(y1bar(1) == 29.0);
    CHECK(y1bar(2) == 28.0);
    CHECK(y1bar(3) == 27.0);
    CHECK(y1bar(4) == 26.0);
    CHECK(y1bar(5) == 25.0);
    CHECK(y1bar(6) == 24.0);
    CHECK(y1bar(7) == 23.0);
    CHECK(y1bar(8) == 22.0);
    CHECK(y1bar(9) == 21.0);
    CHECK(y1bar(10) == 20.0);
    CHECK(y1bar(11) == 19.0);
    CHECK(y1bar(12) == 18.0);
    CHECK(y1bar(13) == 17.0);
    CHECK(y1bar(14) == 16.0);
    CHECK(y1bar(15) == 15.0);
    CHECK(y1bar(16) == 14.0);
    CHECK(y1bar(17) == 13.0);
    CHECK(y1bar(18) == 12.0);
    CHECK(y1bar(19) == 11.0);
    CHECK(y1bar(20) == 10.0);
    CHECK(y1bar(21) == 9.0);
    CHECK(y1bar(22) == 8.0);
    CHECK(y1bar(23) == 7.0);
    CHECK(y1bar(24) == 6.0);
    CHECK(y1bar(25) == 5.0);
    CHECK(y1bar(26) == 4.0);
    CHECK(y1bar(27) == 3.0);
    CHECK(y1bar(28) == 2.0);
    CHECK(y1bar(29) == 1.0);
}

TEST_CASE("CompositeAlgebra get_matrix", "[domain]")
{
    /*!
    * Tests the get_matrix operation.
    */

    using namespace Lielab::domain;

    CompositeAlgebra y1 = CompositeAlgebra({ycn1, yglc1, yglr1, yrn1, yse1, yso1, ysp1, ysu1});
    y1.set_vector({30.0, 29.0, 28.0, 27.0, 26.0, 25.0, 24.0, 23.0, 22.0, 21.0,
                   20.0, 19.0, 18.0, 17.0, 16.0, 15.0, 14.0, 13.0, 12.0, 11.0,
                   10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0});

    const Eigen::MatrixXcd y1hat = y1.get_matrix();

    REQUIRE(y1hat.rows() == 20);
    REQUIRE(y1hat.cols() == 20);
    
    // TODO: Check the 0's
    // TODO: Check each submatrix

}

TEST_CASE("CompositeAlgebra operator()", "[domain]")
{
    using namespace Lielab::domain;

    const CompositeAlgebra y1 = CompositeAlgebra({ycn1, yglc1, yglr1, yrn1, yse1, yso1, ysp1, ysu1});

    // In bounds cn component
    CHECK(y1(0) == 1.0);
    CHECK(y1(1) == 2.0);
    CHECK(y1(2) == 3.0);
    CHECK(y1(3) == 4.0);
    CHECK(y1(-30) == 1.0);
    CHECK(y1(-29) == 2.0);
    CHECK(y1(-28) == 3.0);
    CHECK(y1(-27) == 4.0);

    // In bounds glc component
    CHECK(y1(4) == 5.0);
    CHECK(y1(5) == 6.0);
    CHECK(y1(6) == 7.0);
    CHECK(y1(7) == 8.0);
    CHECK(y1(8) == 9.0);
    CHECK(y1(9) == 10.0);
    CHECK(y1(10) == 11.0);
    CHECK(y1(11) == 12.0);
    CHECK(y1(-26) == 5.0);
    CHECK(y1(-25) == 6.0);
    CHECK(y1(-24) == 7.0);
    CHECK(y1(-23) == 8.0);
    CHECK(y1(-22) == 9.0);
    CHECK(y1(-21) == 10.0);
    CHECK(y1(-20) == 11.0);
    CHECK(y1(-19) == 12.0);

    // In bounds glr component
    CHECK(y1(12) == 13.0);
    CHECK(y1(13) == 14.0);
    CHECK(y1(14) == 15.0);
    CHECK(y1(15) == 16.0);
    CHECK(y1(-18) == 13.0);
    CHECK(y1(-17) == 14.0);
    CHECK(y1(-16) == 15.0);
    CHECK(y1(-15) == 16.0);

    // In bounds rn component
    CHECK(y1(16) == 17.0);
    CHECK(y1(17) == 18.0);
    CHECK(y1(-14) == 17.0);
    CHECK(y1(-13) == 18.0);

    // In bounds se component
    CHECK(y1(18) == 19.0);
    CHECK(y1(19) == 20.0);
    CHECK(y1(20) == 21.0);
    CHECK(y1(-12) == 19.0);
    CHECK(y1(-11) == 20.0);
    CHECK(y1(-10) == 21.0);

    // In bounds so component
    CHECK(y1(21) == 22.0);
    CHECK(y1(22) == 23.0);
    CHECK(y1(23) == 24.0);
    CHECK(y1(-9) == 22.0);
    CHECK(y1(-8) == 23.0);
    CHECK(y1(-7) == 24.0);

    // In bounds sp component
    CHECK(y1(24) == 25.0);
    CHECK(y1(25) == 26.0);
    CHECK(y1(26) == 27.0);
    CHECK(y1(-6) == 25.0);
    CHECK(y1(-5) == 26.0);
    CHECK(y1(-4) == 27.0);

    // In bounds su component
    CHECK(y1(27) == 28.0);
    CHECK(y1(28) == 29.0);
    CHECK(y1(29) == 30.0);
    CHECK(y1(-3) == 28.0);
    CHECK(y1(-2) == 29.0);
    CHECK(y1(-1) == 30.0);

    // Out of bounds
    CHECK(std::isnan(y1(-31)));
    CHECK(std::isnan(y1(30)));

    // In bounds cn component
    CHECK(y1(0, 0) == std::complex<double>(0.0, 0.0));
    CHECK(y1(0, 1) == std::complex<double>(0.0, 0.0));
    CHECK(y1(0, 2) == std::complex<double>(1.0, 2.0));
    CHECK(y1(1, 0) == std::complex<double>(0.0, 0.0));
    CHECK(y1(1, 1) == std::complex<double>(0.0, 0.0));
    CHECK(y1(1, 2) == std::complex<double>(3.0, 4.0));
    CHECK(y1(2, 0) == std::complex<double>(0.0, 0.0));
    CHECK(y1(2, 1) == std::complex<double>(0.0, 0.0));
    CHECK(y1(2, 2) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-20, -20) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-20, -19) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-20, -18) == std::complex<double>(1.0, 2.0));
    CHECK(y1(-19, -20) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-19, -19) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-19, -18) == std::complex<double>(3.0, 4.0));
    CHECK(y1(-18, -20) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-18, -19) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-18, -18) == std::complex<double>(0.0, 0.0));

    // In bounds glc component
    CHECK(y1(3, 3) == std::complex<double>(5.0, 6.0));
    CHECK(y1(3, 4) == std::complex<double>(7.0, 8.0));
    CHECK(y1(4, 3) == std::complex<double>(9.0, 10.0));
    CHECK(y1(4, 4) == std::complex<double>(11.0, 12.0));
    CHECK(y1(-17, -17) == std::complex<double>(5.0, 6.0));
    CHECK(y1(-17, -16) == std::complex<double>(7.0, 8.0));
    CHECK(y1(-16, -17) == std::complex<double>(9.0, 10.0));
    CHECK(y1(-16, -16) == std::complex<double>(11.0, 12.0));

    // In bounds glr component
    CHECK(y1(5, 5) == std::complex<double>(13.0, 0.0));
    CHECK(y1(5, 6) == std::complex<double>(14.0, 0.0));
    CHECK(y1(6, 5) == std::complex<double>(15.0, 0.0));
    CHECK(y1(6, 6) == std::complex<double>(16.0, 0.0));
    CHECK(y1(-15, -15) == std::complex<double>(13.0, 0.0));
    CHECK(y1(-15, -14) == std::complex<double>(14.0, 0.0));
    CHECK(y1(-14, -15) == std::complex<double>(15.0, 0.0));
    CHECK(y1(-14, -14) == std::complex<double>(16.0, 0.0));

    // In bounds rn component
    CHECK(y1(7, 7) == std::complex<double>(0.0, 0.0));
    CHECK(y1(7, 8) == std::complex<double>(0.0, 0.0));
    CHECK(y1(7, 9) == std::complex<double>(17.0, 0.0));
    CHECK(y1(8, 7) == std::complex<double>(0.0, 0.0));
    CHECK(y1(8, 8) == std::complex<double>(0.0, 0.0));
    CHECK(y1(8, 9) == std::complex<double>(18.0, 0.0));
    CHECK(y1(9, 7) == std::complex<double>(0.0, 0.0));
    CHECK(y1(9, 8) == std::complex<double>(0.0, 0.0));
    CHECK(y1(9, 9) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-13, -13) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-13, -12) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-13, -11) == std::complex<double>(17.0, 0.0));
    CHECK(y1(-12, -13) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-12, -12) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-12, -11) == std::complex<double>(18.0, 0.0));
    CHECK(y1(-11, -13) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-11, -12) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-11, -11) == std::complex<double>(0.0, 0.0));

    // In bounds se component
    CHECK(y1(10, 10) == std::complex<double>(0.0, 0.0));
    CHECK(y1(10, 11) == std::complex<double>(-21.0, 0.0));
    CHECK(y1(10, 12) == std::complex<double>(19.0, 0.0));
    CHECK(y1(11, 10) == std::complex<double>(21.0, 0.0));
    CHECK(y1(11, 11) == std::complex<double>(0.0, 0.0));
    CHECK(y1(11, 12) == std::complex<double>(20.0, 0.0));
    CHECK(y1(12, 10) == std::complex<double>(0.0, 0.0));
    CHECK(y1(12, 11) == std::complex<double>(0.0, 0.0));
    CHECK(y1(12, 12) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-10, -10) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-10, -9) == std::complex<double>(-21.0, 0.0));
    CHECK(y1(-10, -8) == std::complex<double>(19.0, 0.0));
    CHECK(y1(-9, -10) == std::complex<double>(21.0, 0.0));
    CHECK(y1(-9, -9) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-9, -8) == std::complex<double>(20.0, 0.0));
    CHECK(y1(-8, -10) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-8, -9) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-8, -8) == std::complex<double>(0.0, 0.0));

    // In bounds so component
    CHECK(y1(13, 13) == std::complex<double>(0.0, 0.0));
    CHECK(y1(13, 14) == std::complex<double>(-24.0, 0.0));
    CHECK(y1(13, 15) == std::complex<double>(23.0, 0.0));
    CHECK(y1(14, 13) == std::complex<double>(24.0, 0.0));
    CHECK(y1(14, 14) == std::complex<double>(0.0, 0.0));
    CHECK(y1(14, 15) == std::complex<double>(-22.0, 0.0));
    CHECK(y1(15, 13) == std::complex<double>(-23.0, 0.0));
    CHECK(y1(15, 14) == std::complex<double>(22.0, 0.0));
    CHECK(y1(15, 15) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-7, -7) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-7, -6) == std::complex<double>(-24.0, 0.0));
    CHECK(y1(-7, -5) == std::complex<double>(23.0, 0.0));
    CHECK(y1(-6, -7) == std::complex<double>(24.0, 0.0));
    CHECK(y1(-6, -6) == std::complex<double>(0.0, 0.0));
    CHECK(y1(-6, -5) == std::complex<double>(-22.0, 0.0));
    CHECK(y1(-5, -7) == std::complex<double>(-23.0, 0.0));
    CHECK(y1(-5, -6) == std::complex<double>(22.0, 0.0));
    CHECK(y1(-5, -5) == std::complex<double>(0.0, 0.0));

    // In bounds sp component
    CHECK(y1(16, 16) == std::complex<double>(25.0, 0.0));
    CHECK(y1(16, 17) == std::complex<double>(26.0, 0.0));
    CHECK(y1(17, 16) == std::complex<double>(27.0, 0.0));
    CHECK(y1(17, 17) == std::complex<double>(-25.0, 0.0));
    CHECK(y1(-4, -4) == std::complex<double>(25.0, 0.0));
    CHECK(y1(-4, -3) == std::complex<double>(26.0, 0.0));
    CHECK(y1(-3, -4) == std::complex<double>(27.0, 0.0));
    CHECK(y1(-3, -3) == std::complex<double>(-25.0, 0.0));

    // In bounds su component
    CHECK(y1(18, 18) == std::complex<double>(0.0, 30.0));
    CHECK(y1(18, 19) == std::complex<double>(-29.0, 28.0));
    CHECK(y1(19, 18) == std::complex<double>(29.0, 28.0));
    CHECK(y1(19, 19) == std::complex<double>(0.0, -30.0));
    CHECK(y1(-2, -2) == std::complex<double>(0.0, 30.0));
    CHECK(y1(-2, -1) == std::complex<double>(-29.0, 28.0));
    CHECK(y1(-1, -2) == std::complex<double>(29.0, 28.0));
    CHECK(y1(-1, -1) == std::complex<double>(0.0, -30.0));

    // In bounds sparse components
    CHECK(y1(0, 3) == std::complex<double>(0.0, 0.0));
    CHECK(y1(0, 4) == std::complex<double>(0.0, 0.0));
    CHECK(y1(1, 3) == std::complex<double>(0.0, 0.0));
    CHECK(y1(1, 4) == std::complex<double>(0.0, 0.0));
    CHECK(y1(2, 3) == std::complex<double>(0.0, 0.0));
    CHECK(y1(2, 4) == std::complex<double>(0.0, 0.0));
    CHECK(y1(3, 0) == std::complex<double>(0.0, 0.0));
    CHECK(y1(4, 0) == std::complex<double>(0.0, 0.0));
    CHECK(y1(3, 1) == std::complex<double>(0.0, 0.0));
    CHECK(y1(4, 1) == std::complex<double>(0.0, 0.0));
    CHECK(y1(3, 2) == std::complex<double>(0.0, 0.0));
    CHECK(y1(4, 2) == std::complex<double>(0.0, 0.0));
    // Probably don't need to check the rest

    // Out of bounds
    CHECK(std::isnan(y1(0, -21).real()));
    CHECK(std::isnan(y1(0, -21).imag()));
    CHECK(std::isnan(y1(-21, 0).real()));
    CHECK(std::isnan(y1(-21, 0).imag()));
    CHECK(std::isnan(y1(-21, -21).real()));
    CHECK(std::isnan(y1(-21, -21).imag()));
    CHECK(std::isnan(y1(0, 20).real()));
    CHECK(std::isnan(y1(0, 20).imag()));
    CHECK(std::isnan(y1(20, 0).real()));
    CHECK(std::isnan(y1(20, 0).imag()));
    CHECK(std::isnan(y1(20, 20).real()));
    CHECK(std::isnan(y1(20, 20).imag()));
}

// TODO: mathops int here

TEST_CASE("CompositeAlgebra math_ops_double", "[domain]")
{
    using namespace Lielab::domain;

    CompositeAlgebra x1 = CompositeAlgebra({ycn1, yglc1, yglr1, yrn1, yse1, yso1, ysp1, ysu1});

    const CompositeAlgebra x1_lm_2 = 2.0*x1;
    const Eigen::VectorXd x1_lm_2bar = x1_lm_2.get_vector();
    REQUIRE(x1_lm_2bar.size() == 30);
    CHECK(x1_lm_2bar(0) == 2.0);
    CHECK(x1_lm_2bar(1) == 4.0);
    CHECK(x1_lm_2bar(2) == 6.0);
    CHECK(x1_lm_2bar(3) == 8.0);
    CHECK(x1_lm_2bar(4) == 10.0);
    CHECK(x1_lm_2bar(5) == 12.0);
    CHECK(x1_lm_2bar(6) == 14.0);
    CHECK(x1_lm_2bar(7) == 16.0);
    CHECK(x1_lm_2bar(8) == 18.0);
    CHECK(x1_lm_2bar(9) == 20.0);
    CHECK(x1_lm_2bar(10) == 22.0);
    CHECK(x1_lm_2bar(11) == 24.0);
    CHECK(x1_lm_2bar(12) == 26.0);
    CHECK(x1_lm_2bar(13) == 28.0);
    CHECK(x1_lm_2bar(14) == 30.0);
    CHECK(x1_lm_2bar(15) == 32.0);
    CHECK(x1_lm_2bar(16) == 34.0);
    CHECK(x1_lm_2bar(17) == 36.0);
    CHECK(x1_lm_2bar(18) == 38.0);
    CHECK(x1_lm_2bar(19) == 40.0);
    CHECK(x1_lm_2bar(20) == 42.0);
    CHECK(x1_lm_2bar(21) == 44.0);
    CHECK(x1_lm_2bar(22) == 46.0);
    CHECK(x1_lm_2bar(23) == 48.0);
    CHECK(x1_lm_2bar(24) == 50.0);
    CHECK(x1_lm_2bar(25) == 52.0);
    CHECK(x1_lm_2bar(26) == 54.0);
    CHECK(x1_lm_2bar(27) == 56.0);
    CHECK(x1_lm_2bar(28) == 58.0);
    CHECK(x1_lm_2bar(29) == 60.0);

    const CompositeAlgebra x1_rm_2 = x1*2.0;
    const Eigen::VectorXd x1_rm_2bar = x1_rm_2.get_vector();
    REQUIRE(x1_rm_2bar.size() == 30);
    CHECK(x1_rm_2bar(0) == 2.0);
    CHECK(x1_rm_2bar(1) == 4.0);
    CHECK(x1_rm_2bar(2) == 6.0);
    CHECK(x1_rm_2bar(3) == 8.0);
    CHECK(x1_rm_2bar(4) == 10.0);
    CHECK(x1_rm_2bar(5) == 12.0);
    CHECK(x1_rm_2bar(6) == 14.0);
    CHECK(x1_rm_2bar(7) == 16.0);
    CHECK(x1_rm_2bar(8) == 18.0);
    CHECK(x1_rm_2bar(9) == 20.0);
    CHECK(x1_rm_2bar(10) == 22.0);
    CHECK(x1_rm_2bar(11) == 24.0);
    CHECK(x1_rm_2bar(12) == 26.0);
    CHECK(x1_rm_2bar(13) == 28.0);
    CHECK(x1_rm_2bar(14) == 30.0);
    CHECK(x1_rm_2bar(15) == 32.0);
    CHECK(x1_rm_2bar(16) == 34.0);
    CHECK(x1_rm_2bar(17) == 36.0);
    CHECK(x1_rm_2bar(18) == 38.0);
    CHECK(x1_rm_2bar(19) == 40.0);
    CHECK(x1_rm_2bar(20) == 42.0);
    CHECK(x1_rm_2bar(21) == 44.0);
    CHECK(x1_rm_2bar(22) == 46.0);
    CHECK(x1_rm_2bar(23) == 48.0);
    CHECK(x1_rm_2bar(24) == 50.0);
    CHECK(x1_rm_2bar(25) == 52.0);
    CHECK(x1_rm_2bar(26) == 54.0);
    CHECK(x1_rm_2bar(27) == 56.0);
    CHECK(x1_rm_2bar(28) == 58.0);
    CHECK(x1_rm_2bar(29) == 60.0);

    x1 *= 2.0;
    const Eigen::VectorXd x1_im_2bar = x1.get_vector();
    REQUIRE(x1_im_2bar.size() == 30);
    CHECK(x1_im_2bar(0) == 2.0);
    CHECK(x1_im_2bar(1) == 4.0);
    CHECK(x1_im_2bar(2) == 6.0);
    CHECK(x1_im_2bar(3) == 8.0);
    CHECK(x1_im_2bar(4) == 10.0);
    CHECK(x1_im_2bar(5) == 12.0);
    CHECK(x1_im_2bar(6) == 14.0);
    CHECK(x1_im_2bar(7) == 16.0);
    CHECK(x1_im_2bar(8) == 18.0);
    CHECK(x1_im_2bar(9) == 20.0);
    CHECK(x1_im_2bar(10) == 22.0);
    CHECK(x1_im_2bar(11) == 24.0);
    CHECK(x1_im_2bar(12) == 26.0);
    CHECK(x1_im_2bar(13) == 28.0);
    CHECK(x1_im_2bar(14) == 30.0);
    CHECK(x1_im_2bar(15) == 32.0);
    CHECK(x1_im_2bar(16) == 34.0);
    CHECK(x1_im_2bar(17) == 36.0);
    CHECK(x1_im_2bar(18) == 38.0);
    CHECK(x1_im_2bar(19) == 40.0);
    CHECK(x1_im_2bar(20) == 42.0);
    CHECK(x1_im_2bar(21) == 44.0);
    CHECK(x1_im_2bar(22) == 46.0);
    CHECK(x1_im_2bar(23) == 48.0);
    CHECK(x1_im_2bar(24) == 50.0);
    CHECK(x1_im_2bar(25) == 52.0);
    CHECK(x1_im_2bar(26) == 54.0);
    CHECK(x1_im_2bar(27) == 56.0);
    CHECK(x1_im_2bar(28) == 58.0);
    CHECK(x1_im_2bar(29) == 60.0);

    x1 = CompositeAlgebra({ycn1, yglc1, yglr1, yrn1, yse1, yso1, ysp1, ysu1});

    const CompositeAlgebra x1_d_2 = x1/2.0;
    const Eigen::VectorXd x1_d_2bar = x1_d_2.get_vector();
    REQUIRE(x1_d_2bar.size() == 30);
    CHECK(x1_d_2bar(0) == 0.5);
    CHECK(x1_d_2bar(1) == 1.0);
    CHECK(x1_d_2bar(2) == 1.5);
    CHECK(x1_d_2bar(3) == 2.0);
    CHECK(x1_d_2bar(4) == 2.5);
    CHECK(x1_d_2bar(5) == 3.0);
    CHECK(x1_d_2bar(6) == 3.5);
    CHECK(x1_d_2bar(7) == 4.0);
    CHECK(x1_d_2bar(8) == 4.5);
    CHECK(x1_d_2bar(9) == 5.0);
    CHECK(x1_d_2bar(10) == 5.5);
    CHECK(x1_d_2bar(11) == 6.0);
    CHECK(x1_d_2bar(12) == 6.5);
    CHECK(x1_d_2bar(13) == 7.0);
    CHECK(x1_d_2bar(14) == 7.5);
    CHECK(x1_d_2bar(15) == 8.0);
    CHECK(x1_d_2bar(16) == 8.5);
    CHECK(x1_d_2bar(17) == 9.0);
    CHECK(x1_d_2bar(18) == 9.5);
    CHECK(x1_d_2bar(19) == 10.0);
    CHECK(x1_d_2bar(20) == 10.5);
    CHECK(x1_d_2bar(21) == 11.0);
    CHECK(x1_d_2bar(22) == 11.5);
    CHECK(x1_d_2bar(23) == 12.0);
    CHECK(x1_d_2bar(24) == 12.5);
    CHECK(x1_d_2bar(25) == 13.0);
    CHECK(x1_d_2bar(26) == 13.5);
    CHECK(x1_d_2bar(27) == 14.0);
    CHECK(x1_d_2bar(28) == 14.5);
    CHECK(x1_d_2bar(29) == 15.0);

    x1 /= 2.0;
    const Eigen::VectorXd x1_id_2bar = x1.get_vector();
    REQUIRE(x1_id_2bar.size() == 30);
    CHECK(x1_id_2bar(0) == 0.5);
    CHECK(x1_id_2bar(1) == 1.0);
    CHECK(x1_id_2bar(2) == 1.5);
    CHECK(x1_id_2bar(3) == 2.0);
    CHECK(x1_id_2bar(4) == 2.5);
    CHECK(x1_id_2bar(5) == 3.0);
    CHECK(x1_id_2bar(6) == 3.5);
    CHECK(x1_id_2bar(7) == 4.0);
    CHECK(x1_id_2bar(8) == 4.5);
    CHECK(x1_id_2bar(9) == 5.0);
    CHECK(x1_id_2bar(10) == 5.5);
    CHECK(x1_id_2bar(11) == 6.0);
    CHECK(x1_id_2bar(12) == 6.5);
    CHECK(x1_id_2bar(13) == 7.0);
    CHECK(x1_id_2bar(14) == 7.5);
    CHECK(x1_id_2bar(15) == 8.0);
    CHECK(x1_id_2bar(16) == 8.5);
    CHECK(x1_id_2bar(17) == 9.0);
    CHECK(x1_id_2bar(18) == 9.5);
    CHECK(x1_id_2bar(19) == 10.0);
    CHECK(x1_id_2bar(20) == 10.5);
    CHECK(x1_id_2bar(21) == 11.0);
    CHECK(x1_id_2bar(22) == 11.5);
    CHECK(x1_id_2bar(23) == 12.0);
    CHECK(x1_id_2bar(24) == 12.5);
    CHECK(x1_id_2bar(25) == 13.0);
    CHECK(x1_id_2bar(26) == 13.5);
    CHECK(x1_id_2bar(27) == 14.0);
    CHECK(x1_id_2bar(28) == 14.5);
    CHECK(x1_id_2bar(29) == 15.0);
}

TEST_CASE("CompositeAlgebra math_ops_CompositeAlgebra", "[domain]")
{
    using namespace Lielab::domain;

    CompositeAlgebra x1 = CompositeAlgebra({ycn1, yglc1, yglr1, yrn1, yse1, yso1, ysp1, ysu1});
    CompositeAlgebra x2 = CompositeAlgebra({ycn1, yglc1, yglr1, yrn1, yse1, yso1, ysp1, ysu1});

    const CompositeAlgebra x1_add_x2 = x1 + x2;
    const Eigen::VectorXd x1_add_x2bar = x1_add_x2.get_vector();
    REQUIRE(x1_add_x2bar.size() == 30);
    CHECK(x1_add_x2bar(0) == 2.0);
    CHECK(x1_add_x2bar(1) == 4.0);
    CHECK(x1_add_x2bar(2) == 6.0);
    CHECK(x1_add_x2bar(3) == 8.0);
    CHECK(x1_add_x2bar(4) == 10.0);
    CHECK(x1_add_x2bar(5) == 12.0);
    CHECK(x1_add_x2bar(6) == 14.0);
    CHECK(x1_add_x2bar(7) == 16.0);
    CHECK(x1_add_x2bar(8) == 18.0);
    CHECK(x1_add_x2bar(9) == 20.0);
    CHECK(x1_add_x2bar(10) == 22.0);
    CHECK(x1_add_x2bar(11) == 24.0);
    CHECK(x1_add_x2bar(12) == 26.0);
    CHECK(x1_add_x2bar(13) == 28.0);
    CHECK(x1_add_x2bar(14) == 30.0);
    CHECK(x1_add_x2bar(15) == 32.0);
    CHECK(x1_add_x2bar(16) == 34.0);
    CHECK(x1_add_x2bar(17) == 36.0);
    CHECK(x1_add_x2bar(18) == 38.0);
    CHECK(x1_add_x2bar(19) == 40.0);
    CHECK(x1_add_x2bar(20) == 42.0);
    CHECK(x1_add_x2bar(21) == 44.0);
    CHECK(x1_add_x2bar(22) == 46.0);
    CHECK(x1_add_x2bar(23) == 48.0);
    CHECK(x1_add_x2bar(24) == 50.0);
    CHECK(x1_add_x2bar(25) == 52.0);
    CHECK(x1_add_x2bar(26) == 54.0);
    CHECK(x1_add_x2bar(27) == 56.0);
    CHECK(x1_add_x2bar(28) == 58.0);
    CHECK(x1_add_x2bar(29) == 60.0);

    x1 += x2;
    const Eigen::VectorXd x1_iadd_x2bar = x1.get_vector();
    REQUIRE(x1_iadd_x2bar.size() == 30);
    CHECK(x1_iadd_x2bar(0) == 2.0);
    CHECK(x1_iadd_x2bar(1) == 4.0);
    CHECK(x1_iadd_x2bar(2) == 6.0);
    CHECK(x1_iadd_x2bar(3) == 8.0);
    CHECK(x1_iadd_x2bar(4) == 10.0);
    CHECK(x1_iadd_x2bar(5) == 12.0);
    CHECK(x1_iadd_x2bar(6) == 14.0);
    CHECK(x1_iadd_x2bar(7) == 16.0);
    CHECK(x1_iadd_x2bar(8) == 18.0);
    CHECK(x1_iadd_x2bar(9) == 20.0);
    CHECK(x1_iadd_x2bar(10) == 22.0);
    CHECK(x1_iadd_x2bar(11) == 24.0);
    CHECK(x1_iadd_x2bar(12) == 26.0);
    CHECK(x1_iadd_x2bar(13) == 28.0);
    CHECK(x1_iadd_x2bar(14) == 30.0);
    CHECK(x1_iadd_x2bar(15) == 32.0);
    CHECK(x1_iadd_x2bar(16) == 34.0);
    CHECK(x1_iadd_x2bar(17) == 36.0);
    CHECK(x1_iadd_x2bar(18) == 38.0);
    CHECK(x1_iadd_x2bar(19) == 40.0);
    CHECK(x1_iadd_x2bar(20) == 42.0);
    CHECK(x1_iadd_x2bar(21) == 44.0);
    CHECK(x1_iadd_x2bar(22) == 46.0);
    CHECK(x1_iadd_x2bar(23) == 48.0);
    CHECK(x1_iadd_x2bar(24) == 50.0);
    CHECK(x1_iadd_x2bar(25) == 52.0);
    CHECK(x1_iadd_x2bar(26) == 54.0);
    CHECK(x1_iadd_x2bar(27) == 56.0);
    CHECK(x1_iadd_x2bar(28) == 58.0);
    CHECK(x1_iadd_x2bar(29) == 60.0);

    x1 = CompositeAlgebra({ycn1, yglc1, yglr1, yrn1, yse1, yso1, ysp1, ysu1});

    const CompositeAlgebra x1_sub_x2 = x1 - x2;
    const Eigen::VectorXd x1_sub_x2bar = x1_sub_x2.get_vector();
    REQUIRE(x1_sub_x2bar.size() == 30);
    CHECK(x1_sub_x2bar(0) == 0.0);
    CHECK(x1_sub_x2bar(1) == 0.0);
    CHECK(x1_sub_x2bar(2) == 0.0);
    CHECK(x1_sub_x2bar(3) == 0.0);
    CHECK(x1_sub_x2bar(4) == 0.0);
    CHECK(x1_sub_x2bar(5) == 0.0);
    CHECK(x1_sub_x2bar(6) == 0.0);
    CHECK(x1_sub_x2bar(7) == 0.0);
    CHECK(x1_sub_x2bar(8) == 0.0);
    CHECK(x1_sub_x2bar(9) == 0.0);
    CHECK(x1_sub_x2bar(10) == 0.0);
    CHECK(x1_sub_x2bar(11) == 0.0);
    CHECK(x1_sub_x2bar(12) == 0.0);
    CHECK(x1_sub_x2bar(13) == 0.0);
    CHECK(x1_sub_x2bar(14) == 0.0);
    CHECK(x1_sub_x2bar(15) == 0.0);
    CHECK(x1_sub_x2bar(16) == 0.0);
    CHECK(x1_sub_x2bar(17) == 0.0);
    CHECK(x1_sub_x2bar(18) == 0.0);
    CHECK(x1_sub_x2bar(19) == 0.0);
    CHECK(x1_sub_x2bar(20) == 0.0);
    CHECK(x1_sub_x2bar(21) == 0.0);
    CHECK(x1_sub_x2bar(22) == 0.0);
    CHECK(x1_sub_x2bar(23) == 0.0);
    CHECK(x1_sub_x2bar(24) == 0.0);
    CHECK(x1_sub_x2bar(25) == 0.0);
    CHECK(x1_sub_x2bar(26) == 0.0);
    CHECK(x1_sub_x2bar(27) == 0.0);
    CHECK(x1_sub_x2bar(28) == 0.0);
    CHECK(x1_sub_x2bar(29) == 0.0);

    x1 -= x2;
    const Eigen::VectorXd x1_isub_x2bar = x1.get_vector();
    REQUIRE(x1_isub_x2bar.size() == 30);
    CHECK(x1_isub_x2bar(0) == 0.0);
    CHECK(x1_isub_x2bar(1) == 0.0);
    CHECK(x1_isub_x2bar(2) == 0.0);
    CHECK(x1_isub_x2bar(3) == 0.0);
    CHECK(x1_isub_x2bar(4) == 0.0);
    CHECK(x1_isub_x2bar(5) == 0.0);
    CHECK(x1_isub_x2bar(6) == 0.0);
    CHECK(x1_isub_x2bar(7) == 0.0);
    CHECK(x1_isub_x2bar(8) == 0.0);
    CHECK(x1_isub_x2bar(9) == 0.0);
    CHECK(x1_isub_x2bar(10) == 0.0);
    CHECK(x1_isub_x2bar(11) == 0.0);
    CHECK(x1_isub_x2bar(12) == 0.0);
    CHECK(x1_isub_x2bar(13) == 0.0);
    CHECK(x1_isub_x2bar(14) == 0.0);
    CHECK(x1_isub_x2bar(15) == 0.0);
    CHECK(x1_isub_x2bar(16) == 0.0);
    CHECK(x1_isub_x2bar(17) == 0.0);
    CHECK(x1_isub_x2bar(18) == 0.0);
    CHECK(x1_isub_x2bar(19) == 0.0);
    CHECK(x1_isub_x2bar(20) == 0.0);
    CHECK(x1_isub_x2bar(21) == 0.0);
    CHECK(x1_isub_x2bar(22) == 0.0);
    CHECK(x1_isub_x2bar(23) == 0.0);
    CHECK(x1_isub_x2bar(24) == 0.0);
    CHECK(x1_isub_x2bar(25) == 0.0);
    CHECK(x1_isub_x2bar(26) == 0.0);
    CHECK(x1_isub_x2bar(27) == 0.0);
    CHECK(x1_isub_x2bar(28) == 0.0);
    CHECK(x1_isub_x2bar(29) == 0.0);

    x1 = CompositeAlgebra({ycn1, yglc1, yglr1, yrn1, yse1, yso1, ysp1, ysu1});

    const CompositeAlgebra x1_unary_sub = (-x1);
    const Eigen::VectorXd x1_usub_bar = x1_unary_sub.get_vector();
    REQUIRE(x1_usub_bar.size() == 30);
    CHECK(x1_usub_bar(0) == -1.0);
    CHECK(x1_usub_bar(1) == -2.0);
    CHECK(x1_usub_bar(2) == -3.0);
    CHECK(x1_usub_bar(3) == -4.0);
    CHECK(x1_usub_bar(4) == -5.0);
    CHECK(x1_usub_bar(5) == -6.0);
    CHECK(x1_usub_bar(6) == -7.0);
    CHECK(x1_usub_bar(7) == -8.0);
    CHECK(x1_usub_bar(8) == -9.0);
    CHECK(x1_usub_bar(9) == -10.0);
    CHECK(x1_usub_bar(10) == -11.0);
    CHECK(x1_usub_bar(11) == -12.0);
    CHECK(x1_usub_bar(12) == -13.0);
    CHECK(x1_usub_bar(13) == -14.0);
    CHECK(x1_usub_bar(14) == -15.0);
    CHECK(x1_usub_bar(15) == -16.0);
    CHECK(x1_usub_bar(16) == -17.0);
    CHECK(x1_usub_bar(17) == -18.0);
    CHECK(x1_usub_bar(18) == -19.0);
    CHECK(x1_usub_bar(19) == -20.0);
    CHECK(x1_usub_bar(20) == -21.0);
    CHECK(x1_usub_bar(21) == -22.0);
    CHECK(x1_usub_bar(22) == -23.0);
    CHECK(x1_usub_bar(23) == -24.0);
    CHECK(x1_usub_bar(24) == -25.0);
    CHECK(x1_usub_bar(25) == -26.0);
    CHECK(x1_usub_bar(26) == -27.0);
    CHECK(x1_usub_bar(27) == -28.0);
    CHECK(x1_usub_bar(28) == -29.0);
    CHECK(x1_usub_bar(29) == -30.0);
}

TEST_CASE("CompositeAlgebra operator[]", "[domain]")
{
    using namespace Lielab::domain;

    const CompositeAlgebra x1 = CompositeAlgebra({ycn1, yglc1, yglr1, yrn1, yse1, yso1, ysp1, ysu1});
    
    // It would be nice if these could be called without std::get<> like:
    // const cn x10 = x1[0];
    const cn x10 = std::get<cn>(x1[0]);
    const glc x11 = std::get<glc>(x1[1]);
    const glr x12 = std::get<glr>(x1[2]);
    const rn x13 = std::get<rn>(x1[3]);
    const se x14 = std::get<se>(x1[4]);
    const so x15 = std::get<so>(x1[5]);
    const sp x16 = std::get<sp>(x1[6]);
    const su x17 = std::get<su>(x1[7]);
    const cn x1m8 = std::get<cn>(x1[-8]);
    const glc x1m7 = std::get<glc>(x1[-7]);
    const glr x1m6 = std::get<glr>(x1[-6]);
    const rn x1m5 = std::get<rn>(x1[-5]);
    const se x1m4 = std::get<se>(x1[-4]);
    const so x1m3 = std::get<so>(x1[-3]);
    const sp x1m2 = std::get<sp>(x1[-2]);
    const su x1m1 = std::get<su>(x1[-1]);

    CHECK(x10.to_string() == "c^2");
    CHECK(x11.to_string() == "gl(2, C)");
    CHECK(x12.to_string() == "gl(2, R)");
    CHECK(x13.to_string() == "r^2");
    CHECK(x14.to_string() == "se(2)");
    CHECK(x15.to_string() == "so(3)");
    CHECK(x16.to_string() == "sp(2, R)");
    CHECK(x17.to_string() == "su(2)");
    CHECK(x1m8.to_string() == "c^2");
    CHECK(x1m7.to_string() == "gl(2, C)");
    CHECK(x1m6.to_string() == "gl(2, R)");
    CHECK(x1m5.to_string() == "r^2");
    CHECK(x1m4.to_string() == "se(2)");
    CHECK(x1m3.to_string() == "so(3)");
    CHECK(x1m2.to_string() == "sp(2, R)");
    CHECK(x1m1.to_string() == "su(2)");

    // Out of bounds
    const glc x18 = std::get<glc>(x1[8]);
    const glc x1m9 = std::get<glc>(x1[-9]);

    CHECK(x18.to_string() == "gl(0, C)");
    CHECK(x1m9.to_string() == "gl(0, C)");
}
