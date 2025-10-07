#include <catch2/catch_all.hpp>
#include <Eigen/Core>

#include <Lielab.hpp>
#include "../test_utils.hpp"

TEST_CASE("concatenate", "[utils]")
{
    using namespace Lielab::utils;

    std::vector<Eigen::VectorXi> li;
    std::vector<Eigen::VectorXd> ld;

    REQUIRE(concatenate<int>({}).size() == 0);
    REQUIRE(concatenate<double>({}).size() == 0);
    REQUIRE(concatenate<int>(li).size() == 0);
    REQUIRE(concatenate<double>(ld).size() == 0);

    Eigen::VectorXi i1 = Eigen::VectorXi::Zero(3);
    i1(0) = 1;
    i1(1) = 2;
    i1(2) = 3;

    Eigen::VectorXd d1 = Eigen::VectorXd::Zero(3);
    d1(0) = 1.5;
    d1(1) = 2.5;
    d1(2) = 3.5;

    const Eigen::VectorXi t1 = concatenate<int>({i1});
    REQUIRE(t1.size() == 3);
    CHECK(t1(0) == 1);
    CHECK(t1(1) == 2);
    CHECK(t1(2) == 3);

    const Eigen::VectorXd t2 = concatenate<double>({d1});
    REQUIRE(t2.size() == 3);
    CHECK(t2(0) == 1.5);
    CHECK(t2(1) == 2.5);
    CHECK(t2(2) == 3.5);

    li.push_back(i1);
    ld.push_back(d1);

    const Eigen::VectorXi t3 = concatenate<int>(li);
    REQUIRE(t3.size() == 3);
    CHECK(t3(0) == 1);
    CHECK(t3(1) == 2);
    CHECK(t3(2) == 3);

    const Eigen::VectorXd t4 = concatenate<double>(ld);
    REQUIRE(t4.size() == 3);
    CHECK(t4(0) == 1.5);
    CHECK(t4(1) == 2.5);
    CHECK(t4(2) == 3.5);

    Eigen::VectorXi i2 = Eigen::VectorXi::Zero(4);
    i2(0) = 4;
    i2(1) = 5;
    i2(2) = 6;
    i2(3) = 7;

    Eigen::VectorXd d2 = Eigen::VectorXd::Zero(4);
    d2(0) = 4.5;
    d2(1) = 5.5;
    d2(2) = 6.5;
    d2(3) = 7.5;

    const Eigen::VectorXi t5 = concatenate<int>({i1, i2});
    REQUIRE(t5.size() == 7);
    CHECK(t5(0) == 1);
    CHECK(t5(1) == 2);
    CHECK(t5(2) == 3);
    CHECK(t5(3) == 4);
    CHECK(t5(4) == 5);
    CHECK(t5(5) == 6);
    CHECK(t5(6) == 7);

    const Eigen::VectorXd t6 = concatenate<double>({d1, d2});
    REQUIRE(t6.size() == 7);
    CHECK(t6(0) == 1.5);
    CHECK(t6(1) == 2.5);
    CHECK(t6(2) == 3.5);
    CHECK(t6(3) == 4.5);
    CHECK(t6(4) == 5.5);
    CHECK(t6(5) == 6.5);
    CHECK(t6(6) == 7.5);

    li.push_back(i2);
    ld.push_back(d2);

    const Eigen::VectorXi t7 = concatenate<int>(li);
    REQUIRE(t7.size() == 7);
    CHECK(t7(0) == 1);
    CHECK(t7(1) == 2);
    CHECK(t7(2) == 3);
    CHECK(t7(3) == 4);
    CHECK(t7(4) == 5);
    CHECK(t7(5) == 6);
    CHECK(t7(6) == 7);

    const Eigen::VectorXd t8 = concatenate<double>(ld);
    REQUIRE(t8.size() == 7);
    CHECK(t8(0) == 1.5);
    CHECK(t8(1) == 2.5);
    CHECK(t8(2) == 3.5);
    CHECK(t8(3) == 4.5);
    CHECK(t8(4) == 5.5);
    CHECK(t8(5) == 6.5);
    CHECK(t8(6) == 7.5);
}

TEST_CASE("arange", "[utils]")
{
    using namespace Lielab::utils;

    REQUIRE(arange<int>(0, 0).size() == 0);
    REQUIRE(arange<double>(0, 0).size() == 0);
    REQUIRE(arange<int>(0).size() == 0);
    REQUIRE(arange<double>(0).size() == 0);

    REQUIRE(arange<int>(1, 0).size() == 0);
    REQUIRE(arange<double>(1, 0).size() == 0);

    const Eigen::VectorXi t1 = arange<int>(0, 1);
    REQUIRE(t1.size() == 1);
    CHECK(t1(0) == 0);
    const Eigen::VectorXd t2 = arange<double>(0.5, 1.5);
    REQUIRE(t2.size() == 1);
    CHECK(t2(0) == 0.5);
    const Eigen::VectorXi t3 = arange<int>(1);
    REQUIRE(t3.size() == 1);
    CHECK(t3(0) == 0);
    const Eigen::VectorXd t4 = arange<double>(0.5);
    REQUIRE(t4.size() == 1);
    CHECK(t4(0) == 0.0);

    const Eigen::VectorXi t5 = arange<int>(2, 6);
    REQUIRE(t5.size() == 4);
    CHECK(t5(0) == 2);
    CHECK(t5(1) == 3);
    CHECK(t5(2) == 4);
    CHECK(t5(3) == 5);
    const Eigen::VectorXd t6 = arange<double>(2.5, 6.5);
    REQUIRE(t6.size() == 4);
    CHECK(t6(0) == 2.5);
    CHECK(t6(1) == 3.5);
    CHECK(t6(2) == 4.5);
    CHECK(t6(3) == 5.5);
    const Eigen::VectorXi t7 = arange<int>(4);
    REQUIRE(t7.size() == 4);
    CHECK(t7(0) == 0);
    CHECK(t7(1) == 1);
    CHECK(t7(2) == 2);
    CHECK(t7(3) == 3);
    const Eigen::VectorXd t8 = arange<double>(3.5);
    REQUIRE(t8.size() == 4);
    CHECK(t8(0) == 0.0);
    CHECK(t8(1) == 1.0);
    CHECK(t8(2) == 2.0);
    CHECK(t8(3) == 3.0);

    const Eigen::VectorXi t9 = arange<int>(-6, -2);
    REQUIRE(t9.size() == 4);
    CHECK(t9(0) == -6);
    CHECK(t9(1) == -5);
    CHECK(t9(2) == -4);
    CHECK(t9(3) == -3);
    const Eigen::VectorXd t10 = arange<double>(-6.5, -2.5);
    REQUIRE(t10.size() == 4);
    CHECK(t10(0) == -6.5);
    CHECK(t10(1) == -5.5);
    CHECK(t10(2) == -4.5);
    CHECK(t10(3) == -3.5);
    const Eigen::VectorXi t11 = arange<int>(-4);
    REQUIRE(t11.size() == 0);
    const Eigen::VectorXd t12 = arange<double>(-3.5);
    REQUIRE(t12.size() == 0);

    const Eigen::VectorXd t13 = arange<double>(2.5, 6.5 + 1e-8);
    REQUIRE(t13.size() == 5);
    CHECK(t13(0) == 2.5);
    CHECK(t13(1) == 3.5);
    CHECK(t13(2) == 4.5);
    CHECK(t13(3) == 5.5);
    CHECK(t13(4) == 6.5);
    const Eigen::VectorXd t14 = arange<double>(4.0 + 1e-8);
    REQUIRE(t14.size() == 5);
    CHECK(t14(0) == 0.0);
    CHECK(t14(1) == 1.0);
    CHECK(t14(2) == 2.0);
    CHECK(t14(3) == 3.0);
    CHECK(t14(4) == 4.0);
}

TEST_CASE("repeat", "[utils]")
{
    using namespace Lielab::utils;

    std::vector<Eigen::VectorXi> li;
    std::vector<Eigen::VectorXd> ld;

    REQUIRE(repeat<int>({}, -1).size() == 0);
    REQUIRE(repeat<double>({}, -1).size() == 0);
    REQUIRE(repeat<int>(Eigen::VectorXi::Zero(0), -1).size() == 0);
    REQUIRE(repeat<double>(Eigen::VectorXd::Zero(0), -1).size() == 0);

    REQUIRE(repeat<int>({}, 0).size() == 0);
    REQUIRE(repeat<double>({}, 0).size() == 0);
    REQUIRE(repeat<int>(Eigen::VectorXi::Zero(0), 0).size() == 0);
    REQUIRE(repeat<double>(Eigen::VectorXd::Zero(0), 0).size() == 0);

    REQUIRE(repeat<int>({}, 1).size() == 0);
    REQUIRE(repeat<double>({}, 1).size() == 0);
    REQUIRE(repeat<int>(Eigen::VectorXi::Zero(0), 1).size() == 0);
    REQUIRE(repeat<double>(Eigen::VectorXd::Zero(0), 1).size() == 0);

    REQUIRE(repeat<int>({}, 2).size() == 0);
    REQUIRE(repeat<double>({}, 2).size() == 0);
    REQUIRE(repeat<int>(Eigen::VectorXi::Zero(0), 2).size() == 0);
    REQUIRE(repeat<double>(Eigen::VectorXd::Zero(0), 2).size() == 0);

    const Eigen::VectorXi t1 = repeat<int>({1, 2, 3}, -1);
    REQUIRE(t1.size() == 0);
    const Eigen::VectorXi t2 = repeat<int>({1, 2, 3}, 0);
    REQUIRE(t2.size() == 0);
    const Eigen::VectorXi t3 = repeat<int>({1, 2, 3}, 1);
    REQUIRE(t3.size() == 3);
    CHECK(t3(0) == 1);
    CHECK(t3(1) == 2);
    CHECK(t3(2) == 3);
    const Eigen::VectorXi t4 = repeat<int>({1, 2, 3}, 2);
    REQUIRE(t4.size() == 6);
    CHECK(t4(0) == 1);
    CHECK(t4(1) == 1);
    CHECK(t4(2) == 2);
    CHECK(t4(3) == 2);
    CHECK(t4(4) == 3);
    CHECK(t4(5) == 3);

    Eigen::VectorXi i1 = Eigen::VectorXi::Zero(3);
    i1(0) = 1;
    i1(1) = 2;
    i1(2) = 3;

    const Eigen::VectorXi t5 = repeat<int>(i1, -1);
    REQUIRE(t5.size() == 0);
    const Eigen::VectorXi t6 = repeat<int>(i1, 0);
    REQUIRE(t6.size() == 0);
    const Eigen::VectorXi t7 = repeat<int>(i1, 1);
    REQUIRE(t7.size() == 3);
    CHECK(t7(0) == 1);
    CHECK(t7(1) == 2);
    CHECK(t7(2) == 3);
    const Eigen::VectorXi t8 = repeat<int>(i1, 2);
    REQUIRE(t8.size() == 6);
    CHECK(t8(0) == 1);
    CHECK(t8(1) == 1);
    CHECK(t8(2) == 2);
    CHECK(t8(3) == 2);
    CHECK(t8(4) == 3);
    CHECK(t8(5) == 3);

    const Eigen::VectorXd t9 = repeat<double>({1.5, 2.5, 3.5}, -1);
    REQUIRE(t9.size() == 0);
    const Eigen::VectorXd t10 = repeat<double>({1.5, 2.5, 3.5}, 0);
    REQUIRE(t10.size() == 0);
    const Eigen::VectorXd t11 = repeat<double>({1.5, 2.5, 3.5}, 1);
    REQUIRE(t11.size() == 3);
    CHECK(t11(0) == 1.5);
    CHECK(t11(1) == 2.5);
    CHECK(t11(2) == 3.5);
    const Eigen::VectorXd t12 = repeat<double>({1.5, 2.5, 3.5}, 2);
    REQUIRE(t12.size() == 6);
    CHECK(t12(0) == 1.5);
    CHECK(t12(1) == 1.5);
    CHECK(t12(2) == 2.5);
    CHECK(t12(3) == 2.5);
    CHECK(t12(4) == 3.5);
    CHECK(t12(5) == 3.5);

    Eigen::VectorXd d1 = Eigen::VectorXd::Zero(3);
    d1(0) = 1.5;
    d1(1) = 2.5;
    d1(2) = 3.5;

    const Eigen::VectorXd t13 = repeat<double>(d1, -1);
    REQUIRE(t13.size() == 0);
    const Eigen::VectorXd t14 = repeat<double>(d1, 0);
    REQUIRE(t14.size() == 0);
    const Eigen::VectorXd t15 = repeat<double>(d1, 1);
    REQUIRE(t15.size() == 3);
    CHECK(t15(0) == 1.5);
    CHECK(t15(1) == 2.5);
    CHECK(t15(2) == 3.5);
    const Eigen::VectorXd t16 = repeat<double>(d1, 2);
    REQUIRE(t16.size() == 6);
    CHECK(t16(0) == 1.5);
    CHECK(t16(1) == 1.5);
    CHECK(t16(2) == 2.5);
    CHECK(t16(3) == 2.5);
    CHECK(t16(4) == 3.5);
    CHECK(t16(5) == 3.5);
}

TEST_CASE("tile", "[utils]")
{
    using namespace Lielab::utils;

    std::vector<Eigen::VectorXi> li;
    std::vector<Eigen::VectorXd> ld;

    REQUIRE(tile<int>({}, -1).size() == 0);
    REQUIRE(tile<double>({}, -1).size() == 0);
    REQUIRE(tile<int>(Eigen::VectorXi::Zero(0), -1).size() == 0);
    REQUIRE(tile<double>(Eigen::VectorXd::Zero(0), -1).size() == 0);

    REQUIRE(tile<int>({}, 0).size() == 0);
    REQUIRE(tile<double>({}, 0).size() == 0);
    REQUIRE(tile<int>(Eigen::VectorXi::Zero(0), 0).size() == 0);
    REQUIRE(tile<double>(Eigen::VectorXd::Zero(0), 0).size() == 0);

    REQUIRE(tile<int>({}, 1).size() == 0);
    REQUIRE(tile<double>({}, 1).size() == 0);
    REQUIRE(tile<int>(Eigen::VectorXi::Zero(0), 1).size() == 0);
    REQUIRE(tile<double>(Eigen::VectorXd::Zero(0), 1).size() == 0);

    REQUIRE(tile<int>({}, 2).size() == 0);
    REQUIRE(tile<double>({}, 2).size() == 0);
    REQUIRE(tile<int>(Eigen::VectorXi::Zero(0), 2).size() == 0);
    REQUIRE(tile<double>(Eigen::VectorXd::Zero(0), 2).size() == 0);

    const Eigen::VectorXi t1 = tile<int>({1, 2, 3}, -1);
    REQUIRE(t1.size() == 0);
    const Eigen::VectorXi t2 = tile<int>({1, 2, 3}, 0);
    REQUIRE(t2.size() == 0);
    const Eigen::VectorXi t3 = tile<int>({1, 2, 3}, 1);
    REQUIRE(t3.size() == 3);
    CHECK(t3(0) == 1);
    CHECK(t3(1) == 2);
    CHECK(t3(2) == 3);
    const Eigen::VectorXi t4 = tile<int>({1, 2, 3}, 2);
    REQUIRE(t4.size() == 6);
    CHECK(t4(0) == 1);
    CHECK(t4(1) == 2);
    CHECK(t4(2) == 3);
    CHECK(t4(3) == 1);
    CHECK(t4(4) == 2);
    CHECK(t4(5) == 3);

    Eigen::VectorXi i1 = Eigen::VectorXi::Zero(3);
    i1(0) = 1;
    i1(1) = 2;
    i1(2) = 3;

    const Eigen::VectorXi t5 = tile<int>(i1, -1);
    REQUIRE(t5.size() == 0);
    const Eigen::VectorXi t6 = tile<int>(i1, 0);
    REQUIRE(t6.size() == 0);
    const Eigen::VectorXi t7 = tile<int>(i1, 1);
    REQUIRE(t7.size() == 3);
    CHECK(t7(0) == 1);
    CHECK(t7(1) == 2);
    CHECK(t7(2) == 3);
    const Eigen::VectorXi t8 = tile<int>(i1, 2);
    REQUIRE(t8.size() == 6);
    CHECK(t8(0) == 1);
    CHECK(t8(1) == 2);
    CHECK(t8(2) == 3);
    CHECK(t8(3) == 1);
    CHECK(t8(4) == 2);
    CHECK(t8(5) == 3);

    const Eigen::VectorXd t9 = tile<double>({1.5, 2.5, 3.5}, -1);
    REQUIRE(t9.size() == 0);
    const Eigen::VectorXd t10 = tile<double>({1.5, 2.5, 3.5}, 0);
    REQUIRE(t10.size() == 0);
    const Eigen::VectorXd t11 = tile<double>({1.5, 2.5, 3.5}, 1);
    REQUIRE(t11.size() == 3);
    CHECK(t11(0) == 1.5);
    CHECK(t11(1) == 2.5);
    CHECK(t11(2) == 3.5);
    const Eigen::VectorXd t12 = tile<double>({1.5, 2.5, 3.5}, 2);
    REQUIRE(t12.size() == 6);
    CHECK(t12(0) == 1.5);
    CHECK(t12(1) == 2.5);
    CHECK(t12(2) == 3.5);
    CHECK(t12(3) == 1.5);
    CHECK(t12(4) == 2.5);
    CHECK(t12(5) == 3.5);

    Eigen::VectorXd d1 = Eigen::VectorXd::Zero(3);
    d1(0) = 1.5;
    d1(1) = 2.5;
    d1(2) = 3.5;

    const Eigen::VectorXd t13 = tile<double>(d1, -1);
    REQUIRE(t13.size() == 0);
    const Eigen::VectorXd t14 = tile<double>(d1, 0);
    REQUIRE(t14.size() == 0);
    const Eigen::VectorXd t15 = tile<double>(d1, 1);
    REQUIRE(t15.size() == 3);
    CHECK(t15(0) == 1.5);
    CHECK(t15(1) == 2.5);
    CHECK(t15(2) == 3.5);
    const Eigen::VectorXd t16 = tile<double>(d1, 2);
    REQUIRE(t16.size() == 6);
    CHECK(t16(0) == 1.5);
    CHECK(t16(1) == 2.5);
    CHECK(t16(2) == 3.5);
    CHECK(t16(3) == 1.5);
    CHECK(t16(4) == 2.5);
    CHECK(t16(5) == 3.5);
}

TEST_CASE("linspace", "[utils]")
{
    using namespace Lielab::utils;

    const Eigen::VectorXd t1 = linspace<double>(1.0, 5.0, -1);
    REQUIRE(t1.size() == 0);

    const Eigen::VectorXd t2 = linspace<double>(1.0, 5.0, 0);
    REQUIRE(t2.size() == 0);

    const Eigen::VectorXd t3 = linspace<double>(1.0, 5.0, 1);
    REQUIRE(t3.size() == 1);
    CHECK(t3(0) == 1.0);
    
    const Eigen::VectorXd t4 = linspace<double>(1.0, 5.0, 2);
    REQUIRE(t4.size() == 2);
    CHECK(t4(0) == 1.0);
    CHECK(t4(1) == 5.0);

    const Eigen::VectorXd t5 = linspace<double>(1.0, 5.0, 3);
    REQUIRE(t5.size() == 3);
    CHECK(t5(0) == 1.0);
    CHECK(t5(1) == 3.0);
    CHECK(t5(2) == 5.0);

    const Eigen::VectorXd t6 = linspace<double>(1.0, 5.0, 5);
    REQUIRE(t6.size() == 5);
    CHECK(t6(0) == 1.0);
    CHECK(t6(1) == 2.0);
    CHECK(t6(2) == 3.0);
    CHECK(t6(3) == 4.0);
    CHECK(t6(4) == 5.0);
}

TEST_CASE("column_stack", "[utils]")
{
    using namespace Lielab::utils;

    std::vector<Eigen::VectorXi> li;
    std::vector<Eigen::VectorXd> ld;

    REQUIRE(column_stack<int>({}).size() == 0);
    REQUIRE(column_stack<double>({}).size() == 0);
    REQUIRE(column_stack<int>(li).size() == 0);
    REQUIRE(column_stack<double>(ld).size() == 0);

    Eigen::VectorXi i1 = Eigen::VectorXi::Zero(3);
    i1(0) = 1;
    i1(1) = 2;
    i1(2) = 3;

    Eigen::VectorXi i2 = Eigen::VectorXi::Zero(2);
    i2(0) = 4;
    i2(1) = 5;

    const Eigen::MatrixXi t1 = column_stack({i1});
    REQUIRE(t1.rows() == 3);
    REQUIRE(t1.cols() == 1);
    CHECK(t1(0, 0) == 1);
    CHECK(t1(1, 0) == 2);
    CHECK(t1(2, 0) == 3);

    const Eigen::MatrixXi t2 = column_stack({i1, i2});
    REQUIRE(t2.rows() == 2);
    REQUIRE(t2.cols() == 2);
    CHECK(t2(0, 0) == 1);
    CHECK(t2(1, 0) == 2);
    CHECK(t2(0, 1) == 4);
    CHECK(t2(1, 1) == 5);

    const Eigen::MatrixXi t3 = column_stack({i1, i2, i1});
    REQUIRE(t3.rows() == 2);
    REQUIRE(t3.cols() == 3);
    CHECK(t3(0, 0) == 1);
    CHECK(t3(1, 0) == 2);
    CHECK(t3(0, 1) == 4);
    CHECK(t3(1, 1) == 5);
    CHECK(t3(0, 2) == 1);
    CHECK(t3(1, 2) == 2);

    li.push_back(i1);
    const Eigen::MatrixXi t4 = column_stack(li);
    REQUIRE(t4.rows() == 3);
    REQUIRE(t4.cols() == 1);
    CHECK(t4(0, 0) == 1);
    CHECK(t4(1, 0) == 2);
    CHECK(t4(2, 0) == 3);

    li.push_back(i2);
    const Eigen::MatrixXi t5 = column_stack(li);
    REQUIRE(t5.rows() == 2);
    REQUIRE(t5.cols() == 2);
    CHECK(t5(0, 0) == 1);
    CHECK(t5(1, 0) == 2);
    CHECK(t5(0, 1) == 4);
    CHECK(t5(1, 1) == 5);

    li.push_back(i1);
    const Eigen::MatrixXi t6 = column_stack(li);
    REQUIRE(t6.rows() == 2);
    REQUIRE(t6.cols() == 3);
    CHECK(t6(0, 0) == 1);
    CHECK(t6(1, 0) == 2);
    CHECK(t6(0, 1) == 4);
    CHECK(t6(1, 1) == 5);
    CHECK(t6(0, 2) == 1);
    CHECK(t6(1, 2) == 2);

    Eigen::VectorXd d1 = Eigen::VectorXd::Zero(3);
    d1(0) = 1.5;
    d1(1) = 2.5;
    d1(2) = 3.5;

    Eigen::VectorXd d2 = Eigen::VectorXd::Zero(2);
    d2(0) = 4.5;
    d2(1) = 5.5;

    const Eigen::MatrixXd t7 = column_stack({d1});
    REQUIRE(t7.rows() == 3);
    REQUIRE(t7.cols() == 1);
    CHECK(t7(0, 0) == 1.5);
    CHECK(t7(1, 0) == 2.5);
    CHECK(t7(2, 0) == 3.5);

    const Eigen::MatrixXd t8 = column_stack({d1, d2});
    REQUIRE(t8.rows() == 2);
    REQUIRE(t8.cols() == 2);
    CHECK(t8(0, 0) == 1.5);
    CHECK(t8(1, 0) == 2.5);
    CHECK(t8(0, 1) == 4.5);
    CHECK(t8(1, 1) == 5.5);

    const Eigen::MatrixXd t9 = column_stack({d1, d2, d1});
    REQUIRE(t9.rows() == 2);
    REQUIRE(t9.cols() == 3);
    CHECK(t9(0, 0) == 1.5);
    CHECK(t9(1, 0) == 2.5);
    CHECK(t9(0, 1) == 4.5);
    CHECK(t9(1, 1) == 5.5);
    CHECK(t9(0, 2) == 1.5);
    CHECK(t9(1, 2) == 2.5);

    ld.push_back(d1);
    const Eigen::MatrixXd t10 = column_stack(ld);
    REQUIRE(t10.rows() == 3);
    REQUIRE(t10.cols() == 1);
    CHECK(t10(0, 0) == 1.5);
    CHECK(t10(1, 0) == 2.5);
    CHECK(t10(2, 0) == 3.5);

    ld.push_back(d2);
    const Eigen::MatrixXd t11 = column_stack(ld);
    REQUIRE(t11.rows() == 2);
    REQUIRE(t11.cols() == 2);
    CHECK(t11(0, 0) == 1.5);
    CHECK(t11(1, 0) == 2.5);
    CHECK(t11(0, 1) == 4.5);
    CHECK(t11(1, 1) == 5.5);

    ld.push_back(d1);
    const Eigen::MatrixXd t12 = column_stack(ld);
    REQUIRE(t12.rows() == 2);
    REQUIRE(t12.cols() == 3);
    CHECK(t12(0, 0) == 1.5);
    CHECK(t12(1, 0) == 2.5);
    CHECK(t12(0, 1) == 4.5);
    CHECK(t12(1, 1) == 5.5);
    CHECK(t12(0, 2) == 1.5);
    CHECK(t12(1, 2) == 2.5);
}
