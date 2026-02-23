#include <Lielab.hpp>

#include <catch2/catch_all.hpp>

#include <cmath>
#include <numbers>

TEST_CASE("Grassmannian to_string", "[domain]")
{
    using Lielab::domain::Grassmannian;

    const Grassmannian xblank = Grassmannian();
    CHECK(xblank.to_string() == "Grassmannian(0, 0, R)");

    const Grassmannian x0 = Grassmannian(0, 0);
    CHECK(x0.to_string() == "Grassmannian(0, 0, R)");
    const Grassmannian x1 = Grassmannian(1, 1);
    CHECK(x1.to_string() == "Grassmannian(1, 1, R)");
    const Grassmannian x10 = Grassmannian(5, 10);
    CHECK(x10.to_string() == "Grassmannian(5, 10, R)");
}

TEST_CASE("Grassmannian main_initializer", "[domain]")
{
    using Lielab::domain::Grassmannian;

    const Grassmannian xblank = Grassmannian();
    CHECK(xblank.get_dimension() == 0);

    const Grassmannian x0 = Grassmannian(0, 0);
    CHECK(x0.get_dimension() == 0);
    const Grassmannian x1 = Grassmannian(1, 1);
    CHECK(x1.get_dimension() == 1);
    const Grassmannian x10 = Grassmannian(5, 10);
    CHECK(x10.get_dimension() == 5);

    CHECK_THROWS(Grassmannian(6, 5));
}

TEST_CASE("Grassmannian get_dimension", "[domain]")
{
    using Lielab::domain::Grassmannian;

    Grassmannian zero(0, 0), one(1, 1), two(2, 2), three(3, 3), four(4, 4), five(5, 5), six(6, 6), seven(7, 7), eight(8, 8);

    CHECK(zero.get_dimension() == 0);
    CHECK(one.get_dimension() == 1);
    CHECK(two.get_dimension() == 2);
    CHECK(three.get_dimension() == 3);
    CHECK(four.get_dimension() == 4);
    CHECK(five.get_dimension() == 5);
    CHECK(six.get_dimension() == 6);
    CHECK(seven.get_dimension() == 7);
    CHECK(eight.get_dimension() == 8);
}

TEST_CASE("Grassmannian serialize/unserialize", "[domain]")
{
    /*!
    * Tests the serialize/unserialize operation.
    */

    using Lielab::domain::Grassmannian;

    Grassmannian x0 = Grassmannian(0, 0);
    x0.unserialize({});
    Eigen::VectorXd x0bar = x0.serialize();

    CHECK(x0bar.size() == 0);

    Grassmannian x1 = Grassmannian(1, 2);
    x1.unserialize({1.0, 2.0, 3.0});
    Eigen::VectorXd x1bar = x1.serialize();

    REQUIRE(x1bar.size() == 2);
    CHECK(x1bar(0) == 1.0);
    CHECK(x1bar(1) == 2.0);

    x1.unserialize({4.0, 5.0});
    x1bar = x1.serialize();

    REQUIRE(x1bar.size() == 2);
    CHECK(x1bar(0) == 4.0);
    CHECK(x1bar(1) == 5.0);

    x1.unserialize({6.0});
    x1bar = x1.serialize();

    REQUIRE(x1bar.size() == 2);
    CHECK(x1bar(0) == 6.0);
    CHECK(x1bar(1) == 5.0);

    Grassmannian x2 = Grassmannian(2, 4);
    x2.unserialize({1.0, 2.0, 3.0});
    Eigen::VectorXd x2bar = x2.serialize();

    REQUIRE(x2bar.size() == 4);
    CHECK(x2bar(0) == 1.0);
    CHECK(x2bar(1) == 2.0);
    CHECK(x2bar(2) == 3.0);
    CHECK(x2bar(3) == 0.0);

    x2.unserialize({4.0, 5.0, 6.0, 7.0});
    x2bar = x2.serialize();

    REQUIRE(x2bar.size() == 4);
    CHECK(x2bar(0) == 4.0);
    CHECK(x2bar(1) == 5.0);
    CHECK(x2bar(2) == 6.0);
    CHECK(x2bar(3) == 7.0);

    x2.unserialize({8.0, 9.0, 10.0, 11.0, 12.0});
    x2bar = x2.serialize();

    REQUIRE(x2bar.size() == 4);
    CHECK(x2bar(0) == 8.0);
    CHECK(x2bar(1) == 9.0);
    CHECK(x2bar(2) == 10.0);
    CHECK(x2bar(3) == 11.0);
}

TEST_CASE("Grassmannian project_point", "[domain]")
{
    using Lielab::domain::Grassmannian;
    using Lielab::utils::to_VectorXd;
    
    const double pi = std::numbers::pi_v<double>;

    Eigen::MatrixXd axes(2, 1);
    axes(0,0) = std::cos(1.0/6.0*pi);
    axes(1,0) = std::sin(1.0/6.0*pi);

    const Grassmannian Gr12 = Grassmannian(to_VectorXd({0.0, 0.0}), axes);

    const Eigen::VectorXd vals1 = to_VectorXd({std::cos(1.0/6.0*pi), std::sin(1.0/6.0*pi)});
    const Eigen::VectorXd proj1 = Gr12.project_point(vals1);

    REQUIRE(proj1.size() == 2);
    CHECK_THAT(proj1(0), Catch::Matchers::WithinAbs(std::cos(1.0/6.0*pi), 1e-15));
    CHECK_THAT(proj1(1), Catch::Matchers::WithinAbs(std::sin(1.0/6.0*pi), 1e-15));

    const Eigen::VectorXd vals2 = to_VectorXd({-std::sin(1.0/6.0*pi), std::cos(1.0/6.0*pi)});
    const Eigen::VectorXd proj2 = Gr12.project_point(vals2);

    REQUIRE(proj2.size() == 2);
    CHECK_THAT(proj2(0), Catch::Matchers::WithinAbs(0.0, 1e-15));
    CHECK_THAT(proj2(1), Catch::Matchers::WithinAbs(0.0, 1e-15));

    const Eigen::VectorXd vals3 = to_VectorXd({std::cos(1.0/6.0*pi), std::sin(1.0/6.0*pi) - 4.0});
    const Eigen::VectorXd proj3 = Gr12.project_point(vals3);

    REQUIRE(proj3.size() == 2);
    CHECK_THAT(proj3(0), Catch::Matchers::WithinAbs(-std::cos(1.0/6.0*pi), 1e-15));
    CHECK_THAT(proj3(1), Catch::Matchers::WithinAbs(-std::sin(1.0/6.0*pi), 1e-15));
}

TEST_CASE("Grassmannian project_vector_onto_tangent_space", "[domain]")
{
    using Lielab::domain::Grassmannian;
    using Lielab::utils::to_VectorXd;
    
    const double pi = std::numbers::pi_v<double>;
    const Eigen::VectorXd vector = to_VectorXd({std::sqrt(2.0)/2.0, std::sqrt(2.0)/2.0});

    Eigen::MatrixXd axes1(2,1);
    axes1(0,0) = 1.0;
    axes1(1,0) = 0.0;
    const Grassmannian Gr121 = Grassmannian::project(to_VectorXd({0.0, 0.0}), axes1);
    const Eigen::VectorXd proj1 = Gr121.project_vector_onto_tangent_space(vector);

    REQUIRE(proj1.size() == 2);
    CHECK_THAT(proj1(0), Catch::Matchers::WithinAbs(std::sqrt(2.0)/2.0, 1e-15));
    CHECK_THAT(proj1(1), Catch::Matchers::WithinAbs(0.0, 1e-15));

    Eigen::MatrixXd axes2(2,1);
    axes2(0,0) = 0.0;
    axes2(1,0) = 1.0;
    const Grassmannian Gr122 = Grassmannian::project(to_VectorXd({0.0, 0.0}), axes2);
    const Eigen::VectorXd proj2 = Gr122.project_vector_onto_tangent_space(vector);

    REQUIRE(proj2.size() == 2);
    CHECK_THAT(proj2(0), Catch::Matchers::WithinAbs(0.0, 1e-15));
    CHECK_THAT(proj2(1), Catch::Matchers::WithinAbs(std::sqrt(2.0)/2.0, 1e-15));

    Eigen::MatrixXd axes3(2,1);
    axes3(0,0) = 1.0;
    axes3(1,0) = 1.0;
    const Grassmannian Gr123 = Grassmannian::project(to_VectorXd({0.0, 0.0}), axes3);
    const Eigen::VectorXd proj3 = Gr123.project_vector_onto_tangent_space(vector);

    REQUIRE(proj3.size() == 2);
    CHECK_THAT(proj3(0), Catch::Matchers::WithinAbs(std::sqrt(2.0)/2.0, 1e-15));
    CHECK_THAT(proj3(1), Catch::Matchers::WithinAbs(std::sqrt(2.0)/2.0, 1e-15));
}

TEST_CASE("Grassmannian project_vector_onto_normal_space", "[domain]")
{
    using Lielab::domain::Grassmannian;
    using Lielab::utils::to_VectorXd;
    
    const double pi = std::numbers::pi_v<double>;
    const Eigen::VectorXd vector = to_VectorXd({std::sqrt(2.0)/2.0, std::sqrt(2.0)/2.0});

    Eigen::MatrixXd axes1(2,1);
    axes1(0,0) = 1.0;
    axes1(1,0) = 0.0;
    const Grassmannian Gr121 = Grassmannian::project(to_VectorXd({0.0, 0.0}), axes1);
    const Eigen::VectorXd proj1 = Gr121.project_vector_onto_normal_space(vector);

    REQUIRE(proj1.size() == 2);
    CHECK_THAT(proj1(0), Catch::Matchers::WithinAbs(0.0, 1e-15));
    CHECK_THAT(proj1(1), Catch::Matchers::WithinAbs(std::sqrt(2.0)/2.0, 1e-15));

    Eigen::MatrixXd axes2(2,1);
    axes2(0,0) = 0.0;
    axes2(1,0) = 1.0;
    const Grassmannian Gr122 = Grassmannian::project(to_VectorXd({0.0, 0.0}), axes2);
    const Eigen::VectorXd proj2 = Gr122.project_vector_onto_normal_space(vector);

    REQUIRE(proj2.size() == 2);
    CHECK_THAT(proj2(0), Catch::Matchers::WithinAbs(std::sqrt(2.0)/2.0, 1e-15));
    CHECK_THAT(proj2(1), Catch::Matchers::WithinAbs(0.0, 1e-15));

    Eigen::MatrixXd axes3(2,1);
    axes3(0,0) = 1.0;
    axes3(1,0) = 1.0;
    const Grassmannian Gr123 = Grassmannian::project(to_VectorXd({0.0, 0.0}), axes3);
    const Eigen::VectorXd proj3 = Gr123.project_vector_onto_normal_space(vector);

    REQUIRE(proj3.size() == 2);
    CHECK_THAT(proj3(0), Catch::Matchers::WithinAbs(0.0, 1e-15));
    CHECK_THAT(proj3(1), Catch::Matchers::WithinAbs(0.0, 1e-15));
}

TEST_CASE("Grassmannian axes_intersection", "[domain]")
{
    using Lielab::domain::Grassmannian;
    using Lielab::utils::to_VectorXd;

    Eigen::MatrixXd axes1(2,1);
    axes1(0,0) = 1.0;
    axes1(1,0) = 0.0;
    const Grassmannian Gr121 = Grassmannian::project(to_VectorXd({0.0, 0.0}), axes1);
    const Eigen::VectorXd intersection1 = Gr121.axes_intersection(to_VectorXd({0.0, 1.0}), to_VectorXd({0.0, -1.0}));

    REQUIRE(intersection1.size() == 2);
    CHECK_THAT(intersection1(0), Catch::Matchers::WithinAbs(0.0, 1e-15));
    CHECK_THAT(intersection1(1), Catch::Matchers::WithinAbs(0.0, 1e-15));

    const Eigen::VectorXd intersection2 = Gr121.axes_intersection(to_VectorXd({0.0, 1.0}), to_VectorXd({1.0, 0.0}));

    REQUIRE(intersection2.size() == 2);
    CHECK(std::isnan(intersection2(0)));
    CHECK(std::isnan(intersection2(1)));

    const Eigen::VectorXd intersection3 = Gr121.axes_intersection(to_VectorXd({0.0, 1.0}), to_VectorXd({1.0, -1.0}));

    REQUIRE(intersection3.size() == 2);
    CHECK_THAT(intersection3(0), Catch::Matchers::WithinAbs(1.0, 1e-15));
    CHECK_THAT(intersection3(1), Catch::Matchers::WithinAbs(0.0, 1e-15));

    Eigen::MatrixXd axes2(2,1);
    axes2(0,0) = 1.0;
    axes2(1,0) = 1.0;
    const Grassmannian Gr122 = Grassmannian::project(to_VectorXd({0.0, 0.0}), axes2);
    const Eigen::VectorXd intersection4 = Gr122.axes_intersection(to_VectorXd({0.0, 1.0}), to_VectorXd({0.0, -1.0}));

    REQUIRE(intersection4.size() == 2);
    CHECK_THAT(intersection4(0), Catch::Matchers::WithinAbs(0.0, 1e-15));
    CHECK_THAT(intersection4(1), Catch::Matchers::WithinAbs(0.0, 1e-15));

    const Eigen::VectorXd intersection5 = Gr122.axes_intersection(to_VectorXd({0.0, 1.0}), to_VectorXd({1.0, 0.0}));

    REQUIRE(intersection5.size() == 2);
    CHECK_THAT(intersection5(0), Catch::Matchers::WithinAbs(1.0, 1e-15));
    CHECK_THAT(intersection5(1), Catch::Matchers::WithinAbs(1.0, 1e-15));

    const Eigen::VectorXd intersection6 = Gr122.axes_intersection(to_VectorXd({0.0, 1.0}), to_VectorXd({1.0, -1.0}));

    REQUIRE(intersection6.size() == 2);
    CHECK_THAT(intersection6(0), Catch::Matchers::WithinAbs(0.5, 1e-15));
    CHECK_THAT(intersection6(1), Catch::Matchers::WithinAbs(0.5, 1e-15));
}
