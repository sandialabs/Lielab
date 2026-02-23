#include <Lielab/domain/liegroups/SO.hpp>
#include <catch2/catch_all.hpp>

#include "../../test_utils.hpp"

#include <cmath>
#include <functional>
#include <numbers>

// using std::numbers::pi;
constexpr double PI = std::numbers::pi_v<double>;
const Lielab::domain::SO DCMId(3);
const Lielab::domain::SO DCMrotx(Lielab::functions::exp(PI/2.0*Lielab::domain::so::basis(0,3)));
const Lielab::domain::SO DCMroty(Lielab::functions::exp(PI/2.0*Lielab::domain::so::basis(1,3)));
const Lielab::domain::SO DCMrotz(Lielab::functions::exp(PI/2.0*Lielab::domain::so::basis(2,3)));
const double some_angle = PI/2.0*5.0/7.0;

TEST_CASE("SO", "[domain]")
{
    /*!
    * Tests SO against well-known identities.
    */

    for (int shape = 2; shape <= 4; shape++)
    {
        const int D = Lielab::domain::so::basis(0, shape).get_dimension();

        // Construct the SO elements
        std::vector<Lielab::domain::SO> elements;
        for (int ii = 0; ii < D; ii++)
        {
            elements.push_back(Lielab::functions::exp(Lielab::domain::so::basis(ii, shape)));
        }

        const Lielab::domain::SO identity(shape);

        is_group<Lielab::domain::SO>(elements, identity);
    }
}

TEST_CASE("from_SU2", "[domain]")
{
    /*!
    * Tests the from_SU2 function
    */

    using Lielab::domain::CompositeGroup;
    using Lielab::domain::SO;
    using Lielab::domain::SU;
    using Lielab::domain::su;
    using Lielab::functions::exp;
    using Lielab::testing::check_almost_equal_nulp;

    // Build 90 degree rotations in x, y, and z
    const su u = su::from_vector({0.0, 0.0, 0.5});
    const su v = su::from_vector({0.0, 0.5, 0.0});
    const su w = su::from_vector({0.5, 0.0, 0.0});
    const SU qx = exp(PI/2.0*u);
    const SU qy = exp(PI/2.0*v);
    const SU qz = exp(PI/2.0*w);

    // Test 90 degree x rotation
    const SO rx = SO::from_SU2(qx);

    CHECK(check_almost_equal_nulp(CompositeGroup{rx}, CompositeGroup{DCMrotx}, 4, true));

    // Test 90 degree y rotation
    const SO ry = SO::from_SU2(qy);

    CHECK(check_almost_equal_nulp(CompositeGroup{ry}, CompositeGroup{DCMroty}, 4, true));

    // Test 90 degree z rotation
    const SO rz = SO::from_SU2(qz);

    CHECK(check_almost_equal_nulp(CompositeGroup{rz}, CompositeGroup{DCMrotz}, 4, true));
}

TEST_CASE("from_eulerangles_body123", "[domain]")
{
    /*!
    * Tests the eanglebody123 function.
    */

    using Lielab::domain::SO;
    using Lielab::domain::CompositeGroup;
    using Lielab::testing::check_almost_equal_nulp;
    
    SO dcm(3);

    // Identity
    dcm = Lielab::domain::SO::from_eulerangles_body123(0.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMId}, CompositeGroup{dcm}, 0));

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_body123(PI/2.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotx}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_body123(0.0, PI/2.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMroty}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_body123(0.0, 0.0, PI/2.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotz}, CompositeGroup{dcm}, 4, true));

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_body123(some_angle, some_angle, some_angle);
    CHECK_THAT(dcm.get_matrix().determinant(), Catch::Matchers::WithinULP(1.0, 4));
    CHECK_THAT((dcm.get_matrix()*dcm.get_matrix().transpose()).trace(), Catch::Matchers::WithinULP(3.0, 4));
}

TEST_CASE("from_eulerangles_body231", "[domain]")
{
    /*!
    * Tests the eanglebody231 function.
    */

    using Lielab::domain::SO;
    using Lielab::domain::CompositeGroup;
    using Lielab::testing::check_almost_equal_nulp;
    
    SO dcm(3);

    // Identity
    dcm = Lielab::domain::SO::from_eulerangles_body231(0.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMId}, CompositeGroup{dcm}, 0));

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_body231(0.0, 0.0, PI/2.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotx}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_body231(PI/2.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMroty}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_body231(0.0, PI/2.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotz}, CompositeGroup{dcm}, 4, true));

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_body231(some_angle, some_angle, some_angle);
    CHECK_THAT(dcm.get_matrix().determinant(), Catch::Matchers::WithinULP(1.0, 4));
    CHECK_THAT((dcm.get_matrix()*dcm.get_matrix().transpose()).trace(), Catch::Matchers::WithinULP(3.0, 4));
}

TEST_CASE("from_eulerangles_body312", "[domain]")
{
    /*!
    * Tests the eanglebody312 function.
    */

    using Lielab::domain::SO;
    using Lielab::domain::CompositeGroup;
    using Lielab::testing::check_almost_equal_nulp;
    
    SO dcm(3);

    // Identity
    dcm = Lielab::domain::SO::from_eulerangles_body312(0.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMId}, CompositeGroup{dcm}, 0));

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_body312(0.0, PI/2.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotx}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_body312(0.0, 0.0, PI/2.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMroty}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_body312(PI/2.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotz}, CompositeGroup{dcm}, 4, true));

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_body312(some_angle, some_angle, some_angle);
    CHECK_THAT(dcm.get_matrix().determinant(), Catch::Matchers::WithinULP(1.0, 4));
    CHECK_THAT((dcm.get_matrix()*dcm.get_matrix().transpose()).trace(), Catch::Matchers::WithinULP(3.0, 4));
}

TEST_CASE("from_eulerangles_body132", "[domain]")
{
    /*!
    * Tests the eanglebody132 function.
    */

    using Lielab::domain::SO;
    using Lielab::domain::CompositeGroup;
    using Lielab::testing::check_almost_equal_nulp;
    
    SO dcm(3);

    // Identity
    dcm = Lielab::domain::SO::from_eulerangles_body132(0.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMId}, CompositeGroup{dcm}, 0));

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_body132(PI/2.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotx}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_body132(0.0, 0.0, PI/2.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMroty}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_body132(0.0, PI/2.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotz}, CompositeGroup{dcm}, 4, true));

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_body132(some_angle, some_angle, some_angle);
    CHECK_THAT(dcm.get_matrix().determinant(), Catch::Matchers::WithinULP(1.0, 4));
    CHECK_THAT((dcm.get_matrix()*dcm.get_matrix().transpose()).trace(), Catch::Matchers::WithinULP(3.0, 4));
}

TEST_CASE("from_eulerangles_body213", "[domain]")
{
    /*!
    * Tests the eanglebody213 function.
    */

    using Lielab::domain::SO;
    using Lielab::domain::CompositeGroup;
    using Lielab::testing::check_almost_equal_nulp;
    
    SO dcm(3);

    // Identity
    dcm = Lielab::domain::SO::from_eulerangles_body213(0.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMId}, CompositeGroup{dcm}, 0));

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_body213(0.0, PI/2.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotx}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_body213(PI/2.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMroty}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_body213(0.0, 0.0, PI/2.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotz}, CompositeGroup{dcm}, 4, true));

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_body213(some_angle, some_angle, some_angle);
    CHECK_THAT(dcm.get_matrix().determinant(), Catch::Matchers::WithinULP(1.0, 4));
    CHECK_THAT((dcm.get_matrix()*dcm.get_matrix().transpose()).trace(), Catch::Matchers::WithinULP(3.0, 4));
}

TEST_CASE("from_eulerangles_body321", "[domain]")
{
    /*!
    * Tests the eanglebody321 function.
    */

    using Lielab::domain::SO;
    using Lielab::domain::CompositeGroup;
    using Lielab::testing::check_almost_equal_nulp;
    
    SO dcm(3);

    // Identity
    dcm = Lielab::domain::SO::from_eulerangles_body321(0.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMId}, CompositeGroup{dcm}, 0));

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_body321(0.0, 0.0, PI/2.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotx}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_body321(0.0, PI/2.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMroty}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_body321(PI/2.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotz}, CompositeGroup{dcm}, 4, true));

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_body321(some_angle, some_angle, some_angle);
    CHECK_THAT(dcm.get_matrix().determinant(), Catch::Matchers::WithinULP(1.0, 4));
    CHECK_THAT((dcm.get_matrix()*dcm.get_matrix().transpose()).trace(), Catch::Matchers::WithinULP(3.0, 4));
}

TEST_CASE("from_eulerangles_body121", "[domain]")
{
    /*!
    * Tests the eanglebody121 function.
    */

    using Lielab::domain::SO;
    using Lielab::domain::CompositeGroup;
    using Lielab::testing::check_almost_equal_nulp;
    
    SO dcm(3);

    // Identity
    dcm = Lielab::domain::SO::from_eulerangles_body121(0.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMId}, CompositeGroup{dcm}, 0));

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_body121(PI/2.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotx}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_body121(0.0, 0.0, PI/2.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotx}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_body121(0.0, PI/2.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMroty}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_body121(PI/2.0, PI/2.0, -PI/2.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotz}, CompositeGroup{dcm}, 4, true));

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_body121(some_angle, some_angle, some_angle);
    CHECK_THAT(dcm.get_matrix().determinant(), Catch::Matchers::WithinULP(1.0, 4));
    CHECK_THAT((dcm.get_matrix()*dcm.get_matrix().transpose()).trace(), Catch::Matchers::WithinULP(3.0, 4));
}

TEST_CASE("from_eulerangles_body131", "[domain]")
{
    /*!
    * Tests the eanglebody131 function.
    */

    using Lielab::domain::SO;
    using Lielab::domain::CompositeGroup;
    using Lielab::testing::check_almost_equal_nulp;
    
    SO dcm(3);

    // Identity
    dcm = Lielab::domain::SO::from_eulerangles_body131(0.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMId}, CompositeGroup{dcm}, 0));

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_body131(PI/2.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotx}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_body131(0.0, 0.0, PI/2.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotx}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_body131(PI/2.0, -PI/2.0, -PI/2.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMroty}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_body131(0.0, PI/2.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotz}, CompositeGroup{dcm}, 4, true));

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_body131(some_angle, some_angle, some_angle);
    CHECK_THAT(dcm.get_matrix().determinant(), Catch::Matchers::WithinULP(1.0, 4));
    CHECK_THAT((dcm.get_matrix()*dcm.get_matrix().transpose()).trace(), Catch::Matchers::WithinULP(3.0, 4));
}

TEST_CASE("from_eulerangles_body212", "[domain]")
{
    /*!
    * Tests the eanglebody212 function.
    */

    using Lielab::domain::SO;
    using Lielab::domain::CompositeGroup;
    using Lielab::testing::check_almost_equal_nulp;
    
    SO dcm(3);

    // Identity
    dcm = Lielab::domain::SO::from_eulerangles_body212(0.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMId}, CompositeGroup{dcm}, 0));

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_body212(0.0, PI/2.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotx}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_body212(PI/2.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMroty}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_body212(0.0, 0.0, PI/2.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMroty}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_body212(PI/2.0, -PI/2.0, -PI/2.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotz}, CompositeGroup{dcm}, 4, true));

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_body212(some_angle, some_angle, some_angle);
    CHECK_THAT(dcm.get_matrix().determinant(), Catch::Matchers::WithinULP(1.0, 4));
    CHECK_THAT((dcm.get_matrix()*dcm.get_matrix().transpose()).trace(), Catch::Matchers::WithinULP(3.0, 4));
}

TEST_CASE("from_eulerangles_body232", "[domain]")
{
    /*!
    * Tests the eanglebody232 function.
    */

    using Lielab::domain::SO;
    using Lielab::domain::CompositeGroup;
    using Lielab::testing::check_almost_equal_nulp;
    
    SO dcm(3);

    // Identity
    dcm = Lielab::domain::SO::from_eulerangles_body232(0.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMId}, CompositeGroup{dcm}, 0));

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_body232(PI/2.0, PI/2.0, -PI/2.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotx}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_body232(PI/2.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMroty}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_body232(0.0, 0.0, PI/2.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMroty}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_body232(0.0, PI/2.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotz}, CompositeGroup{dcm}, 4, true));

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_body232(some_angle, some_angle, some_angle);
    CHECK_THAT(dcm.get_matrix().determinant(), Catch::Matchers::WithinULP(1.0, 4));
    CHECK_THAT((dcm.get_matrix()*dcm.get_matrix().transpose()).trace(), Catch::Matchers::WithinULP(3.0, 4));
}

TEST_CASE("from_eulerangles_body313", "[domain]")
{
    /*!
    * Tests the eanglebody313 function.
    */

    using Lielab::domain::SO;
    using Lielab::domain::CompositeGroup;
    using Lielab::testing::check_almost_equal_nulp;
    
    SO dcm(3);

    // Identity
    dcm = Lielab::domain::SO::from_eulerangles_body313(0.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMId}, CompositeGroup{dcm}, 0));

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_body313(0.0, PI/2.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotx}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_body313(-PI/2.0, -PI/2.0, PI/2.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMroty}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_body313(PI/2.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotz}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_body313(0.0, 0.0, PI/2.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotz}, CompositeGroup{dcm}, 4, true));

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_body313(some_angle, some_angle, some_angle);
    CHECK_THAT(dcm.get_matrix().determinant(), Catch::Matchers::WithinULP(1.0, 4));
    CHECK_THAT((dcm.get_matrix()*dcm.get_matrix().transpose()).trace(), Catch::Matchers::WithinULP(3.0, 4));
}

TEST_CASE("from_eulerangles_body323", "[domain]")
{
    /*!
    * Tests the eanglebody323 function.
    */

    using Lielab::domain::SO;
    using Lielab::domain::CompositeGroup;
    using Lielab::testing::check_almost_equal_nulp;
    
    SO dcm(3);

    // Identity
    dcm = Lielab::domain::SO::from_eulerangles_body323(0.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMId}, CompositeGroup{dcm}, 0));

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_body323(PI/2.0, -PI/2.0, -PI/2.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotx}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_body323(0.0, PI/2.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMroty}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_body323(PI/2.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotz}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_body323(0.0, 0.0, PI/2.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotz}, CompositeGroup{dcm}, 4, true));

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_body323(some_angle, some_angle, some_angle);
    CHECK_THAT(dcm.get_matrix().determinant(), Catch::Matchers::WithinULP(1.0, 4));
    CHECK_THAT((dcm.get_matrix()*dcm.get_matrix().transpose()).trace(), Catch::Matchers::WithinULP(3.0, 4));
}

TEST_CASE("from_eulerangles_space123", "[domain]")
{
    /*!
    * Tests the eanglespace123 function.
    */

    using Lielab::domain::SO;
    using Lielab::domain::CompositeGroup;
    using Lielab::testing::check_almost_equal_nulp;
    
    SO dcm(3);

    // Identity
    dcm = Lielab::domain::SO::from_eulerangles_space123(0.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMId}, CompositeGroup{dcm}, 0));

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_space123(PI/2.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotx}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_space123(0.0, PI/2.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMroty}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_space123(0.0, 0.0, PI/2.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotz}, CompositeGroup{dcm}, 4, true));

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_space123(some_angle, some_angle, some_angle);
    CHECK_THAT(dcm.get_matrix().determinant(), Catch::Matchers::WithinULP(1.0, 4));
    CHECK_THAT((dcm.get_matrix()*dcm.get_matrix().transpose()).trace(), Catch::Matchers::WithinULP(3.0, 4));
}

TEST_CASE("from_eulerangles_space231", "[domain]")
{
    /*!
    * Tests the eanglespace231 function.
    */

    using Lielab::domain::SO;
    using Lielab::domain::CompositeGroup;
    using Lielab::testing::check_almost_equal_nulp;
    
    SO dcm(3);

    // Identity
    dcm = Lielab::domain::SO::from_eulerangles_space231(0.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMId}, CompositeGroup{dcm}, 0));

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_space231(0.0, 0.0, PI/2.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotx}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_space231(PI/2.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMroty}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_space231(0.0, PI/2.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotz}, CompositeGroup{dcm}, 4, true));

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_space231(some_angle, some_angle, some_angle);
    CHECK_THAT(dcm.get_matrix().determinant(), Catch::Matchers::WithinULP(1.0, 4));
    CHECK_THAT((dcm.get_matrix()*dcm.get_matrix().transpose()).trace(), Catch::Matchers::WithinULP(3.0, 4));
}

TEST_CASE("from_eulerangles_space312", "[domain]")
{
    /*!
    * Tests the eanglespace231 function.
    */

    using Lielab::domain::SO;
    using Lielab::domain::CompositeGroup;
    using Lielab::testing::check_almost_equal_nulp;
    
    SO dcm(3);

    // Identity
    dcm = Lielab::domain::SO::from_eulerangles_space312(0.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMId}, CompositeGroup{dcm}, 0));

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_space312(0.0, PI/2.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotx}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_space312(0.0, 0.0, PI/2.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMroty}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_space312(PI/2.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotz}, CompositeGroup{dcm}, 4, true));

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_space312(some_angle, some_angle, some_angle);
    CHECK_THAT(dcm.get_matrix().determinant(), Catch::Matchers::WithinULP(1.0, 4));
    CHECK_THAT((dcm.get_matrix()*dcm.get_matrix().transpose()).trace(), Catch::Matchers::WithinULP(3.0, 4));
}

TEST_CASE("from_eulerangles_space132", "[domain]")
{
    /*!
    * Tests the eanglespace132 function.
    */

    using Lielab::domain::SO;
    using Lielab::domain::CompositeGroup;
    using Lielab::testing::check_almost_equal_nulp;
    
    SO dcm(3);

    // Identity
    dcm = Lielab::domain::SO::from_eulerangles_space132(0.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMId}, CompositeGroup{dcm}, 0));

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_space132(PI/2.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotx}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_space132(0.0, 0.0, PI/2.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMroty}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_space132(0.0, PI/2.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotz}, CompositeGroup{dcm}, 4, true));

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_space132(some_angle, some_angle, some_angle);
    CHECK_THAT(dcm.get_matrix().determinant(), Catch::Matchers::WithinULP(1.0, 4));
    CHECK_THAT((dcm.get_matrix()*dcm.get_matrix().transpose()).trace(), Catch::Matchers::WithinULP(3.0, 4));
}

TEST_CASE("from_eulerangles_space213", "[domain]")
{
    /*!
    * Tests the eanglespace213 function.
    */

    using Lielab::domain::SO;
    using Lielab::domain::CompositeGroup;
    using Lielab::testing::check_almost_equal_nulp;
    
    SO dcm(3);

    // Identity
    dcm = Lielab::domain::SO::from_eulerangles_space213(0.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMId}, CompositeGroup{dcm}, 0));

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_space213(0.0, PI/2.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotx}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_space213(PI/2.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMroty}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_space213(0.0, 0.0, PI/2.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotz}, CompositeGroup{dcm}, 4, true));

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_space213(some_angle, some_angle, some_angle);
    CHECK_THAT(dcm.get_matrix().determinant(), Catch::Matchers::WithinULP(1.0, 4));
    CHECK_THAT((dcm.get_matrix()*dcm.get_matrix().transpose()).trace(), Catch::Matchers::WithinULP(3.0, 4));
}

TEST_CASE("from_eulerangles_space321", "[domain]")
{
    /*!
    * Tests the eanglespace321 function.
    */

    using Lielab::domain::SO;
    using Lielab::domain::CompositeGroup;
    using Lielab::testing::check_almost_equal_nulp;
    
    SO dcm(3);

    // Identity
    dcm = Lielab::domain::SO::from_eulerangles_space321(0.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMId}, CompositeGroup{dcm}, 0));

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_space321(0.0, 0.0, PI/2.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotx}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_space321(0.0, PI/2.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMroty}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_space321(PI/2.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotz}, CompositeGroup{dcm}, 4, true));

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_space321(some_angle, some_angle, some_angle);
    CHECK_THAT(dcm.get_matrix().determinant(), Catch::Matchers::WithinULP(1.0, 4));
    CHECK_THAT((dcm.get_matrix()*dcm.get_matrix().transpose()).trace(), Catch::Matchers::WithinULP(3.0, 4));
}

TEST_CASE("from_eulerangles_space121", "[domain]")
{
    /*!
    * Tests the eanglespace121 function.
    */

    using Lielab::domain::SO;
    using Lielab::domain::CompositeGroup;
    using Lielab::testing::check_almost_equal_nulp;
    
    SO dcm(3);

    // Identity
    dcm = Lielab::domain::SO::from_eulerangles_space121(0.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMId}, CompositeGroup{dcm}, 0));

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_space121(PI/2.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotx}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_space121(0.0, 0.0, PI/2.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotx}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_space121(0.0, PI/2.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMroty}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_space121(PI/2.0, -PI/2.0, -PI/2.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotz}, CompositeGroup{dcm}, 4, true));

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_space121(some_angle, some_angle, some_angle);
    CHECK_THAT(dcm.get_matrix().determinant(), Catch::Matchers::WithinULP(1.0, 4));
    CHECK_THAT((dcm.get_matrix()*dcm.get_matrix().transpose()).trace(), Catch::Matchers::WithinULP(3.0, 4));
}

TEST_CASE("from_eulerangles_space131", "[domain]")
{
    /*!
    * Tests the eanglespace131 function.
    */

    using Lielab::domain::SO;
    using Lielab::domain::CompositeGroup;
    using Lielab::testing::check_almost_equal_nulp;
    
    SO dcm(3);

    // Identity
    dcm = Lielab::domain::SO::from_eulerangles_space131(0.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMId}, CompositeGroup{dcm}, 0));

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_space131(PI/2.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotx}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_space131(0.0, 0.0, PI/2.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotx}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_space131(PI/2.0, PI/2.0, -PI/2.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMroty}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_space131(0.0, PI/2.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotz}, CompositeGroup{dcm}, 4, true));

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_space131(some_angle, some_angle, some_angle);
    CHECK_THAT(dcm.get_matrix().determinant(), Catch::Matchers::WithinULP(1.0, 4));
    CHECK_THAT((dcm.get_matrix()*dcm.get_matrix().transpose()).trace(), Catch::Matchers::WithinULP(3.0, 4));
}

TEST_CASE("from_eulerangles_space212", "[domain]")
{
    /*!
    * Tests the eanglespace212 function.
    */

    using Lielab::domain::SO;
    using Lielab::domain::CompositeGroup;
    using Lielab::testing::check_almost_equal_nulp;
    
    SO dcm(3);

    // Identity
    dcm = Lielab::domain::SO::from_eulerangles_space212(0.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMId}, CompositeGroup{dcm}, 0));

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_space212(0.0, PI/2.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotx}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_space212(PI/2.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMroty}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_space212(0.0, 0.0, PI/2.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMroty}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_space212(PI/2.0, PI/2.0, -PI/2.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotz}, CompositeGroup{dcm}, 4, true));

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_space212(some_angle, some_angle, some_angle);
    CHECK_THAT(dcm.get_matrix().determinant(), Catch::Matchers::WithinULP(1.0, 4));
    CHECK_THAT((dcm.get_matrix()*dcm.get_matrix().transpose()).trace(), Catch::Matchers::WithinULP(3.0, 4));
}

TEST_CASE("from_eulerangles_space232", "[domain]")
{
    /*!
    * Tests the eanglespace232 function.
    */

    using Lielab::domain::SO;
    using Lielab::domain::CompositeGroup;
    using Lielab::testing::check_almost_equal_nulp;
    
    SO dcm(3);

    // Identity
    dcm = Lielab::domain::SO::from_eulerangles_space232(0.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMId}, CompositeGroup{dcm}, 0));

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_space232(PI/2.0, -PI/2.0, -PI/2.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotx}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_space232(PI/2.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMroty}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_space232(0.0, 0.0, PI/2.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMroty}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_space232(0.0, PI/2.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotz}, CompositeGroup{dcm}, 4, true));

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_space232(some_angle, some_angle, some_angle);
    CHECK_THAT(dcm.get_matrix().determinant(), Catch::Matchers::WithinULP(1.0, 4));
    CHECK_THAT((dcm.get_matrix()*dcm.get_matrix().transpose()).trace(), Catch::Matchers::WithinULP(3.0, 4));
}

TEST_CASE("from_eulerangles_space313", "[domain]")
{
    /*!
    * Tests the eanglespace313 function.
    */

    using Lielab::domain::SO;
    using Lielab::domain::CompositeGroup;
    using Lielab::testing::check_almost_equal_nulp;
    
    SO dcm(3);

    // Identity
    dcm = Lielab::domain::SO::from_eulerangles_space313(0.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMId}, CompositeGroup{dcm}, 0));

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_space313(0.0, PI/2.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotx}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_space313(-PI/2.0, PI/2.0, PI/2.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMroty}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_space313(PI/2.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotz}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_space313(0.0, 0.0, PI/2.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotz}, CompositeGroup{dcm}, 4, true));

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_space313(some_angle, some_angle, some_angle);
    CHECK_THAT(dcm.get_matrix().determinant(), Catch::Matchers::WithinULP(1.0, 4));
    CHECK_THAT((dcm.get_matrix()*dcm.get_matrix().transpose()).trace(), Catch::Matchers::WithinULP(3.0, 4));
}

TEST_CASE("from_eulerangles_space323", "[domain]")
{
    /*!
    * Tests the eanglespace323 function.
    */

    using Lielab::domain::SO;
    using Lielab::domain::CompositeGroup;
    using Lielab::testing::check_almost_equal_nulp;
    
    SO dcm(3);

    // Identity
    dcm = Lielab::domain::SO::from_eulerangles_space323(0.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMId}, CompositeGroup{dcm}, 0));

    // Rotate 90 degrees by x-axis
    dcm = Lielab::domain::SO::from_eulerangles_space323(PI/2.0, PI/2.0, -PI/2.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotx}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by y-axis
    dcm = Lielab::domain::SO::from_eulerangles_space323(0.0, PI/2.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMroty}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_space323(PI/2.0, 0.0, 0.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotz}, CompositeGroup{dcm}, 4, true));

    // Rotate 90 degrees by z-axis
    dcm = Lielab::domain::SO::from_eulerangles_space323(0.0, 0.0, PI/2.0);
    CHECK(check_almost_equal_nulp(CompositeGroup{DCMrotz}, CompositeGroup{dcm}, 4, true));

    // Rotate by a random angle
    dcm = Lielab::domain::SO::from_eulerangles_space323(some_angle, some_angle, some_angle);
    CHECK_THAT(dcm.get_matrix().determinant(), Catch::Matchers::WithinULP(1.0, 4));
    CHECK_THAT((dcm.get_matrix()*dcm.get_matrix().transpose()).trace(), Catch::Matchers::WithinULP(3.0, 4));
}

TEST_CASE("to_eulerangles_body123", "[domain]")
{
    /*!
    * Tests the to_eulerangles_body123 function.
    */

    using Lielab::testing::check_almost_equal_nulp;

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_body123());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_body123());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_body123();
    CHECK(check_almost_equal_nulp(ex1theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta3, 0.0, 1, true));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_body123();
    CHECK(check_almost_equal_nulp(ex2theta1, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta3, 0.0, 1, true));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_body123();
    CHECK(check_almost_equal_nulp(ex3theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta3, 0.0, 1, true));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_body123();
    CHECK(check_almost_equal_nulp(ex4theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta3, PI/2.0, 1, true));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_body123(0.5, PI/2.0, 0.5).to_eulerangles_body123();
    CHECK(check_almost_equal_nulp(ex5theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta3, 0.0, 1, true));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_body123(1.0, PI/2.0, 0.0).to_eulerangles_body123();
    CHECK(check_almost_equal_nulp(ex5g1theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta3, 0.0, 1, true));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_body123(1.0, -PI/2.0, 0.0).to_eulerangles_body123();
    CHECK(check_almost_equal_nulp(ex5g2theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta2, -PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta3, 0.0, 1, true));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_body123(some_angle, some_angle, some_angle).to_eulerangles_body123();
    CHECK(check_almost_equal_nulp(ex6theta1, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta2, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta3, some_angle, 1, true));
}

TEST_CASE("to_eulerangles_body231", "[domain]")
{
    /*!
    * Tests the to_eulerangles_body231 function.
    */

    using Lielab::testing::check_almost_equal_nulp;

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_body231());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_body231());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_body231();
    CHECK(check_almost_equal_nulp(ex1theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta3, 0.0, 1, true));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_body231();
    CHECK(check_almost_equal_nulp(ex2theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta3, PI/2.0, 1, true));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_body231();
    CHECK(check_almost_equal_nulp(ex3theta1, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta3, 0.0, 1, true));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_body231();
    CHECK(check_almost_equal_nulp(ex4theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta3, 0.0, 1, true));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_body231(0.5, PI/2.0, 0.5).to_eulerangles_body231();
    CHECK(check_almost_equal_nulp(ex5theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta3, 0.0, 1, true));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_body231(1.0, PI/2.0, 0.0).to_eulerangles_body231();
    CHECK(check_almost_equal_nulp(ex5g1theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta3, 0.0, 1, true));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_body231(1.0, -PI/2.0, 0.0).to_eulerangles_body231();
    CHECK(check_almost_equal_nulp(ex5g2theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta2, -PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta3, 0.0, 1, true));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_body231(some_angle, some_angle, some_angle).to_eulerangles_body231();
    CHECK(check_almost_equal_nulp(ex6theta1, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta2, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta3, some_angle, 1, true));
}

TEST_CASE("to_eulerangles_body312", "[domain]")
{
    /*!
    * Tests the to_eulerangles_body312 function.
    */

    using Lielab::testing::check_almost_equal_nulp;

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_body312());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_body312());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_body312();
    CHECK(check_almost_equal_nulp(ex1theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta3, 0.0, 1, true));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_body312();
    CHECK(check_almost_equal_nulp(ex2theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta3, 0.0, 1, true));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_body312();
    CHECK(check_almost_equal_nulp(ex3theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta3, PI/2.0, 1, true));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_body312();
    CHECK(check_almost_equal_nulp(ex4theta1, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta3, 0.0, 1, true));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_body312(0.5, PI/2.0, 0.5).to_eulerangles_body312();
    CHECK(check_almost_equal_nulp(ex5theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta3, 0.0, 1, true));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_body312(1.0, PI/2.0, 0.0).to_eulerangles_body312();
    CHECK(check_almost_equal_nulp(ex5g1theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta3, 0.0, 1, true));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_body312(1.0, -PI/2.0, 0.0).to_eulerangles_body312();
    CHECK(check_almost_equal_nulp(ex5g2theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta2, -PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta3, 0.0, 1, true));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_body312(some_angle, some_angle, some_angle).to_eulerangles_body312();
    CHECK(check_almost_equal_nulp(ex6theta1, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta2, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta3, some_angle, 1, true));
}

TEST_CASE("to_eulerangles_body132", "[domain]")
{
    /*!
    * Tests the to_eulerangles_body132 function.
    */

    using Lielab::testing::check_almost_equal_nulp;

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_body132());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_body132());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_body132();
    CHECK(check_almost_equal_nulp(ex1theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta3, 0.0, 1, true));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_body132();
    CHECK(check_almost_equal_nulp(ex2theta1, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta3, 0.0, 1, true));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_body132();
    CHECK(check_almost_equal_nulp(ex3theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta3, PI/2.0, 1, true));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_body132();
    CHECK(check_almost_equal_nulp(ex4theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta3, 0.0, 1, true));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_body132(0.5, PI/2.0, 0.5).to_eulerangles_body132();
    CHECK(check_almost_equal_nulp(ex5theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta3, 0.0, 1, true));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_body132(1.0, PI/2.0, 0.0).to_eulerangles_body132();
    CHECK(check_almost_equal_nulp(ex5g1theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta3, 0.0, 1, true));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_body132(1.0, -PI/2.0, 0.0).to_eulerangles_body132();
    CHECK(check_almost_equal_nulp(ex5g2theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta2, -PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta3, 0.0, 1, true));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_body132(some_angle, some_angle, some_angle).to_eulerangles_body132();
    CHECK(check_almost_equal_nulp(ex6theta1, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta2, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta3, some_angle, 1, true));
}

TEST_CASE("to_eulerangles_body213", "[domain]")
{
    /*!
    * Tests the to_eulerangles_body213 function.
    */

    using Lielab::testing::check_almost_equal_nulp;

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_body213());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_body213());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_body213();
    CHECK(check_almost_equal_nulp(ex1theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta3, 0.0, 1, true));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_body213();
    CHECK(check_almost_equal_nulp(ex2theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta3, 0.0, 1, true));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_body213();
    CHECK(check_almost_equal_nulp(ex3theta1, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta3, 0.0, 1, true));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_body213();
    CHECK(check_almost_equal_nulp(ex4theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta3, PI/2.0, 1, true));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_body213(0.5, PI/2.0, 0.5).to_eulerangles_body213();
    CHECK(check_almost_equal_nulp(ex5theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta3, 0.0, 1, true));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_body213(1.0, PI/2.0, 0.0).to_eulerangles_body213();
    CHECK(check_almost_equal_nulp(ex5g1theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta3, 0.0, 1, true));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_body213(1.0, -PI/2.0, 0.0).to_eulerangles_body213();
    CHECK(check_almost_equal_nulp(ex5g2theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta2, -PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta3, 0.0, 1, true));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_body213(some_angle, some_angle, some_angle).to_eulerangles_body213();
    CHECK(check_almost_equal_nulp(ex6theta1, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta2, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta3, some_angle, 1, true));
}

TEST_CASE("to_eulerangles_body321", "[domain]")
{
    /*!
    * Tests the to_eulerangles_body321 function.
    */

    using Lielab::testing::check_almost_equal_nulp;

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_body321());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_body321());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_body321();
    CHECK(check_almost_equal_nulp(ex1theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta3, 0.0, 1, true));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_body321();
    CHECK(check_almost_equal_nulp(ex2theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta3, PI/2.0, 1, true));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_body321();
    CHECK(check_almost_equal_nulp(ex3theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta3, 0.0, 1, true));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_body321();
    CHECK(check_almost_equal_nulp(ex4theta1, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta3, 0.0, 1, true));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_body321(0.5, PI/2.0, 0.5).to_eulerangles_body321();
    CHECK(check_almost_equal_nulp(ex5theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta3, 0.0, 1, true));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_body321(1.0, PI/2.0, 0.0).to_eulerangles_body321();
    CHECK(check_almost_equal_nulp(ex5g1theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta3, 0.0, 1, true));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_body321(1.0, -PI/2.0, 0.0).to_eulerangles_body321();
    CHECK(check_almost_equal_nulp(ex5g2theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta2, -PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta3, 0.0, 1, true));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_body321(some_angle, some_angle, some_angle).to_eulerangles_body321();
    CHECK(check_almost_equal_nulp(ex6theta1, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta2, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta3, some_angle, 1, true));
}

TEST_CASE("to_eulerangles_body121", "[domain]")
{
    /*!
    * Tests the to_eulerangles_body121 function.
    */

    using Lielab::testing::check_almost_equal_nulp;

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_body121());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_body121());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_body121();
    CHECK(check_almost_equal_nulp(ex1theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta3, 0.0, 1, true));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_body121();
    CHECK(check_almost_equal_nulp(ex2theta1, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta3, 0.0, 1, true));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_body121();
    CHECK(check_almost_equal_nulp(ex3theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta3, 0.0, 1, true));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_body121();
    CHECK(check_almost_equal_nulp(ex4theta1, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta3, -PI/2.0, 1, true));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_body121(0.5, 0.0, 0.5).to_eulerangles_body121();
    CHECK(check_almost_equal_nulp(ex5theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta3, 0.0, 1, true));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_body121(1.0, 0.0, 0.0).to_eulerangles_body121();
    CHECK(check_almost_equal_nulp(ex5g1theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta3, 0.0, 1, true));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_body121(1.0, PI, 0.0).to_eulerangles_body121();
    CHECK(check_almost_equal_nulp(ex5g2theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta2, PI, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta3, 0.0, 1, true));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_body121(some_angle, some_angle, some_angle).to_eulerangles_body121();
    CHECK(check_almost_equal_nulp(ex6theta1, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta2, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta3, some_angle, 1, true));
}

TEST_CASE("to_eulerangles_body131", "[domain]")
{
    /*!
    * Tests the to_eulerangles_body131 function.
    */

    using Lielab::testing::check_almost_equal_nulp;

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_body131());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_body131());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_body131();
    CHECK(check_almost_equal_nulp(ex1theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta3, 0.0, 1, true));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_body131();
    CHECK(check_almost_equal_nulp(ex2theta1, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta3, 0.0, 1, true));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_body131();
    CHECK(check_almost_equal_nulp(ex3theta1, -PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta3, PI/2.0, 1, true));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_body131();
    CHECK(check_almost_equal_nulp(ex4theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta3, 0.0, 1, true));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_body131(0.5, 0.0, 0.5).to_eulerangles_body131();
    CHECK(check_almost_equal_nulp(ex5theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta3, 0.0, 1, true));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_body131(1.0, 0.0, 0.0).to_eulerangles_body131();
    CHECK(check_almost_equal_nulp(ex5g1theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta3, 0.0, 1, true));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_body131(1.0, PI, 0.0).to_eulerangles_body131();
    CHECK(check_almost_equal_nulp(ex5g2theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta2, PI, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta3, 0.0, 1, true));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_body131(some_angle, some_angle, some_angle).to_eulerangles_body131();
    CHECK(check_almost_equal_nulp(ex6theta1, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta2, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta3, some_angle, 1, true));
}

TEST_CASE("to_eulerangles_body212", "[domain]")
{
    /*!
    * Tests the to_eulerangles_body212 function.
    */

    using Lielab::testing::check_almost_equal_nulp;

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_body212());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_body212());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_body212();
    CHECK(check_almost_equal_nulp(ex1theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta3, 0.0, 1, true));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_body212();
    CHECK(check_almost_equal_nulp(ex2theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta3, 0.0, 1, true));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_body212();
    CHECK(check_almost_equal_nulp(ex3theta1, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta3, 0.0, 1, true));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_body212();
    CHECK(check_almost_equal_nulp(ex4theta1, -PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta3, PI/2.0, 1, true));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_body212(0.5, 0.0, 0.5).to_eulerangles_body212();
    CHECK(check_almost_equal_nulp(ex5theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta3, 0.0, 1, true));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_body212(1.0, 0.0, 0.0).to_eulerangles_body212();
    CHECK(check_almost_equal_nulp(ex5g1theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta3, 0.0, 1, true));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_body212(1.0, PI, 0.0).to_eulerangles_body212();
    CHECK(check_almost_equal_nulp(ex5g2theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta2, PI, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta3, 0.0, 1, true));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_body212(some_angle, some_angle, some_angle).to_eulerangles_body212();
    CHECK(check_almost_equal_nulp(ex6theta1, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta2, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta3, some_angle, 1, true));
}

TEST_CASE("to_eulerangles_body232", "[domain]")
{
    /*!
    * Tests the to_eulerangles_body232 function.
    */

    using Lielab::testing::check_almost_equal_nulp;

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_body232());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_body232());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_body232();
    CHECK(check_almost_equal_nulp(ex1theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta3, 0.0, 1, true));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_body232();
    CHECK(check_almost_equal_nulp(ex2theta1, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta3, -PI/2.0, 1, true));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_body232();
    CHECK(check_almost_equal_nulp(ex3theta1, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta3, 0.0, 1, true));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_body232();
    CHECK(check_almost_equal_nulp(ex4theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta3, 0.0, 1, true));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_body232(0.5, 0.0, 0.5).to_eulerangles_body232();
    CHECK(check_almost_equal_nulp(ex5theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta3, 0.0, 1, true));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_body232(1.0, 0.0, 0.0).to_eulerangles_body232();
    CHECK(check_almost_equal_nulp(ex5g1theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta3, 0.0, 1, true));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_body232(1.0, PI, 0.0).to_eulerangles_body232();
    CHECK(check_almost_equal_nulp(ex5g2theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta2, PI, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta3, 0.0, 1, true));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_body232(some_angle, some_angle, some_angle).to_eulerangles_body232();
    CHECK(check_almost_equal_nulp(ex6theta1, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta2, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta3, some_angle, 1, true));
}

TEST_CASE("to_eulerangles_body313", "[domain]")
{
    /*!
    * Tests the to_eulerangles_body313 function.
    */

    using Lielab::testing::check_almost_equal_nulp;

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_body313());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_body313());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_body313();
    CHECK(check_almost_equal_nulp(ex1theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta3, 0.0, 1, true));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_body313();
    CHECK(check_almost_equal_nulp(ex2theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta3, 0.0, 1, true));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_body313();
    CHECK(check_almost_equal_nulp(ex3theta1, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta3, -PI/2.0, 1, true));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_body313();
    CHECK(check_almost_equal_nulp(ex4theta1, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta3, 0.0, 1, true));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_body313(0.5, 0.0, 0.5).to_eulerangles_body313();
    CHECK(check_almost_equal_nulp(ex5theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta3, 0.0, 1, true));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_body313(1.0, 0.0, 0.0).to_eulerangles_body313();
    CHECK(check_almost_equal_nulp(ex5g1theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta3, 0.0, 1, true));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_body313(1.0, PI, 0.0).to_eulerangles_body313();
    CHECK(check_almost_equal_nulp(ex5g2theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta2, PI, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta3, 0.0, 1, true));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_body313(some_angle, some_angle, some_angle).to_eulerangles_body313();
    CHECK(check_almost_equal_nulp(ex6theta1, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta2, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta3, some_angle, 1, true));
}

TEST_CASE("to_eulerangles_body323", "[domain]")
{
    /*!
    * Tests the to_eulerangles_body323 function.
    */

    using Lielab::testing::check_almost_equal_nulp;

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_body323());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_body323());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_body323();
    CHECK(check_almost_equal_nulp(ex1theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta3, 0.0, 1, true));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_body323();
    CHECK(check_almost_equal_nulp(ex2theta1, -PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta3, PI/2.0, 1, true));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_body323();
    CHECK(check_almost_equal_nulp(ex3theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta3, 0.0, 1, true));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_body323();
    CHECK(check_almost_equal_nulp(ex4theta1, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta3, 0.0, 1, true));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_body323(0.5, 0.0, 0.5).to_eulerangles_body323();
    CHECK(check_almost_equal_nulp(ex5theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta3, 0.0, 1, true));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_body323(1.0, 0.0, 0.0).to_eulerangles_body323();
    CHECK(check_almost_equal_nulp(ex5g1theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta3, 0.0, 1, true));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_body323(1.0, PI, 0.0).to_eulerangles_body323();
    CHECK(check_almost_equal_nulp(ex5g2theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta2, PI, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta3, 0.0, 1, true));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_body323(some_angle, some_angle, some_angle).to_eulerangles_body323();
    CHECK(check_almost_equal_nulp(ex6theta1, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta2, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta3, some_angle, 1, true));
}

TEST_CASE("to_eulerangles_space123", "[domain]")
{
    /*!
    * Tests the to_eulerangles_space123 function.
    */

    using Lielab::testing::check_almost_equal_nulp;

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_space123());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_space123());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_space123();
    CHECK(check_almost_equal_nulp(ex1theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta3, 0.0, 1, true));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_space123();
    CHECK(check_almost_equal_nulp(ex2theta1, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta3, 0.0, 1, true));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_space123();
    CHECK(check_almost_equal_nulp(ex3theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta3, 0.0, 1, true));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_space123();
    CHECK(check_almost_equal_nulp(ex4theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta3, PI/2.0, 1, true));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_space123(0.5, PI/2.0, 0.5).to_eulerangles_space123();
    CHECK(check_almost_equal_nulp(ex5theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta3, 0.0, 1, true));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_space123(1.0, PI/2.0, 0.0).to_eulerangles_space123();
    CHECK(check_almost_equal_nulp(ex5g1theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta3, 0.0, 1, true));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_space123(1.0, -PI/2.0, 0.0).to_eulerangles_space123();
    CHECK(check_almost_equal_nulp(ex5g2theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta2, -PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta3, 0.0, 1, true));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_space123(some_angle, some_angle, some_angle).to_eulerangles_space123();
    CHECK(check_almost_equal_nulp(ex6theta1, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta2, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta3, some_angle, 1, true));
}


TEST_CASE("to_eulerangles_space231", "[domain]")
{
    /*!
    * Tests the to_eulerangles_space231 function.
    */

    using Lielab::testing::check_almost_equal_nulp;

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_space231());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_space231());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_space231();
    CHECK(check_almost_equal_nulp(ex1theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta3, 0.0, 1, true));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_space231();
    CHECK(check_almost_equal_nulp(ex2theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta3, PI/2.0, 1, true));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_space231();
    CHECK(check_almost_equal_nulp(ex3theta1, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta3, 0.0, 1, true));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_space231();
    CHECK(check_almost_equal_nulp(ex4theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta3, 0.0, 1, true));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_space231(0.5, PI/2.0, 0.5).to_eulerangles_space231();
    CHECK(check_almost_equal_nulp(ex5theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta3, 0.0, 1, true));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_space231(1.0, PI/2.0, 0.0).to_eulerangles_space231();
    CHECK(check_almost_equal_nulp(ex5g1theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta3, 0.0, 1, true));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_space231(1.0, -PI/2.0, 0.0).to_eulerangles_space231();
    CHECK(check_almost_equal_nulp(ex5g2theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta2, -PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta3, 0.0, 1, true));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_space231(some_angle, some_angle, some_angle).to_eulerangles_space231();
    CHECK(check_almost_equal_nulp(ex6theta1, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta2, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta3, some_angle, 1, true));
}


TEST_CASE("to_eulerangles_space312", "[domain]")
{
    /*!
    * Tests the to_eulerangles_space312 function.
    */

    using Lielab::testing::check_almost_equal_nulp;

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_space312());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_space312());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_space312();
    CHECK(check_almost_equal_nulp(ex1theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta3, 0.0, 1, true));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_space312();
    CHECK(check_almost_equal_nulp(ex2theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta3, 0.0, 1, true));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_space312();
    CHECK(check_almost_equal_nulp(ex3theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta3, PI/2.0, 1, true));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_space312();
    CHECK(check_almost_equal_nulp(ex4theta1, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta3, 0.0, 1, true));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_space312(0.5, PI/2.0, 0.5).to_eulerangles_space312();
    CHECK(check_almost_equal_nulp(ex5theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta3, 0.0, 1, true));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_space312(1.0, PI/2.0, 0.0).to_eulerangles_space312();
    CHECK(check_almost_equal_nulp(ex5g1theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta3, 0.0, 1, true));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_space312(1.0, -PI/2.0, 0.0).to_eulerangles_space312();
    CHECK(check_almost_equal_nulp(ex5g2theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta2, -PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta3, 0.0, 1, true));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_space312(some_angle, some_angle, some_angle).to_eulerangles_space312();
    CHECK(check_almost_equal_nulp(ex6theta1, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta2, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta3, some_angle, 1, true));
}

TEST_CASE("to_eulerangles_space132", "[domain]")
{
    /*!
    * Tests the to_eulerangles_space132 function.
    */

    using Lielab::testing::check_almost_equal_nulp;

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_space132());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_space132());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_space132();
    CHECK(check_almost_equal_nulp(ex1theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta3, 0.0, 1, true));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_space132();
    CHECK(check_almost_equal_nulp(ex2theta1, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta3, 0.0, 1, true));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_space132();
    CHECK(check_almost_equal_nulp(ex3theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta3, PI/2.0, 1, true));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_space132();
    CHECK(check_almost_equal_nulp(ex4theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta3, 0.0, 1, true));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_space132(0.5, PI/2.0, 0.5).to_eulerangles_space132();
    CHECK(check_almost_equal_nulp(ex5theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta3, 0.0, 1, true));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_space132(1.0, PI/2.0, 0.0).to_eulerangles_space132();
    CHECK(check_almost_equal_nulp(ex5g1theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta3, 0.0, 1, true));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_space132(1.0, -PI/2.0, 0.0).to_eulerangles_space132();
    CHECK(check_almost_equal_nulp(ex5g2theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta2, -PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta3, 0.0, 1, true));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_space132(some_angle, some_angle, some_angle).to_eulerangles_space132();
    CHECK(check_almost_equal_nulp(ex6theta1, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta2, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta3, some_angle, 1, true));
}

TEST_CASE("to_eulerangles_space213", "[domain]")
{
    /*!
    * Tests the to_eulerangles_space213 function.
    */

    using Lielab::testing::check_almost_equal_nulp;

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_space213());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_space213());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_space213();
    CHECK(check_almost_equal_nulp(ex1theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta3, 0.0, 1, true));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_space213();
    CHECK(check_almost_equal_nulp(ex2theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta3, 0.0, 1, true));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_space213();
    CHECK(check_almost_equal_nulp(ex3theta1, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta3, 0.0, 1, true));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_space213();
    CHECK(check_almost_equal_nulp(ex4theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta3, PI/2.0, 1, true));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_space213(0.5, PI/2.0, 0.5).to_eulerangles_space213();
    CHECK(check_almost_equal_nulp(ex5theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta3, 0.0, 1, true));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_space213(1.0, PI/2.0, 0.0).to_eulerangles_space213();
    CHECK(check_almost_equal_nulp(ex5g1theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta3, 0.0, 1, true));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_space213(1.0, -PI/2.0, 0.0).to_eulerangles_space213();
    CHECK(check_almost_equal_nulp(ex5g2theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta2, -PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta3, 0.0, 1, true));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_space213(some_angle, some_angle, some_angle).to_eulerangles_space213();
    CHECK(check_almost_equal_nulp(ex6theta1, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta2, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta3, some_angle, 1, true));
}

TEST_CASE("to_eulerangles_space321", "[domain]")
{
    /*!
    * Tests the to_eulerangles_space321 function.
    */

    using Lielab::testing::check_almost_equal_nulp;

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_space321());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_space321());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_space321();
    CHECK(check_almost_equal_nulp(ex1theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta3, 0.0, 1, true));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_space321();
    CHECK(check_almost_equal_nulp(ex2theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta3, PI/2.0, 1, true));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_space321();
    CHECK(check_almost_equal_nulp(ex3theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta3, 0.0, 1, true));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_space321();
    CHECK(check_almost_equal_nulp(ex4theta1, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta3, 0.0, 1, true));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_space321(0.5, PI/2.0, 0.5).to_eulerangles_space321();
    CHECK(check_almost_equal_nulp(ex5theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta3, 0.0, 1, true));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_space321(1.0, PI/2.0, 0.0).to_eulerangles_space321();
    CHECK(check_almost_equal_nulp(ex5g1theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta3, 0.0, 1, true));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_space321(1.0, -PI/2.0, 0.0).to_eulerangles_space321();
    CHECK(check_almost_equal_nulp(ex5g2theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta2, -PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta3, 0.0, 1, true));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_space321(some_angle, some_angle, some_angle).to_eulerangles_space321();
    CHECK(check_almost_equal_nulp(ex6theta1, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta2, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta3, some_angle, 1, true));
}

TEST_CASE("to_eulerangles_space121", "[domain]")
{
    /*!
    * Tests the to_eulerangles_space121 function.
    */

    using Lielab::testing::check_almost_equal_nulp;

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_space121());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_space121());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_space121();
    CHECK(check_almost_equal_nulp(ex1theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta3, 0.0, 1, true));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_space121();
    CHECK(check_almost_equal_nulp(ex2theta1, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta3, 0.0, 1, true));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_space121();
    CHECK(check_almost_equal_nulp(ex3theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta3, 0.0, 1, true));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_space121();
    CHECK(check_almost_equal_nulp(ex4theta1, -PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta3, PI/2.0, 1, true));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_space121(0.5, 0.0, 0.5).to_eulerangles_space121();
    CHECK(check_almost_equal_nulp(ex5theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta3, 0.0, 1, true));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_space121(1.0, 0.0, 0.0).to_eulerangles_space121();
    CHECK(check_almost_equal_nulp(ex5g1theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta3, 0.0, 1, true));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_space121(1.0, PI, 0.0).to_eulerangles_space121();
    CHECK(check_almost_equal_nulp(ex5g2theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta2, PI, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta3, 0.0, 1, true));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_space121(some_angle, some_angle, some_angle).to_eulerangles_space121();
    CHECK(check_almost_equal_nulp(ex6theta1, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta2, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta3, some_angle, 1, true));
}

TEST_CASE("to_eulerangles_space131", "[domain]")
{
    /*!
    * Tests the to_eulerangles_space131 function.
    */

    using Lielab::testing::check_almost_equal_nulp;

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_space131());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_space131());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_space131();
    CHECK(check_almost_equal_nulp(ex1theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta3, 0.0, 1, true));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_space131();
    CHECK(check_almost_equal_nulp(ex2theta1, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta3, 0.0, 1, true));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_space131();
    CHECK(check_almost_equal_nulp(ex3theta1, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta3, -PI/2.0, 1, true));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_space131();
    CHECK(check_almost_equal_nulp(ex4theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta3, 0.0, 1, true));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_space131(0.5, 0.0, 0.5).to_eulerangles_space131();
    CHECK(check_almost_equal_nulp(ex5theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta3, 0.0, 1, true));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_space131(1.0, 0.0, 0.0).to_eulerangles_space131();
    CHECK(check_almost_equal_nulp(ex5g1theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta3, 0.0, 1, true));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_space131(1.0, PI, 0.0).to_eulerangles_space131();
    CHECK(check_almost_equal_nulp(ex5g2theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta2, PI, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta3, 0.0, 1, true));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_space131(some_angle, some_angle, some_angle).to_eulerangles_space131();
    CHECK(check_almost_equal_nulp(ex6theta1, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta2, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta3, some_angle, 1, true));
}

TEST_CASE("to_eulerangles_space212", "[domain]")
{
    /*!
    * Tests the to_eulerangles_space212 function.
    */

    using Lielab::testing::check_almost_equal_nulp;

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_space212());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_space212());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_space212();
    CHECK(check_almost_equal_nulp(ex1theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta3, 0.0, 1, true));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_space212();
    CHECK(check_almost_equal_nulp(ex2theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta3, 0.0, 1, true));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_space212();
    CHECK(check_almost_equal_nulp(ex3theta1, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta3, 0.0, 1, true));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_space212();
    CHECK(check_almost_equal_nulp(ex4theta1, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta3, -PI/2.0, 1, true));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_space212(0.5, 0.0, 0.5).to_eulerangles_space212();
    CHECK(check_almost_equal_nulp(ex5theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta3, 0.0, 1, true));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_space212(1.0, 0.0, 0.0).to_eulerangles_space212();
    CHECK(check_almost_equal_nulp(ex5g1theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta3, 0.0, 1, true));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_space212(1.0, PI, 0.0).to_eulerangles_space212();
    CHECK(check_almost_equal_nulp(ex5g2theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta2, PI, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta3, 0.0, 1, true));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_space212(some_angle, some_angle, some_angle).to_eulerangles_space212();
    CHECK(check_almost_equal_nulp(ex6theta1, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta2, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta3, some_angle, 1, true));
}

TEST_CASE("to_eulerangles_space232", "[domain]")
{
    /*!
    * Tests the to_eulerangles_space232 function.
    */

    using Lielab::testing::check_almost_equal_nulp;

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_space232());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_space232());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_space232();
    CHECK(check_almost_equal_nulp(ex1theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta3, 0.0, 1, true));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_space232();
    CHECK(check_almost_equal_nulp(ex2theta1, -PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta3, PI/2.0, 1, true));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_space232();
    CHECK(check_almost_equal_nulp(ex3theta1, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta3, 0.0, 1, true));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_space232();
    CHECK(check_almost_equal_nulp(ex4theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta3, 0.0, 1, true));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_space232(0.5, 0.0, 0.5).to_eulerangles_space232();
    CHECK(check_almost_equal_nulp(ex5theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta3, 0.0, 1, true));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_space232(1.0, 0.0, 0.0).to_eulerangles_space232();
    CHECK(check_almost_equal_nulp(ex5g1theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta3, 0.0, 1, true));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_space232(1.0, PI, 0.0).to_eulerangles_space232();
    CHECK(check_almost_equal_nulp(ex5g2theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta2, PI, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta3, 0.0, 1, true));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_space232(some_angle, some_angle, some_angle).to_eulerangles_space232();
    CHECK(check_almost_equal_nulp(ex6theta1, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta2, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta3, some_angle, 1, true));
}

TEST_CASE("to_eulerangles_space313", "[domain]")
{
    /*!
    * Tests the to_eulerangles_space313 function.
    */

    using Lielab::testing::check_almost_equal_nulp;

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_space313());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_space313());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_space313();
    CHECK(check_almost_equal_nulp(ex1theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta3, 0.0, 1, true));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_space313();
    CHECK(check_almost_equal_nulp(ex2theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta3, 0.0, 1, true));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_space313();
    CHECK(check_almost_equal_nulp(ex3theta1, -PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta3, PI/2.0, 1, true));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_space313();
    CHECK(check_almost_equal_nulp(ex4theta1, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta3, 0.0, 1, true));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_space313(0.5, 0.0, 0.5).to_eulerangles_space313();
    CHECK(check_almost_equal_nulp(ex5theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta3, 0.0, 1, true));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_space313(1.0, 0.0, 0.0).to_eulerangles_space313();
    CHECK(check_almost_equal_nulp(ex5g1theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta3, 0.0, 1, true));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_space313(1.0, PI, 0.0).to_eulerangles_space313();
    CHECK(check_almost_equal_nulp(ex5g2theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta2, PI, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta3, 0.0, 1, true));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_space313(some_angle, some_angle, some_angle).to_eulerangles_space313();
    CHECK(check_almost_equal_nulp(ex6theta1, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta2, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta3, some_angle, 1, true));
}

TEST_CASE("to_eulerangles_space323", "[domain]")
{
    /*!
    * Tests the to_eulerangles_space323 function.
    */

    using Lielab::testing::check_almost_equal_nulp;

    // Improperly sized
    CHECK_THROWS(Lielab::domain::SO(2).to_eulerangles_space323());
    CHECK_THROWS(Lielab::domain::SO(4).to_eulerangles_space323());

    // Identity
    const auto [ex1theta1, ex1theta2, ex1theta3] = DCMId.to_eulerangles_space323();
    CHECK(check_almost_equal_nulp(ex1theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex1theta3, 0.0, 1, true));

    // Rotate 90 degrees by x-axis
    const auto [ex2theta1, ex2theta2, ex2theta3] = DCMrotx.to_eulerangles_space323();
    CHECK(check_almost_equal_nulp(ex2theta1, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex2theta3, -PI/2.0, 1, true));

    // Rotate 90 degrees by y-axis
    const auto [ex3theta1, ex3theta2, ex3theta3] = DCMroty.to_eulerangles_space323();
    CHECK(check_almost_equal_nulp(ex3theta1, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta2, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex3theta3, 0.0, 1, true));

    // Rotate 90 degrees by z-axis
    const auto [ex4theta1, ex4theta2, ex4theta3] = DCMrotz.to_eulerangles_space323();
    CHECK(check_almost_equal_nulp(ex4theta1, PI/2.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex4theta3, 0.0, 1, true));

    // Check gimbal lock condition
    const auto [ex5theta1, ex5theta2, ex5theta3] = Lielab::domain::SO::from_eulerangles_space323(0.5, 0.0, 0.5).to_eulerangles_space323();
    CHECK(check_almost_equal_nulp(ex5theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5theta3, 0.0, 1, true));

    const auto [ex5g1theta1, ex5g1theta2, ex5g1theta3] = Lielab::domain::SO::from_eulerangles_space323(1.0, 0.0, 0.0).to_eulerangles_space323();
    CHECK(check_almost_equal_nulp(ex5g1theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta2, 0.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g1theta3, 0.0, 1, true));

    const auto [ex5g2theta1, ex5g2theta2, ex5g2theta3] = Lielab::domain::SO::from_eulerangles_space323(1.0, PI, 0.0).to_eulerangles_space323();
    CHECK(check_almost_equal_nulp(ex5g2theta1, 1.0, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta2, PI, 1, true));
    CHECK(check_almost_equal_nulp(ex5g2theta3, 0.0, 1, true));

    // Rotate by a random angle (also checks the inverse function)
    const auto [ex6theta1, ex6theta2, ex6theta3] = Lielab::domain::SO::from_eulerangles_space323(some_angle, some_angle, some_angle).to_eulerangles_space323();
    CHECK(check_almost_equal_nulp(ex6theta1, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta2, some_angle, 1, true));
    CHECK(check_almost_equal_nulp(ex6theta3, some_angle, 1, true));
}

