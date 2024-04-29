
// TEST_CASE("vfex1", "[dynamics]")
// {
//     /*!
//     * Tests vfex1.
//     */

//     Lielab::dynamics::vectorfield vf = Lielab::dynamics::vfex1();
//     Lielab::domain::SO y0(3);
//     Lielab::domain::CompositeManifold M0{y0};

//     Lielab::domain::CompositeAlgebra out = vf(2.0, M0);

//     Lielab::domain::so dy = std::get<Lielab::domain::so>(out.space[0]);

//     CHECK(std::abs(dy(0,0) - 0.0) < TOL_FINE);
//     CHECK(std::abs(dy(0,1) - 2.0) < TOL_FINE);
//     CHECK(std::abs(dy(0,2) - 1.0) < TOL_FINE);
//     CHECK(std::abs(dy(1,0) + 2.0) < TOL_FINE);
//     CHECK(std::abs(dy(1,1) - 0.0) < TOL_FINE);
//     CHECK(std::abs(dy(1,2) + 4.0) < TOL_FINE);
//     CHECK(std::abs(dy(2,0) + 1.0) < TOL_FINE);
//     CHECK(std::abs(dy(2,1) - 4.0) < TOL_FINE);
//     CHECK(std::abs(dy(2,2) - 0.0) < TOL_FINE);
// }
