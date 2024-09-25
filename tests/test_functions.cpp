// TODO: pair
// TODO: copair
// TODO: factorial
// TODO: precomputed factorial

TEST_CASE("Ad", "[Ad]")
{
    /*!
    * Tests the Ad function.
    */
    Lielab::domain::so u(3);
    Lielab::domain::so v(3);
    Lielab::domain::so w(3);
    Lielab::domain::so ansso(3);
    Lielab::domain::SO Gso(3);

    Eigen::VectorXd xx(3);
    xx << 1.0, 0.0, 0.0;
    Eigen::VectorXd yy(3);
    yy << 0.0, 1.0, 0.0;
    Eigen::VectorXd zz(3);
    zz << 0.0, 0.0, 1.0;
    Eigen::MatrixXd truthso(3,3);

    u.set_vector(xx);
    v.set_vector(yy);
    w.set_vector(zz);
    Gso = Lielab::functions::exp(v);

    // GuG^-1
    ansso = Lielab::functions::Ad(Gso, u);
    truthso << 0, 0.841470984807896, 0,
              -0.841470984807897, 0, -0.540302305868140,
               0, 0.540302305868140, 0;
    
    assert_matrix(ansso.get_matrix(), truthso);

    // GvG^-1 = v when G = exp(v)
    ansso = Lielab::functions::Ad(Gso, v);
    
    assert_matrix(ansso.get_matrix(), v.get_matrix());

    // GwG^-1
    ansso = Lielab::functions::Ad(Gso, w);
    truthso << 0, -0.540302305868140, 0,
               0.540302305868140, 0, -0.841470984807897,
               0, 0.841470984807897, 0;
    
    assert_matrix(ansso.get_matrix(), truthso);
}

TEST_CASE("cayley1", "[cayley1]")
{
    /*!
    * Tests the cayley1 function
    */

    const Lielab::domain::so rx({1.0, 0.0, 0.0});
    const Lielab::domain::so ry({0.0, 1.0, 0.0});
    Eigen::MatrixXd ans(3, 3);

    // Values calculated by hand
    const Lielab::domain::SO ex1 = Lielab::functions::cayley1(rx);
    ans << 1.0, 0.0, 0.0,
           0.0, 0.6,-0.8,
           0.0, 0.8, 0.6;

    assert_matrix(ex1.get_matrix(), ans);

    const Lielab::domain::SO ex2 = Lielab::functions::cayley1(ry);
    ans << 0.6, 0.0, 0.8,
           0.0, 1.0, 0.0,
           -0.8, 0.0, 0.6;

    assert_matrix(ex2.get_matrix(), ans);
}

TEST_CASE("cayley2", "[functions]")
{
    /*!
    * Tests the cayley2 function.
    */

    const Lielab::domain::so rx({1.0, 0.0, 0.0});
    const Lielab::domain::so ry({0.0, 1.0, 0.0});
    Eigen::MatrixXd ans(3, 3);

    // Values calculated by hand
    const Lielab::domain::SO ex1 = Lielab::functions::cayley2(rx);
    ans << 1.0, 0.0, 0.0,
           0.0, 0.6, -0.8,
           0.0, 0.8, 0.6;
    
    assert_matrix(ex1.get_matrix(), ans);

    const Lielab::domain::SO ex2 = Lielab::functions::cayley2(rx + 2*ry);
    ans << 0.0, 0.0, 1.0,
           0.8, 0.6, 0.0,
           -0.6, 0.8, 0.0;
    
    assert_matrix(ex2.get_matrix(), ans);
}

TEST_CASE("cayley 1 and 2", "[functions]")
{
    /*!
    * Tests cayley1 and cayley2 together with known identities.
    */

    // Identity cayley1 = cayley2 for all basis elements
    const size_t dim = Lielab::domain::so::basis(0,10).get_dimension();
    for (size_t ii = 0; ii < dim; ii++)
    {
        const Lielab::domain::so g = Lielab::domain::so::basis(ii, 10);
        assert_domain(Lielab::functions::cayley1(g), Lielab::functions::cayley2(g));
    }
}

// TODO: commutator

TEST_CASE("Killing", "[functions]")
{
    /*!
    * Tests the Killing function.
    */

    Lielab::domain::so rx({1.0, 0.0, 0.0});
    Lielab::domain::so ry({0.0, 1.0, 0.0});
    Lielab::domain::so rz({0.0, 0.0, 1.0});

    CHECK(std::abs(Lielab::functions::Killing(rx,rx) + 2.0) <= TOL_FINE);
    CHECK(std::abs(Lielab::functions::Killing(rx,ry) + 0.0) <= TOL_FINE);
    CHECK(std::abs(Lielab::functions::Killing(rx,rz) + 0.0) <= TOL_FINE);
    CHECK(std::abs(Lielab::functions::Killing(ry,rx) + 0.0) <= TOL_FINE);
    CHECK(std::abs(Lielab::functions::Killing(ry,ry) + 2.0) <= TOL_FINE);
    CHECK(std::abs(Lielab::functions::Killing(ry,rz) + 0.0) <= TOL_FINE);
    CHECK(std::abs(Lielab::functions::Killing(rz,rx) + 0.0) <= TOL_FINE);
    CHECK(std::abs(Lielab::functions::Killing(rz,ry) + 0.0) <= TOL_FINE);
    CHECK(std::abs(Lielab::functions::Killing(rz,rz) + 2.0) <= TOL_FINE);
}

TEST_CASE("Killingform", "[functions]")
{
    /*!
    * Tests the killing form function.
    */

    Lielab::domain::so rx({1.0, 0.0, 0.0});

    Eigen::MatrixXd K = Lielab::functions::Killingform(rx);
    Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(rx.get_dimension(), rx.get_dimension());

    CHECK(std::abs(K.trace() + 6) <= TOL_FINE);
    assert_matrix(K*K.inverse(), Id);

    Lielab::domain::so so6 = Lielab::domain::so::basis(0,6);

    K = Lielab::functions::Killingform(so6);
    Id = Eigen::MatrixXd::Identity(so6.get_dimension(), so6.get_dimension());

    CHECK(std::abs(K.trace() + 120) <= TOL_FINE);
    assert_matrix(K*K.inverse(), Id);
}

// TODO: ad
// TODO: coAd
// TODO: coad
// TODO: exp
// TODO: log
// TODO: bernoulli

TEST_CASE("dcayley1inv", "[dcayley1inv]")
{
    /*!
    * Tests the inverse of the dcayley1 function
    */
    Lielab::domain::so u(3);
    Lielab::domain::so v(3);
    Lielab::domain::so ansso(3);

    Eigen::VectorXd xx(3);
    xx << 1.0, 0.0, 0.0;
    Eigen::VectorXd yy(3);
    yy << 0.0, 1.0, 0.0;
    Eigen::MatrixXd truthso(3,3);

    u.set_vector(xx);
    v.set_vector(yy);

    ansso = Lielab::functions::dcayley1inv(u, v);
    truthso << 0.0, 0.5, 1.0,
              -0.5, 0.0, 0.0,
              -1.0, 0.0, 0.0;

    assert_matrix(ansso.get_matrix(), truthso);

    ansso = Lielab::functions::dcayley1inv(v, u);
    truthso << 0.0,-0.5, 0.0,
               0.5, 0.0,-1.0,
               0.0, 1.0, 0.0;

    assert_matrix(ansso.get_matrix(), truthso);
}

TEST_CASE("dexp_numerical", "[dexp_numerical]")
{
    /*!
    * Tests the dexp_numerical function.
    */
    Lielab::domain::so u(3);
    Lielab::domain::so v(3);
    Lielab::domain::so ansso(3);

    Eigen::VectorXd xx(3);
    xx << 1.0, 0.0, 0.0;
    Eigen::VectorXd yy(3);
    yy << 0.0, 1.0, 0.0;
    Eigen::MatrixXd truthso(3,3);

    u.set_vector(xx);
    v.set_vector(yy);


    // order = 0
    ansso = Lielab::functions::dexp_numerical(u, v, 0);
    truthso << 0.0, 0.0, 1.0,
               0.0, 0.0, 0.0,
              -1.0, 0.0, 0.0;

    assert_matrix(ansso.get_matrix(), truthso);

    // order = 1
    ansso = Lielab::functions::dexp_numerical(u, v, 1);
    truthso << 0.0, -0.5, 1.0,
               0.5, 0.0, 0.0,
              -1.0, 0.0, 0.0;

    assert_matrix(ansso.get_matrix(), truthso);

    // order = 2
    ansso = Lielab::functions::dexp_numerical(u, v, 2);
    truthso << 0.0, -0.5, 0.833333333333333,
               0.5, 0.0, 0.0,
              -0.833333333333333, 0.0, 0.0;

    assert_matrix(ansso.get_matrix(), truthso);

    // order = 3
    ansso = Lielab::functions::dexp_numerical(u, v, 3);
    truthso << 0.0, -0.458333333333333, 0.833333333333333,
               0.458333333333333, 0.0, 0.0,
              -0.833333333333333, 0.0, 0.0;

    assert_matrix(ansso.get_matrix(), truthso);


    // order = 4
    ansso = Lielab::functions::dexp_numerical(u, v, 4);
    truthso << 0.0, -0.458333333333333, 0.841666666666667,
               0.458333333333333, 0.0, 0.0,
              -0.841666666666667, 0.0, 0.0;

    assert_matrix(ansso.get_matrix(), truthso);

    // order = 8
    ansso = Lielab::functions::dexp_numerical(u, v, 8);
    truthso << 0.0, -0.459697420634921, 0.841471009700176,
               0.459697420634921, 0.0, 0.0,
              -0.841471009700176, 0.0, 0.0;

    assert_matrix(ansso.get_matrix(), truthso);

    Lielab::domain::rn x(4);
    Lielab::domain::rn y(4);
    Lielab::domain::rn ansrn(4);
    Eigen::MatrixXd truthrn(4,4);
    x.set_vector(xx);
    y.set_vector(yy);

    // default order
    ansrn = Lielab::functions::dexp_numerical(x, y);
    truthrn << 0, 0, 0, 0,
               0, 0, 0, 1,
               0, 0, 0, 0,
               0, 0, 0, 0;

    assert_matrix(ansrn.get_matrix(), truthrn);

    // ridiculous order (checks abelian speedhack)
    ansrn = Lielab::functions::dexp_numerical(x, y, 999999999);
    truthrn << 0, 0, 0, 0,
               0, 0, 0, 1,
               0, 0, 0, 0,
               0, 0, 0, 0;

    assert_matrix(ansrn.get_matrix(), truthrn);
}

TEST_CASE("dexpinv_numerical", "[dexpinv_numerical]")
{
    /*!
    * Tests the dexpinv_numerical function.
    */
    Lielab::domain::so u(3);
    Lielab::domain::so v(3);
    Lielab::domain::so ansso(3);

    Eigen::VectorXd xx(3);
    xx << 1.0, 0.0, 0.0;
    Eigen::VectorXd yy(3);
    yy << 0.0, 1.0, 0.0;
    Eigen::MatrixXd truthso(3,3);

    u.set_vector(xx);
    v.set_vector(yy);


    // order = 0
    ansso = Lielab::functions::dexpinv_numerical(u, v, 0);
    truthso << 0.0, 0.0, 1.0,
               0.0, 0.0, 0.0,
              -1.0, 0.0, 0.0;

    assert_matrix(ansso.get_matrix(), truthso);

    // order = 1
    ansso = Lielab::functions::dexpinv_numerical(u, v, 1);
    truthso << 0.0, 0.5, 1.0,
              -0.5, 0.0, 0.0,
              -1.0, 0.0, 0.0;

    assert_matrix(ansso.get_matrix(), truthso);

    // order = 3
    ansso = Lielab::functions::dexpinv_numerical(u, v, 3);
    truthso << 0.0, 0.5, 0.916666666666667,
              -0.5, 0.0, 0.0,
              -0.916666666666667, 0.0, 0.0;

    assert_matrix(ansso.get_matrix(), truthso);

    // order = 11
    ansso = Lielab::functions::dexpinv_numerical(u, v, 11);
    truthso << 0.0, 0.500000000000000, 0.915243861398375,
              -0.500000000000000, 0.0, 0.0,
              -0.915243861398375, 0.0, 0.0;

    assert_matrix(ansso.get_matrix(), truthso);

    Lielab::domain::rn x(4);
    Lielab::domain::rn y(4);
    Lielab::domain::rn ansrn(4);
    Eigen::MatrixXd truthrn(4,4);
    x.set_vector(xx);
    y.set_vector(yy);

    // default order
    ansrn = Lielab::functions::dexpinv_numerical(x, y);
    truthrn << 0, 0, 0, 0,
               0, 0, 0, 1,
               0, 0, 0, 0,
               0, 0, 0, 0;

    assert_matrix(ansrn.get_matrix(), truthrn);

    // ridiculous order (checks abelian speedhack)
    ansrn = Lielab::functions::dexpinv_numerical(x, y, 999999999);
    truthrn << 0, 0, 0, 0,
               0, 0, 0, 1,
               0, 0, 0, 0,
               0, 0, 0, 0;

    assert_matrix(ansrn.get_matrix(), truthrn);
}
