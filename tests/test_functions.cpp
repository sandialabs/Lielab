// TODO: pair
// TODO: copair
// TODO: factorial
// TODO: precomputed factorial

TEST_CASE("Ad", "[Ad]")
{
    /*!
    * Tests the Ad function.
    */
    lielab::domain::so u(3);
    lielab::domain::so v(3);
    lielab::domain::so w(3);
    lielab::domain::so ansso(3);
    lielab::domain::SO Gso(3);

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
    Gso = lielab::functions::exp(v);

    // GuG^-1
    ansso = lielab::functions::Ad(Gso, u);
    truthso << 0, 0.841470984807896, 0,
              -0.841470984807897, 0, -0.540302305868140,
               0, 0.540302305868140, 0;
    
    assert_matrix(ansso.get_ados_representation(), truthso);

    // GvG^-1 = v when G = exp(v)
    ansso = lielab::functions::Ad(Gso, v);
    
    assert_matrix(ansso.get_ados_representation(), v.get_ados_representation());

    // GwG^-1
    ansso = lielab::functions::Ad(Gso, w);
    truthso << 0, -0.540302305868140, 0,
               0.540302305868140, 0, -0.841470984807897,
               0, 0.841470984807897, 0;
    
    assert_matrix(ansso.get_ados_representation(), truthso);
}

TEST_CASE("cayley1", "[cayley1]")
{
    /*!
    * Tests the cayley1 function
    */

    const lielab::domain::so rx({1.0, 0.0, 0.0});
    const lielab::domain::so ry({0.0, 1.0, 0.0});
    Eigen::MatrixXd ans(3, 3);

    // Values calculated by hand
    const lielab::domain::SO ex1 = lielab::functions::cayley1(rx);
    ans << 1.0, 0.0, 0.0,
           0.0, 0.6,-0.8,
           0.0, 0.8, 0.6;

    assert_matrix(ex1.get_ados_representation(), ans);

    const lielab::domain::SO ex2 = lielab::functions::cayley1(ry);
    ans << 0.6, 0.0, 0.8,
           0.0, 1.0, 0.0,
           -0.8, 0.0, 0.6;

    assert_matrix(ex2.get_ados_representation(), ans);
}

TEST_CASE("cayley2", "[functions]")
{
    /*!
    * Tests the cayley2 function.
    */

    const lielab::domain::so rx({1.0, 0.0, 0.0});
    const lielab::domain::so ry({0.0, 1.0, 0.0});
    Eigen::MatrixXd ans(3, 3);

    // Values calculated by hand
    const lielab::domain::SO ex1 = lielab::functions::cayley2(rx);
    ans << 1.0, 0.0, 0.0,
           0.0, 0.6, -0.8,
           0.0, 0.8, 0.6;
    
    assert_matrix(ex1.get_ados_representation(), ans);

    const lielab::domain::SO ex2 = lielab::functions::cayley2(rx + 2*ry);
    ans << 0.0, 0.0, 1.0,
           0.8, 0.6, 0.0,
           -0.6, 0.8, 0.0;
    
    assert_matrix(ex2.get_ados_representation(), ans);
}

TEST_CASE("cayley 1 and 2", "[functions]")
{
    /*!
    * Tests cayley1 and cayley2 together with known identities.
    */

    // Identity cayley1 = cayley2 for all basis elements
    const size_t dim = lielab::domain::so::basis(0,10).get_dimension();
    for (size_t ii = 0; ii < dim; ii++)
    {
        const lielab::domain::so g = lielab::domain::so::basis(ii, 10);
        assert_domain(lielab::functions::cayley1(g), lielab::functions::cayley2(g));
    }
}

// TODO: commutator

TEST_CASE("Killing", "[functions]")
{
    /*!
    * Tests the Killing function.
    */

    lielab::domain::so rx({1.0, 0.0, 0.0});
    lielab::domain::so ry({0.0, 1.0, 0.0});
    lielab::domain::so rz({0.0, 0.0, 1.0});

    CHECK(std::abs(lielab::functions::Killing(rx,rx) + 2.0) <= TOL_FINE);
    CHECK(std::abs(lielab::functions::Killing(rx,ry) + 0.0) <= TOL_FINE);
    CHECK(std::abs(lielab::functions::Killing(rx,rz) + 0.0) <= TOL_FINE);
    CHECK(std::abs(lielab::functions::Killing(ry,rx) + 0.0) <= TOL_FINE);
    CHECK(std::abs(lielab::functions::Killing(ry,ry) + 2.0) <= TOL_FINE);
    CHECK(std::abs(lielab::functions::Killing(ry,rz) + 0.0) <= TOL_FINE);
    CHECK(std::abs(lielab::functions::Killing(rz,rx) + 0.0) <= TOL_FINE);
    CHECK(std::abs(lielab::functions::Killing(rz,ry) + 0.0) <= TOL_FINE);
    CHECK(std::abs(lielab::functions::Killing(rz,rz) + 2.0) <= TOL_FINE);
}

TEST_CASE("Killingform", "[functions]")
{
    /*!
    * Tests the killing form function.
    */

    lielab::domain::so rx({1.0, 0.0, 0.0});

    Eigen::MatrixXd K = lielab::functions::Killingform(rx);
    Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(rx.get_dimension(), rx.get_dimension());

    CHECK(std::abs(K.trace() + 6) <= TOL_FINE);
    assert_matrix(K*K.inverse(), Id);

    lielab::domain::so so6 = lielab::domain::so::basis(0,6);

    K = lielab::functions::Killingform(so6);
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
    lielab::domain::so u(3);
    lielab::domain::so v(3);
    lielab::domain::so ansso(3);

    Eigen::VectorXd xx(3);
    xx << 1.0, 0.0, 0.0;
    Eigen::VectorXd yy(3);
    yy << 0.0, 1.0, 0.0;
    Eigen::MatrixXd truthso(3,3);

    u.set_vector(xx);
    v.set_vector(yy);

    ansso = lielab::functions::dcayley1inv(u, v);
    truthso << 0.0, 0.5, 1.0,
              -0.5, 0.0, 0.0,
              -1.0, 0.0, 0.0;

    assert_matrix(ansso.get_ados_representation(), truthso);

    ansso = lielab::functions::dcayley1inv(v, u);
    truthso << 0.0,-0.5, 0.0,
               0.5, 0.0,-1.0,
               0.0, 1.0, 0.0;

    assert_matrix(ansso.get_ados_representation(), truthso);
}

TEST_CASE("dexp", "[dexp]")
{
    /*!
    * Tests the dexp function.
    */
    lielab::domain::so u(3);
    lielab::domain::so v(3);
    lielab::domain::so ansso(3);

    Eigen::VectorXd xx(3);
    xx << 1.0, 0.0, 0.0;
    Eigen::VectorXd yy(3);
    yy << 0.0, 1.0, 0.0;
    Eigen::MatrixXd truthso(3,3);

    u.set_vector(xx);
    v.set_vector(yy);


    // order = 0
    ansso = lielab::functions::dexp(u, v, 0);
    truthso << 0.0, 0.0, 1.0,
               0.0, 0.0, 0.0,
              -1.0, 0.0, 0.0;

    assert_matrix(ansso.get_ados_representation(), truthso);

    // order = 1
    ansso = lielab::functions::dexp(u, v, 1);
    truthso << 0.0, -0.5, 1.0,
               0.5, 0.0, 0.0,
              -1.0, 0.0, 0.0;

    assert_matrix(ansso.get_ados_representation(), truthso);

    // order = 2
    ansso = lielab::functions::dexp(u, v, 2);
    truthso << 0.0, -0.5, 0.833333333333333,
               0.5, 0.0, 0.0,
              -0.833333333333333, 0.0, 0.0;

    assert_matrix(ansso.get_ados_representation(), truthso);

    // order = 3
    ansso = lielab::functions::dexp(u, v, 3);
    truthso << 0.0, -0.458333333333333, 0.833333333333333,
               0.458333333333333, 0.0, 0.0,
              -0.833333333333333, 0.0, 0.0;

    assert_matrix(ansso.get_ados_representation(), truthso);


    // order = 4
    ansso = lielab::functions::dexp(u, v, 4);
    truthso << 0.0, -0.458333333333333, 0.841666666666667,
               0.458333333333333, 0.0, 0.0,
              -0.841666666666667, 0.0, 0.0;

    assert_matrix(ansso.get_ados_representation(), truthso);

    // order = 8
    ansso = lielab::functions::dexp(u, v, 8);
    truthso << 0.0, -0.459697420634921, 0.841471009700176,
               0.459697420634921, 0.0, 0.0,
              -0.841471009700176, 0.0, 0.0;

    assert_matrix(ansso.get_ados_representation(), truthso);

    lielab::domain::rn x(4);
    lielab::domain::rn y(4);
    lielab::domain::rn ansrn(4);
    Eigen::MatrixXd truthrn(4,4);
    x.set_vector(xx);
    y.set_vector(yy);

    // default order
    ansrn = lielab::functions::dexp(x, y);
    truthrn << 0, 0, 0, 0,
               0, 0, 0, 1,
               0, 0, 0, 0,
               0, 0, 0, 0;

    assert_matrix(ansrn.get_ados_representation(), truthrn);

    // ridiculous order (checks abelian speedhack)
    ansrn = lielab::functions::dexp(x, y, 999999999);
    truthrn << 0, 0, 0, 0,
               0, 0, 0, 1,
               0, 0, 0, 0,
               0, 0, 0, 0;

    assert_matrix(ansrn.get_ados_representation(), truthrn);
}

TEST_CASE("dexpinv", "[dexpinv]")
{
    /*!
    * Tests the dexpinv function.
    */
    lielab::domain::so u(3);
    lielab::domain::so v(3);
    lielab::domain::so ansso(3);

    Eigen::VectorXd xx(3);
    xx << 1.0, 0.0, 0.0;
    Eigen::VectorXd yy(3);
    yy << 0.0, 1.0, 0.0;
    Eigen::MatrixXd truthso(3,3);

    u.set_vector(xx);
    v.set_vector(yy);


    // order = 0
    ansso = lielab::functions::dexpinv(u, v, 0);
    truthso << 0.0, 0.0, 0.0,
               0.0, 0.0, 0.0,
               0.0, 0.0, 0.0;

    assert_matrix(ansso.get_ados_representation(), truthso);

    // order = 1
    ansso = lielab::functions::dexpinv(u, v, 1);
    truthso << 0.0, 0.0, 1.0,
               0.0, 0.0, 0.0,
              -1.0, 0.0, 0.0;

    assert_matrix(ansso.get_ados_representation(), truthso);

    // order = 2
    ansso = lielab::functions::dexpinv(u, v, 2);
    truthso << 0.0, 0.5, 1.0,
              -0.5, 0.0, 0.0,
              -1.0, 0.0, 0.0;

    assert_matrix(ansso.get_ados_representation(), truthso);

    // order = 4
    ansso = lielab::functions::dexpinv(u, v, 4);
    truthso << 0.0, 0.5, 0.916666666666667,
              -0.5, 0.0, 0.0,
              -0.916666666666667, 0.0, 0.0;

    assert_matrix(ansso.get_ados_representation(), truthso);

    // order = 12
    ansso = lielab::functions::dexpinv(u, v, 12);
    truthso << 0.0, 0.500000000000000, 0.915243861398375,
              -0.500000000000000, 0.0, 0.0,
              -0.915243861398375, 0.0, 0.0;

    assert_matrix(ansso.get_ados_representation(), truthso);

    lielab::domain::rn x(4);
    lielab::domain::rn y(4);
    lielab::domain::rn ansrn(4);
    Eigen::MatrixXd truthrn(4,4);
    x.set_vector(xx);
    y.set_vector(yy);

    // default order
    ansrn = lielab::functions::dexpinv(x, y);
    truthrn << 0, 0, 0, 0,
               0, 0, 0, 1,
               0, 0, 0, 0,
               0, 0, 0, 0;

    assert_matrix(ansrn.get_ados_representation(), truthrn);

    // ridiculous order (checks abelian speedhack)
    ansrn = lielab::functions::dexpinv(x, y, 999999999);
    truthrn << 0, 0, 0, 0,
               0, 0, 0, 1,
               0, 0, 0, 0,
               0, 0, 0, 0;

    assert_matrix(ansrn.get_ados_representation(), truthrn);
}
