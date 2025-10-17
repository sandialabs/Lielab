#include <Lielab.hpp>

#include <catch2/catch_all.hpp>

TEST_CASE("test_Coefficients", "[integrate]")
{
    using namespace Lielab::integrate;

    const int start = static_cast<int>(Coefficients::FE1);
    const int last = static_cast<int>(Coefficients::BE1);

    for (int intm = start; intm <= last; intm++)
    {
        const Coefficients method = static_cast<Coefficients>(intm);
        const auto [A, b, bhat, c, e, order, stages, variable, implicit] = get_butcher_tableau(method);

        CHECK(A.rows() == stages);
        CHECK(b.size() == stages);

        if (method != Coefficients::RKV65e)
        {
            CHECK_THAT(b.sum(), Catch::Matchers::WithinULP(1.0, 1));
        }

        if (variable)
        {
            CHECK(bhat.size() == stages);

            if (method != Coefficients::RKV87e)
            {
                CHECK_THAT(bhat.sum(), Catch::Matchers::WithinULP(1.0, 1));
            }
        }
        
        CHECK(c.size() == stages);
        CHECK(e.size() == stages);
    }
}