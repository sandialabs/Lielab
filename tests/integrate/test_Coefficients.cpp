#include <Lielab.hpp>

#include <catch2/catch_all.hpp>

TEST_CASE("test_RungeKuttaCoefficients", "[integrate]")
{
    using namespace Lielab::integrate;

    const int start = static_cast<int>(RungeKuttaCoefficients::FE1);
    const int last = static_cast<int>(RungeKuttaCoefficients::Lobatto3A6);

    for (int intm = start; intm <= last; intm++)
    {
        const RungeKuttaCoefficients method = static_cast<RungeKuttaCoefficients>(intm);
        const auto [A, b, bhat, c, e, order, stages, variable, implicit] = get_butcher_tableau(method);

        // TODO: Check the tableau is fully filled out
        for (int ii = 0; ii < stages; ii++)
        {
            for (int jj = 0; jj < stages; jj++)
            {
                INFO("Checking A for nans of " + std::to_string(static_cast<int>(method)));
                CHECK(!std::isnan(A[ii][jj]));
            }

            INFO("Checking b for nans of " + std::to_string(static_cast<int>(method)));
            CHECK(!std::isnan(b[ii]));

            if (variable)
            {
                INFO("Checking bhat for nans of " + std::to_string(static_cast<int>(method)));
                CHECK(!std::isnan(bhat[ii]));

                INFO("Checking e for nans of " + std::to_string(static_cast<int>(method)));
                CHECK(!std::isnan(e[ii]));
            }

            INFO("Checking c for nans of " + std::to_string(static_cast<int>(method)));
            CHECK(!std::isnan(c[ii]));
        }

        if (method != RungeKuttaCoefficients::RKV65e)
        {
            INFO("Checking sum(b) of " + std::to_string(static_cast<int>(method)));
            // CHECK_THAT(b.sum(), Catch::Matchers::WithinULP(1.0, 1));
        }

        if (variable)
        {
            if (method != RungeKuttaCoefficients::RKV87e)
            {
                INFO("Checking sum(bhat) of " + std::to_string(static_cast<int>(method)));
                // CHECK_THAT(bhat.sum(), Catch::Matchers::WithinULP(1.0, 1));
            }
        }
    }
}