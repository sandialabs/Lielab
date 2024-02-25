#ifndef _LIELAB_H
#define _LIELAB_H

#include <Eigen/Core>
#include <vector>

// using Eigen::VectorXf;

// Global variable definitions.
const std::complex<double> i(0.0,1.0);

namespace lielab
{

    enum class ALGO_STATUS {OK, MAXITER, FINISHED};
    class BaseAlgorithm
    {
        public:

        size_t iterations = 0;
        int max_iterations = -1;
        bool success = false;
        ALGO_STATUS algo_status = ALGO_STATUS::OK;

        BaseAlgorithm()
        {

        }

        void init()
        {
            iterations = 0;
            success = false;
            algo_status = ALGO_STATUS::OK;
        }

        void step()
        {
            iterations++;
            if (max_iterations > 0 && iterations >= max_iterations)
            {
                algo_status = ALGO_STATUS::MAXITER;
            }
        }
    };
    
}

#endif
