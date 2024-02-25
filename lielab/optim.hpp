#ifndef _LIELAB_OPTIM_H
#define _LIELAB_OPTIM_H

#include <functional>

#include <Eigen/Core>
#include <cmath>

#include <iomanip>
#include "utils.hpp"

namespace lielab
{
    namespace optim
    {
        class nlp : public BaseAlgorithm
        {
            public:

            int num_objective_evals = 0;
            int num_jacobian_evals = 0;
            int num_hessian_evals = 0;
            std::string message = "None.";
            double tolerance = 1e-6;

            double val_objective;
            Eigen::VectorXd val_jacobian;
            Eigen::MatrixXd val_hessian;
            
            Eigen::VectorXd x0;
            Eigen::VectorXd lower;
            Eigen::VectorXd upper;

            nlp()
            {

            }

            void init()
            {
                BaseAlgorithm::init();

                num_objective_evals = 0;
                num_jacobian_evals = 0;
                num_hessian_evals = 0;
            }

            void step()
            {
                BaseAlgorithm::step();
            }

        };

        class hnewton : public nlp
        {
            public:

            double fdt = 1e-6;
            lielab::domain::halie lower;
            lielab::domain::halie upper;

            double f;
            double fnext;
            lielab::domain::halie m;
            lielab::domain::halie mnext;

            size_t _dim;

            hnewton()
            {

            }

            void init(const lielab::domain::halie & m0, double fval)
            {
                nlp::init();
                
                this->_dim = m0.get_dimension();

                if (this->_dim != 1)
                {
                    throw lielab::utils::NotImplementedError("hnewton only implemented for dimensions = 1.");
                }

                this->fnext = fval;
                this->mnext = m0;
            }

            void step0()
            {
                this->f = this->fnext;
                this->m = this->mnext;

                lielab::domain::halie dm = 0.0*this->m;
                Eigen::VectorXd _dm(1);
                _dm << this->fdt;
                std::get<lielab::domain::rn>(dm.space[0]).set_vector(_dm);

                this->mnext = this->m + dm;
            }

            void step1()
            {
                const double fprime = (this->fnext - this->f)/(this->fdt);

                lielab::domain::halie dm = 0.0*this->m;
                Eigen::VectorXd _dm(1);
                _dm << -this->f / fprime;
                std::get<lielab::domain::rn>(dm.space[0]).set_vector(_dm);

                this->mnext = this->m + dm;

                dm = this->mnext - this->m;
                lielab::domain::rn err = std::get<lielab::domain::rn>(dm.space[0]);
                Eigen::VectorXd err2 = err._data;

                if (err2.norm() < tolerance)
                {
                    success = true;
                    algo_status = lielab::ALGO_STATUS::FINISHED;
                }
            }

            lielab::domain::halie operator()(const std::function<double(lielab::domain::halie)> & fun, const lielab::domain::halie & m0)
            {
                init(m0, fun(m0));

                while (algo_status == lielab::ALGO_STATUS::OK)
                {
                    step0();
                    this->fnext = fun(this->mnext);
                    step1();
                    this->fnext = fun(this->mnext);
                }

                this->m = this->mnext;

                return this->m;
            }
        };

        class opt_golden : public nlp
        {
            public:

            double tau = (std::sqrt(5.0)-1.0)/2.0;

            double _f1;
            double _f2;
            Eigen::VectorXd _X;
            Eigen::VectorXd _X1;
            Eigen::VectorXd _X2;
            Eigen::VectorXd _A;
            Eigen::VectorXd _B;

            opt_golden() : nlp()
            {

            }

            void init()
            {
                nlp::init();

                _X1 = lower + (1-tau)*(upper - lower);
                _X2 = lower + tau*(upper - lower);
                _A = lower;
                _B = upper;
            }

            void step()
            {
                nlp::step();

                if (_f1 < _f2)
                {
                    val_objective = _f1;
                }
                else
                {
                    val_objective = _f2;
                }
                
                if (_f1 < _f2)
                {
                    _B = _X2;
                    _X2 = _X1;
                    _X1 = _A + (1 - tau)*(_B - _A);
                }
                else
                {
                    _A = _X1;
                    _X1 = _X2;
                    _X2 = _A + tau*(_B - _A);
                }

                if (std::abs((_B - _A).norm()) < tolerance)
                {
                    success = true;
                    algo_status = lielab::ALGO_STATUS::FINISHED;
                }
            }

            Eigen::VectorXd operator()(double (*obj_fun)(Eigen::VectorXd))
            {
                init();

                while (algo_status == lielab::ALGO_STATUS::OK)
                {
                    _f1 = obj_fun(_X1);
                    num_objective_evals++;

                    _f2 = obj_fun(_X2);
                    num_objective_evals++;
                    step();
                }
                
                _f1 = obj_fun(_X1);
                num_objective_evals++;

                _f2 = obj_fun(_X2);
                num_objective_evals++;

                if (_f1 < _f2)
                {
                    val_objective = _f1;
                    _X = _X1;
                }
                else
                {
                    val_objective = _f2;
                    _X = _X2;
                }

                return _X;
            }

        };

        class search_linearx : public nlp
        {
            public:

            double fdx = 1e-6;
            double lower = -9e9;
            double upper = 9e9;

            double _dx;
            double _x;
            double _x1;
            double _x2;
            double _y1;
            double _y2;
            int k;
            int _lo;
            int _hi;

            search_linearx() : nlp()
            {

            }

            void init(double x0)
            {
                nlp::init();
                
                _dx = 0.0;
                _x = x0;
                _x1 = 0.0;
                _x2 = 0.0;

                k = 2;
                _lo = 0;
                _hi = 0;
            }

            double step(double x, double y)
            {
                nlp::step();

                if (k == 2)
                {
                    if (abs(y) < tolerance)
                    {
                        k = 8;
                    }
                    else
                    {
                        _x1 = x;
                        _y1 = y;
                        _dx = fdx + fdx * abs(x);
                        x = x + _dx;
                        k = 3;
                    }
                }
                else if (k == 3)
                {
                    if (std::abs(y) < tolerance)
                    {
                        k = 8;
                    }
                    else
                    {
                        _x2 = x;
                        _y2 = y;

                        if (_y1 == _y2)
                        {
                            k = 9;
                        }
                        else
                        {
                            x = _x1 - _y1*(_x2 - _x1)/(_y2 - _y1);
                            if (x < lower)
                            {
                                x = lower;
                                _lo += 1;
                            }
                            else if (x > upper)
                            {
                                x = upper;
                                _hi += 1;
                            }
                            else
                            {
                                _lo = 0;
                                _hi = 0;
                            }
                            k = 4;
                        }
                    }
                    
                }
                else
                {
                    if (_lo > 2)
                    {
                        k = 9;
                        message = "Lower search bound reached.";
                        algo_status = ALGO_STATUS::FINISHED;
                    }
                    else if (_hi > 2)
                    {
                        k = 9;
                        message = "Upper search bound reached.";
                        algo_status = ALGO_STATUS::FINISHED;
                    }
                    else if (abs(y) < tolerance)
                    {
                        k = 8;
                        message = "Searched reached answer within tolerance.";
                        algo_status = ALGO_STATUS::FINISHED;
                    }
                    else
                    {
                        _x1 = x;
                        _y1 = y;
                        x = x + _dx;
                        k = 3;
                    }
                }

                return x;
            }

            double operator()(std::function<double(double)> obj_fun, double x0)
            {
                init(x0);

                while (algo_status == lielab::ALGO_STATUS::OK)
                {
                    _x = step(_x, obj_fun(_x));
                }

                return _x;
            }
        };

    }
}

#endif
  