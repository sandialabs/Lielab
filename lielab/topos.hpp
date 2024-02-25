#ifndef _LIELAB_TOPOS_H
#define _LIELAB_TOPOS_H

#include <exception>
#include <fstream>
#include <functional>

#include <Eigen/Core>

#include "src/topos/rk_methods.hpp"
#include "src/topos/functions/cayley1.hpp"
#include "src/topos/functions/exp.hpp"

#include "domain"
#include "functions"
#include "optim.hpp"

struct InputError : public std::exception
{
   std::string s;
   InputError(std::string ss) : s(ss) {}
   ~InputError() throw () {}
   const char* what() const throw() { return s.c_str(); }
};

namespace lielab
{
    namespace topos
    {
        enum class COORDINATES {EXPONENTIAL, CAYLEY1};
        enum class RKTYPE {RKTYPE_NONE, RKTYPE_EXPLICIT, RKTYPE_IMPLICIT};

        lielab::domain::halie dcayley1inv(const lielab::domain::halie & a, const lielab::domain::halie & b)
        {
            /*!
            * halie dcayley1inv overload
            */

            lielab::domain::halie out;

            for (int ii = 0; ii < a.space.size(); ii++)
            {
                const size_t ind = a.space[ii].index();
                if (ind == lielab::domain::halie::INDEX_rn)
                {
                    out.space.push_back(lielab::functions::dcayley1inv(std::get<lielab::domain::rn>(a.space[ii]),
                                                                       std::get<lielab::domain::rn>(b.space[ii])));
                }
                else if (ind == lielab::domain::halie::INDEX_so)
                {
                    out.space.push_back(lielab::functions::dcayley1inv(std::get<lielab::domain::so>(a.space[ii]),
                                                                       std::get<lielab::domain::so>(b.space[ii])));
                }
                else if (ind == lielab::domain::halie::INDEX_sp)
                {
                    out.space.push_back(lielab::functions::dcayley1inv(std::get<lielab::domain::sp>(a.space[ii]),
                                                                       std::get<lielab::domain::sp>(b.space[ii])));
                }
                else if (ind == lielab::domain::halie::INDEX_su)
                {
                    out.space.push_back(lielab::functions::dcayley1inv(std::get<lielab::domain::su>(a.space[ii]),
                                                                       std::get<lielab::domain::su>(b.space[ii])));
                }
            }

            return out;
        }

        lielab::domain::halie dexpinv(const lielab::domain::halie & a, const lielab::domain::halie & b, const size_t order = 5)
        {
            /*!
            * halie dexpinv overload
            */

            lielab::domain::halie out;

            for (int ii = 0; ii < a.space.size(); ii++)
            {
                const size_t ind = a.space[ii].index();
                if (ind == lielab::domain::halie::INDEX_gl)
                {
                    out.space.push_back(lielab::functions::dexpinv(std::get<lielab::domain::gl>(a.space[ii]),
                                                                   std::get<lielab::domain::gl>(b.space[ii]), order));
                }
                else if (ind == lielab::domain::halie::INDEX_rn)
                {
                    out.space.push_back(lielab::functions::dexpinv(std::get<lielab::domain::rn>(a.space[ii]),
                                                                   std::get<lielab::domain::rn>(b.space[ii]), order));
                }
                else if (ind == lielab::domain::halie::INDEX_so)
                {
                    out.space.push_back(lielab::functions::dexpinv(std::get<lielab::domain::so>(a.space[ii]),
                                                                   std::get<lielab::domain::so>(b.space[ii]), order));
                }
                else if (ind == lielab::domain::halie::INDEX_sp)
                {
                    out.space.push_back(lielab::functions::dexpinv(std::get<lielab::domain::sp>(a.space[ii]),
                                                                   std::get<lielab::domain::sp>(b.space[ii]), order));
                }
                else if (ind == lielab::domain::halie::INDEX_su)
                {
                    out.space.push_back(lielab::functions::dexpinv(std::get<lielab::domain::su>(a.space[ii]),
                                                                   std::get<lielab::domain::su>(b.space[ii]), order));
                }
            }

            return out;
        }

        /*!
        * The IntegralCurve class. Data structure for data returned by the Flow class.
        */
        class IntegralCurve
        {
            public:

            size_t chunk = 2048;
            size_t length = chunk;
            size_t num_eoms;

            int chs = 0;

            Eigen::VectorXd t;
            Eigen::MatrixXd y;

            IntegralCurve()
            {
                this->num_eoms = 0;
                this->t = Eigen::VectorXd::Zero(chunk);
                this->y = Eigen::MatrixXd::Zero(chunk, 0);
            }

            IntegralCurve(const IntegralCurve & other) :
            chunk(other.chunk),
            length(other.length),
            num_eoms(other.num_eoms),
            chs(other.chs),
            t(other.t),
            y(other.y)
            {
                /*!
                * Copy constructor for IntegralCurve
                */
            }

            IntegralCurve(const size_t num_eoms)
            {
                this->num_eoms = num_eoms;
                this->t = Eigen::VectorXd::Zero(chunk);
                this->y = Eigen::MatrixXd::Zero(chunk, num_eoms);
            }

            IntegralCurve & operator=(const IntegralCurve & other)
            {
                this->num_eoms = other.num_eoms;
                this->chunk = other.chunk;
                this->length = other.length;
                this->chs = other.chs;
                this->t = other.t;
                this->y = other.y;
                return *this;
            }

            void trim_chunk(const size_t last_index)
            {
                Eigen::VectorXd _t(last_index + 1);
                Eigen::MatrixXd _y(last_index + 1, num_eoms);

                _t = this->t.head(last_index + 1);
                _y = this->y.block(0, 0, last_index + 1, num_eoms);

                this->t = _t;
                this->y = _y;
                this->length = this->t.rows();
            }

            void add_chunk()
            {
                Eigen::VectorXd _t(length + chunk);
                Eigen::MatrixXd _y(length + chunk, num_eoms);

                _t.head(length) = this->t;
                _y.block(0, 0, length, num_eoms) = this->y;

                this->chs++;

                this->t = _t;
                this->y = _y;
                this->length = this->t.rows();
            }

            void set_data(const double & t, const Eigen::VectorXd & y, const int & index)
            {
                while (index >= length)
                {
                    this->add_chunk();
                }
                this->y.block(index, 0, 1, num_eoms) = y.transpose();
                this->t(index) = t;
            }

            void write_to_file(const std::string & name)
            {
                // TODO: This should be a function until utils but imports
                //       need to be cleaned up first

                // TODO: Probably just delete this now that we have a complete module

                std::ofstream outputfile(name);
                outputfile << std::setprecision(16);
                size_t n = this->t.rows();
                size_t num_eoms = this->y.cols();

                for (int i = 0; i < n; i++)
                {
                    outputfile << this->t(i) << ',';
                    for (int j = 0; j < num_eoms; j++)
                    {   
                        outputfile << this->y(i,j);
                        
                        if (j < num_eoms - 1)
                        {
                            outputfile << ",";
                        }
                    }
                    outputfile << "\n";
                }
                outputfile.close();
            }
        };

        class TSOutput
        {
            public:

            lielab::domain::hmlie low;
            lielab::domain::hmlie high;
            double error = -1;

            TSOutput() {}

            TSOutput(const lielab::domain::hmlie & low) :
            low(low)
            {

            }

            TSOutput(const lielab::domain::hmlie & low, const lielab::domain::hmlie & high, const double error) :
            low(low), high(high), error(error)
            {

            }
        };

        /*!
         * The TimeStepper class. Numerically integrates a vector field over a single time segment.
         */
        class TimeStepper
        {
            public:

            Eigen::MatrixXd A;
            Eigen::VectorXd B;
            Eigen::VectorXd Bhat;
            Eigen::VectorXd C;
            size_t n;
            size_t order;
            bool variable_step;

            RKTYPE RKT = RKTYPE::RKTYPE_NONE;

            /* from base algorithm */
            size_t iterations = 0;
            int max_iterations = -1;
            bool success = false;
            ALGO_STATUS algo_status = ALGO_STATUS::OK;

            TimeStepper() : TimeStepper(lielab::topos::RKMETHOD::RK45)
            {
                /*!
                * Creates a timestepper with default RK method RK45.
                */
                
            }

            TimeStepper(const RKMETHOD & rkmethod) :
            A(lielab::topos::RKMETHOD_to_A(rkmethod)),
            B(lielab::topos::RKMETHOD_to_B(rkmethod)),
            Bhat(lielab::topos::RKMETHOD_to_Bhat(rkmethod)),
            C(lielab::topos::RKMETHOD_to_C(rkmethod)),
            n(lielab::topos::RKMETHOD_to_n(rkmethod)),
            order(lielab::topos::RKMETHOD_to_order(rkmethod)),
            variable_step(lielab::topos::RKMETHOD_to_variable(rkmethod))
            {
                /*!
                * Creates a timestepper with a built-in RK method specified by RKMETHOD.
                */
                
            }

            void init()
            {
                this->iterations = 0;
                this->success = false;
                this->algo_status = ALGO_STATUS::OK;
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

        /*!
         * A Munthe-Kaas type timestepper.
         *
         * Let \f$M\f$ be a homogeneous manifold with left action \f$L : (G, M) \rightarrow M\f$,
         * let \f$t \in \mathbb{R}\f$ be time, and let \f$f : (t, M) \rightarrow \mathfrak{g}\f$ be a set of first-order ordinary
         * differential equations. A Munthe-Kaas method solves the
         * differential equations by solving the related set of equations:
         * 
         * \f{equation*}{
         * \frac{dU}{dt} = d\phi_U^{-1}(f(t, y))
         * \f}
         * 
         * where \f$ \phi : \mathfrak{g} \rightarrow G \f$ is a coordinate map and
         * \f$ d\phi_U^{-1} : \mathfrak{g} \rightarrow T \mathfrak{g} \f$. Given an \f$s\f$-stage
         * Runge-Kutta method, the iterative Munthe-Kaas method is summarized as
         * 
         * \f{eqnarray*}{
         * U_i &=& dt \sum_{j=1}^{i-1} a_{ij} K_j \\
         * y_i &=& L(\phi(U_i), y_{i-1}) \\
         * K_i &=& d\phi_U^{-1} f(t, y_i) \\
         * i &=& 1, \dots, s \\
         * V &=& dt \sum_{i=1}^s b_i K_i \\
         * y_f &=& L(\phi(V), y_0)
         * \f}
         * 
         * Author / Date: Sparapany / 2020
         */
        class MuntheKaas : public TimeStepper
        {
            public:
            RKTYPE RKT = RKTYPE::RKTYPE_EXPLICIT;

            std::function<lielab::domain::hmlie(const lielab::domain::halie, const lielab::domain::hmlie)> left = &lielab::functions::left_exp_default;
            std::function<lielab::domain::halie(const lielab::domain::halie, const lielab::domain::halie, int)> dphiinv = &lielab::topos::dexpinv;

            double _t0;
            lielab::domain::hmlie _y0;
            double _dt;

            double next_t;
            lielab::domain::hmlie next_y;

            lielab::domain::halie _dy;
            lielab::domain::halie _U;
            std::vector<lielab::domain::halie> _KK;

            MuntheKaas() : TimeStepper()
            {
                /*! \f{equation*}{ () \rightarrow MuntheKaas \f}
                * Instantiates a new MuntheKaas object.
                */

            }

            MuntheKaas(const lielab::topos::RKMETHOD rkmethod) : TimeStepper(rkmethod)
            {
                /*! \f{equation*}{ (RKMETHOD) \rightarrow MuntheKaas \f}
                * Instantiates a new MuntheKaas object with a specific RK method.
                */
            }

            void init(const double t0, const lielab::domain::hmlie & y0, const double dt, const lielab::domain::halie & dy0)
            {
                /*! \f{equation*}{ (\mathbb{R}, M, \mathbb{R}, \mathfrak{g}^{n}) \rightarrow () \f}
                * Initializes the MuntheKaas process.
                */

                TimeStepper::init();

                _t0 = t0;
                _y0 = y0;
                _dt = dt;
                _KK = std::vector<lielab::domain::halie>(n);
                _KK[0] = dy0; // [dy,dy]=0 therefore dexp(dy, dy) = dy
                _dy = dy0;
                _U = 0*dy0;

                // Set algo finish for 1st order methods.
                if (iterations >= n-1)
                {
                    this->algo_status = lielab::ALGO_STATUS::FINISHED;
                }
            }

            void step_0()
            {
                /*!
                * Advances the solution to current iteration and returns pair
                * (next_t, next_y) to be evaluated by the vectorfield.
                */

                _U *= 0;

                for (int jj = 0; jj < iterations+1; jj++)
                {
                    _U += _dt*A(iterations+1, jj)*_KK[jj];
                }

                next_y = this->left(_U, _y0);
                next_t = _t0 + _dt*C(iterations+1);
            }

            void set_dy(const lielab::domain::halie & dy)
            {
                /*!
                * Sets the current vectorfield value. Do this after step_0()
                * and before step_1().
                */
                _dy = dy;
            }

            void step_1()
            {
                /*!
                * Evaluates the \f$ d\phi^{-1}_u(\xi(s,y)) \f$ map.
                */
                _KK[iterations+1] = this->dphiinv(_U, _dy, order=n-1);

                TimeStepper::step();

                if (iterations == n-1)
                {
                    algo_status = lielab::ALGO_STATUS::FINISHED;
                }
            }

            TSOutput postprocess() const
            {
                lielab::domain::halie Ulow = 0*_KK[0];
                lielab::domain::halie Uhigh = 0*_KK[0];

                for (int ii = 0; ii < n; ii++)
                {
                    Ulow += _dt*B(ii)*_KK[ii];
                }

                const lielab::domain::hmlie low = this->left(Ulow, _y0);

                if (variable_step == true)
                {
                    for (int ii = 0; ii < n; ii++)
                    {
                        Uhigh += _dt*Bhat(ii)*_KK[ii];
                    }
                    
                    const lielab::domain::hmlie high = this->left(Uhigh, _y0);

                    const double error = (Uhigh.get_vector() - Ulow.get_vector()).norm();
                    return TSOutput(low, high, error);
                }
                return TSOutput(low);
            }

            TSOutput operator()(const std::function<lielab::domain::halie(double, lielab::domain::hmlie)> & vectorfield, const lielab::domain::hmlie & y0, const double t0, const double dt)
            {
                init(t0, y0, dt, vectorfield(t0, y0));
                
                if (RKT == RKTYPE::RKTYPE_EXPLICIT)
                {
                    while (algo_status == lielab::ALGO_STATUS::OK)
                    {
                        step_0();
                        set_dy(vectorfield(next_t, next_y));
                        step_1();
                    }
                }

                return postprocess();
            }
        };


        /*
        * The Flow class. Numerically calculates the integral curve of a vector field through y0 from t0 to tf.
        */
        class Flow : public BaseAlgorithm
        {
            public:
            double small = 0.5;
            double large = 2.0;
            double pessimist = 0.9;
            double accept = 1.2;
            double tol = 1e-6;
            double dt_min = 1e-16;
            double dt_max = 1.0;

            double dt = 0.02;
            double default_local = 10.0;
            double default_global = 10.0;

            bool variable_time_step = true;

            MuntheKaas stepper;
            lielab::optim::hnewton search;

            int new_exact = 1;
            int num_step = 8;

            double _dt;
            double _dt_temp;
            size_t tind;
            std::vector<double> _tspan;
            double _error_estimate;
            double _dt_new;
            double _eps;
            size_t _num_eoms;
            lielab::domain::hmlie _ynext;

            // Output variables
            IntegralCurve _out;

            void init(const std::vector<double> & tspan, const lielab::domain::hmlie & y0)
            {
                if (tspan.size() < 2)
                {
                    throw InputError("Length of tspan must be 2 (received " + std::to_string(tspan.size()) + ").");
                }

                if ((this->variable_time_step == true) && (this->stepper.variable_step == false))
                {
                    this->variable_time_step = false;
                }

                _dt = dt;
                _dt_temp = std::numeric_limits<double>::quiet_NaN();
                _tspan = tspan;
                tind = 1;
                _num_eoms = y0.serialize().size();

                _eps = std::numeric_limits<double>::epsilon();
                
                _ynext = y0;

                BaseAlgorithm::init();

                _out = IntegralCurve(_num_eoms);
                _out.t(0) = tspan.front();

                Eigen::VectorXd v1(_num_eoms);
                v1 = y0.serialize();

                _out.y.block(0, 0, 1, _num_eoms) = v1.transpose();
            }

            bool step0(const TSOutput & next)
            {
                bool accepted = false;

                if (variable_time_step)
                {
                    _dt_new = new_step_size(next.error, stepper.order);
                    if ((accept*tol - next.error) > 0)
                    {
                        accepted = true;
                        _ynext = next.high;
                    }
                    else
                    {
                        accepted = false;
                        _dt = _dt_new;
                    }
                }
                else
                {
                    accepted = true;
                    _ynext = next.low;
                }

                return accepted;
            }

            void step()
            {
                const size_t sz = _ynext.serialize().size();
                Eigen::VectorXd v1(sz);
                int pos = 0;
                v1 << _ynext.serialize();

                _out.set_data(_out.t(iterations) + _dt, v1, iterations + 1);

                if (variable_time_step == true)
                {
                    _dt = _dt_new;
                }

                // If a temporary time step is stored, reset it.
                if (!std::isnan(_dt_temp))
                {
                    _dt = _dt_temp;
                    _dt_temp = std::numeric_limits<double>::quiet_NaN();
                }

                // Check if the next time step will be over the next time in tspan
                if (_out.t(iterations + 1) + _dt > _tspan[tind])
                {
                    // Check to see if there are more tspan values. If there are, save the current dt
                    if (tind < _tspan.size()-1)
                    {
                        _dt_temp = _dt;
                    }

                    // Set the next dt to end exactly on the next tspan value
                    _dt = _tspan[tind] - _out.t(iterations + 1);

                    // Set up the next tind if there is one
                    if (tind < _tspan.size()-1)
                    {
                        tind += 1;
                    }
                }

                BaseAlgorithm::step();

                // If we're passed the last value in tspan, we're done.
                if (_out.t(iterations) + _dt > _tspan.back())
                {
                    algo_status = lielab::ALGO_STATUS::FINISHED;
                }

                // Stop integration if it grinds too slow.
                // TODO: Remove this? Handle in another way?
                if (_dt < 100*_eps)
                {
                    algo_status = lielab::ALGO_STATUS::FINISHED;
                }

            }

            void stepE(double event_val)
            {
                if (event_val < 0 || std::abs(event_val) < tol)
                {
                    algo_status = lielab::ALGO_STATUS::FINISHED;
                }
            }

            void postprocess()
            {
                _out.trim_chunk(iterations);
            }

            IntegralCurve operator()(std::function<lielab::domain::halie(double, lielab::domain::hmlie)> vectorfield, const std::vector<double> & tspan, const lielab::domain::hmlie & y0)
            {
                /*!
                * The main evaluation method for the Flow class.
                */

                // Initialize the Flow object
                init(tspan, y0);

                // Main loop. Run until algorithm says it's done
                while (algo_status == lielab::ALGO_STATUS::OK)
                {
                    bool accepted = false;

                    // This is the RK loop. It's the error control loop that will
                    // change step size until we have an acceptable step
                    while (accepted == false)
                    {
                        const TSOutput next = stepper(vectorfield, _ynext, _out.t(iterations), _dt);
                        accepted = step0(next);
                    }

                    // We have an acceptable step at this point. Step the solution forward by one
                    step();
                }
                postprocess();

                return IntegralCurve(_out);
            }

            IntegralCurve operator()(std::function<lielab::domain::halie(double, lielab::domain::hmlie)> vectorfield, const std::vector<double> & tspan, const lielab::domain::hmlie & y0, std::function<double(double, lielab::domain::hmlie)> event)
            {
                init(tspan, y0);

                bool accepted = false;

                TSOutput next;

                while (algo_status == lielab::ALGO_STATUS::OK)
                {
                    accepted = false;

                    while (accepted == false)
                    {
                        next = stepper(vectorfield, _ynext, _out.t(iterations), _dt);

                        if (event(_out.t(iterations), _ynext) > 0 && event(_out.t(iterations) + _dt, next.high) <= 0)
                        {
                            search.lower = lielab::domain::halie{lielab::domain::rn{tol}};
                            search.upper = lielab::domain::halie{lielab::domain::rn{_dt}};

                            auto sfun = [&](lielab::domain::halie sdt)
                            {
                                const double _sdt = std::get<lielab::domain::rn>(sdt.space[0])(0);
                                return event(_out.t(iterations) + _sdt, stepper(vectorfield, _ynext, _out.t(iterations), _sdt).high);
                            };

                            const lielab::domain::halie dt_guess = lielab::domain::halie{lielab::domain::rn{_dt/2.0}};
                            const lielab::domain::halie _temp_dt = search(sfun, dt_guess);
                            const double temp_dt = std::get<lielab::domain::rn>(_temp_dt.space[0])(0);
                            _dt = temp_dt;
                            next = stepper(vectorfield, _ynext, _out.t(iterations), _dt);
                        }

                        accepted = step0(next);
                    }
                    step();
                    stepE(event(_out.t(iterations), _ynext));
                }

                postprocess();

                return _out;
            }

            IntegralCurve operator()(std::function<Eigen::VectorXd(double, Eigen::VectorXd)> vectorfield, const std::vector<double> & tspan, const Eigen::VectorXd & y0)
            {
                /*!
                * Overloaded main evaluation method for non-MK integration.
                */

                const lielab::domain::RN Y0(y0);
                const lielab::domain::hmlie M0{Y0};

                std::function<lielab::domain::halie(double, lielab::domain::hmlie)> wrappedfunction = [vectorfield](const double t, const lielab::domain::hmlie & M)
                {
                    const lielab::domain::RN s = std::get<lielab::domain::RN>(M.space[0]);
                    const Eigen::VectorXd states = s.serialize();
                    const Eigen::VectorXd gradient = vectorfield(t, states);
                    const lielab::domain::rn dx(gradient);
                    return lielab::domain::halie{dx};
                };

                return this->operator()(wrappedfunction, tspan, M0);
            }

            IntegralCurve operator()(std::function<Eigen::VectorXd(double, Eigen::VectorXd)> vectorfield, const std::vector<double> & tspan, const Eigen::VectorXd & y0, std::function<double(double, Eigen::VectorXd)> event)
            {
                /*!
                * Overloaded evaluation method for non-MK integration with events.
                */

                const lielab::domain::RN Y0(y0);
                const lielab::domain::hmlie M0{Y0};

                std::function<lielab::domain::halie(double, lielab::domain::hmlie)> wrappedfunction = [vectorfield](const double t, const lielab::domain::hmlie & M)
                {
                    const lielab::domain::RN s = std::get<lielab::domain::RN>(M.space[0]);
                    const Eigen::VectorXd states = s.serialize();
                    const Eigen::VectorXd gradient = vectorfield(t, states);
                    const lielab::domain::rn dx(gradient);
                    return lielab::domain::halie{dx};
                };

                std::function<double(double, lielab::domain::hmlie)> wrappedevent = [event](const double t, const lielab::domain::hmlie & M)
                {
                    const lielab::domain::RN s = std::get<lielab::domain::RN>(M.space[0]);
                    const Eigen::VectorXd states = s.serialize();
                    return event(t, states);
                };

                return this->operator()(wrappedfunction, tspan, M0, wrappedevent);
            }

            double new_step_size(const double error_estimate, const size_t RKOrder) const
            {
                /*!
                * Evaluates, but does not set, a new step size based on order of RK method and error estimated.
                */

                if (!variable_time_step)
                {
                    return this->dt;
                }

                double dt_new = pessimist*std::pow((tol/error_estimate), 1.0/(1.0*RKOrder + 1.0))*_dt;
                dt_new = std::min(dt_new, _dt*large);
                dt_new = std::max(dt_new, _dt*small);
                dt_new = std::min(dt_new, dt_max);
                dt_new = std::max(dt_min, dt_new);
                return dt_new;
            }
        };
    }
}

#endif
