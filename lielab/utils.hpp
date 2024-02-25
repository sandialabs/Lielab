#ifndef _LIELAB_UTILS_H
#define _LIELAB_UTILS_H

#include <exception>
#include <optional>
#include <vector>

#include <Eigen/Core>

namespace lielab
{
    namespace utils
    {
        struct NotImplementedError : public std::exception
        {
            std::string s;
            NotImplementedError(std::string ss) : s(ss) {}
            ~NotImplementedError() throw () {}
            const char* what() const throw() { return s.c_str(); }
        };

        Eigen::VectorXd concatenate(std::vector<Eigen::VectorXd> vlist)
        {
            if (vlist.size() == 0)
            {
                Eigen::VectorXd out(0);
                return out;
            }

            Eigen::VectorXd out = vlist[0];
            for (int ii = 1; ii < vlist.size(); ii++)
            {
                Eigen::VectorXd temp(out.size() + vlist[ii].size());
                temp << out, vlist[ii];
                out = temp;
            }
            return out;
        }

        Eigen::VectorXcd concatenate(std::vector<Eigen::VectorXcd> vlist)
        {
            if (vlist.size() == 0)
            {
                Eigen::VectorXcd out(0);
                return out;
            }

            Eigen::VectorXcd out = vlist[0];
            for (int ii = 1; ii < vlist.size(); ii++)
            {
                Eigen::VectorXcd temp(out.size() + vlist[ii].size());
                temp << out, vlist[ii];
                out = temp;
            }
            return out;
        }

        template <typename T>
        std::optional<T> smartfind(const std::vector<T> list, const std::function<bool(T)> cond)
        {
            const auto it = std::find_if(list.begin(), list.end(), cond);
            if (it == list.end())
            {
                return std::nullopt;
            }
            return *it;
        }

        template <typename T>
        std::optional<size_t> smartfindind(const std::vector<T> list, const std::function<bool(T)> cond)
        {
            const auto it = std::find_if(list.begin(), list.end(), cond);
            if (it == list.end())
            {
                return std::nullopt;
            }
            return std::distance(list.begin(), it);
        }
    }
}

#endif
