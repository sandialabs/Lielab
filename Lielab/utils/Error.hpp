#ifndef LIELAB_UTILS_ERROR_HPP
#define LIELAB_UTILS_ERROR_HPP

#include <exception>
#include <string>

namespace Lielab::utils
{

struct Error : public std::exception
{
   std::string s;
   Error(std::string ss);
   ~Error() noexcept;
   const char* what() const noexcept override;
};

struct NotImplementedError : public std::exception
{
    std::string s;
    NotImplementedError(std::string ss);
    ~NotImplementedError() noexcept;
    const char* what() const noexcept override;
};

struct InputError : public std::exception
{
    std::string s;
    InputError(std::string ss);
    ~InputError() noexcept;
    const char* what() const noexcept override;
};

}

#endif
