#include "Error.hpp"

#include <exception>
#include <string>

namespace Lielab::utils
{

Error::Error(std::string ss) : s(ss) {}
Error::~Error() noexcept {}
const char* Error::what() const noexcept { return s.c_str(); }

NotImplementedError::NotImplementedError(std::string ss) : s(ss) {}
NotImplementedError::~NotImplementedError() noexcept {}
const char* NotImplementedError::what() const noexcept { return s.c_str(); }

InputError::InputError(std::string ss) : s(ss) {}
InputError::~InputError() noexcept {}
const char* InputError::what() const noexcept { return s.c_str(); }

}
