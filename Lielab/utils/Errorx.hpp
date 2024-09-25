#ifndef LIELAB_UTILS_ERRORX_HPP_
#define LIELAB_UTILS_ERRORX_HPP_

struct Errorx : public std::exception
{
   std::string s;
   Errorx(std::string ss) : s(ss) {}
   ~Errorx() throw () {}
   const char* what() const throw() { return s.c_str(); }
};

#endif