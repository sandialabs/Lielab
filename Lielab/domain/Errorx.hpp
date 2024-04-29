#ifndef LIELAB_ERRORX_
#define LIELAB_ERRORX_

struct Errorx : public std::exception
{
   std::string s;
   Errorx(std::string ss) : s(ss) {}
   ~Errorx() throw () {}
   const char* what() const throw() { return s.c_str(); }
};

#endif
