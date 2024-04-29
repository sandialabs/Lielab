#ifndef _LIELAB_DOMAIN_EXTRA_HPP
#define _LIELAB_DOMAIN_EXTRA_HPP

struct IndexError : public std::exception
{
   std::string s;
   IndexError(std::string ss) : s(ss) {}
   ~IndexError() throw () {}
   const char* what() const throw() { return s.c_str(); }
};

struct NotImplementedError : public std::exception
{
   std::string s;
   NotImplementedError(std::string ss) : s(ss) {}
   ~NotImplementedError() throw () {}
   const char* what() const throw() { return s.c_str(); }
};

struct SizeError : public std::exception
{
   std::string s;
   SizeError(std::string ss) : s(ss) {}
   ~SizeError() throw () {}
   const char* what() const throw() { return s.c_str(); }
};

namespace Lielab
{
/*!
 * The domain namespace. Houses all functions and classes related to domains available in lielab.
 * 
 * The naming convention used within the domain namespace is different than the standard naming
 * convention used elsewhere. Typically, classes are in PascalCase and functions are in lowercase.
 * Somewhat deviating from this standard, lie algebras are in lowercase and Lie Groups are in
 * UpperCase. This distinction is made to follow mathematical convections. For example, 
 * su is the special unitary lie algebra and SU is the associated Special Unitary group.
 */

namespace domain
{

}
}

#endif
