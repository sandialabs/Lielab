# Lielab

[![Conan Center](https://img.shields.io/conan/v/lielab)](https://conan.io/center/recipes/lielab)
[![PyPI - Version](https://img.shields.io/pypi/v/lielab)](https://pypi.org/project/lielab/)

C++ code for numerical algorithms and processes defined on Lie-type domains, including Lie algebras and Lie groups. Also fully usable from Python with a set of wrappers.

## User Installation

Source code and prebuilt binaries are available on various artifact hosting websites. Try these pre-built binaries first.

### Python

Install Lielab with pip

```
pip install lielab
```

### C++

Add Lielab to another project by adding it to the conanfile.txt and include with

```
#include <Lielab.hpp>
```

Alternatively, the header files can be included directly into other projects. Eigen must be made available in this case.

## Developer Installation

### Conan

Assuming Python is installed to your system already, the following block will generate the binary file for the Lielab Python bindings

```
conan build . -o with_python=True
```

This will build the pyd file required for Python, and also skip building the pure C++ test cases. The Python wrapper can then be installed with

```
cd python
pip install -e .
```

### CMake

Lielab is template header-only and the only other hard requirement is Eigen. The easiest way to make this available is installing it through Conan. This allows you to build and run Lielab with a minimal set of requirements. But building straight from the CMakeLists is also possible. Pybind11, Catch2, and Eigen3 targets need to be made available. This can be done from another CMakeLists. Alternatively, acquire these codes yourself, add them to the `/include/` folder, and uncomment the following lines from the top level Lielab CMakeLists

```
add_subdirectory(include/Catch2)
add_subdirectory(include/Eigen3)
add_subdirectory(include/pybind11)
```

from the top level Lielab CMakeLists. Then, build the tools as usual with CMake

```
cmake .. -DLIELAB_BUILD_TESTS=False
cmake --build . --target cppLielab
```

This will build the pyd file required for Python, and also skip building the pure C++ test cases. The Python wrapper can then be installed with

```
cd python
pip install -e .
```

## Citation

Find this repo useful?

```
@misc{Lielab,
  author = {Sparapany, Michael J.},
  title = {Lielab: Numerical Lie-theory in C++ and Python},
  year = {2024},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/sandialabs/lielab}}
}
```

## Dev Note

Using the Lielab Pybind11 wrapper to communicate to a heavier weight C++ code also using Lielab can be a very effective way at developing relatively performant C++ code while still retaining the usability of Python. When doing this, error messages thrown can be very unhelpful. The most common one may look like:

```
$ python myexample.py
Traceback (most recent call last):
  File "myexample.py", line <number>, in <module>
outputfromlielabcpp = myprojectusinglielabcpp.solve(inputfromlielabpython)
TypeError: solve(): incompatible function arguments. The following argument types are supported:
    1. (lielabinput: lielab::domain::so) -> [othercode::otherdatatype]

Invoked with: <lielab::domain::so object at <memory location>>
```

This error is thrown when calling some C++ code from a Python package and passing Lielab objects as an argument.

### Fix

Make sure both the C++ binary and the Lielab Python wrapper are compiled with the **same CMake flags**. Often times this error is caused when the Lielab Python wrapper was compiled as a `Release`, or `-DCMAKE_BUILD_TYPE=Release` and the other C++ code was compiled as a `Debug`, or `-DCMAKE_BUILD_TYPE=Debug`. Both need to be set to either `Release` or `Debug`. Only release versions of Lielab are included on PyPi so you may be required to compile your own set of Pybind11 wrappers for easier debugging.

As a rule of thumb, do not use the pre-compiled versions of Python Lielab if you intend to link it to another C++ project. Compile _both_ the C++ project and Python Lielab using the same platform and compiler to minimize issues like this.
