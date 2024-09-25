from conan import ConanFile
from conan.tools.cmake import CMakeToolchain, CMakeDeps, CMake, cmake_layout
from conan.tools.build import check_min_cppstd
from conan.tools.scm import Version
from conan.errors import ConanInvalidConfiguration
import os

required_conan_version = ">=1.52.0"

class LielabConan(ConanFile):
    name = "lielab"
    homepage = "https://github.com/sandialabs/Lielab"
    settings = "os", "compiler", "build_type", "arch"

    no_copy_source = True

    def requirements(self):
        self.requires("eigen/3.4.0")
        self.requires("pybind11/2.12.0")
        self.requires("catch2/3.4.0")

    @property
    def _min_cppstd(self):
        return 20

    @property
    def _compilers_minimum_version(self):
        return {
            "gcc": "11",
            "clang": "12",
            "apple-clang": "13.1",
            "Visual Studio": "17",
            "msvc": "193",
        }

    def validate(self):
        if self.settings.compiler.cppstd:
            check_min_cppstd(self, self._min_cppstd)
        minimum_version = self._compilers_minimum_version.get(str(self.settings.compiler), False)
        if minimum_version and Version(self.settings.compiler.version) < minimum_version:
            raise ConanInvalidConfiguration(f"{self.ref} requires C++{self._min_cppstd}, which your compiler does not support.")

    def build_requirements(self):
        self.tool_requires("cmake/[>=3.23 <4]")

    def layout(self):
        cmake_layout(self)

    def generate(self):
        tc = CMakeToolchain(self)
        tc.variables["LIELAB_INSTALL_LIBRARY"] = True
        tc.variables["LIELAB_BUILD_TESTS"] = False
        tc.variables["LIELAB_BUILD_PYTHON"] = True
        # tc.variables["PYTHON_EXECUTABLE"] = "path to python executable" # Define this to explicitly tell pybind11 which Python to use
        tc.generate()
        cmake = CMakeDeps(self)
        cmake.generate()
    
    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()
        cmake.install()
