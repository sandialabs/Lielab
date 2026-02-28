from conan import ConanFile
from conan.tools.cmake import CMakeToolchain, CMakeDeps, CMake, cmake_layout
from conan.tools.build import check_min_cppstd
from conan.tools.files import copy, get, rmdir
from conan.tools.scm import Version
from conan.errors import ConanInvalidConfiguration
import os

required_conan_version = ">=2.0"

class LielabConan(ConanFile):
    name = "lielab"
    homepage = "https://github.com/sandialabs/Lielab"
    settings = "os", "compiler", "build_type", "arch"
    options = {"fPIC": [True, False],
               "with_assertions": [True, False],
               "with_tests": [True, False],
               "with_coverage": [True, False],
               "with_python": [True, False]}
    
    default_options = {"fPIC": True,
                       "with_assertions": True,
                       "with_tests": False,
                       "with_coverage": False,
                       "with_python": False}

    exports_sources = ("Lielab/*",
                       "CMakeLists.txt",
                       "Lielab.hpp",
                       "LICENSE",
                       "SCR")

    def requirements(self):
        self.requires("eigen/[>=5.0.0 <6]", transitive_headers=True)
        self.requires("fmt/[>=12.0.0 <13]")
        if self.options.get_safe("with_tests"):
            self.requires("catch2/3.11.0")
        if self.options.get_safe("with_python"):
            self.requires("pybind11/[>=3.0.0]")

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

    def config_options(self):
        if self.settings.os == "Windows":
            del self.options.fPIC

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

        if "with_tests" in self.options:
            tc.variables["LIELAB_BUILD_TESTS"] = self.options.with_tests
        if "with_assertions" in self.options:
            tc.variables["LIELAB_WITH_ASSERTIONS"] = self.options.with_assertions
        if "with_coverage" in self.options:
            tc.variables["LIELAB_WITH_COVERAGE"] = self.options.with_coverage
        if "with_python" in self.options:
            tc.variables["LIELAB_BUILD_PYTHON"] = self.options.with_python
        # tc.variables["PYTHON_EXECUTABLE"] = "path to python executable" # Define this to explicitly tell pybind11 which Python to use
        tc.generate()
        cmake = CMakeDeps(self)
        cmake.generate()
    
    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()
    
    def package(self):
        copy(self, "LICENSE", dst=os.path.join(self.package_folder, "licenses"), src=self.source_folder)
        copy(self, "SCR", dst=os.path.join(self.package_folder, "licenses"), src=self.source_folder)
        cmake = CMake(self)
        cmake.configure()
        cmake.install()
        rmdir(self, os.path.join(self.package_folder, "share"))
    
    def package_info(self):
        # self.cpp_info.bindirs = []
        self.cpp_info.libs = ["Lielab"]

        self.cpp_info.set_property("cmake_file_name", "Lielab")
        self.cpp_info.set_property("cmake_target_name", "Lielab::Lielab")

        self.cpp_info.names["cmake_find_package"] = "Lielab"
        self.cpp_info.names["cmake_find_package_multi"] = "Lielab"
