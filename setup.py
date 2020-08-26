import os, re, shutil, sys, sysconfig
import platform
import subprocess

from distutils.version import LooseVersion
from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext

# Project structure and CMake build steps adapted from
# https://www.benjack.io/2018/02/02/python-cpp-revisited.html

CMAKE_COMMAND = os.environ.get("CMAKE_COMMAND", "cmake")


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def run(self):
        try:
            out = subprocess.check_output([CMAKE_COMMAND, "--version"])
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: "
                + ", ".join(e.name for e in self.extensions)
            )

        if platform.system() == "Windows":
            cmake_version = LooseVersion(
                re.search(r"version\s*([\d.]+)", out.decode()).group(1)
            )
            if cmake_version < "3.1.0":
                raise RuntimeError("CMake >= 3.1.0 is required on Windows")

        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        cmake_args = [
            "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=" + extdir,
            "-DPYTHON_EXECUTABLE=" + sys.executable,
        ]

        cfg = "Debug" if self.debug else "Release"
        build_args = ["--config", cfg]

        if platform.system() == "Windows":
            cmake_args += [
                "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_{}={}".format(cfg.upper(), extdir)
            ]
            if sys.maxsize > 2 ** 32:
                cmake_args += ["-A", "x64"]
            build_args += ["--", "/m"]
        else:
            cmake_args += ["-DCMAKE_BUILD_TYPE=" + cfg]
            build_args += ["--", "-j2"]

        env = os.environ.copy()
        env["CXXFLAGS"] = '{} -DVERSION_INFO=\\"{}\\"'.format(
            env.get("CXXFLAGS", ""), self.distribution.get_version()
        )
        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)
        subprocess.check_call(
            [CMAKE_COMMAND, ext.sourcedir] + cmake_args, cwd=self.build_temp, env=env
        )
        subprocess.check_call(
            [CMAKE_COMMAND, "--build", "."] + build_args, cwd=self.build_temp
        )
        print()  # Add an empty line for cleaner output


def file_replace(original_path, match, replace):
    """Replace lines in files by applying regex

    Args:
        original_path (str): Path to file to update in place
        match (str): Matching regex
        replace (str): Replacement regex
    """
    # Adapted from SeqLib python package.
    with open(original_path) as original_file:
        original = original_file.readlines()
    with open(original_path, "w") as replaced_file:
        for line in original:
            replaced_line = re.sub(match, replace, line)
            replaced_file.write(replaced_line)


class SeqLibCMakeBuild(CMakeBuild):
    def run(self):
        # To link into a shared library we need to add the -fPIC flag to SeqLib dependencies
        # before building. Adapted from SeqLib python package.
        cflags_line_re = r"^(CFLAGS\s*=.*)$"
        bwa_makefile_path = os.path.join(
            os.path.dirname(os.path.realpath(__file__)), "lib/seqlib", "bwa", "Makefile"
        )
        file_replace(bwa_makefile_path, cflags_line_re, r"\1 -fPIC -Wno-unused-result")

        super().run()


with open("README.md") as f:
    readme = f.read()

with open("LICENSE") as f:
    license = f.read()

setup(
    name="npsv",
    version="0.1.0",
    description="Non-parametric genotyper for structural variants",
    long_description=readme,
    author="Michael Linderman",
    author_email="mlinderman@middlebury.edu",
    license=license,
    url="https://github.com/mlinderm/npsv",
    scripts=[
        "scripts/synthBAM",
        "scripts/vcf2bed",
        "scripts/vcf2samplot",
        "scripts/svviz22vcf",
        "scripts/padAlleles",
    ],
    entry_points="""
        [console_scripts]
        npsv=npsv.main:main
        npsvg=npsv.npsvg:main
    """,
    packages=find_packages("src"),
    package_dir={"": "src"},
    ext_modules=[CMakeExtension("npsv/npsv")],
    cmdclass=dict(build_ext=SeqLibCMakeBuild),
    zip_safe=False,
    test_suite="tests",
    include_package_data=True,
    data_files=[
        (
            "etc",
            [
                "etc/human_g1k_v37.genome",
                "etc/human_g1k_v37.gaps.bed.gz",
                "etc/human_g1k_v37.gaps.bed.gz.tbi",
            ],
        )
    ],
    classifiers=[
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python :: 3.6",
    ],
)
