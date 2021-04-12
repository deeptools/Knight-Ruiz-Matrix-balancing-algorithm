import os
import re
import sys
import sysconfig
import platform
import subprocess
from pathlib import Path

__version__ = "0.5.0.a2"

from setuptools import setup, Extension, find_packages
from setuptools.command.build_ext import build_ext


class CMakeExtension(Extension):
    def __init__(self, name):
        Extension.__init__(self, name, sources=[])


class CMakeBuild(build_ext):

    user_options = build_ext.user_options + [
        ("cmake-extra-args=", None, "Extra cmake arguments")
    ]

    def initialize_options(self):
        super().initialize_options()
        self.cmake_extra_args = ''

    def finalize_options(self):
        super().finalize_options()

    def run(self):
        try:
            out = subprocess.check_output(["cmake", "--version"])
        except OSError:
            raise RuntimeError(
                "CMake must be installed to build the following extensions: "
                + ", ".join(e.name for e in self.extensions)
            )

        print(f"[INFO] Extra cmake args: {self.cmake_extra_args}")
        build_directory = os.path.abspath(self.build_temp)

        cmake_args = [
            "-DCMAKE_LIBRARY_OUTPUT_DIRECTORY=" + build_directory,
            "-DPython3_EXECUTABLE=" + sys.executable,
        ] + self.cmake_extra_args.split()

        cfg = "Debug" if self.debug else "Release"
        build_args = ["--config", cfg]

        cmake_args += ["-DCMAKE_BUILD_TYPE=" + cfg]

        # Assuming Makefiles
        build_args += ["--", "-j2"]

        self.build_args = build_args

        env = os.environ.copy()
        env["CXXFLAGS"] = '{} -DVERSION_INFO=\\"{}\\"'.format(
            env.get("CXXFLAGS", ""), self.distribution.get_version()
        )

        if not os.path.exists(self.build_temp):
            os.makedirs(self.build_temp)

        # CMakeLists.txt is in the same directory as this setup.py file
        cmake_list_dir = os.path.abspath(os.path.dirname(__file__))
        print(f"\n\n{'='*10} Running CMake prepare {'='*40}")
        print(f"[INFO] Running {' '.join(cmake_args)} in {self.build_temp}.")
        subprocess.check_call(
            ["cmake", cmake_list_dir] + cmake_args, cwd=self.build_temp, env=env
        )

        print(f"\n\n {'='*10} Building extensions {'='*40}")
        cmake_cmd = ["cmake", "--build", "."] + self.build_args
        subprocess.check_call(cmake_cmd, cwd=self.build_temp)

        # Move from build temp to final position
        for ext in self.extensions:
            self.move_output(ext)

    def move_output(self, ext):
        build_temp = Path(self.build_temp).resolve()
        dest_path = Path(self.get_ext_fullpath(ext.name)).resolve()
        source_path = build_temp / self.get_ext_filename(ext.name)
        dest_directory = dest_path.parents[0]
        dest_directory.mkdir(parents=True, exist_ok=True)
        print(f"\n\n[INFO] moving {source_path} to {dest_path}")
        self.copy_file(source_path, dest_path)


ext_modules = [
    CMakeExtension("krbalancing"),
]

setup(
    name="krbalancing",
    description="A c++ extension for python to balance a matrix using KR method",
    version=__version__,
    author="Leily Rabbani",
    author_email="leila.rabbani@gmail.com",
    maintainer="Dilawar Singh",
    maintainer_email="dilawar.s.rajput@gmail.com",
    packages=["krbalancing"],
    package_dir={"krbalancing": "."},
    ext_modules=ext_modules,
    cmdclass=dict(build_ext=CMakeBuild),
    zip_safe=False,
)

## setuptools.setup(
##     name="krbalancing",
##     description="A c++ extension for python to balance a matrix using KR method",
##     version=__version__,
##     author="Leily Rabbani",
##     author_email="leila.rabbani@gmail.com",
##     maintainer="Dilawar Singh",
##     maintainer_email="dilawar.s.rajput@gmail.com",
##     packages=["krbalancing"],
##     package_dir={"krbalancing": "."},
##     package_data={"krbalancing": ["krbalancing*.*"]},
##     has_ext_modules=lambda: True,
## )
