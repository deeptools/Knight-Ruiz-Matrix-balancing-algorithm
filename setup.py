from setuptools.command.build_ext import build_ext
from distutils.core import setup, Extension
import platform
import sys
from sysconfig import get_config_var, get_paths
import logging

__version__ = '0.0.5dev'


def get_include():  # TODO
    info = get_paths()
    Eigen_path = '/'.join(info['include'].split('/')[:-1])
    Eigen_path += '/eigen3'
    return Eigen_path


def __extra_compile_args():
    extra_compile_args = []

    if platform.system() == 'Darwin':
        extra_compile_args = ["-std=c++11"]
    else:
        extra_compile_args = ["-fopenmp", "-std=c++11"]
    return extra_compile_args


def __extra_link_args():
    extra_link_args = []
    if platform.system() != 'Darwin':
        extra_link_args = ["-lgomp", "-lm", "-lrt"]
    return extra_link_args


sources_list = ['src/krbalancing.cpp']

kr_module = Extension('krbalancing',
                      sources=sources_list,
                      include_dirs=[
                          # Path to eigen3 headers
                          get_include()
                      ],
                      extra_link_args=__extra_link_args(),
                      extra_compile_args=__extra_compile_args()
                      )


setup(
    name='krbalancing',
    version=__version__,
    author='Leily Rabbani',
    author_email='leila.rabbani@gmail.com',
    description='A c++ extension for python to balance a matrix using KR method',
    ext_modules=[kr_module],
    install_requires=['pybind11>=2.2'],
#    headers = ['src/krbalancing.hpp']
)
