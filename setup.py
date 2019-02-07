from setuptools.command.build_ext import build_ext
from distutils.core import setup, Extension
import sys
from sysconfig import get_config_var, get_paths
import logging

__version__ = '0.0.1dev'


def get_include():  # TODO
    info = get_paths()
    Eigen_path = '/'.join(info['include'].split('/')[:-1])
    Eigen_path += '/eigen3'
    return Eigen_path


sources_list = ['src/KRBalancing.cpp']

kr_module = Extension('KRBalancing',
                      sources=sources_list,
                      include_dirs=[
                          # Path to eigen3 headers
                          get_include()
                      ],
                      extra_link_args=["-lgomp", "-lm", "-lrt"],
                      extra_compile_args=["-fopenmp", "-std=c++11"]
                      )


setup(
    name='KRBalancing',
    version=__version__,
    author='Leily Rabbani',
    author_email='leila.rabbani@gmail.com',
    description='A c++ extension for python to balance a matrix using KR method',
    ext_modules=[kr_module],
    install_requires=['pybind11>=2.2']
)
