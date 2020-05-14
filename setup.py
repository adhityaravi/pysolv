#! usr/bin/env python 3
# _*_ coding: utf-8 _*_

import setuptools
from numpy.distutils.core import Extension
from numpy.distutils.core import setup

# Extension libraries
f_sor = Extension(
                    name='pysolv.f_lib.f_sor',
                    sources=['pysolv/st_iter_solv/f_sor.f90',
                             'pysolv/tools/f_tools.f90'],
                    extra_f90_compile_args=['-fcheck=all', '-O3'],
                    library_dirs=['/usr/lib'],
                    libraries=['lapack', 'blas']
                 )

f_jacobi = Extension(
                        name='pysolv.f_lib.f_jacobi',
                        sources=['pysolv/st_iter_solv/f_jacobi.f90',
                                 'pysolv/tools/f_tools.f90'],
                        extra_f90_compile_args=['-fcheck=all', '-O3'],
                        library_dirs=['/usr/lib'],
                        libraries=['lapack', 'blas']
                    )

f_cg = Extension(
                    name='pysolv.f_lib.f_cg',
                    sources=['pysolv/krylov_solv/f_cg.f90',
                             'pysolv/tools/f_tools.f90'],
                    extra_f90_compile_args=['-fcheck=all', '-O3'],
                    library_dirs=['/usr/lib'],
                    libraries=['lapack', 'blas']
                )

f_sd = Extension(
                    name='pysolv.f_lib.f_sd',
                    sources=['pysolv/oneD_proj_solv/f_sd.f90',
                             'pysolv/tools/f_tools.f90'],
                    extra_f90_compile_args=['-fcheck=all', '-O3'],
                    library_dirs=['/usr/lib'],
                    libraries=['lapack', 'blas']
                )


# Setup
with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
        name="pysolv",
        version="0.0.1",
        author="Adhitya-Sriram Ravi",
        author_email="adhityashan95@gmail.com",
        description="A collection of linear solvers including stationary iterative solvers like Gauss-Siedel (GS), "
                    "Jacobi (JM), Successive over-relaxation (SOR) and Krylov subspace methods like Conjugate gradient "
                    "(CG), generalized minimal residual (GMRES) etc.",
        long_description=long_description,
        long_description_content_type="text/markdown",
        url="https://framagit.org/adhityaravi/pysolv",
        packages=setuptools.find_packages(),
        ext_modules=[f_sor, f_jacobi, f_cg, f_sd, ],
        classifiers=[
            "Development Status :: 1 - Planning",
            "License :: OSI Approved :: MIT License",
            "Natural Language :: English",
            "Operating System :: OS Independent",
            "Programming Language :: Python :: 3",
            "Programming Language :: Fortran",
            "Topic :: Scientific/Engineering :: Mathematics",
        ],
        python_requires='>=3.6',
        install_requires=['numpy'],
     )
