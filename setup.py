import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
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
    classifiers=[
        "Development Status :: 1 - Planning",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Mathematics",
    ],
    python_requires='>=3.6',
)
