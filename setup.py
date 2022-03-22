from setuptools import setup, find_packages, Extension

setup(
    name='mykmeanssp',
    description="The C kmeans module",
    packages=find_packages(),
    ext_modules=[Extension('mykmeanssp', sources=['kmeans.c'])]
)
