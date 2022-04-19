from setuptools import setup, find_packages, Extension

setup(
    name='spkmeans',
    description="The C kmeans module",
    packages=find_packages(),
    ext_modules=[Extension('spkmeans', sources=['spkmeansmodule.c'])]
)
