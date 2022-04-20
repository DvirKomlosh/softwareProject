from setuptools import setup, find_packages, Extension

setup(
    name='spkmeansmodule',
    description="The C kmeans module",
    packages=find_packages(),
    ext_modules=[Extension('spkmeansmodule', sources=['spkmeansmodule.c'])]
)
