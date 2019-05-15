from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize



extensions = [
    Extension(
        'cython_collision',
        ['cython_collision.pyx']
        ),
    Extension(
        'cython_force_calc',
        ['cython_force_calc.pyx']
        ),
    ]


setup(name='solar_system',
         version='0.1.0',
         author='Zach Bateman',
         author_email='zkbateman@gmail.com',
         description='A solar system formation simulation modeling many particles under the force of gravity.',
         packages=find_packages(),
         ext_modules = cythonize(extensions)
         )
