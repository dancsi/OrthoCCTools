#!python

from distutils.core import setup, Extension

setup(name='fastscoreCC',
      version='0.1.0',
      description='Fast scoring function',
      author='Ajasja Ljubetic',
      ext_modules=[Extension('_fastscoreCC', ['Interaction.cpp', 'fastscoreCC.i'],
                             swig_opts=['-c++'], extra_compile_args=['-std=c++11', '-D_hypot=hypot'])
                  ],
     py_modules=['fastscoreCC']
     #package_data={'': ['fastscoreCC.py']}
     )