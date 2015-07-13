#!python

from distutils.core import setup, Extension

setup(name='fast',
      version='0.1.0',
      description='Fast scoring function',
      author='Ajasja Ljubetic',
      ext_modules=[Extension('_fastscoreCC', ['fastscoreCC.i', 'Interaction.cpp'],
                             swig_opts=['-c++'])],
    #  py_modules=['fastscoreCC.py']
     )