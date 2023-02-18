from setuptools import setup, find_packages


setup(
    name='SingleCellGGM',
    version='1',
    license='BSD-3 Cluse',
    author="MaShisongLab",
    packages=find_packages('src'),
    package_dir={'': 'src'},
    url='https://github.com/MaShisongLab/SingleCellGGM',
    keywords='SingleCellGGM',
    install_requires=[
          'NumPy','pandas'
      ],

)
