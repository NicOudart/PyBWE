from setuptools import setup,find_packages

setup(name='PyBWE',
version='1.0',
description='Python tools for Bandwidth Extrapolation (BWE) of radar data',
long_description=open("requirements.txt").read(),
long_description_content_type='text/markdown',
url='#',
author='Nicolas OUDART',
author_email='nicolas.oudart@latmos.ipsl.fr',
packages=find_packages('src'),
install_requires=open("requirements.txt").read(),
package_dir={ '' : 'src' },
python_requires='>=3.7',
zip_safe=False,
classifiers=[
    'Programming Language :: Python :: 3',
    'License :: OSI Approved :: MIT License',
    'Operating System :: OS Independent',
    'Intended Audience :: Science/Research'
    'Topic :: Scientific/Engineering :: Astronomy'])