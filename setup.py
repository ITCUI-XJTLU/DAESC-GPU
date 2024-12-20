from setuptools import setup, find_packages

# 读取长描述
try:
    with open('README.md', 'r', encoding='utf-8') as fh:
        long_description = fh.read()
except FileNotFoundError:
    long_description = "A Python package for allele specific gene expression"

setup(
    name='DAESC-GPU',  
    version='1.0', 
    author='Tengfei Cui',      
    author_email='tfcui23@uw.edu',  
    description='A Python package for allele specific gene expression',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/ITCUI-XJTLU/DAESC-GPU',  
    packages=find_packages(),  
    install_requires=[
        'numpy',
        'cupy',
        'pandas',
        'matplotlib',
        'scipy',
        'statsmodels',
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.10',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    python_requires='>=3.10', 
)
