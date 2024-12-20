from setuptools import setup, find_packages

setup(
    name='DAESC-GPU',  
    version='1.0', 
    author='Tengfei Cui',      
    author_email='tfcui23@uw.edu',  
    description='A Python package for allele specific gene expresion',
    long_description=open('README.md').read(), 
    long_description_content_type='text/markdown',
    url='https://github.com/ITCUI-XJTLU/DAESC-GPU',  
    packages=find_packages(),  
    install_requires=[
        'numpy', 
        'cupy',
        'time',
        'math',
        'pandas',
        'csv', 
        'matplotlib',
        'scipy',
        'gc',
        'inspect',
        'numpy', 
        'pandas', 
        'statsmodels' 
        'scipy' 
        'time'
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.10', 
)
