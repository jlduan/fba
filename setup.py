from setuptools import setup, find_packages

with open('README.md', 'r') as fh:
    long_description = fh.read()

setup(
    name='fba',
    version='0.0.8',
    author='JD',
    description='Tools for feature barcoding analyses',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/jlduan/fba',
    packages=find_packages(exclude=('docs', 'tests')),
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX :: Linux',
        'Operating System :: Unix',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    python_requires='>=3.6',
    install_requires=[
        'regex',
        'polyleven>=0.5',
        'dnaio',
        'pysam>=0.14.0',
        'pandas',
        'numpy',
        'umi_tools>=1.0.0',
        'scipy',
        'scikit-learn',
        'statsmodels',
        'pyclustering',
        'hdbscan',
        'matplotlib>=3.3',
        'umap-learn'
    ],
    entry_points={
        'console_scripts': ['fba=fba.__main__:main'],
    }
)
