from setuptools import setup, find_packages

setup(
    name='hamburger',
    version='0.2.1',
    author='David Williams',
    author_email='dwilliams001@dundee.ac.uk',
    description='A tool to extract and analyse contiguous sets of genes in bacterial genomes using HMMs.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/djw533/hamburger',
    packages=find_packages(
        #exclude these as packages - but import the data within these directories below
        exclude=["hamburger.models.T6SS",
                 "hamburger.models.full",
                 "hamburger.models.individual",
                 "hamburger.r_scripts",
                 "hamburger.t6ss_reference_set"]),    
    include_package_data=True,
    package_data={
        # include everything under hamburger/models
        'hamburger': ['models/*/*.hmm','r_scripts/*.R','t6ss_reference_set/*']
    },
    install_requires=[
        'biopython>1.80',
        'tqdm',
        'pandas'
    ],
    entry_points={
        'console_scripts': [
            'hamburger=hamburger.cli:main',  # Assumes hamburger.py has a main() function
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.7',
)
