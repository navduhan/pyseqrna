import setuptools 
setuptools.setup(
    name='pySeqRNA',
    version='0.1',
    description ="pySeqRNA a python based RNA analysis package",
    url="#",
    author="Naveen Duhan",
    author_email = "naveen.duhan@usu.edu",
    license='MIT',
    packages = setuptools.find_packages(),
    include_package_data=True,
    entry_points={
            'console_scripts': [
                    'pyseqrna = pyseqrna.__main__:main',
            ]
    },

 classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix",
        
        # Indicate who your project is intended for
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',


        # Specify the Python versions you support here. In particular, ensure
        # that you indicate whether you support Python 2, Python 3 or both.
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7'
    ],
    package_data={'': ['param/*.ini']},
)