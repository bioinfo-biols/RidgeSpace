from setuptools import setup

setup(
    name='RidgeSpace',
    version='0.1.0',
    author='Ruiqiao He',
    author_email='ruiqiaohe@gmail.com',
    packages=['RidgeSpace'],
    license="GPL",
    url='http://pypi.python.org/pypi/RidgeSpace/',
    description='Unravelling three-dimensionally dynamics of spatial multi-modal data',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    install_requires=[
        "matplotlib>=3.5.0",
        "numpy",
        "pandas",
        "scipy",
    ],
    python_requires='>=3.7.1',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        "Programming Language :: Python :: 3",
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        "Programming Language :: Python :: 3.9",
        'Topic :: Scientific/Engineering :: Bio-Informatics',
      ],
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'RidgeSpace=RidgeSpace.main:RidgeSpace_command',
        ]
    }
)
