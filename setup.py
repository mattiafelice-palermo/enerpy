from setuptools import setup, find_packages

setup(
    name='enerpy',
    version='0.1.0',  # Starting with an initial version
    packages=find_packages(),
    description='A Python tool for exploring bonded interactions in Molecular Mechanics force fields',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='Mattia Felice Palermo',
    author_email='mattiafelice.palermo@gmail.com',
    url='https://github.com/mattiafelice-palermo/enerpy',
    install_requires=[
        'numpy',  # Assuming numpy might be used, adjust as needed
        'scipy',  # Assuming scipy might be used for calculations
        'matplotlib',
    ],
    entry_points={
        'console_scripts': [
            'enerpy = enerpy.main:main',  # Point to the main function in enerpy.main
        ],
    },
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
    python_requires='>=3.7',
)
