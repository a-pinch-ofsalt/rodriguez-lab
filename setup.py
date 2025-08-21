from setuptools import setup, find_packages

setup(
    name='franklin-lab-notebook',
    version='0.1.0',
    description="Data & Results from Franklin's lab notebook.",
    author='Franklin Hiciano',
    author_email='fhiciano5@gmail.com',
    url='https://github.com/a-pinch-ofsalt/rodriguez-lab',
    packages=find_packages(),
    install_requires=[
        'requests',
        'pandas',
        'biopython',
        'receptor-utils'
    ],
    python_requires='>=3.8',
)
