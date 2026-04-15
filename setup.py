from setuptools import setup, find_packages

setup(
    name='xrd_profile',
    version='0.1.0',
    packages=find_packages(),
    install_requires=['numpy', 'scipy', 'matplotlib'],
    python_requires='>=3.8',
)
