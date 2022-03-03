from setuptools import setup, find_packages

setup(
    name='hubmap_notebook_helper_functions',
    version='0.1.0',
    description='Utilities for using HuBMAP APIs and visualizing their data in Jupyter notebooks',
    url='https://github.com/hubmapconsortium/cross-dataset-common',
    author='Sean Donahue',
    author_email='seandona@andrew.cmu.edu',
    license='MIT',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Programming Language :: Python :: 3 :: Only',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    packages=find_packages(),
    install_requires=[
        'hubmap_api_py_client',
        'matplotlib',
        'pandas',
        'requests>=2.22.0',
        'seaborn',
    ],
    python_requires='>=3.6',
)
