from setuptools import setup, find_packages

setup(
    name='ncbi-cluster-tracker',
    version='0.1.0',
    packages=find_packages(),
    py_modules=['main'],
    entry_points={
        'console_scripts': [
            'ncbi-cluster-tracker = main:main'
        ]
    },
    include_package_data=True,
)
