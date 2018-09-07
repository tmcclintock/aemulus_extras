from setuptools import setup,find_packages

setup(
    name='aemulus_extras',
    version='0.1',
    author='Tom McClintock',
    description='Analytic quantities for the Aemulus simulations.',
    py_modules=['aemulus_extras'],
    include_package_data=True,
    # adding packages
    packages=find_packages(),
    install_requires=['numpy',],
    setup_requires=['pytest_runner'],
    tests_require=['pytest']
)
