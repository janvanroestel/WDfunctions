import setuptools

setuptools.setup(
    name="WDfunctions",
    version="0.1.0",
    url="https://github.com/janvanroestel/WDfunctions",
    author="Jan van Roestel",
    author_email="jcjvanroestel@gmail.com",
    description="A collection of functions related to white dwarf stars",
    long_description=open('README.md').read(),
    packages=setuptools.find_packages(),
    install_requires=['setuptools-git'],
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
    ],
    include_package_data=True,
    package_data={'': ['data/WDtracks/*','WDfunctions/data/WDtracks/*']},
)
