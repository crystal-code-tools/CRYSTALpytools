import setuptools

long_description = 'This repository contains functions to be used with the\
 <a href="https://www.crystal.unito.it/index.php">CRYSTAL code</a>.'

setuptools.setup(
    name="crystal_functions",
    version="2022.6.16",
    author_email="crystalcodetools@gmail.com",
    description="Functions to be used with the CRYSTAL code.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/crystal-code-tools/crystal_functions",
    project_urls={
        "Bug Tracker": "https://github.com/crystal-code-tools/crystal_functions/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=setuptools.find_packages(include=['crystal_functions', 'crystal_functions.*']),
    #python_requires=">=3.8",
    install_requires=[
	"pymatgen",
	"mendeleev",
	"ase"
    ]
)
