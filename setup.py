import setuptools
from datetime import date

long_description = 'Python tools for the <a href="https://www.crystal.unito.it/index.php">CRYSTAL code</a> developed and mantained by the CRYSTAL code developers'

setuptools.setup(
    name="CRYSTALpytools",        
    version=date.today().strftime('%Y.%m.%d'),
    author_email="crystalcodetools@gmail.com",
    description="Python tools for the CRYSTAL code developed and mantained by the CRYSTAL code developers.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/crystal-code-tools/CRYSTALpytools",
    project_urls={
        "Bug Tracker": "https://github.com/crystal-code-tools/CRYSTALpytools/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=setuptools.find_packages(include=['CRYSTALpytools', 'CRYSTALpytools.*']),
    python_requires=">=3.9",
    install_requires=[
	"numpy<2.0",
	"sympy",
	"scipy",
	"matplotlib",
	"pandas",
	"PyYAML",
	"mendeleev>=0.14.0",
    "pymatgen>=2022.7.25",
    "ase>=3.22.1",
    "basis_set_exchange>=0.9.1"
    ]
)
