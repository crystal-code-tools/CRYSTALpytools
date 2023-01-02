import setuptools

long_description = 'Python tools for the <a href="https://www.crystal.unito.it/index.php">CRYSTAL code</a> developed and mantained by the CRYSTAL code developers'

setuptools.setup(
    name="CRYSTALpytools",        
    version="2023.01.03",
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
    #python_requires=">=3.8",
    install_requires=[
	"pymatgen",
	"mendeleev",
	"ase"
    ]
)
