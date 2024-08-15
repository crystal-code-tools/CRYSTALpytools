# CRYSTALpytools

This repository contains functions that allow the user to access the [CRYSTAL code](https://www.crystal.unito.it/index.php), input, output and execution from a python infrastructure, such as Jupyter Notebooks. Although they can achieve this goal on their own, they achieve their full potential when used together with [pymatgen](https://pymatgen.org/index.html). In this latter scenario, the CRYSTALpytools may be used as a layer between CRYSTAL and pymatgen.

In January 2022 the first stable version (v2022.1.10) was released.

## Documentation

The **[documentation site](https://crystal-code-tools.github.io/CRYSTALpytools/)** of CRYSTALpytoosl explains the usage of submodules and classes. The latest released version of CRYSTALpytools is available via [PyPI](https://pypi.org/project/CRYSTALpytools/#history).

## Structure

The CRYSTALpytools module aims at providing the user a python interface to the CRYSTAL code. The central data structure, called 'Crystal_object' is created by the `crystal_io` by parsing CRYSTAL input/output files. The flowchart below is aimed at showing how different parts of the module interact with the Crystal_objects.

![crystal_object](docs_source/_static/crystal_object.png)

## Installation

### Create a conda/anaconda environment

This step is not mandatory, but it makes using CRYSTALpytools very smooth. It is, therefore, very recommended. If you are new to anaconda, please follow [these steps](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) to install it on your computer.

Create a new conda environment:

``` console
conda create --name crystal python=3.9
```

In the line above, “crystal” is the name of the environment and can be set to any you like. The “python=3.9” ensures that the suitable python distribution is installed.

Activate the conda environment:

``` console
conda activate crystal
```

### Install CRYSTALpytools

The CRYSTALpytools package can be installed from pip. Pip is a package-management system written in Python and is used to install and manage software packages (called modules in python).

``` console
pip install --upgrade CRYSTALpytools
```

Windows users might need to install windows-curses. This can be done by using:

``` console
pip install windows-curses
```

To check that CRYSTALpytools was install please type

``` console
conda list
```

This will return a list of all the modules installed in the environment. Here there should be CRYSTALpytools. If this was not the case, something went wrong during the installation. Please check the location of the environment that is being displayed. This appears at the beginning of the `conda list` command. The most common mistake at this stage is that the environment was not activated as described above.


Please note that pip will only install the functions and not the example notebooks. This decision was taken in order to reduce the volume of data transferred when installing. If you are interested in the example notebooks please read the section below.

### Set the path to runcry and runprop

If you intend to run CRYSTAL on the machine where you are running the CRYSTALpytools, the path to your local runcry amd runprop needs to be specified. To do so, please run the `set_runcry_path` and `set_runprop_path` functions:

``` console python
>>> from CRYSTALpytools.execute import set_runcry_path, set_runprop_path
>>> set_runcry_path('path_to_your_runcry')
>>> set_runprop_path('path_to_your_runcry')
```

### Testing

To test the CRYSTALpytools please run the example notebooks in the [example folder](examples/).

## Examples and Test Cases

Notebooks in the [example folder](examples/). are the testing examples used during the development of CRYSTALpytools. They are designed to be self-explanatory and can be used as (basic) tutorials.

## Tutorials

Tutorials designed as workflows can be found in the [tutorial folder](tutorial/). *IN developing*

