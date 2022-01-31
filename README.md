# crystal_functions
This repository contains functions to be used with the
<a href="https://www.crystal.unito.it/index.php">CRYSTAL code</a>.

This is a work in progress, if you would like to contribute please get in touch with camino.bruno@gmail.com.

## Installation
A list of the dependecies used by crystal_functions can be found in this file.

### conda
To install crystal_functions using conda, please use:
```console
conda install -c crystal-python-tools crystal_functions
```

### pip
To install crystal_functions using pip, please use:
```console
pip install crystal_functions
```
Please note that both conda and pip will only install the functions and not the example notebooks. This decision was taken in order to reduce the volume of data transferred when installing. If you are interested in the example notebooks please read the section below.

## Examples
Each function is documented in Jupyter Notebooks that can be found in the  [example folder](example/). There is one notebook per function file (e.g. the functions contained in file_read_write.py are explained in the example/file_read_write.ipynb notebook).

## Usage
The functions are divided into files depending on their ultimate goal. For example, all the i/o functions are saved in crystal_functions/file_read_write.py. To access them, please use:

```console
from crystal_functions.file_read_write import Crystal_output

Crystal_output('output_name.out')
```
Each individual function contains either 'crystal' or 'cry' in its name. This was chosen, despite making the names of the functions longer, in order to avoid ambiguity. This means that when calling a function, you will know that it refers to a crystal_functions function and not, for example, a pymatgen one with a similar name.
