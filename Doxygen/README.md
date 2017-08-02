# Overview

This directory contain the source files for the PICSAR documentation. The documentation can be built (in HTML format) using Doxygen.

# Building the documentation

In order to build the documentation, you need to first install Doxygen and pandoc. 

Doxygen can be installed under Linux Debian by using
```
apt-get install doxygen
```
and under MacOSX by using MacPorts
```
sudo port install doxygen
```

Pandoc can be installed with anaconda
```
conda install pandoc
```

Then the documentation can be built by typing
```
make html
```


# Visualizing the documentation

Once built, the documentation can be visualized by open the file `html/index.html` with a web browser.




