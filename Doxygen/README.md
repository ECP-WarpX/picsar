# Overview

This directory contain the source files for the PICSAR documentation. The documentation can be built (in HTML format) using Doxygen.

# Building the documentation

In order to build the documentation, you need to first install. Doxygen can be installed under Linux Debian by using
```
apt-get install doxygen
```
and under MacOSX by using MacPorts
```
sudo port install doxygen
```

The documentation can then be built by typing, from the current directory:
```
doxygen Doxyfile
```

# Visualizing the documentation

Once built, the documentation can be visualized by open the file `html/index.html` with a web browser.




