PyCORN
======

A script to extract data from UNICORN result (.res) files and plot them.

![With fractions - yay!](https://github.com/pyahmed/PyCORN/blob/dev/samples/sample1_2009Jun16no001_plot.jpg)

Description: 

A script to extract data from .res (results) files generated by UNICORN Chromatography software supplied with ÄKTA Systems. This script will find all data blocks, extract the data and write out csv-files. If you have matplotlib installed it will also plot all the curves including fractions if present. Plots can be saved in any format supported by matplotlib (default is pdf).

Limitations:
- See https://github.com/pyahmed/PyCORN/issues

Requirements:
- Python 2.7 or 3.x (Tested on Windows 7 / Mac OSX)
- optional: working matplotlib installation (for plotting)

Usage:
- See docs/USAGE.txt

License:
- GPLv2 see docs/LICENSE.txt
