# VIPER

**VIPER** is an open-source software tool for the analysis of captured velocity time series. VIPER stand for **V**elocity time-series **i**nspection, **p**ost-processing and **er**ror filtering tool and contains a collection of most common, best-practice methods used for velocity time series evaluation.


## System requirements

When using from MATLAB, VIPER requires pwelch.m for computing Welch's spectrum:
- either pwelch.m provided by MathWorks (e.g. [Signal Processing Toolbox](https://www.mathworks.com/products/signal.html)) 
- or a 3rd-patry pwelch.m function (e.g. pwelch.m written by Peter V. Lanspeary)

When using as executable, **VIPER** requires installed [MATLAB Runtime](https://www.mathworks.com/products/compiler/matlab-runtime.html) of version specified in the release.


## Applied methods

The algorithms implemented in **VIPER** are described in the this book, which we kindly ask you to cite when publishing results obtained by **VIPER**:

Sokoray-Varga, Béla (2022): Inspection, post-processing and error filtering of captured velocity time series. Karlsruhe: Bundesanstalt für Wasserbau.
Available at: https://hdl.handle.net/20.500.11970/108631


## License 

**VIPER** is distributed by the [Federal Waterways Engineering and Research Institute](http://www.baw.de/) 
and is freely available and open source, licensed under the [GNU General Public License 3](https://www.gnu.org/licenses/gpl.html). 
See [LICENSE.txt](LICENSE.txt) for details.

**VIPER** uses the following third-party function, which is distributed under its own terms. See [3RD-PARTY.txt](3RD-PARTY.txt) for details.

- [pwelch.m](pwelch.m) by Peter V. Lanspeary: GNU General Public License 2
