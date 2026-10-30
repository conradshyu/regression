# Polynomial Regression for Free Energy Estimates

Written by [Conrad Shyu](mailto:conradshyu@hotmail.com)<br>
*release 1, October 1, 2009*<br>
*release 2, March 11, 2014*<br>
*last updated on October 30, 2026*

## Synopsis

This project contains a C++ program that fits a polynomial to thermodynamic integration data and estimates the
free energy difference by integrating the fitted polynomial. It also reports a trapezoidal-rule estimate and
includes an example data file. The implementation uses the C++ Standard Library and requires a C++ compiler; no
additional libraries are required.

To compile the program, simply type:

`make`

The make program will invoke the Makefile and compiles the source code. Alternatively, to compile the files
manually, issue the command:

`g++ -I. driver.cpp regress.cpp -o regress -lm`

The project files are:

| Filename | Description |
| --- | --- |
| `driver.cpp` | the user interface program |
| `regress.cpp` | the main implementation of polynomial regression |
| `regress.h` | the header file for polynomial regression |
| `example.csv` | example file of simulations data from thermodynamic integration |
| `Makefile` | the makefile to compile the source code |
| `.gitignore` | ignores compiler outputs and intermediate files |
| `LICENSE` | GNU General Public License, version 3 or later |
| `README.md` | this file |

The program fits the data by linear least squares using powers of the lambda values. Its estimate function writes
values of the fitted polynomial for plotting.

The polynomial regression program takes four command-line parameters: `input_file`, `degree` (optional),
`plot_file` (optional), and `data_points` (optional). Only one parameter, the input file, is required to run the
software. The followings detail the format and meaning of each command-line parameter.

## `input_file`

The input file must be plain text, with one data row per line. The first column contains lambda values and the
second contains the corresponding ensemble averages (dG/dlambda) from the simulations. Separate the columns with
a comma, semicolon, or whitespace. Lines whose first character is `#` are ignored.
For example, the followings list the contents of the example input file:

```
0.0,  51.49866347
0.1,  23.92508775
0.2,  10.35390700
0.3,   2.58426990
0.4,  -2.18351656
0.5,  -5.41745387
0.6,  -7.62452181
0.7,  -9.25455804
0.8, -10.45592989
0.9, -11.39244138
1.0, -12.12433704
```

To run the analysis program with the required parameter, type:

`regress example.csv`

The program prints the fitted polynomial coefficients and estimates the free energy difference by integrating that
polynomial and by the trapezoidal rule. For example, the example file should generate the following output:

```
Degree, Coefficients
     0,      51.49866347
     1,    -440.06771988
     2,    2724.49944200
     3,  -16355.52325204
     4,   74074.45611540
     5, -226406.51705359
     6,  454849.19406749
     7, -590543.95494381
     8,  476225.01322752
     9, -216633.75369268
    10,   42443.03080908

Free energy difference
Regression: 0.65644834
Trapezoid: 1.02220063
```

## `degree`

This parameter specifies the degree of polynomial for the regression model. This parameter is optional. If omitted,
the program uses the smaller of the number of input rows minus one and 10. Explicit degrees above 10 are also
capped at 10. For example, to run the analysis with a polynomial of degree eight, type

`regress example.csv 8`

Alternatively, to use the default degree, type

`regress example.csv`

## `plot_file`

This parameter specifies the name of the file that contains estimates from the fitted polynomial. This parameter is
optional. If it is omitted, no plot data will be generated. The output covers lambda values from 0 to 1, regardless
of the range of lambda values in the input file. For example, to generate plot data, issue the command:

`regress example.csv 10 plot.csv`

The file, example.csv, is the input file that contains the thermodynamic integration data. The numerical value,
10, specifies the degree of polynomial. To generate a plot, it is necessary to specify the degree of polynomial
for the regression model. The file, plot.csv, is the output file that contains the estimates from the
interpolating polynomial.

## `data_points`

This parameter sets the number of equal intervals between lambda 0 and 1 in the plot output. It is optional; when
omitted, the number of intervals is one less than the number of input rows, so the output has the same number of
rows as the input. The output includes both endpoints, so a value of 100 produces 101 rows. For example, to
generate a file of plot data with 100 intervals, issue the command:

`regress example.csv 10 plot.csv 100`

## Reference

C. Shyu and F.M. Ytreberg (2009). Reducing the bias and uncertainty of free energy estimates by using regression
to fit thermodynamic integration data. *Journal of Computational Chemistry*, 30(14), 2297-2304.

## License
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any
later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
details.

You should have received a copy of the GNU General Public License along with this program. If not, see
<http://www.gnu.org/licenses/>.
