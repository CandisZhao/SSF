

# Protein Complexes Identification with FWER Control

## Introduction

SSF algorithm can be divided into two key modules: a seed expansion procedure for significant protein complexes
search and a set cover strategy for redundancy elimination. In the search stage of SSF, we adopt a local search
procedure to find protein complexes in an iterative manner.  In the redundancy elimination stage of SSF, we employ 
a set cover formulation to remove some redundant protein complexes.

## Usage

Implementation of SSF is available in two different languages;
in C++ and in Python 2.7.

* Import the C++ code into Miscrosoft Visual Studio, and then load
some datasets to test in the file conf.h. By the way, the
significance level alpha also could be user-specific. After
that complie and run the code. Finally, the result will be obtained.

* As to SSF in python, users could perform the following command in the command line:

  * python upper_bound_test.py input_file.txt output_file.txt alpha

  And then, the result will be reported in the output_file.txt, and the corresponding
assessment measure will be layout on the terminal.

In addition, we also provide the five datasets and three reference datasets.
We adopt the "oNMI" to calculate the NMI between the predicted protein complex
and the true protein complex in the reference set. And the file of match_standalone.py
is used to assess the algorithm's performance.
