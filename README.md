# xtpurt
_Stata implementation of modern panel unit root tests for heteroskedastic panels_

The Stata command `xtpurt` implements the heteroskedasticity-robust panel unit root tests (PURTs) suggested in

- Herwartz, H. and F. Siedenburg. 2008. Homogenous panel unit root tests under cross sectional dependence: Finite sample modifications and the wild bootstrap. _Computational Statistics and Data Analysis_ 53(1): 137-150.

- Demetrescu, M. and C. Hanck. 2012. A simple nonstationary-volatility robust panel unit root test. _Economics Letters_ 117(2): 10-13.

- Herwartz, H., S. Maxand and Y. M. Walle. 2017. Heteroskedasticity-robust unit root testing for trending panels. _Center for European, Governance and Economic Development Research discussion papers 314_.

While the former two tests are robust to time varying volatility when the data contain an intercept only, the latter test is unique in the sense that it is asymptotically pivotal for trending heteroskedastic panels. Moreover, `xtpurt` incorporates lag order selection, prewhitening and detrending procedures to account for serial correlation and trending data.

_The vignette paper for `xtpurt` is under review at the Stata Journal. Hence, we cannot publish it here. As soon as possible, we will provide references to the Stata Journal article._

## Installation

Simply download the contents of the program folder _/ado_ to your hardisk. The folder contains two files:

- xtpurt.ado : the Stata command `xtpurt`.
- xtpurt.sthlp : the corresponding Stata help file, e.g. `help xtpurt`.

Basically, you have three options how to employ additional ado files:

1. Same directory: Copy both files into your project folder, i.e. typically where your do-file is located.

2. Specific directory: You might have or want to create a specific command library folder where you keep your collection of Stata extensions. Copy both files in this folder, e.g. _C:\\mystata\\ado\\_. In Stata you import all functions by means of `adopath + "C:\mystata\ado"`.

3. Personal directory: You can install `xtpurt` system wide by copying both files to the  [personal ado directory](http://www.stata.com/support/faqs/programming/personal-ado-directory/) at _C:\\ado\\personal\\_.

## Examples

We provide an example application on price levels and inflation in _/examples_.

1. Download the folder _/examples_.

2. Start Stata and change the working directory to your _/examples_ copy.

3. Run _prices.do_ and observe the results.

4. Consider the help file, `help xtpurt`, for further options.
