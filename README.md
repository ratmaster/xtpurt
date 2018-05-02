# xtpurt

Current version: 1.1.0 (24 August 2017)

[![Package-License](https://img.shields.io/github/license/mashape/apistatus.svg)](https://opensource.org/licenses/MIT)

_Stata implementation of modern panel unit root tests for heteroskedastic panels_

The Stata command `xtpurt` implements the heteroskedasticity-robust panel unit root tests (PURTs) suggested in

- Herwartz, H. and F. Siedenburg. 2008. Homogenous panel unit root tests under cross sectional dependence: Finite sample modifications and the wild bootstrap. _Computational Statistics and Data Analysis_ 53(1): 137-150.

- Demetrescu, M. and C. Hanck. 2012. A simple nonstationary-volatility robust panel unit root test. _Economics Letters_ 117(2): 10-13.

- Herwartz, H., S. Maxand and Y. M. Walle. 2017. Heteroskedasticity-robust unit root testing for trending panels. _Center for European, Governance and Economic Development Research discussion papers 314_.

While the former two tests are robust to time-varying volatility when the data contain an intercept only, the latter one is unique in the sense that it is asymptotically pivotal for trending heteroskedastic panels. Moreover, `xtpurt` incorporates lag order selection, prewhitening and detrending procedures to account for serial correlation and trending data.

## Release

The software package and vignette paper for `xtpurt` were published in the Stata Journal ([direct link](https://www.stata-journal.com/article.html?article=st0519)), and the package is now available directly within Stata and can be downloaded and installed automatically. When employing the implemented tests of `xtpurt` in research, we would appreciate a reference to our software as suggested below:

 > Herwartz, H., Maxand, S., Raters, F. H. C., and Walle, Y. M. (2018). Panel unit root tests for heteroskedastic panels. The Stata Journal, 18(1), 184-196.

## Installation

Just download the contents of the program folder _/ado_ to your harddisk. The folder contains two files:

- xtpurt.ado : the Stata command `xtpurt`.
- xtpurt.sthlp : the corresponding Stata help file, e.g., `help xtpurt`.

Basically, you have three options how to employ additional ado files:

1. The same directory: Copy both files into your project folder, i.e., typically at the location of your do-file.

2. A specific directory: You might have or want to create a specific command library folder where you keep your collection of Stata extensions. Copy both files in this folder, e.g., _C:\\mystata\\ado\\_. In Stata, you import all functions through `adopath + "C:\mystata\ado"`.

3. The personal directory: You can install `xtpurt` system-wide by copying both files to the  [personal ado directory](http://www.stata.com/support/faqs/programming/personal-ado-directory/) at _C:\\ado\\personal\\_.

### Update

Re-download the files from the folder _/ado_ to your harddisk and replace the existing.

Afterward, restart Stata or run the command `discard` for refreshing your ado libraries.

## Examples

We provide an example application on price levels and inflation in _/examples_.

1. Download the folder _/examples_.
2. Start Stata and change the working directory to your _/examples_ copy.
3. Run _prices.do_ and observe the results.
4. Consider the help file, `help xtpurt`, for further options.
