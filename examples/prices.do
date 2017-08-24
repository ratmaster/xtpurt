// replication of results

version 11


// define environment

clear all
set more off
set mem 1000m
set maxvar 30000


// load quarterly OECD data

use prices
xtset


// compute log prices and differences

gen lprices = log(prices)
gen dlprices = D.lprices


// run unit root tests

xtpurt lprices, test(all) trend maxlags(9)
xtpurt dlprices, test(all) maxlags(9)
xtpurt dlprices, test(all) trend maxlags(9)
