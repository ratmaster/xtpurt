// replication of results

version 11


// define environment

clear all
set more off
set mem 1000m
set maxvar 30000


// run example

use prices

xtset

xtpurt lprices, test(all) trend maxlags(9)

xtpurt dlprices, test(all) maxlags(9)

xtpurt dlprices, test(all) trend maxlags(9)
