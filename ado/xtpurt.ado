*! version 1.1.0  24aug2017
//******************************************************************************
// xtpurt: A Comprehensive Testing Framework for Panel Unit Root Tests
// -------------------------------------------------------------------
//
// This Stata package provides easy access to the following panel unit root
// tests (PURTs):
//
//
// - Herwartz and Siedenburg (2008)
// - Demetrescou and Hanck (2012)
// - Herwartz, Maxand and Walle (2017)
//
// For further information, please type "help xtpurt".
//
//
// Implementation:
//
// Fabian H. C. Raters
// Yabibal M. Walle
// Department of Economics
// University of Goettingen
//
// For user convenience, some code pieces are quite similar to xtunitroot().
//
//
// (C) 2017 by the authors.
//******************************************************************************

program xtpurt, rclass sortpreserve

version 11

syntax varname(ts) [if] [in], [ TEST(string)     ///
                                Trend            ///
                                NOCONstant       ///
                                LAGSEL(string)   ///
                                MAXlags(string)  ///
                                LAGSEXport(name) ///
                                LAGS(string)     ///
                              ]

tsunab ylevel : `varlist'
marksample touse

if "`test'" != "" {
	if "`test'" != "hs" && "`test'" != "dh" && "`test'" != "hmw"  ///
			      && "`test'" != "all" {
		di as error "panel unit root test {bf:`test'} is unknown"
		exit 198
	}
}
else {
	// set default tests
	if "`trend'" != "" {
		local test = "hmw"
	}
	else {
		local test = "hs"
	}
}

// errors and warnings for combinations of deterministics and tests
local warnings = ""
if "`trend'" != "" {
	if "`noconstant'" != "" {
		di in smcl as err "cannot specify both {bf:noconstant} and {bf:trend}"
		exit 198
	}
	
	// deterministics for lag order selection and prewhitening
	local noconst = ""

	if "`test'" == "all" {
		local warnings = "t_hs and t_dh are not robust to heteroskedasticity with {bf:trend}"
	}
	else if "`test'" != "hmw" {
		local warnings = "t_`test' is not robust to heteroskedasticity with {bf:trend}"
	}
}
else {
	// deterministics for lag order selection and prewhitening
	local noconst = "noconstant"

	if "`test'" == "hmw" {
		di in smcl as err "cannot apply panel unit root test {bf:hmw} without {bf:trend}"
		exit 198
	}
}

if "`lagsel'" != "" {
	if "`lags'" != "" {
		di as error "cannot specify both {bf:lags()} and {bf:lagsel()}"
		exit 198
	}
	if "`lagsel'" == "aic" {
		local lagcriterion = "AIC"
	}
	else if "`lagsel'" == "bic" {
		local lagcriterion = "BIC"
	}
	else if "`lagsel'" == "hqic" {
		local lagcriterion = "HQIC"
	}
	else {
		di as error "lag selection method {bf:`lagsel'} is unknown"
		exit 198
	}
}
else {
	if "`lags'" != "" {
		local lagsel = "fix"
		local lagcriterion = "FIX"
	}
	else {
		local lagsel = "bic"
		local lagcriterion = "BIC"
	}
}

if "`maxlags'" != "" {
	if "`lags'" != "" {
		di as error "cannot specify both {bf:lags()} and {bf:maxlags()}"
		exit 198
	}
	if `maxlags' < 2 {
		di as error "{bf:maxlags()} must be at least 2"
		exit 198
	}
}
else {
	local maxlags = "4"
}


_xt, trequired
capture xtset
tempname usrdelta
local panelvar      `r(panelvar)'
local timevar       `r(timevar)'
scalar `usrdelta' = `r(tdelta)'

// check for balanced sample without missing data
qui _xtstrbal `panelvar' `timevar' `touse'
if r(strbal) == "no" {
	di as error "test requires strongly balanced data"
	exit 498
}


// work with data
preserve
quietly {
drop if `touse' == 0

// do not allow gaps
tempname timediff
by `panelvar': gen `timediff' = D.`timevar'
capture assert `timediff' == `usrdelta' if `timevar' != `timevar'[1]
if c(rc) {
	di as error "test does not allow for gaps in data"
	exit 498
}

// generate differences
tempname ydiff
gen double `ydiff' = D.`ylevel'

// store panel levels
qui levelsof `panelvar', local(ggs)

// check for constant values
foreach i of local ggs {
	summarize `ydiff' if `panelvar' == `i'
	if r(min) == r(max) {
		di as error "the first differences of {bf:`ylevel'} for" ///
		" {bf:`panelvar'}=`i' are constant"
		exit 498
	}
}


// determine lagorder
tempvar lagorder
if "`lagsel'" != "fix" {
	// receive optimal lag-order
	gen int `lagorder' = .
	tempname llz lorder
	scalar `llz' = .
	if "`noconst'" != "" {
		tempvar ydiff2
		gen double `ydiff2' = `ydiff'^2
		by `panelvar': replace `ydiff2' = . if _n <= `maxlags' + 1
	}
	foreach i of local ggs {
		// compute information criteria for zero lags and noconstant
		if "`noconst'" != "" {
			summarize `ydiff2' if `panelvar' == `i', meanonly
			scalar `llz' = ln(r(mean)) + ln(2 * _pi) + 1
		}
	
		varsoc `ydiff' if `panelvar' == `i', maxlag(`maxlags') `noconst'
		mata getLagorder("`lagsel'", "`llz'", "`lorder'")
		replace `lagorder' = `lorder' if `panelvar' == `i'
	}
}
else {
	gen int `lagorder' = `lags'
	summarize `lagorder', meanonly
	if r(min) < 0 {
		di as error "{bf:lags()} must be nonnegative"
		exit 198
	}
}
mkmat `lagorder' if `timevar' == `timevar'[1], matrix(lom)
mat colnames lom = "Lags"
mat rownames lom = `ggs'


// prewhitening
tempvar ylevelpre ydiffpre
gen double `ylevelpre' = `ylevel'
gen double `ydiffpre' = `ydiff'
tempname beta
local row = 1
foreach i of local ggs {
	local lo = lom[`row', 1]
	if `lo' > 0 {
		regress `ydiff' L(1/`lo').`ydiff' if `panelvar' == `i', `noconst'
		matrix `beta' = e(b)
		forvalues j = 1/`lo' {
			replace `ylevelpre' = `ylevelpre' - `beta'[1,`j'] * L`j'.`ylevel' ///
						if `panelvar' == `i'
			replace `ydiffpre' = `ydiffpre' - `beta'[1,`j'] * L`j'.`ydiff' ///
						if `panelvar' == `i'
		}
	}
	local row = `row' + 1
}
replace `ylevel' = `ylevelpre'
replace `ydiff' = `ydiffpre'
drop `ylevelpre' `ydiffpre'


// re-balance panel
tempname maxlag minlag
summarize `lagorder', meanonly
scalar `maxlag' = r(max)
scalar `minlag' = r(min)
by `panelvar': drop if _n <= `maxlag' + 1


if "`trend'" != "" {
	// recursive detrending according to Chang (2002)
	by `panelvar': replace `ylevel' = `ylevel' + (2 / _n) * sum(`ylevel') ///
				- 6 / ((_n + 1) * _n) * sum(_n * `ylevel')
	
	// first two differences might not be exactly zero
	// which affects DH via the sign() transformations
	replace `ylevel' = 0 if `timevar' == `timevar'[1] | `timevar' == `timevar'[2]
	
	foreach i of local ggs {
		summarize `ydiff' if `panelvar' == `i', meanonly
		replace `ydiff' = `ydiff' - r(mean) if `panelvar' == `i'
	}
}
else {
	// demeaning by subtracting the first observation
	if "`noconstant'" == "" {
		tempname firstelem
		by `panelvar': gen double `firstelem' = `ylevel'[1]
		by `panelvar': replace `ylevel' = `ylevel' - `firstelem'
	}
}

}

// apply testing method
tempname tau
if "`test'" == "hmw" || ("`test'" == "all" && "`trend'" != "") {
	// applying HMW test
	mata getHMW("`timevar'", "`ylevel'", "`ydiff'", "`tau'")
	return scalar t_hmw_p = normal(`tau')
	return scalar t_hmw = `tau'
	local testname = "Herwartz et al. (2017)"
}
if "`test'" == "dh" || "`test'" == "all" {
	// applying DH test
	mata getHSDH("`timevar'", "`ylevel'", "`ydiff'", "`tau'", 1)
	return scalar t_dh_p = normal(`tau')
	return scalar t_dh = `tau'
	local testname = "Demetrescou and Hanck (2012)"
}
if "`test'" == "hs" || "`test'" == "all" {
	// applying HS test
	mata getHSDH("`timevar'", "`ylevel'", "`ydiff'", "`tau'", 0)
	return scalar t_hs_p = normal(`tau')
	return scalar t_hs = `tau'
	local testname = "Herwartz and Siedenburg (2008)"
}
if "`test'" == "all" {
	local testname = "All methods"
}


// clean up after computations
restore

// optional lagorder export
if ("`lagsexport'" != "") {
	qui gen int `lagsexport' = .
	local row = 1
	foreach i of local ggs {
		qui replace `lagsexport' = lom[`row', 1] if `panelvar' == `i'
		local row = `row' + 1
	}
}


// define return values
qui xtsum `panelvar' if `touse'
return scalar N_r  = r(N) / r(n) - `maxlag' - 1
return scalar N_t  = r(N) / r(n)
return scalar N_g  = r(n)
return scalar N    = r(N)
return matrix lags = lom
return local warnings = "`warnings'"
return local lagsel = "`lagsel'"
if "`trend'`noconstant'" != "" {
	return local determ "`trend'`noconstant'"
}
else {
	return local determ "constant"
}
return local test = "`test'"


// display output table
di
di as text "`testname' unit-root test for " as res "`ylevel'"
local linelen = 21 + length("`testname'") + length("`ylevel'")
di in smcl as text "{hline `linelen'}"
di as text "Ho: Panels contain unit roots"	///
	_col(45) "Number of panels"		///
	_col(63) "="				///
	_col(64) as res %7.0g return(N_g)
di as text "Ha: Panels are stationary"		///
	_col(45) "Number of periods"		///
	_col(63) "="				///
	_col(64) as res %7.0g return(N_t)
di as text 					///
	_col(45) "After rebalancing"		///
	_col(63) "="				///
	_col(64) as res %7.0g return(N_r)
di
di as text "Constant:   " _c
if "`noconstant'" != "" {
	di as res "Not included" _c
}
else {
	di as res "Included" _c
}
di as text _col(45) "Prewhitening: " _c
di as res "`lagcriterion'"
di as text "Time trend: " _c
if "`trend'" != "" {
	di as res "Included" _c
}
else {
	di as res "Not included" _c
}
di as text _col(45) "Lag orders: " _c
di as res  "  min=" as res %1.0f `minlag' " max=" as res %1.0f `maxlag'
di in smcl as text "{hline 78}"
di as text _col(2) "Name" _col(21) "Statistic" _col(36) "p-value"
di in smcl as text "{hline 78}"
if "`test'" == "hs" || "`test'" == "all" {
	di as text _col(2) "t_hs"			///
		as res _col(21) %9.4f return(t_hs)	///
		as res _col(37) %6.4f return(t_hs_p)
}
if "`test'" == "dh" || "`test'" == "all" {
	di as text _col(2) "t_dh"			///
		as res _col(21) %9.4f return(t_dh)	///
		as res _col(37) %6.4f return(t_dh_p)
}
if "`test'" == "hmw" || ("`test'" == "all" && "`trend'" != "") {
	di as text _col(2) "t_hmw"			///
		as res _col(21) %9.4f return(t_hmw)	///
		as res _col(37) %6.4f return(t_hmw_p)
}
if "`warnings'" != "" {
	di
	di in smcl as text "{hline 6}"
	di as text "Caution:    {it:`warnings'}"
}
di in smcl as text "{hline 78}"

end


mata
// define lag order subfunction
void getLagorder(string scalar lagsel, string scalar llzn, string scalar lordern) {
	X = st_matrix("r(stats)")
	z = st_numscalar(llzn)
	if (lagsel == "aic") {
		a = 7
	}
	else if (lagsel == "bic") {
		a = 9
	}
	else {
		a = 8
	}
	mv = colmin(X[, a])
	if (z == . || mv < z) {
		y = select(X[, 1], X[, a] :== mv)
	}
	else {
		y = 0
	}
	st_numscalar(lordern, y)
}


// define HS and DH test
void getHSDH(string scalar timevarn, string scalar yleveln, string scalar ydiffn, string scalar taun, real scalar dh) {
	timevar = st_data(., timevarn)
	T = rows(uniqrows(timevar))

	ylevel = st_data(., yleveln)
	ydiff = st_data(., ydiffn)
	yLevel = colshape(ylevel, T)'
	if (dh == 1) {
		yLevel = sign(yLevel)
	}
	yDiff = colshape(ydiff, T)'

	// compute HS or DH test statistic
	x = colshape(yLevel[1..T-1, .]', 1)' * colshape(yDiff[2..T, .]', 1)
	z = 0
	for (t = 2; t <= T; t++) {
		z = z + yLevel[t - 1, .] * yDiff[t, .]' * yDiff[t, .] * yLevel[t - 1, .]'
	}
	tau = x / sqrt(z)

	// return estimated tau
	st_numscalar(taun, tau)
}


// define HMW test
void getHMW(string scalar timevarn, string scalar yleveln, string scalar ydiffn, string scalar taun) {
	timevar = st_data(., timevarn)
	T = rows(uniqrows(timevar))

	ylevel = st_data(., yleveln)
	ydiff = st_data(., ydiffn)
	yLevel = colshape(ylevel, T)'
	yDiff = colshape(ydiff, T)'
	yDiff2 = yDiff * yDiff'
	yDiff22 = yDiff2 :* yDiff2

	x = colshape(yLevel[1..T-1, .]', 1)' * colshape(yDiff[2..T, .]', 1)

	a = J(T, T, 0)
	nu = J(T, 1, 0)
	for (t = 2; t <= T; t++) {
		v = 0
		for (i = 1; i <= t - 1; i++) {
			a[i, t - 1] = 1 + 2 * (t - i) / (t - 1) ///
					- 3 * (1 - ((i - 1) * i) / ((t - 1) * t))
			v = v + (a[i, t - 1] / T) * yDiff2[i, i]
		}
		nu[t] = -v
	}
	ab = a / T
	at = a - ab


	// compute s
	c1 = 0
	for (i = 1; i <= T - 1; i++) {
		for (j = i + 1; j <= T - 1; j++) {
			for (s = i + 1; s <= T; s++) {
				for (t = j + 1; t <= T; t++) {
					c1 = c1 + ab[i, s - 1] * ab[j, t - 1] * yDiff22[i, j]
				}
			}
		}
	}
	c1 = 2 * c1

	c2 = 0
	for (i = 1; i <= T - 1; i++) {
		for (s = i + 1; s <= T; s++) {
			for (t = i + 1; t <= T; t++) {
				c2 = c2 + at[i, s - 1] * ab[i, t - 1] * yDiff22[i, s]
			}
		}
	}
	c2 = 2 * c2

	c3 = 0
	for (i = 1; i <= T - 1; i++) {
		for (t = i + 1; t <= T; t++) {
			c3 = c3 + at[i, t - 1]^2 * yDiff22[i, t]
		}
	}

	c4 = 0
	for (i = 1; i <= T - 1; i++) {
		for (t = i + 1; t <= T; t++) {
			z = 0
			for (j = 1; j <= T; j++) {
				if (j != i && j != t) {
					z = z + yDiff22[i, j]
				}
			}
			c4 = c4 + z * ab[i, t - 1]^2
		}
	}

	c5 = 0
	for (i = 1; i <= T - 2; i++) {
		for (s = i + 1; s <= T - 1; s++) {
			for (t = s + 1; t <= T; t++) {
				z = 0
				for (j = 1; j <= T; j++) {
					if (j != i && j != t && j != s) {
						z = z + yDiff22[i, j]
					}
				}
				c5 = c5 + z * ab[i, t - 1] * ab[i, s - 1]
			}
		}
	}
	c5 = 2 * c5

	// compute HMW test statistic
	tau = (x - sum(nu)) / sqrt(c1 - c2 + c3 + c4 + c5)

	// return estimated tau
	st_numscalar(taun, tau)
}
end
