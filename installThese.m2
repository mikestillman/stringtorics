-- To install StringTorics:
-- Do these lines in order.
restart
uninstallAllPackages()
restart
installPackage "IntegerEquivalences"

-- Install PALPInterface BEFORE StringTorics: some StringTorics doc examples
-- (e.g. label(CalabiYauInToric)) need PALPInterface to be loadable, otherwise
-- installing StringTorics reports "file not found on path: PALPInterface.m2".
-- Loading PALPInterface does not require the PALP binaries: an executable is
-- looked up only when a function actually needs one.  Installing it does, since
-- its documentation examples run PALP.  Put the palp executables (poly.x,
-- cws.x, nef.x, and the -11d.x variants) on your PATH, or add their directory
-- to programPaths.  On macOS:  brew install macaulay2/tap/palp
restart
installPackage "PALPInterface"
restart
installPackage "StringTorics"
check "IntegerEquivalences" -- test # 0 takes 14 sec on my M4 MBP
check "PALPInterface" -- 19 tests, all pass.
check "StringTorics"


-- At this point, GV invariants functions won't work yet, but everything else should.
-- To get GV invariants going, compile the code in ComputeGV
-- and place computeGV on your PATH.

-- to install DanilovKhovanskii (warning: although all examples and tests install and check,
-- this package has not been well tested!  Also the interface will change, perhaps
-- even this month at ESI)
restart
uninstallPackage "DanilovKhovanskii"
restart
installPackage "DanilovKhovanskii"
restart
check "DanilovKhovanskii"
