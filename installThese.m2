-- To install StringTorics:
-- Do these lines in order.
restart
uninstallAllPackages()
restart
installPackage "IntegerEquivalences"

-- Install PALPInterface BEFORE StringTorics: some StringTorics doc examples
-- (e.g. label(CalabiYauInToric)) need PALPInterface to be loadable, otherwise
-- installing StringTorics reports "file not found on path: PALPInterface.m2".
-- Note: to actually USE the PALP functions you must also install PALP on your
-- computer (place the palp executables, e.g. poly.x, cws.x, on your PATH).
-- Installing/loading the package itself does not require the binaries.
-- (PALPInterface has no documentation or tests at all.)
restart
installPackage "PALPInterface"
restart
installPackage "StringTorics"
check "IntegerEquivalences"
check "StringTorics"
check "PALPInterface"

-- At this point, GV invariants functions won't work yet, but everything else should.
-- To get GV invariants going, compile the code in ComputeGV
-- and place computeGV on your PATH.

-- to install DanilovKhovanskii
restart
uninstallPackage "DanilovKhovanskii"
restart
installPackage "DanilovKhovanskii"
restart
check "DanilovKhovanskii"  -- one error (hash table has wrong class) 3/5/26.
