# stringtorics

The package **StringTorics** is ready to use, even though it is being actively worked on.
If you run across bugs, please let me know (mes15@cornell.edu)

The package **IntegerEquivalences** is usable.

The package **PALPInterface** is in its infancy.  It doesn't interface to too many
routines in the PALP package, although this month at ESI, I hope that will change!

The package **DanilovKhovanskii** is still being debugged and cleaned up.  Do **not**
trust results from this package without checking them!  However, all documentation
examples and all tests install and check cleanly.

# Installing StringTorics

First install the latest Macaulay2 (version 1.26.06, not 1.26.05!).  
For instructions, see https://github.com/Macaulay2/M2/wiki

Second, download, or checkout, the repository, and change directoy into this directory.
  https://github.com/mikestillman/stringtorics
via
```m2
git clone https://github.com/mikestillman/stringtorics
cd stringtorics

```

At this point, run M2, and run the following commands to install, check, or view the documentation.
    You may follow these instructions replacing `StringTorics` with `IntegerEquivalences`, or `PALPInterfac`.
    (or even `DanilovKhovanskii`, but recall that it hasn't been checked carefull for correctness yet!).

```
M2
```

Now (inside running M2), install StringTorics as follows

```m2
installPackage "StringTorics"
```

You can check that it is correct with 
```m2
check needsPackage "StringTorics"
```

You then can use it (and `IntegerEquivalences`, and `PALPInterface`)
by doing (inside M2, as above):

```m2
needsPackage "StringTorics"
```

For examples and documentation (which we are improving during this month!)
please try (inside M2, as above):

```m2
viewHelp "StringTorics"
```

or 

```m2
help "StringTorics"
```
