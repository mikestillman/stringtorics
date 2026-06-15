# stringtorics

The package **StringTorics** is ready to use, even though it is being actively worked on.
If you run across bugs, 

The package **IntegerEquivalences** is pretty much usable, although one example and test fail.

The package **DanilovKhovanskii** is still being debugged and cleaned up.  Do **not**
trust results from this package without checking them!

The package **PALPInterface** is in its infancy.  It doesn't interface to too many
routines in the PALP package, although this month at ESI, I hope that will change!

# Installing StringTorics

First install the latest Macaulay2 (version 1.26.06, not 1.26.05!).  
For instructions, see https://github.com/Macaulay2/M2/wiki

Second, download, or checkout, the repository
  https://github.com/mikestillman/stringtorics
via
  `git clone https://github.com/mikestillman/stringtorics`

Second, change into the stringtorics directory, and run M2.
```
cd stringtorics`
`m2`
```

Now install StringTorics as follows

```m2
installPackage "StringTorics"
```

You can check that it is correct with 
```m2
check needsPackage "StringTorics"
```

You then can use it (and `IntegerEquivalences`, and `PALPInterface`)
by doing (inside M2):

```m2
needsPackage "StringTorics"
```

For examples and documentation (which we are improving during this month!)
please try

```m2
viewHelp "StringTorics"
```

or 

```m2
help "StringTorics"
```
