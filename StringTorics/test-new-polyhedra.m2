restart
uninstallPackage "CohomCalg"
uninstallPackage "Topcom"
uninstallPackage "ReflexivePolytopesDB"
uninstallPackage "AbstractToricVarieties"
uninstallPackage "StringTorics"

restart
check "CohomCalg"
check "Topcom"
check "ReflexivePolytopesDB"
check "AbstractToricVarieties"
check "StringTorics" --  currently (16 May 2018), one failure, due to changes in triangulation code.




path = prepend("../m2-code/", path)

uninstallPackage "Polyhedra"
needsPackage("Polyhedra", FileName=>"~/src/M2-lkastner/M2/Macaulay2/packages/Polyhedra.m2")
restart
installPackage("Polyhedra", FileName=>"~/src/M2-lkastner/M2/Macaulay2/packages/Polyhedra.m2")

needsPackage "StringTorics"
peek loadedFiles

