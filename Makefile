#You may change M2 to point to the executable M2 on your machine, e.g. M2=/opt/homebrew/bin/M2
M2=M2

install:
	$(M2) -e 'uninstallAllPackages()' -e 'exit 0'
	$(M2) -e 'elapsedTime installPackage "PALPInterface"' -e '-exit 0' # 
	$(M2) -e 'elapsedTime installPackage "IntegerEquivalences"' -e '-exit 0'
	$(M2) -e 'elapsedTime installPackage "DanilovKhovanskii"' -e '-exit 0'
	$(M2) -e 'elapsedTime installPackage "StringTorics"' -e '-exit 0'

check:
	$(M2) -e 'elapsedTime check "PALPInterface"' -e '-exit 0'
	$(M2) -e 'elapsedTime check "IntegerEquivalences"' -e '-exit 0'
	$(M2) -e 'elapsedTime check "DanilovKhovanskii"' -e '-exit 0'
	$(M2) -e 'elapsedTime check "StringTorics"' -e '-exit 0'

clean:
	rm -rf foo-normalForm-foo foo run-palp*

