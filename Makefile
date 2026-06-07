M2=~/src/M2-worktree/computegv/M2/BUILD/mike/builds.tmp/arm64-appleclang/M2

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
	rm -rf foo-normalForm-foo
