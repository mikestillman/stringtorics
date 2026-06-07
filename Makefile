M2=~/src/M2-worktree/computegv/M2/BUILD/mike/builds.tmp/cmake-appleclang/M2

install:
	$(M2) -e 'uninstallAllPackages()' -e 'exit 0'
	$(M2) -e 'installPackage "PALPInterface"' -e '-exit 0'
	$(M2) -e 'installPackage "IntegerEquivalences"' -e '-exit 0'
	$(M2) -e 'installPackage "DanilovKhovanskii"' -e '-exit 0'
	$(M2) -e 'installPackage "StringTorics"' -e '-exit 0'

check:
	$(M2) -e 'check "PALPInterface"' -e '-exit 0'
	$(M2) -e 'check "IntegerEquivalences"' -e '-exit 0'
	$(M2) -e 'check "DanilovKhovanskii"' -e '-exit 0'
	$(M2) -e 'check "StringTorics"' -e '-exit 0'

clean:
	rm -rf foo-normalForm-foo
