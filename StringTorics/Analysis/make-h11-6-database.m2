-- h11=6 database construction, 20 Aug 2023.
-- We have already constructed (taking several hours on 'heaviside'
-- the polytopes database:
--   This was done by creatng several databases, then combining them.
-- Here we add in the CY3s, but to different files

  restart
  debug needsPackage "StringTorics"
  DB6 = "../Databases/cys-ntfe-h11-6.dbm"
  addToCYDatabase("Foo6/h11-6-test.dbm", DB6, {0,1})
  
  

  topes = kreuzerSkarke(5, Limit => 20000); -- 4990 of these
  assert(#topes == 4990)
  elapsedTime createCYDatabase(DB5, topes)
  -- Let's find which are not favorable, not torsion free.
  -- This took 6300 seconds (a bit less than 2 hours).
  topes2 = drop(topes,4873)
  elapsedTime createCYDatabase(DB5, topes2)
  -- OK made it to here.
    
end--
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=0' -e 'hi=100' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=101' -e 'hi=200' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=201' -e 'hi=300' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=301' -e 'hi=400' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=401' -e 'hi=500' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=501' -e 'hi=600' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=601' -e 'hi=700' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=701' -e 'hi=800' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=801' -e 'hi=900' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=901' -e 'hi=1000' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=1001' -e 'hi=1100' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=1101' -e 'hi=1200' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=1201' -e 'hi=1300' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=1301' -e 'hi=1400' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=1401' -e 'hi=1500' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=1501' -e 'hi=1600' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=1601' -e 'hi=1700' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=1701' -e 'hi=1800' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=1801' -e 'hi=1900' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=1901' -e 'hi=2000' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=2001' -e 'hi=2100' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=2101' -e 'hi=2200' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=2201' -e 'hi=2300' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=2301' -e 'hi=2400' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=2401' -e 'hi=2500' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=2501' -e 'hi=2600' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=2601' -e 'hi=2700' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=2701' -e 'hi=2800' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=2801' -e 'hi=2900' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=2901' -e 'hi=3000' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=3001' -e 'hi=3100' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=3101' -e 'hi=3200' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=3201' -e 'hi=3300' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=3301' -e 'hi=3400' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=3401' -e 'hi=3500' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=3501' -e 'hi=3600' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=3601' -e 'hi=3700' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=3701' -e 'hi=3800' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=3801' -e 'hi=3900' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=3901' -e 'hi=4000' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=4001' -e 'hi=4100' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=4101' -e 'hi=4200' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=4201' -e 'hi=4300' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=4301' -e 'hi=4400' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=4401' -e 'hi=4500' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=4501' -e 'hi=4600' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=4601' -e 'hi=4700' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=4701' -e 'hi=4800' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=4801' -e 'hi=4900' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=4901' -e 'hi=5000' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=5001' -e 'hi=5100' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=5101' -e 'hi=5200' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=5201' -e 'hi=5300' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=5301' -e 'hi=5400' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=5401' -e 'hi=5500' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=5501' -e 'hi=5600' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=5601' -e 'hi=5700' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=5701' -e 'hi=5800' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=5801' -e 'hi=5900' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=5901' -e 'hi=6000' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=6001' -e 'hi=6100' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=6101' -e 'hi=6200' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=6201' -e 'hi=6300' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=6301' -e 'hi=6400' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=6401' -e 'hi=6500' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=6501' -e 'hi=6600' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=6601' -e 'hi=6700' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=6701' -e 'hi=6800' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=6801' -e 'hi=6900' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=6901' -e 'hi=7000' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=7001' -e 'hi=7100' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=7101' -e 'hi=7200' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=7201' -e 'hi=7300' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=7301' -e 'hi=7400' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=7401' -e 'hi=7500' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=7501' -e 'hi=7600' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=7601' -e 'hi=7700' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=7701' -e 'hi=7800' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=7801' -e 'hi=7900' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=7901' -e 'hi=8000' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=8001' -e 'hi=8100' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=8101' -e 'hi=8200' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=8201' -e 'hi=8300' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=8301' -e 'hi=8400' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=8401' -e 'hi=8500' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=8501' -e 'hi=8600' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=8601' -e 'hi=8700' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=8701' -e 'hi=8800' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=8801' -e 'hi=8900' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=8901' -e 'hi=9000' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=9001' -e 'hi=9100' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=9101' -e 'hi=9200' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=9201' -e 'hi=9300' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=9301' -e 'hi=9400' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=9401' -e 'hi=9500' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=9501' -e 'hi=9600' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=9601' -e 'hi=9700' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=9701' -e 'hi=9800' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=9801' -e 'hi=9900' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=9901' -e 'hi=10000' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=10001' -e 'hi=10100' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=10101' -e 'hi=10200' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=10201' -e 'hi=10300' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=10301' -e 'hi=10400' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=10401' -e 'hi=10500' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=10501' -e 'hi=10600' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=10601' -e 'hi=10700' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=10701' -e 'hi=10800' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=10801' -e 'hi=10900' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=10901' -e 'hi=11000' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=11001' -e 'hi=11100' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=11101' -e 'hi=11200' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=11201' -e 'hi=11300' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=11301' -e 'hi=11400' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=11401' -e 'hi=11500' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=11501' -e 'hi=11600' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=11601' -e 'hi=11700' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=11701' -e 'hi=11800' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=11801' -e 'hi=11900' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=11901' -e 'hi=12000' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=12001' -e 'hi=12100' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=12101' -e 'hi=12200' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=12201' -e 'hi=12300' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=12301' -e 'hi=12400' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=12401' -e 'hi=12500' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=12501' -e 'hi=12600' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=12601' -e 'hi=12700' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=12701' -e 'hi=12800' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=12801' -e 'hi=12900' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=12901' -e 'hi=13000' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=13001' -e 'hi=13100' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=13101' -e 'hi=13200' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=13201' -e 'hi=13300' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=13301' -e 'hi=13400' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=13401' -e 'hi=13500' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=13501' -e 'hi=13600' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=13601' -e 'hi=13700' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=13701' -e 'hi=13800' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=13801' -e 'hi=13900' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=13901' -e 'hi=14000' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=14001' -e 'hi=14100' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=14101' -e 'hi=14200' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=14201' -e 'hi=14300' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=14301' -e 'hi=14400' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=14401' -e 'hi=14500' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=14501' -e 'hi=14600' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=14601' -e 'hi=14700' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=14701' -e 'hi=14800' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=14801' -e 'hi=14900' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=14901' -e 'hi=15000' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=15001' -e 'hi=15100' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=15101' -e 'hi=15200' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=15201' -e 'hi=15300' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=15301' -e 'hi=15400' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=15401' -e 'hi=15500' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=15501' -e 'hi=15600' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=15601' -e 'hi=15700' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=15701' -e 'hi=15800' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=15801' -e 'hi=15900' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=15901' -e 'hi=16000' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=16001' -e 'hi=16100' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=16101' -e 'hi=16200' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=16201' -e 'hi=16300' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=16301' -e 'hi=16400' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=16401' -e 'hi=16500' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=16501' -e 'hi=16600' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=16601' -e 'hi=16700' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=16701' -e 'hi=16800' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=16801' -e 'hi=16900' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=16901' -e 'hi=17000' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=17001' -e 'hi=17100' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0' &

///
M2 --stop -e 'needsPackage "StringTorics"' -e 'lo=0' -e 'hi=2' -e 'addToCYDatabase("Foo6/range-"|lo|"-"|hi|".dbm", "../Databases/cys-ntfe-h11-6.dbm", splice{lo..hi})' -e 'exit 0'

lo = 0;
hi = -1;
for i from 1 to 171 do (
    lo = hi + 1;
    hi = 100 * i;
    m2str = "M2 --stop -e 'needsPackage \"StringTorics\"' -e 'lo="|lo|"' -e 'hi="|hi|"' -e 'addToCYDatabase(\"Foo6/range-\"|lo|\"-\"|hi|\".dbm\", \"../Databases/cys-ntfe-h11-6.dbm\", splice{lo..hi})' -e 'exit 0' &";
    print m2str;
    )
///
