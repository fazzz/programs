#######################################################
#                                                     #
#                                                     #
#                                                     #
# M    M     A     K  KK   EEEEE  FFFFF IIIII   L     #
# MM  MM    A A    K K     EE     F       I     L   E #
# M MM M   AAAAA   KKK     E E    FFF     I     L     #
# M    M  AA   AA  K  KK   EEEEE  F     IIIII   LLLL  #
#                                                     #
#                                                     #
#                                                     #
#######################################################


make = make  --no-print-directory
makeclean = make clean --no-print-directory
maketest = make test --no-print-directory

default : massMD pre_MD
	$(make) -C massMD/src;
	$(make) -C pre_MD/src
	echo "Compile is Completed!"
	date >> history.txt
	mv history.txt ../man/
	cd ../man;tail -n 1 history.txt

.PHONY: clean
clean:
	$(makeclean) -C massMD/src
	$(makeclean) -C pre_MD/src

.PHONY: test
test:
	$(maketest) -C massMD/src


















