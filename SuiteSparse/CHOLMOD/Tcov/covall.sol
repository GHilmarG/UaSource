#!/bin/csh
	tcov -x cm.profile z*.c >& /dev/null
	echo -n "statments not yet tested: "
	./covs > covs.out
	grep "#####" *tcov | wc -l
