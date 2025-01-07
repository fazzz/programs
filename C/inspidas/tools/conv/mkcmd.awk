{
	print "cut -c1-72 <",$1," >tmp;sed -e 's/ *$//' <tmp"," >",$2
	print "sed -e 's/^      INCLUDE/CMSP  INCLUDE/' <",$2,">tmp"
	print "mv tmp ",$2
}
