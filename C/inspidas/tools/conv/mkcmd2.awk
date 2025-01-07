{
	print "sed -e 's/ *$//' <",$1," >tmp"
	print "mv tmp ",$1
}
