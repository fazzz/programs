#CalcMean.rb

def CalcMean(inputfile,valuetype,refnum)
	num=0
	ave_value_p=0.0
	linenow=[]

	inputfile.each do |line|
		linenow[num]=line.split(' ')
		if num==refnum then
			ave_value_p = linenow[num][valuetype].to_f
		elsif num > refnum
			ave_value_p = (linenow[num][valuetype].to_f+ave_value_p*(num-refnum))/(num-refnum+1)
		end
		p ave_value_p
		num = num+1
	end
	inputfile.close
end

if (argc=ARGV.length) < 3
  print "USAGE: ./CalcMean inputfile valuetype refnum > outputfile"
  exit
end

inputfile = open(ARGV.shift)
valuetype = (ARGV.shift).to_i
refnum = (ARGV.shift).to_i
CalcMean(inputfile,valuetype,refnum)
