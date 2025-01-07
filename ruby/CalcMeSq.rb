#CalcMeSq.rb

def CalcMean(inputfile,valuetype,refnum)
	num=0
	mesqdv_value_p=[]
	ref_value=0.0
	linenow=[]

	num=0
	inputfile.each do |line|
		linenow[num]=line.split(/\s+/)
		if num==refnum then
			ref_value = linenow[num][valuetype].to_f
			mesqdv_value_p[num] = 0.0
			p mesqdv_value_p[num]
		elsif num > refnum then
			mesqdv_value_p[num] = ((num-refnum-1)*mesqdv_value_p[num-1]+(linenow[num][valuetype].to_f-ref_value)**2)/(num-refnum)
			p mesqdv_value_p[num]
		end
		num = num+1
	end
	inputfile.close
end

inputfile = open(ARGV.shift)
valuetype = (ARGV.shift).to_i
refnum = (ARGV.shift).to_i
CalcMean(inputfile,valuetype,refnum)
