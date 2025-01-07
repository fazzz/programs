#CalcRtMeSq.rb
#usage: ruby CalcRtMeSq.rb inputfile valuetype refnum
#last modified 3III10

def CalcRtMeSq(inputfile,valuetype,refnum)
        rt_me_sq_value=0.0
        squre_ave_value=0.0
        ave_value=0.0
	linenow=[]

	num=0
	inputfile.each do |line|
		linenow[num]=line.split(' ')
		if num==refnum then
                        ave_value = linenow[num][valuetype].to_f
                        squre_ave_value = (linenow[num][valuetype].to_f)**2
			rt_me_sq_value = 0.0
                        p rt_me_sq_value
		elsif num > refnum then
			ave_value = ((num-refnum)*ave_value+linenow[num][valuetype].to_f)/(num-refnum+1)
			squre_ave_value = ((num-refnum)*squre_ave_value+(linenow[num][valuetype].to_f)**2)/(num-refnum+1)
			rt_me_sq_value = ((squre_ave_value-ave_value**2)/(ave_value**2))**0.5 
                        p rt_me_sq_value
		end
		num = num+1
	end
	inputfile.close
end

inputfile = open(ARGV.shift)
valuetype = (ARGV.shift).to_i
refnum = (ARGV.shift).to_i
CalcRtMeSq(inputfile,valuetype,refnum)
