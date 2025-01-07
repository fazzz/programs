#test
inputfile = open(ARGV.shift)
num=0
deltakineticenergy = 0.0
kineticenergy_now=[]
linenow=[]
inputfile.each do |line|
	if num==0 then
		lineini=line.split(/\s+/)
		kineticenergy_now[0]=lineini[7].to_f
	elsif
		linenow[num]=line.split(/\s+/)
		kineticenergy_now[num] = linenow[num][1]
		deltakineticenergy=deltakineticenergy+(kineticenergy_now[0]-(kineticenergy_now[num]).to_f)**2
	end
	num = num+1
end
inputfile.close
output = open("deltaKineticenergy.txt","a")
output <<  (deltakineticenergy / num)**0.5
output.close
