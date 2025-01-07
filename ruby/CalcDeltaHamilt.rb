#test
inputfile = open(ARGV.shift)
deltat = ARGV.shift
#p inputfile
num=0
deltahamiltonian = 0.0
hamiltonian_now=[]
linenow=[]
inputfile.each do |line|
	if num==0 then
		lineini=line.split(/\s+/)
		hamiltonian_now[0]=lineini[1].to_f
#		p line
#		p hamiltonianinitial
	elsif
		linenow[num]=line.split(/\s+/)
#		p line
#		p linenow
		hamiltonian_now[num] = linenow[num][1]
#		p hamiltonian_now[num].to_f
#		p hamiltonianinitial.to_f-(hamiltonian_now[num]).to_f
#		p hamiltonianinitial
		deltahamiltonian=deltahamiltonian+(hamiltonian_now[0]-(hamiltonian_now[num]).to_f)**2
	end
	num = num+1
end
inputfile.close
output = open("deltahamiltonian.txt","a")
output << deltat
output << " "
output << (deltahamiltonian / num)**0.5
output << " "
output.close
