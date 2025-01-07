#CalcDistribution.rb
#ruby profile.txt 0.0 600.0 1.0 7

def CalcDistribution(inputfile,mini,maxm,step,valuetype)
  linenow=[]
  histgram=[]

  num=0
  inputfile.each do |line|
    linenow[num]=line.split(/\s+/)
    if (linenow[num][valuetype].to_f<mini)
	if histgram[0]==nil
	  histgram[0]=0
	else
          histgram[0]=histgram[0]+1
	end
#	p 0
#	p linenow[num][valuetype].to_f
    elsif (linenow[num][valuetype].to_f>maxm)
	n=((maxm-mini)/step).to_i
	if histgram[n]==nil
	  histgram[n]=0
	else
	  histgram[n]=histgram[n]+1
	end
#	p n
#	p linenow[num][valuetype].to_f
    else
	n=((linenow[num][valuetype].to_f-mini)/step).to_i
	if histgram[n]==nil
	  histgram[n]=0
	else
	  histgram[n]=histgram[n]+1
	end
#	p n
#	p linenow[num][valuetype].to_f
    end
    num=num+1
  end
  inputfile.close
#  p num

#  p histgram
	i=0
	histgram.each do |line|
	  hist=(line.to_f/num.to_f).to_f
#	  puts hist
	  print mini+i*step
	  print " "
	  print hist
	  print "\n"
	  i=i+1
	end
end

inputfile=open(ARGV.shift)
mini=(ARGV.shift).to_f
maxm=(ARGV.shift).to_f
step=(ARGV.shift).to_f
valuetype=(ARGV.shift).to_i

#p maxm
#p mini
#p step
#p valuetype

CalcDistribution(inputfile,mini,maxm,step,valuetype)

