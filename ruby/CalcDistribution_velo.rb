#CalcDistribution.rb
#ruby profile.txt 0.0 600.0 1.0 7

def CalcDistribution_velo(inputfile,max,step,valuetype,atom,numatom)
  linenow=[]
  histgram=[]

  num=0
  inputfile.each do |line|
    linenow[num]=line.split(/\s+/)
    if ((num%(numatom+2))==atom)
       if((linenow[num][valuetype].to_f).abs>maxm)
#	n=((maxm)/step).to_i
#	if histgram[n]==nil
#	  histgram[n]=0
#	else
#	  histgram[n]=histgram[n]+1
#	end
#	p n
#	p linenow[num][valuetype].to_f
       else
         n=(((linenow[num][valuetype].to_f)+max)/step).to_i
#         p linenow[num][valuetype].to_f
#         p n
	if histgram[n]==nil
	  histgram[n]=0
	else
	  histgram[n]=histgram[n]+1
	end
#	p n
#	p linenow[num][valuetype].to_f
#    end
    end
    num=num+1
  end
  inputfile.close
#  p num

#  p histgram
	i=0
	histgram.each do |line|
          if (line==nil)
            line=0
          end
	  hist=(line.to_f/num.to_f).to_f
#	  puts hist
#          if (i<max)
            print -max+i*step            
#          end
#          else 
#            print -max+i*step            
#          end
	  print " "
	  print hist
	  print "\n"
	  i=i+1
	end
end

inputfile=open(ARGV.shift)
maxm=(ARGV.shift).to_f
step=(ARGV.shift).to_f
valuetype=(ARGV.shift).to_i
atom=(ARGV.shift).to_i
numatom=(ARGV.shift).to_i

#p maxm
#p mini
#p step
#p valuetype

CalcDistribution_velo(inputfile,max,step,valuetype,atom,numatom)

