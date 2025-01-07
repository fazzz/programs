#CalcRamChad.rb
#ruby CalcRamChad.rb DihedData.txt 0.01

def CalcRamChad(inputfile,step)
  linenow=[]
  histgram=[]

  pi=3.14
  num=0
  sum=0
  numhist = (pi/step).to_i*2
  sum = numhist**2
  i=0
  while i < numhist
    i=i+1
    j=0
    histgram[i]=[]
    while j < numhist
      j=j+1
      histgram[i][j]=0.0
    end
  end

  inputfile.each do |line|
    linenow[num]=line.split(" ")
#    linenow[num].shift
    cum=((linenow[num][0].to_f)/step).to_i
    row=((linenow[num][1].to_f)/step).to_i
    cum = (pi/step.to_f).to_i+cum
    row = (pi/step.to_f).to_i+row
#    p cum; p row;
#    histgram
    if histgram[cum]==nil
      histgram[cum]=[]
      histgram[cum][row]=1
      sum=sum+1
    else
      if histgram[cum][row]==nil
        histgram[cum][row]=1
        sum=sum+1
      else
        histgram[cum][row]=histgram[cum][row]+1
        sum=sum+1
      end
    end
    num=num+1
  end
  inputfile.close

#  p histgram

#  numhist = histgram[10].size



  i=0
  while i < numhist
#    psi=i*step
    j=0
    if histgram[i] != nil
      while j < numhist
#        phi=j*step
        if histgram[i][j] != nil
          hist=((histgram[i][j].to_f)/sum).to_f
          if (hist!=0.0)
            hist=Math.log(hist)
          else
            his=0.0
          end
          psi=i.to_f*step-pi
          phi=j.to_f*step-pi
          printf("%5f %5f %f\n",psi,phi,hist)
        end
        j=j+1
      end
      printf("\n")
    end
    i=i+1
  end
end

inputfile=open(ARGV.shift)
step=(ARGV.shift).to_f

CalcRamChad(inputfile,step)
