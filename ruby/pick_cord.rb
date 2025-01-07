#!~/ruby/bin/ruby

inputfile = open(ARGV.shift)
numatom = (ARGV.shift).to_i

#############################################
i=0

coord=[]
thisline=[]
num=0
j=1
inputfile.each do |line|
  thisline = line.split(' ')
  coord << thisline
#  p coord
  coord.flatten!

  if coord.size==numatom*3
    num=num+1
#    p num
    judge=num%100
#    p judge

    out = File.open(name,"w")
    if judge==0
      i=0
      out.printf("ALA\n")
      out.printf("    %d\n",numatom)
      while i < numatom
        out.printf("%12.7f%12.7f%12.7f",coord[i*3],coord[i*3+1],coord[i*3+2])
        i=i+1
        if (i%2==0)
          out.printf("\n")
        end
      end
      printf("\n")
    out.close
    j=j+1
   end
     coord.clear
  end
end

#############################################
