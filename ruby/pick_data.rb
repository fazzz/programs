inputfile = open(ARGV.shift)
numatom = (ARGV.shift).to_i

#############################################
i=0

coord=[]
thisline=[]
num=0

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
    if judge==0
      i=0
      while i < numatom
        printf("%f %f %f\n",coord[i*3],coord[i*3+1],coord[i*3+2])
        i=i+1
      end
      printf("\n")
    end
    coord.clear
  end
end

#############################################
