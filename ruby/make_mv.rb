#!/home/yamamori/ruby/bin

inpfile=open(ARGV.shift)
totalnatom=(ARGV.shift).to_i
nameatom=["C ","H ","H ","H ","C ","O ","N ","H ","C ","H ","C ","H ","H ","H ","C ","O ","N ","H ","C ","H ","H ","H "]

f=0
natom=0
num=0
printf("MODEL\n")
inpfile.each do |line|
  coord=line.split(' ')
  x=coord[0];y=coord[1];z=coord[2];
  if x!=nil 
    natom=(num)%(totalnatom)
    printf("ATOM%#7d  %-2s           %10.6f%8.3f%8.3f\n",natom+1,nameatom[natom],x,y,z)
    num=num+1
    if (natom==totalnatom-1)
      printf("ENDMOD\n")
      printf("MODEL\n")
  end
  end

end
