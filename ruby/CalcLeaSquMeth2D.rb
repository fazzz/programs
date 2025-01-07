#CalcLeaSquMeth2D.rb
#ruby CalcRamChad.rb DihedData.txt 0.01

class Profile
  def initialize
    @data=[]
  end

  def inputdata(x,y)
    @data << [x,y]
  end

  def calcleasqumeth2d()
    sumx=0.0;sumy=0.0;sumx2=0.0;sumxy=0.0;
    n=0
    @data.each do |xy|
      sumx  = sumx+xy[0].to_f
      sumy  = sumy+xy[1].to_f
      sumx2 = sumx2+(xy[0].to_f)**2
      sumxy = sumxy+(xy[0].to_f)*(xy[1].to_f)
      n = n+1
    end
    f=sumx2*n-(sumx)**2
    a=(n*sumxy-sumx*sumy)/f
    b=(-sumx*sumxy+sumx2*sumy)/f
    p a
    p b
  end
end

inputfile = open(ARGV.shift)

thisline=[]
prof=Profile.new
inputfile.each do |line|
  thisline = line.split(' ')
  prof.inputdata(thisline[0],thisline[1])
end

#p prof

prof.calcleasqumeth2d()


