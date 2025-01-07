# -*- coding: cp932 -*-
#CalcDihed.rb

class Peptide
  def initialize
    @coord=[]
  end

  def inputcoord(x,y,z)
    @coord << [x,y,z]
  end

  def clear
    @coord.clear
  end

  #二面角(i,j,k,l)を計算する
  def dihed(i,j,k,l)
    dist_ij=dist(i,j)
    dist_jk=dist(j,k)
    dist_kl=dist(k,l)

    angl_ijk=angle(i,j,k)
    angl_jkl=angle(j,k,l)

    v_ji=[];v_jk=[];v_kj=[];v_kl=[];
    [0,1,2].each do |alpha|
      v_ji[alpha]=((@coord[i][alpha]).to_f-(@coord[j][alpha]).to_f)/dist_ij
      v_jk[alpha]=((@coord[k][alpha]).to_f-(@coord[j][alpha]).to_f)/dist_jk
      v_kj[alpha]=-v_jk[alpha]
      v_kl[alpha]=((@coord[l][alpha]).to_f-(@coord[k][alpha]).to_f)/dist_kl  
    end

    v_jijk=[]
    v_klkj=[]
    [0,1,2].each do |alpha|
      v_jijk[alpha]=(otpd(v_ji,v_jk))[alpha]/(angl_ijk[1])
      v_klkj[alpha]=(otpd(v_kl,v_kj))[alpha]/(angl_jkl[1])
    end
    costheta=inpd(v_jijk,v_klkj)
#    p costheta
    if costheta>=1.0
      costheta=1
    end
    theta=Math.acos(costheta)

 #   p v_jijk
 #   abs=0.0
 #   [0,1,2].each  do |a|
 #     abs = abs+v_jijk[a]**2
 #   end
 #   p abs
 #   p v_klkj
    v_judge=otpd(v_jijk,v_klkj)
    abs=0.0
    [0,1,2].each  do |a|
      abs = abs+v_judge[a]**2
    end
    abs = abs**0.5
    cos_judge=inpd(v_jk,v_judge)/abs
#    p cos_judge
#    theta_judge=Math.acos(cos_judge)
    if cos_judge <=0.0 then
      return theta
    else
      theta=2.0*(Math.acos(-1.0))-theta
    end
    return theta
  end

  #結合長(i,j)を計算する
  def dist(i,j)
    distd=0.0
    [0,1,2].each do |alpha|
      distd=distd+((@coord[i][alpha]).to_f-(@coord[j][alpha]).to_f)**2
    end
    dist=distd**0.5
  end

  #結合角(i,j,k)を計算する
  def angle(i,j,k)
    dist_ij=dist(i,j)
    dist_kj=dist(j,k)

    v_ji=[]
    v_jk=[]
    [0,1,2].each do |alpha|
      v_ji[alpha]=((@coord[i][alpha]).to_f-(@coord[j][alpha]).to_f)/dist_ij
      v_jk[alpha]=((@coord[k][alpha]).to_f-(@coord[j][alpha]).to_f)/dist_kj
    end

    cos=inpd(v_ji,v_jk)
    sin=(1.0-cos**2)**0.5
#    p Math.acos(cos)*180/3.14
#    p Math.asin(sin)*180/3.14
    return [cos,sin]
  end
end

def inpd(v1,v2)
  innerproduct=0.0
  [0,1,2].each do |i|
    innerproduct=innerproduct+v1[i]*v2[i]
  end
  return innerproduct
end

def otpd(v1,v2)
  outerproduct=[]
  outerproduct[0]=v1[1]*v2[2]-v1[2]*v2[1]
  outerproduct[1]=v1[2]*v2[0]-v1[0]*v2[2]
  outerproduct[2]=v1[0]*v2[1]-v1[1]*v2[0]
  return outerproduct
end

#############################################
if (argc=ARGV.length) < 3
  print "USAGE: ./CalcDihed.rb inputfile(trj) numatom i j k l  > outputfile"
  exit
end

inputfile = open(ARGV.shift)
numatom = (ARGV.shift).to_i
atomi=(ARGV.shift).to_i
atomj=(ARGV.shift).to_i
atomk=(ARGV.shift).to_i
atoml=(ARGV.shift).to_i
#atompairfile = open(ARGV.shift)

i=0
#atompair=[]
#atompairfile.each do |line|
#  atompair[i]=line.split(/\s+/)
#  i=i+1
#end

#p atompair

#phi=psi=0.0
pi=Math.acos(-1.0)
coord=[]
thisline=[]
j=0
num=0
pep=Peptide.new
#inputfile.gets
inputfile.each do |line|
  num=num+1
  thisline = line.split(' ')
#  thisline.shift
  coord << thisline
  coord.flatten!
#  coord.slice!(/\s+/)
#  p coord.size
  if coord.size==numatom*3
    i=0
    while i < numatom
      pep.inputcoord(coord[i*3],coord[i*3+1],coord[i*3+2])
      i=i+1
    end
#    p pep
    phi=pep.dihed(atomi,atomj,atomk,atoml)  
    if (phi>pi)              
      phi=phi-2.0*pi         
    end                      
    ############################
    # phi=pep.dihed(4,6,8,14)  #
    # if (phi>pi)              #
    #   phi=phi-2.0*pi         #
    # end                      #
    # psi=pep.dihed(6,8,14,16) #
    # if (psi>pi)              #
    #   psi=psi-2.0*pi         #
    # end                      #
    ############################
#    p j
#   p pep
    #############################
    # printf("%f %f\n",phi,psi) #
    #############################
    printf("%f\n",phi)
    pep.clear
    coord.clear
#    break
  end
#    p coord
#    break
#    i=0
#    break
#     atompair.each do |pair|
#      i=pair[0];j=pair[1];k=pair[2];l=pair[3];
#     end
#     i=0
#  end
end
#############################################
