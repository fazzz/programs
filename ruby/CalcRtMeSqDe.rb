#CalcRtMeSqDe.rb
#usage: ruby CalcRtMeSq.rb inputfile valuetype refnum
#last modified 3III10

def CalcRtMeSq(inputfile,valuetype,refnum)
  me_sq_de_value=0.0
  ref_value=0.0
  linenow=[]
  
  num=0
  inputfile.each do |line|
    linenow[num]=line.split(' ')
    if num==refnum then
      ref_value = linenow[num][valuetype].to_f
      me_sq_de_value = 0.0
      p me_sq_de_value
    elsif num > refnum then
      me_sq_de_value = ((num-refnum)*me_sq_de_value+(linenow[num][valuetype].to_f-ref_value)**2)/(num-refnum+1)
      p (me_sq_de_value**0.5)/(ref_value**2)**0.5
    end
    num = num+1
  end
  inputfile.close
end

if (argc=ARGV.length)<3 then
  print "USAGE: ./CalcRtMeSqDe.rb inputfile valuetype refnum"
  exit
end

inputfile = open(ARGV.shift)
valuetype = (ARGV.shift).to_i
refnum = (ARGV.shift).to_i
CalcRtMeSq(inputfile,valuetype,refnum)
