#!/home/tsukasa/.pyenv/shims/python
# -*- coding: utf-8 -*-
from modeller import *
from modeller.automodel import *	# Load the automodel class
from modeller.parallel import *

import shutil


import re
import os
import random
import sys
import shutil
import fcntl
import datetime
import locale

from mymodel import MyModel
#class MyModel(automodel):
#	def add_sec_str(self,fname):
#		self.myvar_secstr = fname;
#	def add_contact(self,fname):
#		self.myvar_contact = fname;
#	
#	def special_restraints(self, aln):
#		rsr = self.restraints
#		at = self.atoms
#		
#		
#		
#		if hasattr(self, 'myvar_secstr'):
#			tt = "";
#			
#			f = open(self.myvar_secstr, 'r')
#			for line in f:
#				l = re.findall('>', line)
#				if len(l) == 0:
#					tt += line;
#			f.close()
#
#			tt = re.sub(r'[\r\n]', '', tt)
#			char_list = list(tt+"#")
#			istart = 1;
#			prevchar = char_list[0]
#			for i in range(1, len(char_list)):
#				if char_list[i] != prevchar :
#					if prevchar == "H":
#						rsr.add(secondary_structure.alpha(self.residue_range(str(istart)+':', str(i)+':')))
#					if prevchar == "E":
#						rsr.add(secondary_structure.strand(self.residue_range(str(istart)+':', str(i)+':')))
#					
#					
#					prevchar = char_list[i];
#					istart = i+1;
#		
#		
#		if hasattr(self, 'myvar_contact'):
#			tt = "";
#			
#			f = open(self.myvar_contact, 'r')
#			pat = re.compile(r"^[\s]*([0-9]+)[\s]+([0-9]+)");
#			
#			for line in f:
#				
#				match = pat.search(line);
#				if match is not None:
#					rsr.add(forms.gaussian(group=physical.xy_distance,
#						feature=features.distance(at["CA:"+match.group(1)], at["CA:"+match.group(2)]),
#						mean=8, stdev=2))
#			f.close()
		
		







currentdate = datetime.datetime.today()
datecode = currentdate.strftime("%m%d%H%M");

## modeller_run.py pir_file pdb_directory tmp_directory secondarystructure_file(optional) contact_file(optional)
## ss_file should be like "HHHHEEEEEECCCCCEEEEECCCCEEEE"  and  ">" line like fasta can be permitted
## contact_file residue_number<tab/space>residue_number<tab/space>restraint


fcount = 0;
flimit = 1;

file = os.path.abspath(sys.argv[1]);
# modeldir = os.path.abspath(sys.argv[2]);
outdir = sys.argv[2];
strfile = "";
contactfile = "";
if len(sys.argv) > 4:
	strfile = sys.argv[3];

if len(sys.argv) > 5:
	contactfile = sys.argv[4];

aln= file;
cname = 'test';
seqname=[]

log.verbose();

if not os.path.exists(outdir):
	os.mkdir(outdir)

os.chdir(outdir);		   #change current directory
#shutil.copyfile(file,outdir)


env = environ(rand_seed=int(os.getpid())+int(random.random()*10000))
aaln = alignment(env)

for line in open(aln, 'r'):
	tm = re.search('^sequence\:([^\:]+)\:', line)
	if tm:
		cname = tm.group(1);

dirlist = [];
# dirlist.append(modeldir);
#dirlist.append("./worktmp2/");
dirlist.append("/work2/tsukasa/pdb/");
dirlist.append("/work2/tsukasa/pdb_obsolete");
env.io.atom_files_directory = dirlist

#print(env.io.atom_files_directory)

for line in open(aln, 'r'):
	tm = re.search('>[^;]+;[\s]*([^\s]+)', line)
	if tm:
		ss =  tm.group(1)
		if ss != cname:
			seqname.append(ss);
	
	tm = re.search('structure[XNM]?\:([^\:]+)', line)
	if tm:
		ss =  tm.group(1)
		mdl = model(env,file=ss, model_format="PDB_OR_MMCIF")
		aaln.append_model(mdl,align_codes=ss)


#aaln.write(file='debug.seq')

if len(seqname) == 1:
	a = MyModel(env, alnfile=aln,knowns=seqname, sequence=cname,assess_methods=(assess.DOPE, assess.GA341))
else:
	a = MyModel(env, alnfile=aln,knowns=seqname, sequence=cname,assess_methods=(assess.DOPE, assess.GA341))




if len(strfile) > 0:
	a.add_sec_str(strfile);


if len(contactfile) > 0:
	a.add_contact(contactfile);




a.starting_model = 1
a.ending_model = 5
#a.md_level = refine.slow
#j=job()
#for _ in range(5):
#    j.append(local_slave())
#a.use_parallel_job(j)
a.make()
