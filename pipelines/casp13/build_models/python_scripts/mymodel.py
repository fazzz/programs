from modeller import *
from modeller.automodel import *        # Load the automodel class
from modeller.parallel import *
class MyModel(automodel):
        def add_sec_str(self,fname):
                self.myvar_secstr = fname;
        def add_contact(self,fname):
                self.myvar_contact = fname;

        def special_restraints(self, aln):
                rsr = self.restraints
                at = self.atoms



                if hasattr(self, 'myvar_secstr'):
                        tt = "";

                        f = open(self.myvar_secstr, 'r')
                        for line in f:
                                l = re.findall('>', line)
                                if len(l) == 0:
                                        tt += line;
                        f.close()

                        tt = re.sub(r'[\r\n]', '', tt)
                        char_list = list(tt+"#")
                        istart = 1;
                        prevchar = char_list[0]
                        for i in range(1, len(char_list)):
                                if char_list[i] != prevchar :
                                        if prevchar == "H":
                                                rsr.add(secondary_structure.alpha(self.residue_range(str(istart)+':', str(i)+':')))
                                        if prevchar == "E":
                                                rsr.add(secondary_structure.strand(self.residue_range(str(istart)+':', str(i)+':')))


                                        prevchar = char_list[i];
                                        istart = i+1;


                if hasattr(self, 'myvar_contact'):
                        tt = "";

                        f = open(self.myvar_contact, 'r')
                        pat = re.compile(r"^[\s]*([0-9]+)[\s]+([0-9]+)");

                        for line in f:

                                match = pat.search(line);
                                if match is not None:
                                        rsr.add(forms.gaussian(group=physical.xy_distance,
                                                feature=features.distance(at["CA:"+match.group(1)], at["CA:"+match.group(2)]),
                                                mean=8, stdev=2))
                        f.close()

