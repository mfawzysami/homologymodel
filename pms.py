from Bio.Blast import NCBIWWW, NCBIXML
import os, sys, argparse, StringIO
from Bio.SeqIO.PirIO import PirWriter
import Bio.SeqIO
import urllib2
from modeller import *
from Bio.Alphabet import generic_protein
from modeller.automodel import *



class ProteinModeling(object):
    def __init__(self):
        self.db_pir_name = "db.pir"
        self.db_bin_name = "db"

    def parse(self):
        self.parser = argparse.ArgumentParser(description="Protein Homology Modelling Service.")
        self.parser.add_argument("-i", "--input", help="Query protein Fasta file")
        self.parser.add_argument("-o", "--output", default=os.curdir,
                                 help="Output Directory where all downloaded files are stored alongside "
                                      "predicted query protein model(s). Defaults to current working directory")
        self.parser.add_argument("-r", "--results",
                                 help="Blast results XML file if exists for the query sequence. if not provided, the software will perform blast search against the query protein sequence")
        self.parser.add_argument("-d", "--nodownloads", nargs='?', const=True, default=False,
                                 help="If specified PMS will not download PDB files but it will search the output directory for these files during downstream analysis")
        self.parser.add_argument("-s","--min",default=30,help="Sequences below this length will be ignored during the modelling process. Defaults : 30 residues",type=int)
        self.parser.add_argument("-m","--max",default=4000,help="Sequences above this length will be ignored during the modeling process. Defaults: 4000 residues",type=int)
        self.parser.add_argument("-c","--clean",default=True,help="Clean sequences from all non-standard protein residues. Defaults: True",nargs="?",const=True)
        self.parser.add_argument("-t","--iterations",type=int,help="Number of Search iterations",default=1)
        self.parser.add_argument('-k',"--check",default=False,help="Check the modeling profile for deviations . Defaults: False",nargs="?",const=True)
        self.parser.add_argument('-e',"--evalue",default=0.01,type=float,help="Include sequences whose E-value is larger or equal to this provided value. Defaults: 0.01")
        self.parser.add_argument('-n','--count',type=int,help="Number of blast search results to use for model building.",default=100)
        self.parser.add_argument('-b','--build',const=True,nargs="?",default=False,help="Instruct PMS to build a homology model for a specific PDB")
        self.parser.add_argument("-u","--take",default=None,type=str,help="The specific PDB ID to use to create a homology model for the unknown query protein")
        self.parser.add_argument("-z","--models",default=5,type=int,help="Number of predicted models to generate for the unknown protein")
        self.args = self.parser.parse_args()
        if len(sys.argv) <= 1:
            self.parser.print_help()
            return
        self.predict()

    def convert_to_pir(self,fasta):
        buff = StringIO.StringIO("")
        fasta_fp = Bio.SeqIO.parse(fasta,'fasta',alphabet=generic_protein)
        buff = StringIO.StringIO("")
        pir_value = ""
        for record in fasta_fp:
            pdb_id = record.name.split(':')[0]
            buff.write(">P1;{0}\r".format(pdb_id))
            buff.write("sequence:{0}:::::::0.00: 0.00\r".format(pdb_id))
            buff.write(record.seq.upper()+"*\r")
        pir_value = buff.getvalue()
        return pir_value




    def predict(self):
        # step perform blast search
        try:
            if not os.path.exists(self.args.output):
                os.mkdir(self.args.output)
            self.read_input()
            blast_results = self.perform_blast_search()
            self.create_modeller_db(blast_results)
            self.model(blast_results)
            self.align_and_build()
            self.make_model()
        except Exception as e:
            print(e.message)
            #self.parser.print_help()

    def read_input(self):
        if not self.args.input:
            print("You have to provide query protein sequence in Fasta Format.")

            raise Exception()
        self.input_fasta = Bio.SeqIO.read(open(self.args.input, 'r'), "fasta")
        pir_value = self.convert_to_pir(self.args.input)
        if pir_value:
            base_file, _ = os.path.splitext(os.path.basename(self.args.input))
            pir_output = open(os.path.join(self.args.output,"{0}.pir".format(base_file)),'w')
            pir_output.write(pir_value)
            pir_output.close()

        #convert the fasta into PIR format for later usage
        # pir_file = StringIO.StringIO("")
        # n = Bio.SeqIO.convert(self.args.input,"fasta",pir_file,"pir",alphabet=generic_protein)
        # if n > 0:
        #     base_file , _ = os.path.splitext(os.path.basename(self.args.input))
        #     pir_output = open(os.path.join(self.args.output,"{0}.pir".format(base_file)),'w')
        #     pir_output.write(pir_file.getvalue())
        #     pir_output.close()
        #     pir_file.close()



    def perform_blast_search(self):
        if self.args.results:
            blast_results = NCBIXML.read(open(self.args.results, 'r'))
        else:
            print("Performing Blast Search. Please Wait.....")
            results = NCBIWWW.qblast("blastp", "pdb", self.input_fasta.seq)
            xml_contents = results.read()
            results.close()
            with open(os.path.join(self.args.output, "blast_results.xml"), 'w') as output:
                output.write(xml_contents)
            blast_results = NCBIXML.read(StringIO.StringIO(xml_contents))
            if self.args.count > len(blast_results.alignments):
                self.args.count = len(blast_results.alignments)
        if not self.args.nodownloads:
            self.download_pdbs(blast_results)

        return blast_results

    def download_pdbs(self, blast_results):
        download_url = "https://files.rcsb.org/download/{0}.pdb"
        print("Downloading PDBs. Please Wait...")
        if not blast_results:
            print("There are not blast results. Please try again. Aborting....")
            raise Exception()
        for hit in blast_results.alignments[:self.args.count]:
            splitted_tags = hit.hit_id.split('|')
            pdb_id = splitted_tags[1]
            print("Downloading: {0}".format(pdb_id))
            response = urllib2.urlopen(download_url.format(pdb_id))
            output_file = os.path.join(self.args.output, "{0}.pdb".format(pdb_id))
            file = open(output_file, 'w')
            file.write(response.read())
            file.close()
            response.close()

    def create_modeller_db(self, blast_results):
        print("Creating MODELLER database. Please Wait....")
        db_loc = os.path.join(self.args.output,"db.fa")
        db = open(db_loc,'w')

        for hit in blast_results.alignments[:self.args.count]:
            for hsp in hit.hsps:
                db.write(">{0}\r".format(hit.title))
                db.write(hsp.sbjct)
                db.write('\r')
        db.close()
        pir_value = self.convert_to_pir(db_loc)
        if pir_value:
            with open(os.path.join(self.args.output, self.db_pir_name), 'w') as output:
                 output.write(pir_value)
        env = environ()
        sdb = sequence_db(env)
        sdb.read(seq_database_file=os.path.join(self.args.output,self.db_pir_name),seq_database_format="PIR",chains_list='ALL',minmax_db_seq_len=(self.args.min,self.args.max),
                 clean_sequences=self.args.clean)
        sdb.write(seq_database_file=os.path.join(self.args.output,"{0}.bin".format(self.db_bin_name)),seq_database_format="BINARY",chains_list='ALL')
        sdb.read(seq_database_file=os.path.join(self.args.output,"{0}.bin".format(self.db_bin_name)),seq_database_format="BINARY",chains_list='ALL')
        aln = alignment(env)
        base_name , _ = os.path.splitext(os.path.basename(self.args.input))
        aln.append(file=os.path.join(self.args.output,"{0}.pir".format(base_name)),alignment_format='PIR',align_codes="ALL")
        prf = aln.to_profile()
        prf.build(sdb, matrix_offset=-450, rr_file=os.path.join(os.curdir,'blosum62.sim.mat'),
                  gap_penalties_1d=(-500, -50), n_prof_iterations=self.args.iterations,
                  check_profile=self.args.check, max_aln_evalue=self.args.evalue,score_statistics=False)
        prf.write(file=os.path.join(self.args.output,"build_profile.prf"),profile_format='TEXT')
        aln = prf.to_alignment()
        aln.write(file=os.path.join(self.args.output,"build_profile.pir"),alignment_format="PIR")

    def model(self,blast_results):
        print("Building Protein models. Please Wait...")
        env = environ()
        aln = alignment(env)
        for hit in blast_results.alignments[:self.args.count]:
            id = hit.hit_id.split('|')[1]
            m = model(env,file=open(os.path.join(self.args.output,"{0}.pdb".format(id)),'r'))
            aln.append_model(m,atom_files=str(id),align_codes=str(id))

        aln.malign(rr_file=os.path.join(os.curdir,"as1.sim.mat"))
        aln.malign3d()
        aln.compare_structures()
        aln.id_table(matrix_file=os.path.join(self.args.output,"family.mat"))
        env.dendrogram(matrix_file=os.path.join(self.args.output,"family.mat"),cluster_cut=1.0)
        print ("Protein models have been built successfully.")

    def align_and_build(self):
        if self.args.build and not self.args.take:
            raise Exception("Please specify the best PDB model to use to create a homology model, Give a PDB ID.")
        if not self.args.build:
            return
        env = environ()
        aln = alignment(env)
        mdl = model(env, file=open(os.path.join(self.args.output,self.args.take+".pdb"),'r'))
        aln.append_model(mdl, align_codes=str(self.args.take), atom_files=str("{0}.pdb".format(self.args.take)))
        base_name, _ = os.path.splitext(os.path.basename(self.args.input))
        aln.append(file=os.path.join(self.args.output,"{0}.pir".format(base_name)), align_codes=base_name)
        aln.align2d()
        aligned_file = os.path.join(self.args.output,"{0}-{1}".format(base_name,self.args.take))
        aln.write(file=str(aligned_file+".pir"), alignment_format='PIR')
        aln.write(file=str(aligned_file+".pap"), alignment_format='PAP')

    def make_model(self):
        if self.args.build and not self.args.take:
            raise Exception("Please specify the best PDB model to use to create a homology model, Give a PDB ID.")
        if not self.args.build:
            return
        base_name, _ = os.path.splitext(os.path.basename(self.args.input))
        aligned_file = os.path.join(self.args.output, "{0}-{1}".format(base_name, self.args.take))
        env = environ()
        a = automodel(env, alnfile=str(aligned_file+".pir"),
                      knowns=self.args.take, sequence=str(base_name),
                      assess_methods=(assess.DOPE,
                                      # soap_protein_od.Scorer(),
                                      assess.GA341))
        a.starting_model = 1
        a.ending_model = self.args.models
        a.make()



def main():
    modeller = ProteinModeling()
    modeller.parse()


if __name__ == '__main__':
    main()
