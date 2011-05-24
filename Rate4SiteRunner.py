'''
Run rate4site on given alignment and a tree. The tree is given as a filename,
the aln as a set of records in fasta format
'''
from subprocess import Popen,PIPE
import dendropy
from Bio import SeqIO
import sys,os,os.path
import tempfile
from Bio import AlignIO
def narrow_tree(tree,aln):
    nodes=[]
    for record in  aln:
        nodes.append(record.id)  
    nd = tree.nodes()
    for n in nd:
        if n.is_leaf() and (n.taxon is not None) and (str(n.taxon) not in nodes):
            tree2=dendropy.treemanip.prune_subtree(tree,n)
    return tree

class Rate4Site(object):
    def __init__(self,aln,treef,cmd=None):
        if os.path.isfile(aln):
            self.alnfile=aln
            self.aln = AlignIO.read(open(self.alnfile),'fasta')
        else:
            self.aln=aln
            self.alnfile = tempfile.NamedTemporaryFile()
            AlignIO.write(aln,self.alnfile,'fasta')
            self.alnfile.flush()
        if not cmd:
            import sys
            if sys.maxint==9223372036854775807: #64 bit
                cmd='rate4site64'
            else:
                cmd='rate4site'
        if isinstance(treef,dendropy.Tree):
            parent_tree=treef
        elif os.path.isfile(treef):
            parent_tree=dendropy.Tree.get_from_path(treef,'newick')
        self.tree = narrow_tree(parent_tree,self.aln)
        self.treefile = tempfile.NamedTemporaryFile()
#        self.tree.write(self.treefile,'newick',internal_labels=False)
        self.treefile.write(self.tree.as_string('newick',internal_labels=False)[5:])
        self.treefile.flush()
        self.cmd=cmd

    def runRate(self,refseq=None,outname='outfile_rate4site'):
        model=' -mn -bn -Im '
        aln = ' -s '+self.alnfile
        tree= ' -t '+self.treefile.name
        outf =' -o '+outname
        runcmd = self.cmd+model+aln+tree+outf
        if refseq:
            runcmd+=' -a '+refseq
        pipe = Popen(runcmd,shell=True, stdout=PIPE, stderr=PIPE).stdout
        rates=[]
        for line in pipe:
            if 'rate of pos:' in line:
                rates.append(float(line.strip().split()[-1]))
        return rates
                

        
        
        
    


