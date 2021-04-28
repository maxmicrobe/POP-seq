
def write_full_seq(gbk_file, out_path):
	''' makes a fasta file that contains 
	all the nucleotide sequence from genbank file
	'''
	from Bio import SeqIO
	rec = SeqIO.read(gbk_file,'gb')
	seq = str(rec.seq)
	header = rec.name
	with open(out_path, 'w') as f:
		f.write('>%s\n%s'%(header,seq))
def get_full_info_df(gbk_file_path, records='genome'):

    '''gets all available information pesent in 
    a genbank file depending if that information 
    is there or not
    Parameters:
    records='genome' or "plasmid"
    '''
    import pandas as pd
    from Bio import SeqIO
    from Bio.Seq import Seq

    aa=[]
    nt=[]
    start=[]
    stop=[]
    locus_tag=[]
    old_locus_tag=[]
    gene=[]
    func=[]
    strand=[]
    accessions=[]

    for rec in (SeqIO.parse(gbk_file_path,'genbank')):
        accessions.append(rec.id)
    print (accessions[0])
    for rec in (SeqIO.parse(gbk_file_path,'genbank')):
        seq = rec.seq
        if accessions[0] and records=='genome':

            for feat in rec.features:



                if feat.type =='CDS' or feat.type=='rRNA' or feat.type=='tRNA':
                    #try:    aa.append(feat.qualifiers['translation'][0])
                    #except: aa.append('na')

                    try:    start.append(feat.location.start)
                    except: start.append('na')

                    try:    stop.append(feat.location.end)
                    except: stop.append('na')

                    try:    strand.append(feat.strand)
                    except: strand.append('na')

                    try:    func.append(feat.qualifiers['product'][0])
                    except: func.append('na')

                    try:    locus_tag.append(feat.qualifiers['locus_tag'][0])
                    except: locus_tag.append('na')

                    try:    old_locus_tag.append(feat.qualifiers['old_locus_tag'][0])
                    except: old_locus_tag.append('na')

                    try:    gene.append(feat.qualifiers['gene'][0])
                    except: gene.append('na')

                    try: nt.append(str(feat.location.extract(seq)))
                    except: nt.append('na')
                    
                    dna = Seq(str(feat.location.extract(seq)))
                    try: aa.append(str(dna.translate()))
                    except: aa.append('na')


            gene_df = pd.DataFrame({'locus_tag':locus_tag,'old_locus_tag':old_locus_tag,'gene':gene,'function':func,'start':start,'stop':stop,
                                    'strand':strand,'nts':nt,'aa':aa},
                                columns=['locus_tag','old_locus_tag','gene','function','start','stop','strand','nts','aa'])
            gene_df['gene_length']=gene_df['stop']-gene_df['start']

            return gene_df

        elif records== "plasmid" and rec.id == accessions[1]:
            aa=[]
            nt=[]
            start=[]
            stop=[]
            locus_tag=[]
            old_locus_tag=[]
            gene=[]
            func=[]
            strand=[]
            for feat in rec.features:



                if feat.type =='CDS' or feat.type=='rRNA' or feat.type=='tRNA':
                    #try:    aa.append(feat.qualifiers['translation'][0])
                    #except: aa.append('na')

                    try:    start.append(feat.location.start)
                    except: start.append('na')

                    try:    stop.append(feat.location.end)
                    except: stop.append('na')

                    try:    strand.append(feat.strand)
                    except: strand.append('na')

                    try:    func.append(feat.qualifiers['product'][0])
                    except: func.append('na')

                    try:    locus_tag.append(feat.qualifiers['locus_tag'][0])
                    except: locus_tag.append('na')

                    try:    old_locus_tag.append(feat.qualifiers['old_locus_tag'][0])
                    except: old_locus_tag.append('na')

                    try:    gene.append(feat.qualifiers['gene'][0])
                    except: gene.append('na')

                    try: nt.append(str(feat.location.extract(seq)))
                    except: nt.append('na')
                    
                    dna = Seq(str(feat.location.extract(seq)))
                    try: aa.append(str(dna.translate()))
                    except: aa.append('na')
            plasmid_df = pd.DataFrame({'locus_tag':locus_tag,'old_locus_tag':old_locus_tag,'gene':gene,'function':func,'start':start,'stop':stop,
                                    'strand':strand,'nts':nt,'aa':aa},
                                columns=['locus_tag','old_locus_tag','gene','function','start','stop','strand','nts','aa'])
            plasmid_df['gene_length']=plasmid_df['stop']-plasmid_df['start']
            print len(plasmid_df)
            return plasmid_df
        
        
        elif records== "plasmid2" and rec.id == accessions[2]:
            aa=[]
            nt=[]
            start=[]
            stop=[]
            locus_tag=[]
            old_locus_tag=[]
            gene=[]
            func=[]
            strand=[]
            for feat in rec.features:



                if feat.type =='CDS' or feat.type=='rRNA' or feat.type=='tRNA':
                    #try:    aa.append(feat.qualifiers['translation'][0])
                    #except: aa.append('na')

                    try:    start.append(feat.location.start)
                    except: start.append('na')

                    try:    stop.append(feat.location.end)
                    except: stop.append('na')

                    try:    strand.append(feat.strand)
                    except: strand.append('na')

                    try:    func.append(feat.qualifiers['product'][0])
                    except: func.append('na')

                    try:    locus_tag.append(feat.qualifiers['locus_tag'][0])
                    except: locus_tag.append('na')

                    try:    old_locus_tag.append(feat.qualifiers['old_locus_tag'][0])
                    except: old_locus_tag.append('na')

                    try:    gene.append(feat.qualifiers['gene'][0])
                    except: gene.append('na')

                    try: nt.append(str(feat.location.extract(seq)))
                    except: nt.append('na')
                    
                    dna = Seq(str(feat.location.extract(seq)))
                    try: aa.append(str(dna.translate()))
                    except: aa.append('na')


            plasmid_df2 = pd.DataFrame({'locus_tag':locus_tag,'old_locus_tag':old_locus_tag,'gene':gene,'function':func,'start':start,'stop':stop,
                                    'strand':strand,'nts':nt,'aa':aa},
                                columns=['locus_tag','old_locus_tag','gene','function','start','stop','strand','nts','aa'])
            plasmid_df2['gene_length']=plasmid_df2['stop']-plasmid_df2['start']
            print len(plasmid_df2)
            return plasmid_df2

    

def write_fna_file(gbk_file,out_path):
	'''writes a fasta file that has all the nucleotide sequences of CDS, rRNA and tRNA'''
	from Bio import SeqIO

	nt=[]
	locus_tag=[]


	rec = SeqIO.read(gbk_file,'gb')
	seq = rec.seq

	for feat in rec.features:
		if feat.type =='CDS' or feat.type=='rRNA' or feat.type=='tRNA':
			try: nt.append(str(feat.location.extract(seq)))
			except: nt.append('na')
			try:    locus_tag.append(feat.qualifiers['locus_tag'][0])
			except: locus_tag.append('na')

	with open(out_path,'w') as f:
		for i,j in zip(locus_tag,nt):
				f.write('>%s\n%s\n'%(i,j))



def write_faa_file(gbk_file,out_path):
	from Bio import SeqIO

	cds=[]
	locus_tag=[]


	rec = SeqIO.read(gbk_file,'gb')
	seq = rec.seq

	for feat in rec.features:
		if feat.type =='CDS':
			try: cds.append(feat.qualifiers['translation'][0])
			except: cds.append('na')
			try:    locus_tag.append(feat.qualifiers['locus_tag'][0])
			except: locus_tag.append('na')

	with open(out_path,'w') as f:
		for i,j in zip(locus_tag,cds):
				f.write('>%s\n%s\n'%(i,j))

def method_list():
	func_list=['write_full_seq', 'get_full_info_df', 'write_fna_file', 'write_faa_file', 'method_list',
               'make_saf_file', 'make_promoter_saf_file','AT_content_wig','GC_content_wig', 
               'get_codon_adaptation_index','get_rbs_sites', 'get_upstream_full_info_df', 'make_saf_file_defined_bits']
	for i in func_list:
		print '%s\n'%i


def make_saf_file(gbk_infile,saf_outfile):
    '''takes gbk file and makes saf file. User to specify the path to both gbk and saf
the saf file is required for FeatureCounts. The saf file only has CDS [no t/rRNA]'''
   
    import pandas as pd
    from Bio import SeqIO
    infile = SeqIO.read(gbk_infile,'gb')
    genes =[]
    strand = []
    start = []
    stop = []
    accession=infile.annotations['accessions'][0]

    for feature in infile.features:
        if feature.type == 'CDS':  #Only protein-coding genes
            genes.append(feature.qualifiers['locus_tag'][0])
            if feature.strand == 1:
                strand.append("+")
                start.append(feature.location.start.real+1)
                stop.append(feature.location.end.real)
            elif feature.strand == -1:
                strand.append("-")
                start.append(feature.location.start.real+1)
                stop.append(feature.location.end.real)
    gene_df = pd.DataFrame({"GeneID": genes,  'Chr':accession ,"Start": start, "Stop": stop,'Strand': strand},
                              columns = ['GeneID', 'Chr', 'Start','Stop', 'Strand'])
    gene_df.to_csv(saf_outfile, sep='\t', index=False)


def make_promoter_saf_file(gbk_infile,saf_outfile,upstream_length=100):
    '''takes gbk file and makes saf file of promoter region upstream of each gene. 
    User to specify the path to both gbk and saf
the saf file is required for FeatureCounts. The saf file only has CDS [no t/rRNA]'''
   
    import pandas as pd
    from Bio import SeqIO
    genes =[]
    strand = []
    start = []
    stop = []
    accession=[]
    for infile in SeqIO.parse(open(gbk_infile,"r"), "gb") :

        for feature in infile.features:
            if feature.type == 'CDS':  #Only protein-coding genes
                accession.append(infile.id)
                genes.append(feature.qualifiers['locus_tag'][0])


                if feature.strand == 1:
                    strand.append("+")
                    stop.append(feature.location.start.real+1)
                    start.append(feature.location.start.real-99)
                elif feature.strand == -1:
                    strand.append("-")
                    stop.append(feature.location.start.real+1)
                    start.append(feature.location.start.real-99)
    gene_df = pd.DataFrame({"GeneID": genes,  'Chr':accession ,"Start": start, "Stop": stop,'Strand': strand},
                           columns = ['GeneID', 'Chr', 'Start','Stop', 'Strand'])
    gene_df.to_csv(saf_outfile, sep='\t', index=False)



def make_saf_file_defined_bits(gbk_infile, saf_outfile, bitsize=15):
    '''takes gbk file and makes saf file but in defined lengths. User to specify the length of each bit, and the path to both gbk and saf
the saf file is required for FeatureCounts'''

    import pandas as pd
    from Bio import SeqIO
    infile = SeqIO.read(gbk_infile, 'gb')
    features = []
    strand = []
    start = []
    stop = []
    accession = infile.annotations['accessions'][0]
    sequence = str(infile.seq)

    for i in list(range(1, len(sequence), bitsize)):
        st = i
        sp = i+ bitsize
        feature= "%s_%s"%(st,sp)
        dna = '+'
        
        features.append(feature)
        strand.append(dna)
        start.append(st)
        stop.append(sp)      

    
    gene_df = pd.DataFrame({"GeneID": features,  'Chr': accession, "Start": start, "Stop": stop, 'Strand': strand},
                           columns=['GeneID', 'Chr', 'Start', 'Stop', 'Strand'])
    
    gene_df.to_csv(saf_outfile, sep='\t', index=False)


def GC_content_wig(fasta_path, wig_output, smooth=10001):
    '''calculates smoothed GC content of the whole genome. The smoothing has to be 1001, 10001 etc. 
    the wig output should be in a "path/file.wig" format at the end'''
    from scipy import signal
    with open(fasta_path, 'r') as infile:
        x = infile.read().splitlines()
        fil = x[1:2]
	accession = x[0][1:]

    gc_cont=[]
    window = 100
    for i in range(len(fil[0])):
        ss = fil[0][i:window+i+1]
        gg =(ss.count('G'))
        cc =(ss.count('C'))
        gc =gg+cc/float(window)
        gc_cont.append(gc)


    gc_cont = signal.savgol_filter(gc_cont, smooth, 3)

    header1 = 'track\ttype=wiggle_0\tname=GC_content\tgraphType=points\tvisibility=full\tcolor=168,130,88\n'
    header2='fixedStep\tchrom=%s\tstart=1\tstep=1\tspan=1\n'%accession
    with open (wig_output, 'w+') as f:
        f.write(header1)
        f.write(header2)
        for i in gc_cont:
            f.write('%0.1f\n'%i)


def AT_content_wig(fasta_path, wig_output, smooth=10001):
    '''calculates smoothed GC content of the whole genome. The smoothing has to be 1001, 10001 etc. 
    the wig output should be in a "path/file.wig" format at the end'''
    from scipy import signal
    with open(fasta_path, 'r') as infile:
        x = infile.read().splitlines()
        fil = x[1:2]
	accession = x[0][1:]

    at_cont=[]
    window = 100
    for i in range(len(fil[0])):
        ss = fil[0][i:window+i+1]
        tt =(ss.count('T'))
        aa =(ss.count('A'))
        at =tt+aa/float(window)
        at_cont.append(at)


    at_cont = signal.savgol_filter(at_cont, smooth, 3)

    header1 = 'track\ttype=wiggle_0\tname=AT_content\tgraphType=points\tvisibility=full\tcolor=168,130,88\n'
    header2='fixedStep\tchrom=%s\tstart=1\tstep=1\tspan=1\n'%accession
    with open (wig_output, 'w+') as f:
        f.write(header1)
        f.write(header2)
        for i in at_cont:
            f.write('%0.1f\n'%i)


def get_codon_adaptation_index(fasta): 
    '''returns a dataframe that calculates codon adaptation for gene-DNA fasta file'''
    index=[]
    gene_ids=[]
    import pandas as pd
    from Bio import SeqIO
    from Bio.SeqUtils import CodonUsage
    from Bio.SeqUtils.CodonUsage import CodonAdaptationIndex
    for record in SeqIO.parse(fasta, "fasta"):
        sequence=str(record.seq)
        
        try:
            Seq_Cai = CodonUsage.CodonAdaptationIndex()
            cai=Seq_Cai.cai_for_gene(sequence)
            index.append(cai)
            gene_ids.append(record.id)
        except:
            index.append('problem')
            gene_ids.append(record.id)
    df = pd.DataFrame({'gene_id':gene_ids, 'cai':index}, columns=['gene_id','cai'])
    return df


def get_rbs_sites(gbk, temp_output,antiSD_seq='ACCTCCTT' ):
    '''
    returns a dataframe that contains the delta G values and the ribosome bindng sites. 
    The input is a genbank file
    The anti Shine Dalgarno sequences can be set according to the organism
    '''
    from Bio import SeqIO
    from Bio import Seq
    import pandas as pd
    import numpy as np
    import subprocess
    
    infile = next(SeqIO.parse(gbk,'gb')) # I use next in order to avoid errors from multiple record genbank files.
    
    genes =[]
    name = []
    product = []
    func = []
    strand = []
    start = []
    stop = []
    aaseq = []
    cds_seq = []
    real_st=[]
    real_sp=[]

    SD_start=[]
    SD_end=[]
    SD_dg=[]
    SD_seq=[]
    genes2=[]

        #genome_seq_df = pd.DataFrame({'sequence':list(str(infile.seq))},index=range(1,len(str(infile.seq))+1))
    for feature in infile.features:
        if feature.type == 'CDS' and 'product' in feature.qualifiers:  #Only cares for coding sequences which are not pseudogenes
            genes.append(feature.qualifiers['locus_tag'][0])
            try: name.append(feature.qualifiers['gene'][0])
            except: name.append('')
            product.append(feature.qualifiers['product'][0])
            #cds_seq.append(feature.location.extract(infile.seq).tostring())
            if 'function' in feature.qualifiers:                       #not all genes have known functions
                func.append(feature.qualifiers['function'][0])
            else:
                func.append("N/A")
            try:
                aaseq.append(feature.qualifiers['translation'][0])
            except:
                aaseq.append("N/A")
            if feature.strand == 1:
                strand.append("plus")
                real_st.append(feature.location.start.real+1)
                real_sp.append(feature.location.end.real)
                st=feature.location.start.real-25
                sp=feature.location.start.real
                start.append(st)
                stop.append(sp)
                cds_seq.append(str(infile.seq[st:sp]))
            elif feature.strand == -1:
                real_st.append(feature.location.start.real+1)
                real_sp.append(feature.location.end.real)
                strand.append("minus")
                st=feature.location.end.real
                sp=feature.location.end.real+25
                start.append(st)
                stop.append(sp)
                my_seq = str(infile.seq[st:sp])
                my_seq = Seq.reverse_complement(my_seq)
                cds_seq.append(my_seq)
    gene_df = pd.DataFrame({"gene_id": genes,  "product": product, "strand": strand, "start": start, "stop": stop, 
                            "cds_seq":cds_seq,"aaseq": aaseq, "real_st" :real_st, "real_sp": real_sp},
                              columns = ["gene_id",  "product", "strand", "start", "stop", "cds_seq","aaseq",'real_st', 'real_sp'])




    for index, row in gene_df.iterrows():

        seq = row['cds_seq']
        fromm=row['start']
        gene_id=row['gene_id']
        scores=[]
        position_from=[]
        position_to=[]

        #for i in range(0,(len(seq)+1)-10): # windows of 10 nts step 1 to cover all seq
        for i in range(0,16):
            newseq= seq[i:i+10]
            fromm=fromm + 1
            f=open('%s/temp'%temp_output, 'w') #write each window seq
            writing= ">{}\n{}&{}\n".format(fromm,antiSD_seq,newseq)
            f.write(writing)
            f.close()


            #will have to install The ViennaRNA Package. RNAsubopt is part of the package (version 2.3.1)
            '''RNAsubopt reads RNA sequences from stdin and (in the default -e mode) calculates all
            suboptimal secondary structures within a user defined energy range above the
            minimum free energy (mfe). It prints the suboptimal structures in dot-bracket
            notation followed by the energy in kcal/mol to stdout. Be careful, the number
            of structures returned grows exponentially with both sequence length and energy
            range.
            '''
            command = 'RNAsubopt -e 1 -s --noconv < /home/m/Desktop/temp'#version 2.3.1
            output=subprocess.check_output(command, shell = True)
                ####dealing with the RNAsubopt results####
            array = output.split('\n')

            scores.append(array[1]) #element 1 contains the scores
            coordinates=(array[0][1:-6])
            position_from.append(coordinates)
            position_to.append(int(coordinates)+9)


        dg=[]
        seq_ar=[]

        for j in scores: #this calculates the minimum delta G and the element position of the minimum deltaG
            deltaG=j.split()[1]
            seq1 = j.split()[0]
            seq_ar.append(seq1)
            dg.append(np.float64(deltaG))
        final_dg=np.min(dg) #this calculates the minimum delta G and the element position of the minimum deltaG
        final_element=np.argmin(dg) # the element position of the minimum deltaG
        final_from=position_from[final_element]
        final_to = position_to[final_element]
        seq2=seq_ar[final_element]
        seq3=seq2[9:]
        SD_start.append(final_from)
        SD_end.append(final_to)
        SD_dg.append(final_dg)
        SD_seq.append(seq3)
        genes2.append(gene_id)

    RBS_df =pd.DataFrame({'gene_id':genes2,
                         'RBS_start':SD_start,
                         'RBS_end':SD_end,
                         'RBS_dg':SD_dg,
                         'RBS_seq':SD_seq},
                        columns=['gene_id','RBS_start', 'RBS_end', 'RBS_dg', 'RBS_seq'])
    return RBS_df


def get_upstream_full_info_df(gbk_file_path, upstream=25):

    '''gets all available information pesent in

    a genbank file depending if that information

    is there or not

    '''

    import pandas as pd

    from Bio import SeqIO

    from Bio import Seq

 

    aa=[]

    nt=[]

    start=[]

    stop=[]

    locus_tag=[]

    old_locus_tag=[]

    gene=[]

    func=[]

    strand=[]

    upstream_seq=[]

 

    rec = next(SeqIO.parse(gbk_file_path,'gb'))

    seq = rec.seq

 

    for feat in rec.features:

        if feat.type =='CDS' or feat.type=='rRNA' or feat.type=='tRNA':

            try:    aa.append(feat.qualifiers['translation'][0])

            except: aa.append('na')

 

            try:    start.append(feat.location.start)

            except: start.append('na')

 

            try:    stop.append(feat.location.end)

            except: stop.append('na')

 

            try:    strand.append(feat.strand)

            except: strand.append('na')

 

            try:    func.append(feat.qualifiers['product'][0])

            except: func.append('na')

 

            try:    locus_tag.append(feat.qualifiers['locus_tag'][0])

            except: locus_tag.append('na')

 

            try:    old_locus_tag.append(feat.qualifiers['old_locus_tag'][0])

            except: old_locus_tag.append('na')

 

            try:    gene.append(feat.qualifiers['gene'][0])

            except: gene.append('na')

 

            try: nt.append(str(feat.location.extract(seq)))

            except: nt.append('na')

           

    

    

    ############################### This part takes the upstream sequences ##########################

       

            if feat.strand == 1:

                try:

                    st=feat.location.start.real - upstream

                    sp=feat.location.start.real

                    upstream_seq.append(str(seq[st:sp]))

                except:

                    upstream_seq.append('na')

               

            elif feat.strand == -1:

               

                try:

                    st=feat.location.end.real

                    sp=feat.location.end.real + upstream

                    my_seq = seq[st:sp]

                    my_seq1 = Seq.reverse_complement(my_seq)

                    upstream_seq.append(str(my_seq1))

                except:

                    upstream_seq.append('na')

 

    

    ################################ This part saves and returns dataframe #############################

           

            

    gene_df = pd.DataFrame({'locus_tag':locus_tag,'old_locus_tag':old_locus_tag,'gene':gene,'function':func,'start':start,'stop':stop,

                            'strand':strand,'nts':nt,'aa':aa,'upstream_seq':upstream_seq},

                        columns=['locus_tag','old_locus_tag','gene','function','start','stop','strand','nts','aa','upstream_seq'])

    gene_df['gene_length']=gene_df['stop']-gene_df['start']

    return gene_df
