'''
Created on 05.06.2014

@author: david
'''


class Feature:
    '''
    handle one feature (line) of a GTF file
    '''
    
    def __init__(self):
        self.chr = "."
        self.source = "."
        self.featureType = "."
        self.frame = "."
        self.start = 0
        self.end = 0
        self.score = "."
        self.strand = True
        self.frame = "."
        self.geneId = None
        self.transcriptId =  None
        self.attributes = {}
    
    
    def readline(self,line):
        '''
        process one line of the gtf file
        <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
        :param line: one line of a gtf file
        '''

        line = line.rstrip().split("\t")
        
        try:
            self.chr = line[0]
            self.source = line[1]
            self.featureType = line[2]
            self.start = int(line[3])
            self.end = int(line[4])
            self.score = line[5]
            self.strand = line[6] in ['1','+'] # 
            self.frame = line[7]
        except ValueError:
            raise ValueError("Error in line '%s'" % " ".join(line))
        
        self.start=int(self.start)
        self.end=int(self.end)
        
        #parse attributes
        attributes = line[8]
        #trim comments
        attributes=attributes[:attributes.find("#")].rstrip()
        values = map(lambda x: x.strip(), attributes.split(";")[:-1])
        
        for info in values:
            info = map( lambda x: x.strip(), info.split(" "))
            name, value=info[0], info[1].replace("\"","")

            if name == "gene_id":
                self.geneId = value
            elif name == "transcript_id":
                self.transcriptId = value
            elif name == "gene_biotype":
                self.source=value
            elif name == "gene_name":
                self.attributes[name]=value    
            else:
                try:
                    value=float(value)
                    value=int(value)
                except ValueError:
                    pass
                except TypeError:
                    pass    
                self.attributes[name]=value

        if not self.geneId:
            raise TypeError("no gene_id found in line %s" % " ".join(line))
        if not self.transcriptId and self.featureType != "gene":  #added to support GRCH38
            raise TypeError("no transcript_id found in line %s" % " ".join(line))

         
def iterator(infile):

    while 1:
        line = infile.readline()
        if not line: raise StopIteration
        if line.startswith("#"): continue #skip comments
        #added to handle GRCH38 which contains features for genes, transcripts and UTR's which have to be skipped
        #TODO: change this to handle UTR's more precisely
        if line.split("\t")[2] not in ("CDS","exon","start_codon","stop_codon"): continue
        gtf = Feature()
        gtf.readline(line)
        yield gtf


