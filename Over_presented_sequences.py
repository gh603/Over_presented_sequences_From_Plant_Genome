"""
This program is for cis-element analysis based on binomial test.
The user can upload a list of targeted gene and over-presented cis-elements
will be returned to users
"""
from Bio import SeqIO
import scipy.stats as stat
import math
import matplotlib.pyplot as plt
import os
print(os.getcwd())
os.chdir('/Users/gh603/Desktop/Python_project')

def main():
    def program_intro():
        """
        A short description of this test.
        """
        print('='*80)
        print('This program is for identification of over-presented cis-elements\n'+'from the promoter regions of given gene list.\n')
        print('Following species genomes are supported by this program:\n')
        print('1.Zea mays'+' '*10+'2.Arabidopsis thaliana'+' '*10+'3.Vitis vinifera')
        print('='*80)
    
    def parameter():
        """
        Function to ask user specify parameters for test.
        """
        species_ind = eval(input('Please select your targeted species:\n'+'1 for Zea mays; 2 for Arabidposis thaliana; 3 for Vitis vinifera\n'+':'))#ask users selection on organisms
        pro_length = eval(input('Input length of promoter for analysis: '))#ask user to input promoter length
        query_list = input('Input list name that contains your target genes: ')#as user to specify the list containing target genes.
        threshold = eval(input('Input your significant threshold for test: '))
        
        print()
        print('Cis-element analysis specification:')
        print('Your selected organism is: '+ all_species[species_ind-1])
        print('Your selected promoter length is: '+str(pro_length)+'bp')
        print('Your significant threshold for test is: '+str(threshold))
        print('='*80)
        return(species_ind,pro_length,query_list,threshold)
        
    def file_read(name,sep):
        """
        file_read function to automatically read file
        """
        f = open(name)
        file = []
        for line in f:
            if sep==None:
                line = line.rstrip().split()
            elif sep =='\n':
                line = line.rstrip()
            else:
                line = line.rstrip().split(sep)
            file.append(line)
        return(file)
    
    def locusExtract(annotation):
        """
        function to extract the locus information of genes from annotation file
        """
        annotation_file = file_read(annotation,sep=None)#read annotation file in python
        coord_gene = dict()
        for line in annotation_file:
            if line[0].startswith('#'):#any line starts with '#' will be header line which is useless
                continue
            elif line[2] == 'gene':
                gene_info = line[8]#get the 8 column for gene information
                name_ind = gene_info.find('Name')#the index of 'Name' followed by te gene id is returned
                gene_id = gene_info[name_ind+5:]#retrieve the gene ID
                chromosome = line[0]#retrieve the chromosome information
                strand = line[6]
                start = line[3]#retrive the start position
                coord_gene[gene_id] = [chromosome,strand,start]#store gene position to the dictionary.
            else:
                continue 
        print('Coordinates extraction of genes are completed!')
        return(coord_gene)
    
    def promoterExtract(reference, coordinates, length):
        """
        function to extract promoter from genome sequences file. The reference genome sequences, coordinates for each gene in genome, and 
        expected length of promoter regions shoud be specified.
        """
        reference = SeqIO.index(all_geno['Zea mays'],'fasta')#read reference file using biopython packages
        chr_seq = dict()#create an empty dictionary to store the chromosome sequences
        for key in reference.keys():#for loop to append chromosome sequences in to dictionary.
            if reference[key].description.split()[1] == 'dna:chromosome':#sequences followed by 'dna:chromosome mark' will be the whole chromosome sequences
                chr_name = 'chr'+reference[key].id#get chromosome ID
                chr_seq[chr_name]= reference[key].seq#append chromosome sequences in dictionary chr_seq and indexed by its corresponding chromosome ID
        
        all_promoter = dict()#empty dictionary to store the promoter of each gene
        for key in coordinates.keys():
            gene = coordinates[key]
            chr = 'chr'+gene[0]#the first element is the chromosome ID where genes locate
            start_pos = int(gene[2])-length-1#start position of promoter regions
            
            if start_pos <0:#start position cannot be negative
                start_pos = 0
            
            stop_pos = int(gene[2])-1#stop position of promoter region
            if chr not in chr_seq.keys():#only get genes located in chromosome regions provided by reference sequences.
                continue
            else:
                prom_seq = chr_seq[chr][start_pos:stop_pos]#retrieve the promoter sequences using start_pos and stop_pos
                all_promoter[key] = prom_seq#append gene sequences to dictionary 'all_promoter' and indexed by its gene ID
        print('Promoter extraction of genes are completed!')
        return(all_promoter)
    
    def cisExtract(cis_file):
        """
        By providing cis-element file name, functionally characterized cis-element and its corresponding sequences will be extracted
        """
        cis = open(cis_file)#read cis-file
        cis_info = dict()#empty dictionary to store cis-element information
        line = cis.readline()#read the first line of cis-file
        while line:#while loop to keep line will not be empty
            if line[0].startswith('>'):#the first element of each line that startswith '>' is the name of cis-element
                cis_id = line.split()[0]#get the name of of cis-element
                line = cis.readline()#read next line followed by the name of cis-element
                cis_seq = line.split()[0]#the sequence is always the first element.
                cis_info[cis_id]=cis_seq#append the sequences of cis-element into cis_info dictionary indexed by its name
            else:
                line = cis.readline()
        print('Extraction of cis-element names and sequences are completed!')
        return(cis_info)
    
    def cisCount(promoter, cis_file):
        """
        Function to count the occurence of cis-element by providing promoter sequences of genes and cis-element sequences file.
        """
        cis_count = dict()#empty dict to store occurrence of each cis-element
        for cis_key in cis_file.keys():#for loop to get the sequences of each cis-element
            target_cis = cis_file[cis_key]
            count = 0#set the initial count as 0
            for pro_key in promoter.keys():#for loop to count of each cis-element in the given promoter file
                target_pro = promoter[pro_key]#read each individual promoters
                count += target_pro.count(target_cis)#get the count and add to the total count
            cis_count[cis_key] = count#append the total count of each cis-element in the 'cis_count' dict indexed by the cis name
        print('Count of cis-element frequency completed!')
        return(cis_count)
            
    def subset(all_file, subset_index):
        """
        Function to sebset all_file by providing a list of subset_index. All_file must be a dictionary.
        """
        subset_file = dict()
        for index in subset_index:
            if index in list(all_file.keys()):
                subset_file[index] = all_file[index]
            else:
                continue
        return(subset_file)
    
    def dictCount(target_dict):
        """
        Function to get the sum of all the element in a give dictionary
        """
        total = 0
        for value in target_dict.values():
            total+=value
        return(total)
    
    def dictFreq(target_dict):
        """
        Function to calculate the frequency of each element in a given dictionary
        """
        total = dictCount(target_dict)
        freq = dict()
        for key, value in target_dict.items():
            freq[key] = value/total
        return(freq)
        
    def my_binom_test(example_count,cis_count):
        """
        Function to run a binomial test to identify the over-presented or down-presented cis-element
        """
        background_freq = dictFreq(cis_count)#calculate backgroud frequence
        example_total = dictCount(example_count)#get the total cis-element event found in the promoter regions of genes in given list
        example_freq = dict()
        for key,value in example_count.items():
            pvalue = stat.binom_test(value,n=example_total,p=background_freq[key])#binomial test using background frequency as expected frequency.
            example_freq[key] = pvalue
        return(example_freq)
    
    def findSig(test_output,alpha):
        """
        Function to find significant element in given binomial test output file.
        Significant element must have a p-value less than specified alpha.
        """
        sig_cis = dict()#empty dictionary to store significant element
        for key, value in test_output.items():
            if value< alpha:#compare each p-value with alpha.
                sig_cis[key] = value
            else:
                continue
        return(sig_cis)
    
    def my_bar_graph(target_pvalue,threshold):
        """
        Function to plot the binomial results output.
        """
        log_pvalue=dict()
        for key,value in target_pvalue.items():
            log_pvalue[key] = -math.log(value,10)
        plt.bar(range(len(log_pvalue)),sorted(log_pvalue.values(),reverse=True))
        plt.axhline(y = -math.log(threshold,10), color = 'r')
        plt.show()
            
    def writeCis(output,cis,sig_cis):
        """
        Function to write significant cis-element into a local file. One file containing sequences of each cis-element and the significant cis-element file should be included.Output name should be specified.
        """
        w = open(str(output)+'.txt','w')
        line = 'Name'+'\t'+'Sequence'+'\t'+'p-value'+'\n'
        for key in sig_cis.keys():
            w.write(line)
            line = key+'\t'+cis[key]+'\t'+str(sig_cis[key])+'\n'
        w.close()
        
    class species:
        """
        Class to create an object which contains annotation, genome sequences, gene coordinates, promoter sequecesn with given length, counts of cis-elements. It also enables users to subset whole genome data by providing a gene list. If provided with interested genes, binomial test can be conducted by user specifying threshold. Write to local disk file is also supported.
        """
        def __init__(self,species_name,pro_length,cis):
            self.anno = all_anno[species_name]
            self.geno = all_geno[species_name]
            self.length = pro_length
            self.coord = locusExtract(self.anno)#gene coordinates
            self.all_promoter = promoterExtract(self.geno,self.coord,self.length)#extract promoter
            self.cis = cis
            self.cis_count = cisCount(self.all_promoter,cis)#count cis-elements at whole genome level
        
        def queryExtract(self,query):
            self.query = file_read(query,'\n')
            self.query_coord = subset(self.coord,self.query)#coordinates of target genes
            self.query_promoter = subset(self.all_promoter, self.query)#promoter of target_genes
            self.query_cis_count = cisCount(self.query_promoter, self.cis)#count of cis-elements in promoters of target genes.
        
        def sig_test(self,threshold):
            """
            Binomial test to compare target gene with whole gene level
            """
            self.output = my_binom_test(self.query_cis_count, self.cis_count)
            self.sig = findSig(self.output, threshold)
        
        def sig_write(self, file_name):
            writeCis(file_name,self.cis,self.sig)
            print('Significant cis-elements are saved as '+file_name+'.txt')
            print('='*80)
    
    def query(query_list,threshold):
        """
        query function to carry out statistical test on query genes.
        """
        #extract information regarding uploaded target genes
        bg.queryExtract(query_list)
        #carry out statistical test
        bg.sig_test(threshold)
        #plot distribution of significant cis-elment
        askGraph = input("Type 'y' to show cis-element distribution plot, 'n' to skip: ")
        if askGraph.lower() == 'y':
            my_bar_graph(bg.output,threshold)
        #write to file
        askWrite = input("Type 'y' to write file to local disk, 'n' to quit: ")
        if askWrite.lower()=='y':
            name = input('Name your output file: ')
            bg.sig_write(name)   
    
    #Prepare the annotation file name and genome sequences file name.
    all_anno = dict()
    all_anno['Zea mays'] = 'Zmays_284_6a.gene'
    
    all_geno = dict()
    all_geno['Zea mays'] = 'Zea_mays.sample.fa'
    
    #extract cis information(name,sequences) from local disk
    cis = cisExtract('place.fasta')
    all_species = ['Zea mays', 'Arabidopsis thaliana', 'Vitis vinifera']
    
    #print introduction page
    program_intro()
    
    #specify parameters for this test
    species_ind, pro_length, query_list,threshold = parameter()
    species_name = all_species[species_ind-1]
    
    #extract information regarding selected organism
    bg = species(species_name,pro_length,cis)
    #query function
    query(query_list,threshold)
    
    #while function to continously ask user's selection.
    switch = True
    while switch:
        decision = input("Type 'n' to start a new analysis, 'q' to quit: ")#ask user whether they want to quit the test
        if decision == 'n':#'n' to start a new job
            species_ind, pro_length, query_list = parameter()
            species_name = all_species[species_ind-1]
            if all_anno[species_name] == bg.anno and pro_length == self.length:
                query()
            else:
                bg = species(species_name,pro_length,cis)
                query(query_list,threshold)
        elif decision == 'q':#'q' to quit job
            switch = False
        else:#keep in the loop
            continue
    
    print('Have a good day!')
    
main()
        

            
    

