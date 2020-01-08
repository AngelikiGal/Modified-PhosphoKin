#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
# 1st Function: Tool for searching kinase binding sites and potential phosphorylated residues in a protein. 
# 2nd Funtion: This scripts also categorizes the phosphorylated residues according to their location relatively to the proteins' known active sites (inside of, close to within six residues proximity and outside of active sites).
# 3rd Function: This script also links kinases to active sites according to the phosphorylated residues they are responsible for and shows their activity inside of, close to and outside of active sites.
# 4th Function: This script finds the possibly responsible kinases for the experimentally observed phosphorylated residues in a protein according to kinases' binding motifs.
# 5th Function: This script identifies the missense point mutations that have been found in cancer patients and possibly disturb the predicted phosphorylation of the protein by either not allowing kinase binding or the phosphorylation according to kinases' motifs.

# Author: Panagiota-Aggeliki Galliou * (ag.gal.work@gmail.com))
# Author: Kleio-Maria Verrou (kleioverrou@yahoo.com)
# * to whom correspondance should be addressed.
## Copyright: © Copyright 2018+, Panagiota-Aggeliki Galliou, All rights reserved.
##License: This script follows the BSD License.


### Note: the kinase motifs are traslated to regular expressions
#### Input: a Sequence.txt file, an Active_sites.txt, an exp_phospho.txt file, a Motifs.txt file and a Mutation.xlxs file.
#### OUTPUT: Five output files; one for the 1st function (.txt file), one for the 2nd function (.txt file), one for the 3rd function (.txt file), one for the 4th function (.txt file) and one for the 5th function (.xls file).

# Edit: 29/6/2019 by Panagiota-Aggeliki Galliou (ag.gal.work@gmail.com)
## changes: -included important residues for kinases motifs, -included recognizing missense point mutations from a file and matching them to the proteins' active sites, kinases' binding sites and phosphorylations, identifying those mutations that possibly disturb the phoshporylation of the protein according to kinases' motifs either by changing an important residue for kinase binding or by changing the residue to be phosphorylated. 
# Edit: 10/7/2019 by Panagiota-Aggeliki Galliou (ag.gal.work@gmail.com)
## changes: -included the Tyrosine Kinases and Tyrosines phosphorylated residues. So the find_phospho_indexes_in_a_motif(motif) function was changed to identify pY as phosphorylated residues in order to recognise the motifs of Tyrosine Kinases.
# Last Edit: 08/01/2020 by Panagiota-Aggeliki Galliou (ag.gal.work@gmail.com)
## changes: -changed the comments and the prints to contain more details and to be more specific,
"""
import re
import time
import sys
import os
import xlrd
import xlwt
import xlsxwriter


# ~ classes ~ #
## ~ ~ for input entries ~ ~ ##
class protein:
    def __init__(self):
        self.name = ""  #name of the protein
        self.type = ""  #Kinase
        self.motifs = [] # list of motifs the protein has. it takes motif objects

class motif:
    def __init__(self):
        self.regex = ""  #regular expression of the motif
        self.pr=[]   #it's an array. it holds numbers that show from the beginning of the motif in how many residues the kinase will phosphorylate. E.g +0 means the kinase will phosphorylate the first recognized residue in the motif. +1 the second secognized residue from the start of the motif etc. So, the motif is [E/D]..pS the value would be +3. 
        self.ir=[]   #it's an array. it  holds numbers that show the important residues in the motif (the residues that are required in the protein sequence for the kinase to bind and phosphorylate there. so basically, anything that is not an X (any aminoacid)). The important residues are presented with numbers, which show from the beggining of the motif (indexed as 0) the position of the important residue. So in a motif XXApS the values would be +3,+4.
        
        
## ~ ~ for result entries ~ ~ ##      
class slim:
    def __init__(self):
        self.start = -1
        self.end = -1
        self.sequence  = -1
        self.phospho_res=[]  #each element of the array shows the position of the residue that will be phosphorylated according to the motif
        self.important_res=[] #the element of the array shows the position of each important residue according to the motif
        
class result_motif:
    def __init__(self):
        self.regex = ""
        self.Slims = [] # the slims found for this motif. It takes slim objects.
        self.M_slims=0 #counts the slims found for this motif. It is the len(Slims). The 0 value it has by default means that the Slims array is empty.
        self.M_p_residues=[] #it has all the phosphorylated residues for this motif.
        self.M_i_residues=[] #it has all the important residues for this motif.
        
class result_protein:
    def __init__(self):
        self.name = ""
        self.type="" #Kinase
        self.Motifs = []
        self.total_slims= 0 #counts all the slims found for a kinase.
        self.total_p_residues=[] #has all the phosphorylated residues for this kinase.
        self.p_residues_inside=[] #has the residues in total_p_residues that are inside of active sites 
        self.p_residues_close=[] #has the residues in total_p_residues that are close within six residues proximity to active sites
        self.p_residues_outside=[] #has the residues in total_p_residues that are outside of active sites
        self.total_i_residues=[] #has all the importabt residues for this kinase.
        self.i_residues_inside=[] #has the residues in total_i_residues that are inside of active sites 
        self.i_residues_close=[] #has the residues in total_i_residues that are close within six residues proximity to active sites
        self.i_residues_outside=[] #has the residues in total_i_residues that are outside of active sites
                
class active_site:
    def __init__(self):
        self.start=0 #start of active site. Shows the position of the residue (eg 160).
        self.end=0   #end of active site. Shows the position of the residue (eg 160).

class result_exp_phosph:
    def __init__(self):
        self.residue="" #is the residue that have been experimentally found to be phosphorylated
        self.responsible_kinases=[] #The kinases that were predicted as possibly responsible for the phosphorylation of this residue. Takes result_kinase objects.


class AMutation(): 
	def __init__(self, SampleID, CancerType, Mutation, FuncImpact, MutationType):
	    self.SampleID = SampleID 
	    self.CancerType = CancerType
	    self.Mutation = Mutation
	    self.FuncImpact = FuncImpact
	    self.MutationType = MutationType

class result_mutation():
   def __init__(self):
        self.SampleID = ""
        self.CancerType = ""
        self.Mutation = ""
        self.FuncImpact = ""
        self.MutationType = ""
        self.locationAS=""  #location to active sites. Either inside, close  or outside
        self.AS = "" #which active site. In case the variable locationAS is inside or close this variable will have a value. in case lthe value of ocationAS is outside, the value of AS will be "-".
        self.dPlaces=[] #a list holding all the motifs, slims, and kinases that this mutation disturbs. 

class dPlace:
    def __init__(self):
        self.residue="" #either "phosphorylated residue" or "important residue" depending on whether the residue that the mutation distrubs is an important or a phosphorylated residue
        self.kinase="" #the kinase that is disturbed by the mutation
        
# ~ END classes ~ #
        
        
# ~ functions ~ # 

def get_sequence (filename):
    ## ~ ~ ~ Opens the Protein.txt file and gets protein sequence ~ ~ ~ ##
    seq=""
    protein_file= open(filename) # here you can put any protein sequence in fasta format
    if protein_file:
        print("Input file ",filename," opened.\n")
        print("\n------------------------------------------------------------------\n")
        for line in protein_file:
            line=line.strip()
            if ">" in line:
                continue
            else:
                seq+=line        
    else:
        raise IOError("Sequence file ",filename," could not open")
    
    protein_file.close()    
    return seq

def get_exp_phosph_residues (filename):
    exp_phosph=[] #an array to hold all the experimentally observed phosphorylations of the protein in a string format
    string="" #buffer string
    Exp_phosph_file= open(filename) # here you can put any txt file with the experimentally observed phosphorylation with a format: e.g. S13, T43, S145
    if Exp_phosph_file:
        print("Input file ",filename," opened.\n")
        for line in Exp_phosph_file:
            line=line.strip()
            string+=line #get everything written in the file in a single string     
    else:
        raise IOError("Sequence file ",filename," could not open")
    
    Exp_phosph_file.close()   
    
    ## edit the string to get the containing experimentally observed phosphorylated residues
    buff=string.split(",")
    
    ## removing the whitespace characters from the start and end of the string and specifically the " " from the start of the string
    for i in range(0,len(buff)):
        buff[i]=buff[i].strip()
        exp_phosph.append(buff[i])

    return exp_phosph

def get_active_sites(filename):

    active_sites=[] #an array to hold all the active_site objects
    string="" #buffer string
    AA_file= open(filename) # here you can put any txt file with the active sites of the protein with the following format: 15-22, 45-53, 63-74
    if AA_file:
        print("Input file ",filename," opened.\n")
        for line in AA_file:
            line=line.strip()
            string+=line #get everything written in the file in a single string     
    else:
        raise IOError("Sequence file ",filename," could not open")
    
    AA_file.close()   
    
    ## edit the string to get the containing active sites
    buff=string.split(",")

    ## removing the whitespace characters from the start and end of the string and specifically the " " from the start of the string
    for i in buff:
        i=i.strip()
        start,end=i.split("-")
        active_sites.append(create_active_site(int(start),int(end)))   
    return active_sites

def get_motifs (filename):
    kinases=[] #an array to hold all the protein objects
    buff=[] #buffer contains in each position the motifs of each kinase
    M_file= open(filename) # here you can put any txt file with the active sites of the protein with the following format: 15-22, 45-53, 63-74
    if M_file:
        print("Input file ",filename," opened.\n")
        for line in M_file:
            line=line.strip()
            if (len(line)>0): #in case there are unnecessary new lines (enters) in the txt file
                buff.append(line)
    else:
        raise IOError("Sequence file ",filename," could not open")
    
    M_file.close()   
    
    ## edit each element in the buffer to get the motifs of each kinase
    for i in buff:
       motifs=[] #an array to hold all the motif objects
       protein, protein_motifs=i.split(":") #to split the protein name from the motifs of each kinase
       protein=protein.strip()
       protein=protein.split(" ")

       if (len(protein)<1):
           raise Exception("Protein name or protein type was not given correctly in ",i)
       
       protein_name=protein[0].strip()
       protein_type=protein[1].strip()
       
       protein_motifs=protein_motifs.split("-") #to split the motifs of each kinase. protein_motifs acts as a buffer array that contains all the motifs of each kinase
              
       for i in protein_motifs:
           i=i.strip()
           regex,pr,ir=find_phospho_indexes_in_a_motif(i)
           motifs.append(create_motif(regex,pr,ir))
       kinases.append(create_protein(protein_name,protein_type,motifs))    
    return kinases


def find_phospho_indexes_in_a_motif(motif):
    #reads a motif and finds the residues to be phosphorylated. Also converts the motif into a regular expression.
    nl = len(motif)
    ir = []   #array that holds the positions for the residues indicated by the motif to be phosphorylated. # list of indexes of importand residues (including phospho..) in the motif
    aS = []   #array that holds the possible residues for each position in the protein.
    ni = -1   #shows the current aminoacid of the aS element
    ii = 0    #shows the current character of the input string
    pr = []   # list of indexes of phosphorylated residues in the motif

    while ii < nl-1:  #check the whole string until one character before the last character
        c = motif[ii]    #current character
        cn = motif[ii+1] #next character

        if (c == "["): #if "["
            ts = c
            c = cn
            ii += 1
            while c != "]":    # [pS*/pT*] #scan the string until the "]"
                ts += c    # save the characters before "]" into the ts string
                ii += 1
                c = motif[ii]
            ts += c       #save the "]" into the string
            aS.append(ts) #append it to the aS

            if "*" in ts: #if "*" in the ts just move on
                ni += 1
            elif ("p" in ts): #if "p" in the ts it means the position will be phosphorylated so append to pr
                ni += 1
                pr.append(ni)
            else:
                ni+=1   #just move on
            ir.append(ni)
#            ni -= 1

        elif (c == "p"):   #if "p" search the next character
            if cn in "STY":  #if "S" or "T" or "Y"
                ni += 1
                ii += 1
                ts = c+cn   #save the characters into a ts sting
                ir.append(ni)  #keep it in importand residues
                if (ii<nl-1) and (motif[ii+1]=="*"):  #-- check if "*" which means that the residues should already be phosphoprylated for the kinase to bind, move on.
                    ts +="*"
                    ii +=1
                else:
                    pr.append(ni)  #If not "*" then its the residue to be phosphorylated so append in pr
                aS.append(ts)
        else:   #else if not "[" or  "p" then it is an amino acid
            ni += 1
            if not ('X' in c):    # if it has a specific aminoacid then add it to ir
                ir.append(ni)
            aS.append(c)  #just append to aS
        ii += 1

    if ii<nl:  #check the last character of the string. this condition is true if the last character is "X" or an amino acid
        c = motif[ii]
        aS.append(c)

        if not ('X' in c):    # if it has a specific aminoacid then add it to ir
            ni +=1
            ir.append(ni)

    #convert the array into a string for output and better handling
    sret = ""
    for ci in aS:
        sret += ci

    #replace certain characters of the string to convert the motif language into regular expressions.
    sret = sret.replace("X", ".")
    sret = sret.replace("/", ",")
    sret = sret.replace("pS*", "S")
    sret = sret.replace("pT*", "T")
    sret = sret.replace("pY*", "Y")
    sret = sret.replace("pS", "S")
    sret = sret.replace("pT", "T")
    sret = sret.replace("pY", "Y")


    if len(pr) < 1:
        raise Exception("A motif should show at least one residue to be phosphorylated.")
    if len(ir) < 1:
        raise Exception("A motif should show at least one important residue for binding.")

    spr = []
    for elem in pr:
        spr.append(str(elem))
    sir = []
    for elem in ir:
        sir.append(str(elem))

    return sret, spr, sir

def create_motif(regex,pr,ir):
    #creates a motif object with values the given arguments.
    #the pi arguments are for the phospho_index value of the motif object and the ir arguments are for the important_res values of the motif object.
    m=motif()
    m.regex=regex
    m.pr=pr
    m.ir=ir
    
    return m
 
def create_protein (name,the_type,motif_array):
    #creates a protein objects with values the given arguments.

    if (len(motif_array)>=1):    # a protein must have at least one motif.
        #create a protein element with a name, a type and the array with its motifs
        p=protein()
        p.name=name
        p.type=the_type
        p.motifs=motif_array
        #return the protein element
    else:
        raise Exception("A protein should have at least one motif")
    return p


def create_slim(seq,st,en,pr,ir):
    #pi argument corresponds to phosphorylated index [counting from the beginning of the slim found]. 
    #ir argument corresponds to the important index [counting fromt the beginning of the slim found].
    s=slim()
    s.sequence=seq  #shows the sequence of the slim
    s.start=st     #shows the location of residue in the start of the slim (found by re.finditer())
    s.end=en       #shows the location of residue in the end of the slim (found by re.finditer())
    
    for pi in pr:
        if(is_correct(s,int(pi))):
            s.phospho_res.append(s.start+int(pi)) 
    
    for ii in ir:
        if(is_correct(s,int(ii))):
            s.important_res.append(s.start+int(ii)) 
    return s

def is_correct(s, index):
    #checks if the slim's residue that the index points to fall out of the slim's limits.
    status=""
    if (s.start<=s.start+index<=s.end):
        status="TRUE"
    else: #out of slim limits
        status="FALSE"
        print("Index ", index, "out of slim limits for ", s.sequence, "!")
    return status

def create_result_motif(regex,array):
    rm=result_motif()
    rm.regex=regex
    rm.Slims=array
    rm.M_slims=len(rm.Slims)
    
    for i in rm.Slims:
        for a in i.phospho_res:
            rm.M_p_residues.append(a)
    for i in rm.Slims:   
        for b in i.important_res:
            rm.M_i_residues.append(b)
            
    return rm
                
def create_result_protein(protein, motif_array):
    rp=result_protein()
    rp.name=protein.name
    rp.type=protein.type
    rp.Motifs=motif_array
    
    for i in rp.Motifs:
        rp.total_slims=rp.total_slims+i.M_slims
        for a in i.M_p_residues:
            rp.total_p_residues.append(a)
 
        for b in i.M_i_residues:
            rp.total_i_residues.append(b)
            
    #convert to set to keep only the unique entries and then to list again for easier manipulation.    
    rp.total_p_residues=set(rp.total_p_residues)
    rp.total_p_residues=list(rp.total_p_residues)
    #sorting the list
    rp.total_p_residues.sort()
    #convert to set to keep only the unique entries and then to list again for easier manipulation.    
    rp.total_i_residues=set(rp.total_i_residues)
    rp.total_i_residues=list(rp.total_i_residues)
    #sorting the list
    rp.total_i_residues.sort()
    #print(rp.total_i_residues)
    
    return rp           
                
def create_active_site(st,en):
    aa=active_site()
    aa.start=st
    aa.end=en
    return aa

def create_result_exp_phosph(residue,responsible_kinases):
    rexp=result_exp_phosph()
    rexp.residue=residue
    rexp.responsible_kinases=responsible_kinases
    return rexp


def create_AMutation(SampleID, CancerType, Mutation, FuncImpact, MutationType):
    mutation_entry = AMutation(SampleID, CancerType, Mutation, FuncImpact, MutationType)
    return mutation_entry

def create_result_mutation(SampleID, CancerType, Mutation, FuncImpact, MutationType, locationAS, AS, dPlaces):
    rm = result_mutation()
    rm.SampleID=SampleID
    rm.CancerType=CancerType
    rm.FuncImpact=FuncImpact
    rm.MutationType=MutationType
    rm.Mutation=Mutation
    rm.locationAS=locationAS
    rm.AS=AS
    rm.dPlaces=dPlaces
    return rm

def create_dPlace(kinase, residue):
    dp=dPlace()
    dp.kinase=kinase
    dp.residue=residue
    return dp
    
    
def find_motifs_for_a_kinase(protein,sequence):
    #finds all the motifs for a protein in a sequence
    slim_array=[] #all the slims found for a motif
    motifs_array=[] #all the motifs found for a protein
    
    for a in protein.motifs: #for each motif in this protein
        slim_buffer=re.findall(a.regex, sequence) #find the regular expression in the given sequence.  
                        
        slim_buffer=set(slim_buffer) #convert it to a set to keep only the unique entries
        slim_buffer=list(slim_buffer) # convert it to a list for easier manipulation
                    
        for x in slim_buffer:                
            for m in (re.finditer(x,sequence)):
                st=m.start()+1  #start of the slim. +1 because list start from 0 but residues from 1.
                en=m.end()+0    #end if the slim.
            
                new_slim=create_slim(x,st,en,a.pr,a.ir) #creates an new slim object
                slim_array.append(new_slim)
       
        slim_buffer=[] #clear slim_buffer
        new_result_motif=create_result_motif(a.regex, slim_array) 
        motifs_array.append(new_result_motif) 
        slim_array=[] #clear slim_array
       
    
    new_result_protein=create_result_protein(protein,motifs_array)
    motifs_array=[] #clear motif_array
    
    return new_result_protein

def find_motifs_for_all_proteins(proteins_array, sequence):
    #finds all motifs of all proteins that are included in a protein_array in a sequence. Returns a Results array with all the results
    Results=[] #has all the result_protein objects.
    
    for i in proteins_array:
        Results.append(find_motifs_for_a_kinase(i,sequence))
       
    return Results
   
def merge_phospho_residues(Results) :
    #merges the total_p_residues arrays of all result_protein objects that are included in the Results array on a single list. Returns that list.
    All_phospho_residues=[]
    for a in Results:
        for b in a.total_p_residues:
            All_phospho_residues.append(b)
    #convert to set to keep only the unique entries and then to list for easier manipulation
    All_phospho_residues=set(All_phospho_residues)
    All_phospho_residues=list(All_phospho_residues)
    All_phospho_residues.sort() #sorting
    return All_phospho_residues

def merge_important_residues(Results) :
    #merges the total_p_residues arrays of all result_protein objects that are included in the Results array on a single list. Returns that list.
    All_important_residues=[]
    for a in Results:
        for b in a.total_i_residues:
            All_important_residues.append(b)
    #convert to set to keep only the unique entries and then to list for easier manipulation
    All_important_residues=set(All_important_residues)
    All_important_residues=list(All_important_residues)
    All_important_residues.sort() #sorting
    return All_important_residues
    
def assign_residues(array,sequence):
    #get a array with the position of residues and returns an array with the residues and their position (e.g. S4)
    return_array=[]
    for i in range (0, len(array)):
        index=array[i]
        residue=sequence[index-1]
        return_array.append(str(residue)+str(index))
    return return_array

def de_assign_residues (array):
    #separates the position from the residue. It takes an array with value eg S231 and return 231.
    return_array=[]
    for i in array:
        return_array.append(int(i[1:]))
    return return_array

def print_results(Results,protein_name,protein_sequence,inside,close,outside,Results_Exp_phosphp,Exp_Inside,Exp_Close,Exp_Outside):   
    
    # ~ ~ for the first output file ~ ~ 
    outputname=protein_name+"-Kinases_Motifs.txt"
    
    with open (outputname, "w") as output:
        if output: 
            print("\nFile ",outputname, "created.\n")
            
            underscore="_________________________________________________________________________________"
            dashes="-------------------------------------------------------------"
            
            All_phospho_residues=merge_phospho_residues(Results)  
            All_important_residues=merge_important_residues(Results)  
            
            All_phospho_residues_with_AA=assign_residues(All_phospho_residues,protein_sequence)
                        
            output.write("{:s}".format(underscore))
            output.write("\nAll phosphorylated residues predicted for {:s}: (Total={:d})\n".format(protein_name,len(All_phospho_residues)))
            output.write(str(All_phospho_residues))

            #output.write("{:d}\n".format(All_phospho_residues))
            output.write("\n{:s}\n".format("or (with assigned residues:)"))
            for i in All_phospho_residues_with_AA: 
                output.write("{:s}, ".format(i))

            output.write("\n{:s}\n".format(underscore))
            
            output.write("{:s}".format(underscore))
            output.write("\nAll Important residues for the predicted phosphorylation of {:s}: (Total={:d})\n".format(protein_name,len(All_important_residues)))
            output.write(str(All_important_residues))
            output.write("\n{:s}".format(underscore))
            
            output.write("\nKinase's Binding Motifs in {:s}:\n".format(protein_name))
            
            for a in Results:
                motif_count=0
                output.write("\n{:s}".format(dashes))
                output.write("\n{:s} {:s} (Total motifs found={:d})\n".format(a.name,a.type, a.total_slims))
                output.write("\n{:5s}{:s} (Total={:d})\n".format(" ","Phosphorylated Residues:",len(a.total_p_residues)))
                ouf=" ["
                for i in range(0,len(a.total_p_residues)):
                    ouf=ouf+str(a.total_p_residues[i])
                    if (i+1<len(a.total_p_residues)):
                        ouf=ouf+", "
                ouf=ouf+"]"    
                output.write("{:5s}{:5s}\n".format(" ",ouf))
                
                output.write("\n{:5s}{:s} (Total={:d})\n".format(" ","Important Residues:",len(a.total_i_residues)))
                ouf=" ["
                for i in range(0,len(a.total_i_residues)):
                    ouf=ouf+str(a.total_i_residues[i])
                    if (i+1<len(a.total_i_residues)):
                        ouf=ouf+", "
                ouf=ouf+"]"    
                output.write("{:5s}{:5s}\n".format(" ",ouf))
                
                
                
                output.write("\n\n{:5s} Motifs:\n".format(" "))
                for b in a.Motifs:
                    output.write("\n{:8s}{:3d}) {:s} (Total={:d})\n".format(" ",motif_count+1,b.regex,b.M_slims))
                    motif_count=motif_count+1
                    slim_count=0
                        
                    for c in b.Slims: 
                        ouf1=""
                        for i in range(0,len(c.phospho_res)):
                            ouf1=ouf1+str(c.phospho_res[i])
                            if (i+1<len(c.phospho_res)):
                                ouf1=ouf1+", "
                                                            
                        ouf2=""
                        for i in range(0,len(c.important_res)):
                            ouf2=ouf2+str(c.important_res[i])
                            if (i+1<len(c.important_res)):
                                ouf2=ouf2+", "
                        
                        output.write("{:11s}{:3d}) {:5s} {:4d} - {:4d}: [Phosphorylated residue(s)= {:4s}] [Important residue(s)= {:4s}]\n".format(" ",slim_count+1,c.sequence,c.start,c.end,ouf1,ouf2))    
                        slim_count=slim_count+1    
                    output.write("\n{:13s}Motif's phosphorylated residues: (Total={:d})\n".format(" ",len(b.M_p_residues)))
                    ouf=" ["
                    for i in range(0,len(b.M_p_residues)):
                        ouf=ouf+str(b.M_p_residues[i])
                        if (i+1<len(b.M_p_residues)):
                            ouf=ouf+", "
                    ouf=ouf+"]"    
                    output.write("{:13s}{:s}\n".format(" ",ouf))
                    
                    output.write("\n{:13s}Motif's important residues: (Total={:d})\n".format(" ",len(b.M_i_residues)))
                    ouf=" ["
                    for i in range(0,len(b.M_i_residues)):
                        ouf=ouf+str(b.M_i_residues[i])
                        if (i+1<len(b.M_i_residues)):
                            ouf=ouf+", "
                    ouf=ouf+"]"    
                    output.write("{:13s}{:s}\n".format(" ",ouf))
                    
                #motif_count=0 #set again in zero to count the next protein's motifs
                output.write("{:s}\n".format(dashes))
            
            print("Results written on ",outputname," file.")      
        else:
            raise IOError("Output file ",outputname," not created")
    output.close()
    
    
    # ~ ~ for the second output file ~ ~ 
    outputname=protein_name+"-Categorization_of_phosphorylations_comparatively_to_active_sites.txt"
    
    with open (outputname, "w") as output:
        if output: 
            print("\nFile ",outputname, "created.\n")
            
            output.write("{:s}-Categorization of phosphorylated residues relatively to active sites \n(inside of, close to within six residues proximity and outside of active sites)\n".format(protein_name))
            output.write("\n{:s}\n{:s}\n".format("Experimentally observed phosphorylated residues:",underscore))
            output.write("Phosphorylated residues Inside Active sites: (Total={:d})\n".format(len(Exp_Inside)))
            for i in Exp_Inside:
                output.write("{:s}, ".format(i))
            output.write("\n{:s}\n\nPhosphorylated residues Close within six residues proximity to Active sites: (Total={:d})\n".format(dashes,len(Exp_Close)))
            for i in Exp_Close:
                output.write("{:s}, ".format(i))
            output.write("\n{:s}\n\nPhosphorylated residues Outside of Active sites: (Total={:d})\n".format(dashes,len(Exp_Outside)))
            for i in Exp_Outside:
                output.write("{:s}, ".format(i))
            output.write("\n{:s}\n\n".format(underscore))    
                   
            output.write("\n\n\n{:s}\n{:s}\n".format("Predicted observed phosphorylated residues:",underscore))   
            output.write("\nAll phosphorylated residues: (Total={:d})\n".format(len(All_phospho_residues)))
            for i in All_phospho_residues_with_AA:
                output.write("{:s}, ".format(i))
            output.write("\n\n{:s}\nPhosphorylated residues Inside Active sites: (Total={:d})\n".format(dashes,len(inside)))
            for i in assign_residues(inside,protein_sequence):
                output.write("{:s}, ".format(i))
            output.write("\n{:s}\n\nPhosphorylated residued Close within six residues proximity to Active sites: (Total={:d})\n".format(dashes,len(close)))
            for i in assign_residues(close,protein_sequence):
                output.write("{:s}, ".format(i))
            output.write("\n{:s}\n\nPhosphorylated residued Outside of Active sites: (Total={:d})\n".format(dashes,len(outside)))
            for i in assign_residues(outside,protein_sequence):
                output.write("{:s}, ".format(i))
            output.write("\n{:s}".format(underscore))   
            print("Results written on ",outputname," file.") 
        else:
            raise IOError("Output file ",outputname," not created")
    output.close()
    
    
    # ~ ~ for the third output file ~ ~ 
    outputname=protein_name+"-Link_ActiveSite_with_kinases.txt"
    
    with open (outputname, "w") as output:
        if output: 
            print("\nFile ",outputname, "created.\n")
            
            output.write("{:s}\n".format(dashes))
            for i in Results:
                output.write("{:s} {:s}:\n".format(i.name,i.type))
                output.write("\n{:5s}Phosphorylated residues Inside of Active sites: (Total={:d})\n{:5s}".format(" ",len(i.p_residues_inside)," "))
                for a in assign_residues(i.p_residues_inside,protein_sequence):
                    output.write("{:s}, ".format(a))
                output.write("\n\n{:5s}Phosphorylated residues Close to Active sites (within 6 residues proximity): (Total={:d})\n{:5s}".format(" ",len(i.p_residues_close)," "))
                for a in assign_residues(i.p_residues_close,protein_sequence):
                    output.write("{:s}, ".format(a))
                output.write("\n\n{:5s}Phosphorylated residues Outside of Active site: (Total={:d})\n{:5s}".format(" ",len(i.p_residues_outside)," "))
                for a in assign_residues(i.p_residues_outside,protein_sequence):
                    output.write("{:s}, ".format(a))
                output.write("\n{:s}\n\n".format(dashes))
            print("Results written on ",outputname," file.") 
        else:
            raise IOError("Output file ",outputname," not created")
    output.close()    
    
    # ~ ~ for the fourth output file ~ ~ 
    outputname=protein_name+"-Possily_responsible_kinases.txt"
    with open (outputname, "w") as output:
        if output: 
            print("\nFile ",outputname, "created.\n")
            
            for i in Results_Exp_phosphp:
                output.write("{:s}:\n".format(i.residue))
                for a in i.responsible_kinases:
                    output.write("{:5s}{:s}\n".format("",a.name))
                output.write("\n")
            print("Results written on ",outputname," file.") 
        else:
            raise IOError("Output file ",outputname," not created")
    output.close()  
    return

def find_phosphorylations_inside(array):
    #finds the residues in the array argument that are inside of active sites appends them in a list and returns that list.
    inside_list=[]
    for i in array:
        for x in Active_sites:
             if x.start<=i<=x.end: # residue is inside active sites
                 inside_list.append(i)
    return inside_list

def find_phosphorylations_close(array, inside_list):
    #finds the residues in the array argument that are close within six residues proximiry to active sites appends them in a list and returns that list.
    close_list=[]
    for i in array:
        if (i not in inside_list): #if residue is not inside then is either close to or outside of active sites
            for x in Active_sites:
                if (x.start-6<=i<x.start or x.end<i<=x.end+6):  #residue close within 6 residues proximity to active sites
                    close_list.append(i)
    return close_list
   
def find_phosphorylations_outside (array, inside_list, close_list):
    #finds the residues in the array argument that are outside of active sites appends them in a list and returns that list.
    outside_list=[]
    for i in array:
        if ((i in inside_list) or (i in close_list)): 
            pass
        else: #residue not inside neither close to active sites, therefore residue is outside active sites
            outside_list.append(i)
    return outside_list
             
def categorize_phosphorylations_relatively_to_active_sites(all_phosphorylations):
    #categorizes the residues in all_phosphorylations in three lists (inside, close and outside) accordind to their location relatively to active sites(inside active sites, close within 6 residues proximity to active sites and outside of active sites)
    Inside=find_phosphorylations_inside(all_phosphorylations)
    Close=find_phosphorylations_close(all_phosphorylations, Inside)
    Outside=find_phosphorylations_outside(all_phosphorylations,Inside,Close)
    
    return Inside,Close,Outside
 
def link_kinases_with_active_sites(Results,inside,close,outside):
    #Splits the total_p_residues of each result_kinase in the Results in three lists in the result_kinase that until now were empty.
    #The lists are p_residues_inside, p_residues_close and p_residues_outside and take the residues of the kinase that are inside, close and outside of active sites,respectively.
    for i in Results:
        for x in i.total_p_residues:
            if(x in inside):
                i.p_residues_inside.append(x)
            else: #x not inside so it is either close or outside of active sites
                if (x in close):
                    i.p_residues_close.append(x)
                else: #x not close so it is outside of active sites
                    if (x in outside):
                        i.p_residues_outside.append(x)
                    else: 
                        print("Residue ",x," is not inside nor close nor outside of active sites!" )
    return Results #returns again the Results but now the p_residues_inside, p_residues_close and p_residues_outside of each result_kinase in the Result are not empty but have values.

def find_responsible_kinases(Exp_phosph,Results):
    #finds the possibly responsible kinases for the experimentally observed phosphorylations according to the predicted phosphorylated residues (total_p_residues attribute) of each result_protein object in Results.
    Responsible_kinases=[] #buffer that contains all the kinases that could phosphorylate each residue. It takes result_kinase objects.
    Results_Exp_phosph=[] #Has all the result_exp_phosph objects
    
    for i in Exp_phosph:
        residue_num=i[1:]
        for a in Results:
           if(len(a.total_p_residues)>0): #if the kinase is predicted to phosphorylate residues in the protein
               if(int(residue_num) in a.total_p_residues):
                   Responsible_kinases.append(a)
        Results_Exp_phosph.append(create_result_exp_phosph(i,Responsible_kinases))
        Responsible_kinases=[] #empty list for the next Exp_phosph
    return Results_Exp_phosph

def get_Mutations(TableFile):
    #it parses an excel (.xlsx) file with mutations mined from CBioPortal. It only gets the point mutations.
    Mutations=[] #List holding all the mutation objects
    try:    
        wb=xlrd.open_workbook(TableFile)
        sheet=wb.sheet_by_index(0)
        
        try:
            for i in range(sheet.nrows):            
                if ("Missense_Mutation" in sheet.row_values(i)): #only getting the Missense mutations and they result in a different amino acid. If a mutation is nonsense will result in the same amino acid that will not disturb the binding of kinases in the protein.
                    MutationType="Missense_Mutation" 
                   
                    if("deleterious" in sheet.cell_value(i,5)): #the column 5 should have this information in each row in the excel file
                        FuncImpact="deleterious"
                    if ("tolerated" in sheet.cell_value(i,5)): #the column 5 should have this information in each row in the excel file
                        FuncImpact="tolerated"
                                       
                    sampleID=str(sheet.cell_value(i,1)).strip()
                    cancerType=str(sheet.cell_value(i,2)).strip()
                    mutation= str(sheet.cell_value(i,3)).strip()
                   
                    Mutation=create_AMutation(SampleID=sampleID, CancerType=cancerType, Mutation=mutation, FuncImpact=FuncImpact, MutationType=MutationType)
    
                    Mutations.append(Mutation)
    
            #print(len(Mutations))
        except:
            print("\nSomething went wrong and could get the mutations.. Check the format of values of the excel file")
            sys.exit() #telling to the user were the problem is, and quiting from the program
    except:
        print("\nSomething went wrong and could not open file... Check the format of the excel file with the mutations")
        sys.exit() #telling to the user were the problem is, and quiting from the program

    return Mutations   

def assignResidueToAS(my_mutation):
    #assigns a mutation to an active site if it is inside or close to it, and if it is outside of all active sites it marks the mutation as outside of active sites
    locationToAS=""
    activeSite=""

    #put the mutation in an array to give it as an argument to the function categorize_phosphorylations_relatively_to_active_sites
    my_mutationArray=[]
    my_mutationArray.append(int(my_mutation))
    inside, close, outside= categorize_phosphorylations_relatively_to_active_sites(my_mutationArray)
    #because the my_mutationArray has only one value, only one of the arrays inside close, outside will have a length 1 and the rest two will have a length 0. The array that has the value 1 shows the location of the mutation relatively to active sites
    if (len(outside)>0): #mutation is outside
        locationToAS="Outside"
        activeSite="-"
    elif (len(close)>0):
        locationToAS="Close" #mutation is close
        for i in Active_sites:
            if((int(i.start)-6<=int(my_mutation)<int(i.start)) or (int(i.end)<int(my_mutation)<=int(i.end)+6)):
                activeSite= str(i.start) + "-" + str(i.end)
                break
    elif (len(inside)>0): #mutation is inside
        locationToAS="Inside"
        for i in Active_sites:
            if(int(i.start)<=int(my_mutation)<=int(i.end)): #mutation is inside an active site
                activeSite= str(i.start) + "-" + str(i.end)
                break
    else:
        print("Something went wrong. A mutation should either be inside, close to or outside of an active site.")
    return locationToAS,activeSite   

def mutationExists(mutation_entry, array):
    #finds if a mutation entry already exists in an array
    status="does not exist"
    for i in array: # two mutation entries are considered the same if they have the same Mutation value, the same SampleId value and the same Cancer type value.
        if (i.Mutation==mutation_entry.Mutation):
            if(i.SampleID==mutation_entry.SampleID):
               if(i.CancerType==mutation_entry.CancerType): 
                   status="exists"
    return status

def only_keep_unique_mutations (Result_Mutations):
    #it keeps only the unique entries of the Result_Mutations array and returns a new array with only the unique mutations
    new_array=[]

    i=0
    new_array.append(Result_Mutations[i])

    for x in range(i+1,len(Result_Mutations)):
        if(mutationExists(Result_Mutations[x],new_array)=="does not exist"):
            new_array.append(Result_Mutations[x])
                         
    return new_array
    
def find_Mutations(Results,MutationList,protein_name):
    
    All_important_residues=merge_important_residues(Results)    
    Result_Mutations=[]
    
    for i in MutationList:
        my_mutation=str(i.Mutation[1:len(i.Mutation)-1])
    
        dResidue=""
        dKinase=""
        locationAS=""
        AS=""      
        dPlaces=[]
        
        for x in All_important_residues:          
            #if mutation is in All_important_residues [important residues include phosphorylated residues]
            if (str(my_mutation)==str(x)):
                #then search if for which kinase it is an important residue.
                for a in Results:                     
                    for y in a.total_i_residues:
                        if(str(my_mutation)==str(y)):
                            dKinase=a.name
                            dResidue="important residue"
                            #search if mutation is in total_p_residues because then it would be disturbing a phosphorylated residue
                            for y in a.total_p_residues:
                                if(str(my_mutation)==str(y)):
                                    dResidue="phosphorylated residue"
        
                            dPlaces.append(create_dPlace(dKinase,dResidue))
                            break

                locationAS,AS=assignResidueToAS(my_mutation)
                Result_Mutations.append(create_result_mutation(i.SampleID, i.CancerType, i.Mutation, i.FuncImpact, i.MutationType, locationAS, AS, dPlaces))             
                
    Result_Mutations=only_keep_unique_mutations(Result_Mutations)
    print_mutations_Results(Result_Mutations,protein_name)                        
             
    return

def print_mutations_Results(Result_Mutations,protein_name):
    # ~ ~ for the fifth output file ~ ~ 
    outputname=protein_name+"-Mutations_disturbing_Phosphorylations.xls"
      
    try:    
        wb=xlwt.Workbook(encoding = 'ascii')
    
        sheet=wb.add_sheet("Missense Mutations")
            
        sheet.write(0,0,"Sample ID")
        sheet.write(0,1,"Cancer Type")
        sheet.write(0,2,"Mutation")
        sheet.write(0,3,"Amino Acid Change")
        sheet.write(0,4,"Functional Impact")
        sheet.write(0,5,"Location to AS")
        sheet.write(0,6,"Active Site (AS)")
        sheet.write(0,7,"No. disturbed Kinases")
        sheet.write(0,8,"disturbed Kinases")
        sheet.write(0,9,"disturbed Residue")
        
        rowindex=1
        for i in Result_Mutations:
            sheet.write(rowindex,0, i.SampleID)
            sheet.write(rowindex,1, i.CancerType)
            mutation=str(i.Mutation[1:len(i.Mutation)-1])
            sheet.write(rowindex,2, mutation)
            aaChange=str(i.Mutation[0]+" -> "+i.Mutation[-1])
            sheet.write(rowindex,3, aaChange)
            sheet.write(rowindex,4, i.FuncImpact)
            sheet.write(rowindex,5, i.locationAS)
            sheet.write(rowindex,6, i.AS)
            sheet.write(rowindex,7, len(i.dPlaces))
            
            increase=0
            for x in i.dPlaces:
                sheet.write(rowindex,8+increase, x.kinase)
                sheet.write(rowindex,9+increase, x.residue)
                increase=increase+2        
            
            rowindex=rowindex+1
            
        wb.save(outputname)
        print("Results written on ",outputname," file.") 
        
       
    except:
        print("\nSomething went wrong and could write result mutations in an excel file...")
        sys.exit() #telling to the user were the problem is, and quiting from the program   

    return

def search_in_protein(protein_sequence):
    #one function that does it all. Get the sequence, finds the motifs of all kinases in the sequence, categorizes the phosphorylated residues according to their position relatively to active sites, links the kinases with active sites, finds the possibly responsible kinases for the experimentally observed phosphorylations and sents all the results from wrtting in output files.
    seq=get_sequence(protein_sequence)
    protein_name=protein_sequence.split(".",1)
    protein_name=str(protein_name[0])
    
    Results=find_motifs_for_all_proteins(All_proteins,seq)
    Inside,Close,Outside=categorize_phosphorylations_relatively_to_active_sites(merge_phospho_residues(Results)) 
    #in the "new" Results the p_residues_inside,p_residues_close and p_residues_outside arguments are not empty lists but contain the kinase's residues that are inside, close and outside of active sites, correspondingly.
    Results=link_kinases_with_active_sites(Results,Inside,Close,Outside)
    Results_Exp_phosphp=find_responsible_kinases(Exp_Phosph,Results)
    Exp_Inside, Exp_Close, Exp_Outside=categorize_phosphorylations_relatively_to_active_sites(de_assign_residues(Exp_Phosph))
    #assign residues to them again for writting
    Exp_Inside=assign_residues(Exp_Inside,seq)
    Exp_Close=assign_residues(Exp_Close,seq)
    Exp_Outside=assign_residues(Exp_Outside,seq)
    #sent it all for writting on a file
    print_results(Results,protein_name,seq,Inside,Close,Outside,Results_Exp_phosphp,Exp_Inside,Exp_Close,Exp_Outside)
    
    
    find_Mutations(Results,MutationList,protein_name)
    return 



# ~ END functions ~ # 

# ~ Main ~ #

##Inform user for the origin, function and requirments of the code

print("\nThis code is a modification of the PhosphoKin tool \n(Galliou, P.A., Verrou, K.M., 2019. An in silico method for studying the phosphorylation in association to active sites. Aristotle Biomedical Journal 1, 48–59). \n")
print("This code was implemented by Galliou Panagiota-Angeliki (email: ag.gal.work@gmail.com).\n")
print("Briefly, the PhosphoKin tool predicts phosphorylation sites and phosphorylated residues in a protein sequence based on kinase binding motifs.\n")
print("This code, modified the PhosphoKin tool to read cancer mutations found in the protein sequence and identify those deleterious missense point mutations that could obstruct the predicted phosphorylation of the protein according to the kinases' binding motifs, by either occuring directly on a predicted phosphorylated residue or  occuring on required residues for kinase binding in a phosphorylation site.")
print("\nThe code takes 5 files as input: \n 1) A file with the active sites of the protein.\n 2) A file with the experimentally observed phosphorylated residues in the protein.\n 3) A file with the motif(s) of kinases(s) for which the user wants to find binding sites in the protein sequence.\n 4) A file with the sequence of the protein in a fasta format.\n 5) An xlsx file (excel) with the mutations of the protein that have been found in all types of cancer as downloaded from CBioPortal.\n")
print("The code creates 5 output files: \n 1) A file that shows the phosphorylation sites and phosphorylated residues for each given motif of each given protein.\n 2) A file that categorizes the phosphorylated residues in the protein according to their position relatively to the active sites of the protein.\n 3) A file that links the given kinases with the active sites of the protein.\n 4) A file that shows the possibly responsible kinase(s) for each experimentally observed phosphorylated residue in the protein.\n 5) An excel file with those missense point mutations found in the protein in cancer patients that could obstruct the protein's predicted phosphorylation either by occuring on the very exact amino acid that the kinase could phosphorylate or by occuring on the required surrounding residues for kinase binding to the phosphorylation site, according to kinases' motifs.\n")

print("The name of the above 5 output files are:\n 1) LAMC1-Kinases_Motifs.txt \n 2) LAMC1-Categorization_of_phosphorylations_comparatively_to_active_sites.txt \n 3) LAMC1-Link_ActiveSite_with_kinases.txt \n 4) LAMC1-Possily_responsible_kinases.txt \n 5) LAMC1-Mutations_disturbing_Phosphorylations.xls\n")
print("The 5 input files used by this code are uploaded in the following link:\nhttps://www.dropbox.com/sh/5p0b4fkri1uisbl/AABm9MqYvhaxEinaPKVTgfMGa?dl=0\n\nFor details on the format of the first 4 files, please refer to the manuscript of PhosphoKin Tool \n(Galliou, P.A., Verrou, K.M., 2019. An in silico method for studying the phosphorylation in association to active sites. Aristotle Biomedical Journal 1, 48–59).\n\nThe 5th input file should be formated as the LAMC1-Mutations.xlxs file in the above link \n(https://www.dropbox.com/s/e2paiz1jvhc8t07/LAMC1-Mutations.xlsx?dl=0).\n")
print("The input files should be in the same directory as the code upon running. \nThe output files are created in the directory of the code upon running.\n\n")
print("\n------------------------------------------------------------------\n\n")  
"""
## Take the names of the required files as input from the user.  
print("[ATTENTION: Please read the manual to make sure what format the required files should have.]\n")


### ~ ~ ~ Asking from user files and opening each file ~ ~ ~ ###
print ("\nNow, I would like you to give me the 5 files.\n ")



###Asking for the active sites
active_sites_file=input("\nGive the Active file, please: ")
if active_sites_file == 'help':
	print("\nThe file must have the following format:\n Start_of_active_site_1 - End_of_active_site_1,  Start_of_active_site_2 - End_of_active_site_2, etc \n e.g. 13-22, 150-162, 1147-1458, etc\n")
	active_sites_file=input("\nGive the Active Sites file, please: ")
try: 
    #active_sites_file = "Active_sites.txt"
    Active_sites=get_active_sites(active_sites_file)  #Holds all the active site objects
except:
	print("\nSomething went wrong... Check the format of the Active Sites file")
	sys.exit() #telling to the user were the problem is, and quiting from the program
    
###Asking for the experimentally observed phosphorylated residues 
exp_phospho_file=input("\nGive the Experimentally observed phosphorylated residues in the protein file, please: ")
if exp_phospho_file == 'help':
	print("\nThe file must have the following format:\n phosphorylated_residue_1, phosphorylated_residue_2, phosphorylated_residue_3, etc..\n e.g S4, S156, T445, etc \n ")
	exp_phospho_file=input("\nGive the Experimentally observed phosphorylated residues in the protein file, please: ")
try: 
    Exp_Phosph=get_exp_phosph_residues(exp_phospho_file) #Holds all the experimentallt observed phosphorylated residues as strings
except:
	print("\nSomething went wrong... Check the format of the Experimentally observed phosphorylated residues in the protein file")
	sys.exit() #telling to the user were the problem is, and quiting from the program

###Asking for the Motifs
motifs_file=input("\nGive the Motifs file, please: ")
if motifs_file == 'help':
	print("\nThe file must have the following format: \n Name_of_kinase_1[space]protein type: motif_1, motif_2, motif_3, etc [enter] Name_of_kinase_2[space]protein type: motif_1,motif_2, etc [enter]\n H1K Kinase:[pS/pT]P[R/K]-[pS/pT]PX[R/K]-[R/K][pS/pT]P \nATM Kinase:pSQ-[P/L/I/M]X[L/I/E/D]pSQ-LpSQE\n")
	motifs_file=input("\nGive the Motifs file, please: ")
try: 
    All_proteins=[] # Holds all the protein objects.
    All_proteins=get_motifs(motifs_file)
except:
	print("\nSomething went wrong... Check the format of the Motifs file")
	sys.exit() #telling to the user were the problem is, and quiting from the program

###Asking for the Mutations
mutations_file=input("\nGive the Mutations file, please: ")
if mutations_file == 'help':
	print("\nThe file must have the format as downloaded by CBioPortal \n [The columns should be in the following order: Sample ID, Cancer Type,	Protein Change,	Annotation,	Functional Impact,	Mutation Type,	Copy, COSMIC, MS, VS, Center, Chromosome, Start Pos, End Pos, Ref, Var, Allele Freq (T), Allele Freq (N), Variant Reads, Ref Reads,	Variant Reads (N),	Ref Reads (N),	# Mut in Sample] ")          
	mutations_file=input("\nGive the Mutations file, please: ")
try: 
    MutationList=[] # Holds all the protein objects.
    MutationList=get_Mutations(mutations_file)
except:
	print("\nSomething went wrong... Check the format of the Mutations file")
	sys.exit() #telling to the user were the problem is, and quiting from the program 

###Asking for the Seqeunce 
sequence_file=input("\nGive the Protein Sequence file, please: ")
if sequence_file == 'help':
	print("\nThe file must have a fasta format (https://en.wikipedia.org/wiki/FASTA_format).\n ")
	sequence_file=input("\nGive the Protein Sequence file, please: ")
try: 
    search_in_protein(sequence_file) 
except:
	print("\nSomething went wrong... Check the format of the Sequence file")
	sys.exit() #telling to the user were the problem is, and quiting from the program

print("\n\n------------------------------------------------------------------\n\n")
"""
##Taking the input files.
Active_sites=get_active_sites("LAMC1-Active_Sites.txt") #Takes the file with the active sites as input
Exp_Phosph=get_exp_phosph_residues("LAMC1-Experimentally_Observed_Phosphorylations.txt") #Takes the file with the experimentally observed phosphorylations as input
All_proteins=[] # Holds all the protein objects.
All_proteins=get_motifs("LAMC1-Kinase_Motifs.txt") #Takes the file with the kinases' recognition motifs as input
MutationList=get_Mutations("LAMC1-Mutations.xlsx") #Takes the file with the cancer mutations in the protein as input
search_in_protein("LAMC1.txt") #Takes the file with the protein sequence as input
# ~ END Main ~ #
    