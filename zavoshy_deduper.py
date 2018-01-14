#!/usr/bin/python


import argparse

parser = argparse.ArgumentParser(description="A program to remove PCR duplicates from SAM file.")
parser.add_argument("-f", "--file", help='absolute file path to sam file',required=True,type= str)
parser.add_argument("-p", "--paired", help='if data is paired, use this flag with no arguements to indicate paired end data.',required=False,action='store_true')
parser.add_argument("-u", "--umi", help='file with known umis. No option for randomers at this time',required=False,type= str)

args = parser.parse_args()

file=args.file
pair=args.paired
umi_file=args.umi

output=file.split('/')
output=output[len(output)-1]
if pair == True:
    print('paired end capability not available at this time. Run as single end')
    exit()

umi_dict={}
tup_dict={}

def bit_checker(bit):
    '''Takes the bitwise flag and checks for strandedness. Assumes read is mapped, assumes single ended, and that map is sorted. returns "+" or "-" depending on strand''' 
    
   # if (bit & 4) == 4:
        #print('Error: unmapped read detected')
    strand= "+"
    if (bit & 16) == 16: #check to see if the bitwise flag is 16. If it is, then it is in the reverse orientation
        strand= "-"
    #if (bit & 1) == 1: #checking to see if the data is paired or not
        #forwrd='first'
        #if (bit & 128) == 128:
            #forwrd='second'
        #return strand,forwrd
    
    return strand
    
def umi_fun(pos):
    '''checks if the umi is one of the known umis in the provided UMI list. There is no option for randomers at this time.'''
    pos=pos.split(':')
    #print(pos)
    umi=pos[7]
    return umi
def soft_clip(cigar,pos):
    '''checks for softclipping and then edits as necessary'''
    import re
    if re.match( r'([0-9]+)S([0-9]+)M*', cigar):
        match=re.match( r'([0-9]+)S([0-9]+)M*', cigar)
        adj_start=int(pos)- int(match.group(1))
        return adj_start
    else:
        return pos
        
with open(umi_file, 'r') as umis, open(file, 'r') as file:

    #n=0
    for thing in umis:
        thing=thing.strip('\n')

        umi_dict[thing]=0
        #tup_dict[n]=thing
        #n+=1



    for line in file:
        if line.startswith('@'):
            line=line.strip('\n')
            file= open(output+'_deduped', 'a+')
            file.write(line)
            file.close
        else:
            newline=line.strip('\n').split('\t')

            desc=newline[0] #this takes the description string which will be converted with the umi function
            umi=umi_fun(desc)
            pos=newline[3] #looks at the starting position
            chrom=newline[2] #this takes the third column which contains the chromosome this read is on
            bit=bit_checker(int(newline[1]))
            cigar=soft_clip(newline[5],pos)


            info=(chrom,cigar,bit,umi)
            #print(info)

            if umi in umi_dict:
                #print(tup_dict)
                #if tup_dict[info]==umi:
                if info in tup_dict:
                    
                    file= open(output+'_deduped_discard', 'a+')

                    file.write(line)
                    file.close()
                    

                else:

                    tup_dict[info]=0

                    file= open(output+'_deduped', 'a')

                    file.write(line)

                    file.close()
            else:
                file= open(output+'_deduped_discard', 'a+')

                file.write(line)
                file.close()
                
print(umi_dict)
print(len(tup_dict))