'''Scripts have been written to run on Python 2.7'''

def convert_chr(str):
    '''Converts chromosome to correct format'''
    if len(str)==2:
        chromosome='Pf3D7_'+str+'_v3'
    elif len(str)==1:
        chromosome='Pf3D7_0'+str+'_v3'
    else:
        chromosome=str
    return chromosome

def convert_bedfile(location_filename):
    '''Convert DNAscent files to bedfiles with chromosome IDs converted to PF3D7_01_v3 ... PF3D7_14_v3.
    Complete path to the file to be converted is required ("location_filename").'''
    b=open(location_filename,'r')
    a=b.read()
    a=a.split('\n')
    a=filter(lambda x:x!='',a)
    a=filter(lambda x:x[0]!='#',a)
    a=map(lambda x:x.split(' '),a)
    convert=map(lambda x:[convert_chr(x[0]),int(x[1]),int(x[2]),x[3],x[7]],a)
    convert=sorted(convert,key=lambda x:(x[0],x[1]))
    convert='\n'.join(map(lambda x:'\t'.join(map(lambda y:str(y),x)),convert))
    return convert

def merge_bedfiles(filename_list):
    '''This function can be used to merge bedfiles of runs from the same timepoint.
    Comma separated list of the complete path for each file to be merged is required ("filename_list").
    Prints out a summary of the total number of reads per run.
    This returns a bedfile containing all reads sorted by chromosome number and coordinate.'''
    merged_bed=[]
    filename_list=filename_list.split(',')
    filename_list=filter(lambda x:x!='',filename_list)
    for filename in filename_list:
        b=open(filename,'r')
        a=b.read()
        if '\t' in a:
            a=a.split('\n')
            a=filter(lambda x:x!='',a)
            a=map(lambda x:x.split('\t'),a)
        else:
            a=a.split('\n')
            a=filter(lambda x:x!='',a)
            a=map(lambda x:x.split(' '),a)
        print 'Run'+filename[filename.rfind('Run')+3:-12]+'\t'+str(len(a))
        bedfile=map(lambda x:[convert_chr(x[0]),int(x[1]),int(x[2]),x[3],x[4]],a)
        merged_bed+=bedfile
    merged_bed=sorted(merged_bed,key=lambda x:(x[0],x[1]))
    merged_bed='\n'.join(map(lambda x:'\t'.join(map(lambda y:str(y),x)),merged_bed))
    return merged_bed

def GO_term(location_filename):
    '''This summarises GO terms obtained from PlasmoDB (Aurrecoechea et al., 2009) by total count.
    Complete path to tab separated *.txt file downloaded from PlasmoDB with headers. First two columns should contain "Gene ID" and "source_id", and succeeding columns should be the GO terms (e.g. Computed GO Components, Curated GO Processes).'''
    b=open(location_filename,'r')
    a=b.read()
    a=a.split('\n')
    if '' in a:
        a.remove('')
    a=map(lambda x:x.split('\t'),a)
    for column in a[0][2:]:
        print column
        index=a[0].index(column)
        group=map(lambda x:x[index],a[1:])
        group='\t'.join(group)
        group=group.replace(';','\t')
        group=group.split('\t')
        grp_set=list(set(group))
        grp_set=sorted(grp_set,key=lambda x:(group.count(x)),reverse=True)
        for i in grp_set:
            print i+'\t'+str(group.count(i))
        print '-------'
