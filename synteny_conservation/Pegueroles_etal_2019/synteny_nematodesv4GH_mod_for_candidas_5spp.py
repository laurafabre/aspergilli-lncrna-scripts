
import sys,re

def recodeList(file1, file2, string):#gene order list for spX, orthology file, ex: 'calb'
    #returns a list with the gene order of a given sp were geneIDs were renamed to myID 
    #to allow direct comparisson between species;
    #genes with no homology are also included!
    #ex: ['myID1', 'WBGene000xxx' 'myID2', 'myID3', 'myID4', 'myID5','XLOC_001892', 'myID6', 'myID7', ...]  
    mylist = []
    matched=0; nonmatched=0; found=0; lncRNA=0
    for row1 in open(file1, 'r'):  
        row1 = row1.rstrip()
        # ~ if re.search(r'MSTRG', row1):
        if row1.startswith("MSTRG"):
            #print(row1)
            mylist.append(row1)
            found=1
            lncRNA=lncRNA+1
        for row2 in open(file2 , 'r').readlines():  
            row2= row2.rstrip().split("\t")
            if string=='afla':
                if row1== row2[0]:
                    mylist.append(row2[4])
                    matched=matched+1
                    found=1
                    break
            if string=='afum':
                if row1== row2[1]:
                    mylist.append(row2[4])
                    matched=matched+1
                    found=1
                    break
            if string=='anid':
                if row1== row2[2]:
                    mylist.append(row2[4])
                    matched=matched+1
                    found=1
                    break
            ##    
            if string=='anig':
                if row1== row2[3]:
                    mylist.append(row2[4])
                    matched=matched+1
                    found=1
                    break
                    
        if found==0:
            nonmatched=nonmatched+1
            if nonMatch=='yes':
                mylist.append(row1)
        found=0
    print "number of matched IDs for",string,": " ,matched
    print "number of NON matched IDs for",string,": " ,nonmatched
    print "number of candidate lncRNA",string,": " ,lncRNA,"\n" 
    return mylist

def dictionaryOfClusters(myidx,mylist):#myidx=positions in myList for each lncRNA; myList=renamed gene order
    #to create a dict; for each lncRNA (key) stores the number of nearby renamed geneID indicated by the user  
    mydict={}
    for idx in myidx:
        key=mylist[idx]
        val={'left':[], 'right':[], 'all':[]}
        if not key in mydict:
            mydict[key]=val
        i=1
        while i <=genesNearby:
            try:
                if idx-i >=0:                 
                    mydict[key]['left'].append(mylist[idx-i])
                    mydict[key]['all'].append(mylist[idx-i])
                mydict[key]['right'].append(mylist[idx+i])                                   
                mydict[key]['all'].append(mylist[idx+i])
            except IndexError:
                pass
            i=i+1
     
    return mydict
    
def comparingDict(sp1, sp2):
    #compares dictionaries from dictionaryOfClusters() for two species; 
    #if the number of shared genes is >= minOverlap, lncRNA are stored in myHomologs list
    if sp1=='afla':
        dict1= dictionaryOfClusters(afla_idx,aflaList)
    if sp1=='afum':
        dict1= dictionaryOfClusters(afum_idx,afumList)
    if sp1=='anid':
        dict1= dictionaryOfClusters(anid_idx,anidList)
    if sp1=='anig':
        dict1= dictionaryOfClusters(anig_idx,anigList)
       
        
    if sp2=='afla':
        dict2= dictionaryOfClusters(afla_idx,aflaList)
    if sp2=='afum':
        dict2= dictionaryOfClusters(afum_idx,afumList)
    if sp2=='anid':
        dict2= dictionaryOfClusters(anid_idx,anidList)
    if sp2=='anig':
        dict2= dictionaryOfClusters(anig_idx,anigList)
    myHomologs=[]
    homologFound='false'
    
    
    for key1, val1 in dict1.iteritems():
        #print key1, val1
        for key2, val2 in dict2.iteritems():
            if len(set(dict1[key1]['all']).intersection(dict2[key2]['all'])) >=minOverlap:
            #at least one of the two lncRNA share genes in the left and the right side    
            #forces that both lncRNA share genes in the left and the right side    
                if len(set(dict1[key1]['right']).intersection(dict2[key2]['right'])) >=minSideOverlap:
                    if len(set(dict1[key1]['left']).intersection(dict2[key2]['left'])) >=minSideOverlap: 
                        homologFound='true'
                if len(set(dict1[key1]['right']).intersection(dict2[key2]['left'])) >=minSideOverlap:
                    if len(set(dict1[key1]['left']).intersection(dict2[key2]['right'])) >=minSideOverlap:
                        homologFound='true'
                if homologFound=='true':             
                        mytup=(key1,key2)
                        myHomologs.append(mytup)
                        homologFound='false'
    return myHomologs

#####################################

in1= sys.argv[1] #gene order list for sp1
in2= sys.argv[2] #gene order list for sp2
in3= sys.argv[3] #gene order list for sp3
in4= sys.argv[4] #gene order list for sp4
in6= sys.argv[5] #orthology file
out= open(sys.argv[6] , 'w') #out file
genesNearby= int(sys.argv[7])
minOverlap= int(sys.argv[8])
minSideOverlap= int(sys.argv[9])
nonMatch = sys.argv[10]
temp= open('temp', 'w')

aflaList=recodeList(in1, in6, 'afla')
afumList=recodeList(in2, in6, 'afum')
anidList=recodeList(in3, in6, 'anid')
anigList=recodeList(in4, in6, 'anig')


for x in aflaList:
    temp.write("%s\tafla\n" %(x))       
for x in afumList :
    temp.write("%s\tafum\n" %(x))    
for x in anidList:
    temp.write("%s\tanid\n" %(x))
for x in anigList:
    temp.write("%s\tanig\n" %(x))    


afla_idx = [i for i, item in enumerate(aflaList) if item.startswith('MSTRG')]
#list containing the positions in calbList for each lncRNA
afum_idx = [i for i, item in enumerate(afumList) if item.startswith('MSTRG')]
anid_idx = [i for i, item in enumerate(anidList) if item.startswith('MSTRG')]
anig_idx = [i for i, item in enumerate(anigList) if item.startswith('MSTRG')]



for x in comparingDict('afla','afum'):
    out.write('afla\tafum\t%s\n'% ('\t'.join(x)))  
for x in comparingDict('afla','anid'):
    out.write('afla\tanid\t%s\n'% ('\t'.join(x)))
for x in comparingDict('afla','anig'):
    out.write('afla\tanig\t%s\n'% ('\t'.join(x)))
for x in comparingDict('afum','anid'):
    out.write('afum\tanid\t%s\n'% ('\t'.join(x)))
for x in comparingDict('afum','anig'):
    out.write('afum\tanig\t%s\n'% ('\t'.join(x)))
for x in comparingDict('anid','anig'):     
    out.write('anid\tanig\t%s\n'% ('\t'.join(x)))    

    
    
for x in comparingDict('afum','afla'):
    out.write('afum\tafla\t%s\n'% ('\t'.join(x)))  
for x in comparingDict('anid','afla'):
    out.write('anid\tafla\t%s\n'% ('\t'.join(x)))
for x in comparingDict('anig','afla'):
    out.write('anig\tafla\t%s\n'% ('\t'.join(x)))
for x in comparingDict('anid','afum'):
    out.write('anid\tafum\t%s\n'% ('\t'.join(x)))
for x in comparingDict('anig','afum'):
    out.write('anig\tafum\t%s\n'% ('\t'.join(x)))
for x in comparingDict('anig','anid'):     
    out.write('anig\tanid\t%s\n'% ('\t'.join(x)))    

    
