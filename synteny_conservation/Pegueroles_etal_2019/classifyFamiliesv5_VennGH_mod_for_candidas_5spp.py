import sys,os,re

def renameGenes(in1):#ex: 4spv4_6cluster3min.out
    #renames genes adding the species; ex: XLOC_010119 -> XLOC_010119anid
    #creates a dictionary with sp1 (key) pointing to sp2 (val)
    dictTemp={}
    for row in open(in1, 'r').readlines():
        row = row.rstrip().split('\t')
        key= row[2]#+row[0]
        val=[]
        if not key in dictTemp:
            dictTemp[key]=val
        dictTemp[key].append(row[2]) #+row[0])
        dictTemp[key].append(row[3]) #+row[1])
    print(dictTemp)
    return dictTemp
    
##################################################
afla= sys.argv[1] #number of expressed (protein coding/lncRNA) genes
afum= sys.argv[2] #number of expressed (protein coding/lncRNA) genes
anid= sys.argv[3] #number of expressed (protein coding/lncRNA) genes
anig=sys.argv[4]
#cglab= sys.argv[5] #number of expressed (protein coding/lncRNA) genes
in1= sys.argv[5] #out file from synteny_nematodesv4GH.py; ex: 4spv4_6cluster3min.out
out= open(sys.argv[6], 'w') #.fam
out2= open(sys.argv[7], 'w') #.txt
out3= open(sys.argv[8], 'w') #number of families
out4= open(sys.argv[9], 'w') #number of genes

#to group into families
dictTemp1= renameGenes(in1)

dictTemp2= dictTemp1
print(dictTemp1)
for key1, val1 in dictTemp1.iteritems():
    for key2, val2 in dictTemp2.iteritems():
        if key1 in val2:
            val1=list(set(val1+val2))
            dictTemp1[key1]=val1

dictTemp2= dictTemp1
for key1, val1 in dictTemp1.iteritems():
    for key2, val2 in dictTemp2.iteritems():
        if key1 in val2:
            val1=list(set(val1+val2))
            dictTemp1[key1]=val1
           
mylist=[val for val in dictTemp1.values()]
mylistSorted=[]
for x in mylist:
    mylistSorted.append(sorted(x)) 
uniq_mylist = [list(t) for t in set(sorted(map(tuple, mylistSorted)))]
     
print(uniq_mylist)
dictFam={}
#dictFam structure: fam -> {sp1:[,], sp2:[,], sp3:[,], sp4:[,]}
i=1
for x in uniq_mylist:
    key=i
    val=x
    dictFam[key]=val
    i=i+1

aflam=0; afumm=0; anidm=0; anigm=0;  #number of classified genes
for key, val in dictFam.iteritems():
    for x in val:
        out.write("fam%s\t%s\n" % (key, x))
        if re.search("afla", x):
            aflam=aflam+1
        if re.search("afum", x):
            afumm=afumm+1
        if re.search("anid", x):
            anidm=anidm+1
        if re.search("anig", x):
            anigm=anigm+1

#to count the number of species for gene family        
flafum=0; flanid=0;flanig=0; fumnid=0; fumnig=0; nidnig=0
flafumnid=0; flafumnig=0; flanidnig=0;fumnidnig=0
flafumnidnig=0


genesflafum=0; genesflanid=0; genesflanig=0; genesfumnid=0; genesfumnig=0; genesnidnig=0
genesflafumnid=0; genesflafumnig=0; genesflanidnig=0; genesfumnidnig=0
genesflafumnidnig=0

flaF=0; nidF=0; fumF=0; nigF=0
flaG=0; nidG=0; fumG=0; nigG=0

for key, val in dictFam.iteritems():
    mySp=[]
    for x in val:
        x=x.split('|a')
        sp =x[1]
        mySp.append(sp)
    out2.write("fam%s\t%s\n" % (key,'\t'.join(list(set(mySp)))))
    if 'fla' in mySp and 'nid' in mySp and 'fum' in mySp and 'nig' in mySp:
        flafumnidnig +=1
        genesflafumnidnig = genesflafumnidnig+len(mySp)
    if 'fla' in mySp and 'nid' in mySp and 'fum' in mySp and not 'nig' in mySp:
        flafumnid +=1
        genesflafumnid =genesflafumnid+len(mySp)
    if 'fla' in mySp and 'nig' in mySp and 'fum' in mySp and not 'nid' in mySp:
        flafumnig +=1
        genesflafumnig = genesflafumnig+len(mySp)
    if 'fla' in mySp and 'nig' in mySp and 'nid' in mySp and not 'fum' in mySp:
        flanidnig +=1
        genesflanidnig = genesflanidnig+len(mySp)
    if 'fum' in mySp and 'nig' in mySp and 'nid' in mySp and not 'fla' in mySp:
        fumnidnig +=1
        genesfumnidnig = genesfumnidnig+len(mySp)
    if 'fla' in mySp and 'fum' in mySp and not 'nig' in mySp and not 'nid' in mySp :
        flafum +=1
        genesflafum = genesflafum+len(mySp)
    if 'fla' in mySp and 'nid' in mySp and not 'fum' in mySp and not 'nig' in mySp:
        flanid +=1
        genesflanid = genesflanid+len(mySp)
    if 'fla' in mySp and 'nig' in mySp and not 'nid' in mySp and not 'fum' in mySp:
        flanig +=1
        genesflanig = genesflanig+len(mySp)
    if 'fum' in mySp and 'nid' in mySp and not 'fla' in mySp and not 'nig' in mySp:
        fumnid +=1
        genesfumnid= genesfumnid+len(mySp)
    if 'fum' in mySp and 'nig' in mySp and not 'nid' in mySp and not 'fla' in mySp:
        fumnig +=1
        genesfumnig = genesfumnig+len(mySp)
    if 'nig' in mySp and 'nid' in mySp and not 'fla' in mySp and not 'fum' in mySp:
        nidnig +=1
        genesnidnig = genesnidnig+len(mySp)
    if 'fla' in mySp:
        flaF +=1; flaG += len(mySp)
    if 'nid' in mySp:
        nidF +=1; nidG += len(mySp)
    if 'fum' in mySp:
        fumF +=1; fumG += len(mySp)
    if 'nig' in mySp:
        nigF +=1; nigG += len(mySp)
        
print 'flafumnid',flafumnid
print 'flafumnig',flafumnig
print 'flanidnig', flanidnig
print 'fumnidnig',fumnidnig
print 'flafum',flafum
print 'flanid',flanid
print 'flanig',flanig
print 'fumnid',fumnid 
print 'fumnig',fumnig
print 'nidnig',nidnig
print 'flafumnidnig',flafumnidnig

print afla, aflam, flaF
print afum, afumm, fumF
out3.write('rm(list=ls())\n')
out3.write('library(VennDiagram)\n')
out3.write('setwd("%s")\n' % os.getcwd())
out3.write('fla <-%s\n' % (int(afla)-(aflam -flaF)))
out3.write('fum <-%s\n' % (int(afum)-(afumm-fumF)))
out3.write('nid <-%s\n' % (int(anid)-(anidm-nidF)))
out3.write('nig <-%s\n' % (int(anig)-(anigm-nigF)))
out3.write('overlapflafum <- %s\n' % (flafum + flafumnid + flafumnig + flafumnidnig))
out3.write('overlapflanid <- %s\n' % (flanid + flanidnig + flafumnid + flafumnidnig))
out3.write('overlapflanig <- %s\n' % (flanig + flafumnig + flanidnig + flafumnidnig))
out3.write('overlapfumnid <- %s\n' % (fumnid + flafumnid + fumnidnig + flafumnidnig))
out3.write('overlapfumnig <- %s\n' % (fumnig + flafumnig + fumnidnig + flafumnidnig))
out3.write('overlapnidnig <- %s\n' % (nidnig + flanidnig + fumnidnig + flafumnidnig))
out3.write('overlapflafumnid <- %s\n' % (flafumnid + flafumnidnig))
out3.write('overlapflafumnig <- %s\n' % (flafumnig + flafumnidnig))
out3.write('overlapflanidnig <- %s\n' % (flanidnig + flafumnidnig))
out3.write('overlapfumnidnig <- %s\n' % (fumnidnig + flafumnidnig))
out3.write('overlapflafumnidnig <- %s\n' % (flafumnidnig))
out3.write('draw.quad.venn(area1=fla, area2=fum, area3=nid, area4=nig,\n n12=overlapflafum, n13=overlapflanid, n14=overlapflanig,\n n23=overlapfumnid, n24=overlapfumnig, n34=overlapnidnig,\n n123=overlapflafumnid, n124=overlapflafumnig, n134=overlapflanidnig,\n n234=overlapfumnidnig, n1234=overlapflafumnidnig,\n fill = c("skyblue", "pink1", "mediumorchid", "orange"))\n')

print afla,aflam,flaG
print afum,afumm,fumG
out4.write('rm(list=ls())\n')
out4.write('library(VennDiagram)\n')
out4.write('setwd("%s")\n' % os.getcwd())
out4.write('fla <-%s\n' % (int(afla)-(aflam -flaG)))
out4.write('fum <-%s\n' % (int(afum)-(afumm -fumG)))
out4.write('nid <-%s\n' % (int(anid)-(anidm -nidG)))
out4.write('nig <-%s\n' % (int(anig)-(anigm -nigG)))
out4.write('overlapflafum <- %s\n' % (genesflafum + genesflafumnid + genesflafumnig + genesflafumnidnig))
out4.write('overlapflanid <- %s\n' % (genesflanid + genesflanidnig + genesflafumnid + genesflafumnidnig))
out4.write('overlapflanig <- %s\n' % (genesflanig + genesflafumnig + genesflanidnig + genesflafumnidnig))
out4.write('overlapfumnid <- %s\n' % (genesfumnid + genesflafumnid + genesfumnidnig + genesflafumnidnig))
out4.write('overlapfumnig <- %s\n' % (genesfumnig + genesflafumnig + genesfumnidnig + genesflafumnidnig))
out4.write('overlapnidnig <- %s\n' % (genesnidnig + genesflanidnig + genesfumnidnig + genesflafumnidnig))
out4.write('overlapflafumnid <- %s\n' % (genesflafumnid + genesflafumnidnig))
out4.write('overlapflafumnig <- %s\n' % (genesflafumnig + genesflafumnidnig))
out4.write('overlapflanidnig <- %s\n' % (genesflanidnig + genesflafumnidnig))
out4.write('overlapfumnidnig <- %s\n' % (genesfumnidnig + genesflafumnidnig))
out4.write('overlapflafumnidnig <- %s\n' % (genesflafumnidnig))
out4.write('draw.quad.venn(area1=fla, area2=fum, area3=nid, area4=nig,\n n12=overlapflafum, n13=overlapflanid, n14=overlapflanig,\n n23=overlapfumnid, n24=overlapfumnig, n34=overlapnidnig,\n n123=overlapflafumnid, n124=overlapflafumnig, n134=overlapflanidnig,\n n234=overlapfumnidnig, n1234=overlapflafumnidnig,\n fill = c("skyblue", "pink1", "mediumorchid", "orange"))\n')