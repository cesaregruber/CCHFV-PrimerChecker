

from Bio import SeqIO

#from Bio.Seq import Seq
#from Bio.Alphabet import IUPAC
#from Bio.SeqUtils import nt_search

from os import listdir  #, chdir
from os.path import join #,  isfile #, isdir,


from xlwt import Workbook
import xlwt

import re
import itertools as iter

#from numpy import mean, median

# funzioni scritte da me :
import TrovaMatch 
import IUPACfunc
import CodeMatch
import PrintMatch







#~ MAIN ~ #

# ######## INPUTS : ###############################################################

#path generico di input:
path=".."   # ci vado a cercare ....


#dizionari delle sequenze allineate:
tuttef=path+"/inputfiles/CCHFVsequencesAligned.fasta"
newworldvirus=SeqIO.to_dict(SeqIO.parse(tuttef, "fasta"))

# AGGIUSTAMENTO PER DEYDE  e gli altri :
#DEYDE: parte alla pos  4 nt  e finisce alla pos -1
# ATKINSON : parte alla pos 4 nt
# KAMBOJ: parte alla pos 35 nt

worldvirus= {}
for k in newworldvirus :
    myseq = newworldvirus[k].seq 
    #if myseq[:35] != '-'*35 :
    worldvirus[k] = newworldvirus[k]



#dizionari delle sequenze NON ALLINEATE:
tuttefna=path+"/inputfiles/CCHFVsequences.fasta"
worldvirusNOal=SeqIO.to_dict(SeqIO.parse(tuttefna, "fasta"))
worldvirusNOaln={}
for k in worldvirusNOal :
    if k in worldvirus :
        worldvirusNOaln[k] = worldvirusNOal[k]



#sequenza di riferimento NCBI :
myref = path+"/inputfiles/NC_005302_S_Aligned.fasta"
refvirus = SeqIO.to_dict(SeqIO.parse(myref, "fasta"))



#  #1:  ... la lista dei files con gli headers dei clades : 
#   !!!! ATTENZIONE: questa poi la uso per aprire e leggere tutti i file, utilizzando  join(path,clade) !!!!!!!!! 
clades=[f for f in listdir(path+"/inputfiles") if 'clade' in f and 'txt' in f and '~' not in f and 'clades' not in f and "SENZA_SPAGNA" not in f]
# poi la uso per :
    #2 : dizionario dei dizionari: sequenze allineate divise per clades
    #3: lista degli headers dei clades

#dizionario dei primers:
primf = open("../inputfiles/Primers_table.csv", 'r')

prim = path+"/Primers_sequences.fa"  # nome dell' OUTPUT 'primfa'


# ######## OUTPUTS : ###############################################################

#dizionario dei primers:
primfa = open(prim, 'w')    

# tabella dei risultati :
outf = open(join(path,'PrimersMetches.csv'), 'w')                  #schema dei primer

outfasta=open(join(path,  'Primers_Sequences_aligned.fasta'), 'w')    # FASTA dei primer (per controllo)


myworkbook = join(path,  'PrimersMetches_RESULTS.xls')            #file exel dei primer FORMATTATO !!!


# ##############################################################################


#info sulle codifiche di ambiguita' :
ambidic ,  ambinonambidic ,  invambidic ,  complambidic ,  ambikeys  ,  ambiset  = IUPACfunc.iupacdic()



#dizionari delle sequenze allineate:
wordvirus = {}
for k in worldvirus :
    newk = k.split('|')[0]
    wordvirus[newk] = worldvirus[k]
    wordvirus[newk].seq = worldvirus[k].seq.upper()
    

#dizionari delle sequenze NON ALLINEATE:
wordvirusNOaln = {}
for k in worldvirusNOaln :
    newk = k.split('|')[0]
    wordvirusNOaln[newk] = worldvirusNOaln[k]
    wordvirusNOaln[newk].seq = worldvirusNOaln[k].seq.upper()

#dizionari dei gap aggiunti nell'allineamento:
wordvirusgaps = {}
for vir in wordvirus :
    alns = str(wordvirus[vir].seq)
    gapoints = [i for i, letter in enumerate(alns) if letter == '-']
    wordvirusgaps[vir] = gapoints

# sequenza di riferimento NCBI ( con e senza gap !!)
RefSeq = ''
RefSeqNoGaps = ''
refname = ''
for vir in refvirus :
    if 'NC_' in vir :
        RefSeq = str(refvirus[vir].seq)
        RefSeqNoGaps = str(refvirus[vir].seq).replace('-', '')
        refname = vir



#dizionario dei dizionari: sequenze allineate divise per clades

#clades: lista dei files con gli headers dei clades

cladevirus = {}
cladegaps = {}
cladeNOaln = {}


for clade in clades :
    k= clade[:-4]
    cladevirus[k] = {} 
    cladegaps[k] = {}
    cladeNOaln[k] = {}
    
    mycladef = open(    join(path+'/inputfiles', clade), 'r'    )
    mycladelines = mycladef.readlines()
    mycladeheads = [l.split('\n')[0]  for l in mycladelines if len(l)> 2 ]
    
    for myh in mycladeheads :
        for wh in wordvirus.keys() :
            if myh in wh :
                cladevirus[k][myh] = wordvirus[wh]
                cladegaps[k][myh] = wordvirusgaps[wh]
                cladeNOaln[k][myh] = wordvirusNOaln[wh]

cladeviruslist =   sorted(cladevirus.keys())

Nclades = len(cladeviruslist)




#dizionario dei primers:

primerstype = {} # dizionario di dizionari PER TROVARE LE SEQUENZE IDENTICHE alle singole sonde

primersgoals = {} # dizionario di dizionari PER REGISTRARE IL NUMERO DI MATCH E MISMATCH delle sequenze con le singole sonde

for l in primf :
    line = l.replace('"', '')
    if line[:6] != "AUTHOR" :
        ls = line.split(',')
        name = '|'.join( [ls[0],  ls[1], ls[2],  ls[5] ,  ls[4] ] )
        primfa.write('>' + name + '\n' )
        seq = ls[6]#[:-1]
        #print seq
        primfa.write(seq+ '\n' )
        
        primerstype[ls[0]] = { 'type': ls[7]  }


primfa.write('\n' )

primf.close()
primfa.close()



primersdic=SeqIO.to_dict(SeqIO.parse(prim, "fasta"))

primerskeys = []
for record in SeqIO.parse(prim, "fasta"):
    primerskeys.append(record.id)



    

#schema dei primer
# FASTA dei primer (per controllo)

book = Workbook(style_compression=2)                                                            #file exel dei primer FORMATTATO !!!
sheet1 = book.add_sheet('Sheet 1')
#larghezza prima colonna:newprimer
numchars=0
lenprime = 0
for p in primerskeys :
    l = len(str(primersdic[p].seq))
    if len(p) > numchars :
        numchars = len(p)
    if l > lenprime :
        lenprime = l


# larghezza delle colonne:
sheet1.col(0).width = 256 * (numchars + 1)
sheet1.col(1).width = 256 * 7
sheet1.col(2).width = 256 * 5
#sheet1.col(3).width = 256 * 7

for  co in range(3, lenprime+200) :
    sheet1.col(co).width = 256*3


#stile della tabella:
style = xlwt.easyxf('font: bold 1, color black; ')# pattern: pattern solid;')
    


strainordercontrol = 0      #tutti gli strain sono nello stesso ordine (in tutti i clade )




#cerco le sequenze dei primers nei clades:


#p = primersdic['Yapa_et_al.2005|CCRealP1|Primer_S']

lunghezzatotalesequenze = len(str(wordvirus[wordvirus.keys()[0]].seq))

book_col = 0
book_row = 0

strainlistINclades={}

DELTA = 10     # quanta sequenza stampo prima e dopo il primer  ################
LASTIMPORTANT = 5 # numero di nucleotidi importanti alla fine del primer, per non avere problemi con la trascrittasi

primersorder = []

for p in primerskeys :
    
    primer = str(primersdic[p].seq)
    pname = p.split('|')[0]
    primersorder.append(pname)
    primerstype[pname][p] = {}
    revertito = 0
    outf.write('\n\n'+p+ '\t')
    
    book_row += 2
    book_col = 0
    sheet1.write(book_row, book_col, p,style)
    book_col +=1
    
    [ TOTALESATTI ,  init,  fine , primergaps,  trueinit ,  truefine,  gappedprimer,  qualiesatti]  = TrovaMatch.trova_match(primer, wordvirus,  wordvirusgaps,  wordvirusNOaln, RefSeq,  RefSeqNoGaps  )
    
    
    if TOTALESATTI == 0 :
        primer = str(primersdic[p].seq)
        # provo a vedere se e' complement- reverse :
        primer = IUPACfunc.revcomp(primer, complambidic)
        
        [ TOTALESATTI ,  init,  fine , primergaps,  trueinit ,  truefine,  gappedprimer,  qualiesatti] = TrovaMatch.trova_match(primer, wordvirus,  wordvirusgaps,  wordvirusNOaln, RefSeq,  RefSeqNoGaps  )
        
        if TOTALESATTI != 0 :
            print "primer trovato: era revertito!!!!!"
            revertito = 1


    print  p
    print primer
    print   TOTALESATTI ,  init,  fine , primergaps,  trueinit ,  truefine
    print gappedprimer
    print "# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------  #"



    if TOTALESATTI != 0 or revertito == 1 :
        
        primerstype[pname][p]['all'] = qualiesatti
        
        firstORlast = ''
        
        if 'primer_rev' in p :
            firstORlast = 'first'
        elif 'primer_for' in p :
            firstORlast = 'last'
        elif 'probe' in p :
            firstORlast = 'none'
        else :
            print "------------------------------"
            print "ERROR IN FIRST or LAST"
            print p
            print "------------------------------"
        
        [newtarget ,  newcounts , expandedtarget , dstart, dend, strainlist ,  numvarlist, numvarlast] =  CodeMatch.codifica_match( gappedprimer, wordvirus, init, fine ,  ambidic,  invambidic,  DELTA, firstORlast , LASTIMPORTANT )    # codifica_match( gappedprimer, myclade, init, fine ,  ambidic)
        
        
        primerstype[pname][p]['all_MismatchCounts'] = numvarlist
        
        primerstype[pname][p]['all_StrainsOrder'] = strainlist
        
        primerstype[pname][p]['all_NOTMATCHlastFIVE'] = numvarlast
        
        
        #controllo dell'ordine degli strains con il primer precedente:
        if "all" in strainlistINclades :
            print "-----------------control strain order---------------------------------"
            if strainlistINclades["all"] == strainlist :
                print  "strain order OK"
            else :
                print "!!!!!!!!!!!!!!!! ERROR  ' !!!! problem in strain order !!!!!!!!!!!!!!!"
                strainordercontrol += 1        # <- conto il numero degli errori nell'ordine degli strain
        else :
            strainlistINclades["all"] = strainlist
        
        # COMINCIO A SCRIVERE :
        outfasta.write('>'+p+'\n'+'-'*init + gappedprimer+'-'*(lunghezzatotalesequenze-fine) +'\n')
        
        outf.write(str(init) + '\t')
        outf.write(str(fine) +'\tprimer\t')
        outf.write( '\t'.join([c for c in gappedprimer]) + '\n')
        
        #sheet1.write(book_row, book_col, 'primer', style)
        book_col += dstart+1
        
        style = xlwt.easyxf('font: bold 1, color black; alignment: horizontal right;  ')# pattern: pattern solid;')
        sheet1.write(book_row, book_col, str(trueinit), style)
        style = xlwt.easyxf('font: bold 1, color black; ')# pattern: pattern solid;')
        book_col += 1
    
        oldpattern = xlwt.Pattern()
        
        for c in gappedprimer :
            style.pattern.pattern_fore_colour = 7
            style = xlwt.easyxf('pattern: pattern solid; font: bold 1, color black; alignment: horizontal center; ')
            sheet1.write(book_row, book_col, c, style)
            style.pattern.pattern_fore_colour = 7
            book_col += 1
        
        style.pattern = oldpattern
        style = xlwt.easyxf('font: bold 1, color black; ')
        
        #book_col += 2
        style = xlwt.easyxf('font: bold 1, color black; alignment: horizontal left;  ')# pattern: pattern solid;')
        sheet1.write(book_row, book_col, str(truefine), style)
        style = xlwt.easyxf('font: bold 1, color black; ')# pattern: pattern solid;')
        book_col += 1
        
        book_row +=1
    
        freqmatches = round( float(TOTALESATTI)/ len(wordvirus.keys()) ,  2)            #frequenza di volte che il primer "meccia" esattamente con le sequenze target
        outf.write('\t\t\tworld\t'+  '\t'.join([c for c in newtarget ]) + '\n')
    
        #ripeto il gioco per tutti i clades :
        for clade in cladeviruslist :
            
            TOTALESATTI = 0
            TOTALRILEVABILI = 0
            qualiesatti = []
            qualirilevabili = []           # rilevabile se e' uguale negli ultimi 15 nucleotidi e nei restanti e' almeno uguale nei due terzi (come da UCSC in-silico PCR e da: 
                                                    #Yu B.; Zhang C. (2011). "In Silico PCR Analysis". In Silico Tools for Gene Discovery. Methods in Molecular Biology. 760. pp. 91-107. DOI:10.1007/978-1-61779-176-5_6     
                                                    # NON COMPLETATO !!!!
            
            myclade = cladevirus[clade]
            mygaps = cladegaps[clade]
            myNOaln = cladeNOaln[clade]
            
            # TRADUCO GLI AMBIGUITY CODES  
            s = str(gappedprimer)         # ricopio il primer com' e' 
            for a in ambikeys:
                if a in s :
                    s = s.replace(a,  '[' +ambidic[a] + ']' )            # <- 
            
            # QUALI ESATTI :            ####################################
#            myvir = []  # [nome, sequenza]
#            myseq = ''
#            for vir in myclade :
#                myseq =  str(myclade[vir].seq)
#                if re.search( s,  myseq ):     # <- "s" e' la sequenza del primer tradotta con gli ambiguity codes
#                    TOTALESATTI += 1
#                    qualiesatti.append(vir)

            # QUALI RILEVABILI  : #########################################   NON l'HO FINITO !!!!
            #myseq = ''
            #for vir in myNOaln :
            #    myseq =  str(myNOaln[vir].seq)[init:fine]
            #    if re.search( s,  myseq ):     # <- "s" e' la sequenza del primer tradotta con gli ambiguity codes
            #        TOTALRILEVABILI += 1
            #        qualirilevabili.append(vir)
            
            # primerstype[pname][p][clade] = qualiesatti   <---vecchio !!
            
            [newtarget ,  newcounts , expandedtarget , dstart, dend, strainlist ,  numvarlist, numvarlast]  =  CodeMatch.codifica_match( gappedprimer, myclade, init, fine ,  ambidic,  invambidic,  DELTA, firstORlast , LASTIMPORTANT) 
            
            primerstype[pname][p]['MismatchCounts_'+ clade] = numvarlist
        
            primerstype[pname][p]['StrainsOrder_'    + clade] = strainlist
            
            primerstype[pname][p]['NOTMATCHlastFIVE_'    + clade] = numvarlast
            
            qualiesatti = []
            for i in range(len(strainlist)) :
                if numvarlist[i]  == 0 :
                    qualiesatti.append(strainlist[i])
            
            primerstype[pname][p][clade] = qualiesatti 
            
            
            
            #controllo del nuovo ordine degli strain rispetto al precedente assay
            if clade in strainlistINclades :
                print "-----------------control strain order---------------------------------"
                if strainlistINclades[clade] == strainlist :
                    print  "strain order OK"
                else :
                    print "!!!!!!!!!!!!!!!! ERROR  ' !!!! problem in strain order !!!!!!!!!!!!!!!"
                    strainordercontrol += 1        # <- conto il numero degli errori nell'ordine degli strain
            else :
                strainlistINclades[clade] = strainlist
            
            
            outf.write('\t\t\t'+clade+ '\t'+  '\t'.join([c for c in newtarget]) + '\n')
            
            book_col = 1
            sheet1.write(book_row, book_col, clade, style)
            book_col += 2
            
            [style, book_row, book_col ] =  PrintMatch.stampa_match(sheet1, book_row, book_col , style , newtarget ,   newcounts , expandedtarget, dstart, dend )
            print newtarget
            print newcounts
        
            
            #freqmatches = round( float(TOTALESATTI)/ len(myclade.keys()) ,  2)           #frequenza di volte che il primer matcha esattamente con le sequenze target
            
            freqmatches = sum(newcounts)/len(newcounts)             #frequenza di volte che i nucleotidi del primer matchano esattamente con le sequenze target
    
            book_col = 120
            sheet1.write(book_row, book_col, str(freqmatches), style)
            book_col += 1
            

            num_NOTMATCH_lastFIVE = sum(numvarlast)            #NUMERO di volte che i nucleotidi del primer NON matcha esattamente con le sequenze target NEGLI ULTIMI 5 NUCLEOTIDI
        
            book_col += 3
            sheet1.write(book_row, book_col, str(num_NOTMATCH_lastFIVE), style)
            book_col += 1                        

            
            book_row +=1
        




primersorder = list(set(primersorder))

all_dgenset = {}  #dizionaro degli assay, con dentro il dizionario dei clade, con il nome delle sequenze "perfectly marching" ( match al 100% )

for pname in primersorder : #primerstype :
    
    all_dgenset[pname] = {}
    
    book_col  = 0
    book_row += 4
            
    sheet1.write(book_row, book_col, pname, style)
    book_col += 1
    sheet1.write(book_row, book_col, primerstype[pname]['type'] , style)
    book_row += 1
    
    #        '###########################################################################################'
    # MODIFICATO A OTTOBRE: STAMPO SOLO GLI ASSAY CON PRIMER MULTIPLI !!!!
    #print '###########################################################################################'
    #print pname ,  primerstype[pname]['type'] 
    
    # CASO SEMPLICE: NON CI SONO SONDE MULTIPLE (c'e' da considerare un elemento nel dizionario in piu', che corrisponde alla voce 'type' )
    if ( ( primerstype[pname]['type'] == 'Nested' and  len(primerstype[pname].keys()) <= 5 ) 
           or (  primerstype[pname]['type'] == 'PCR' and  len(primerstype[pname].keys()) <= 3 ) 
           or  (  primerstype[pname]['type'] == 'RealTime' and  len(primerstype[pname].keys()) <= 4 )  
           or (  primerstype[pname]['type'] == 'LAMP' and  len(primerstype[pname].keys()) <= 9 )               ):
        
        Nseq = {}
        dgenset = {}
        numismatch = {}    #numero dei mismatch fra primer e clade
        totmismatch = {} #somma dei mismatch di TUTTE le sonde di un assay 
        numLASTmism = {} # numero dei mismatch negli ultimi LASTIMPORTANT nucleotidi 
        totLASTmismatch = {} #somma dei mismatch di TUTTE le sonde di un assay negli ultimi LASTIMPORTANT nucleotidi 
        
        for clade in cladeviruslist :
            dgenset[clade] = set(strainlistINclades[clade])    # set dei genomi appartenenti al clade  ------> sara' il set dei genomi esatti !!!!!
            Nseq[clade] = len( strainlistINclades[clade])      # numero dei genomi appartenenti al clade
            numismatch[clade] = 0                                            #numero dei mismatch fra primer e clade
            totmismatch[clade] = [0]*Nseq[clade]          #somma dei mismatch di TUTTE le sonde di un assay 
            numLASTmism[clade] = 0     #numero dei mismatch fra primer e clade  negli ultimi LASTIMPORTANT nucleotidi 
            totLASTmismatch[clade] = [0]*Nseq[clade]          #somma dei mismatch di TUTTE le sonde di un assay negli ultimi LASTIMPORTANT nucleotidi 
            
        
        for primer in primerstype[pname] :
            if primer != 'type' :
                for clade in primerstype[pname][primer] :
                    if clade[:5] == 'clade' :         # se e' effetteivamente un clade
                        dgenset[clade] = dgenset[clade] & set( primerstype[pname][primer][clade] )
                    elif clade[:15] == 'MismatchCounts_' :        #se invece e' il conteggio dei mismatches :
                        numismatch[ clade[15:] ] = max(  [numismatch[ clade[15:]] ]  +primerstype[pname][primer][clade]  )       # appendo tutti i nuovi valori di mismatch dei singoli primer/probe, e prendo solo il valore massimo 
                        
                        totmismatch[ clade[15:] ] = [ totmismatch[ clade[15:] ][i] + primerstype[pname][primer][clade][i] for i in range(Nseq[ clade[15:] ]) ]      # sommo i mismatch degli assays, sequenza per sequenza
                    
                    elif clade[:17] == 'NOTMATCHlastFIVE_' :        #se invece e'  il numero di mismatch negli ultimi LASTIMPORTANT nucleotidi 
                        numLASTmism[ clade[17:] ] =  max(  [numLASTmism[ clade[17:]] ]  +primerstype[pname][primer][clade]  )              # appendo tutti i nuovi valori di mismatch dei singoli primer/probe, e prendo solo il valore massimo 
                        
                        totLASTmismatch[ clade[17:] ] = [ totLASTmismatch[ clade[17:] ][i] + primerstype[pname][primer][clade][i] for i in range(Nseq[ clade[17:] ]) ]      # sommo i mismatch degli assays, sequenza per sequenza
                    
                        
        for clade in cladeviruslist :
            freqmatches = round(len(dgenset[clade])*100.0 / Nseq[clade] , 1)
            
            #print clade ,  freqmatches ,  "massimo numero mismatch =" ,  numismatch[clade] ,  totmismatch[clade]
            
            #book_row += 1
            book_col = 1
            sheet1.write(book_row, book_col, clade, style)
            book_col += 1
            sheet1.write(book_row, book_col, (str(freqmatches) + '%' ), style)
            book_col += 1
            sheet1.write(book_row, book_col, ( '('+ str(Nseq[clade] )+')'  ), style)
            book_col += 1
            sheet1.write(book_row, book_col,   str(numismatch[clade] ) , style)        # stampo  il numero massimo di mismatch trovati per ciascuna sequenza, per ciascun primer/probe
            book_col += 1
            sheet1.write(book_row, book_col,   str(max(totmismatch[clade] )) , style)     # stampo solo il numero massimo di mismatch trovati per ciascuna sequenza
            book_col += 1
            sheet1.write(book_row, book_col,   str(numLASTmism[clade] ) , style)        # stampo  il numero massimo di mismatch negli ultimi LASTIMPORTANT nucleotidi trovati per ciascuna sequenza, per ciascun primer/probe
            book_col += 1
            sheet1.write(book_row, book_col,   str(max(totLASTmismatch[clade] )) , style)     # stampo solo il numero massimo di mismatch negli ultimi LASTIMPORTANT nucleotidi trovati per ciascuna sequenza

            book_row +=1

            
            all_dgenset[pname][clade] = {}
            all_dgenset[pname][clade]["SetGeniEsatti"] = dgenset[clade]
            all_dgenset[pname][clade]["MismatchPerSequenza"] = totmismatch[clade]            
            all_dgenset[pname][clade]["MismatchUltimiNucleotidi"] = totLASTmismatch[clade]
        
#        # SALVO UN CONTROLLO CON UN ASSAY A CASO
#        if pname == 'Koehler_er_al.2018' :
#            ko_numismatch = numismatch
#            ko_totmismatch = totmismatch
#            ko_dgenset = dgenset
    
    # CASO COMPLICATO :  CI SONO SONDE MULTIPLE !!! 
    else  :
        
        print '###########################################################################################'
        print pname ,  primerstype[pname]['type'] 
        
        allfors = {}
        allrevs = {}
        allprobes = {}
        
        allformism = {}         #  per ogni forward, il valore di TUTTI i mismatch con le sequenze di quel clade
        allrevmism = {}         # per ogni reverse, il valore di TUTTI i mismatch con le sequenze di quel clade
        allprobemism = {}           #per ogni probe, il valore di TUTTI i mismatch con le sequenze di quel clade
        
        allforLASTmism = {}         #  per ogni forward, il valore di TUTTI i mismatch con le sequenze di quel clade negli ultimi LASTIMPORTANT nucleotidi 
        allrevLASTmism = {}         # per ogni reverse, il valore di TUTTI i mismatch con le sequenze di quel clade negli ultimi LASTIMPORTANT nucleotidi 
        allprobeLASTmism = {}           #per ogni probe, il valore di TUTTI i mismatch con le sequenze di quel clade negli ultimi LASTIMPORTANT nucleotidi 
        
        allcombmism  = {}           #per ogni combinazione (forward+probe+reverse), il valore di  TUTTI i mismatch con le sequenze di quel clade
        allcombLASTmism  = {}           #per ogni combinazione (forward+probe+reverse), il valore di  TUTTI i mismatch  negli ultimi LASTIMPORTANT nucleotidi con le sequenze di quel clade
        
        
        # NOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO    totcombmism  = {}       #per ogni combinazione (forward+probe+reverse), il valore di TUTTI i mismatch con le sequenze di quel clade
        
        numismatch = {}    #numero dei mismatch fra primer e clade
        totmismatch = {} #somma dei mismatch di TUTTE le sonde di un assay
                
        numLASTmismatch = {}    #numero dei mismatch fra primer e clade negli ultimi LASTIMPORTANT nucleotidi
        totLASTmismatch = {} #somma dei mismatch di TUTTE le sonde di un assay negli ultimi LASTIMPORTANT nucleotidi
        
        Nseq = {}
        dgenset = {}
        for clade in cladeviruslist :
            dgenset[clade] = set(cladevirus[clade].keys())
            Nseq[clade] = len( cladevirus[clade].keys() )
            
            allfors[clade] = {}
            allrevs[clade] = {}
            allprobes[clade] = {}
            
            allformism[clade] = {}                    #  per ogni forward, il valore di TUTTI i mismatch con le sequenze di quel clade
            allrevmism[clade] = {}                  # per ogni reverse, il valore di TUTTI i mismatch con le sequenze di quel clade
            allprobemism[clade] = {}                #per ogni probe, il valore di TUTTI i mismatch con le sequenze di quel clade
            
            allcombmism[clade] = {}              #per ogni combinazione (forward+probe+reverse), il valore di TUTTI i mismatch con le sequenze di quel clade
                        
            allforLASTmism[clade] = {}                    #  per ogni forward, il valore di TUTTI i mismatch con le sequenze di quel clade negli ultimi LASTIMPORTANT nucleotidi 
            allrevLASTmism[clade] = {}                  # per ogni reverse, il valore di TUTTI i mismatch con le sequenze di quel clade negli ultimi LASTIMPORTANT nucleotidi 
            allprobeLASTmism[clade] = {}                #per ogni probe, il valore di TUTTI i mismatch con le sequenze di quel clade negli ultimi LASTIMPORTANT nucleotidi 
            
            allcombLASTmism[clade] = {}              #per ogni combinazione (forward+probe+reverse), il valore di TUTTI i mismatch negli ultimi LASTIMPORTANT nucleotidi con le sequenze di quel clade
            
            #NOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO      totcombmism[clade]  = {}       #per ogni combinazione (forward+probe+reverse), il valore di TUTTI i mismatch con le sequenze di quel clade

        
        for primer in primerstype[pname] :
            if primer != 'type' :
                for clade in primerstype[pname][primer] :
                    if clade[:5] == 'clade' :         # se e' effetteivamente un clade
                        if 'primer_for' in primer :
                            allfors[clade][primer] = set(primerstype[pname][primer][clade])          # set delle sequenze identiche ai forward, clade per clade 
                        elif 'primer_rev' in primer :
                            allrevs[clade][primer] = set(primerstype[pname][primer][clade])          # set delle sequenze identiche ai reverse, clade per clade 
                        elif 'probe' in primer :
                            allprobes[clade][primer] = set(primerstype[pname][primer][clade])       # set delle sequenze identiche ai probe, clade per clade 
                        else :
                            print '*************************************************PRIMER NAME ERROR*******************************************************************'
                            print "ERRORE IN:"
                            print primer
                            print '*************************************************************************************************************************************'
                    
                    elif clade[:15] == 'MismatchCounts_' :        #se invece e' il conteggio dei mismatches :
                        if 'primer_for' in primer :
                            allformism[ clade[15:] ][primer] = primerstype[pname][primer][clade]               # prendo, per ogni forward, il valore di TUTTI i mismatch con le sequenze di quel clade
                        elif 'primer_rev' in primer :
                            allrevmism[ clade[15:] ][primer] = primerstype[pname][primer][clade]                # prendo, per ogni reverse, il valore di TUTTI i mismatch con le sequenze di quel clade
                        elif 'probe' in primer :
                            allprobemism[ clade[15:] ][primer] = primerstype[pname][primer][clade]            # prendo, per ogni probe, il valore di TUTTI i mismatch con le sequenze di quel clade
                        else :
                            print '*************************************************PRIMER NAME ERROR*******************************************************************'
                            print "ERRORE IN:"
                            print primer
                            print '*************************************************************************************************************************************'                        
                    
                    elif clade[:17] == 'NOTMATCHlastFIVE_' :        #se invece e' il conteggio dei mismatches  negli ultimi LASTIMPORTANT nucleotidi  :
                        if 'primer_for' in primer :
                            allforLASTmism[ clade[17:] ][primer] = primerstype[pname][primer][clade]               # prendo, per ogni forward, il valore di TUTTI i mismatch con le sequenze di quel clade negli ultimi LASTIMPORTANT nucleotidi
                        elif 'primer_rev' in primer :
                            allrevLASTmism[ clade[17:] ][primer] = primerstype[pname][primer][clade]                # prendo, per ogni reverse, il valore di TUTTI i mismatch con le sequenze di quel clade negli ultimi LASTIMPORTANT nucleotidi
                        elif 'probe' in primer :
                            allprobeLASTmism[ clade[17:] ][primer] = primerstype[pname][primer][clade]            # prendo, per ogni probe, il valore di TUTTI i mismatch con le sequenze di quel clade negli ultimi LASTIMPORTANT nucleotidi
                        else :
                            print '*************************************************PRIMER NAME ERROR*******************************************************************'
                            print "ERRORE IN:"
                            print primer
                            print '*************************************************************************************************************************************'                        
                    
        for clade in cladeviruslist :
            if allprobes[clade] == {}  :   # se non c'e' il probe 
                allprobes[clade]['fintoprobe'] = set(cladevirus[clade].keys())          # faccio finta che ce ne sia uno con tutte le sequenze  possibili
                allprobemism[clade]['fintoprobe'] = [0]*Nseq[clade]
                allprobeLASTmism[clade]['fintoprobe'] = [0]*Nseq[clade]
            
            allcombinations = [] # lista di tutte le combinazioni possibili 
            
            
            for fors in allfors[clade] :
                for revs in allrevs[clade] :
                    for probs in allprobes[clade] :
                        stessigenomi = allfors[clade][fors] & allrevs[clade][revs] & allprobes[clade][probs]        # set dei genomi identici condivisi
                        mycomb = [fors, revs, probs]           # combinazione relativa
                        
                        
                        #per ogni combinazione (forward+probe+reverse), la SOMMA di TUTTI i mismatch con le sequenze di quel clade:
                        allcombmism[clade][fors, revs, probs] =   [ sum(mism) for mism in zip ( allformism[ clade ][fors]    ,   allrevmism[ clade ][revs] ,   allprobemism[ clade ][probs]   )  ]
                        #per ogni combinazione (forward+probe+reverse), la SOMMA di TUTTI i mismatch con le sequenze di quel clade negli ultimi LASTIMPORTANT nucleotidi:
                        allcombLASTmism[clade][fors, revs, probs] =   [ sum(mism) for mism in zip ( allforLASTmism[ clade ][fors]    ,   allrevLASTmism[ clade ][revs] ,   allprobeLASTmism[ clade ][probs]   )  ]
                        
                        
                        #totcombmism[clade][fors, revs, probs]  =  [ totmismatch[ clade[15:] ][i] + primerstype[pname][primer][clade][i] for i in range(Nseq[ clade[15:] ]) ]   # sommo i mismatch degli assays, sequenza per sequenza
                    
                        mincomb = min( allcombmism[clade][fors, revs, probs] )                                                       #numero minimo di mismatch per la combinazione
                        maxcomb = max( allcombmism[clade][fors, revs, probs]  )                                                     #numero massimo di mismatch per la combinazione
                            
                        minLASTcomb = min( allcombLASTmism[clade][fors, revs, probs] )                                                       #numero minimo di mismatch per la combinazione negli ultimi LASTIMPORTANT nucleotidi
                        maxLASTcomb = max( allcombLASTmism[clade][fors, revs, probs]  )                                                     #numero massimo di mismatch per la combinazione negli ultimi LASTIMPORTANT nucleotidi
                        
                        
                        #meancomb = mean( allcombmism[clade][fors, revs, probs] )                                              #numero medio di mismatch per la combinazione
                        #meadiancomb = median( allcombmism[clade][fors, revs, probs] )                                        #valore mediano di mismatch per la combinazione
                        
                        
                        mytuple = (len(stessigenomi),  mycomb ,  stessigenomi,  mincomb,  maxcomb, allcombmism[clade][fors, revs, probs]  ,  minLASTcomb,  maxLASTcomb ,  allcombLASTmism[clade][fors, revs, probs]   )             
                        allcombinations.append(mytuple)
                        
           
            
            mismfors = min( [ max( allformism[clade][p] )   for p in allformism[clade].keys() ] )              #prendo il forward con il numero mininmo di mismatches massimi
            mismrevs = min( [  max( allrevmism[clade][p] )   for p in allrevmism[clade].keys() ] )            #prendo il reverse con il numero mininmo di mismatches massimi
            mismprobes = min( [ max( allprobemism[clade][p] )   for p in allprobemism[clade].keys() ] )            #prendo il probe con il numero mininmo di mismatches massimi
 
            numismatch[clade] = max(mismfors, mismrevs ,  mismprobes)        # il numero di mismatch riportato sara' quello del peggiore fra forward,reverse e probes
            
            
            mismLASTfors = min( [ max( allforLASTmism[clade][p] )   for p in allforLASTmism[clade].keys() ] )              #prendo il forward con il numero mininmo di mismatches massimi
            mismLASTrevs = min( [  max( allrevLASTmism[clade][p] )   for p in allrevLASTmism[clade].keys() ] )            #prendo il reverse con il numero mininmo di mismatches massimi
            mismLASTprobes = min( [ max( allprobeLASTmism[clade][p] )   for p in allprobeLASTmism[clade].keys() ] )            #prendo il probe con il numero mininmo di mismatches massimi
 
            numLASTmismatch[clade] = max(mismLASTfors, mismLASTrevs ,  mismLASTprobes)        # il numero di mismatch riportato sara' quello del peggiore fra forward,reverse e probes
      
            
            # cerco se ci sono combinazioni che beccano al 100% almeno un genoma :
            
            classifica_esatti = sorted(allcombinations)[::-1]
            
            thewinner_exact = classifica_esatti[0]          # la migliore combinazione e' quella che becca almeno una sequenza al 100%
            for esatta in classifica_esatti :
                if esatta[0] == thewinner_exact[0] :      # se la combinazione successiva ha lo stesso numero di beccate 
                    if   (  esatta[7]   < thewinner_exact[7]    )  or   (  esatta[7]   <= thewinner_exact[7]    and  esatta[4]  < thewinner_exact[4]  )   or ( esatta[7]   <= thewinner_exact[7]   and  esatta[4]  <= thewinner_exact[4]  and esatta[3] < thewinner_exact[3] )  :
                        thewinner_exact = esatta                 #  allora prendo quella col minor numero di mismatch massimi per sequenza
                

            
            setesatti = set()                     # IMPORTANTE :  SET COMPLETO DEI GENOMI BECCATI AL 100% DA UNA QUALUNQUE FRA TUTTE LE COMBINAZIONI DI PRIMERS/PROBE
            
            maxmism = 1000          # numero massimo di mismatch che l'assay puo' incontrare (quindi sara' il minimo fra tutti i massimi )
            minmism = 100               # numero minimo di mismatch che l'assay puo' incontrare (quindi sara' il minimo fra tutti i minimi )
            thewinner_mism = 0            # combinazione migliore
            
            mismperseq = [1000 ]*Nseq[clade]        #numero di mismatch per sequenza (scelto come il minimo fra tutte le combinazioni )
            mismLASTperseq = [6]*Nseq[clade]        #numero di mismatch per sequenza (scelto come il minimo fra tutte le combinazioni )
            
            for el in allcombinations :            # cerco in tutte le combinazioni   
                
                for seq in range(Nseq[clade] ) :         # NOTABENE : tratto i mimsatch negli ultimi 5 e i mismatch totali come parametri INDIPENDENTI (semplificazione)          
                    combmism = el[5]
                    mismperseq[seq] = min(mismperseq[seq]  , combmism[seq] )      # per ogni sequenza, trovo la combinazione che la becca meglio !!!
                    
                    combLASTmism = el[8]
                    mismLASTperseq[seq] = min(mismLASTperseq[seq]  , combLASTmism[seq] )      # per ogni sequenza, trovo la combinazione che la becca meglio NEGLI ULTIMI NUCLEOTIDI !!!
                
                if el[0] > 0  :                                           # se ci sono combinazioni che beccano al 100% almeno un genoma 
                    setesatti = setesatti | el[2]               # registro il nome del genoma "beccato"                     <-----  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FURBATA !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
                if el[4]  < maxmism  or (el[4]  <= maxmism  and el[3] < minmism ) :              # se il numero di mismatch massimo e' migliore (cioe' minore ) del precedente registrato --  o il mismatch massimo e' uguale ma il minimo e' migliore:
                    maxmism = el[4]        #registro il nuovo massimo
                    minmism = el[3]
                    thewinner_mism = el          # la migliore combinazione e' quella con il numero minimo di mismatch
                    
                elif (el[4]  == maxmism  and el[3] == minmism ) and (el[0] > thewinner_mism[0] ) :  # oppure se i mismatch sono gli stessi, ma il numero di beccate al 100% e' maggiore:
                    thewinner_mism = el          # la migliore combinazione e' quella con il numero minimo di mismatch , e poi con il numero di 100% maggiore
            
                
            
            totesatti = len(setesatti)  
            
            freqmatches = round( totesatti*100.0 / Nseq[clade] ,  1 )
            
            totmismatch[clade] = maxmism
            
            
            all_dgenset[pname][clade] = {}
            all_dgenset[pname][clade]["SetGeniEsatti"] = setesatti
            all_dgenset[pname][clade]["MismatchPerSequenza"] = mismperseq
            all_dgenset[pname][clade]["MismatchUltimiNucleotidi"] = mismLASTperseq
            
            

#            if "Sas" in pname  and clade == "clade4":
#                print "----------------------------------------------------------------------------------------------"
#                print "------------------------------------SAS---------------------------------------------------"                
#                print classifica_esatti
#                print "----------------------------------------------------------------------------------------------"
#                print setesatti
#                print totesatti
#                print freqmatches
#                print "------------------------------------SAS---------------------------------------------------"       
#                print "----------------------------------------------------------------------------------------------"            

            if thewinner_exact != thewinner_mism :      # se le combinazioni sono diverse
                if (thewinner_exact[0] != thewinner_mism[0] )   or  (thewinner_exact[3] != thewinner_mism[3])  or  (thewinner_exact[4] != thewinner_mism[4]) :     #..e se il numero di esatti e' diverso O se il numero di mismatch e' diverso
                    print "Risultati discordanti!!!"
                    print thewinner_exact
                    print thewinner_mism
                    print "QUELLO SCELTO E' SEMPRE IL MIGLIORE FRA I MISMATCH, OVVERO:"
                    print thewinner_mism
            
            print clade ,  freqmatches ,    "massimo numero mismatch =" ,  numismatch[clade] ,  totmismatch[clade] 

            #book_row += 1
            book_col = 1
            sheet1.write(book_row, book_col, clade, style)
            book_col += 1
            sheet1.write(book_row, book_col, (str(freqmatches) + '%' ), style)
            book_col += 1
            sheet1.write(book_row, book_col, ( '('+ str(Nseq[clade] )+')'  ), style)
            book_col += 1
            
            for el in thewinner_mism[1] :
                if el != 'fintoprobe' :
                    sheet1.write(book_row, book_col, el, style)
                    book_col += 1
            
            book_col += 1
            book_col += 1
            sheet1.write(book_row, book_col,   str(numismatch[clade] ) , style)
            book_col += 1
            sheet1.write(book_row, book_col,   str(totmismatch[clade] ) , style)
            book_row +=1


        #for clade in cladeviruslist :
        #   print clade ,  len(dgenset[clade])*100.0 / Nseq[clade]




mystring = 'MismatchCounts_'

totmism = []

for pname in primerstype :
    for sonda in  primerstype[pname] :
        if sonda != 'type' :
            for clade in cladeviruslist :
                if mystring+clade in primerstype[pname][sonda].keys() :
                    #print sonda , clade , primerstype[pname][sonda][mystring+clade]
                    totmism += primerstype[pname][sonda][mystring+clade]
                else :
                    print sonda ,  'CONTEGGIO VARIANTI NON PERVENUTO !!!'
                



#FIND BEST COMBINATION OF ASSAYS :

ptotest=list(primersorder)   # lista dei primer da testare, da cui elimino quelli con problemi per la target region :
ptotest.remove('Deyde_et_al.2006')
ptotest.remove('Atkinson_et_al.2012')
ptotest.remove('Bonney_et_al.2017')


Na = len(ptotest)   # number of assays to test
Ncombs = 5  #max numerosity of combinations to test

#prima di tutto, vedo se c'e' un assay che li becca tutti :
cladeGiaPerfetti = []
# ...e faccio una classifica dei migliori assay :
primerparade = {}


# aggiungo il clade "TUTTICLADE"
for pname in primerstype : 
    all_dgenset[pname]["tutticlade"] = {} 
    allsetesatti = set()
    allmismatches = []
    allastmismatches = []
    for clade in cladeviruslist :
        allsetesatti= all_dgenset[pname][clade]['SetGeniEsatti'] | allsetesatti
        allmismatches = allmismatches + all_dgenset[pname][clade]["MismatchPerSequenza"]
        allastmismatches = allastmismatches + all_dgenset[pname][clade]["MismatchUltimiNucleotidi"]
        
    all_dgenset[pname]["tutticlade"]['SetGeniEsatti'] = allsetesatti
    all_dgenset[pname]["tutticlade"]["MismatchPerSequenza"] = allmismatches
    all_dgenset[pname]["tutticlade"]["MismatchUltimiNucleotidi"] = allastmismatches


# e aggiungo TUTTICLADE alla lista dei clade :
cladeviruslist.append("tutticlade")
Nseq["tutticlade"] = len(newworldvirus)

for clade in  cladeviruslist :
    primerparade[clade] = []
    for pname in ptotest :
        misms = all_dgenset[pname][clade][ 'MismatchPerSequenza']       #lista del numeto di mismatch del primer con le sequenze in quel clade
        
        mismLAST = all_dgenset[pname][clade][ 'MismatchUltimiNucleotidi']       #lista del numeto di mismatch NEGLI ULTIMI NUCLEOTIDI del primer con le sequenze in quel clade
        
        NONesatti =  Nseq[clade] -  len( all_dgenset[pname][clade]['SetGeniEsatti'])        # numero di sequenze NON beccate = numero di sequenze totali - numero di sequenze esatte 
        
        if max(misms) == 0 :
            cladeGiaPerfetti.append( clade )
            cladeGiaPerfetti.append( pname )
        
        primerparade[clade].append( (NONesatti, max(mismLAST) , max(misms),  pname) )
    
    primerparade[clade].sort()        #  <- ordino prima per numero di sequenze (non) beccate , poi per il numero di mismatch massimi
        

#poi lavoro sui clade non perfetti :

bestnameLAST  = {}
bestnameNUM = {}
bestnamePERF  = {}
bestnameSUM = {}

bestLASTmism = {}
bestNUMmism = {}
bestPERFmatch = {}
bestSUMmism = {}

genomibeccati  = {}
genomiesatti = {}


for clade in  [ cladeviruslist[9] ] :
    print clade
    if clade not in cladeGiaPerfetti:
        bestnameLAST[clade] = [ (primerparade[clade][0][3], )  ] # nome del miglior assay
        bestnameNUM[clade] = [ (primerparade[clade][0][3], )  ] # nome del miglior assay
        bestnamePERF[clade] = [ (primerparade[clade][0][3], )  ] # nome del miglior assay
        bestnameSUM[clade] = [ (primerparade[clade][0][3], )  ] # nome del miglior assay
        
        bestLASTmism[clade] = [ all_dgenset[ bestnameLAST[clade][0][0] ][clade]['MismatchUltimiNucleotidi'] ]
        maxlast = max( all_dgenset[ bestnameLAST[clade][0][0] ][clade]['MismatchUltimiNucleotidi']  )
        
        bestNUMmism[clade] = [ all_dgenset[ bestnameLAST[clade][0][0] ][clade]['MismatchPerSequenza'] ]
        maxnum = max( all_dgenset[ bestnameLAST[clade][0][0] ][clade]['MismatchPerSequenza']  )
        
        bestSUMmism[clade] =  [ all_dgenset[ bestnameLAST[clade][0][0] ][clade]['MismatchPerSequenza'] ]
        sumnum = sum( all_dgenset[ bestnameLAST[clade][0][0] ][clade]['MismatchPerSequenza']  )
        
        bestPERFmatch[clade] =  all_dgenset[bestnameLAST[clade][0][0] ][clade]['SetGeniEsatti']  
        
        for k in range(1, Ncombs+1) :
            print k
            for comb in iter.combinations(ptotest, k) :
                    lastmism  = [6] *Nseq[clade]
                    numsism =  [100]*Nseq[clade]
                    exgen = set()
                    
                    
                    for assay in comb :
                        exgen = exgen | all_dgenset[assay][ clade ]['SetGeniEsatti']
                        
                        for n in range(Nseq[clade]) :
                            lastmism[n] = min(lastmism[n] , all_dgenset[assay][clade]['MismatchUltimiNucleotidi'][n]  )
                                            
                            numsism[n] = min( numsism[n] ,  all_dgenset[assay][clade]['MismatchPerSequenza'][n]   )
                    
#                    if clade == "clade1" :
#                        print comb
#                        print lastmism
#                        print numsism
#                        print exgen
                    
                    # TEST 1 : comparo il numero di mismatches negli ultimi 5 nucleotidi ##########################################
                    if  max(lastmism) < maxlast   :
                        bestnameLAST[clade] = [ comb ]
                        bestLASTmism[clade] = [ lastmism ]
                        maxlast = max(lastmism) 
                        
                               #stessto numero mismatch                                       # scelgo la combinazione piu' corta
#                    elif  max(lastmism) == maxlast : # and len(comb) < len(bestnameLAST[clade][len(bestnameLAST[clade]) -1])  :
#                        bestnameLAST[clade] = [ comb ]
#                        bestLASTmism[clade] = [ lastmism ]
#                        maxlast = max(lastmism) 
                                            
                             #stessto numero mismatch       #no tutti zeri                      # scelgo la combinazione piu' corta
                    elif  max(lastmism) == maxlast and comb not in bestnameLAST[clade]  : # and maxlast != 0 : # and len(comb) == len(bestnameLAST[clade][len(bestnameLAST[clade]) -1])  :
                        bestnameLAST[clade].append( comb )
                        bestLASTmism[clade].append(lastmism) 
                    
#                    elif  max(lastmism) == maxlast  and maxlast == 0  : # and len(comb) == len(bestnameLAST[clade][len(bestnameLAST[clade]) -1])  :
#                        bestnameLAST[clade].append( comb )
#                        # bestLASTmism[clade].append(lastmism) <- se sono tutti zero, mi basta una lista sola !
                    
                    
                    # TEST 2 : comparo il numero di mismatches nell'intera regione target ##########################################
                    if  max(numsism) < maxnum    :
                        bestnameNUM[clade] = [ comb ]
                        bestNUMmism[clade] = [ numsism ]
                        maxnum = max(numsism ) 

#                    elif  max(numsism) == maxnum : #  and  len(comb) < len(bestnameNUM[clade][len(bestnameNUM[clade]) -1])   :
#                        bestnameNUM[clade] = [ comb ]
#                        bestNUMmism[clade] = [ numsism ]
#                        maxnum = max(numsism) 
                        
                    elif  max(numsism) == maxnum and comb not in bestnameNUM[clade]   :  #and maxnum != 0 : # and  len(comb) == len(bestnameNUM[clade][len(bestnameNUM[clade]) -1])     :
                        bestnameNUM[clade].append( comb )
                        bestNUMmism[clade].append(numsism) 
                    
#                    elif  max(numsism) == maxnum  and maxnum == 0  : # and  len(comb) == len(bestnameNUM[clade][len(bestnameNUM[clade]) -1])  :
#                        bestnameNUM[clade].append( comb )
#                        # bestNUMmism[clade].append(numsism) <- se sono tutti zero, mi basta una lista sola !
                        
                    
                    # TEST 3 : comparo il numero di sequenze matchate al 100% ##########################################                    
                    if  len(exgen) > len(bestPERFmatch[clade] )    :
                        bestnamePERF[clade] = [ comb ]
                        bestPERFmatch[clade] = exgen
                        
#                    elif  len(exgen) ==  len(bestPERFmatch[clade] ) : #  and  len(comb) < len(bestnamePERF[clade][len(bestnamePERF[clade]) -1]) :
#                        bestnamePERF[clade] = [ comb ]
#                        bestPERFmatch[clade] = exgen
                        
                    elif  len(exgen) == len(bestPERFmatch[clade] )  and comb not in bestnamePERF[clade] : #  and  len(comb) == len(bestnamePERF[clade][len(bestnamePERF[clade]) -1]) :
                        bestnamePERF[clade].append( comb )
                        
                        
                    # TEST 4 : comparo il numero di mismatches nell'intera regione target ##########################################
                    if  sum(numsism) < sumnum    :
                        bestnameSUM[clade] = [ comb ]
                        bestSUMmism[clade] = [ numsism ]
                        sumnum = sum(numsism )
                                               
                    elif  sum(numsism) == sumnum  and comb not in bestnameSUM[clade]  :  #and maxnum != 0 : # and  len(comb) == len(bestnameNUM[clade][len(bestnameNUM[clade]) -1])     :
                        bestnameSUM[clade].append( comb )
                        bestSUMmism[clade].append(numsism) 
                        

bests = {}
for clade in [ cladeviruslist[9] ] :
    if clade not in cladeGiaPerfetti:
        bests[clade] =   list(  (( set(bestnameSUM[clade]) & set(bestnameNUM[clade]) ) &   set(bestnameLAST[clade])    )   & set(bestnamePERF[clade]) )
        
        bests[clade].sort(key=len)
        
        print '---------------------------------------------------------------------------------------------------------------'
        print clade
        print("Il miglioor assay per il "+clade+ " e' :")
        bestassay=primerparade[clade][0][3]
        print bestassay
        print "MismatchPerSequenza:"
        print all_dgenset[bestassay][clade][ 'MismatchPerSequenza']
        print max(all_dgenset[bestassay][clade][ 'MismatchPerSequenza'])
        print "MismatchUltimiNucleotidi:"
        print all_dgenset[bestassay][clade][ 'MismatchUltimiNucleotidi']   
        print max(all_dgenset[bestassay][clade][ 'MismatchUltimiNucleotidi']   )
        print "SetGeniEsatti :"
        print  all_dgenset[bestassay][clade]['SetGeniEsatti']
        print round( len(all_dgenset[bestassay][clade]['SetGeniEsatti'])*100.0 / Nseq[clade] ,  1 )
        
        for i in range(5) :
            if len(list(bests[clade])) > i :
                print("Una delle migliori combinazioni e' ")
                mybestcomb = list(bests[clade])[i]
                print mybestcomb
                
                qualeSUM = bestnameSUM[clade].index(mybestcomb)
                mynummism = bestSUMmism[clade][qualeSUM]
                
                qualeLAST = bestnameLAST[clade].index(mybestcomb)            
                mylastmism   =      bestLASTmism[clade][qualeLAST]
                
                qualePERF = bestnamePERF[clade].index(mybestcomb)
                myexgen=   bestPERFmatch[clade]
                
                print "MismatchPerSequenza:"
                print mynummism
                print max(mynummism)
                print "MismatchUltimiNucleotidi:"
                print mylastmism
                print  max(mylastmism)
                print "SetGeniEsatti :"
                print  myexgen
                print round( len(myexgen)*100.0 / Nseq[clade] ,  1 )


print "QUINDI:"

for clade in [ cladeviruslist[9]  ]:
    if clade not in cladeGiaPerfetti:
        print clade
        print bests[clade][:3]



book.save(myworkbook)             #file exel dei primer FORMATTATO !!!

outfasta.close()
outf.close()



#
#for clade in  cladeviruslist :
#    if clade not in cladeGiaPerfetti:
#        
#        if  primerparade[clade][1][0]  <  Nseq[clade]   :    # se ci sono almeno 2 assays che beccano qualche sequenza al 100% , vedo se combinando gli assay migliora la situazione
#            
#            #seqesatte = Nseq[clade]  - primerparade[clade][0][0]      # il numero di sequenze esatte
#            
#            bestname = primerparade[clade][0][2]  # nome del miglior assay
#            
#            bestcombination[clade] = [  bestname  ]   # come minimo, il miglior modo per beccare una sequenza e' usare il miglior assay !
#            
#            bestcombination["NumeroGenomiMancanti_"+clade] =  primerparade[clade][0][0] 
#            
#            genomiesatti[clade] = all_dgenset[ bestname ][ clade ]['SetGeniEsatti']    # quali sequenze becca il migliore ?
#            
#            genomibeccati[clade] = list(all_dgenset[ bestname ][ clade ]['SetGeniEsatti'] )       # creo una lista, per poter contare gli stessi genomi piu' volte
#            
#            # Ora aggiungo gli assay, uno ad un (se ottengo un numero piu' alto di sequenze esatte )
#            for tupla in  primerparade[clade] :        
#                
#                addassayname = tupla[2]
#                
#                addgenomiesatti =  all_dgenset[addassayname][ clade ]['SetGeniEsatti']  |  genomiesatti[clade] 
#                
#                if len(addgenomiesatti ) > len(genomiesatti[clade]) :
#                    
#                    bestcombination[clade].append(addassayname)
#                    genomiesatti[clade] = addgenomiesatti
#                    
#                    genomibeccati[clade] = genomibeccati[clade] + list( all_dgenset[addassayname][ clade ]['SetGeniEsatti']   )            # qui conto tutte le volte che ho trovato un genoma
#                    
#            
#            bestcombination["NumeroGenomiMancanti_"+clade] =  Nseq[clade] -  len(genomiesatti[clade] )
#
#print "BEST COMBINATIONS:"
#for clade in cladeviruslist:
#    if clade in bestcombination.keys() :
#        print clade 
#        print bestcombination[clade]
#        print bestcombination['NumeroGenomiMancanti_'+clade]
#
#
#
#
## Ora controllo di non aver agginto assay inutilmente :
#            
#best_smart_combination = {}
#for k in bestcombination :
#    best_smart_combination[k] = bestcombination[k]
#                
#for clade in  cladeviruslist:
#    if clade in bestcombination.keys() :               
#            
#        for myassay in  bestcombination[clade] :        
#                
#            removeassayname = myassay   
#            
#            removegenomes =  all_dgenset[removeassayname][ clade ]['SetGeniEsatti']   
#            
#            removegenomibeccati = list(genomibeccati[clade])
#            
#            for genome in removegenomes  :                    
#                removegenomibeccati.remove(genome)                 # rimuovo ad uno ad uno i genomi beccati dell'ultimo assay aggiunto
#            
#            if set(removegenomibeccati)  == genomiesatti[clade]     :   # se questo non cambia il numero totale dei genomi esatti
#                best_smart_combination[clade].remove(myassay)
#            
#            
#print "BEST SMART COMBINATIONS:"
#for clade in cladeviruslist:
#    if clade in best_smart_combination.keys() :
#        print clade 
#        print best_smart_combination[clade]
#        print best_smart_combination['NumeroGenomiMancanti_'+clade]
#        
#        
#        










#PER fare l'istogramma in R:
#totmism1 <- c(3,3,2, 2, 2, 2, 2, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 2, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 1, 3, 3, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 2, 2, 2, 3, 2, 2, 1, 4, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 3, 2, 3, 2, 3, 3, 3, 3, 3, 3, 2, 2, 3, 3, 4, 3, 2, 2, 3, 3, 2, 3, 2, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 2, 2, 2, 2, 2, 2, 3, 2, 3, 2, 2, 2, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 2, 4, 4, 3, 2, 3, 3, 2, 0, 0, 0, 0, 0, 1, 0, 2, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 3, 0, 0, 0, 0, 0, 2, 3, 2, 3, 2, 2, 2, 2, 3, 3, 3, 2, 2, 3, 1, 1, 3, 1, 3, 2, 3, 3, 1, 3, 3, 1, 4, 3, 1, 2, 3, 2, 1, 3, 2, 3, 3, 2, 1, 3, 2, 1, 2, 3, 4, 3, 3, 3, 3, 3, 4, 3, 3, 3, 3, 3, 3, 4, 4, 3, 3, 3, 3, 3, 3, 3, 4, 3, 3, 4, 3, 3, 3, 3, 3, 3, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 7, 4, 4, 7, 1, 1, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 3, 0, 0, 0, 0, 0, 1, 1, 1, 2, 1, 1, 6, 4, 2, 2, 2, 3, 1, 2, 3, 3, 3, 3, 3, 1, 1, 1, 3, 1, 2, 3, 2, 2, 2, 3, 2, 1, 3, 3, 4, 2, 1, 1, 3, 2, 1, 3, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 3, 2, 2, 2, 2, 1, 1, 2, 2, 2, 3, 2, 2, 2, 2, 2, 3, 2, 2, 3, 2, 2, 2, 2, 2, 5, 3, 3, 4, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 4, 2, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 2, 1, 0, 1, 1, 1, 2, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 3, 3, 3, 2, 1, 2, 0, 1, 0, 0, 2, 0, 0, 13, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 1,0)
#              
#totmism2 <- c(0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 2, 2, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 2, 0, 1, 2, 2, 0, 1, 1, 0, 1, 0, 2, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 2, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 1, 2, 1, 1, 2, 0, 0, 0, 0, 1, 1, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 1, 1, 2, 1, 1, 2, 1, 2, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 2, 2, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 2, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 3, 3, 4, 1, 1, 2, 1, 2, 2, 1, 3, 1, 3, 2, 2, 2, 2, 1, 3, 2, 1, 2, 2, 3, 2, 1, 3, 2, 3, 2, 3, 1, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 2, 2, 2, 1, 2, 2, 2, 1, 2, 1, 1, 1, 1, 2, 1, 2, 2, 2, 2, 2, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 7, 7, 6, 6, 7, 3, 3, 3, 3, 3, 3, 4, 3, 3, 3, 3, 3, 5, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 4, 3, 3, 5, 3, 3, 3, 4, 3, 5, 4, 4, 4, 4, 4, 4, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 5, 4, 4, 4, 4, 4, 4, 4, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 4, 4, 4, 4, 4, 4, 4, 4, 2, 0, 0, 2, 4, 4, 3, 4, 3, 3, 4, 4, 2, 4, 5, 6, 5, 4, 4, 4, 5, 2, 4, 5, 4, 5, 4, 4, 6, 4, 6, 4, 4, 5, 6, 5, 5, 4, 3, 3, 3, 3, 3, 4, 5, 5, 3, 2, 2, 5, 3, 2, 5, 4, 2, 5, 2, 3, 3, 3, 5, 3, 2, 5, 5, 2, 5, 5, 3, 3, 5, 2, 5, 2, 4, 3, 5, 3, 3, 5, 3, 3, 5, 5, 5, 5, 5, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 4, 4, 4, 5, 5, 5, 5, 5, 4, 5, 5, 5, 5, 0, 0, 5, 2, 2, 5, 5, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 5, 5, 4)
#              
#totmism3 <- c(4, 4, 6, 4, 4, 6, 6, 5, 6, 5, 3, 4, 3, 5, 3, 4, 6, 5, 4, 4, 5, 4, 4, 6, 5, 5, 4, 4, 4, 6, 4, 4, 6, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 4, 6, 4, 4, 5, 4, 3, 3, 4, 4, 4, 6, 4, 4, 4, 4, 4, 5, 4, 5, 4, 4, 4, 4, 4, 4, 4, 0, 0, 3, 4, 3, 5, 6, 5, 5, 6, 6, 4, 4, 5, 5, 5, 5, 6, 5, 4, 4, 5, 5, 5, 5, 5, 4, 4, 5, 4, 5, 6, 5, 4, 5, 4, 5, 5, 6, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 5, 4, 6, 6, 4, 6, 5, 5, 5, 6, 5, 6, 6, 6, 6, 6, 6, 6, 6, 5, 6, 6, 6, 6, 6, 5, 6, 5, 5, 6, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 5, 6, 6, 6, 6, 6, 6, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 5, 6, 6, 6, 6, 6, 6, 6, 6, 5, 0, 0, 5, 2, 2, 2, 2, 2, 2, 2, 0, 0, 20, 0, 0, 0, 0, 5, 0, 0, 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 4, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 3, 1, 1, 4, 2, 3, 2, 3, 2, 2, 3, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 2, 1, 0, 0, 1, 2, 1, 1, 3, 0, 0, 13, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 18, 0, 0, 21, 21, 1, 21, 9, 0, 0, 0, 15, 0, 0, 0, 13, 8, 0, 21, 21, 0, 0, 0, 0, 0, 0, 13, 0, 13, 0, 13, 9, 2, 9, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 21, 0, 18, 18, 0, 0, 0, 0, 0, 18, 0, 4, 0)
#              
#totmism4 <- c( 0, 8, 18, 17, 0, 1, 1, 0, 18, 0, 9, 18, 0, 0, 18, 0, 1, 0, 1, 18, 0, 15, 21, 15, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 7, 0, 1, 15, 0, 0, 0, 0, 1, 0, 0, 0, 1, 18, 21, 0, 1, 18, 0, 1, 16, 1, 18, 1, 0, 1, 21, 0, 0, 21, 0, 0, 0, 21, 0, 0, 21, 21, 0, 21, 0, 0, 0, 0, 21, 0, 0, 0, 0, 0, 0, 21, 21, 21, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 21, 0, 21, 21, 0, 0, 0, 0, 0, 21, 0, 21, 0, 0, 21, 21, 21, 0, 0, 0, 0, 21, 0, 21, 21, 0, 0, 21, 0, 2, 0, 2, 14, 0, 18, 2, 21, 2, 0, 2, 0, 1, 0, 2, 1, 0, 1, 3, 9, 0, 2, 21, 0, 0, 0, 0, 2, 21, 0, 0, 2, 18, 2, 2, 2, 21, 0, 2, 21, 2, 14, 2, 0, 2, 21, 0, 0, 21, 0, 0, 0, 19, 0, 0, 19, 19, 0, 19, 0, 0, 0, 0, 19, 0, 0, 0, 0, 0, 0, 19, 19, 19, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 19, 0, 19, 19, 0, 0, 0, 0, 0, 19, 0, 19, 0, 0, 19, 19, 19, 0, 0, 0, 0, 19, 0, 19, 19, 0, 0, 19, 0, 2, 0, 2, 14, 0, 18, 2, 19, 2, 0, 2, 0, 1, 0, 2, 1, 0, 1, 3, 9, 0, 2, 19, 0, 0, 0, 0, 2, 19, 0, 0, 2, 18, 2, 2, 2, 19, 0, 2, 19, 2, 14, 2, 0, 2, 19, 0, 0, 19, 0, 0, 1, 26, 1, 1, 9, 3, 1, 26, 1, 1, 1, 1, 9, 2, 2, 1, 1, 1, 1, 26, 10, 2, 1, 1, 1, 1, 1, 1, 2, 2, 2, 1, 2, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 26, 1, 4, 2, 1, 1, 1, 2, 1, 18, 1, 2, 1, 1, 4, 3, 2, 1, 1, 1, 1, 4, 1, 14, 2, 1, 1, 4, 1, 1, 4, 2, 2, 3, 2, 1, 20, 1, 1, 2, 1, 2, 1, 2, 1, 3, 1, 2, 1, 1, 1, 18, 1, 1, 1, 1, 2, 26, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 25, 1, 2, 2, 3, 2, 3, 2, 2, 5, 3, 3, 2, 3, 2, 2, 3, 1, 1, 21, 0, 0, 0, 1, 2, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 2, 0, 1, 0, 0, 1, 2, 2, 2, 2, 2, 1, 3, 4, 2, 2, 2, 4, 1, 2, 3, 4, 2, 3, 3, 1, 2, 2, 3, 2, 2, 3, 2, 2, 3, 3, 2, 2, 3, 2, 4, 2, 2, 2, 3, 2, 2, 3, 1, 2, 2, 3, 2, 1, 3, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 1, 3, 2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 3, 2, 2, 2, 1, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 4, 1, 4, 2, 4, 3, 1, 2, 2, 2, 2, 1, 3, 1, 1, 1, 1, 1, 1, 3, 3, 1, 1, 2, 2, 2, 1, 3, 1, 1, 1, 4, 1, 1, 2, 2, 1, 3, 1, 2, 4, 1, 3, 1, 1, 0, 1, 1, 1, 1)
#
#totmism5 <- c( 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 2, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 3, 1, 1, 3, 3, 3, 3, 3, 2, 2, 3, 3, 3, 3, 2, 3, 2, 2, 2, 2, 3, 3, 3, 2, 3, 3, 2, 3, 2, 3, 3, 3, 2, 2, 2, 2, 2, 2, 3, 2, 2, 3, 2, 3, 3, 3, 4, 2, 4, 4, 3, 3, 4, 2, 2, 2, 2, 4, 2, 3, 4, 2, 4, 3, 4, 2, 2, 3, 3, 2, 4, 2, 2, 4, 3, 2, 4, 2, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 2, 3, 3, 2, 3, 3, 2, 2, 2, 2, 2, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 2, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 2, 0, 1, 1, 2, 0, 1, 2, 2, 1, 2, 1, 2, 0, 0, 2, 0, 1, 2, 0, 1, 2, 2, 0, 1, 2, 1, 2, 1, 0, 1, 2, 0, 1, 2, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 1, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 1, 3, 3, 3, 4, 3, 3, 4, 2, 4, 2, 1, 2, 2, 2, 3, 1, 2, 4, 2, 3, 2, 1, 2, 2, 1, 1, 1, 2, 3, 2, 1, 2, 2, 2, 4, 3, 4, 3, 4, 2, 3, 2, 3, 2, 1, 0, 3, 3, 0, 0, 4, 0, 3, 4, 3, 3, 0, 3, 1, 0, 3, 2, 2, 0, 3, 4, 0, 3, 2, 2, 5, 4, 0, 3, 4, 0, 4, 3, 3, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 3, 3, 4, 4, 4, 3, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 0, 0, 4, 2, 2, 3, 2, 3, 3, 2, 2, 2, 2, 1, 1, 1, 2, 2, 3, 1, 2, 1, 1, 2, 1, 0, 2, 2, 3, 2, 3, 2, 1, 2, 1, 1, 1, 2, 2, 2, 2, 2, 2, 1, 0, 2, 2, 2, 1, 2, 3, 1, 1, 2, 1, 2, 2, 2, 2, 1, 2, 2, 1, 3, 2, 1, 1, 2, 2, 1, 2, 0, 2, 2, 2, 1, 2, 2, 1, 2, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 7, 7, 6, 1, 1, 2, 1, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0, 3, 0, 2, 0, 0, 0, 2, 1, 2, 2, 2, 2, 1, 1, 2, 3, 3, 1, 2, 2, 1, 1, 2, 1, 1, 2, 2, 1, 1, 1, 3, 1, 2, 3, 1, 1, 1, 2, 1, 2, 1, 3, 2, 2, 1, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2, 2, 3, 2, 2, 2, 2, 2)
#
#totmism6 <- c( 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 3, 1, 1, 3, 2, 2, 3, 3, 3, 3, 3, 0, 0, 18, 1, 1, 1, 0, 5, 0, 1, 0, 0, 2, 0, 1, 0, 0, 2, 0, 0, 0, 4, 1, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 3, 2, 1, 2, 2, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 3, 2, 2, 3, 2, 2, 2, 2, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2, 3, 1, 1, 2, 3, 2, 4, 1, 1, 5, 2, 2, 3, 6, 3, 3, 6, 1, 1, 23, 1, 1, 1, 1, 1, 2, 3, 1, 2, 1, 1, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 2, 3, 1, 0, 1, 0, 1, 0, 0, 0, 2, 1, 2, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 2, 1, 0, 2, 1, 1, 1, 0, 1, 1, 2, 2, 2, 0, 1, 1, 0, 2, 0, 1, 0, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 15, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 4, 4, 5, 3, 3, 3, 19, 3, 3, 3, 1, 1, 19, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 19, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 2, 0, 2, 0, 1, 2, 1, 2, 0, 0, 2, 0, 0, 2, 2, 0, 19, 0, 0, 0, 1, 2, 1, 0, 2, 6, 0, 2, 2, 1, 1, 3, 1, 1, 0, 0, 0, 2, 2, 2, 2, 0, 1, 0, 3, 3, 5, 4, 4, 5, 4, 3, 12, 3, 3, 4, 3, 4, 3, 4, 3, 5, 3, 4, 3, 3, 3, 9, 3, 3, 3, 3, 4, 19, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 17, 3, 4, 4, 3, 4, 5, 3, 3, 5, 2, 4, 0, 0, 0, 0, 0, 4, 4, 5, 5, 5, 5, 4, 1, 3, 4, 4, 5, 5, 5, 6, 4, 5, 3, 3, 4, 3, 1, 5, 3, 4, 4, 4, 2, 2, 2, 2, 2, 3, 1, 2, 2, 2, 2, 1, 3, 2, 1, 1, 2, 1, 2, 2, 2, 2, 1, 2, 2, 1, 3, 2, 1, 1, 2, 2, 1, 2, 2, 2, 2, 2, 1, 2, 2, 1, 2, 2, 3, 2, 1, 2, 2, 2, 2, 2, 2, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 1, 1, 2, 2, 3, 2, 3, 4, 4, 3, 3, 5, 2, 2, 2, 2, 2, 4, 4, 5, 5, 5, 5, 4, 2, 3, 4, 4, 5, 5, 5, 4, 4, 5, 3, 3, 4, 3, 2, 5, 3, 4, 4, 4, 0, 1, 0, 1, 0, 1, 2, 1, 1, 1, 1, 2, 1, 1, 2, 2, 1, 2, 1, 0, 0, 1, 2, 1, 1, 2, 2, 1, 2, 2, 1, 0, 2, 1, 1, 1, 1, 0, 2, 1, 0, 2, 0, 1, 2, 1, 2, 1, 1, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2, 2, 1, 1, 2)
#
#totmism7 <- c( 1, 2, 5, 5, 2, 3, 3, 4, 4, 4, 4, 4, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0, 2, 1, 2, 1, 2, 1, 0, 1, 1, 2, 2, 1, 1, 2, 1, 1, 2, 1, 2, 3, 1, 1, 1, 1, 2, 1, 1, 2, 1, 1, 1, 2, 1, 2, 1, 2, 1, 2, 1, 1, 2, 1, 2, 2, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 1, 1, 2, 1, 1, 1, 1, 1, 5, 5, 6, 6, 6, 6, 5, 0, 4, 5, 5, 6, 6, 6, 6, 5, 6, 4, 4, 5, 4, 0, 6, 4, 5, 5, 5, 2, 1, 2, 1, 2, 3, 0, 1, 1, 1, 1, 0, 3, 1, 0, 0, 1, 0, 1, 2, 2, 1, 0, 1, 1, 0, 2, 1, 0, 0, 1, 2, 0, 1, 1, 1, 1, 2, 0, 1, 2, 0, 2, 1, 2, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 2, 1, 3, 2, 2, 2, 0, 1, 2, 2, 2, 2, 2, 6, 6, 7, 7, 7, 7, 6, 1, 5, 6, 6, 7, 7, 7, 7, 6, 7, 5, 5, 6, 5, 1, 7, 5, 6, 6, 6, 3, 2, 3, 2, 3, 4, 1, 2, 2, 2, 2, 1, 4, 2, 1, 1, 2, 1, 2, 3, 3, 2, 1, 2, 2, 1, 1, 2, 1, 1, 2, 3, 1, 2, 2, 2, 2, 3, 1, 2, 3, 1, 3, 2, 1, 2, 1, 2, 2, 2, 2, 2, 2, 1, 1, 1, 2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 1, 1, 2, 2, 3, 2, 4, 2, 2, 3, 5, 5, 4, 5, 4, 4, 5, 3, 4, 2, 0, 0, 0, 2, 2, 2, 2, 4, 2, 0, 2, 0, 2, 2, 2, 2, 1, 2, 2, 0, 2, 2, 2, 2, 1, 1, 1, 1, 1, 2, 2, 1, 1, 2, 2, 2, 2, 3, 1, 1, 4, 1, 2, 2, 2, 1, 1, 1, 2, 1, 1, 3, 3, 1, 1, 1, 1, 3, 1, 2, 2, 1, 1, 1, 1, 1, 2, 1, 1, 5, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 3, 5, 4, 4, 4, 4, 3, 4, 3, 3, 3, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 4, 4, 4, 4, 3, 4, 5, 4, 4, 5, 5, 3, 4, 4, 3, 3, 3, 3, 3, 2, 2, 1, 1, 1, 1, 1, 3, 1, 1, 2, 1, 1, 1, 1, 2, 1, 1, 1, 0, 1, 3, 1, 1, 1, 1, 1, 3, 2, 3, 2, 3, 2, 1, 2, 2, 3, 3, 2, 2, 3, 2, 2, 3, 2, 3, 4, 2, 2, 2, 2, 3, 2, 2, 3, 2, 2, 2, 3, 2, 3, 2, 3, 2, 3, 2, 2, 3, 2, 3, 3, 2, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 4, 2, 5, 3, 2, 3, 2, 2, 3, 3, 2, 2, 2, 2, 2, 2, 2)
#
#totmism8 <- c(2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 2, 1, 2, 2, 2, 0, 2, 2, 2, 3, 1, 3, 1, 3, 2, 0, 1, 1, 2, 2, 2, 2, 3, 1, 1, 4, 1, 2, 2, 2, 1, 1, 1, 2, 1, 1, 3, 1, 1, 1, 3, 1, 3, 1, 2, 2, 3, 1, 1, 3, 1, 2, 1, 1, 3, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 1, 3, 2, 2, 2, 2, 3, 2, 1, 3, 3, 2, 2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 2, 1, 2, 3, 2, 4, 5, 5, 3, 6, 6, 4, 4, 4, 4, 4, 0, 2, 1, 1, 1, 1, 0, 5, 1, 0, 2, 1, 1, 1, 2, 0, 1, 3, 1, 2, 1, 5, 1, 3, 0, 0, 0, 4, 4, 4, 4, 4, 3, 5, 4, 4, 4, 4, 5, 3, 4, 5, 5, 4, 5, 4, 4, 4, 4, 5, 4, 4, 5, 5, 4, 5, 5, 4, 4, 5, 4, 4, 4, 4, 4, 5, 4, 4, 5, 4, 4, 5, 4, 5, 4, 4, 4, 4, 4, 4, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 4, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 4, 5, 5, 4, 4, 3, 4, 3, 6, 6, 3, 4, 4, 3, 4, 3, 3, 4, 2, 3, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 2, 0, 2, 0, 2, 1, 1, 0, 0, 1, 1, 1, 1, 2, 0, 0, 3, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 2, 2, 0, 0, 2, 0, 2, 0, 1, 1, 2, 0, 0, 2, 0, 1, 0, 0, 4, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 4, 3, 3, 3, 3, 2, 3, 2, 2, 4, 3, 3, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 3, 3, 3, 3, 2, 3, 4, 3, 3, 6, 6, 2, 1, 3, 1, 1, 1, 1, 1, 5, 5, 6, 6, 6, 6, 5, 0, 4, 5, 5, 6, 6, 6, 6, 5, 6, 4, 4, 5, 4, 0, 6, 4, 5, 5, 5, 2, 1, 2, 1, 2, 3, 0, 1, 1, 1, 1, 0, 3, 1, 0, 0, 1, 0, 1, 2, 2, 1, 0, 1, 1, 0, 2, 1, 0, 0, 1, 2, 0, 1, 1, 1, 1, 2, 0, 1, 2, 0, 2, 1, 2, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 1, 2, 1, 3, 3, 3, 2, 2, 3, 2, 2, 2, 2, 2, 4, 4, 5, 5, 5, 5, 4, 1, 3, 4, 4, 5, 5, 5, 5, 4, 5, 3, 3, 4, 3, 1, 5, 3, 4, 4, 4, 1, 0, 1, 0, 1, 2, 1, 0, 0, 0, 0, 1, 2, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 2, 3, 3, 1, 7, 7, 5, 5, 5, 5, 5, 1, 1, 0, 0, 0, 0, 1, 6, 2, 1, 1, 0, 0, 0, 1, 1, 0, 2, 2, 1, 2, 6, 0, 2, 1, 1, 1, 5, 5)
#              
#totmism9 <- c( 5, 5, 5, 4, 6, 5, 5, 5, 5, 6, 4, 5, 6, 6, 5, 6, 5, 5, 5, 5, 6, 5, 5, 6, 6, 5, 6, 6, 5, 5, 6, 5, 5, 5, 5, 5, 6, 5, 5, 6, 5, 5, 6, 5, 6, 5, 5, 5, 5, 5, 5, 6, 6, 6, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 5, 6, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 5, 6, 6, 5, 5, 4, 5, 4, 7, 7, 4, 2, 1, 3, 3, 3, 3, 3, 5, 5, 6, 6, 6, 6, 5, 2, 4, 5, 5, 6, 6, 6, 6, 5, 6, 4, 4, 5, 4, 2, 6, 4, 5, 5, 5, 4, 3, 4, 3, 4, 5, 2, 3, 3, 3, 3, 2, 5, 3, 2, 2, 3, 2, 3, 4, 4, 3, 2, 3, 3, 2, 3, 3, 2, 2, 3, 4, 2, 3, 3, 3, 3, 4, 2, 3, 4, 2, 4, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 2, 2, 3, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 2, 2, 3, 3, 4, 3, 3, 0, 0, 2, 5, 4, 3, 3, 3, 3, 3, 1, 1, 2, 2, 2, 2, 1, 4, 0, 1, 1, 2, 2, 2, 3, 1, 2, 2, 0, 1, 0, 4, 2, 2, 1, 1, 1, 3, 3, 3, 3, 3, 4, 4, 3, 3, 3, 3, 4, 4, 3, 4, 4, 3, 4, 3, 3, 3, 3, 4, 3, 3, 4, 4, 3, 4, 4, 3, 3, 4, 3, 3, 3, 3, 3, 4, 3, 3, 4, 3, 3, 4, 3, 4, 3, 3, 3, 3, 3, 3, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 3, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 3, 4, 4, 3, 3, 4, 3, 2, 4, 4, 2, 1, 0, 3, 3, 3, 3, 3, 5, 5, 6, 6, 6, 6, 5, 2, 4, 5, 5, 6, 6, 6, 6, 5, 6, 4, 4, 5, 4, 2, 6, 4, 5, 5, 5, 4, 3, 4, 3, 4, 5, 2, 3, 3, 3, 3, 2, 5, 3, 2, 2, 3, 2, 3, 4, 4, 3, 2, 3, 3, 2, 2, 3, 2, 2, 3, 4, 2, 3, 3, 3, 3, 4, 2, 3, 4, 2, 4, 3, 2, 3, 2, 3, 3, 3, 3, 3, 3, 2, 2, 2, 3, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 2, 2, 3, 3, 4, 3, 3, 1, 1, 2, 3, 4, 3, 3, 3, 3, 3, 3, 5, 4, 4, 4, 4, 3, 2, 4, 3, 5, 4, 4, 4, 4, 3, 4, 4, 4, 5, 4, 2, 4, 4, 3, 3, 3, 2, 1, 2, 1, 2, 1, 2, 1, 1, 1, 1, 2, 1, 1, 2, 2, 1, 2, 1, 2, 2, 1, 2, 1, 1, 2, 2, 1, 2, 2, 1, 2, 2, 1, 1, 1, 1, 2, 2, 1, 2, 2, 2, 1, 2, 1, 2, 1, 1, 1, 1, 1, 1, 2, 2, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2, 2, 1, 1, 0, 1, 3, 4, 4, 2, 2, 2, 5, 5, 5, 5, 5, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 3, 2, 3, 2, 3, 2, 1, 0, 2, 3, 3, 0, 2, 3, 0, 0, 3, 0, 3, 4, 2, 2)
#
#totmism10 <- c(0, 2, 3, 0, 2, 3, 0, 0, 2, 3, 0, 3, 0, 3, 2, 3, 0, 2, 3, 0, 3, 3, 2, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 4, 0, 6, 6, 4, 4, 4, 4, 4, 2, 0, 1, 1, 1, 1, 2, 5, 1, 2, 0, 1, 1, 1, 2, 2, 1, 1, 1, 0, 1, 5, 1, 1, 2, 2, 2, 4, 4, 4, 4, 4, 5, 5, 4, 4, 4, 4, 5, 5, 4, 5, 5, 4, 5, 4, 4, 4, 4, 5, 4, 4, 5, 5, 4, 5, 5, 4, 4, 5, 4, 4, 4, 4, 4, 5, 4, 4, 5, 4, 4, 5, 4, 5, 4, 4, 4, 4, 4, 4, 5, 5, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 4, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 4, 5, 5, 4, 4, 5, 4, 3, 6, 6, 3, 1, 2, 1, 2, 1, 1, 2, 2, 1, 2, 2, 2, 2, 2, 0, 2, 3, 1, 2, 2, 2, 3, 2, 2, 2, 2, 3, 2, 0, 2, 2, 3, 3, 2, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 3, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 4, 3, 2, 3, 2, 3, 4, 2, 2, 3, 3, 3, 3, 3, 2, 2, 3, 3, 3, 2, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 2, 8, 8, 1, 1, 1, 2, 2, 2, 2, 3, 0, 0, 1, 1, 1, 1, 0, 4, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 4, 1, 1, 0, 0, 0, 3, 2, 3, 3, 3, 3, 2, 2, 3, 3, 3, 2, 3, 3, 2, 2, 3, 2, 2, 3, 3, 2, 2, 2, 3, 2, 3, 4, 2, 2, 2, 3, 1, 3, 2, 3, 3, 3, 2, 3, 3, 2, 3, 3, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 3, 4, 4, 3, 3, 3, 3, 2, 3, 3, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 4, 4, 4, 3, 4, 4, 3, 4, 3, 3, 4, 4, 3, 3, 3, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 2, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 2, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 2, 2, 3, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 4, 2, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 0)
#               
#totmism11 <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 2, 0, 1, 1, 0, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 2, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 3, 3, 1, 1, 1, 1, 1, 2, 3, 2, 2, 2, 2, 3, 1, 2, 2, 3, 2, 2, 2, 2, 2, 2, 3, 2, 3, 2, 0, 2, 3, 2, 2, 2, 2, 1, 2, 0, 2, 2, 3, 2, 0, 0, 0, 2, 2, 0, 2, 2, 1, 2, 1, 1, 2, 1, 1, 1, 0, 2, 1, 0, 2, 1, 0, 2, 2, 1, 2, 0, 2, 2, 2, 0, 2, 2, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 2, 0, 0, 1, 0, 1, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 2, 4, 4, 2, 4, 4, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2, 1, 2, 1, 2, 1, 1, 1, 1, 0, 1, 0, 1, 0, 3, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 2, 0, 0, 1, 0, 0, 1, 2, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 3, 3, 1, 3, 3, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2, 1, 2, 2, 3, 3, 1, 0, 0, 2, 2, 1, 2, 2, 0, 2, 1, 2, 2, 2, 2, 2, 0, 2, 1, 0, 2, 2, 2, 2, 2, 0, 3, 0, 2, 2, 2, 1, 2, 2, 2, 1, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4)
#
#totmism12 <- c( 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 2, 4, 4, 2, 1, 1, 2, 3, 2, 2, 3, 3, 5, 3, 3, 3, 3, 3, 3, 3, 3, 5, 2, 3, 3, 3, 3, 4, 4, 3, 4, 3, 3, 3, 4, 3, 3, 3, 2, 1, 2, 1, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 2, 1, 1, 2, 1, 2, 2, 1, 1, 2, 1, 0, 2, 1, 0, 1, 0, 0, 0, 0, 0, 1, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 3, 1, 1, 2, 3, 3, 3, 3, 3, 3, 3, 2, 1, 2, 3, 2, 2, 2, 1, 2, 2, 1, 2, 2, 2, 2, 2, 2, 1, 2, 1, 2, 1, 2, 1, 2, 2, 2, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 2, 4, 4, 1, 1, 1, 1, 1, 2, 1, 2, 1, 1, 2, 1, 2, 1, 1, 1, 2, 2, 2, 1, 1, 2, 1, 1, 2, 1, 2, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 3, 2, 2, 2, 2, 2, 3, 2, 2, 2, 2, 3, 3, 2, 3, 2, 2, 2, 2, 3, 3, 3, 2, 2, 3, 2, 3, 3, 1, 2, 2, 3, 2, 2, 3, 2, 2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2, 2, 4, 3, 3, 2, 2, 2, 2, 2, 2, 1, 3, 2, 2, 2, 2, 1, 2, 2, 1, 2, 2, 3, 2, 3, 3, 3, 2, 3, 2, 1, 1, 3, 2, 2, 2, 1, 3, 1, 4, 1, 0, 1, 2, 3, 2, 2, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 3, 1, 3, 2, 1, 0, 2, 2, 1, 3, 1, 1, 1, 2, 2, 1, 1, 1, 3, 1, 1, 0, 4, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 0, 1, 2, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 1, 1, 2, 2, 2, 2, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 2, 1, 3, 1, 1, 1, 2, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 2, 0, 0, 0, 2, 1, 0, 2, 2, 0, 2, 0, 1, 1, 0, 2, 0, 0, 2, 1, 0, 3, 2, 0, 1, 2, 0, 2, 0, 1, 1, 2, 0, 1, 2, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#
#totmism13 <- c( 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 2, 3, 3, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 1, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 2, 1, 3, 1, 1, 2, 3, 1, 2, 2, 1, 2, 2, 3, 3, 2, 2, 2, 1, 2, 0, 1, 2, 2, 3, 3, 2, 1, 1, 1, 4, 3, 2, 3, 3, 2, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 3, 3, 3, 3, 3, 3, 3, 5, 5, 5, 4, 4, 4, 5, 4, 5, 5, 5, 5, 4, 5, 4, 5, 5, 4, 5, 3, 5, 4, 4, 4, 5, 5, 5, 2, 2, 2, 2, 2, 1, 2, 1, 2, 1, 1, 3, 1, 1, 3, 3, 1, 3, 1, 2, 1, 2, 3, 2, 1, 3, 1, 1, 3, 3, 2, 2, 3, 1, 1, 1, 1, 2, 3, 2, 2, 3, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 6, 6, 3, 4, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 2, 3, 2, 2, 3, 4, 2, 2, 2, 3, 3, 3, 3, 3, 4, 2, 3, 3, 2, 2, 3, 3, 2, 3, 5, 2, 3, 2, 3, 3, 3, 3, 3, 2, 3, 3, 2, 3, 3, 2, 3, 3, 2, 3, 2, 4, 3, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 4, 3, 3, 4, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 3, 1, 2, 1, 2, 2, 2, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 1, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 1, 2, 3, 2, 2, 1, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 2, 3, 3, 2, 2, 2, 2, 2, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 3, 2, 3, 3, 3, 2, 1, 2, 2, 2, 2, 1, 2, 1, 1, 1, 1, 1, 1, 2, 2, 2, 1, 2, 2, 1, 1, 2, 2, 1, 2, 3, 1, 1, 2, 2, 3, 3, 1, 2, 3, 1, 2, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 2, 3, 4, 3, 3, 4, 3, 3, 3, 3, 3)
#
#totmism14 <- c(1, 1, 0, 0, 0, 0, 1, 2, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 1, 2, 3, 2, 3, 2, 2, 1, 2, 3, 4, 4, 3, 2, 3, 2, 3, 4, 2, 3, 2, 2, 3, 3, 3, 4, 2, 3, 4, 2, 3, 3, 2, 2, 3, 2, 4, 2, 2, 2, 3, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 1, 0, 2, 1, 1, 1, 0, 5, 0, 0, 0, 1, 2, 2, 1, 0, 2, 0, 0, 1, 0, 5, 1, 0, 1, 0, 0, 3, 4, 3, 4, 3, 3, 2, 4, 4, 2, 2, 3, 3, 2, 2, 2, 2, 2, 2, 3, 3, 4, 3, 4, 2, 2, 3, 2, 3, 3, 3, 3, 2, 2, 4, 2, 3, 3, 2, 4, 3, 2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 4, 4, 4, 4, 4, 4, 4, 5, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 5, 4, 4, 4, 4, 4, 4, 4, 4, 5, 3, 3, 5, 1, 1, 3, 3, 3, 3, 3, 0, 2, 0, 0, 0, 0, 0, 2, 0, 0, 2, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 2, 0, 1, 0, 0, 0, 3, 1, 3, 1, 3, 3, 2, 2, 1, 0, 0, 1, 2, 0, 1, 1, 0, 1, 2, 1, 1, 1, 1, 1, 0, 1, 1, 0, 2, 1, 1, 3, 1, 0, 2, 0, 3, 3, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 5, 4, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 2, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 2, 1, 1, 1, 0, 1, 2, 0, 2, 0, 2, 2, 1, 1, 0, 2, 2, 1, 3, 2, 1, 1, 2, 1, 2, 3, 2, 0, 1, 0, 2, 1, 2, 2, 3, 1, 0, 1, 1, 2, 1, 2, 2, 2, 1, 0, 2, 1, 3, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 3, 1, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 2, 3, 2, 2, 3, 3, 3, 3, 3, 0, 1, 0, 0, 0, 0, 0, 4, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 3, 1, 0, 0, 0, 0, 4, 4, 4, 4, 4, 4, 3, 3, 4, 3, 3, 3, 4, 2, 3, 3, 3, 3, 4, 4, 4, 4, 3, 4, 3, 3, 4, 3, 3, 3, 4, 4, 3, 3, 3, 3, 4, 4, 3, 4, 4, 3, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 3, 3, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 2, 2, 2, 1, 1, 1, 2, 2, 0, 2, 1, 1, 1, 1, 1, 1, 1, 2, 1, 2)
#               
#totmism15 <- c( 2, 2, 1, 1, 0, 1, 3, 3, 3, 3, 3, 3, 2, 2, 3, 3, 3, 2, 3, 3, 2, 2, 3, 2, 4, 4, 3, 3, 2, 3, 3, 2, 3, 3, 2, 2, 3, 3, 2, 3, 2, 3, 3, 3, 2, 3, 3, 2, 3, 3, 3, 3, 4, 3, 4, 3, 3, 3, 4, 4, 4, 4, 4, 4, 3, 3, 3, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 3, 3, 4, 3, 3, 5, 5, 4, 3, 3, 2, 2, 2, 2, 2, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 2, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 1, 4, 4, 3, 2, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 4, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 2, 3, 2, 3, 2, 2, 3, 3, 3, 3, 3, 4, 2, 4, 4, 4, 3, 4, 3, 3, 3, 3, 4, 3, 3, 4, 3, 3, 3, 4, 3, 2, 4, 3, 3, 3, 3, 2, 4, 3, 2, 4, 3, 3, 3, 1, 2, 1, 1, 1, 3, 1, 2, 2, 2, 2, 1, 0, 1, 1, 1, 2, 2, 2, 1, 0, 1, 2, 2, 1, 2, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 4, 2, 2, 4, 4, 3, 2, 3, 2, 2, 3, 2, 3, 4, 3, 3, 3, 4, 2, 3, 3, 3, 3, 3, 4, 3, 4, 4, 2, 3, 2, 3, 2, 3, 2, 3, 3, 3, 2, 3, 2, 5, 2, 2, 2, 1, 5, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 3, 2, 3, 2, 2, 2, 2, 2, 2, 3, 2, 2, 2, 1, 2, 3, 2, 2, 5, 2, 2, 2, 3, 2, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 2, 4, 4, 2)
#
#
# totmism <- c(totmism1,totmism2,totmism3,totmism4,totmism5,totmism6,totmism7,totmism8,totmism9,totmism10,totmism11,totmism12,totmism13,totmism14,totmism15)
# hist(totmism)




#
#
#
#pname = 'Deyde_et_al.2006'
#for primer in primerstype[pname] :
#        if primer != 'type' :
#            for clade in primerstype[pname][primer] :
#                if clade[:15] == 'MismatchCounts_' : 
#                    print primer , primerstype[pname][primer][clade] 
#
