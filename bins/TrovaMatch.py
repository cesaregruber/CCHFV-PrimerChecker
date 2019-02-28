

import re
import IUPACfunc



# ###################################################################################################################################
#  TROVO I GENOMI DEL CLADE CHE MATCHANO AL 100 PER 100 con le sequenze NON allineate 


def trova_match(primer  , alnvirs,  gapsvirs,  simplevirs, RefSeq,  RefSeqNoGaps  ): #myclade,  primer,  s,  init ,  fine,  ambidic,  ambiset ) :
    
    # mi serve: 
   #1) la sequenza del primer       (STRINGA)
   #2) le sequenze dei virus allineati  (DIZIONARIO)
   #3) le posizioni dei gap negli allineamenti  (DIZIONARIO)  
   #4) le sequenze dei virus non allineati  (DIZIONARIO)
   #5) l'id della sequenza che si vuole usare come riferimento

    ambidic ,  ambinonambidic ,  invambidic ,  complambidic ,  ambikeys  ,  ambiset  = IUPACfunc.iupacdic()
    
    TOTALESATTI = 0
    primergaps = []
    init = 0
    fine = 0
    trueinit = 0 
    truefine=0
    gappedprimer = ''           # <- primer con i gap dentro
    trovato  = 0
    
    # TRADUCO GLI AMBIGUITY CODES  
    s = str(primer)         # ricopio il primer com' e' 
    for a in ambikeys:
        if a in s :
            s = s.replace(a,  '[' +ambidic[a] + ']' )            # <- 
    
    myvir = []  # [nome, sequenza]
    myseq = ''
    qualiesatti = []   # nome di tutte le sequenze esatte
    for vir in simplevirs :
        myseq =  str(simplevirs[vir].seq)
        if re.search( s,  myseq ):     # <- "s" e' la sequenza del primer tradotta con gli ambiguity codes
            TOTALESATTI += 1
            myvir = [vir, myseq]
            trovato = 1
            qualiesatti.append(vir)
            
    
    if TOTALESATTI == 0 :  # se non l'ho trovata
        #  #giochetto: provo a sostituire tutti i nucleotidi con gli ambiguity codes della corrispettiva transizione (purina = R, piramidina = Y)
        newprimer =  ''
        for i in range(len(primer)) :
            b = primer[i]
            if b == 'A' or b == 'G' :                  #purina
                newprimer = newprimer +'R'
            elif b == 'C' or b == 'T' :             #piramidina
                newprimer = newprimer + 'Y'
            
        # TRADUCO GLI AMBIGUITY CODES  
        s = str(newprimer)         # ricopio il primer com' e' 
        for a in ambikeys:
            if a in s :
                s = s.replace(a,  '[' +ambidic[a] + ']' )            # <- 
        
        myvir = []  # [nome, sequenza]
        myseq = ''
        for vir in simplevirs :
            myseq =  str(simplevirs[vir].seq)
            if re.search( s,  myseq ):     # <- "s" e' la sequenza del primer tradotta con gli ambiguity codes
                trovato = 1
                myvir = [vir, myseq]
        
        if trovato == 0 :
        
            print("il primer  non sembra coincidere con nessuna sequenza in database - neanche sostituendo gli ambiguity codes!")
            
    if trovato == 1 :

        m = re.search(s, myvir[1] )        
        
        init =  m.start()       # <- inizio e fine nella sequenza SEMPLICE (non allineata)
        fine = m.end()
        
        #alnseq = str(alnvirs[myvir[0]].seq)
        
        # cerco la posizione nell'allineamento : SE CI SONO GAPS segno la posizione del gap nel primer
        gaplist = gapsvirs[myvir[0]]
        
        
        
        # ORA init e fine diventano le posizioni nell'allineamento !
        for i in range(len(gaplist)) :
           if gaplist[i] <= init :
               init += 1
               fine += 1
           elif init < gaplist[i] < fine : 
               fine += 1
               puntogap = gaplist[i] - init
               primergaps.append(puntogap)
        
        if primergaps != [] :

            p = 0
            for i in range(fine - init) :
                if i in primergaps :
                    gappedprimer = gappedprimer + '-' 
                else :
                    gappedprimer = gappedprimer + primer[p]
                    p += 1
        
        else :
            gappedprimer = primer
        #
        
        

        
        primematch = RefSeq[init:fine].replace('-', '')    #sequenza che matcha col primer all'interno del riferimento
        resto =  RefSeqNoGaps.split(primematch)   #parti nel riferimento (senza gaps) che sono a sinistra e a destra del match
        
        trueinit = len(resto[0])+1
        truefine = len(RefSeqNoGaps) - len(resto[1])
            #
    
    return TOTALESATTI ,  init,  fine , primergaps,  trueinit ,  truefine,  gappedprimer , qualiesatti
