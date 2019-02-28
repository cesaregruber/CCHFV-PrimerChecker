


import IUPACfunc
from collections import Counter








# ###################################################################################################################################

def codifica_match( gappedprimer, myclade, init, fine ,  ambidic, invambidic, DELTA, firstORlast , LASTIMPORTANT  ) :

    #info sulle codifiche di ambiguita' :
    ambidic ,  ambinonambidic ,  invambidic ,  complambidic ,  ambikeys  ,  ambiset  = IUPACfunc.iupacdic()
    
    
    primer = [c for c in gappedprimer]        # primer   in formato lista con i gap 
    
    target = []
    expandedtarget = []   # salvo la sequenza con anche 3 nucleotidi prima e 3 dopo 
    counts = []       # frequenza dei nucleotidi nelle sequenze target
    
    expstart = init-DELTA if init >= DELTA else 0 
    expend = fine+DELTA if len(myclade[myclade.keys()[0]].seq) - fine >= DELTA else len(myclade[myclade.keys()[0]].seq) 
    
    dstart = init - expstart     # differenza di posizione fra inizio sequenza e inizio primer
    dend = expend - fine
    
    
    for pos in range(expstart, expend) :   # salvo la sequenza con anche 3 nucleotidi prima e 3 dopo 
        
        isgap = 0
        targetc = ''
        
        charlist = [ myclade[ky].seq[pos] for ky in myclade.keys() ] 
        charset = set(charlist)
        
        if '-' in charset :      # se c' e' un gap ...
            charset.remove('-') # ...lo tolgo .....
            isgap = 1
        
        if len(charset) > 0 :
            character = ''.join(charset)
            for c in  (ambiset & charset) :
                character = character.replace(c, ambidic[c])
            
            character = ''.join(sorted(set(character)))
            targetc = invambidic[character]
        
        if isgap == 1 :    # ...per poi rimetterlo.....
            targetc = '-' + targetc
        
        expandedtarget.append(targetc)
        
        if pos in range(init, fine) :            
            target.append(targetc)
            counts.append(dict(Counter(charlist)))    #conto le varianti
        
    newtarget = ['666']*len(target)
    
    #conto le volte che ci azzecco, base per base  :
    #countarget = []
    tot_sequenze = len(charlist)
    newcounts = [100.0]*len(counts)        
    
    for i in range(len(target)) :
        # # CASO NO GAP  ################################################################################
        if '-' not in target[i]  and '-' not in primer[i]  :     #se non ci sono gap in quella posizione
            # Casi Buoni
            if target[i] == primer[i]  :                                         # uguali normali
                newtarget[i] =  '.' 
            elif (primer[i] in ambidic)  and (target[i] not in ambidic) and  (target[i] in ambidic[primer[i]] )  :                  #        target semplice e primer ambigui
                newtarget[i] =  '.' 
            elif (primer[i] in ambidic ) and (target[i] in ambidic ) and ( ( set(ambidic[target[i]]) | set(ambidic[primer[i]]) ) == set(ambidic[primer[i]]) ) :          #           target ambiguo e primer ambiguo
                newtarget[i] =  '.' 
                
            # Casi No Buoni :
    
            elif  (primer[i] not in ambidic )  and (target[i] in ambidic)  :      # target ambiguo e primer semplice 
                myc = 0
                for b in ambidic[target[i]] :
                    if (b == primer[i] ) and (b in counts[i] ) :  # <- mi tengo conservativo: se c'e' una sequenza con ambiguita', considero come se certamente la sequenza non contenga lo stesso nucleotide del primer
                        myc += counts[i][b]
                newcounts[i] = (myc*100.0)/tot_sequenze
                newtarget[i] = target[i] 
                
            
            elif ( primer[i] in ambidic ) and ( target[i] in ambidic ) and ( ( set(ambidic[target[i]]) | set(ambidic[primer[i]]) ) != set(ambidic[primer[i]]) ) :        #           target ambiguo e primer ambiguo
                myc = 0  #conto il numero delle volte che ci ho azzeccato
                #myb =''          # controllo quali caratteri non sono nell'ambiguita' del primer 
                for b in ambidic[target[i]] :
                    if b in ambidic[primer[i]]  and b in counts[i] : 
                        myc += counts[i][b]
                if set(ambidic) & set(counts[i]) : # se pero' almeno una sequenza ha una ambiguita' (quindi in counts c'e' una ambiguita' ) 
                    for amb in set(ambidic) & set(counts[i]) :     #per ogni ambiguita' trovata 
                        if ambidic[amb] in ambidic[primer[i]]  :    #se l'ambiguita e' compresa in quella del primer :
                            myc += counts[i][amb]
                    #else :
                    #    myb = myb + b
                newcounts[i] = (myc*100.0)/tot_sequenze
                newtarget[i] = target[i] #invambidic[myb]
    
            else :                                                                                                               #           target semplice e primer qualunque, ma non "meccia"
                newtarget[i] =  target[i]
                newcounts[i] = 0.0
        
        
        
        # # CASO GAP 1 ################################################################################        
        elif '-' in target[i]  and '-'  not in  primer[i]  :        # se c'e' il gap nelle sequenze target MA NON NEL PRIMER
            
            #prima analizzo tutto quello senza gap
            target[i] = target[i][1:]

            # Casi Buoni
            if target[i] == primer[i] and  (primer[i] not in ambidic ) and (target[i] not in ambidic ):                                         # uguali normali
                newtarget[i] =  '.' 
                myc = counts[i][target[i]]
                newcounts[i] = (myc*100.0)/tot_sequenze
                
            elif target[i] != primer[i] and (primer[i] in ambidic)  and (target[i] not in ambidic) and  (target[i] in ambidic[primer[i]] )  :                  #        target semplice e primer ambigui
                newtarget[i] =  '.' 
                myc = counts[i][target[i]]
                newcounts[i] = (myc*100.0)/tot_sequenze
                
            elif target[i] != primer[i] and (primer[i] in ambidic ) and (target[i] in ambidic ) and  ( ( set(ambidic[target[i]]) | set(ambidic[primer[i]]) ) == set(ambidic[primer[i]]) ):          #           target ambiguo e primer ambiguo
                myc = 0  #conto il numero delle volte che ci ho azzeccato
                myb ='.'          # controllo quali caratteri non sono nell'ambiguita' del primer 
                for b in ambidic[target[i]] :
                    if b in ambidic[primer[i]]  and b in counts[i] : 
                        myc += counts[i][b]
                if set(ambidic) & set(counts[i]) : # se pero' almeno una sequenza ha una ambiguita' (quindi in counts c'e' una ambiguita' ) 
                    for amb in set(ambidic) & set(counts[i]) :     #per ogni ambiguita' trovata 
                        if ambidic[amb] in ambidic[primer[i]]  :    #se l'ambiguita e' compresa in quella del primer :
                            myc += counts[i][amb]
                
                newcounts[i] = (myc*100.0)/tot_sequenze
                newtarget[i] = myb
                
            # Casi No Buoni :
    
            elif  target[i] != primer[i] and (primer[i] not in ambidic )  and (target[i] in ambidic)  :      # target ambiguo e primer semplice 
                newtarget[i] = target[i] 
                myc = 0
                for b in ambidic[target[i]] :
                    if (b == primer[i] ) and (b in counts[i] ) :  # <- mi tengo conservativo: se c'e' una sequenza con ambiguita', considero come se certamente la sequenza non contenga lo stesso nucleotide del primer
                        myc += counts[i][b]
                newcounts[i] = (myc*100.0)/tot_sequenze
                newtarget[i] = target[i] 
                
            
            elif target[i] != primer[i] and ( primer[i] in ambidic ) and ( target[i] in ambidic ) and (  ( set(ambidic[target[i]]) | set(ambidic[primer[i]])  ) != set(ambidic[primer[i]])   ) :        #           target ambiguo e primer ambiguo
                myc = 0  #conto il numero delle volte che ci ho azzeccato
                #myb =''          # controllo quali caratteri non sono nell'ambiguita' del primer 
                for b in ambidic[target[i]] :
                    if b in ambidic[primer[i]]  and b in counts[i] : 
                        myc += counts[i][b]
                if set(ambidic) & set(counts[i]) : # se pero' almeno una sequenza ha una ambiguita' (quindi in counts c'e' una ambiguita' ) 
                    for amb in set(ambidic) & set(counts[i]) :     #per ogni ambiguita' trovata 
                        if ambidic[amb] in ambidic[primer[i]]  :    #se l'ambiguita e' compresa in quella del primer :
                            myc += counts[i][amb]
                            
                    #else :
                    #    myb = myb + b
                newcounts[i] = (myc*100.0)/tot_sequenze
                newtarget[i] = target[i] #invambidic[myb]
    
            else :                                                                                                               #           target semplice e primer qualunque, ma non "meccia"
                newtarget[i] = target[i]
                newcounts[i] = 0.0
            
            #poi aggiungo una 'X' che significa delezione
            newtarget[i] = 'X' if newtarget[i]  == '.' else newtarget[i] + 'X'    #          aggiungo un una 'X' che significa DELEZIONE

        # # CASO GAP 2 ################################################################################
        elif not '-' in target[i]  and '-'  in  primer[i]  :        # se NON c'e' il gap nelle sequenze target MA c'e' NEL PRIMER
        
            newtarget[i] = 'i' + target[i] + 'i'         # aggiungo una 'i' che significa inserzione
            newcounts[i] = 0.0
            

        # # CASO GAP 3 ################################################################################        
        elif '-' in target[i]  and '-'  in  primer[i]  :        # se c'e' il gap nelle sequenze del target e del primer
            
            #conto quanti gap ci sono nelle sequenze target:
            myc = counts[i]['-']
            newcounts[i] = (myc*100.0)/tot_sequenze
            
            # e salvo la stringa con l'inserzione :
            newtarget[i] = 'i' + target[i][1:] + 'i'  if len(target[i]) > 1 else '.'
            
    #confronto la stringa "gappedprimer" con le sequenze str(myclade[k].seq) per vedere quante varianti ci sono fra il primer e OGNI sequenza :
    strainlist = []     #lista delle sequenze analizzate
    numvarlist = []     #lista del numero di varianti fra sequenza e primer
    numvarlast = []  #lista del numero di varianti fra sequenza e primer NEGLI ULTIMI "LASTIMPORTANT" NUCLEOTIDI
    for  strain in myclade :     #ciclo su tutti gli strain
        mys =  str(myclade[strain].seq)[init:fine]
        #if 10*'-' not in mys : <- ELIMINATO!   # NON CONSIDERO QUELLE SEQUENZE CHE NON SONO STATE SEQUENZIATE IN QUELLA REGIONE !!!
        strainlist.append(strain)
        numvarlist.append(0)
        numvarlast.append(0)
        for i in range(len(gappedprimer) ): #ciclo su tutti i nucleotidi
            schar = ambinonambidic[mys[i]]                    #nucleotide nello strain !!!mi tengo conservativo: se nella sequenza c'e' un ambiguita', allora anche il primer deve avere quell'ambiguita'
            gchar = ambinonambidic[gappedprimer[i]]     #nucleotide nel primer
            if schar not in gchar :         # se il primer non "match"a  con la sequenza
                numvarlist[len(numvarlist) -1 ] += 1    # conto una mutazione in piu'
                if (i > len(gappedprimer)  - (LASTIMPORTANT + 1) ) and firstORlast == 'last'   :     # se mi trovo negli ultimi 5 nucleotidi e il primer e' forward
                    numvarlast[len(numvarlast) -1 ] += 1 # conto una mutazione in piu' anche negli ultimi 5 nucleotidi
                elif ( i < LASTIMPORTANT )  and firstORlast == 'first'   :     # se mi trovo negli ultimi 5 nucleotidi e il primer e' reverse
                    numvarlast[len(numvarlast) -1 ] += 1 # conto una mutazione in piu' anche negli ultimi 5 nucleotidi
        
    
    
    return newtarget ,  newcounts , expandedtarget , dstart, dend ,  strainlist ,  numvarlist , numvarlast
    
    
    
    
    
    
    
    
    
    
    
    
