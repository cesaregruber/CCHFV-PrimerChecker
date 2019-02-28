




# ###################################################################################################################################
# fa il reverse complement di una sequenza :


def revcomp(seq, complambidic ) :
    newseq = ''
    for c in seq[::-1] :
        newseq = newseq + complambidic[c]
    
    return newseq






# DIZIONARIO DEGLI AMBIGUITY CODES #####

# http://pythonforbiologists.com/index.php/introduction-to-python-for-biologists/regular-expressions/



def iupacdic():
    
    iname3 = "AmbiguityCodes.csv"
    
    ambidic = {}         #dizionario dei significati degli ambiguity codes          (es: ambidic['Y'] = 'CT' )
    ambinonambidic = {}  #dizionario dei significati degli ambiguity codes + tutti i caratteri 'normali'  + il cararrere 'gap'        (es: ambinonambidic['Y'] = 'CT' ; ambinonambidic['C'] = 'C' ; ambinonambidic['-'] = '-' )
    invambidic = {}     # dizionario inverso degli ambiguity codes                  (es: invambidic['AG'] = 'R' ;  invambidic['T'] = 'T' )
    complambidic = {}   # dizionario dei complementary reverses degli ambiguity codes   (es: complambidic['R'] = 'Y' ;  complambidic['C'] = 'G' )
    
    
    fi = open(iname3, 'r')
    #fi.readline()
    for line in fi :
        if 'IUPAC' not in line : 
            ls = line.split('\t')
            k = ls[0]
            means = ls[1].split(' or ')
            ambidic[k] =  ''.join(sorted(means))  #means
            ambinonambidic[k] = str(ambidic[k]) 
            invambidic[ambidic[k] ] = k
            complambidic[k] = ls[2][:-1]
    
    
    ambinonambidic['-'] = '-'
    
    
    
    ambikeys = ambidic.keys()
    fi.close()
    
    
    ambiset = set(ambidic.keys())
    bases =  ['A', 'C','G', 'T'] 
    combases = bases[::-1]
    for x in range(len(bases)) :
        invambidic[bases[x]] = bases[x]
        complambidic[bases[x]] = combases[x]
        ambinonambidic[bases[x]] = bases[x]
        
    return ambidic ,  ambinonambidic ,  invambidic ,  complambidic ,  ambikeys  ,  ambiset
        
    



# ###############################


