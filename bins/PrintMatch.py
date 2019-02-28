


import xlwt



# ###################################################################################################################################
# STAMPA SU FILE EXCELL I RISULTATI DEI PRIMER
# ###################################################################################################################################


def stampa_match(sheet, book_row, book_col , style , newtarget ,   newcounts  ,  expandedtarget,  dstart, dend) :
    #scrivo i risultati dei matches, colorandoli a differenza di quanto sono buoni :
    oldpattern = xlwt.Pattern()
    
    style0 = xlwt.easyxf('font: color black; alignment: horizontal center; ')      # (in font)  bold 1, 
    
    
    style1 = xlwt.easyxf('pattern: pattern solid; font: bold 1, color black; alignment: horizontal center; ') 
    style1.pattern.pattern_fore_colour = 50

    style2 = xlwt.easyxf('pattern: pattern solid; font: bold 1, color black; alignment: horizontal center; ') 
    style2.pattern.pattern_fore_colour = 51

    style3 = xlwt.easyxf('pattern: pattern solid; font: bold 1, color black; alignment: horizontal center; ') 
    style3.pattern.pattern_fore_colour = 52

    style4 = xlwt.easyxf('pattern: pattern solid; font: bold 1, color black; alignment: horizontal center; ') 
    style4.pattern.pattern_fore_colour = 53

    style5 = xlwt.easyxf('pattern: pattern solid; font: bold 1, color black; alignment: horizontal center; ') 
    style5.pattern.pattern_fore_colour = 2
    
    
    for j in range(len(expandedtarget)) :
        d = expandedtarget[j]
        e = d.replace('-', '')
        f = e.replace('X', '')
        q = f.replace('i', '')
        
        if  dstart <= j <  len(expandedtarget) - dend :
            i = j - dstart
            c = newtarget[i]
            
#            g = c.replace('-', '')
#            h = g.replace('X', '')
#            
#            m = h if h != '' else q
            
            if c == '.' :
                sheet.write(book_row, book_col, q, style0)
                book_col += 1
            else :
                if  85 <= newcounts[i] <= 100 :
                    sheet.write(book_row, book_col, q, style1) 
                    book_col += 1
    
                    
                elif  70 <= newcounts[i] < 85 :
                    sheet.write(book_row, book_col, q, style2)
                    book_col += 1
    
                
                elif  55 <= newcounts[i] < 70 :
                    sheet.write(book_row, book_col, q, style3)
                    book_col += 1
    
                
                elif  40 <= newcounts[i] < 55 :  
                    sheet.write(book_row, book_col, q, style4)
                    book_col += 1
    
                
                elif  0 <= newcounts[i] < 45 :
                    sheet.write(book_row, book_col, q, style5)
                    book_col += 1
        else :
    
            sheet.write(book_row, book_col, q, style0)
            book_col += 1
    
    style.pattern = oldpattern
    style = xlwt.easyxf('font: bold 1, color black; ')

    
    return style, book_row, book_col 

# ###################################################################################################################################
