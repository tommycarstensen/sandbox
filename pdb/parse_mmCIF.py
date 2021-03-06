## avoid doing slow len(l)!!!
## rewrite this script at some stage making use of regular expressions...

d_characters = {
    ## key = initiation
    ## value = termination
    "'":"'",
    '"':'"',
    '(':')',
    ';':';',
    }

bool_testmode = False

def read_lines(pdb):

##    fd = open('/data/mmCIF/%s/%s.cif' %(pdb[1:3],pdb,),'r')
##    fd = open('E:/mmCIF/%s/%s.cif' %(pdb[1:3],pdb,),'r')
##    fd = open('/media/Elements/mmCIF/%s/%s.cif' %(pdb[1:3],pdb,),'r')
##    fd = open('/home/tc/Downloads/mmCIF/%s/%s.cif' %(pdb[1:3],pdb,),'r')
##    fd = open('/media/Tommy/mmCIF/%s/%s.cif' %(pdb[1:3],pdb,),'r')
    fd = open('/media/WDMyBook1TB/2TB/mmCIF/%s/%s.cif' %(pdb[1:3],pdb,),'r')
    lines = fd.readlines()
    fd.close()

    return lines

def main(
    pdb,
    lines = None,
    bool_testmode = False,
    d_breaks = {}, ## break if e.g. _exptl.method == SOLUTION NMR
    l_data_categories_break = [], ## break if e.g. coordinates(_atom_site)
    l_data_categories_skip = [], ## skip if e.g. coordinates (_atom_site) and anisotropic factors (_atom_site_anisotrop)
    l_data_categories = [], ## parse selected data categories
    d_breaks_negation = {}, ## break if e.g. _pdbx_struct_assembly.oligomeric_details != monomeric
    ):

    if len(l_data_categories) > 0:
        for k in d_breaks.keys():
            l_data_categories += [k[:k.index('.')]]
        for k in d_breaks_negation.keys():
            l_data_categories += [k[:k.index('.')]]

    if not lines:
        lines = read_lines(pdb)

    d = main_loop(
        pdb,lines,
        d_breaks = d_breaks,
        d_breaks_negation = d_breaks_negation,
        l_data_categories_break = l_data_categories_break,
        l_data_categories_skip = l_data_categories_skip,
        l_data_categories = l_data_categories,
        )

    return d


def main_loop(
    pdb,lines,
    bool_testmode = False,
    d_breaks = {}, ## break if e.g. _exptl.method = SOLUTION NMR
    l_data_categories_break = [], ## break if e.g. coordinates(_atom_site)
    l_data_categories_skip = [], ## skip if e.g. coordinates (_atom_site) or anisotropic factors (_atom_site_anisotrop)
    l_data_categories = [],
    d_breaks_negation = {}, ## break if e.g. _pdbx_struct_assembly.oligomeric_details != monomeric
    ):

    d = {}

    i2 = 2
    i3 = i2
    for i2 in range(2,len(lines)):

        ## skip parsed data items
        if i2 < i3:
            continue

        ##
        ## loop over items
        ##
        if lines[i2].strip() == 'loop_':
            bool_loop = True
            ## loop over items
            l_items = []
            for i3 in range(i2+1,len(lines)):
                if lines[i3][0] == '_':
                    item = lines[i3].strip()
                    l_items += [item]
##                            if item in d.keys():
##                                print item
##                                stop
                    data_category = item[:item.index('.')]
                    if l_data_categories == None or l_data_categories == [] or data_category in l_data_categories:
                        d[item] = []
                        d[data_category] = []
                else:
                    break
            continue

        ##
        ## not a loop
        ##
        elif lines[i2][0] == '_':
            bool_loop = False
            item = lines[i2].split()[0]
            try:
                data_category = item[:item.index('.')]
            except:
                print pdb
                print item
                print lines[i2-1]
                print lines[i2]
                print lines[i2+1]
                stop
            ## only item on current line; not value
            if lines[i2].strip() == item:
                ## multiple lines
                if lines[i2+1][0] == ';':
                    value = ''
                    l_characters = []
                    for i3 in range(i2+1,len(lines)):
                        value += lines[i3].strip()
                        if lines[i3][0] == ';':
                            ## multiline value initiation
                            if l_characters == []:
                                l_characters += [';']
                                continue
                            ## multiline value termination
                            else:
                                l_characters = []
                                break
                ## single line
                else:
                    i3 = i2+1
                    value = lines[i3].strip()
            ## item and value on current line
            else:
                i3 = i2
                value = lines[i3][len(item):].strip()

            if (
                (value[0] == "'" and value[-1] == "'")
                or
                (value[0] == ";" and value[-1] == ";")
                ):
                value = value[1:-1]

            if l_data_categories and l_data_categories != [] and not data_category in l_data_categories:
                continue

            l_items = [item]
            l_values = [value]
            for index in range(len(l_items)):
                item = l_items[index]
                value = l_values[index]
                d[item] = [value]
            i3 += 1
            if d_breaks:
                if item in d_breaks.keys():
                    if (
                        value == d_breaks[item]
                        or
                        ## careful here, because "2" in "22"
                        value in d_breaks[item]
                        or
                        ## careful here, because "2" in "22"
                        d_breaks[item] in value
                        ):
                        return d
            if d_breaks_negation:
                if item in d_breaks_negation.keys():
                    if value != d_breaks_negation[item] and value not in d_breaks_negation[item]:
                        return d


        ##
        ## end of data category
        ##
##                elif lines[i2][0].strip() == '#':
        elif lines[i2][0] == '#':
            l_items = []

            ## if previous data category in l_data_categories_break
            ## then break after parsing that data category
            if data_category in l_data_categories_break:
                return d

            continue

        ##
        ## loop over values
        ##
        else:

            if l_data_categories and l_data_categories != [] and data_category not in l_data_categories:
                continue
            if data_category in l_data_categories_skip:
                continue
##            if data_category in l_data_categories_break:
##                return d

            ## loop over values
            l_values = []
            l_characters = []
            value = ''
            
            for i3 in range(i2,len(lines)):

##                print pdb, data_category, i3

                ## end of loop (data category)
                if lines[i3][0] == '#':
                    break

                ## last data item in data category
##                elif len(l_values) > len(l_items):
##                    stopstop
                elif len(l_values) == len(l_items):
                    break
##                elif pdb == '1aa1' and data_category == '_atom_site':
##                    print 'a', lines[i3]

                line = lines[i3].rstrip()
                l = lines[i3].rstrip().split("'")

                ## data item across multiple lines or with spaces
                if (
                    '"' in lines[i3] ## spaces
                    or
                    "'" in lines[i3] ## spaces
                    or
                    "(" in lines[i3]
                    or
                    ';' == lines[i3][0] ## multiline
                    or
                    l_characters != [] ## 1abw
                    ): ## 1wt1, 1a0a
##                            print lines[i3]
                    l_characters, value, l_values = special_character(
                        line,l_characters,lines,i3,l_values,value,
                        )

                ##  no continuation across lines
                elif len(l) == 1:
##                else:
                    l_values += lines[i3].split()

                else: ## get rid of this else and replace elif above with else!!!
                    print lines[i3]
                    print l
                    print l_items
                    stop3
                    for s in l:
                        if s == '': ## first character is '
                            continue
                        if s[0] == ' ' or s[-1] == ' ':
                            l_values += s.split()
                        else:
                            l_values += [s]

##                        if len(l_values) == len(l_items):
##                            break

            if bool_testmode == True:
                if len(l_items) != len(l_values):
                    print
                    print lines[i2]
                    print lines[i3]
                    print
                    print bool_loop
                    print
                    print len(l_items), 'items', l_items
                    print len(l_values), 'values', l_values
                    print
                    for i in range(min(len(l_items),len(l_values))):
                        print l_items[i], l_values[i]
                    print pdb
                    stop_diff_len

            for index in range(len(l_items)):
                item = l_items[index]
                try:
                    value = l_values[index]
                except:
                    print index, l_values
                    print l_items
                    stop

                if bool_testmode == True:
                    if (
                        (value[0] == "'" and value[-1] == "'")
                        or
                        (value[0] == ";" and value[-1] == ";")
                        ):
                        print value
                        stop
                        value = value[1:-1]

                d[item] += [value]

                if d_breaks:
                    if item in d_breaks.keys():
                        if (
                            value == d_breaks[item]
                            or
                            ## careful here, because "2" in "22"
                            value in d_breaks[item]
                            or
                            ## careful here, because "2" in "22"
                            d_breaks[item] in value
                            ):
                            return d
                if d_breaks_negation:
                    if item in d_breaks_negation.keys():
                        if value != d_breaks_negation[item] and value not in d_breaks_negation[item]:
                            return d

            data_category = l_items[0][:l_items[0].index('.')]
            d[data_category] += [l_values]

    return d


def special_character(line,l_characters,lines,i3,l_values,value,):
    
    i_character2 = 0
    for i_character1 in range(len(line)):

        if i_character1 < i_character2:
            continue

        ##
        ## space delimeter
        ##
        if line[i_character1] == ' ':
            continue

        ##
        ## multi-line value
        ##
        if line[i_character1] == ';':

            if bool_testmode == True:
                if i_character1 != 0:
                    print line
                    stop

            ## semicolon end
            if l_characters == [';']:
                if bool_testmode == True:
                    if line.strip() != ';':
                        print line
                        stop
                l_characters = []
                value += line[1:].strip()
                l_values += [value]
                value = ''
                break ## next line

            ## semicolon start
            elif l_characters == []:
                l_characters += [line[i_character1]]
                value += line[1:].strip()
                break ## next line

            else:
                print l_characters
                print line
                stop

        ##
        ## multi-line continuation
        ##
##        elif l_characters == [';'] and ';' not in line:
        elif l_characters == [';'] and line.strip() != ';': ## 1m9e
            value += line.strip()
            break ## next line

        ##
        ## quoted value
        ##
        elif line[i_character1] in ['"',"'",'(',]:

##                                    if len(l_characters) > 0 and line[i_character1] == d_characters[line[i_character1]]:
##                                        l_characters = []
##                                        print line
##                                        stop
##                                    elif len(l_characters) > 0:
##                                        stop
##                                    else:
            l_characters += [line[i_character1]]

            if bool_testmode == True:
                if lines[i3].strip() == ';':
                    print lines[i3-1]
                    print l_characters
                    print value
                    stop
                if len(l_characters) > 1:
                    print l_characters
                    stop

            ## quoted value
            for i_character2 in range(i_character1+1,len(line)):
                if line[i_character2] == d_characters[line[i_character1]]:
                    ## end of special character
                    if len(l_characters) == 1:
                        bool_break = False
                        if (
                            i_character2 == len(line)-1
                            or
                            line[i_character2+1] == ' '
                            ):
                            bool_break = True
                        if bool_break == True:
                            i_character2 += 1
                            value = line[i_character1+1:i_character2-1]
                            if bool_testmode == True:
                                if value.strip() == '' or value.strip() == ';':
                                    print value
                                    print line
                                    stop
                            l_values += [value]
                            l_characters = []
                            value = ''
                            break
                        else: ## e.g. 1y59 _chem_comp.pdbx_synonyms
                            continue
                    ## initiation of special character
                    elif len(l_characters) == 0:
                        l_characters = l_characters[:-1]
                elif bool_testmode == True and line[i_character2] == line[i_character1]:
                    l_characters += [line[i_character2]]
                    print line[i_character2]
                    print line
                    stop7

        ##
        ## non-quoted value
        ##
        else:
            if i_character1+1 == len(line):
                if bool_testmode == True:
                    if line[i_character1] == ' ':
                        stop
                value = line[i_character1]
                l_values += [value]
                value = ''
            for i_character2 in range(i_character1+1,len(line)):
                if line[i_character2] == ' ':
                    value = line[i_character1:i_character2]
                    l_values += [value]
                    value = ''
                    break
                elif i_character2+1 == len(line):
                    value = line[i_character1:i_character2+1]
                    l_values += [value]
                    value = ''
                    break

            ## break loop over i_character1
            if i_character2+1 == len(line):
                break

    return l_characters, value, l_values


if __name__ == '__main__':
    import time
    t1 = time.clock()
##    d = main('1y59',)
    d = main('1aon',)
    t2 = time.clock()
    print t2-t1
