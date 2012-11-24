## http://docs.python.org/howto/regex.html
## metacharacters
## . ^ $ * + ? { } [ ] \ | ( )
## . any character
## ^
## $
## * ## multiple
## +
## ?
## {
## }
## [ list
## ] list
## \
## |
## ( group
## ) group

def main():

    import urllib, re, os

    if not os.path.isfile('alice.txt'):
        urllib.urlretrieve('http://www.upriss.org.uk/python/alice.txt','alice.txt',)

    fd = open('alice.txt','r')
    lines = fd.readlines()
    fd.close()

    pattern = r"the "
    prog = keyword = re.compile(pattern, re.I)
    for line in lines:
        result = keyword.search(line)
##        match = re.search(r"the ", line, re.I)
        if result:
            print result.group(), ':', line
##        if match:
##            print result.group(), ':', line

    ## 2.1 Retrieve lines that have two consecutive o's.
    keyword = re.compile(r'oo')
    for line in lines:
        result = keyword.search(line)
        if result:
            print '2.1', line,

    ##2.2 Retrieve lines that contain a three letter string consisting of "s", then any character, then "e", such as "she".
    keyword = re.compile(r's.e')
    for line in lines:
        match = result = keyword.search(line)
        if result:
            print '2.2', match.group(), line,
    
    ##2.3 Retrieve lines with a three letter word that starts with s and ends with e.
    keyword = re.compile(r'\bs.e\b')
    for line in lines:
        match = result = keyword.search(line)
        if result:
            print '2.3', match.group(), line,

    ##2.4 Retrieve lines that contain a word of any length that starts with s and ends with e. Modify this so that the word has at least four characters.
    keyword = re.compile(r'\bs\w*e\b')
    for line in lines:
        match = result = keyword.search(line)
        if result:
            print '2.4', match.group(), ':', line,
    keyword = re.compile(r'\bs\w\w+e\b')
    for line in lines:
        match = result = keyword.search(line)
        if result:
            print '2.4b', match.group(), ':', line,
    
    ##2.5 Retrieve lines that start with a. Retrieve lines that start with a and end with n.
    keyword = re.compile(r'^a')
    for line in lines:
        match = result = keyword.search(line)
        if result:
            print '2.5a', match.group(), ':', line,
    keyword = re.compile(r'^a.*n$')
    for line in lines:
        match = result = keyword.search(line)
        if result:
            print '2.5b', match.group(), ':', line,

    ##2.6 Retrieve blank lines. Think of at least two ways of doing this.
    keyword = re.compile(r'^$')
    for line in lines:
        match = result = keyword.search(line)
        if result:
            print '2.6', match.group(), ':', line,

    ##2.7 Retrieve lines that do not contain the blank space character.
    keyword = re.compile(r' ')
    for line in lines:
        match = result = keyword.search(line)
        if not match:
            print '2.7', line,

    ##2.8 Retrieve lines that contain more than one blank space character.
    keyword = re.compile(r' .* ')
    for line in lines:
        match = result = keyword.search(line)
        if not match:
            print '2.8', line,
    ##
    ##3 Add a few lines with numbers etc. to the end of the alice.txt file so that you can search for the following regular expressions:
    ##
    ##3.1 an odd digit followed by an even digit (eg. 12 or 74) 
    keyword = re.compile(r'[13579][24680]')
    for line in lines:
        match = result = keyword.search(line)
        if match:
            print '3.1', match.group(), ':', line,
    ##3.2 a letter followed by a non-letter followed by a number 
    keyword = re.compile(r'[a-zA-Z][\W0-9_]\d')
    keyword = re.compile(r'[a-zA-Z][^a-zA-Z]\d')
    for line in lines:
        match = result = keyword.search(line)
        if match:
            print '3.2', match.group(), ':', line,
    ##3.3 a word that starts with an upper case letter 
    keyword = re.compile(r'\b[A-Z]\w*\b')
    for line in lines:
        match = result = keyword.search(line)
        if match:
            print '3.3', match.group(), ':', line,

    ##3.4 the word "yes" in any combination of upper and lower cases letters 
    keyword = re.compile(r'\byes\b', re.I)
    for line in lines:
        match = result = keyword.search(line)
        if match:
            print '3.4', match.group(), ':', line,
    ##3.5 one or more times the word "the" 
    keyword = re.compile(r'(the )+')
    for line in lines:
        match = result = keyword.search(line)
        if match:
            print '3.5', match.group(), ':', line,
    ##3.6 a date in the form of one or two digits, a dot, one or two digits, a dot, two digits 
    keyword = re.compile(r'\d\d?\.\d\d?\.\d\d')
    for line in lines:
        match = result = keyword.search(line)
        if match:
            print '3.6', match.group(), ':', line,
    ##3.7 a punctuation mark 
    keyword = re.compile(r'[\.,!\?\!:;]')
    for line in lines:
        match = result = keyword.search(line)
        if match:
            print '3.7', match.group(), ':', line,
    ##4.1 Write a script that asks users for their name, address and phone number. Test each input for accuracy, for example, there should be no letters in a phone number. A phone number should have a certain length. An address should have a certain format, etc. Ask the user to repeat the input in case your script identfies it as incorrect.
    ##
    ##4.2 Concerning your projects: what kind of checking is needed to ensure that users fill in the forms in a sensible manner? Make certain that your form can handle all kinds of input. For example, users can have several first names, middle initials, several last names (which may or may not be hyphenated).

    ##2.2 Optional: Print all lines in the alice.txt file so that the first and
    ##the last character in each line are switched.
    # compiling the regular expression:
    keyword = re.compile(r"(.)(.*)(.)")
    for line in lines:
        result = keyword.search(line)
        if result:
            print '2.2', result.group(3) + result.group(2) + result.group(1)    

    ##2.3 Print all lines in the alice.txt file that contain two double
    ##characters.
    keyword = re.compile(r"(.)\1(.*)(.)\3")
##    keyword = re.compile(r"(.)\1")
    for line in lines:
        result = keyword.search(line)
        if result:
            print '2.3', result.groups(), ':', line

    ##3.1 all upper case A by lower case a. 
    keyword = re.compile(r'A')
    for line in lines:
        result = keyword.search(line)
        if result:
            print '3.1', keyword.sub('a',line)

    ##3.2 Delete all words with more than 3 characters.
    keyword = re.compile(r'\b\w\w\w+\b')
    for line in lines:
        result = keyword.search(line)
        if result:
            print '3.2', keyword.sub('',line)

    ##3.3 Print two blank space characters after the "." at the end of a sentence.
    keyword = re.compile(r'\.')
    for line in lines:
        result = keyword.search(line)
        if result:
            print '3.3', keyword.sub('.  ',line)

    ##3.4 Replace single quotes (' or `) by double quotes. 
    keyword = re.compile(r'[\'`]')
    for line in lines:
        result = keyword.search(line)
        if result:
            print '3.4', keyword.sub('"',line)

        
    return

if __name__ == '__main__':
    main()
