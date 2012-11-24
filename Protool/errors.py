#
# $Id: errors.py 3233 2008-02-06 15:25:33Z nielsen $
#
# This file is part of the Protool package
#
# (c) Jens Erik Nielsen, http://www.cmbi.kun.nl/gv/nielsen
##
# This file contains error messages and a few flags that control
# the overall behaviour of Protool
#
#



# Invalid atoms - f.ex. side chain atoms for a phi/psi evaluation
#
class ProtoolError(Exception):
    
    def __init__(self,txt=None):
        if txt:
            self.txt=txt
        else:
            self.txt='Protool error'
        return
        
    def __str__(self):
        return repr(self.txt)
        

class InvalidAtomError(ProtoolError):
    def __init__(self):
        self.txt='Invalid atoms for this operation'
        return

class HydrogenInTorsionError(ProtoolError):
    def __init__(self):
        self.txt='Hydrogen atom in torsion angle'
        return
                        

class NotAnAminoAcidError(ProtoolError):
    def __init__(self,residue=''):
        if residue=='':
            self.txt='Non-amino acid residue'
        else:
            self.txt='%s is not an amino acid' %residue
        return

        
# 'Not Found' Errors..
class AtomNotFoundError(ProtoolError):
    def __init__(self,uniqueid=''):
        self.txt='Atom Not Found: %s' %uniqueid
        self.id=uniqueid
        
class ResidueNotFoundError(ProtoolError):
    def __init(self,uniqueid):
        self.txt=='Residue Not Found: %s' %uniqueid
        self.id=uniqueid

        
ResidueNotFound='Residue Not Found'
NotFoundError='Entity Not Found'

# N- and C-terminals
class Nterm(ProtoolError):
    def __init__(self):
        self.txt='N-terminal Residue'
        
class Cterm(ProtoolError):
    def __init__(self):
        self.txt='C-terminal Residue'
        
# Incomplete Errors
#
IncompletePositionError='Not all coordinates were found this atom'

#
# I/O errors
EndOfFileError='End Of File'
FileNotFoundError='File Not Found'

#
# Flags
#
silent=1
