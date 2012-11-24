'''this script has been run and errors fixed'''

import os

pdbpath = '/oxygenase_local/data/pdb/'

d_nucleotide = {}
d_protein = {}
d_hetatm = {}

l_residues = ['ALA', 'CYS', 'GLU', 'ASP', 'GLY', 'PHE', 'ILE', 'HIS', 'LYS', 'MET', 'LEU', 'ASN', 'GLN', 'PRO', 'SER', 'ARG', 'THR', 'TRP', 'VAL', 'TYR']
l_nucleotides = [
    'A','C','G','U','I', ## ribonucleotides
    'DA','DC','DG','DT','DI', ## deoxyribonucleotides
    'N', ## wild card
    ]

d_protein = {'1smj': {'D 460 ': {'HEM': ' ', 'SER': ' '}}, '1smi': {'B 463 ': {'HEM': ' ', 'THR': ' '}}, '1sml': {'A 269 ': {'ZN': ' ', 'ARG': ' '}}, '1lq8': {'G  14 ': {'ASP': ' ', 'NDG': ' '}}, '2idb': {'C 503 ': {'HIS': ' ', 'EDO': ' '}}, '1mcx': {'A   8 ': {'CA': ' ', 'LEU': ' '}}, '1j8u': {'A 425 ': {'ASP': ' ', 'FE2': ' '}}, '1j8t': {'A 425 ': {'ASP': ' ', 'FE2': ' '}}, '1ywh': {'O 313 ': {'BMA': ' ', 'THR': ' '}}, '2atf': {'A   1 ': {'NI': ' ', 'SER': ' '}}, '1f9c': {'A 373 ': {'MN': ' ', 'ARG': ' '}}, '1yq4': {'A   2 ': {'THR': ' ', '3NP': ' '}}, '2hhm': {'B   2 ': {'SO4': ' ', 'ALA': ' '}}, '1eh8': {'A 200 ': {'GLY': ' ', 'ZN': ' '}}, '1dao': {'H 340 ': {'FAB': ' ', 'LEU': ' '}}, '1eh6': {'A 200 ': {'GLY': ' ', 'ZN': ' '}}, '1ju3': {'A 575 ': {'PRO': ' ', 'PBC': ' '}}, '1vyw': {'C 300 ': {'292': ' ', 'ARG': ' '}}, '2i9z': {'A 102 ': {'EDO': ' ', 'ALA': ' '}}, '1hro': {'A   1 ': {'HEM': ' ', 'GLY': ' '}}, '2ovw': {'D 400 ': {'NAG': ' ', 'VAL': ' '}}, '2fer': {'A 417 ': {'MNR': ' ', 'HIS': ' '}}, '2feu': {'B 417 ': {'MNR': ' ', 'HIS': ' '}}, '1glb': {'G 500 ': {'GOL': ' ', 'ASP': ' '}}, '1qf8': {'B 182 ': {'GLN': ' ', 'HOH': ' '}}, '1piy': {'B 352 ': {'FE': ' ', 'GLU': ' '}}, '1piz': {'B 349 ': {'FE': ' ', 'GLN': ' '}}, '1mqf': {'A   1 ': {'HEM': ' ', 'MET': ' '}}, '1beb': {'B   1 ': {'SO4': ' ', 'LEU': ' '}}, '1a7k': {'D 360 ': {'MET': ' ', 'PO4': ' '}}, '1c09': {'C  54 ': {'FE': ' ', 'GLU': ' '}}, '1uiv': {'B1500 ': {'NI': ' ', 'VAL': ' '}}, '2fe6': {'A 417 ': {'MNR': ' ', 'HIS': ' '}}, '1qfk': {'L 152 ': {'GLC': ' ', 'ARG': ' '}}, '1ds6': {'A 192 ': {'MG': ' ', 'LEU': ' '}}, '1sfo': {'C 302 ': {'ZN': ' ', 'ASN': ' '}}, '1lm6': {'A 203 ': {'FE': ' ', 'GLU': ' '}}, '6gss': {'B   1 ': {'PRO': ' ', 'MES': ' '}}, '1ymm': {'A 188 ': {'NAG': ' ', 'GLU': ' '}}, '2exr': {'A   1 ': {'SER': ' ', 'FAD': ' '}}, '2ay1': {'B  30 ': {'GLN': ' ', 'HOH': ' '}}, '2ay2': {'B  30 ': {'GLN': ' ', 'HOH': ' '}}, '2hg0': {'A 403 ': {'FUL': ' ', 'HIS': ' '}}, '2ay4': {'B  30 ': {'GLN': ' ', 'HOH': ' '}}, '2ay5': {'B  30 ': {'GLN': ' ', 'HOH': ' '}}, '2ay6': {'B  30 ': {'GLN': ' ', 'HOH': ' '}}, '2ay7': {'B  30 ': {'GLN': ' ', 'HOH': ' '}}, '1ddo': {'H 341 ': {'ITR': ' ', 'THR': ' '}}, '2ay9': {'B  30 ': {'GLN': ' ', 'HOH': ' '}}, '1y6w': {'A 148 ': {'CA': ' ', 'LYS': ' '}}, '1wiy': {'B 389 ': {'HEM': ' ', 'ALA': ' '}}, '1lth': {'R   1 ': {'OXM': ' ', 'ALA': ' '}}, '1xre': {'B   2 ': {'MN': ' ', 'SER': ' '}}, '4cox': {'D 601 ': {'HEM': ' ', 'ARG': ' '}}, '2paj': {'A 485 ': {'ZN': ' ', 'GLU': ' '}}, '10gs': {'B   1 ': {'PRO': ' ', 'MES': ' '}}, '1y67': {'D 220 ': {'FE': ' ', 'ALA': ' '}}, '16gs': {'B   1 ': {'PRO': ' ', 'MES': ' '}}, '20gs': {'B   1 ': {'PRO': ' ', 'MES': ' '}}, '3sqc': {'C   9 ': {'PRO': ' ', 'HOH': ' '}}, '4sdh': {'B   1 ': {'HEM': ' ', 'PRO': ' '}}, '1mu5': {'A   2 ': {'CA': ' ', 'SER': ' '}}, '1b2y': {'A   1 ': {'PCA': ' ', 'GLN': ' '}}, '2sqc': {'B   3 ': {'C8E': ' ', 'GLU': ' '}}, '1oas': {'B 317 ': {'GLU': ' ', 'PLP': ' '}}, '2gsm': {'A   6 ': {'ILE': ' ', 'HOH': ' '}}, '1muc': {'B 373 ': {'MN': ' ', 'ARG': ' '}}, '1ggy': {'A 730 ': {'YB': ' ', 'SER': ' '}}, '1yro': {'D 125 ': {'MET': ' ', 'MN': ' '}}, '2gss': {'B   1 ': {'PRO': ' ', 'MES': ' '}}, '2g73': {'A 181 ': {'MG': ' ', 'LEU': ' '}}, '2g74': {'A 181 ': {'MG': ' ', 'LEU': ' '}}, '2gsz': {'E 362 ': {'ILE': ' ', 'SO4': ' '}}, '1kif': {'H 341 ': {'BEZ': ' ', 'THR': ' '}}, '1jpz': {'B 460 ': {'HEM': ' ', 'SER': ' '}}, '1mlw': {'A 400 ': {'FE': ' ', 'LYS': ' '}}, '1xvd': {'B   4 ': {'FE': ' ', 'SER': ' '}}, '1xvg': {'B   4 ': {'FE': ' ', 'SER': ' '}}, '2fta': {'D   1 ': {'CU': ' ', 'ALA': ' '}}, '2ogj': {'A 401 ': {'ZN': ' ', 'ARG': ' '}}, '5cox': {'D 601 ': {'HEM': ' ', 'ARG': ' '}}, '2eun': {'A   1 ': {'HEM': ' ', 'MET': ' '}}, '2euo': {'A   1 ': {'HEM': ' ', 'MET': ' '}}, '2eup': {'A   1 ': {'HEM': ' ', 'MET': ' '}}, '2euq': {'A   1 ': {'HEM': ' ', 'MET': ' '}}, '2eur': {'A   1 ': {'HEM': ' ', 'MET': ' '}}, '2eus': {'A   1 ': {'HEM': ' ', 'MET': ' '}}, '2eut': {'A   1 ': {'HEM': ' ', 'MET': ' '}}, '2euu': {'A   1 ': {'HEM': ' ', 'MET': ' '}}, '1azx': {'L   1 ': {'NTP': ' ', 'HIS': ' '}}, '2bc2': {'B   2 ': {'GLN': ' ', 'HOH': ' '}}, '2yu9': {'C 302 ': {'ZN': ' ', 'ASN': ' '}}, '2p0a': {'B 401 ': {'GLY': ' ', 'ANP': ' '}}, '3pgh': {'D 601 ': {'HEM': ' ', 'ARG': ' '}}, '1nm0': {'A   1 ': {'HEM': ' ', 'MET': ' '}}, '3pgm': {'B 233 ': {'SO4': ' ', 'ALA': ' '}}, '1eh7': {'A 200 ': {'GLY': ' ', 'ZN': ' '}}, '3sdh': {'B   1 ': {'HEM': ' ', 'PRO': ' '}}, '1fta': {'D   4 ': {'AMP': ' ', 'ALA': ' '}}, '4blc': {'D   1 ': {'HEM': ' ', 'ALA': ' '}}, '2nvx': {'C 302 ': {'ZN': ' ', 'ASN': ' '}}, '2nvv': {'F 501 ': {'ZN': ' ', 'HIS': ' '}}, '2nvt': {'C 302 ': {'ZN': ' ', 'ASN': ' '}}, '2nvr': {'C 506 ': {'K': ' ', 'LEU': ' '}}, '2nvq': {'C 302 ': {'ZN': ' ', 'ASN': ' '}}, '4gss': {'B   1 ': {'GTX': ' ', 'PRO': ' '}}, '1hh8': {'A 201 ': {'FLC': ' ', 'LYS': ' '}}, '3sdp': {'B 191 ': {'FE': ' ', 'LYS': ' '}}, '1tv3': {'A   1 ': {'MET': ' ', 'BG4': ' '}}, '1tv2': {'A   1 ': {'MET': ' ', 'BG5': ' '}}, '2ani': {'A 320 ': {'FE': ' ', 'GLU': ' '}}, '1gcw': {'D 135 ': {'CMO': ' ', 'TYR': ' '}}, '1bd0': {'A   1 ': {'IN5': ' ', 'MET': ' '}}, '17gs': {'B   1 ': {'GTX': ' ', 'PRO': ' '}}, '1ggz': {'A 148 ': {'CA': ' ', 'LYS': ' '}}, '2iga': {'B   1 ': {'MET': ' ', 'CA': ' '}}, '1qnl': {'A   1 ': {'MET': ' ', 'BMD': ' '}}, '2a8a': {'A 439 ': {'LYS': ' ', 'CD': ' '}}, '2aa2': {'A 984 ': {'LYS': ' ', 'BOG': ' '}}, '2a8h': {'B 477 ': {'ZN': ' ', 'VAL': ' '}}, '2h6a': {'B 318 ': {'ZN': ' ', 'ALA': ' '}}, '2gn7': {'B 241 ': {'ALA': ' ', 'MAN': ' '}}, '1bqh': {'I   3 ': {'NAG': ' ', 'GLN': ' '}}, '1e0e': {'B  47 ': {'GLY': ' ', 'ZN': ' '}}, '1lld': {'B   1 ': {'NAD': ' ', 'ALA': ' '}}, '1zt4': {'A   1 ': {'ALA': ' ', 'AGH': ' '}}, '2axr': {'A 492 ': {'HIS': ' ', 'NAG': ' '}}, '1bc2': {'B   3 ': {'SO4': ' ', 'LYS': ' '}}, '1afs': {'B 321 ': {'ASP': ' ', 'TES': ' '}}, '1ids': {'D 200 ': {'FE': ' ', 'GLN': ' '}}, '1nik': {'C 302 ': {'ZN': ' ', 'ASN': ' '}}, '1l0v': {'M 601 ': {'GLY': ' ', 'FAD': ' '}}, '2gnd': {'B 241 ': {'ALA': ' ', 'MAN': ' '}}, '2gnm': {'A 242 ': {'GLN': ' ', 'MAN': ' '}}, '2min': {'C 492 ': {'CA': ' ', 'ALA': ' '}}, '1jaf': {'B 129 ': {'HEM': ' ', 'LYS': ' '}}, '1rqq': {'F 101 ': {'LYS': ' ', '112': ' '}}, '2ccy': {'B   1 ': {'HEM': ' ', 'GLN': ' '}}, '1fwq': {'A   1 ': {'MET': ' ', 'ZN': ' '}}, '1cx2': {'D 601 ': {'HEM': ' ', 'ARG': ' '}}, '1qj3': {'B 429 ': {'GLN': ' ', 'PLP': ' '}}, '1sky': {'B   1 ': {'MET': ' ', 'SO4': ' '}}, '1fpi': {'B   8 ': {'AMP': ' ', 'THR': ' '}}, '1fpk': {'B   6 ': {'TL': ' ', 'PHE': ' '}}, '1fpj': {'B   8 ': {'AMP': ' ', 'THR': ' '}}, '1fpl': {'B   8 ': {'AMP': ' ', 'THR': ' '}}, '1xu3': {'B   4 ': {'FE': ' ', 'SER': ' '}}, '1fpe': {'B   4 ': {'AMP': ' ', 'ALA': ' '}}, '1fpd': {'B   4 ': {'AMP': ' ', 'ALA': ' '}}, '1fpg': {'B   4 ': {'AMP': ' ', 'ALA': ' '}}, '1fpf': {'B   4 ': {'AMP': ' ', 'ALA': ' '}}, '8cat': {'B   2 ': {'NDP': ' ', 'ASP': ' '}}, '1jwb': {'B   1 ': {'MET': ' ', 'ZN': ' '}}, '2hsi': {'A   2 ': {'ZN': ' ', 'PRO': ' '}}, '2r2f': {'B 288 ': {'FEO': ' ', 'PRO': ' '}}, '1fpp': {'A   1 ': {'MET': ' ', 'FPP': ' '}}, '1g20': {'C 492 ': {'CA': ' ', 'ALA': ' '}}, '2nvz': {'C 302 ': {'ZN': ' ', 'ASN': ' '}}, '1ime': {'B   2 ': {'CA': ' ', 'ALA': ' '}}, '1xuq': {'A   2 ': {'MN': ' ', 'ALA': ' '}}, '2hmj': {'A 449 ': {'FES': ' ', 'ARG': ' '}}, '2hmm': {'A 449 ': {'FES': ' ', 'ARG': ' '}}, '2hml': {'A 449 ': {'FES': ' ', 'ARG': ' '}}, '1imc': {'A   4 ': {'PRO': ' ', 'CL': ' '}}, '1imb': {'B   4 ': {'PRO': ' ', 'LIP': ' '}}, '2ij4': {'B 463 ': {'HEM': ' ', 'THR': ' '}}, '1hkc': {'A 917 ': {'PO4': ' ', 'SER': ' '}}, '1q55': {'D 611 ': {'CA': ' ', 'LEU': ' '}}, '1tlg': {'A   1 ': {'MET': ' ', 'GAL': ' '}}, '2pms': {'B 337 ': {'SO4': ' ', 'GLU': ' '}}, '2ns9': {'B   1 ': {'PO4': ' ', 'LEU': ' '}}, '1r5u': {'C 302 ': {'ZN': ' ', 'ASN': ' '}}, '1ksa': {'B   8 ': {'VAL': ' ', 'HOH': ' '}}, '2ox4': {'H 401 ': {'MG': ' ', 'HIS': ' '}}, '1rgq': {'B 189 ': {'AKP': ' ', 'HIS': ' '}}, '1amu': {'B 551 ': {'AMP': ' ', 'TRP': ' '}}, '2qgp': {'C 101 ': {'ZN': ' ', 'LYS': ' '}}, '1q5a': {'B 611 ': {'CA': ' ', 'LEU': ' '}}, '1q5c': {'D 611 ': {'CA': ' ', 'LEU': ' '}}, '1q5b': {'C 611 ': {'CA': ' ', 'LEU': ' '}}, '1t1i': {'A 358 ': {'CA': ' ', 'SER': ' '}}, '7cat': {'B   2 ': {'NDP': ' ', 'ASP': ' '}}, '1sdd': {'B2183 ': {'TYR': ' ', 'CU': ' '}}, '1t1g': {'A 358 ': {'CA': ' ', 'SER': ' '}}, '2p8e': {'B 301 ': {'MG': ' ', 'HIS': ' '}}, '1bls': {'B 361 ': {'IPP': ' ', 'GLN': ' '}}, '2aio': {'A 312 ': {'ZN': ' ', 'ALA': ' '}}, '1ke1': {'A   1 ': {'CA': ' ', 'ALA': ' '}}, '1pyg': {'D 841 ': {'AMP': ' ', 'ILE': ' '}}, '2iid': {'D 490 ': {'ILE': ' ', 'FAD': ' '}}, '1k9i': {'A 403 ': {'CA': ' ', 'PRO': ' '}}, '1tqn': {'A 500 ': {'HEM': ' ', 'VAL': ' '}}, '1hkb': {'B 917 ': {'GLC': ' ', 'SER': ' '}}, '1n5k': {'B 213 ': {'MG': ' ', 'PRO': ' '}}, '2h8h': {'A 530 ': {'GLY': ' ', 'H8H': ' '}}, '2h34': {'A   2 ': {'NA': ' ', 'ASP': ' '}}, '1n5l': {'B 210 ': {'MG': ' ', 'PRO': ' '}}, '1jk7': {'A 303 ': {'LYS': ' ', 'BME': ' '}}, '1tqs': {'A   1 ': {'NAG': ' ', 'ARG': ' '}}, '1tqt': {'A   1 ': {'NAG': ' ', 'ARG': ' '}}, '1tqu': {'A   1 ': {'NDG': ' ', 'ARG': ' '}}, '1tqv': {'A   1 ': {'NDG': ' ', 'ARG': ' '}}, '1tqw': {'A   1 ': {'NAG': ' ', 'ARG': ' '}}, '1f7r': {'A 135 ': {'MG': ' ', 'SER': ' '}}, '1f7n': {'A 120 ': {'MG': ' ', 'SER': ' '}}, '2ofp': {'B 301 ': {'NAP': ' ', 'ARG': ' '}}, '1nlt': {'A 350 ': {'ZN': ' ', 'GLU': ' '}}, '1n4q': {'L  16 ': {'ZN': ' ', 'GLU': ' '}}, '1dcm': {'B 124 ': {'MG': ' ', 'ALA': ' '}}, '2gmr': {'M 303 ': {'BCL': ' ', 'MET': ' '}}, '2i8e': {'A 101 ': {'TRP': ' ', 'IOD': ' '}}, '1x9u': {'B 105 ': {'GLY': ' ', 'CU': ' '}}, '1byf': {'B 125 ': {'GOL': ' ', 'ASP': ' '}}, '1x9r': {'B 105 ': {'GLY': ' ', 'CU': ' '}}, '2nw9': {'A 301 ': {'HIS': ' ', 'FT6': ' '}}, '2nw8': {'B 305 ': {'MN': ' ', 'HIS': ' '}}, '1gnw': {'B   1 ': {'GTX': ' ', 'ALA': ' '}}, '1kwp': {'B  14 ': {'HG': ' ', 'PHE': ' '}}, '1azs': {'C 395 ': {'GSP': ' ', 'GLY': ' '}}, '1azt': {'B 395 ': {'GSP': ' ', 'GLY': ' '}}, '1r9t': {'C 302 ': {'ZN': ' ', 'ASN': ' '}}, '14gs': {'B   1 ': {'PRO': ' ', 'MES': ' '}}, '1mwo': {'A 435 ': {'GLY': ' ', 'ZN': ' '}}, '1e2y': {'J 180 ': {'GLY': ' ', 'CL': ' '}}, '1rly': {'A  60 ': {'ZN': ' ', 'ALA': ' '}}, '6cox': {'B 601 ': {'HEM': ' ', 'ARG': ' '}}, '2a5k': {'B 305 ': {'AZP': ' ', 'PHE': ' '}}, '3bls': {'B 361 ': {'APB': ' ', 'GLN': ' '}}, '1gn3': {'B 200 ': {'FE': ' ', 'GLN': ' '}}, '1kw0': {'A 425 ': {'ASP': ' ', 'FE2': ' '}}, '1nqg': {'A   5 ': {'CA': ' ', 'PRO': ' '}}, '1nqe': {'A   1 ': {'MG': ' ', 'GLN': ' '}}, '2np5': {'B 201 ': {'LMT': ' ', 'LEU': ' '}}, '1nqh': {'A   4 ': {'CA': ' ', 'SER': ' '}}, '2hkk': {'A 260 ': {'ZN': ' ', 'LYS': ' '}}, '2are': {'B 251 ': {'GLU': ' ', 'MAN': ' '}}, '1ybv': {'B   2 ': {'BEA': ' ', 'ALA': ' '}}, '1zr6': {'A 492 ': {'HIS': ' ', 'NAG': ' '}}, '2p27': {'A 300 ': {'GLY': ' ', 'MG': ' '}}, '2hk6': {'A   1 ': {'MET': ' ', 'HOH': ' '}}, '1ieb': {'B   3 ': {'NAG': ' ', 'SER': ' '}}, '1fat': {'D 236 ': {'CA': ' ', 'THR': ' '}}, '1gxf': {'B 492 ': {'LEU': ' ', 'FAD': ' '}}, '2h7s': {'C   1 ': {'HEM': ' ', 'THR': ' '}}, '2h7r': {'A   1 ': {'HEM': ' ', 'THR': ' '}}, '2h7q': {'A   1 ': {'HEM': ' ', 'THR': ' '}}, '1fah': {'B 471 ': {'HOH': ' ', 'ARG': ' '}}, '1fag': {'D 465 ': {'PAM': ' ', 'GLN': ' '}}, '1oqj': {'B 179 ': {'ZN': ' ', 'ARG': ' '}}, '1oqm': {'D 125 ': {'MET': ' ', 'MN': ' '}}, '1oqc': {'D   1 ': {'MET': ' ', 'FAD': ' '}}, '3gss': {'B   1 ': {'PRO': ' ', 'MES': ' '}}, '1hrp': {'B   1 ': {'NAG': ' ', 'SER': ' '}}, '2ar3': {'C 160 ': {'ZN': ' ', 'ALA': ' '}}, '3min': {'C 492 ': {'CA': ' ', 'ALA': ' '}}, '1e8a': {'B  90 ': {'CA': ' ', 'LYS': ' '}}, '1uaz': {'B 250 ': {'RET': ' ', 'ALA': ' '}}, '1dps': {'L  12 ': {'NA': ' ', 'THR': ' '}}, '1q2l': {'A 961 ': {'ZN': ' ', 'ASN': ' '}}, '1odb': {'F  91 ': {'GLU': ' ', 'CU': ' '}}, '1zoa': {'B 461 ': {'HEM': ' ', 'PRO': ' '}}, '1jv4': {'A 160 ': {'ALA': ' ', 'CD': ' '}}, '2b76': {'M 601 ': {'GLY': ' ', 'FAD': ' '}}, '1zo9': {'B 461 ': {'HEM': ' ', 'PRO': ' '}}, '2b8k': {'J  66 ': {'ZN': ' ', 'LEU': ' '}}, '1zo4': {'B 461 ': {'HEM': ' ', 'PRO': ' '}}, '1e6v': {'D 553 ': {'F43': ' ', 'GLU': ' '}}, '2alu': {'A 687 ': {'CO3': ' ', 'LEU': ' '}}, '2alt': {'A 687 ': {'CO3': ' ', 'LEU': ' '}}, '1w0e': {'A 501 ': {'HEM': ' ', 'SER': ' '}}, '1w0f': {'A 501 ': {'HEM': ' ', 'SER': ' '}}, '1w0g': {'A 501 ': {'HEM': ' ', 'SER': ' '}}, '2iep': {'B 212 ': {'HOH': ' ', 'ALA': ' '}}, '1y0j': {'A 239 ': {'ZN': ' ', 'ARG': ' '}}, '1r33': {'A1045 ': {'NAG': ' ', 'SER': ' '}}, '1tyl': {'D   2 ': {'VAL': ' ', 'CL': ' '}}, '1tym': {'D   2 ': {'VAL': ' ', 'CL': ' '}}, '1kfy': {'M 601 ': {'GLY': ' ', 'FAD': ' '}}, '1dbr': {'D   4 ': {'MG': ' ', 'LYS': ' '}}, '1nkx': {'A 686 ': {'FE': ' ', 'PHE': ' '}}, '2hfz': {'A 902 ': {'ZN': ' ', 'ASP': ' '}}, '2p69': {'A 300 ': {'HIS': ' ', 'PLP': ' '}}, '9gss': {'B   1 ': {'PRO': ' ', 'MES': ' '}}, '1cdl': {'D   4 ': {'CA': ' ', 'LEU': ' '}}, '2q41': {'D 334 ': {'ASN': ' ', 'HOH': ' '}}, '1ovx': {'B  49 ': {'ZN': ' ', 'GLU': ' '}}, '1tfz': {'A 409 ': {'FE': ' ', 'SER': ' '}}, '1tfp': {'B   1 ': {'SO4': ' ', 'VAL': ' '}}, '2ay3': {'B  30 ': {'GLN': ' ', 'HOH': ' '}}, '1kf6': {'M 601 ': {'GLY': ' ', 'FAD': ' '}}, '1lr0': {'A 129 ': {'ZN': ' ', 'HIS': ' '}}, '2cwg': {'E   4 ': {'SIA': ' ', 'ALA': ' '}}, '2ay8': {'B  30 ': {'GLN': ' ', 'HOH': ' '}}, '11gs': {'B   1 ': {'PRO': ' ', 'GTT': ' '}}, '1lrm': {'A 425 ': {'ASP': ' ', 'FE': ' '}}, '2iss': {'C 282 ': {'5RP': ' ', 'GLU': ' '}}, '1frf': {'L 552 ': {'MG': ' ', 'GLU': ' '}}, '5gss': {'B   1 ': {'PRO': ' ', 'MES': ' '}}, '2e2i': {'C 302 ': {'ZN': ' ', 'ASN': ' '}}, '2e2h': {'C 302 ': {'ZN': ' ', 'ASN': ' '}}, '2e2j': {'C 302 ': {'ZN': ' ', 'ASN': ' '}}, '1khv': {'B   2 ': {'LU': ' ', 'SER': ' '}}, '2h98': {'A   2 ': {'GLU': ' ', 'CL': ' '}}, '1sie': {'F   4 ': {'SIA': ' ', 'ARG': ' '}}, '1sid': {'F   3 ': {'GLC': ' ', 'LYS': ' '}}, '1nyk': {'B 204 ': {'GLY': ' ', 'FES': ' '}}, '1pxx': {'A 601 ': {'HEM': ' ', 'ARG': ' '}}, '2gug': {'C 375 ': {'GLY': ' ', 'FMT': ' '}}, '1pxe': {'A   1 ': {'MET': ' ', 'ZN': ' '}}, '1yu1': {'A 381 ': {'GLU': ' ', 'MMC': ' '}}, '8ruc': {'G   1 ': {'MET': ' ', 'CAP': ' '}}, '1r9s': {'C 302 ': {'ZN': ' ', 'ASN': ' '}}, '1a25': {'B 289 ': {'CA': ' ', 'PRO': ' '}}, '1hef': {'I 204 ': {'OHE': ' ', 'PHE': ' '}}, '1heg': {'I 204 ': {'GLY': ' ', 'OHE': ' '}}, '1si1': {'A 317 ': {'FE': ' ', 'SER': ' '}}, '1si0': {'A 320 ': {'FE': ' ', 'LYS': ' '}}, '1imd': {'B   4 ': {'PRO': ' ', 'MN': ' '}}, '1jw9': {'B   1 ': {'MET': ' ', 'ZN': ' '}}, '1pqv': {'C 300 ': {'ZN': ' ', 'TYR': ' '}}, '2h9q': {'D   4 ': {'CCU': ' ', 'ARG': ' '}}, '1n4p': {'L  16 ': {'ZN': ' ', 'GLU': ' '}}, '1n4s': {'L  16 ': {'ZN': ' ', 'GLU': ' '}}, '1n4r': {'L  17 ': {'ZN': ' ', 'ARG': ' '}}, '1q74': {'D 300 ': {'ZN': ' ', 'ALA': ' '}}, '2h4k': {'A   1 ': {'MET': ' ', '509': ' '}}, '1gyl': {'A 360 ': {'FMN': ' ', 'GLY': ' '}}, '2h4g': {'A   1 ': {'MET': ' ', '694': ' '}}, '1lnq': {'H   8 ': {'CA': ' ', 'ILE': ' '}}, '2oyc': {'A 301 ': {'NA': ' ', 'HIS': ' '}}, '2hmo': {'A 449 ': {'FES': ' ', 'ARG': ' '}}, '2hmn': {'A 449 ': {'FES': ' ', 'ARG': ' '}}, '1hw7': {'A 240 ': {'ZN': ' ', 'ALA': ' '}}, '1jni': {'A 111 ': {'HEM': ' ', 'ILE': ' '}}, '1lgc': {'H 513 ': {'FUC': ' ', 'GLN': ' '}}, '1ady': {'D   1 ': {'MET': ' ', 'HAM': ' '}}, '1pj0': {'B 358 ': {'FE': ' ', 'VAL': ' '}}, '18gs': {'B   1 ': {'PRO': ' ', 'MES': ' '}}, '1ima': {'B   4 ': {'PRO': ' ', 'IPD': ' '}}, '1n97': {'B 389 ': {'HEM': ' ', 'ALA': ' '}}, '1bav': {'D 309 ': {'TYR': ' ', 'BIP': ' '}}, '1tkw': {'B 251 ': {'VAL': ' ', 'HEC': ' '}}, '1mz4': {'A 135 ': {'HEM': ' ', 'VAL': ' '}}, '4sbv': {'C   1 ': {'CA': ' ', 'ALA': ' '}}}
for pdb in d_protein.keys():
    for ID in d_protein[pdb].keys():
        hetIDs = ','.join(d_protein[pdb][ID].keys())
        print pdb, ID, hetIDs
stop

subdirs = os.listdir(pdbpath)
subdirs.sort()
for subdir in subdirs:
##    if subdir < 'y0':
##        continue
    print subdirs.index(subdir), len(subdirs), subdir
    files = os.listdir(pdbpath+subdir)
    for file in files:
        pdb = file[3:7]
        fd = open('%s%s/%s' %(pdbpath,subdir,file),'r')
        lines = fd.readlines()
        fd.close()
        d_IDs = {}
        parse = False
        for line in lines:

##            if line[:10] == 'REMARK 470':
##                if line.strip() == 'REMARK 470   M RES CSSEQI  ATOMS':
##                    parse = True
##                    continue
##                if parse == True:
##                    res_name = line[15:18].strip()
##                    chain = line[19]
##                    res_no = line[20:24]
##                    iCode = line[24]
##                    altloc = ' '
##                    ID = chain+res_no+iCode
##                    if ID not in d_IDs.keys():
##                        d_IDs[ID] = {'res_name':res_name,'altloc':altloc}
            
##            if line[:10] == 'REMARK 465':
##                if line.strip() == 'REMARK 465   M RES C SSSEQI':
##                    parse = True
##                    continue
##                if parse == True:
##                    res_name = line[15:18].strip()
##                    chain = line[19]
##                    res_no = line[22:26]
##                    iCode = line[26]
##                    altloc = ' '
##                    ID = chain+res_no+iCode
##                    if ID not in d_IDs.keys():
##                        d_IDs[ID] = {'res_name':res_name,'altloc':altloc}

            if line[:6].strip() in ['ATOM','HETATM']:
                res_name = line[17:20].strip()
                chain = line[21]
                res_no = line[22:26]
                iCode = line[26]
                altloc = line[16]
                ID = chain+res_no+iCode
                if ID not in d_IDs.keys():
                    d_IDs[ID] = {'res_name':res_name,'altloc':altloc}
##                    None
                elif res_name != d_IDs[ID]['res_name']:
##                    if line[:6].strip() == 'ATOM': ## remark 470 only!!!
##                        del d_IDs[ID] ## remark 470 only!!!
##                        continue ## remark 470 only!!!
                    if altloc != d_IDs[ID]['altloc']:
                        if pdb not in [
                            ## remark 465
                            '19gs','2are','1bls','1btw','2g73','1gar','1euw',
                            '1tv4','1n5j','2nw8',
                            ## remark 470
                            '',
                            ## same ID but different compounds (remediation)
                             ## all
                            '427d','1a48','1ajq','1bch','1blc',
                            ## multiplication of atoms necessary
                            '1k4q',
                             ## water
                            '1chw','1zeh','1lpb','1nhk',
                            ## non-chain hetero atoms
                            '2glz', ## ZN/NI
                            ## other remediation errors
                            '1epn','1h9h',
                            ## check again
                            '1aw8',
                            '1gai','1gah', ## MAN/SER/THR
                            ## no error
                            '354d','476d','2ci1','1din','1ysl','1g7b',
                            '2fi5','2ftm', ## ASP/IAS
                            ## stereoisomers
                            '2jge', ## SGR/SGX
                            '2jgj', ## SVV/SVW
                            ## hydroxyl group
                            '2aog','2bwx',
                            ## carboxyl group
                            '1blc','1rqh','1rr2',
                            '1k56','1k55', ## LYS/KCX
                            ## methyl group
                            '2aod','2aoc',
                            '1ob7', ## AIB/DIV
                            ## formyl group
                            '2ok6', ## TQQ/1TQ
                            ## acetyl group
                            '1xpk','1xpm','1xpl',
                            ## phosphate group
                            '2uv2',
                            ## ASP/THR
                            '2lhb',
                            ## SER/GLY
                            '1eis','1enm','1en2',
                            ## SER/PRO
                            '1cbn','1ejg','1jxt','1jxw','1jxx','1jxy','1jxu',
                            ## TRP/TYR
                            '1alx','1w5u',
                            ## TRP/TYR/PHE
                            '2izq',
                            ## ILE/VAL
                            '1alz','1al4','2hmz','2hmq',
                            ## MET/VAL
                            '1fh2','1eta',
                            ## TYR/HIS
                            '2ic4',
                            ## MSE/MET
                            '2q6u',
                            ]:
                            print pdb
                            print '"'+ID+'"'
                            print 'altloc', altloc
                            print 'res_name', res_name
                            print d_IDs[ID]
                            print line
                            stop
                    elif altloc == ' ':
                        if (res_name in l_nucleotides and d_IDs[ID]['res_name'] not in l_nucleotides+l_residues) or (d_IDs[ID]['res_name'] in l_nucleotides and res_name not in l_nucleotides+l_residues):
                            d_nucleotide[pdb] = {ID:{res_name:altloc,d_IDs[ID]['res_name']:d_IDs[ID]['altloc']}}
                        elif (res_name in l_residues and d_IDs[ID]['res_name'] not in l_nucleotides+l_residues) or (d_IDs[ID]['res_name'] in l_residues and res_name not in l_nucleotides+l_residues):
                            d_protein[pdb] = {ID:{res_name:altloc,d_IDs[ID]['res_name']:d_IDs[ID]['altloc']}}
                        elif res_name not in l_nucleotides+l_residues and d_IDs[ID]['res_name'] not in l_nucleotides+l_residues:
                            d_hetatm[pdb] = {ID:{res_name:altloc,d_IDs[ID]['res_name']:d_IDs[ID]['altloc']}}
                        else:
                            if pdb in [
                                ## coordiante section
                                '1h9h',
                                ## remark 465
                                '363d','1hef','1heg','2iid','2p83','1ofr','1og0','2nw8',
                                ## remark 470
                                '2aew','2eue',
                                ]:
                                continue
                            print pdb
                            print line
                            print res_name, '"'+ID+'"', altloc, d_IDs[ID]
                            notexpected
                    else:
                        if pdb in ['1gai','1gah',]:
                            continue
                        print pdb
                        print line
                        print res_name, '"'+ID+'"', altloc, d_IDs[ID]
                        notexpected


print d_nucleotide
print d_protein
print d_hetatm
