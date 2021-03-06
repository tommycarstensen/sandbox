## use normal modes to rank already identified pockets

import os

dn = 'POCASA'

l = os.listdir(dn)

l_chainIDs = [
    '2lipA', '1xqwA', '1bu7A', '1bzcA', '1esoA', '1tmlA', '1a8qA', '1jh6A',
    '1a4lA', '1i78A', '1smlA', '1qe3A', '1d3gA', '1akoA', '1rddA', '1p1xA',
    '1ak0A', '1ptdA', '2aceA', '1j00A', '1gojA', '1lbaA', '1bolA', '1vzxA',
    '1l9xA', '1og1A', '2cpoA', '1ru4A', '1ljlA', '1mrqA', '1w2nA', '1cz1A',
    '1mugA', '135lA', '1o8aA', '1rtuA', '1uchA', '1b6gA', '1chdA', '1pjaA',
    '1oxaA', '1r44A', '1w0hA', '1expA', '1bt1A', '1i9aA', '1im5A', '1pp4A',
    '2rnfA', '2ebnA', '1d1qA', '1ehyA', '1geqA', '1ca2A', '1gnsA', '1eh5A',
    '1l6pA', '1r4zA', '1a2tA', '1di1A', '2ayhA', '1astA', '1cm0A', '1cv2A',
    '1aj0A', '1dveA', '1fobA', '1ga8A', '1un1A', '1btlA', '2pecA', '1fcqA',
    '1czfA', '1thgA', '1booA', '1iu4A', '1bqcA', '206lA', '1cdeA', '1snzA',
    '1gq8A', '1aqlA', '1ps1A', '1s95A', '1pylA', '1ra2A', '1b6bA', '1pntA',
    '1e1aA', '2f9rA', '1v04A', '2nlrA', '1n29A', '1pbgA', '5cpaA', '1agmA', '1byaA', '1r76A', '1u5uA','1vidA', '1h4gA', '1akdA', '1fy2A', '1xqdA', '1d6oA', '1qv0A', '1qjeA', '1fvaA','1bp2A', '1ah7A', '2pthA', '2engA', '2acyA', '1qazA', '2a0nA', '1dl2A', '1gp5A','1onrA', '1cwyA', '1pudA', '1bs9A', '1dinA', '1xyzA', '1bwlA', '1eugA', '1idjA','1g24A', '1oygA', '1hzfA', '9papA', '1eb6A', '1ghsA', '1rbnA', '1bixA', '1bs4A','1celA', '1hkaA', '1b02A', '1qibA', '1u3fA', '1agyA', '1zioA', '1pa9A', '2tpsA','2plcA', '1qk2A', '1j53A', '1m21A',
    ]

##for fn in l:
for chainID in l_chainIDs:
    if os.path.isfile('%4s_%1s_probe.pdb' %(chainID[:4],chainID[-1],)):
        continue
    fn = chainID[:4]+'_Pocket_DepthCenters.pdb'
    print chainID
    os.system(
        'python /home/people/tc/svn/GoodVibes/goodvibes_ligand.py --pdb %s --chain %s --dist_max 6 --dist_min 3 --import_probe_coords %s/%s --multiple_modes' %(
            chainID[:4],chainID[-1],
            dn,fn,
            ))
##    print
##    stop

##for chainID in l_chainIDs:
