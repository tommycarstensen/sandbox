import os
import exclude_redundancies
import sys
sys.path.append('/home/people/tc/svn/tc_sandbox/pdb/')
import parse_mmCIF

def main():

    ## read CSA database
    fd = open('/local/tc/goodvibes/CSA_2_2_12.dat','r')
    lines = fd.readlines()
    fd.close()

    l_chainIDs = []
    d_catres = {}
    for i_line in range(1,len(lines)):
##        if i_line > 1000:
##            break
        line = lines[i_line]
        if line[:4] != line[-5:].strip():
            continue

        l = line.strip().split(',')

        pdbID = l[0]
    ##    site_number = l[1]
        chainID = l[3]
        res_name = l[2]
        res_no = int(l[4])
        ref = l[7]
        ## ion
        if chainID == '':
            continue
        if res_name == '':
            continue

        if not pdbID+chainID in d_catres.keys():
            d_catres[pdbID+chainID] = {'res_names':[],'res_nos':[],'refs':[]}
        d_catres[pdbID+chainID]['res_names'] += [res_name]
        d_catres[pdbID+chainID]['res_nos'] += [res_no]
        d_catres[pdbID+chainID]['refs'] += [ref]

    l_chainIDs = d_catres.keys()

##    l_exclude = []
##    for chainID in l_chainIDs:
##        refs = d_catres[chainID]['refs']
##
##        if chainID[:4] in refs:
##            continue
####        if set([chainID[:4]]) == set(refs):
####            continue
##
##        res_names_ref = []
##        for ref in refs:
##            if chainID[:4] == ref:
##                continue
##            res_names_ref += d_catres[ref]['res_names']
####        if len( set( set(d_catres[chainID]['res_names']) ^ set(res_names_ref) ) & set(['ALA','CYS','ASP',]) ) > 0:
##        if len( set( set(d_catres[chainID]['res_names']) ^ set(res_names_ref) ) - set(['MG','CA',]) ) > 0:
##            for ref in set(refs):
##                print ref, d_catres[ref]['res_names']
##            print chainID, d_catres[chainID]['res_names']
##            print set( set(d_catres[chainID]['res_names']) ^ set(res_names_ref) ) - set(['MG','CA',])
##            print
##            l_exclude += [chainID]
##    stopstop

##    print len(l_chainIDs)
##    print
##    l_chainIDs = exclude(l_chainIDs)
##    print 
##    print l_chainIDs
##    print len(l_chainIDs)
##    fd = open('tmp.txt','w')
##    fd.write(str(l_chainIDs))
##    fd.close()
    l_chainIDs = ['2lipA', '1xqwA', '1bu7A', '1bzcA', '1esoA', '1tmlA', '1a8qA', '1jh6A', '1a4lA', '1i78A', '1smlA', '1qe3A', '1d3gA', '1akoA', '1rddA', '1p1xA', '1ak0A', '1ptdA', '2aceA', '1j00A', '1gojA', '1lbaA', '1bolA', '1vzxA', '1l9xA', '1og1A', '2cpoA', '1ru4A', '1ljlA', '1mrqA', '1w2nA', '1cz1A', '1mugA', '135lA', '1o8aA', '1rtuA', '1uchA', '1b6gA', '1chdA', '1pjaA', '1oxaA', '1r44A', '1w0hA', '1expA', '1bt1A', '1i9aA', '1im5A', '1pp4A', '2rnfA', '2ebnA', '1d1qA', '1ehyA', '1geqA', '1ca2A', '1gnsA', '1eh5A', '1l6pA', '1r4zA', '1a2tA', '1di1A', '2ayhA', '1astA', '1cm0A', '1cv2A', '1aj0A', '1dveA', '1fobA', '1ga8A', '1un1A', '1btlA', '2pecA', '1fcqA', '1czfA', '1thgA', '1booA', '1iu4A', '1bqcA', '206lA', '1cdeA', '1snzA', '1gq8A', '1aqlA', '1ps1A', '1s95A', '1pylA', '1ra2A', '1b6bA', '1pntA', '1e1aA', '2f9rA', '1v04A', '2nlrA', '1n29A', '1pbgA', '5cpaA', '1agmA', '1byaA', '1r76A', '1u5uA', '1vidA', '1h4gA', '1akdA', '1fy2A', '1xqdA', '1d6oA', '1qv0A', '1qjeA', '1fvaA', '1bp2A', '1ah7A', '2pthA', '2engA', '2acyA', '1qazA', '2a0nA', '1dl2A', '1gp5A', '1onrA', '1cwyA', '1pudA', '1bs9A', '1dinA', '1xyzA', '1bwlA', '1eugA', '1idjA', '1g24A', '1oygA', '1hzfA', '9papA', '1eb6A', '1ghsA', '1rbnA', '1bixA', '1bs4A', '1celA', '1hkaA', '1b02A', '1qibA', '1u3fA', '1agyA', '1zioA', '1pa9A', '2tpsA', '2plcA', '1qk2A', '1j53A', '1m21A']

    ##
    ## add probe atoms and calculate normal modes
    ##
    for chainID in l_chainIDs:

        if os.path.isfile('%s_%s_probe.pdb' %(chainID[:4],chainID[-1],)):
            continue

        fd = open('%s_%s_probe.pdb' %(chainID[:4],chainID[-1],),'w')
        fd.write(' ')
        fd.close()

        print chainID

        os.system('python /home/people/tc/svn/GoodVibes/goodvibes_ligand.py --pdb %s --chain %s --dist_max 6 --dist_min 3' %(chainID[:4],chainID[-1],))

    return


def exclude(l_chainIDs):

    ##
    ## exclude obsolete structures and theoretical structures
    ##
    print 'obsolete/theoretical'
    print len(l_chainIDs)
    l_exclude = []
    for chainID in l_chainIDs:
        if not os.path.isfile('/data/mmCIF/%s/%s.cif' %(chainID[1:3],chainID[0:4],)):
            l_exclude += [chainID]
    for chainID in l_exclude:
        l_chainIDs.remove(chainID)
    print len(l_chainIDs)
    print

    ##
    ## exclude multidomain structures
    ##
    print 'multidomain'
    print len(l_chainIDs)
    fd = open('../CathDomall','r')
    lines = fd.readlines()
    fd.close()
    l_single_domain_chains = []
    for line in lines:
        chainID = line[:5]
        if not chainID in l_chainIDs:
            continue
        n_domains = int(line[7:9])
        if n_domains == 1:
            l_single_domain_chains += [chainID]
    l_chainIDs = list( set(l_chainIDs) & set(l_single_domain_chains) )
    print len(l_chainIDs)
    print

    ##
    ## exclude multichain biological units
    ## exclude non-x-ray structures
    ##
    print 'multichain'
    print len(l_chainIDs)
    l_exclude = []
    l_pdbs_parsed = []
    d_resolutions = {}
    for i_chainID in range(len(l_chainIDs)):
        chainID = l_chainIDs[i_chainID]
        print i_chainID, len(l_chainIDs), chainID
        pdbID = chainID[:4]
        if pdbID in l_pdbs_parsed:
            continue
        d_mmCIF = parse_mmCIF.main(pdbID)

        l_pdbs_parsed += [pdbID]
          
        if d_mmCIF['_exptl.method'] != ['X-RAY DIFFRACTION']:
            l_exclude += [pdbID]
            continue

        try:
            l_oligomeric_counts = d_mmCIF['_pdbx_struct_assembly.oligomeric_count']
        except:
            print chainID
            continue
        if l_oligomeric_counts != len(l_oligomeric_counts)*['1']:
            l_exclude += [pdbID]

        try:
            d_resolutions[pdbID] = float(''.join(d_mmCIF['_refine_hist.d_res_high']))
        except:
            print chainID
            stop

    for chainID in list(l_chainIDs):
        if chainID[:4] in l_exclude:
            l_chainIDs.remove(chainID)
    print len(l_chainIDs)
    print

    ##
    ## exclude redundancies
    ##
    print 'redunant'
    print len(l_chainIDs)
    fd = open('../bc-50.out','r')
    lines = fd.readlines()
    fd.close()
    d = {}
    for i_line in range(len(lines)):
        line = lines[i_line]
        l_cluster = line.split()
        for i in range(len(l_cluster)):
            l_cluster[i] = l_cluster[i][:4].lower()+l_cluster[i][-1]
        l = list( set(l_cluster) & set(l_chainIDs) )
        if len(l) > 1:
            max_resolution = ['',None,]
            l.sort()
            for chainID in l:
                pdbID = chainID[:4]
                resolution = d_resolutions[pdbID]
                if resolution < max_resolution[0]:
                    max_resolution = [resolution,chainID,]
            for chainID in l:
                if chainID != max_resolution[1]:
                    l_chainIDs.remove(chainID)
    print len(l_chainIDs)
    print

    return l_chainIDs


if __name__ == '__main__':
    main()
