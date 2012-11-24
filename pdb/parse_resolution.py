## built-in
import os, sys
import parse_mmCIF

def main():

    d = {}

    if os.path.isfile('db_resolution.txt'):
        
        fd = open('db_resolution.txt','r')
        lines = fd.readlines()
        fd.close()

        for line in lines:
            l = line.strip().split()
            pdb = l[0]
            v = l[1]
            d[pdb] = v

    path = '/media/WDMyBook1TB/2TB/mmCIF'
    l_dns = os.listdir(path)
    l_dns.sort()

    lines_out = []

    for i in range(len(l_dns)):
        dn = l_dns[i]

        if dn < sys.argv[-1]:
            continue

        if not os.path.isdir('%s/%s' %(path,dn)):
            continue

        print '%s/%s %s' %(i+1,len(l_dns), dn)
        l_fns = os.listdir('%s/%s' %(path,dn))
        l_fns.sort()
        for fn in l_fns:
            if fn[-3:] == '.gz':
                continue

            pdb = fn[0:4]

            if pdb in d.keys():
                continue

            print pdb

            fd = open('%s/%s/%s' %(path, dn, fn), 'r')
            lines = fd.readlines()
            fd.close()

            d_mmCIF = parse_mmCIF.main(
                pdb,lines,
                l_data_categories = [
                    '_refine',
                    '_refine_hist',
                    ], ## parse selected data categories
                l_data_categories_break = [
                    '_refine',
##                    '_refine_hist',
                    ],
                d_breaks_negation = {
                    ## break if not x-ray diffraction
                    '_exptl.method':'X-RAY DIFFRACTION',
                    }
                )

            if d_mmCIF['_exptl.method'] != ['X-RAY DIFFRACTION']:
                continue

            resolution = d_mmCIF['_refine.ls_d_res_high']

            line = '%s %s\n' %(pdb,resolution,)
            lines_out += [line]

            fd = open('db_resolution.txt','a')
            fd.write(line)
            fd.close()

            d[pdb] = resolution

    ##
    ## write to file
    ##
    lines_out = []
    for pdb,resolution in d.items():
        line = '%s %s\n' %(pdb,resolution,)
        lines_out += [line]
    fd = open('db_resolution.txt','w')
    fd.writelines(lines_out)
    fd.close()

    d = {}
    fd = open('db_resolution.txt','r')
    lines = fd.readlines()
    fd.close()

    lines_out = []
    for line in lines:
        resolution = line.strip().split()[1][2:-2]
        if resolution == '.':
            continue
        resolution = float(resolution)
        resolution = round(resolution,2)
        if not resolution in d.keys():
            d[resolution] = 0
        d[resolution] += 1
        lines_out += ['%s\n' %(resolution,)]
    fd = open('histogram_resolution.txt','w')
    fd.writelines(lines_out)
    fd.close()
    stop

    lines_out = []
    l_resolutions = d.keys()
    l_resolutions.sort()
##    for resolution,count in d.items():
    for resolution in l_resolutions:
        count = d[resolution]
        lines_out += ['%s %s\n' %(resolution,count,)]
    fd = open('histogram_resolution.txt','w')
    fd.writelines(lines_out)
    fd.close()

    return

if __name__ == '__main__':
    main()
