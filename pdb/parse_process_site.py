s = '''
1nzk
2o2n
2qm8
1sg2
1txu
1v58
1w5q
1xo5
1yto
1z2v
1z4s
1gaz
1ilz
1im0
1qa7
2zz0
2zzb
2zzc
1gam
1oc9
1r1z
1abs
1b8c
1b8r
2bqk
2i3v
2i3w
2mbw
2mgk
2mgl
2mgm
2p9k
1tye
1f1f
1wsz
1wt0
1wt1
1wt2
1wt3
1uxw
1e7n
1uoj
6yas
185d
2a97
2bhg
2bws
1qng
2c7d
2w63
1gnz
2jh2
1y1v
1p63
1usm
2c24
1gmv
1o6v
1uwg
2wb1
2jlp
1k4r
2vlc
2voy
2vri
2vwe
2w39
1yj9
1zbq
1e0z
1e3a
1c46
3ena
3etm
3evm
3evo
445d
1zzq
1e0p
1gwd
1ll0
1a37
1bs1
1vqw
1zrm
2hub
2zyp
3e5r
3e6a
3fyj
3guk
1h7r
2vze
3c7w
3c7y
3c7z
3c81
3c8q
3c8s
3cdq
3cdr
3cdv
1p2l
1pqm
1i01
3cdo
1e3h
1e3p
3f9l
3fad
3fi5
3htd
1l35
1lw9
2rb1
1swz
1xep
1xep
1qgc
2f1o
1fft
1fiw
2hoc
1i2a
1lsh
1ml7
1p6k
2psw
1sx2
2aig
3aig
1dyf
1gwd
3ht8
1l03
1l11
1l26
2q0m
1rcm
3e3d
3e3s
3e3t
2pgv
2pgy
1zcs
2hu3
3hwl
3l64
150l
1aki
5lyt
6lyt
1lsa
1lsb
1lsc
1lsd
1lse
1lsf
3lyt
4lyt
5lyt
6lyt
1pqj
3c80
1pqk
1t8g
2a0s
2a5f
2a9f
2ae8
2avy
2aw4
2aw7
2awb
1c5v
3dpr
2hkv
1q9i
1rvy
1rz8
1t98
1tg6
1tk3
1tkr
1txo
1u5v
1u63
1u6m
1u8e
1ut6
1va0
1vbn
1vf6
1vgp
1x6x
1x6y
1xdn
1xdp
1xfp
1xgf
1xv9
1xvp
1y19
1y63
1ylr
1ynl
2zah
1zix
1zjn
2g4p
2g4q
1lzb
1lzc
1lze
1lzg
3fw3
139l
1ckg
1fiy
1pkj
1pkk
1pkh
2vw4
2vw6
2vw7
3b95
3bpd
2bt2
1e2v
1e3u
1ha5
151l
1ppk
1ay6
1b5g
2z18
3g3x
3ht6
2rb0
1wtm
1wtn
1h87
1dpw
1dpx
1flq
1flu
1flw
1fly
2ihl
1lza
1lzb
1lzc
1lzd
1lze
1lzg
1v7s
1v7t
200l
226l
256l
260l
195l
196l
197l
198l
199l
1ks3
1kw5
1kw7
1ky0
1ky1
1l0j
1l0k
1l35
1l55
1l57
1l59
1l61
1l62
1l63
1lpy
2lzm
4lzm
5lzm
6lzm
1pqj
1xei
1xej
1xek
3f8v
3f9l
3fa0
3fad
1bql
1dkj
1dkk
1fn5
1acx
1azu
2hco
151l
1l35
155l
1c69
1c6a
1c6b
1hsw
1hsx
2bpu
2cgi
2d4i
3f8v
3f9l
3fa0
3fad
1l56
1l60
2oe4
1ovj
1p36
211d
1zcs
455d
1cwu
5dnb
3f29
1ggd
1inh
1kkr
1pyz
1qff
1ike
8lyz
1lks
2vb1
1xgr
2cgi
3d0y
1gk1
1thu
1h04
1ak9
1l9w
1y1a
1lks
1znm
1zy8
2beq
130d
2bvr
2by3
2wcl
2whr
2bvl
1d2u
2vbz
2wgm
1ztz
2wo8
1h0h
1o7t
2v8u
403d
1d3s
1dbw
1dck
1eqd
1eyq
1q2c
2uwu
2v1r
2gn5
1q9u
382d
2b2k
1f1h
1mgp
1n0r
4tln
3tmn
2hpl
1lu0
3pca
3pcj
3pck
3pcl
3pcm
7taa
302d
1a5k
1b0d
3hja
1mks
1mku
1tal
1wv7
3edz
3f0e
2fdn
2fe4
3gnq
1jzh
3e9i
1n7q
1njt
1nju
'''

l_pdbs = s.split()

import urllib2


print 'PDB\tSITE\tDEPOSITION\tREVISION\n'

for pdb in l_pdbs:

    print
    urllines = urllib2.urlopen('http://www.pdb.org/pdb/files/%s.cif?headerOnly=YES' %(pdb))
    lines = urllines.readlines()
    for i in range(len(lines)):
        line = lines[i]
        item = line.split()[0]
        if item == '_pdbx_database_status.process_site':
            print pdb, '\t',
            print line.split()[1], '\t',
            break
