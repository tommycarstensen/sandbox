## built in
import os, random
## add on
import numpy
import sys
sys.path.append('/home/tc/svn/tc_sandbox/misc')
import gnuplot, statistics, ols

bool_normalize = False

l_columns_html = [
    'gif1','gif2','pdb1', 'pdb2', 'bm1', 'bm2',
    'rmsd', 'mutations', 'chains', 'residues', 'coordinates',
    'pH1', 'pH2', 'T1', 'T2', 'res1', 'res2','spacegroup1', 'spacegroup2',
    'REMARK465', 'REMARK470','transformations',
    'title1','title2','hetIDs1', 'hetIDs2',
    ]

## HEWL
s_exclude = '2VB1, 3LZT, 1IEE, 4LZT, 3AJN, 1V7S, 1V7T, 2D4J, 2D4K, 2D4I, 1WTN, 2Z18, 2Z12, 3AGI, 2HUB, 2Z19, 1QIO, 2HU3, 2ZYP, 3D9A, 2XBR, 2D6B, 2BLY, 2BLX, 2XBS, 1WTM, 193L, 3P64, 2CGI, 3N9E, 2BPU, 3N9A, 2IHL, 194L, 3A8Z, 3P66, 2HTX, 2PC2, 2X0A, 3EXD, 3N9C, 2FBB, 3A90, 1LZB, 3AGH, 1RJC, 1SQ2, 2F2N, 1LKR, 3A92, 1VDQ, 1VAU, 2W1L, 2C8P, 1HF4, 3E3D, 2C8O, 1XFP, 3A95, 1LZA, 1A2Y, 2HU1, 1ZVH, 3A91, 1AKI, 3AGG, 3A93, 2F30, 1IR8, 3KAM, 1YIL, 3A94, 1GPQ, 1T3P, 1VDS, 1H6M, 1HEL, 1LZ8, 2ZQ3, 1DPW, 3P4Z, 3QY4, 2XJW, 1LSM, 2A7D, 1B2K, 2F4G, 1VAT, 1VDT, 1DPX, 1HEM, 1LCN, 1VDP, 2B5Z, 1W6Z, 3A34, 1LZC, 3P68, 1KXX, 1ZVY, 1HEN, 1UIH, 1RFP, 2W1X, 2W1Y, 1LZE, 1HEP, 1BGI, 1IOS, 1UIB, 1H87, 1IOT, 3IJU, 3EMS, 2WAR, 1HEQ, 1N4F, 1LSE, 1LZG, 1LMA, 1HEO, 1IOR, 1GWD, 1HER, 1LYS, 1YQV, 1LSC, 1YIK, 1T6V, 3OK0, 1FLQ, 1IOQ, 1FN5, 1HC0, 1LSB, 1LSD, 1LZD, 1UIA, 2H9K, 3RT5, 1FLW, 2YVB, 1FLU, 1F10, 2H9J, 1BWI, 2LZT, 3OJP, 1KXY, 1BWJ, 132L, 2AUB, 1LSA, 2W1M, 2Q0M, 1LZ9, 2G4P, 1LSN, 3A96, 2G4Q, 1LSF, 1FLY, 1G7J, 1BVX, 1AT5, 2A7F, 1BWH, 1AT6, 1VFB, 3A6C, 2DQD, 1DKK, 2I25, 3SP3, 1G7I, 3A67, 5LYM, 1RI8, 1UIF, 1J1X, 2DQC, 1J1O, 1J1P, 1YKZ, 1AZF, 3RNX, 3IJV, 3A6B, 2DQJ, 1KIQ, 1LJF, 1LJH, 2XTH, 1IR9, 6LYT, 2F4A, 1VED, 1IR7, 2LYM, 3LYM, 1UC0, 1JJ3, 5LYT, 4LYT, 3RZ4, 3RW8, 1UIG, 1YKX, 1UID, 1UIC, 1UIE, 1UUZ, 1HEW, 1KXW, 1B0D, 1YKY, 1RCM, 1JJ1, 2LYO, 1NBY, 3LYO, 1F0W, 1G7H, 2I6Z, 1JIY, 1JIS, 1HSX, 1JIT, 1LYO, 1KIR, 2DQE, 1JJ0, 1G7M, 2ZQ4, 1UCO, 1LPI, 3M18, 1LJ4, 1LJE, 3LYT, 1DKJ, 1QTK, 1C10, 1LJ3, 1LJG, 1NDG, 2CDS, 3A3R, 1YL1, 3AW6, 1LJI, 1LJJ, 1KIP, 4LYM, 3AW7, 2EIZ, 1NBZ, 1PS5, 4LYO, 1Z55, 2YBI, 2YDG, 1LJK, 1HSW, 1YL0, 2YBH, 1UA6, 1IC7, 2YBJ, 1JTT, 1P2C, 2YBL, 1G7L, 1XEI, 1JPO, 3A3Q, 2YBM, 1ZV5, 1DQJ, 1XEJ, 2DQI, 2YBN, 2EKS, 1IC5, 1NDM, 3P65, 1XGP, 3T6U, 3M3U, 1LZT, 1XGU, 1XGQ, 2ZNX, 1XEK, 2EPE, 3G3A, 2DQG, 2DQH, 1XGR, 2D91, 1XGT, 1IC4, 1MLC, 1FDL, 1C08, 2YSS, 2I26, 1BQL, 1JTO, 1MEL, 2DQF, 2ZNW, 3G3B, 9LYZ, 1BVK, 3F6Z, 1ZMY, 3HFM, 1BHZ, 2HSO, 2HS7, 2HS9, 2A6U, 1SF4, 1SF6, 1SF7, 1SFB, 1SFG, 1GXX, 1GXV, 1JA2, 1JA4, 1JA6, 1JA7, 1IO5, 1E8L, 1LZN, 1LKS, 2IFF, 1LZH, 2LZH, 8LYZ, 7LYZ, 1LYZ, 2LYZ, 3LYZ, 4LYZ, 5LYZ, 6LYZ'
## T4L
s_exclude += '1SX7,1SX2,1SWY,1SWZ,3F8V,3FA0,3F9L,3FAD,3DKE,3HH6,3HH5,3HH3,3HH4,3HTG,3HUK,2RBO,2RBN,2Q9D,2IGC,2NTG,3SBB,3HTD,3HUA,3HUQ,2RBR,3HU9,3GUI,2RB2,2RBP,1P36,1LW9,1PQM,1P7S,3K2R,1P2L,1P6Y,2RBS,3FI5,1XEP,1P3N,1P37,3GUP,2OU9,3HT6,119L,3GUN,1P2R,161L,1L58,1L92,217L,1L83,1L62,1L19,3HT8,1ZUR,1L47,1PQI,244L,1L46,1L17,1L18,1L23,1L33,3LZM,131L,221L,1L66,1L60,1L24,1L52,1P64,241L,237L,1L74,219L,1G0M,128L,110L,1L65,1L30,1L22,1LYH,139L,130L,1L45,1L36,1L03,3C8Q,1L63,129L,1L68,1L44,4LZM,1PQD,1L48,1L32,3DMV,1G07,3GUJ,3CDT,165L,2RAZ,166L,1P46,1L90,211L,138L,1L10,1L59,1L41,1L06,1L13,1L15,1KNI,1L12,1L14,212L,3C7Z,247L,240L,243L,3C82,1L07,1L16,1L94,1KW5,1L29,163L,118L,1D9W,1L09,159L,1L93,1L91,2A4T,1L11,156L,122L,108L,1G0Q,164L,123L,2O79,173L,1L02,1L04,1L05,232L,1L26,1L01,158L,160L,114L,2RBQ,3DNA,162L,1G0L,1DYE,115L,206L,1L08,1ZYT,102L,238L,126L,111L,1LYE,1LYG,7LZM,3HT7,188L,1L38,3CDQ,1PQO,1CVK,1L42,1L43,1L35,1L27,1DYB,181L,1L61,127L,1LPY,1LYJ,5LZM,125L,1QSB,255L,239L,245L,246L,250L,229L,6LZM,157L,2RB1,2LZM,1G0J,1G0P,112L,1LWG,242L,184L,120L,113L,1L98,180L,256L,1LYF,1CV3,186L,220L,3C8S,1C69,1QT7,1NHB,107L,1L86,1L87,1L80,137L,3DN8,1L49,1L37,1L39,1L40,3CDR,1QT5,1QUD,182L,1L88,1L20,185L,187L,1G0K,1CU2,1QT3,1G06,1QS9,3C8R,1G1W,258L,1L56,1L31,226L,109L,1CV4,2CUU,183L,3CDV,3C7W,1CV5,199L,155L,2HUL,3HTF,224L,249L,1L51,1L25,1L50,146L,228L,3DN2,3HU8,1QSQ,1L67,1L28,3DN3,1L72,1L73,3DMX,1G0G,1L54,1D2W,1T8G,1L71,1KW7,235L,3DN6,3G3X,1ZWN,1L00,1CUP,3L2X,3DN4,254L,222L,2RAY,1L21,1CV6,1DYA,3HWL,248L,236L,225L,2OU8,260L,195L,1LGW,1QTB,1L84,1L69,3DN1,1L53,1CX7,1L55,1L57,1G1V,1C6Q,3DN0,1C6H,1C6I,1QT6,2NTH,3C83,1C6L,1QUH,1L79,1L76,223L,148L,1L75,1EPY,1QT8,233L,192L,1DYF,1C6J,210L,1LLH,1LI3,234L,1L89,1C6E,1C6K,1QUO,3HTB,1L70,1I6S,257L,1C6G,3C81,1OWY,1LGX,1L34,2O4W,1L64,2OTY,227L,3CDO,1KY0,1C6P,141L,1LYI,1L0J,2RB0,1L99,103L,214L,1D3J,1PQJ,1LGU,230L,198L,143L,191L,1P56,1B6I,2L78,1CX6,1OWZ,1C6F,1L0K,142L,1T6H,3L64,1QUG,1QTH,1L85,200L,2OEA,172L,1C6C,253L,1L95,1L81,1C60,1C6T,1QTZ,259L,2OE9,145L,147L,1CUQ,1D3N,1TLA,1L96,1L97,1LI2,1C64,152L,1C61,1LI6,1C65,1C6D,190L,1D3F,1CU5,1L77,1C63,1OVH,144L,3G3V,140L,1DYC,3HT9,1KY1,2OU0,1D2Y,1CV1,3DMZ,1PQK,2HUK,1DYD,201L,1CV0,2B6T,2OE7,1CU3,1C6M,1D3M,3GUK,215L,1QT4,252L,197L,2OE4,2B6X,1CU6,1OV7,1C6A,2B75,1LYD,1JTM,205L,1DYG,1KS3,2B72,1CTW,3C80,1OVJ,1L82,3C7Y,2B74,175L,2B73,218L,2B6W,1OVK,1OV5,1C6B,213L,1SSW,216L,1CU0,1C66,196L,167L,2Q9E,176L,1QS5,2B6Z,1C6N,151L,3G3W,1C62,2B70,2OTZ,1C67,1LWK,2B6Y,1QTC,1QTV,3GUO,150L,3GUL,1SSY,1T8F,174L,251L,1C68,1QTD,3GUM,231L,3SB5,3SB9,1JTN,189L,2QAR,149L,209L,171L,2HUM,104L,3SB7,178L,169L,177L,170L,3SB8,1P5C,168L,1JQU,3SB6,3SBA,2B7X,3JR6,2LC9,2LCB'

def analyze_rmsd():

    print 'analyze rmsd'

    ## discrete (not continuous) data
    d_data_ratio_discrete = {
        'mutations':{},
        'chains':{},
##        'residues':{}, ## use average radius of gyration instead
##        'coordinates':{}, ## use average radius of gyration instead
        }
    ## continuous data
    d_data_ratio_continuous = {
        'pH':{'single':[],'min':[],'max':[],'average':[],'difference':[],}, ## discrete?
        'T':{'single':[],'min':[],'max':[],'average':[],'difference':[],},
        'res':{'single':[],'min':[],'max':[],'average':[],'difference':[],}, ## discrete?
        'radiusGyration':{'single':[],'min':[],'max':[],'average':[],'difference':[],},
        'MV':{'single':[],'min':[],'max':[],'average':[],'difference':[],},
        }
    ## nominal (categorical) data - binary in this case
    d_data_nominal = {
##        'spacegroup':{'identical':{},'different':{}},
        'spacegroup_categorical':{'identical':{},'different':{}},
        'spacegroup_binary':{'identical':[],'different':[]},
        'hetIDs':{'identical':[],'different':[]},
        'REMARK465':{'True':[],'False':[]},
        'REMARK470':{'True':[],'False':[]},
        'transformations':{'True':[],'False':[]},
        'authors':{'identical':[],'different':[]},
##        'exptl_crystal_grow':{'identical':[],'different':[]},
        'pH':{'identical':[],'different':[]},
        'T':{'identical':[],'different':[]},
        'cryo':{'True':[],'False':[]},
        }
    ##
    d_data_nominal_ANOVA = {}
    for k1 in d_data_nominal.keys():
        for k2 in d_data_nominal.keys():
            if k1 == k2:
                continue
            if 'transformations' in [k1,k2,]:
                continue
            k = k1+'+'+k2
            d_data_nominal_ANOVA[k] = {}
            for k1sub in d_data_nominal[k1].keys():
                for k2sub in d_data_nominal[k2].keys():
                    ksub = k1sub+'+'+k2sub
                    d_data_nominal_ANOVA[k][ksub] = []
    ## combinations... subsets of data
    d_combined = {
        'chains_no_transformations':{}, ## no transformations not default!!! this checks effect of transformations!!!
##        'residues_no_transformations':{},
##        'mutations_no_transformations':{},
##        'residues_1_chain':{},
##        'mutations_1_chain':{}, ## also same space group - same SG and 1 chain already default
##        'pHdiff_nomutation':{}, ## 0 mutations, 1 chain and identical space groups by default
        }
    ## spacegroup vs spacegroup (average RMSD)
    d_spacegroups_data = {}
    ## cryo/noncryo vs cryo/noncryo (average RMSD)
    d_temperature_data = {'cryo':[],'non-cryo':[],'mixed':[],}

    ## ordinary least squares regression
    d_ols = {
        'rmsd':[],'Tdiff':[],'pHdiff':[],'resmax':[],
        'radiusGyration':[],'MV':[],
        }

    n_lines = 1+len(l_columns_html)+1

    htmls = os.listdir('htm/')
    htmls.sort()
    lines = []
    d_pdbs_skip = {}
##    l_gnuplot = []

    ##
    ## parse radius of gyration
    ##
    fd = open('/home/tc/svn/tc_sandbox/pdb/radius_of_gyration.txt','r')
    lines_rg = fd.readlines()
    fd.close()
    d_rg = {}
    for line in lines_rg:
        l = line.split()
        pdb = l[0]
        rg = l[1]
        d_rg[pdb] = rg

    ##
    ## parse Matthews Coefficient
    ##
    fd = open('/home/tc/svn/tc_sandbox/pdb/db_MatthewsCoefficient.txt','r')
    lines = fd.readlines()
    fd.close()
    d_MV = {}
    for line in lines:
        l = line.split()
        pdb = l[0]
        MV = l[1]
        d_MV[pdb] = MV

    ##
    ## parse authors
    ##
    fd = open('/home/tc/svn/tc_sandbox/pdb/db_authors.txt','r')
    lines_authors = fd.readlines()
    fd.close()
    d_authors = {}
    for line in lines_authors:
        l = line.split()
        pdb = l[0]
        s_authors = ''.join(l[1:])
        l_authors = s_authors.split(';')
        d_authors[pdb] = l_authors

    ##
    ## parse growth conditions
    ##
##    fd = open('/home/tc/svn/tc_sandbox/pdb/db_exptl_crystal_grow.txt','r')
    fd = open('/home/tc/svn/tc_sandbox/pdb/db_exptl_crystal_grow_comp.txt','r')
    lines = fd.readlines()
    fd.close()
    d_grow = {}
    for line in lines:
        pdb = line[:4]
##        s = line[5:]
        set_compounds = set(eval(line[5:].upper()))
##        d_grow[pdb] = s
        d_grow[pdb] = set_compounds

    ##
    ## parse htmls
    ##
    for i in range(len(htmls)):
        html = htmls[i]
        pdb = html[:4]
        if i % 200 == 0:
            print '%4.1f%%' %(100*float(i)/len(htmls))
##        if i > 1000: ## tmp!!!!!!
##            break
        fd = open('htm/%s' %(html),'r')
        lines = fd.readlines()
        fd.close()
        if lines[0][:8] == '<script ':
            lines = lines[2+n_lines:-1]
        else:
            lines = lines[1+n_lines:-1]

        ##
        ## parse data
        ##
        (
            d_data_ratio_discrete,d_data_ratio_continuous,d_data_nominal,
            d_pdbs_skip,
            d_combined,
##            l_gnuplot,
            d_spacegroups_data,
            d_ols,
            d_temperature_data,
            d_data_nominal_ANOVA,
            ) = parse_html(
                lines,n_lines,d_pdbs_skip,
                d_data_ratio_discrete,d_data_ratio_continuous,d_data_nominal,
                d_combined,
                html,
##                l_gnuplot,
                d_rg, d_MV, d_grow, ## pre-parsed values
                d_spacegroups_data,
                d_ols,
                d_temperature_data,
                d_authors, ## pre-parsed values
                d_data_nominal_ANOVA,
                )

##    fd = open('data.gnu','w')
##    fd.writelines(l_gnuplot)
##    fd.close()

    ordinary_linear_regression(d_ols)

    for k in ['non-cryo','cryo','mixed',]:
        print k, len(d_temperature_data[k]), sum(d_temperature_data[k])/len(d_temperature_data[k])

    ##
    ## plot continuous (non-discrete) data
    ##
    plot_continuous(d_data_ratio_continuous)

    ##
    ## plot combined discrete data
    ##
    plot_subsets(d_combined)

    ##
    ## plot discrete (non-continuous) data
    ##
    plot_discrete(d_data_ratio_discrete)

    ##
    ## plot spacegroups against each other
    ##
    plot_spacegroup_vs_spacegroup_vs_rmsd(d_spacegroups_data,d_data_nominal['spacegroup_categorical'],)

    ##
    ## nominal/categorical data
    ## binary: hetID difference, remarks present, transformations, space group difference (True/False)
    ##
    plot_nominal(d_data_nominal,d_data_nominal_ANOVA,)

    print
    print 'n', 'pH', 'diff', len(d_data_ratio_continuous['pH']['difference'])
    print 'n', 'T', 'diff', len(d_data_ratio_continuous['T']['difference'])
    print 'n', 'res', 'max', len(d_data_ratio_continuous['res']['max'])
    print 'n', 'radiusGyration', 'average', len(d_data_ratio_continuous['radiusGyration']['average'])
    print 'n', 'MV', 'average', len(d_data_ratio_continuous['MV']['average'])

    return


def ordinary_linear_regression(d_ols):

    print
    print '**** ordinary least squares linear regression'

    ## ordinary least squares regression
    ## http://www.velocityreviews.com/forums/t357533-scipy-i-need-an-example-of-use-of-linalg-lstsq.html
    y = numpy.array(d_ols['rmsd'])

    x1 = numpy.array(d_ols['pHdiff'])
    x2 = numpy.array(d_ols['radiusGyration'])
    x3 = numpy.array(d_ols['resmax'])
    x4 = numpy.array(d_ols['Tdiff'])
    x5 = numpy.array(d_ols['MV'])

    A = numpy.column_stack([
        x1, x2,
        x3,
        x4,
##        x5,
        numpy.ones(len(y), float)
        ])
    parameters, residuals, rank, singularvalues = numpy.linalg.lstsq(A,y)
    l_parameters = [
        'pHdiff', 'radiusGyration',
        'resmax',
        'Tdiff',
##        'MV',
        ]
    for i_parameter in range(len(l_parameters)):
        print l_parameters[i_parameter], parameters[i_parameter]
    print 'residuals', residuals ## sum of the residuals: sum((p[0] * A[:,0] + p[1] * A[:,1] - y)**2)
    print 'rank', rank ## number of column vectors
    print singularvalues

    ##
    ## OrdLS
    ##
    print
    print 'n', len(x1), len(x2), len(x3), len(x4), len(x5)
    x = numpy.column_stack([
        x1, x2,
        x3,
        x4,
##        x5,
        ])
    mymodel = ols.ols(y,x,'rmsd',l_parameters)
    print 'probabilities', mymodel.p
    print mymodel.summary()

    lines = []
    for i in range(len(y)):
        lines += [
            '%s %s\n' %(
                y[i],
                parameters[0]*x1[i]+parameters[1]*x2[i]+parameters[2]*x3[i]+parameters[3]*x4[i]+parameters[4],
                )
            ]
    fd = open('tmp.txt','w')
    fd.writelines(lines)
    fd.close()

    return


def plot_spacegroup_vs_spacegroup_vs_rmsd(d_spacegroups_data,d_spacegroups_data_1d,):

    ##
    ## exclude low frequency spacegroups
    ##
    l_spacegroups = []
    for spacegroup in d_spacegroups_data.keys():

        ## dont plot rarely occuring space groups
        if not spacegroup in d_spacegroups_data_1d['identical'].keys():
            continue
        if not spacegroup in d_spacegroups_data_1d['different'].keys():
            continue
        ## 200 incl ...
        ## 150 incl P 63, P 43 21 2
        ## 100 incl ...
        ## 50 incl ...
        if len(d_spacegroups_data_1d['identical'][spacegroup]) < 50:
            continue
        if len(d_spacegroups_data_1d['different'][spacegroup]) < 50:
            continue

        l_spacegroups += [spacegroup]

    l_spacegroups.sort()

    ##
    ## loop over high frequency space groups
    ##
    d_xtics_contour = {}
    l_gnuplotdata_contour = []
    for i_spacegroup1 in range(len(l_spacegroups)):
        spacegroup1 = l_spacegroups[i_spacegroup1]

        ## set xtic
        d_xtics_contour['%-10s' %(spacegroup1)] = float(i_spacegroup1+.5)

        for i_spacegroup2 in range(len(l_spacegroups)):
            spacegroup2 = l_spacegroups[i_spacegroup2]
            if spacegroup1 not in d_spacegroups_data.keys():
                rmsd_average = 0.0
            elif spacegroup2 not in d_spacegroups_data[spacegroup1].keys():
                rmsd_average = 0.0
            elif spacegroup1 == spacegroup2:
##                rmsd_average = 0
                l_rmsds = d_spacegroups_data[spacegroup1][spacegroup2]
                rmsd_average = sum(l_rmsds)/len(l_rmsds)
            else:
                l_rmsds = d_spacegroups_data[spacegroup1][spacegroup2]
                rmsd_average = sum(l_rmsds)/len(l_rmsds)
            if len(l_rmsds) < 10 or rmsd_average == 0.0:
                rmsd_average = 0.25

            l_gnuplotdata_contour += [
                '%s %s %s\n' %(i_spacegroup1, i_spacegroup2, rmsd_average,),
                ]

        l_gnuplotdata_contour += [
            '%s %s %s\n' %(i_spacegroup1, i_spacegroup2+1, rmsd_average,),
            ]
        l_gnuplotdata_contour += ['\n']

    for i_spacegroup2 in range(len(l_spacegroups)):
        l_gnuplotdata_contour += [
            '%s %s %s\n' %(i_spacegroup1+1, i_spacegroup2, rmsd_average,),
            ]
    l_gnuplotdata_contour += [
        '%s %s %s\n' %(i_spacegroup1+1, i_spacegroup2+1, rmsd_average,),
        ]
    l_gnuplotdata_contour += ['\n']

    prefix_gnuplot = 'ps/spacegroups'
    gnuplot.contour_plot(
        prefix_gnuplot+'_contour', l_gnuplotdata_contour,
        title = 'RMSD as a function of space group differences',
        d_xtics = d_xtics_contour,
        d_ytics = d_xtics_contour,
        zlabel = 'average RMSD',
        x2 = len(l_spacegroups),
        y2 = len(l_spacegroups),
##        z2=2.0,
        )

    return


def write_html(d_rmsd, d_header, prefix):

    ## keys=htmlkeys, values=htmlcolumns (table headings)
    d_columns_headers = {
        'REMARK465':'missing<br>residues',
        'REMARK470':'missing<br>atoms',
        'transformations':'transformations',
        'identical authors':'identical<br>authors',
        'gif1':'gif1','gif2':'gif2',
        'pdb1':'pdb1', 'pdb2':'pdb2',
        'bm1':'bm1', 'bm2':'bm2',
        'rmsd':'<a href="http://en.wikipedia.org/wiki/Protein_structural_alignment">rmsd</a>',
        'chains':'chains (biounit)', 'residues':'residues (biounit)', 'coordinates':'coords (biounit)',
        'pH1':'pH1', 'pH2':'pH2', 'T1':'T1', 'T2':'T2',
        'res1':'<a href="http://alpha.rcsb.org/pdb/statistics/histogram.do?mdcat=refine&mditem=ls_d_res_high&minLabel=0&maxLabel=5&numOfbars=10">resolution1</a>',
        'res2':'<a href="http://alpha.rcsb.org/pdb/statistics/histogram.do?mdcat=refine&mditem=ls_d_res_high&minLabel=0&maxLabel=5&numOfbars=10">resolution2</a>',
        'spacegroup1':'<a href="http://alpha.rcsb.org/pdb/statistics/histogram.do?mdcat=symmetry&mditem=space_group_name_H-M&numOfbars=200">space1</a>',
        'spacegroup2':'<a href="http://alpha.rcsb.org/pdb/statistics/histogram.do?mdcat=symmetry&mditem=space_group_name_H-M&numOfbars=200">space2</a>',
        'title1':'title1','title2':'title2',
        'hetIDs1':'hetIDs1', 'hetIDs2':'hetIDs2',
        'mutations':'mutations (biounit)',
        }

    ## initiate html lines
    th = '<tr>\n'
    for column in l_columns_html:
        th += '<td>%s</td>\n' %(d_columns_headers[column])
    th += '</tr>\n'

    l_tr = ''

    d_html = {}

    ##
    ## loop over pdbs and biomolecules
    ##
    for pdb1 in d_rmsd:

        ## parse physiochemical properties (pH and temperature)
        if 'REMARK200' in d_header[pdb1].keys():
            T1 = d_header[pdb1]['REMARK200']['TEMPERATURE']
            pH1 = d_header[pdb1]['REMARK200']['PH']
        else:
            T1 = 'N/A'
            pH1 = 'N/A'
        try:
            T1 = '%5.1f' %(float(T1))
        except:
            T1 = '%s' %(T1.rjust(5))
        try:
            pH1 = '%4.1f' %(float(pH1))
        except:
            pH1 = '%s' %(pH1.rjust(4))
        ## parse x-ray space group
        spacegroup1 = d_header[pdb1]['CRYST1'].rjust(10)
        ## parse hetIDs
        hetIDs1 = set()
        ## parse resolution
        res1 = d_header[pdb1]['REMARK2']
        try:
            res1 = '%5.2f' %(float(res1))
        except:
            res1 = '%s' %(res1.rjust(5))
        ## parse hetIDs
        for chain in d_header[pdb1]['HET'].keys():
            for res_no in d_header[pdb1]['HET'][chain].keys():
                for iCode in d_header[pdb1]['HET'][chain][res_no].keys():
                    hetIDs1 |= d_header[pdb1]['HET'][chain][res_no][iCode]
        hetIDs1 = list(hetIDs1)

        for bm1 in d_rmsd[pdb1].keys():

            for pdb2 in d_rmsd[pdb1][bm1]:

                ## parse physiochemical properties (pH and temperature)
                if 'REMARK200' in d_header[pdb2].keys():
                    T2 = d_header[pdb2]['REMARK200']['TEMPERATURE']
                    pH2 = d_header[pdb2]['REMARK200']['PH']
                else:
                    T2 = 'N/A'
                    pH2 = 'N/A'
                try:
                    T2 = '%5.1f' %(float(T2))
                except:
                    T2 = '%s' %(T2.rjust(5))
                try:
                    pH2 = '%4.1f' %(float(pH2))
                except:
                    pH2 = '%s' %(pH2.rjust(4))
                ## parse x-ray space group
                spacegroup2 = d_header[pdb2]['CRYST1'].rjust(10)
                ## parse hetIDs
                hetIDs2 = set()
                ## parse resolution
                res2 = d_header[pdb2]['REMARK2']
                try:
                    res2 = '%5.2f' %(float(res2))
                except:
                    res2 = '%s' %(res2.rjust(5))
                ## parse hetIDs
                for chain in d_header[pdb2]['HET'].keys():
                    for res_no in d_header[pdb2]['HET'][chain].keys():
                        for iCode in d_header[pdb2]['HET'][chain][res_no].keys():
                            hetIDs2 |= d_header[pdb2]['HET'][chain][res_no][iCode]
                hetIDs2 = list(hetIDs2)

                for bm2 in d_rmsd[pdb1][bm1][pdb2].keys():

                    ## parse rmsd and related data
                    rmsd = '%5.2f' %(d_rmsd[pdb1][bm1][pdb2][bm2]['rmsd'])
                    n_chains = '%3i' %(d_rmsd[pdb1][bm1][pdb2][bm2]['chains'])
                    n_residues = '%4i' %(d_rmsd[pdb1][bm1][pdb2][bm2]['residues'])
                    n_coordinates = '%5i' %(d_rmsd[pdb1][bm1][pdb2][bm2]['coordinates'])
                    mutations = '%2i' %(d_rmsd[pdb1][bm1][pdb2][bm2]['mutations'])
                    if 'REMARK465' in d_header[pdb1].keys() or 'REMARK465' in d_header[pdb2].keys():
                        remark465 = True
                    else:
                        remark465 = False
                    if 'REMARK470' in d_header[pdb1].keys() or 'REMARK470' in d_header[pdb2].keys():
                        remark470 = True
                    else:
                        remark470 = False
                    transformations = d_rmsd[pdb1][bm1][pdb2][bm2]['transformations']
                    bool_identical_authors = d_rmsd[pdb1][bm1][pdb2][bm2]['bool_identical_authors']

                    ## write data to dictionary
                    htmlhetIDs1 = ''
                    for i in range(len(hetIDs1)):
                        hetID = hetIDs1[i]
                        htmlhetIDs1 += '<a href="http://ligand-depot.rcsb.org/pub/%s/%s.html">%s</a>, ' %(hetID,hetID,hetID)
                        if (i+1) % 3 == 0:
                            htmlhetIDs1 += '<br>'
                    htmlhetIDs1 = htmlhetIDs1[:-2]
                    htmlhetIDs2 = ''
                    for i in range(len(hetIDs2)):
                        hetID = hetIDs2[i]
                        htmlhetIDs2 += '<a href="http://ligand-depot.rcsb.org/pub/%s/%s.html">%s</a>, ' %(hetID,hetID,hetID)
                        if (i+1) % 3 == 0:
                            htmlhetIDs2 += '<br>'
                    htmlhetIDs2 = htmlhetIDs2[:-2]

                    prefix1 = '%s%s%s%s' %(pdb1, str(bm1).zfill(2), pdb2, str(bm2).zfill(2))
                    prefix2 = '%s%s%s%s' %(pdb2, str(bm2).zfill(2), pdb1, str(bm1).zfill(2))
                    d_columns_data = {
                        'transformations':'<td>%s' %(transformations),
                        'identical authors':'<td>%s' %(bool_identical_authors),
                        'REMARK465':'<td>%s' %(remark465),
                        'REMARK470':'<td>%s' %(remark470),
                        'title1':'<td style="font-size:70%%">%s' %(d_header[pdb1]['TITLE']),
                        'title2':'<td style="font-size:70%%">%s' %(d_header[pdb2]['TITLE']),
                        'gif1':'<td><a href="../pdb/%s/%s.pdb"><img src="../gif/%s/%s.gif"></a>' %(pdb1[1],prefix1,pdb1[1],prefix1),
                        'gif2':'<td><a href="../pdb/%s/%s.pdb"><img src="../gif/%s/%s.gif"></a>' %(pdb2[1],prefix2,pdb2[1],prefix2),
                        'pdb1':'<td><a href="../htm/%s.htm">%s</a>' %(pdb1,pdb1),
                        'pdb2':'<td><a href="../htm/%s.htm">%s</a>' %(pdb2,pdb2),
                        'bm1':'<td style="text-align: right">%s' %(bm1),
                        'bm2':'<td style="text-align: right">%s' %(bm2),
                        'pH1':'<td style="text-align: right">%s' %(pH1),
                        'pH2':'<td style="text-align: right">%s' %(pH2),
                        'T1':'<td style="text-align: right">%s' %(T1),
                        'T2':'<td style="text-align: right">%s' %(T2),
                        'res1':'<td style="text-align: right">%s' %(res1),
                        'res2':'<td style="text-align: right">%s' %(res2),
                        'spacegroup1':'<td nowrap>%s' %(spacegroup1),
                        'spacegroup2':'<td nowrap>%s' %(spacegroup2),
                        'hetIDs1':'<td style="font-size:80%%" nowrap>%s' %(htmlhetIDs1),
                        'hetIDs2':'<td style="font-size:80%%" nowrap>%s' %(htmlhetIDs2),
                        'rmsd':'<td style="text-align: right">%s' %(rmsd),
                        'chains':'<td style="text-align: right">%s' %(n_chains),
                        'residues':'<td style="text-align: right">%s' %(n_residues),
                        'coordinates':'<td style="text-align: right">%s' %(n_coordinates),
                        'mutations':'<td style="text-align: right">%s' %(mutations)
                        }

                    d_columns_data_txt = {
                        'transformations':transformations,
                        'bool_identical_authors':bool_identical_authors,
                        'REMARK465':remark465,
                        'REMARK470':remark470,
                        'pH1':pH1,
                        'pH2':pH2,
                        'T1':T1,
                        'T2':T2,
                        'res1':res1,
                        'res2':res2,
                        'spacegroup1':spacegroup1,
                        'spacegroup2':spacegroup2,
                        'hetIDs1':','.join(hetIDs1),
                        'hetIDs2':','.join(hetIDs2),
                        'rmsd':rmsd,
                        'chains':n_chains,
                        'residues':n_residues,
                        'coordinates':n_coordinates,
                        'mutations':mutations,
                        }

                    lines_txt = []
                    for key in d_columns_data_txt.keys():
                        value = d_columns_data_txt[key]
                        lines_txt += ['%15s %s\n' %(key, value)]
                    for fd in [
                        'txt/%s/%s.txt' %(pdb1[1],prefix1),
                        'txt/%s/%s.txt' %(pdb2[1],prefix2),
                        ]:
                        fd = open(fd,'w')
                        fd.writelines(lines_txt)
                        fd.close()
                    
                    ## write data to html lines
                    tr = '<tr>\n'
                    for column in l_columns_html:
                        tr += '%s</td>\n' %(d_columns_data[column])
                    tr += '</tr>\n'

                    l_tr += tr

                    if not pdb1 in d_html.keys():
                        d_html[pdb1] = ''
                    d_html[pdb1] += tr
                    if not pdb2 in d_html.keys():
                        d_html[pdb2] = ''
                    d_html[pdb2] += tr

    ## write html to global file
    if prefix != 'quickrmsd':
        path = 'tmp/'
        file = '%s.htm' %(prefix)
        append_table_rows(path,file,l_tr,th,d_rmsd,)

    ## write html to local file
    if prefix == 'quickrmsd':
        path = 'htm/'
        for pdb in d_html.keys():
            file = '%s.htm' %(pdb)
            l_tr = d_html[pdb]
            append_table_rows(path,file,l_tr,th,d_rmsd,)
    
    return


def append_table_rows(path, file,l_tr,th,d_rmsd):

    if os.path.isfile('%s%s' %(path,file)):
        fd = open('%s%s' %(path, file),'r')
        lines_prev1 = fd.readlines()
        fd.close()
        for i in range(len(lines_prev1)-1,-1,-1,):
            if lines_prev1[i] in ['<table border="1" class="sortable" id="sortable">\n','<table border="1">\n',]:
                lines_prev1 = lines_prev1[i+1:-1]
                break
            
        n_lines = 1+len(l_columns_html)+1
        lines_prev2 = []

        for i_row in range(1,len(lines_prev1)/(n_lines),1,):
            ## parse pdbs
            pdb_lines = lines_prev1[i_row*n_lines+3:i_row*n_lines+5]
            bm_lines = lines_prev1[i_row*n_lines+5:i_row*n_lines+7]
            l_pdbs = []
            l_bms = []
            for j in range(2):
                pdb_line = pdb_lines[j]
                bm_line = bm_lines[j]
                try:
                    pdb_index2 = pdb_line.index('</a>')
                except:
                    print lines_prev1[i_row*n_lines:(i_row+1)*n_lines]
                    print lines_prev1[0*n_lines:1*n_lines]
                    print lines_prev1[0]
                    print lines_prev1[-1]
                    print path,file, pdb_line
                    print len(lines_prev1), n_lines
                    stop
                pdb_index1 = pdb_line[:pdb_index2].rindex('>')+1
                pdb = pdb_line[pdb_index1:pdb_index2]
                bm_index2 = bm_line.index('</td>')
                bm_index1 = bm_line[:bm_index2].rindex('>')+1
                bm = bm_line[bm_index1:bm_index2]
                l_pdbs += [pdb]
                l_bms += [bm]
            pdb1 = l_pdbs[0]
            pdb2 = l_pdbs[1]
            try:
                bm1 = int(l_bms[0])
            except:
                print lines_prev1[i_row*n_lines:(i_row+1)*n_lines]
                print l_bms
                print bm_lines
                print file
                stop
            bm2 = int(l_bms[1])
            bool_already_in_html = False
            if pdb1 in d_rmsd.keys():
                if bm1 in d_rmsd[pdb1].keys():
                    if pdb2 in d_rmsd[pdb1][bm1].keys():
                        if bm2 in d_rmsd[pdb1][bm1][pdb2].keys():
                            bool_already_in_html = True
            if pdb2 in d_rmsd.keys():
                if bm2 in d_rmsd[pdb2].keys():
                    if pdb1 in d_rmsd[pdb2][bm2].keys():
                        if bm1 in d_rmsd[pdb2][bm2][pdb1].keys():
                            bool_already_in_html = True
            if bool_already_in_html == False:
                lines_prev2 += lines_prev1[i_row*n_lines:(i_row+1)*n_lines]

        html = []
        ## javascript
        html += ['<script type="text/javascript" src="../sortable.js"></script>\n']
        ## table init
        html += ['<table border="1" class="sortable" id="sortable">\n']
        ## table header
        html += lines_prev1[:n_lines]
        ## previous lines
        html += lines_prev2
        ## new/updated lines
        html += [l_tr]
        ## table term
        html += ['</table>\n']
    else:
        html = ['<table border="1">\n'+th+l_tr+'</table>\n']

    fd = open('%s%s'%(path,file),'w')
    fd.writelines(html)
    fd.close()

    return


def dataset_description():

    print 'make sure nothing is written to html files before continuing! continue?'
    s = raw_input()
    if s != 'y':
        return


    ##
    ## table of single point mutations
    ##
    print 'generating table with statistics for single point mutations'
    lines = []
    for s in '0123456789abcddefghijklmnopqrstuvwxyz':
        fd = open('single_point_mutations/%s.txt' %(s),'r')
        lines += fd.readlines()
        fd.close()
    d_counts = {}
    d_counts_wt = {}
    d_counts_mut = {}
    count = 0
    ## temp!!! tmp!!!
    for res in ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR',]: ## temp!!! tmp!!!
        d_counts_wt[res] = 0 ## temp!!! tmp!!!
        d_counts[res] = 0 ## temp!!! tmp!!!
        d_counts_mut[res] = 0 ## temp!!! tmp!!!

    l_pdbs = []
    for line in lines:
        l = line.split()

        res_name1 = l[8]
        res_name2 = l[9]
        if res_name1 == 'N/A':
            continue
        if res_name2 == 'N/A':
            continue

        ## avoid duplicate entries (that shouldnt have been there in the first place!!!)
        pdb1 = l[0]
        pdb2 = l[1]
        if [pdb1,pdb2,] in l_pdbs:
            continue
        else:
            l_pdbs += [[pdb1,pdb2,]]

        mutation = res_name1+res_name2
        if not mutation in d_counts.keys():
            d_counts[mutation] = 0
        d_counts[mutation] += 1
        count += 1

        if not res_name1 in d_counts_wt.keys():
            d_counts_wt[res_name1] = 0
        d_counts_wt[res_name1] += 1
        if not res_name2 in d_counts_mut.keys():
            d_counts_mut[res_name2] = 0
        d_counts_mut[res_name2] += 1

    print 'total number of single point mutations:', count
    s = '   \t'
    for res in ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR',]:
        s += '%5s\t' %(mutation)
    s += 'total\n'
    for res1 in ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR',]:
        s += '%5s\t' %(res1)
        for res2 in ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR',]:
            if res1 == res2:
                s += '\t'
            elif not res1+res2 in d_counts.keys():
                s += '\t'
            else:
                s += '%5i\t' %(d_counts[res1+res2])
        s += '%5i' %(d_counts_wt[res1])
        s += '\n'
    s += 'total\t'
    for res2 in ['ALA','CYS','ASP','GLU','PHE','GLY','HIS','ILE','LYS','LEU','MET','ASN','PRO','GLN','ARG','SER','THR','VAL','TRP','TYR',]:
        s += '%5i\t' %(d_counts_mut[res2])
    s += '%5i' %(count)
    s += '\n'
    fd = open('mutation_table.txt','w')
    fd.write(s)
    fd.close()
##        stop


    ##
    ## statistics
    ##
    print 'generating statistics'

    htmls = os.listdir('htm/')
    n_lines = 1+len(l_columns_html)+1
    d_values = {
        'rmsd':[],
        'chains':{},
        'residues':[],
        'coordinates':[],
        'mutations':{},
        }
    d_parameters = {
        'rmsd':7,
        'mutations':8,
        'chains':9,
        'residues':10,
        'coordinates':11,
        }
    for j in range(len(htmls)):
        if j % 100 == 0:
            print '%4.1f%%' %(100*float(j)/len(htmls))
        html = htmls[j]
        pdb = html[:4]
        fd = open('htm/%s' %(html),'r')
        lines = fd.readlines()
        fd.close()
        if lines[0][:8] == '<script ':
            lines = lines[2+n_lines:-1]
        else:
            lines = lines[1+n_lines:-1]
        for i in range(len(lines)/(n_lines)):
            for parameter in ['mutations','residues','coordinates','chains','rmsd',]:
                i_line = d_parameters[parameter]
                s_line = lines[i*n_lines+i_line]
                index2 = s_line.rindex('</td>')
                index1 = s_line[:index2].rindex('>')+1
                value = s_line[index1:index2].strip()
                if parameter == 'rmsd':
                    if float(value) < 10.:
                        d_values['rmsd'] += [float(value)]
                    else:
                        break
                if parameter == 'chains':
                    chains = int(value)
                    if not int(value) in d_values['chains'].keys():
                        d_values['chains'][chains] = 0
                    d_values['chains'][chains] += 1
                    if chains == 1:
                        if not mutations in d_values['mutations'].keys():
                            d_values['mutations'][mutations] = 0
                            if mutations >= 60:
                                print html, htmls[j-1]
                                stop
                        d_values['mutations'][mutations] += 1
                        d_values['residues'] += [int(residues)]
                        d_values['coordinates'] += [coordinates]
                if parameter == 'mutations':
                    mutations = int(value)
                if parameter == 'residues':
                    residues = int(value)
                if parameter == 'coordinates':
                    coordinates = int(value)
    print 'avg rmsd', sum(d_values['rmsd'])/len(d_values['rmsd'])
    print 'avg residues', sum(d_values['residues'])/len(d_values['residues'])
    print 'avg coordinates', sum(d_values['coordinates'])/len(d_values['coordinates'])
    for k,v in d_values['chains'].items():
        print 'chains', k,v
    for k,v in d_values['mutations'].items():
        print 'mutations', k,v


    ## number of paired PDBs
    l_fn = os.listdir('%s/htm' %(self.topdir))
    print 'number of paired *pdb files*:', len(l_fn)

    ## protein with most partners
    l_fn = os.listdir('%s/htm' %(self.topdir))
    size_max = [0,'N/A',]
    for fn in l_fn:
        if os.path.getsize('%s/htm/%s' %(self.topdir,fn,)) > size_max[0]:
            size_max = [os.path.getsize('%s/htm/%s' %(self.topdir,fn,)),fn,]
    fd = open('%s/htm/%s' %(self.topdir,size_max[1],),'r')
    lines = fd.readlines()
    fd.close()
    i = -1
    for line in lines:
        if '<tr>' in line:
            i += 1
    print '%s has %i partners (incl. structures with multiple biounits)' %(size_max[1],i,)

    ## number of paired biomolecules
    n_pairs = 0
    n_biomolecules = 0
    for s in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789':
        l_fn = os.listdir('%s/pdb/%s' %(self.topdir,s.lower(),))
        n_pairs += len(l_fn)
        set_biomolecules = set()
        for s in l_fn:
            set_biomolecules |= set([s[:6]])
            set_biomolecules |= set([s[6:12]])
        n_biomolecules += len(set_biomolecules)
    print 'number of *biomolecules*', n_biomolecules
    print 'number of paired *biomolecules*', n_pairs

    return


def plot_subsets(d_combined):

    d_x1 = {
        'chains_no_transformations':{
            'x1':1.0,
##            'x2':12.0,
            'y1':0,
            'y2':2.0,
            'xlabel':'chains',
            'step_y':0.05,
            },
        }

    print

    for parameter in d_combined.keys():

##        prefix_gnuplot = 'ps/subset_%s' %(parameter)
        prefix_gnuplot = 'ps/discrete_%s' %(parameter)

        ## non-discrete (continuous) data
##        s_averages = ''
        gnuplotdata_fit = ''
        gnuplotdata_errorbars = ''
        l_fit = [[],[],]
        for value in d_combined[parameter].keys():
            l_rmsds = d_combined[parameter][value]
            if len(l_rmsds) == 1:
                print parameter, value, 'only 1 rmsd value',
                stderr = 0
                average = l_rmsds[0]
            else:
                average, stderr = statistics.do_stderr(l_rmsds)
            gnuplotdata_errorbars += '%s %s %s\n' %(value, average, stderr)
            for rmsd in d_combined[parameter][value]:
                gnuplotdata_fit += '%s %s\n' %(value, rmsd)
                l_fit[0] += [value]
                l_fit[1] += [rmsd]

        ## discrete data
        step_x = 1
        step_y = d_x1[parameter]['step_y']

        ## RMSD as a function of the number of mutations
##        if parameter == 'mutations_1_chain':

        l_discrete = []
        for value in d_combined[parameter].keys():
            for rmsd in d_combined[parameter][value]:
                l_discrete += [[value,rmsd,]]

        l_gnuplotdata_contour = continuous2discrete(l_discrete,step_x,step_y,)

##            fd = open('%s_contour.gnuplotdata' %(prefix_gnuplot),'w')
##            fd.writelines(l_gnuplotdata_contour)
##            fd.close()

        fd = open('%s_scatter_fit.gnuplotdata' %(prefix_gnuplot),'w')
        fd.write(gnuplotdata_fit)
        fd.close()

        fd = open('%s_scatter.gnuplotdata' %(prefix_gnuplot),'w')
        fd.write(gnuplotdata_errorbars)
        fd.close()

##        fd = open('ps/%s_averages.gnuplotdata' %(parameter),'w')
##        fd.write(s_averages)
##        fd.close()

        d_xtics = {}
        for value in d_combined[parameter].keys():
            d_xtics[value] = value+0.5

        ##
        ## plot 2d contour (discrete)
        ##
        gnuplot.contour_plot(
            prefix_gnuplot+'_contour', l_gnuplotdata_contour,
            xlabel = d_x1[parameter]['xlabel'], ylabel = 'RMSD',
            title = 'RMSD as a function of %s' %(d_x1[parameter]['xlabel'],),
##            x2 = d_x1[parameter]['x2']+1,
            x2 = max(d_xtics.keys())+1,
            y2 = d_x1[parameter]['y2'],
            x1 = d_x1[parameter]['x1'],
            d_xtics = d_xtics,
            )

        a,b,r,p = statistics.do_regression(l_fit[0],l_fit[1],)
        print '******** r', parameter, r, p

        ##
        ## plot 2d scatter
        ##
        gnuplot.scatter_plot_2d(
            prefix_gnuplot+'_scatter',
            regression=True, regression_data=prefix_gnuplot+'_scatter_fit',
            regression_title = 'slope = %s, r = %s' %(round(b,3), round(r,2)),
            errorbars=True,
            xlabel = d_x1[parameter]['xlabel'],
##                averages = 'ps/%s_averages.gnuplotdata' %(parameter)
            xmin = d_x1[parameter]['x1']-1, ## extend so errorbar doesnt overlap with frame
            xmax = max(d_xtics.keys())+1, ## extend so errorbar doesnt overlap with frame
            ymin = d_x1[parameter]['y1'],
            ymax = d_x1[parameter]['y2']
            )
            
    return


def plot_nominal(d_data_nominal,d_data_nominal_ANOVA,):

    print

    ## no need to plot binary categorical data!!!

    for parameter in d_data_nominal.keys():

        if parameter == 'spacegroup_categorical':
            continue

        l_keys = d_data_nominal[parameter].keys()
        l1 = d_data_nominal[parameter][l_keys[0]]
        l2 = d_data_nominal[parameter][l_keys[1]]
        n = min(len(l1),len(l2))

        l_t = []
        l_p = []
        for N in range(1000):
            l_l_rmsds = [[],[],]
            for i_l_l_rmsds in range(2):
                l_rmsds = [list(l1),list(l2),][i_l_l_rmsds]
                for x in range(n):
                    i_rmsd = int(random.random()*len(l_rmsds))
                    l_l_rmsds[i_l_l_rmsds] += [l_rmsds.pop(i_rmsd)]
            mean1,mean2,stderr,t,p = statistics.twosamplettest(l_l_rmsds[0],l_l_rmsds[1],verbose=False,)
            l_t += [t]
            l_p += [p]
        l_t.sort()
        l_p.sort()
        print
        print parameter, l_keys[0], l_keys[1]
        print l_t[0], l_t[100], l_t[500], l_t[-101], l_t[-1]
        print l_p[0], l_p[100], l_p[500], l_p[-101], l_p[-1]

##            gnuplotdata = ''
##            for rmsd in l1:
##                gnuplotdata += '1 %s\n' %(rmsd)
##            for rmsd in l2:
##                gnuplotdata += '2 %s\n' %(rmsd)
##            prefix_gnuplot = 'ps/%s' %(parameter)
##            fd = open('%s.gnuplotdata' %(prefix_gnuplot),'w')
##            fd.write(gnuplotdata)
##            fd.close()
##            d_xtics = {d_nominal[parameter][0]:1.,d_nominal[parameter][1]:2.}
##            gnuplot.scatter_plot_2d(prefix_gnuplot, d_xtics=d_xtics, xlabel = parameter)

        print
        print parameter, l_keys[0], l_keys[1]
        statistics.twosamplettest(l1,l2)

##    ##
##    ## two-way ANOVA
##    ##
##    d_table = {}
##    for k in d_data_nominal_ANOVA.keys():
##        if 'spacegroup' in k:
##            continue
##
##        l_l_rmsds = []
##        print
##        print '############'
##        print k
##        l_ksub = d_data_nominal_ANOVA[k].keys()
##        l_ksub.sort() ## sort to submit correct matrix for statistical analysis
##        for ksub in l_ksub:
##            l_rmsds = d_data_nominal_ANOVA[k][ksub]
##            if len(l_rmsds) == 0:
##                break
##            l_l_rmsds += [l_rmsds]
##            print ksub, len(l_rmsds), sum(l_rmsds)/len(l_rmsds)
##        if len(l_rmsds) == 0:
##            continue
##
##        ## make lists of equal length for *balanced* ANOVA design
##        n = min(len(l_l_rmsds[0]),len(l_l_rmsds[1]),len(l_l_rmsds[2]),len(l_l_rmsds[3]))
##
##        d_sampling = {'rows':[],'cols':[],'rxc':[],}
##        for i_sample in range(1000):
##            l_l_rmsds_random = [[],[],[],[],]
##            ## randomly populate groups
##            for i_llrmsd in range(4):
##                l_rmsds = list(l_l_rmsds[i_llrmsd])
##                for x in range(n):
##                    i_rmsd = int(random.random()*len(l_rmsds))
##                    l_l_rmsds_random[i_llrmsd] += [l_rmsds.pop(i_rmsd)]
##
##            p_rows, p_cols, p_interaction, l_means, l_variances = statistics.twofactor_anova(
##                l_l_rmsds_random[0],l_l_rmsds_random[1],
##                l_l_rmsds_random[2],l_l_rmsds_random[3],
##                verbose = False,
##                )
####            if p_interaction == None:
####                continue
##
##            d_sampling['rows'] += [p_rows]
##            d_sampling['cols'] += [p_cols]
##            d_sampling['rxc'] += [p_interaction]
##
##        print 'n', n
##        print 'means', l_means
##        print 'variances', l_variances
##        for k_stat in ['rows','cols','rxc']:
##            d_sampling[k_stat].sort()
##            print k_stat,
##            print 'min', round(d_sampling[k_stat][100],4),
##            print 'max', round(d_sampling[k_stat][900],4),
##            print 'median', round(d_sampling[k_stat][500],4)
##
##        ## append to table
##        k1 = k[k.index('+')+1:]
##        k2 = k[:k.index('+')]
##        if not k1 in d_table.keys():
##            d_table[k1] = {}
##        d_table[k1][k2] = {
##            'n':n,
##            'intmin':d_sampling['rxc'][100],
##            'intmax':d_sampling['rxc'][900],
##            'intmed':d_sampling['rxc'][500],
##            'rowmin':d_sampling['rows'][100],
##            'rowmax':d_sampling['rows'][900],
##            'rowmed':d_sampling['rows'][500],
##            'colmin':d_sampling['cols'][100],
##            'colmax':d_sampling['cols'][900],
##            'colmed':d_sampling['cols'][500],
##            }
##
####        print 'p_rows', round(p_rows,4)
####        print 'p_cols', round(p_cols,4)
####        print 'p_interaction', round(p_interaction,4)
##
##    ## sort keys like in Word document
##    l_keys = ['pH','T','authors','hetIDs','cryo','REMARK465','REMARK470',]
##    lines = []
##    for k1 in l_keys:
##        line1 = ''
##        line2 = ''
##        line3 = ''
##        line4 = ''
##        for k2 in l_keys:
##            if k1 != k2:
##                line1 += '%i\t%.2f-%.2f\t' %(d_table[k1][k2]['n'],abs(d_table[k1][k2]['intmin']),abs(d_table[k1][k2]['intmax']))
##                line2 += '\t%.2f\t' %(abs(d_table[k1][k2]['intmed']),)
##                line3 += '%.2f-%.2f\t%.2f-%.2f\t' %(abs(d_table[k1][k2]['rowmin']),abs(d_table[k1][k2]['rowmax']),abs(d_table[k1][k2]['colmin']),abs(d_table[k1][k2]['colmax']))
##                line4 += '%.2f\t%.2f\t' %(abs(d_table[k1][k2]['rowmed']),abs(d_table[k1][k2]['colmed']))
##            else:
##                line1 += '\t\t'
##                line2 += '\t\t'
##                line3 += '\t\t'
##                line4 += '\t\t'
##                
##        line1 = line1[:-1]+'\n'
##        line2 = line2[:-1]+'\n'
##        line3 = line3[:-1]+'\n'
##        line4 = line4[:-1]+'\n'
##        lines += [line1,line2,line3,line4,]
##
##    fd = open('ANOVA_table.txt','w')
##    fd.writelines(lines)
##    fd.close()

    ##
    ## space groups
    ##
    plot_spacegroup_vs_rmsd_vs_count(d_data_nominal)

    return


def plot_spacegroup_vs_rmsd_vs_count(d_data_nominal,):

    print

    ## statistics
    l1 = []
    l2 = []

##    set_spacegroups = set()
##    set_spacegroups |= set(d_data_nominal['spacegroup_categorical']['identical'].keys())
##    set_spacegroups |= set(d_data_nominal['spacegroup']['different'].keys())
##    l_spacegroups = list(set_spacegroups)
##    l_spacegroups.sort()

##    ## temporary associations...
##    d_spacegroups = {}
##    d_spacegroups_reverse = {}
##    for i in range(len(l_spacegroups)):
##        spacegroup = l_spacegroups[i]
##        d_spacegroups[spacegroup] = float(i)
##        d_spacegroups_reverse[i] = spacegroup

    ## identical, different
    for type in ['identical','different',]:
        ## discrete
        d_discrete = {}
        rmsd_discrete_max = 0
        step_y = 0.1
        for spacegroup in d_data_nominal['spacegroup_categorical'][type].keys():

            ## exclude low frequency space groups
##            if spacegroup in [
##                'H 3','P 1','P 31','P 32','P 42 2 2','P 6','P 63 2 2',
##                ]:
##                continue
            if not spacegroup in d_data_nominal['spacegroup_categorical']['identical'].keys():
                continue
            if not spacegroup in d_data_nominal['spacegroup_categorical']['different'].keys():
                continue
            if (
                len(d_data_nominal['spacegroup_categorical']['identical'][spacegroup]) < 50
                or
                len(d_data_nominal['spacegroup_categorical']['different'][spacegroup]) < 50
                ):
                continue

##            xval = d_spacegroups[spacegroup]
            ## data for plot
            for rmsd in d_data_nominal['spacegroup_categorical'][type][spacegroup]:

                ## discrete
                rmsd_discrete = rmsd-rmsd%step_y

                ## append count
                if rmsd_discrete > rmsd_discrete_max:
                    rmsd_discrete_max = rmsd_discrete
##                if not value_discrete in d_discrete.keys():
##                    d_discrete[value_discrete] = {}
##                if not rmsd_discrete in d_discrete[value_discrete].keys():
##                    d_discrete[value_discrete][rmsd_discrete] = 0
##                d_discrete[value_discrete][rmsd_discrete] += 1
                if not spacegroup in d_discrete.keys():
                    d_discrete[spacegroup] = {}
                if not rmsd_discrete in d_discrete[spacegroup].keys():
                    d_discrete[spacegroup][rmsd_discrete] = 0
                d_discrete[spacegroup][rmsd_discrete] += 1

            ## data for statistics
            if type == 'identical':
                l1 += d_data_nominal['spacegroup_categorical'][type][spacegroup]
            elif type == 'different':
                l2 += d_data_nominal['spacegroup_categorical'][type][spacegroup]

        ## new associations of xtics and space groups after deletion of selected low frequency space groups
        l_spacegroups = d_discrete.keys()
        l_spacegroups.sort()
        d_xtics_contour = {}
##        d_xtics_scatter = {}
        for i_spacegroup in range(len(l_spacegroups)):
            spacegroup = l_spacegroups[i_spacegroup]
            ## set xtic
            d_xtics_contour['%-10s' %(spacegroup)] = float(i_spacegroup+.5)
##            d_xtics_scatter[spacegroup] = float(i_spacegroup)

        if bool_normalize == True:
            d_discrete = normalize_counts(d_discrete)

        l_gnuplotdata_contour = []
##            for value_discrete in range(0,int(max(d_discrete.keys())/step_x)+1,):
##        for value_discrete in range(0,int((len(l_spacegroups)+1)/step_x)+1,):

        ##
        ## write data to file
        ##
        for i_spacegroup in range(len(l_spacegroups)):

##            value_discrete *= step_x
            value_discrete = i_spacegroup
            spacegroup = l_spacegroups[i_spacegroup]

            for rmsd_discrete in range(0,int(rmsd_discrete_max/step_y),):

                rmsd_discrete *= step_y

                count = 0
                if spacegroup in d_discrete.keys():
                    if rmsd_discrete in d_discrete[spacegroup].keys():
                        count = d_discrete[spacegroup][rmsd_discrete]

                l_gnuplotdata_contour += [
                    '%s %s %s\n' %(value_discrete, rmsd_discrete, count,),
                    ]

            l_gnuplotdata_contour += ['\n']

        prefix_gnuplot = 'ps/nominal_spacegroups_%s' %(type)
        print 'plotting', prefix_gnuplot

##        ## plot 2d scatter (non-discrete)
##        gnuplot.scatter_plot_2d(
##            prefix_gnuplot,
##            d_xtics = d_xtics_scatter,
##            xlabel='spacegroups',
##            )

        ## plot 2d contour (discrete)
        gnuplot.contour_plot(
            prefix_gnuplot, l_gnuplotdata_contour,
##            xlabel = 'spacegroups',
            ylabel = 'RMSD',
            title = 'RMSD as a function of %s (%s)' %('spacegroups',type,),
            d_xtics = d_xtics_contour,
            x2 = len(l_spacegroups)-1,
            y2 = 4.0,
            )

    print
    print 'finished plot', prefix_gnuplot

    ## statistics
    print 'calculating statistics', prefix_gnuplot
    statistics.twosamplettest(l1,l2)

    return


def normalize_counts(d_discrete):

    ## normalize for each x-axis value
    for value_discrete in d_discrete.keys():
        count_x = sum(d_discrete[value_discrete].values())
        max_count_x = max(d_discrete[value_discrete].values())
        for rmsd_discrete in d_discrete[value_discrete].keys():
            count_xy = d_discrete[value_discrete][rmsd_discrete]

##            ## get rid of contour point if few data points at contour point or contour column
##            if count_xy < 10 or len(d_discrete[value_discrete].keys()) < 10:
##                count_xy = 0
##            ## convert count to percentage
##            else:
####                    count_xy = int(100*count_xy/count_x)
##                count_xy = int(100*count_xy/max_count_x)
            count_xy = int(100*count_xy/max_count_x)

            d_discrete[value_discrete][rmsd_discrete] = count_xy

    return d_discrete


def plot_discrete(d_data_ratio_discrete,):

    print

    d_xmax = {
        'chains':{
##            'xmax':12.0,
            'xmax':None,
            'xmin':1.,'step_y':0.05,'ymax':2.0,
            },
        'mutations':{
##            'xmax':10.0,
            'xmax':None,
            'xmin':0.,'step_y':0.02,
##            'ymax':0.5, ## if incl. T4L and HEWL
            'ymax':1.0, ## if excl. T4L and HEWL
            },
        }

    for parameter in d_data_ratio_discrete.keys():

        ##
        ## scatter
        ##
        gnuplotdata_errorbars = ''
        gnuplotdata_fit = ''
        l_fit = [[],[],]
        for value in d_data_ratio_discrete[parameter].keys():
            l_rmsds = d_data_ratio_discrete[parameter][value]
            if len(l_rmsds) == 1:
                print parameter, value, 'only 1 rmsd value',
                stderr = 0
                average = l_rmsds[0]
            else:
                average, stderr = statistics.do_stderr(l_rmsds)
            gnuplotdata_errorbars += '%s %s %s\n' %(value, average, stderr)
            for rmsd in d_data_ratio_discrete[parameter][value]:
                gnuplotdata_fit += '%s %s\n' %(value, rmsd)
                l_fit[0] += [value]
                l_fit[1] += [rmsd]

        prefix_gnuplot = 'ps/discrete_%s' %(parameter)
        fd = open('%s_scatter.gnuplotdata' %(prefix_gnuplot),'w')
        fd.write(gnuplotdata_errorbars)
        fd.close()
        fd = open('%s_scatter_fit.gnuplotdata' %(prefix_gnuplot),'w')
        fd.write(gnuplotdata_fit)
        fd.close()

        print 'plotting scatter', prefix_gnuplot

        a,b,r,p = statistics.do_regression(l_fit[0],l_fit[1],)
        print '******** r', parameter, r, p

        l_averages = [[],[]]
        for value_discrete, l_rmsds in d_data_ratio_discrete[parameter].items():
##                if len(l_rmsds) < 1:
##                    continue
            rmsd_average = sum(l_rmsds)/len(l_rmsds)
            l_averages[0] += [value_discrete]
            l_averages[1] += [rmsd_average]
        a,b,r,p = statistics.do_regression(l_averages[0],l_averages[1],)
        print '1******** r', parameter, type, r, p

        l_averages = [[],[]]
        for value_discrete, l_rmsds in d_data_ratio_discrete[parameter].items():
            if len(l_rmsds) < 2:
                continue
            rmsd_average = sum(l_rmsds)/len(l_rmsds)
            l_averages[0] += [value_discrete]
            l_averages[1] += [rmsd_average]
        a,b,r,p = statistics.do_regression(l_averages[0],l_averages[1],)
        print '2******** r', parameter, type, r, p

        l_averages = [[],[]]
        for value_discrete, l_rmsds in d_data_ratio_discrete[parameter].items():
            if len(l_rmsds) < 3:
                continue
            rmsd_average = sum(l_rmsds)/len(l_rmsds)
            l_averages[0] += [value_discrete]
            l_averages[1] += [rmsd_average]
        a,b,r,p = statistics.do_regression(l_averages[0],l_averages[1],)
        print '3******** r', parameter, type, r, p

        l_averages = [[],[]]
        for value_discrete, l_rmsds in d_data_ratio_discrete[parameter].items():
            if len(l_rmsds) < 10:
                continue
            rmsd_average = sum(l_rmsds)/len(l_rmsds)
            l_averages[0] += [value_discrete]
            l_averages[1] += [rmsd_average]
        a,b,r,p = statistics.do_regression(l_averages[0],l_averages[1],)
        print '10******** r', parameter, type, r, p

        l_averages = [[],[]]
        for value_discrete, l_rmsds in d_data_ratio_discrete[parameter].items():
            if len(l_rmsds) < 100:
                continue
            rmsd_average = sum(l_rmsds)/len(l_rmsds)
            l_averages[0] += [value_discrete]
            l_averages[1] += [rmsd_average]
        a,b,r,p = statistics.do_regression(l_averages[0],l_averages[1],)
        print '100******** r', parameter, type, r, p

        ##
        ## scatter
        ##
        gnuplot.scatter_plot_2d(
            prefix_gnuplot+'_scatter', regression=True,
            errorbars = True,
            xlabel = parameter,
            ymin = 0.,
            ymax = d_xmax[parameter]['ymax'],
            xmin = d_xmax[parameter]['xmin']-1,
##            xmax = d_xmax[parameter]['xmax']+1,
            xmax = max(d_data_ratio_discrete[parameter].keys())+1,
            regression_data=prefix_gnuplot+'_scatter_fit',
            regression_title = 'slope = %s, r = %s' %(round(b,3), round(r,2)),
            )

        ##
        ##
        ##
        l_discrete = []
        for value in d_data_ratio_discrete[parameter].keys():
            for rmsd in d_data_ratio_discrete[parameter][value]:
                l_discrete += [[value,rmsd,]]
        step_x = 1.
        step_y = d_xmax[parameter]['step_y']
        l_gnuplotdata_contour = continuous2discrete(l_discrete,step_x,step_y,)

        d_xtics = {}
        for value in d_data_ratio_discrete[parameter].keys():
            d_xtics[value] = value+0.5

        zlabel = 'count'
        if bool_normalize == True:
            zlabel = 'normalized frequency (%)'
        gnuplot.contour_plot(
            prefix_gnuplot+'_contour', l_gnuplotdata_contour,
            xlabel = parameter,
            ylabel = 'RMSD',
            title = 'RMSD as a function of %s' %(parameter,),
##            y2 = d_steps[parameter]['y2'],
##            x2 = d_steps[parameter]['x2'], x1 = d_steps[parameter]['x1'],
##            zlabel = 'count',
            zlabel = zlabel,
            x1 = d_xmax[parameter]['xmin'],
##            x2 = d_xmax[parameter]['xmax']+1,
            x2 = max(d_xtics.keys())+1,
            y2 = d_xmax[parameter]['ymax'],
            d_xtics = d_xtics,
            )

    return


def continuous2discrete(l_continuous,step_x,step_y,):

    d_discrete, max_value_discrete_y = continuous2discrete_loopandcount(l_continuous,step_x,step_y,)

    if bool_normalize == True:
        d_discrete = normalize_counts(d_discrete)

    l_gnuplotdata_contour = continuous2discrete_loopandfillmatrixoflines(
        d_discrete, step_x,step_y,
        max_value_discrete_y,
        )

    return l_gnuplotdata_contour


def continuous2discrete_loopandfillmatrixoflines(
    d_discrete,step_x,step_y,
    max_value_discrete_y,
    ):

    l_gnuplotdata_contour = []
    for value_discrete_x in range(0,int(max(d_discrete.keys())/step_x)+1,):
        value_discrete_x *= step_x
        for value_discrete_y in range(0,int(max_value_discrete_y/step_y),):
            value_discrete_y *= step_y

            count = 0
            if value_discrete_x in d_discrete.keys():
                if value_discrete_y in d_discrete[value_discrete_x].keys():
                    count = d_discrete[value_discrete_x][value_discrete_y]

            l_gnuplotdata_contour += ['%s %s %s\n' %(value_discrete_x, value_discrete_y, count,)]

        l_gnuplotdata_contour += ['%s %s %s\n' %(value_discrete_x, value_discrete_y+step_y, count,)]
        l_gnuplotdata_contour += ['\n']

    for value_discrete_y in range(0,int(max_value_discrete_y/step_y),):
        l_gnuplotdata_contour += ['%s %s %s\n' %(value_discrete_x+step_x, value_discrete_y, count,)]
    l_gnuplotdata_contour += ['%s %s %s\n' %(value_discrete_x+step_x, value_discrete_y+step_y, count,)]
    l_gnuplotdata_contour += ['\n']

    return l_gnuplotdata_contour


def continuous2discrete_loopandcount(l_continuous,step_x,step_y,):

    max_value_discrete_y = 0
    d_discrete = {}
    for i in range(len(l_continuous)):

        value_continuous_x = l_continuous[i][0]
        value_continuous_y = l_continuous[i][1]
        value_discrete_x = continuous2discrete_value(value_continuous_x,step_x,)
        value_discrete_y = continuous2discrete_value(value_continuous_y,step_y,)

        if not value_discrete_x in d_discrete.keys():
            d_discrete[value_discrete_x] = {}
        if not value_discrete_y in d_discrete[value_discrete_x].keys():
            d_discrete[value_discrete_x][value_discrete_y] = 0
        d_discrete[value_discrete_x][value_discrete_y] += 1

        if value_discrete_y > max_value_discrete_y:
            max_value_discrete_y = value_discrete_y

    return d_discrete, max_value_discrete_y


def continuous2discrete_value(value_continuous,step,):

    value_discrete = value_continuous-value_continuous%step

    return value_discrete


def plot_continuous(d_data_ratio_continuous,):

    print

    d_steps = {
        ## 0.1 steps
        'pH':{
            'xlabel':'pH',
            'xstep': 0.5,'ystep':0.05,
            'x1': 0.0,'x2':  4.0,
            'y2':1.0,
            },
        ## 0.01 (minor) and 0.05 (major1) and 0.1 (major2) steps
        'res':{
            'xlabel':'resolution ({\305})',
            'xstep': 0.1,'ystep':0.05,
            'x1': 0.9,'x2':  3.0,
            'y2':1.5,
            },
        ## continuous
        'T':{
            'xlabel':'T (K)',
            'xstep':10.0,'ystep':0.05,
            'x1': 0.0,'x2':250.0,
            'y2':1.0,
            },
        ## continuous
        'radiusGyration':{
            'xlabel':'r_G ({\305})',
            'xstep': 1.0,'ystep':0.05 ,
            'x1':10.0,'x2': 25.0,
            'y2':2.0,
            },
        ## continuous
        'MV':{
            'xlabel':'MV ({\305}^3 / Da)',
            'xstep': 0.01, ## usually significance of 2 decimals...
##            'xstep': 0.005,
            'ystep':0.05, ## smooth gradient (blue) and square blocks when 0.01 xstep
##            'ystep':0.02,
            'x1':0.0,'x2': 0.25,
            'y2':1.0,
            },
        }

    for parameter in d_data_ratio_continuous.keys():

##        ## skip plot
##        if parameter == 'res':
##            continue

        step_x = d_steps[parameter]['xstep']
        step_y = d_steps[parameter]['ystep']
        for type in d_data_ratio_continuous[parameter]:

            import collections
            values = [
                d_data_ratio_continuous[parameter][type][i][0] for i in range(len(
                    d_data_ratio_continuous[parameter][type]
                    ))
                ]
            print
            print parameter, type, collections.Counter(values).most_common(5), len(values)

            ## skip plot if not "difference"
            if type in ['single','min','max','average',]:
                if parameter == 'res' and type == 'max':
                    pass
##                elif parameter == 'T' and type == 'max':
##                    pass
                elif parameter == 'radiusGyration' and type == 'average':
                    pass
                else:
                    continue
            if parameter in ['res','radiusGyration',] and type == 'difference':
                continue

            gnuplotdata_fit = ''
            l_fit = [[],[],]
            d_discrete_count = {}
            d_discrete = {}
            d_grouped = {}
            rmsd_discrete_max = 0
            for values in d_data_ratio_continuous[parameter][type]:
                ## parse x- and y-axis values
                value = values[0]
                rmsd = values[1]

##                if parameter == 'res' and type == 'max' and value in [1,2,3,]: ## tmp. exclude integer values
##                    continue
##                if parameter == 'pH' and type == 'difference' and value in [0,]: ## tmp. exclude integer values
##                    continue
##                if parameter == 'T' and type == 'difference' and value in [0,]: ## tmp. exclude integer values
##                    continue

                ## append continuous values
                gnuplotdata_fit += '%s %s\n' %(value, rmsd)
                l_fit[0] += [value]
                l_fit[1] += [rmsd]

                ## convert continuous values to discrete values
                value_discrete = value-value%step_x
                rmsd_discrete = rmsd-rmsd%step_y

                ## find max
                if rmsd_discrete > rmsd_discrete_max:
                    rmsd_discrete_max = rmsd_discrete

                ## append count of discrete values
                if not value_discrete in d_discrete_count.keys():
                    d_discrete_count[value_discrete] = {}
                    d_discrete[value_discrete] = []
                    d_grouped[value_discrete] = []
                if not rmsd_discrete in d_discrete_count[value_discrete].keys():
                    d_discrete_count[value_discrete][rmsd_discrete] = 0
                d_discrete_count[value_discrete][rmsd_discrete] += 1
                d_discrete[value_discrete] += [rmsd]
                d_grouped[value_discrete] += [[value,rmsd,]]

            gnuplotdata_errorbars = ''
            for value_discrete in d_discrete.keys():
                l_rmsds = d_discrete[value_discrete]
                if len(l_rmsds) == 1:
                    print parameter, value_discrete, 'only 1 rmsd value',
                    stderr = 0
                    average = l_rmsds[0]
                else:
                    average, stderr = statistics.do_stderr(l_rmsds)
                gnuplotdata_errorbars += '%s %s %s\n' %(value_discrete, average, stderr)

            if bool_normalize == True:
                d_discrete_count = normalize_counts(d_discrete_count)

            l_gnuplotdata_contour = continuous2discrete_loopandfillmatrixoflines(
                d_discrete_count,step_x,step_y,rmsd_discrete_max,
                )

            a,b,r,p = statistics.do_regression(l_fit[0],l_fit[1],)
            print '******** r', parameter, type, r, p

            l_averages = [[],[]]
            for value_discrete, l_rmsds in d_discrete.items():
                if len(l_rmsds) < 1:
                    continue
                rmsd_average = sum(l_rmsds)/len(l_rmsds)
                l_averages[0] += [value_discrete]
                l_averages[1] += [rmsd_average]
            a,b,r,p = statistics.do_regression(l_averages[0],l_averages[1],)
            print '1******** r', parameter, type, r, p

            l_averages = [[],[]]
            for value_discrete, l_rmsds in d_discrete.items():
                if len(l_rmsds) < 2:
                    continue
                rmsd_average = sum(l_rmsds)/len(l_rmsds)
                l_averages[0] += [value_discrete]
                l_averages[1] += [rmsd_average]
            a,b,r,p = statistics.do_regression(l_averages[0],l_averages[1],)
            print '2******** r', parameter, type, r, p

            l_averages = [[],[]]
            for value_discrete, l_rmsds in d_discrete.items():
                if len(l_rmsds) < 3:
                    continue
                rmsd_average = sum(l_rmsds)/len(l_rmsds)
                l_averages[0] += [value_discrete]
                l_averages[1] += [rmsd_average]
            a,b,r,p = statistics.do_regression(l_averages[0],l_averages[1],)
            print '3******** r', parameter, type, r, p

            l_averages = [[],[]]
            for value_discrete, l_rmsds in d_discrete.items():
                if len(l_rmsds) < 10:
                    continue
                rmsd_average = sum(l_rmsds)/len(l_rmsds)
                l_averages[0] += [value_discrete]
                l_averages[1] += [rmsd_average]
            a,b,r,p = statistics.do_regression(l_averages[0],l_averages[1],)
            print '10******** r', parameter, type, r, p

            l_averages = [[],[]]
            for value_discrete, l_rmsds in d_discrete.items():
                if len(l_rmsds) < 100:
                    continue
                print '[[[[[', parameter, value_discrete
                rmsd_average = sum(l_rmsds)/len(l_rmsds)
                l_averages[0] += [value_discrete]
                l_averages[1] += [rmsd_average]
            a,b,r,p = statistics.do_regression(l_averages[0],l_averages[1],)
            print '100******** r', parameter, type, r, p

##            ## normalize before regression
##            print d_grouped.keys()
##            l_a = []
##            l_b = []
##            l_r = []
##            l_p = []
##            n = 1
##            for x in range(1000):
##                l_fit_normalized = [[],[],]
##                for value_discrete in d_grouped.keys():
##                    if len(d_grouped[value_discrete]) < n:
##                        continue
##                    if x == 0:
##                        print 'xxx', value_discrete
##                    l = list(d_grouped[value_discrete])
##                    for xx in range(n):
##                        i_random = int(random.random()*len(l))
##                        value, rmsd = l.pop(i_random)
##                        l_fit_normalized[0] += [value]
##                        l_fit_normalized[1] += [rmsd]
##                a,b,r,p = statistics.do_regression(l_fit_normalized[0],l_fit_normalized[1],verbose=False,)
##                l_a += [a]
##                l_b += [b]
##                l_r += [r]
##                l_p += [p]
##            l_a.sort()
##            l_b.sort()
##            l_r.sort()
##            l_p.sort()
##            print parameter
##            print 'a', l_a[0], l_a[100], l_a[500], l_a[-100], l_a[-1]
##            print 'slope', l_b[0], l_b[100], l_b[500], l_b[-100], l_b[-1]
##            print 'r', l_r[0], l_r[100], l_r[500], l_r[-100], l_r[-1]
##            print 'p', l_p[0], l_p[100], l_p[500], l_p[-100], l_p[-1]
##            if parameter == 'MV' and type == 'difference':
##                stop

            ## write scatter plot data
            prefix_gnuplot = 'ps/continuous_%s_%s' %(parameter,type,)
            fd = open('%s_scatter_fit.gnuplotdata' %(prefix_gnuplot),'w')
            fd.write(gnuplotdata_fit)
            fd.close()
            fd = open('%s_scatter.gnuplotdata' %(prefix_gnuplot),'w')
            fd.write(gnuplotdata_errorbars)
            fd.close()

            ## variables shared between plots
            xlabel = d_steps[parameter]['xlabel']+' '+type

            ## linear regression
            if type in ['difference','single',] or parameter == 'res' or parameter == 'radiusGyration':
                regression = True
            else:
                regression = False

            ##
            ## scatter
            ##
            gnuplot.scatter_plot_2d(
                prefix_gnuplot+'_scatter',
                regression=regression,
                regression_data=prefix_gnuplot+'_scatter_fit',
                regression_title = 'slope = %s, r = %s' %(round(b,3), round(r,2)),
                errorbars=True, ## only errorbars if discrete values
                xlabel = xlabel,
                xmin = -step_x,
##                xmax = max(d_discrete.keys())+step_x,
                xmax = d_steps[parameter]['x2'],
                ymin = 0,
                ymax = d_steps[parameter]['y2'],
                )

            ##
            ## contour
            ##
            zlabel = 'count'
            if bool_normalize == True:
                zlabel = 'normalized frequency (%)'
            gnuplot.contour_plot(
                prefix_gnuplot+'_contour', l_gnuplotdata_contour,
                xlabel = xlabel,
                ylabel = 'RMSD',
                title = 'RMSD as a function of %s' %(xlabel,),
                y1 = 0,
                y2 = d_steps[parameter]['y2'],
                x2 = d_steps[parameter]['x2'],
                x1 = d_steps[parameter]['x1'],
                zlabel = zlabel,
                )

            print
            print 'finished', prefix_gnuplot

    return


def parse_html(
    lines,n_lines,d_pdbs_skip,
    d_data_ratio_discrete,d_data_ratio_continuous,d_data_nominal,
    d_combined,
    html,
    d_rg,d_MV,d_grow,
##    l_gnuplot,
    d_spacegroups_data,
    d_ols,
    d_temperature_data,
    d_authors,
    d_data_nominal_ANOVA,
    ):

    ## loop over table rows
    for i in range(len(lines)/(n_lines)):

        s_gnuplot = ''
        d_nominal = {}

        ## parse pdbs
        pdb1_line = lines[i*n_lines+3]
        pdb2_line = lines[i*n_lines+4]
        pdbs = [pdb1_line,pdb2_line]
        for j in range(2):
            pdb_line = pdbs[j]
            try:
                pdb_index2 = pdb_line.index('</a>')
            except:
                print html
                print pdbs
                print pdb_line
                print lines[i*n_lines:i*(n_lines+1)]
                stop
            pdb_index1 = pdb_line[:pdb_index2].rindex('>')+1
            s_pdb = pdb_line[pdb_index1:pdb_index2]
            pdbs[j] = s_pdb
        pdb1 = pdbs[0]
        pdb2 = pdbs[1]

        if pdb1 < pdb2:
            if not pdb1 in d_pdbs_skip.keys():
                d_pdbs_skip[pdb1] = []
            if pdb2 in d_pdbs_skip[pdb1]:
                continue
            else:
                d_pdbs_skip[pdb1] += [pdb2]
        elif pdb2 < pdb1:
            if not pdb2 in d_pdbs_skip.keys():
                d_pdbs_skip[pdb2] = []
            if pdb1 in d_pdbs_skip[pdb2]:
                continue
            else:
                d_pdbs_skip[pdb2] += [pdb1]

        ## general exclusion
        ## exclude frequently occuring T4L and HEWL (Jens' suggestion to clean/modify dataset)
        if pdb1.upper() in s_exclude or pdb2.upper() in s_exclude:
            continue

        errorpdbs = set([
####                    ## error
####                    '2bfk','2bfl',
####                    ## less than 50 residues due to mutations
####                    '1p7i','1p7j',
##                ## less than 50 residues due to REMARK465 records
##                '2yxq','2yxr','1bx7','1bx8',
##                ## less than 50 residues due to non-overlapping residues
##                '1ft8','1koh',
##                ##
##                '1znb','1zlx','1i3h','1j5o','1f33','1qvc',
##                '1ypu','1uc7','1jr8','1g6l','1ixc','1xri',
##                '1xfi','2q40','2q47','2qd0','1m70',
##                ##
##                '1f1g','1f18','1abs','1dxd','1dvb',
##                '1dmm','1dmn','1dmq',
##                '195l','196l','200l','197l','199l','198l',
##                ## unknown whether error or not
##                '1za1','1qtv','1khh','1cvm','1ez9','1ulx','2der','1x98',

##                ## incorrect (celsius)
##                '1ade','1c0e',

##                ## remediation change
##                '1dos',
##                ## high pH
##                '3c90',
####                                '1fgn',
            ])

        if pdb1 in errorpdbs or pdb2 in errorpdbs:
            continue

        ##
        ## parse non-hyperlinkreferenced data
        ##
        d_parameters = {
            'rmsd':7,
            'mutations':8,
            'chains':9,
            'residues':10,
            'coordinates':11,
            'pH1':12,
            'pH2':13,
            'T1':14,
            'T2':15,
            'res1':16,
            'res2':17,
            'spacegroup1':18,
            'spacegroup2':19,
            'REMARK465':20,
            'REMARK470':21,
            'transformations':22,
            }
        for parameter in d_parameters.keys():
            i_line = d_parameters[parameter]
            s_line = lines[i*n_lines+i_line]
            index2 = s_line.rindex('</td>')
            index1 = s_line[:index2].rindex('>')+1
            value = s_line[index1:index2].strip()
            if parameter in [
                'rmsd',
##                    'pH1','pH2','T1','T2','res1','res2',
                ]:
                try:
                    value = float(value)
                except:
                    print parameter, value
                    print pdb1, pdb2
                    print html
                    stop
            if parameter in ['mutations','chains','residues','coordinates',]:
                try:
                    value = int(value)
                except:
                    print parameter, value
                    print pdb1,pdb2
                    print i
                    print lines[i*n_lines+i_line-1:i*n_lines+i_line+2]
                    print pdb1_line
                    print pdb2_line
                    print html
                    stop
            d_parameters[parameter] = value
        d_parameters['pdb1'] = pdb1
        d_parameters['pdb2'] = pdb2

        ## radius of gyration
##        if pdb1 in d_rg.keys():
        if pdb1 in d_rg.keys() and pdb2 in d_rg.keys():
            rg1 = d_rg[pdb1]
            rg2 = d_rg[pdb2]
            d_parameters['radiusGyration1'] = rg1
            d_parameters['radiusGyration2'] = rg2
##            d_parameters['radiusGyration'] = (float(rg1)+float(rg2)))/2.

        ## append Matthews coefficient, VM (not MV...)
        if pdb1 in d_MV.keys() and pdb2 in d_MV.keys():
            MV1 = d_MV[pdb1]
            MV2 = d_MV[pdb2]
            
##            ## adjust Matthews coefficient for resolution...
##            MV1 = str(float(d_MV[pdb1])+0.584872*(2.-float(d_parameters['res1'])))
##            MV2 = str(float(d_MV[pdb2])+0.584872*(2.-float(d_parameters['res2'])))
            
            d_parameters['MV1'] = MV1
            d_parameters['MV2'] = MV2

        ## skip if virus (might exclude other proteins...)
        if int(d_parameters['chains']) % 60 == 0:
            continue

        ## check resolution
        if len(set([pdb1,pdb2]) & set(['2c32','1ye1','2zqp',])) == 0:
            if d_parameters['res1'] != 'N/A':
                if float(d_parameters['res1']) > 5.:
                    print pdb1,pdb2
                    stop_resolution
            if d_parameters['res2'] != 'N/A':
                if float(d_parameters['res2']) > 5.:
                    print pdb2,pdb1
                    stop_resolution
        if d_parameters['res1'] < '0.5':
            print pdb1, 'res', d_parameters['res1'], pdb2
            continue
        if d_parameters['res2'] < '0.5':
            print pdb2, 'res', d_parameters['res2'], pdb1
            continue

        if float(d_parameters['residues']) < 50:
            errorpdbs = set([
####                    ## error
####                    '2bfk','2bfl',
####                    ## less than 50 residues due to mutations
####                    '1p7i','1p7j',
                ## less than 50 residues due to REMARK465 records
                '2yxq','2yxr','1bx7','1bx8',
                ## less than 50 residues due to non-overlapping residues
                '1ft8','1koh',
                ])
            if len(set([pdb1,pdb2])-errorpdbs) > 0:
                print set([pdb1,pdb2])-errorpdbs
                print pdb1, pdb2
                print d_parameters['residues']
                stop1 ## temp!!!

        ##
        ## parse hyperlinkreferenced data
        ##
        d_hetIDs = {
            'hetIDs1':25,
            'hetIDs2':26,
            }
        
        for hetIDs in d_hetIDs.keys():
            i_line = d_hetIDs[hetIDs]
            s_line = lines[i*n_lines+i_line]
            d_hetIDs[hetIDs] = []
            try:
                index3 = s_line.index('</td>')
            except:
                print pdb1, pdb2
                print s_line
                stop
            if s_line[index3-4:index3] == '</a>':
                index2 = 0
                while index2 != index3-4:
                    index2 += 1
                    index2 += s_line[index2:].index('</a>')
                    index1 = s_line[:index2].rindex('>')+1
                    hetID = s_line[index1:index2].strip()
                    d_hetIDs[hetIDs] += [hetID]

        ## set rmsd
        rmsd = d_parameters['rmsd']

        ##
        ## rmsd outlier investigation
        ##
        if float(rmsd) > 99.: ## temp!!! due to incorrect transformations...
            print rmsd, pdb1,pdb2
            continue
        if float(rmsd) > 4. and d_parameters['spacegroup1'] == d_parameters['spacegroup2']:
##            print html, pdb1, pdb2, rmsd
##            fd = open('case_large_rmsd_identical_spacegroup.txt','a')
##            fd.write("'%s','%s', %f\n" %(pdb1,pdb2,rmsd,))
##            fd.close()
            if d_parameters['chains'] == 1 and d_parameters['mutations'] == 0:
                fd = open('case_large_rmsd_identical_spacegroup_single_chain_no_mutations.txt','a')
                fd.write("'%s','%s', %f\n" %(pdb1,pdb2,rmsd,))
                fd.close()
                print pdb1, pdb2, 'rmsd', rmsd
                print 'stop_interesting__identical_SG_and_1chain_and_0mutations_but_large_RMSD'
##                stop_interesting_identical_SG_single_chain_but_large_RMSD_WTF

        ##
        ## cryo or not?
        ##
        if d_parameters['T1'] not in ['N/A','NULL',] and d_parameters['T2'] not in ['N/A','NULL',]:
            T1 = convert_str_value_to_float_value('T', d_parameters['T1'], d_parameters,)
            T2 = convert_str_value_to_float_value('T', d_parameters['T2'], d_parameters,)

            max_temp_cryo = 243. ## 2zsq 140 2c02 203 1ekv 223 1rdr 243
            min_temp_noncryo = 253. ## 1k7c 263 1g6n 253
            min_temp_cryo = 10. ## 2j8t 10 (just a check...)
            if T1 >= min_temp_noncryo and T2 >= min_temp_noncryo:
                k_temp = 'non-cryo'
            elif T1 >= min_temp_cryo and T1 <= max_temp_cryo and T2 >= min_temp_noncryo:
                k_temp = 'mixed'
            elif T2 >= min_temp_cryo and T2 <= max_temp_cryo and T1 >= min_temp_noncryo:
                k_temp = 'mixed'
            elif T1 >= min_temp_cryo and T1 <= max_temp_cryo and T2 >= min_temp_cryo and T2 <= max_temp_cryo:
                k_temp = 'cryo'
            elif T1 == False or T2 == False:
                k_temp = None
            else:
                print pdb1, pdb2
                print T1, T2
                stop
        else:
            k_temp = None

        ##
        ## regression
        ##
        if (
            d_parameters['mutations'] == 0 and d_parameters['chains'] == 1
            and
            d_parameters['spacegroup1'] == d_parameters['spacegroup2']
            and
            d_parameters['pH1'] != 'NULL' and d_parameters['pH2'] != 'NULL'
            and
            d_parameters['T1'] != 'NULL' and d_parameters['T2'] != 'NULL'
            and
            'radiusGyration1' in d_parameters.keys() and 'radiusGyration2' in d_parameters.keys()
            and
            'MV1' in d_parameters.keys() and 'MV2' in d_parameters.keys()
            ):

            pH1 = convert_str_value_to_float_value('pH', d_parameters['pH1'], d_parameters,)
            pH2 = convert_str_value_to_float_value('pH', d_parameters['pH2'], d_parameters,)
            T1 = convert_str_value_to_float_value('T', d_parameters['T1'], d_parameters,)
            T2 = convert_str_value_to_float_value('T', d_parameters['T2'], d_parameters,)
            res1 = convert_str_value_to_float_value('res', d_parameters['res1'], d_parameters,)
            res2 = convert_str_value_to_float_value('res', d_parameters['res2'], d_parameters,)
            rG_average = (
                float(d_parameters['radiusGyration1'])+float(d_parameters['radiusGyration2'])
                )/2.
            MV_average = (
                float(d_parameters['MV1'])+float(d_parameters['MV2'])
                )/2.

##            if 'MV1' in d_parameters.keys() and 'MV2' in d_parameters.keys(): ## tmp!!!
##                if abs(float(d_parameters['MV1'])-float(d_parameters['MV2'])) > 1.0: ## tmp!!!
##                    print '**************', pdb1, pdb2, MV1, MV2
            
            if k_temp != None:
                d_temperature_data[k_temp] += [rmsd]

            if k_temp == 'cryo':
                pH_diff = abs(float(pH1)-float(pH2))
                res_max = max(float(res1),float(res2))
                T_diff = abs(float(T1)-float(T2))
                if pH1 != False and pH2 != False and T1 != False and T2 != False:
                    d_ols['rmsd'] += [d_parameters['rmsd']]
                    d_ols['radiusGyration'] += [rG_average]
                    d_ols['pHdiff'] += [pH_diff]
                    d_ols['resmax'] += [res_max]
                    d_ols['Tdiff'] += [T_diff]
                    d_ols['MV'] += [MV_average]

        ######
        ## nominal scale data
        ######

        ##
        ## spacegroup
        ##
        spacegroup1 = d_parameters['spacegroup1']
        spacegroup2 = d_parameters['spacegroup2']
        if d_parameters['mutations'] == 0 and d_parameters['chains'] == 1:
            if spacegroup1 == spacegroup2:
                if not spacegroup1 in d_data_nominal['spacegroup_categorical']['identical'].keys():
                    d_data_nominal['spacegroup_categorical']['identical'][spacegroup1] = []
                d_data_nominal['spacegroup_categorical']['identical'][spacegroup1] += [rmsd]
                d_data_nominal['spacegroup_binary']['identical'] += [rmsd]
                d_nominal['spacegroup_binary'] = 'identical'
            else:
                ## add sg1
                if not spacegroup1 in d_data_nominal['spacegroup_categorical']['different'].keys():
                    d_data_nominal['spacegroup_categorical']['different'][spacegroup1] = []
                d_data_nominal['spacegroup_categorical']['different'][spacegroup1] += [rmsd]
                ## add sg2
                if not spacegroup2 in d_data_nominal['spacegroup_categorical']['different'].keys():
                    d_data_nominal['spacegroup_categorical']['different'][spacegroup2] = []
                d_data_nominal['spacegroup_categorical']['different'][spacegroup2] += [rmsd]
                ##
                d_data_nominal['spacegroup_binary']['different'] += [rmsd]
                d_nominal['spacegroup_binary'] = 'different'

        ##
        ## nominal) hetIDs
        ##
##        difference_hetID = abs(len(d_hetIDs['hetIDs2'])-len(d_hetIDs['hetIDs1']))
        set_hetID_diff = set(d_hetIDs['hetIDs1']) ^ set(d_hetIDs['hetIDs2'])
        if (
            d_parameters['mutations'] == 0 and d_parameters['chains'] == 1
            and
            spacegroup1 == spacegroup2
            ):
##            if difference_hetID > 0:
            if len(set_hetID_diff) > 0:
                d_data_nominal['hetIDs']['different'] += [rmsd]
                d_nominal['hetIDs'] = 'different'
            else:
                d_data_nominal['hetIDs']['identical'] += [rmsd]
                d_nominal['hetIDs'] = 'identical'

        ##
        ## nominal) authors
        ##
        bool_authors_identical = False
        if (
            d_parameters['mutations'] == 0 and d_parameters['chains'] == 1
            and
            spacegroup1 == spacegroup2
            ):
            if (
                pdb1 in d_authors.keys()
                and
                pdb2 in d_authors.keys()
                ):
                l_authors1 = d_authors[pdb1]
                l_authors2 = d_authors[pdb2]
                if len(set(l_authors1) & set(l_authors2)) > 0:
                    d_data_nominal['authors']['identical'] += [rmsd]
                    bool_authors_identical = True
                    d_nominal['authors'] = 'identical'
                else:
                    d_data_nominal['authors']['different'] += [rmsd]
                    d_nominal['authors'] = 'different'

##            ##
##            ## nominal) growth conditions
##            ##
##            if (
##                pdb1 in d_grow.keys() and pdb2 in d_grow.keys()
##                and
##                len(d_grow[pdb1]) > 0 and len(d_grow[pdb2]) > 0
##                ):
##                set_grow_comp_diff = set(d_grow[pdb1]) ^ set(d_grow[pdb1])
##                if len(set_grow_comp_diff) == 0:
##                    d_data_nominal['exptl_crystal_grow']['identical'] += [rmsd]
##                    d_nominal['exptl_crystal_grow'] = 'identical'
##                else:
##                    d_data_nominal['exptl_crystal_grow']['different'] += [rmsd]
##                    d_nominal['exptl_crystal_grow'] = 'different'

##            if pdb1 in d_grow.keys() and pdb2 in d_grow.keys():
##                grow1 = d_grow[pdb1].upper().strip()
##                grow2 = d_grow[pdb2].upper().strip()
##                if grow1 not in ['','?',] and grow2 not in ['','?',]:
##                    l_grow = [grow1,grow2]
##                    for i_grow in range(2):
##                        grow = l_grow[i_grow]
##
##                        ## remove trailing punctuation for reasons of comparison
##                        if grow[-1] == '.':
##                            grow = grow[:-1]
##                        grow = grow.replace('AT ROOM TEMPERATURE','')
##                        grow = grow.replace(', ROOM TEMPERATURE','')
##                        ## remove pH information given in _exptl_crystal_grow.pH
##                        while (
##                            grow[:3] == 'PH '
##                            or
##                            'PH=' in grow
##                            or
##                            ', PH ' in grow
##                            or
##                            ',PH ' in grow
##                            or
##                            ' PH ' in grow
##                            or
##                            'TEMPERATURE ' in grow
##                            ):
##                            
##                            if grow[:3] == 'PH ':
##                                index1 = 0
##                                add = 3
##                            ## no space after
##                            elif 'PH=' in grow:
##                                index1 = grow.index('PH=')
##                                add = 0
##                            elif 'TEMPERATURE ' in grow:
##                                index1 = grow.index('TEMPERATURE ')
##                                add = len('TEMPERATURE ')
##                            ## space after due to parenthesis
##                            elif ' PH (' in grow:
##                                index1 = grow.index(' PH (')
##                                add = 4
##                            elif ',PH ' in grow:
##                                index1 = grow.index(',PH ')
##                                if ',PH (' in grow:
##                                    print grow
##                                    stop_tmp
##                                add = 0
##                            else:
##                                index1 = grow.index(' PH ')
##                                add = 4
##
##                            if ',' in grow[index1+add+1:]:
##                                index2 = index1+grow[index1+add+1:].index(',')+add+1
##                            elif ' ' in grow[index1+add+1:]:
##                                index2 = index1+grow[index1+add+1:].index(' ')+add+1
##                            else:
##                                index2 = len(grow)
##    ##                        print
##    ##                        print grow
##                            grow = grow[:index1]+grow[index2:]
##
##                        ## not important
##                        grow = grow.replace('VAPOR DIFFUSION','')
##                        grow = grow.replace('VAPOUR DIFFUSION','')
##                        grow = grow.replace('HANGING DROP','')
##                        grow = grow.replace('SITTING DROP','')
##                        grow = grow.replace('BUFFER','')
##                        ##
##                        grow = grow.replace('-','') ## TRI-SODIUM v TRISODIUM
##                        grow = grow.replace(',',' ') ## comma separator vs space separation
##                        grow = grow.replace(', AND ',', ') ##
####                        grow = grow.replace('.','')
##                        grow = grow.replace('000,','K,') ## PEG 8000 v PEG 8K
##                        grow = grow.replace('000 ','K ') ## PEG 8000 v PEG 8K
##                        grow = grow.replace('POLY ETHYLENE GLYCOL','PEG') ## PEG4000 v PEG 4000
##                        grow = grow.replace('PEG ','PEG') ## PEG4000 v PEG 4000
##    ##                        print grow
##                        while '  ' in grow:
##                            grow = grow.replace('  ',' ') ## double spacing introduced by replacements
##
##                        l_grow[i_grow] = grow
##
##                    grow1 = l_grow[0]
##                    grow2 = l_grow[1]
##
##                    if grow1 == grow2:
##                        d_data_nominal['exptl_crystal_grow']['identical'] += [rmsd]
##                        d_nominal['exptl_crystal_grow'] = 'identical'
##                    else:
##                        d_data_nominal['exptl_crystal_grow']['different'] += [rmsd]
##                        d_nominal['exptl_crystal_grow'] = 'different'
##                        if (
##                            not grow1 in grow2
##                            and
##                            not grow2 in grow1
##                            ):
##                            print
##                            print '**********'
##                            print pdb1, grow1
##                            print pdb2, grow2
##                            print
##                            print d_grow[pdb1].strip()
##                            print
##                            print d_grow[pdb2].strip()
                        

        ##
        ## nominal) remarks, transformations
        ##
        if (
            d_parameters['mutations'] == 0
            and
            spacegroup1 == spacegroup2
            ):
            for parameter in ['REMARK465','REMARK470','transformations']:

##                if parameter == 'transformations' and d_parameters['chains'] != 2:
                if parameter == 'transformations' and d_parameters['chains'] == 1:
                    continue
                if parameter != 'transformations' and d_parameters['chains'] != 1:
                    continue
                
                if d_parameters[parameter] == 'True':
                    d_data_nominal[parameter]['True'] += [rmsd]
                    if parameter == 'transformations':
                        binary_transformations = 1
                    d_nominal[parameter] = 'True'
                elif d_parameters[parameter] == 'False':
                    d_data_nominal[parameter]['False'] += [rmsd]
                    if parameter == 'transformations':
                        binary_transformations = 0
                    d_nominal[parameter] = 'False'

        ## nominal pH diff, T diff
        if (
            d_parameters['mutations'] == 0
            and
            d_parameters['chains'] == 1
            and
            spacegroup1 == spacegroup2
            ):
            if d_parameters['pH1'] != 'NULL':
                if d_parameters['pH1'] == d_parameters['pH2']:
                    d_data_nominal['pH']['identical'] += [rmsd]
                    d_nominal['pH'] = 'identical'
                else:
                    d_data_nominal['pH']['different'] += [rmsd]
                    d_nominal['pH'] = 'different'
            if d_parameters['T1'] != 'NULL':
                if d_parameters['T1'] == d_parameters['T2']:
                    d_data_nominal['T']['identical'] += [rmsd]
                    d_nominal['T'] = 'identical'
                else:
                    d_data_nominal['T']['different'] += [rmsd]
                    d_nominal['T'] = 'different'
            if k_temp == 'cryo':
                d_data_nominal['cryo']['True'] += [rmsd]
                d_nominal['cryo'] = 'True'
            elif k_temp == 'non-cryo':
                d_data_nominal['cryo']['False'] += [rmsd]
                d_nominal['cryo'] = 'False'
            ## ignore mixed (cryo and non-cryo)
            else:
                pass

        ##
        ## ANOVA
        ##
        if (
            d_parameters['mutations'] == 0 and d_parameters['chains'] == 1
            and
            spacegroup1 == spacegroup2
            ):
            for k1,k1sub in d_nominal.items():
                for k2,k2sub in d_nominal.items():
                    if k1 == k2:
                        continue
                    k = k1+'+'+k2
                    ksub = k1sub+'+'+k2sub
                    d_data_nominal_ANOVA[k][ksub] += [rmsd]

        ######
        ## continuous ratio scale data
        ######
        if (
            d_parameters['spacegroup1'] == d_parameters['spacegroup2']
            and
            d_parameters['chains'] == 1
            and
            d_parameters['mutations'] == 0
####            and k_temp == 'cryo' ## tmp!!!
##            and len(set_hetID_diff) == 0 ## tmp!!!
##            and bool_authors_identical == True ## tmp!!!
##            and d_parameters['REMARK465'] == 'False' ## tmp!!!
##            and d_parameters['REMARK470'] == 'False' ## tmp!!!
            ):
            (
                d_data_ratio_continuous, d_diff, d_max,
                ) = parse_continuous(
                    d_data_ratio_continuous, d_parameters,
                    )
##            if (
##                d_data_ratio_continuous['res']['max'][-1][0] > 2.9
##                and
##                d_data_ratio_continuous['res']['max'][-1][0] < 3.1
##                and
##                d_parameters['rmsd'] < 0.05
##                ):
##                print pdb1, pdb2
##                print d_parameters['rmsd']
##                stop

        ######
        ## discrete ratio scale data
        ######
        for parameter in d_data_ratio_discrete.keys():

            if d_parameters['spacegroup1'] != d_parameters['spacegroup2']:
                continue

            if d_parameters['chains'] > 1:
                if parameter == 'chains':
                    pass
                elif parameter == 'mutations':
                    continue
                else:
                    print parameter
                    stop

            if d_parameters['mutations'] > 0:
                if parameter == 'chains':
                    continue
                elif parameter == 'mutations':
                    pass
                else:
                    print parameter
                    stop
            value = d_parameters[parameter]
            if value not in d_data_ratio_discrete[parameter].keys():
                d_data_ratio_discrete[parameter][value] = []
            d_data_ratio_discrete[parameter][value] += [rmsd]

        ######
        ## combined data
        ######

        ## RMSD as a function of the number of chains (when not transformed)
        if (
            d_parameters['mutations'] == 0
            and
            d_parameters['spacegroup1'] == d_parameters['spacegroup2']
            and
            d_parameters['transformations'] == 'False'
            ):
            for s in [
                'chains_no_transformations',
                ]:
                value = d_parameters[s[:s.index('_')]]
                if value not in d_combined[s].keys():
                    d_combined[s][value] = []
                d_combined[s][value] += [rmsd]

        ##
        ## space groups
        ##
        if (
            d_parameters['mutations'] == 0 and d_parameters['chains'] == 1
##            and
##            spacegroup1 != spacegroup2
            ):
            if not spacegroup1 in d_spacegroups_data.keys():
                d_spacegroups_data[spacegroup1] = {}
            if not spacegroup2 in d_spacegroups_data.keys():
                d_spacegroups_data[spacegroup2] = {}
            if not spacegroup2 in d_spacegroups_data[spacegroup1].keys():
                d_spacegroups_data[spacegroup1][spacegroup2] = []
            if not spacegroup1 in d_spacegroups_data[spacegroup2].keys():
                d_spacegroups_data[spacegroup2][spacegroup1] = []
            d_spacegroups_data[spacegroup1][spacegroup2] += [rmsd]
            d_spacegroups_data[spacegroup2][spacegroup1] += [rmsd]

## some attempt to fit to all parameters?
##        if 'pH' in d_diff.keys() and 'T' in d_diff.keys():
##            s_gnuplot = '%f %i %i %i %f %f %f %i %i %f\n' %(
##                rmsd,
##                int(d_parameters['mutations']),
##                int(d_parameters['chains']),
##                int(d_parameters['residues']),
##                d_diff['pH'],
##                d_diff['T'],
##                d_max['res'],
##                binary_transformations,
##                (
##                     0.59483772597790929*int(d_parameters['mutations'])
##                    -0.53773619667629347*int(d_parameters['chains'])
##                    +0.24537238074043699*d_max['res']
##                    -0.084648260889715471*d_diff['pH']
##                    -0.046597220190475144*binary_transformations
##                    -0.0011915915999097186*d_diff['T']
##                    -0.00024641340425029173*int(d_parameters['residues'])
##                    )
##                )
##            l_gnuplot += [s_gnuplot]

    return (
        d_data_ratio_discrete,
        d_data_ratio_continuous,
        d_data_nominal,
        d_pdbs_skip,
        d_combined,
##        l_gnuplot,
        d_spacegroups_data,
        d_ols,
        d_temperature_data,
        d_data_nominal_ANOVA,
        )


def convert_str_value_to_float_value(
    parameter, value,
    d_parameters,
    ):

    if '-' in value: ## e.g. 4lzt pH
        separator = '-'
    if ';' in value:
        separator = ';'
    if '-' in value or ';' in value:
        l_values = value.split(separator)
        while 'NULL' in l_values:
                l_values.remove('NULL')
        while ' NULL' in l_values:
                l_values.remove(' NULL')
        value = 0
        try:
            for s_value in l_values:
                value += float(s_value)
        except:
            print pdb1, pdb2, s_value, l_values
            stop
        if len(l_values) == 0:
            return False
##            break
        else:
            value /= len(l_values)

    try:
        value = float(value)
    except:
        print value
        print no, pdb1, pdb2
        stop

    if (value < 50 or value > 303) and parameter == 'T':
        errorpdbs = set([
            ## high temperature
            '1fah','1q0o',
            ## unknown whether error or not
            '1za1','1qtv','1khh','1cvm','1ez9','1ulx',
            '2der','1x98','1a3d','2v19','2zvz',
            '1f1g','1f18','1dxd','1dvb',
            '1dmm','1dmn','1dmq',
            '195l','196l','200l','197l','199l','198l',
            ## correct (cryogen/cryo-cooled)
            '2fbb','1j3y','1j41','1ajg','1ajh','1a6k','1c1g','2i16','2qxw','2j7n','2j8t',
            '1abs',
            ## correct (cryogen/cryo-cooled)
            '2fbb','1j3y','1j41','1ajg','1ajh','1a6k','1c1g','2i16','2qxw',
            '2j7n','3ecl','3e4n',
##                ## incorrect (celsius)
##                '1ade','1c0e',
##                                ## incorrect (celsius)
##                                '1zin',
            ])
        if len(set([d_parameters['pdb1'],d_parameters['pdb2'],])-errorpdbs) == 2:
            print value
            print pdb1, pdb2
            stop2temperature ## temp!!!
    if (value < 1 or value > 12) and value != 'NULL' and parameter == 'pH':
        errorpdbs = set([
##                                ## high pH
##                                '3c90',
####                                '1fgn',
            ])
        if len(set([pdb1,pdb2])-errorpdbs) == 2:
            print pdb1, pdb2
            print 'pH', value
            print d_parameters['pH1'], d_parameters['pH2']
            stop3pH ## temp!!!
    if value > 11.5 and parameter == 'res':
        if len(set([pdb1,pdb2])-set(['1c1g','2tma',])) == 2:
            print pdb1, pdb2
            print d_parameters['res1'], d_parameters['res2']
            stop4resolution ## temp!!

    return value


def parse_continuous(d_data_ratio_continuous,d_parameters):

    d_diff = {}
    d_max = {}
    for parameter in d_data_ratio_continuous.keys():

        if parameter == 'radiusGyration' and 'radiusGyration1' not in d_parameters.keys():
            continue
        if parameter == 'MV' and 'MV1' not in d_parameters.keys():
            continue

        values = []
        for no in ['1','2']:
            value = d_parameters[parameter+no]
            if value not in ['NULL','N/A']:
                value = convert_str_value_to_float_value(parameter,value,d_parameters,)
                if value == False:
                    continue
                values += [value]
                d_data_ratio_continuous[parameter]['single'] += [[value,d_parameters['rmsd']]]

        if len(values) == 2:

            ## check resolution
            if parameter == 'res' and abs(values[0] - values[1]) > 4.61:
                print abs(values[0] - values[1])
                print pdb1,pdb2
                stop

            if parameter == 'pH' and abs(values[0] - values[1]) > 9.8:
                print abs(values[0] - values[1])
                print pdb1,pdb2
                stop

            d_data_ratio_continuous[parameter]['max'] += [[max(values),d_parameters['rmsd']]]
            d_data_ratio_continuous[parameter]['min'] += [[min(values),d_parameters['rmsd']]]
            d_data_ratio_continuous[parameter]['average'] += [[sum(values)/2.,d_parameters['rmsd']]]
            d_data_ratio_continuous[parameter]['difference'] += [[abs(values[1]-values[0]),d_parameters['rmsd']]]
            d_diff[parameter] = abs(values[1]-values[0])
            d_max[parameter] = max(values)

    return d_data_ratio_continuous, d_diff, d_max
