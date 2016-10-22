## DONT USE REUTERS TICKERS!!!
## Add Greece, Isreael (ta100)

## Largest economies:
## Add S&P1200 global (S&P500+S&P700 international)
## EU (add S&P 350!!!)
## US (Russell, S&P)
## China (SSE50, only interested in large cap)
## Japan (no growth, Nikkei225)
## Germany (DAX30 + ?)
## UK
## France
## Brazil (Bovespa 50+/-)
## Italy (MIB30)
## Portugal (PSI20)
## India (Sensex30)
## Russia (dodgy)
## Canada (TSX60)
## Australia (ASX50)
## South Korea
## Spain (IBEX35)
## Mexico
## Indonesia
## Netherlands (AEX25)
## Turkey
## Switzerland
## Nigeria
## Sweden
## Poland
## Argentina
## Belgium
## Taiwan
## Norway (OBX25)
## Austria
## Thailand
## Colombia
## South Africa
## Denmark (OMXC20)

def main():

    import urllib

    ##    index = []
    ##    for i in range(0,2017,50):
    ##        url = 'http://finance.yahoo.com/d/quotes.csv?s=@%'+'5ENYA&f=sl1d1t1c1ohgv&e=.csv&h=%s' %(i)
    ##        urllines = urllib2.urlopen(url)
    ##        lines = urllines.readlines()
    ##        for line in lines:
    ##            if '"' not in line:
    ##                continue
    ##            ticker = line.split(',')[0][1:-1]
    ##            if ticker not in index:
    ##                index.append(ticker)
    ##    index.sort()
    ##    print index
    ##    stop

##        sp500 = []
##        for i in range(0,500,50):
##            url = 'http://download.finance.yahoo.com/d/quotes.csv?s=@%s5EGSPC&f=sl1d1t1c1ohgv&e=.csv&h=%s' %('%', i)
##            urllines = urllib2.urlopen(url)
##            lines = urllines.readlines()
##            for line in lines:
##                ticker = line.split(',')[0][1:-1]
##                if ticker == '':
##                    continue
##                if ticker not in sp500:
##                    sp500 += [ticker]
##        for ticker in sp500:
##            print "%s," %(ticker),
##        print len(sp500)
##        stop

    ## set tickers
    fd = open('tickers/index_ixic.txt'); IXIC = fd.read().strip().split(','); fd.close()
    fd = open('tickers/index_nyse.txt'); NYSE = fd.read().strip().split(','); fd.close()

    dji30 = [
        'GE','XOM','PG','DD','UTX','AA',
        'MMM','IBM','MRK','AXP','MCD','BA','KO',
        'CAT','DIS','JPM','HPQ','JNJ','WMT','HD','INTC','MSFT','T',
        'PFE','VZ','BAC','CVX','KFT','CSCO','TRV',
        ]
    ## http://www2.standardandpoors.com/portal/site/sp/en/us/page.topic/indices_500/2,3,2,2,0,0,0,0,0,3,1,0,0,0,0,0.html
    fd = open('tickers/index_sp500.txt'); sp500 = fd.read().strip().split(','); fd.close()#; sp500.remove('')
    
    ## http://www.nasdaq.com/indexshares/historical_data.stm
    ## http://dynamic.nasdaq.com/dynamic/nasdaq100_activity.stm
    nasdaq100 = [
        'ATVI','ADBE','AKAM','ALXN','GOOG','GOOGL','AMZN','AAL','AMGN','ADI','AAPL','AMAT','ADSK','ADP','BIDU','BBBY','BIIB','BMRN','AVGO','CA','CELG','CERN','CHTR','CHKP','CSCO','CTXS','CTSH','CMCSA','COST','CSX','CTRP','DISCA','DISCK','DISH','DLTR','EBAY','EA','ENDP','EXPE','ESRX','FB','FAST','FISV','GILD','HSIC','ILMN','INCY','INTC','INTU','ISRG','JD','LRCX','LBTYA','LBTYK','LVNTA','LMCA','LMCK','LLTC','MAR','MAT','MXIM','MU','MSFT','MDLZ','MNST','MYL','NTAP','NTES','NFLX','NCLH','NVDA','NXPI','ORLY','PCAR','PAYX','PYPL','QCOM','REGN','ROST','SBAC','STX','SIRI','SWKS','SBUX','SRCL','SYMC','TMUS','TSLA','TXN','KHC','PCLN','TSCO','TRIP','FOX','FOXA','ULTA','VRSK','VRTX','VIAB','VOD','WBA','WDC','WFM','XLNX','YHOO',        
        ]

    ## http://www.nasdaq.com/screening/company-list.aspx
    with open('tickers/nasdaq_nyse_amex.txt') as f:
        nasdaq_nyse_amex = [line.strip() for line in f]

    ## http://www.russell.com/indexes/membership/default.asp
    fd = open('tickers/index_russell2000.txt'); russell2000 = fd.read().strip().split(','); fd.close()
    fd = open('tickers/index_russell1000.txt'); russell1000 = fd.read().strip().split(','); fd.close()
    fd = open('tickers/index_russell3000.txt'); russell3000 = fd.read().strip().split(','); fd.close()

    dju15 = ['AEP', 'AES', 'CNP', 'D', 'DUK', 'ED', 'EIX', 'EXC', 'FE', 'NI', 'PCG', 'PEG', 'SO', 'FPL', 'WMB']
    djt20 = ['AXB', 'AMR', 'CAL', 'CHRW', 'CNW', 'CSX', 'EXPD', 'FDX', 'GMT', 'JBHT', 'JBLU', 'LSTR', 'LUV', 'NSC', 'OSG', 'R', 'UNP', 'UPS', 'YRCW']
    djc65 = ['AA','MO','AEP','AIG','T','CNP','C','KO','ED','CAL','D','DUK','DD','EIX','EXC','XOM','FE','GE','HPQ','HON','IBM','JBLU','JPM','MCD','MRK','NI','PFE','PCG','PEG','R','SO','UNP','VZ','WMB','MMM','AES','AMR','CAT','CSX','GMT','NSC','OSG','FPL','UTX','DIS','YRCW','AXB','AXP','BA','CHRW','CNW','EXPD','FDX','HD','INTC','JBHT','JNJ','LSTR','MSFT','PG','LUV','UPS','WMT']

    ##
    ## Asia
    ##
    ## Japan/Tokyo ^N225 $US:N225
    ## http://www.nni.nikkei.co.jp/CF/FR/MKJ/list_nikkei_constituents.cfm
    ## https://en.wikipedia.org/wiki/Nikkei_225
    nikkei225 = []
    with open('tickers/nikkei225.txt') as f:
        for line in f:
            nikkei225.append('JP:'+line.rstrip().split('\t')[1])

    ## Korea/Seoul KOSPI KRX 100
    ## http://eng.krx.co.kr/m2/m2_2/m2_2_3/JHPENG02002_03.jsp
    kospi = [
        '005930.KS', ## Samsung Electronics
        '005935.KS', ## Samsung Electronics
        '005380.KS', ## Hyundai Motor
        '005490.KS', ## POSCO
        '012330.KS', ## Hyundai Mobis
        '055550.KS', ## Shinhan Financial Group
        '009540.KS', ## Hyundai Heavy Industries (HHI)
        '051910.KS', ##
        '032830.KS', ##
        '105560.KS', ## KB Financial Group (KB Kookmin Bank)
        '015760.KS', ## Korea Electric Power (KEPCO)
        '017670.KS', ## SK Telecom (SKM.N)
        '096770.KS', ## SK Innovation
        '053000.KS', ## Woori Finance Holdings
        '034220.KS', ## LG Display
        ## Chaebols?
        '003550.KS', ## LG Group
        '034730.KS', ## SK Group --> SK Holdings Co Ltd
        '000660.KS',  # SK Hynix
        '035420.KS',  # Naver
        '028260.KS',  # Cheil Industries
        ]
    ## China/Shanghai SSE50
    ## https://en.wikipedia.org/wiki/SSE_Composite_Index
    sse50 = []
    with open('tickers/sse50.txt') as f:
        for line in f:
            sse50.append(line.rstrip().split('\t')[1]+'.SS')

    ## China/Hong Kong HSI45
    hsi45 = [
        '0001.HK','0002.HK',
        '0005.HK','0011.HK','0023.HK','0388.HK','0939.HK','1398.HK','2318.HK','2388.HK',
        '2628.HK', ## China Life Insurance Company
        '3328.HK','3988.HK',
        '0002.HK','0003.HK','0006.HK',
        '0001.HK','0012.HK','0016.HK','0083.HK','0101.HK','0688.HK',
        '0004.HK',
        '0013.HK', ## Hutchison Whampoa
        '0017.HK','0019.HK',
        '0066.HK', ## MTR Corporation
        '0144.HK','0267.HK','0291.HK','0293.HK','0330.HK',
        '0386.HK', ## Sinopec
        '0494.HK','0551.HK',
        '0700.HK', ## Tencent Holdings Ltd
        '0762.HK','0857.HK','0883.HK',
        '0941.HK', ## China Mobile
        '1088.HK','1199.HK','2038.HK','2600.HK',
        '1833.HK',
        ]
    ## China/Taiwan (http://www.twse.com.tw/en/statistics/statistics_list.php?tm=01&stm=004)
    taiex = [
        '2498.TW', ## HTC Corporation
        ]
    ## Indonesia/Jakarta ^JKSE Jakarta Islamic Index
    jii30 = [
        'BUMI.IJ', ## Bumi Resources
        'TLKM.IJ', ## Telekomunikasi Indonesia
        'UNTR.IJ', ## United Tractors
        'SMGR.IJ', ## Semen Gresik
        'LPKR.IJ', ## Lippo Karawaci
        'MNCN.IJ', ## Media Nusantara Citra
        'UNVR.IJ', ## Unilever Indonesia
        'INTP.IJ', ## Indocement Tunggal Prakarsa
        'ISAT.IJ', ## Indostat
        ## not complete
        ]
    ## Pakistan/Karachi
    kse30 = [
        'OGDC.PA',
        'MCB.PA',
        'NBP.PA',
        'PPL.PA', ## Pakistan Petroleum
        'SCBPL.PA', ## Standard Chartered Bank
        'UBL.PA',
        'PSO.PA', ## Pakistan State Oil
        'ABL.PA', ## Allied Bank Limited
        'KAPCO.PA', ## Kot Addu Power Company
        'BOP.PA', ## Bank of Punjab
        'BAFL.PA', ## Bank Alfalah
        'SNGP.PA', ## Sui Northern Gas Pipelines
        'HUBC.PA', ## Hub Power Company
        'DAWH.PA', ## Dawood Hercules Chemicals Limited
        'HMB.PA', ## Habib Metropolitan Bank
        'LUCK.PA', ## Lucky Cement
        ## not complete
        ]

    ## Singapore/Singapore City ^STI Strait Times Index
    ## https://en.wikipedia.org/wiki/Straits_Times_Index
    sti30 = []
    with open('tickers/sti30.txt') as f:
        for line in f:
            sti30.append(line.rstrip().split()[1]+'.SI')


    ## Thailand
    ## https://en.wikipedia.org/wiki/SET50_Index_and_SET100_Index
    set50 = []
    with open('tickers/set50.txt') as f:
        for line in f:
            set50.append(line.rstrip().split('\t')[0]+'.SET')

    ## Mumbai/India Bombay Sensex
    sensex30 = [
        'ACC.BO','BHEL.BO',
        'DLF.BO',
        'GRASIM.BO', ## Grasim Industries
        'HINDALCO.BO', ## Hindalco Industries
        'INFY.BO',
        'ITC.BO',
        'NTPC.BO', 'ONGC.BO','RCOM.BO',
        'TCS.BO',
        ]
    ## Kuala Lumpur/Malaysia ^KLSE
    # http://www.ftse.com/objects/csv_to_table.jsp?infoCode=bm1c&theseFilters=&csvAll=&theseColumns=MSwyLDMsMjA=&theseTitles=&tableTitle=FTSE%20Bursa%20Malaysia%20KLCI&dl=&p_encoded=1
    klci = [
        'IJM.KL',
        'AMM.KL', ## 1015 AMMB HOLDINGS BHD
        'CIMB.KL',  # 1023 CIMB Group Holdings
        'RHBC', ## 1066
        'MAYBANK.KL',  # 1155 Malayan Banking Berhad
        '1295.KL',
        'BJTOTO.KL',  # 1562
        'IGB.KL',  # '1597.KL',
        '1619.KL','1643.KL','1651.KL','1783.KL',
        'BURSA.KL', ## 1818.KL, BURSA MALAYSIA
        '1961.KL','2003.KL','2011.KL','2194.KL','2267.KL','2283.KL','2356.KL','2445.KL','2488.KL','2771.KL','2836.KL','2887.KL','3034.KL',
##        '3158.KL',
##        '3182.KL',
        '3492.KL',
        '3794.KL','3816.KL',
        '3859.KL','3867.KL','3905.KL','4006.KL','4065.KL',
        'ROTH.KL', ## 4162 British American Tobacco
        '4197.KL','4235.KL','4243.KL','4324.KL','4405.KL','4502.KL','4588.KL','4634.KL','4677.KL','4863.KL','4898.KL','5005.KL','5012.KL','5014.KL','5052.KL','5053.KL','5060.KL','5076.KL','5077.KL','5089.KL','5097.KL','5099.KL','5103.KL','5122.KL','5142.KL','5185.KL','5231.KL','5266.KL','5304.KL','5347.KL','5398.KL',
        'PETGAS.KL', ## 6033 Petronas Gas
        '5703.KL','5819.KL','5843.KL',
        'PETDAG.KL', ## 5681 Petronas Dagangan
        '6084.KL','6165.KL','6289.KL','6327.KL','6521.KL','6556.KL','6645.KL','6807.KL','6866.KL',
        'AXIATA.KL', ## 6888 Axiata Group
        'DIGI.KL', ## 6947 Digi.com
        '7100.KL','7108.KL',
        'TPGC.KL', # 7113 Top Glove Corporation
        '7158.KL','7164.KL',
        'DIAL.KL', # 7277 Dialog Group
        '8575.KL','8583.KL','8664.KL','8893.KL','9059.KL',
        'WCT.KL', ## 9679 WCT BHD
        ]

    ##
    ## Europe
    ##
    ## UK/London FTSE100
    ftse100 = [
        'GB:AAL','GB:ABF','GB:ANTO','GB:AV.','GB:AZN','GB:BARC','GB:BATS',
        'GB:BBY','GB:BLND','GB:BLT','GB:BNZL',
        'GB:CCL','GB:CNA','GB:CNE','GB:COB','GB:CPG','GB:CPI',
        'GB:DGE','GB:DRX','GB:EMG','GB:EXPN','GB:FGP',
        'GB:GFS',
        'GB:GSK', ## GlaxoSmithKline
        'GB:HMSO','GB:HOME',
        'GB:HSBA', ## HSBC
        'GB:IAP','GB:IHG','GB:JMAT',
        'GB:KAZ','GB:KAZ','GB:KGF','GB:KGF','GB:LAND','GB:LAND','GB:LGEN',
        'GB:LGEN',
        'GB:LLOY','GB:LSE','GB:MKS','GB:MRW',
        'GB:NXT','GB:OML','GB:PNN','GB:PRU','GB:PSON','GB:RB.',
        'GB:RDSa','GB:RDSb','GB:REL',
        'GB:RIO',
        'GB:RRS',
        'GB:RSA','GB:SAB','GB:SBRY','GB:SDR','GB:SDRC','GB:SGE','GB:SHP',
        'GB:SMIN','GB:SRP','GB:SSE','GB:STAN','GB:SVT',
        'GB:TATE','GB:TLW',
        'GB:TSCO',
        'GB:ULVR',
        'GB:VED',
        'GB:VOD', ## Vodafone
        'GB:WOS','GB:WPP','GB:WTB',
        'GB:ADN','GB:ADM','GB:AAL','GB:ANTO','GB:ARM',
        'GB:ABF','GB:AZN','GB:AV.','GB:BAB','GB:BARC',
        'GB:BLT','GB:BP.','GB:BLND','GB:BNZL',
        'GB:BRBY','GB:CPI','GB:CNA','GB:CCH','GB:CPG','GB:CRH',
        'GB:CRDA','GB:DGE',
        'GB:EVR','GB:EXPN','GB:FRES','GB:GFS',
        'GB:GKN','GB:GSK','GB:GLEN','GB:HMSO','GB:HSBA','GB:IMI',
        'GB:IHG','GB:IAG','GB:ITRK','GB:ITV','GB:SBRY','GB:JMAT',
        'GB:KGF','GB:LAND','GB:LGEN','GB:LLOY','GB:MKS','GB:MGGT','GB:MRO',
        'GB:MRW','GB:NXT','GB:OML',
        'GB:PFC','GB:PRU',
        'GB:RRS','GB:RB.','GB:REL',
        'GB:RIO',
        'GB:RBS','GB:RDSA','GB:RSA','GB:SAB','GB:SGE','GB:SDR','GB:SRP',
        'GB:SVT',
        'GB:SNN','GB:SMIN','GB:SSE','GB:STAN',
        'GB:TATE',
        'GB:TSCO',  # Tesco PLC
        'GB:TLW','GB:ULVR',
        'GB:VED',
        'GB:VOD','GB:WEIR','GB:WTB','GB:WOS','GB:WPP',
        ]
    ftse250 = [
        'GB:ADN','GB:ACA',
        'GB:APF','GB:ASHM','GB:AHT','GB:ATK',
        'GB:AVV',
        'GB:AA.','ACA.L','ADN.L','AGR.L','ALD.L',
        'AMFW.L', 'ASHM.L', 'ATK.L','AUTO.L','AVV.L',
        'BAG.L','BBA.L', 'BBY.L', 'BGEO.L','BME.L',
        'GB:BOK','BOY.L','BRSN.L',
        'BTG.L','BVIC.L',
        'BVS.L',
        'BWNG.L','BWY.L','BYG.L','CAPC.L','CARD.L','CBG.L','CCC.L','CEY.L',
        'CHOO.L','CIR.L','CKN.L','CLI.L','CLLN.L','CNE.L','COB.L',
        'CRDA.L','CRST.L',
        'CWD.L','CWK.L','DCG.L',
        'DEB.L',
        'DFS.L','DJAN.L','DLN.L','DNLM.L','DPH.L','DRX.L',
        'DTY.L','ECM.L','ELM.L','EMG.L','ERM.L','ESNT.L',
        'ESUR.L','ETO.L','EVR.L','FDSA.L','FGP.L',
        'GFRD.L','GFS.L','GFTU.L','GNC.L','GNK.L',
        'GNS.L','GOG.L','GPOR.L','GRG.L','GRI.L',
        'HAS.L','HFD.L',
        'HGG.L',
        'GB:HIK','HOME.L','HSTG.L',
        'HSTN.L','HSV.L','HSX.L', 'HWDN.L','IAP.L','IBST.L','ICP.L',
        'IGG.L','IMI.L','INCH.L','INDV.L','INVP.L','IPF.L','IPO.L',
        'IRV.L',
        'GB:JD.','JDW.L','GB:JE.','JLG.L','JLT.L',
        'JUP.L','KAZ.L','KIE.L','KLR.L',
        'LAD.L','LMP.L','LOOK.L',
        'LRE.L','MARS.L',
        'MGAM.L','MGGT.L',
        'MLC.L','MONY.L',
        'MRO.L','MSLH.L','MTO.L','NCC.L','NEX.L','NMC.L','OCDO.L',
        'OPHR.L',
        'OSB.L','PAG.L','PAY.L','PDG.L','PETS.L','PFC.L','PHNX.L',
        'PLP.L','PNN.L','POLY.L','PTEC.L','PZC.L',
        'GB:QQ.','RAT.L','RDW.L','RMV.L','RNK.L',
        'RPC.L',
        'RTO.L','SAFE.L','SAGA.L',
        'SCT.L',
        'SGC.L','SGP.L','SGRO.L',
        'SHAW.L','SHB.L','SHI.L','SMDS.L','SMIN.L','SMP.L',
        'SMWH.L','SNR.L','SOPH.L','SPD.L','SPI.L','SPX.L','SRP.L',
        'SSPG.L', 'SVS.L','SXS.L','TALK.L','TATE.L','TCG.L',
        'TEP.L','TLPR.L',
        'TLW.L',
        'UBM.L','UDG.L','ULE.L','UTG.L','VCT.L','VEC.L','VED.L',
        'VSVS.L','WEIR.L','GB:WG.','WIZZ.L','WKP.L','WMH.L',
        'ZPLA.L',
        ]
    ## London Stock Exchange (LSE)
    with open('tickers/england.txt') as f:
        LSE = ['GB:'+line.strip() for line in f]

    ## France/Paris ^FCHI $US:PARI
    ## http://uk.finance.yahoo.com/q/cp?s=%5EFCHI
    ## https://en.wikipedia.org/wiki/CAC_40
    cac40 = []
    with open('tickers/cac40.txt') as f:
        for line in f:
            cac40.append('FR:'+line.rstrip().split()[0])

    ## Portugal/Lisbon
    ## https://en.wikipedia.org/wiki/PSI-20
    psi20 = [
        'ALTR.LS',
        'BCP.LS',
        'BPI.LS',
        'CTT.LS',
        'EDPR.LS',
        'EDP.LS',
        'GALP.LS',
        'IPR.LS',
        'JMT.LS',
        'EGL.LS',
        'NOS.LS',
        'NVG.LS',
        'PHR.LS',
        'RENE.LS',
        'SEM.LS',
        'SON.LS',
        'TDSA.LS',
        ]

    cacnext20 = [
        'FR:AF',
        'FR:AKE',
        'FR:BVI',
        'FR:CO',
        'FR:CGG',
        'FR:DSY',
        'FR:EDEN',
        'FR:GET',
        'FR:GTO',
        'FR:RMS',
        'FR:LI',
        'FR:MMB',
        'FR:KN',
        'FR:SESG',
        'FR:SCR',
        'FR:SW',
        'FR:SEV',
        'FR:HO',
        'FR:FR',
        'FR:MF',
        ]

    ## http://www.nasdaq.com/screening/companies-by-region.aspx?region=Europe 
    with open('tickers/europe.txt') as f:
        europe = [line.strip() for line in f]

    ## Germany/XETRA ^GDAXI http://en.wikipedia.org/wiki/DAX
    dax30 = [
        'DE:ADS',
        'DE:ALV', ## Allianz
        'DE:BAS', ## BASF
        'DE:BAY', ## Bayer
        'DE:BEI', ## Beiersdorf (NIVEA, Labello)
        'DE:BMW',
        'DE:CBK',
        'DE:DAI',
        'DE:DBK',
        'DE:DB1',
        'DE:LHA',
        'DE:DPW',
        'DE:DTE', ## Deutsche Telekom
        'DE:EOA', ## E.ON
        'DE:FRE', ## Fresenius
        'DE:FME', ## Fresenius Medical Care
        'DE:HEI', ## HeidelbergCement
        'DE:HEN3',  # Henkel
        'DE:IFX',
        'DE:SDF', ## fertilizers
        'DE:LIN',
        'DE:MAN',
        'DE:MRK', ## "German" Merck
##        'DE:MRC', ## "German" Merck
        'DE:MEO',
        'DE:MUV', ## Munich Re
        'DE:RWE',
        'DE:SAP', ## SAP
        'DE:SIE', ## Siemens
        'DE:TKA', ## ThyssenKrupp AG
##        'DE:VOW',
        'DE:VOW3', ## Volkswagen Group
        ]
    mdax50 = [
        'DE:ARL',
        'DE:SPR',
        'DE:BYW6',
        'DE:GBF',
        'DE:BNR',
        'DE:CLS',
        'DE:CON',
        'DE:DEQ',
        'DE:AIR',
        'DE:FIE',
        'DE:FRA',
        'DE:G1A',
        'DE:GXI',
        'DE:GIL', ## Gildemeister AG
        'DE:HNR',
        'DE:HDD',
        'DE:HHFA',
        'DE:HOT',
        'DE:BOS',
        'DE:KCO',
        'DE:KRN',
        'DE:LXS',
        'DE:LEO',
        'DE:NDA',
        'DE:PRA',
        'DE:PSM',
        'DE:PUM',
        'DE:RAA',
        'DE:RHM',
        'DE:RHK',
        'DE:SGL',
        'DE:SZG', ## Salzgitter AG
        'DE:SAZ',
        'DE:SZU',
        'DE:SY',
        'DE:TUI', ## TUI AG, MDAX demotion
        'DE:VOS',
        'DE:WCH',
        'DE:WIN',
##        'DE:CONG', ## MDAX demotion
##        'DE:TUIGn', ## MDAX demotion
        ]
    ## Spain/Madrid ^IBEX
    ## https://en.wikipedia.org/wiki/IBEX_35
    ibex35 = [
        'ES:ABE', ## Abertis
        'ES:ABG',
        'ES:ACS', ## ACS
        'ES:ACX', ## Acerinox
        'ES:ANA',
        'ES:BBVA',
        'ES:BKT',
        'ES:BME',
        'ES:EBRO',
        'ES:ELE',
        'ES:ENG',
        'ES:FCC',
        'ES:FER', ## Ferrovial
        'ES:GAM',
        'ES:GAS', ## Gas Natural
        'ES:GRF',
        'ES:IBE',
        'ES:IDR',
        'ES:ITX', ## Inditex
        'ES:MAP',
        'ES:OHL',
        'ES:POP',
        'ES:REE',
        'ES:REP',
        'ES:SAB',
        'ES:SAN',
        'ES:TEF',
        'ES:TL5',
        'ES:TRE',
        'ES:DIA',
        ]
    with open('tickers/ibex35.txt') as f:
        for line in f:
            ibex35.append('ES:'+line.rstrip().split('\t')[1])

    ## Italy/Milano ^MIB30 .MDD
    ## https://en.wikipedia.org/wiki/FTSE_MIB
    mib30 = [
        'IT:ENI', ## Eni
        'IT:ENEL',
        'IT:ISP',
        'IT:LUX',
        'IT:G',
        'IT:ATL',
        'IT:SRG',
        'IT:UCG',
        'IT:TEN',
        'IT:TIT',
        'IT:TRN',
        'IT:FCA',
        'IT:PST',
        'IT:CNHI',
        'IT:EXO',
        'IT:RACE',
        'IT:FNC',
        'IT:US',
        'IT:MONC',
        'IT:SFER',
        'IT:SPM',
        'IT:ERG','IT:ISP',
        'IT:LUX','IT:MS',
        'IT:PLT','IT:SRG','IT:STM','IT:TRN',
        'IT:FCA',  # Fiat Chrysler Automobiles NV
        'IT:SPM',
        'IT:YNAP',
        'IT:CPR',  # Campari not actually a part of MIB30...
        ]
    ## Australia/AX ^AXJO
    asx50 = [
        ## ASX20
        'AU:MQG',
        ## ASX50
        'AU:AGL',
        'AU:AMC','AU:AMP','AU:ALL','AU:ASX','AU:ANZ','AU:BHP',
        'AU:BSL',
        'AU:BXB',
        'AU:CBA', ## Commonwealth Bank
        'AU:CSL', ## CSL Limited
        'AU:FXJ','AU:GMG','AU:GPT','AU:IAG','AU:LLC',
        'AU:MGR','AU:NAB','AU:NCM','AU:ORI','AU:ORG','AU:QAN',
        'AU:QBE', ## QBE Insurance
        'AU:RIO',
        'AU:STO','AU:SGP','AU:SUN','AU:TAH','AU:TLS',
        'AU:TCL',
        'AU:WES', ## Westfarmers
        'AU:WBC',
        'AU:WPL', ## Woodside Petroleum
        'AU:WOW',
        ]

    ## Belgium/Brussels BEL20
    ## https://en.wikipedia.org/wiki/BEL20
    bel20 = [
        'BE:ABI',
#        'BE:AD',  # Ahold Delhaize
        'BE:ELI',
        'BE:FORB',  # Ageas (BE:AGS)
        'BE:KBC',
        'BE:SOLB',
        'BE:ONTEX',
        'BE:GBLB',
        'BE:UMI',
        'BE:PROX',
        'BE:COLR',
        'BE:THR',
        'BE:ACKB',
        'BE:ACKB','BE:AGFB','BE:BEKB',
        'BE:PROX', ## Belgacom -> Proximus
        'BE:COFB','BE:COLR',
        'BE:DEXB','BE:FORB','BE:GBLB','BE:KBC',
        'BE:SOLB',
        'FR:ENGI', ## GDF Suez -> Engie
        'BE:UMI',
        ]

    ## Canada/Toronto TSX60
    ## https://en.wikipedia.org/wiki/S%26P/TSX_60
    tsx60 = [
        'CA:AEM','CA:AGU','CA:ARX','CA:BMO','CA:BNS','CA:ABX',
        'CA:BCE','CA:VRX','CA:CCO','CA:CM',
        'CA:CNR','CA:CNQ','CA:CP','CA:ENB','CA:ECA',
        'CA:GIL','CA:G','CA:HSE','CA:IMO',
        'CA:K','CA:L','CA:LUN','CA:MGa',
        'CA:MFC', 'CA:NA',
        'CA:SU', ## Suncor Energy
        'CA:POT','CA:BB','CA:RCIb',
        'CA:RY',
        'CA:SJR/B', ## Shaw Communications
        'CA:SNC','CA:SLF',
        'CA:SU','CA:T','CA:TRI',
        'CA:TD','CA:TA','CA:TRP','CA:WN',
        'CA:YRI',
        ]
    with open('tickers/tsx60.txt') as f:
        for line in f:
            tsx60.append('CA:'+line.rstrip().split('\t')[0])

    ## Ireland/Dublin ISEQ20
    ## https://en.wikipedia.org/wiki/ISEQ_20
    iseq20 = [
        'RY4B.IR',
        'CRG.IR',
        'KRZ.IR',
        'BIR.IR',
        'SK3.IR',
        'YZA.IR',
        'KRX.IR',
        'GL9.IR',
        'GCC.IR',
        'GN1.IR',
        'GCC.IR',
        'KRZ.IR', ## Kerry Group
        'KRX.IR',
        ]

    ## OMX Small Cap
    ## http://www.nasdaqomxnordic.com/shares/listed-companies/nordic-small-cap
    with open('tickers/omx_nordic_small_cap.txt') as f:
        omx_small = [line.strip() for line in f]

    ## OMX Medium Cap
    ## http://www.nasdaqomxnordic.com/shares/listed-companies/nordic-mid-cap
    with open('tickers/omx_nordic_medium_cap.txt') as f:
        omx_medium = [line.strip() for line in f]

    ## Sweden/Stockholm OMXS30
    omxs30 = [
        'SE:ABB','SE:ALFA','SE:ASSAb','SE:ATCOa',
        'SE:AZN','SE:BOL','SE:ELUXb','SE:ENRO','SE:ERICb','SE:HMb',
        'SE:INVEb','SE:LUPE','SE:NDA','SE:SAND','SE:SCAb',
        'SE:SEBa','SE:SECUb','SE:SHBa','SE:SKAb',
        'SE:SKF-B','SE:SSABa','SE:SWEDa','SE:SWMA','SE:TEL2b',
        'SE:VOLVb',
        ]
    ## Denmark/Copenhagen OMXC20
    omxc20 = [
        'CARLb.CO', 'CHR.CO', 'COLOb.CO','DANSKE.CO', 'DSV.CO',
        'DNORD.CO',
        'FLS.CO','GEN.CO','LUN.CO','MAERSKb.CO',
        'NDA.CO','NKT.CO','NOVOb.CO','NZYMb.CO','SYDB.CO',
        'TOP.CO','TRYG.CO','VWS.CO','WDH.CO',
        ]
    ## Norway/Oslo OBX25
    obx25 = [
        'AKVER.OL',
        'SUBC.OL',
        'AKER.OL','AKSO.OL','DNBNOR.OL','DNO.OL',
        'FOE.OL','FRO.OL','MHG.OL','NHY.OL',
        'NSG.OL','ORK.OL','PGS.OL',
        'SDRL.OL','SEVAN.OL',
        'STL.OL',
        'STB.OL',
        'TEL.OL','TGS.OL',
        'YAR.OL',
        'GOGL.OL',
        'REC.OL',
        'SCHA.OL',
        ]
    ## Finland/Helsinki OMX Helsinki 25
    omxh25 = [
        'CGCBV.HE','ELISA.HE','FUM1V.HE','KCR1V.HE','KESBV.HE',
        'KNEBV.HE','MEO1V.HE','NES1V.HE','NOK1V.HE',
        'NRE1V.HE','OTE1V.HE','OUT1V.HE','RMR1V.HE',
        'SAA1V.HE','SAMAS.HE','STERV.HE','TIE1V.HE',
        'UNR1V.HE','UPM1V.HE','WRT1V.HE','YTY1V.HE',
        ]
    ## Iceland/Reykjavik OMX Iceland 6 (http://www.nasdaqomxnordic.com/index/index_info?Instrument=ICEIS0000018919)
    omxi6 = [
        'ATLA.IC', ## Atlantic Petroleum
        'ICEAIR.IC', ## Icelandair Group
        'MARL.IC', ## Marel (Food Systems)
        'OSSRu.IC', ## Ossur
        'FOAIR.IC', ## Atlantic Airways
        'BNORDIK.IC',  # BankNordik
##        'HFEIM.IC',
        ## Foroya Banki...
        ]
    ## OMX Nordic 40
    ## http://www.nasdaqomxnordic.com/index/index_info?Instrument=SE0001809476
    omxn40 = [
        'SE:ABB','SE:ALFA', 'ALV', 'SE:ASSAb', 'SE:ATCOa', 'SE:ATCOb',
        'SE:AZN', 'CARLb.CO', 'COLOb.CO', 'DANSKE.CO', 'DSV.CO',
        'ELISA.HE',
        'SE:ELUXb',
        'SE:ERICb', 'FUM1V.HE',  # Fortum
        'GEN.CO', 'SE:HEXAb',
        'SE:HMb','SE:INVEb',
        'KNEBV.HE',
        'MAERSKb.CO','MEO1V.HE','SE:NDA','NES1V.HE','NOK1V.HE',
        'NOVOb.CO','NZYMb.CO','PNDORA.CO', 'SAMAS.HE','SE:SAND','SE:SCAb',
        'SE:SEBa','SE:SECUb','SE:SHBa','SE:SKAb','SE:SKF-B',
        'SE:SSABa','STERV.HE','SE:SWEDa','SE:SWMA',
        'SE:TELIA', 'SE:TEL2b',
        'SE:TELIA','UPM1V.HE','SE:VOLVb','VWS.CO','WRT1V.HE',
        ]
    ## OMX Baltic 10 (http://www.nasdaqomxnordic.com/index/index_info?Instrument=RISSE0001850033)
    ## http://www.nasdaqomxnordic.com/index/index_info?Instrument=SE0001850033
    omxb10 = [
##        'ETLAT.TL', ## AS Eesti Telekom
##        'IVL1L.VL', ## Invalda
##        'SAB1L.VL', ## AB Siauliu Bankas
        'APG1L.VL', ## Apranga
        'BLT1T.TL', ## AS Baltika
        'KNF1L.VL', ## Klaipedos Nafta
        'OEG1T.TL', ## Olympic Entertainment Group AS
        'TAL1T.TL', ## Tallink Grupp
        'TKM1T.TL', ## Tallinna Kaubamaja AS
        'MRK1T.TL',
        'OLF1R.RG',  # Olainfarm AS
        'SAB1L.VL',
        'SFG1T.TL',
        'TEO1L.VL',
        'TVEAT.TL',
        ]

    ## Netherlands/Amsterdam AEX
    ## https://en.wikipedia.org/wiki/AEX_index
    aex25 = [
        'NL:AKZA','NL:ASML',
        'NL:DSM',
        'NL:INGA',
        'NL:KPN',
        'NL:RAND',
        'NL:RDSa',
        'NL:SBMO',
        'NL:TOM2','NL:UNA','NL:MT',
        'NL:GTO', 'NL:ATC', 'NL:WKL', 'NL:PHIA',
        ]

    ## Austria/Vienna ATX20
    atx20 = [
        'TKA.VI',
        'MMK.VI', ## MAYR-MELNHOF KARTON AG
        'VOES.VI',
        'EVN.VI', ## EVN AG
        'VIE.VI',
        'OMV.VI','RHI.VI','WBS.VI','EBS.VI',
        'ANDR.VI','ATEC.VI','BWIN.VI',
        'POST.VI','RIBH.VI','UNIQ.VI',
        ]
    ## Switzerland/Zurich SMI20
    smi20 = [
        'ABBN.VX','ATLN.VX','ADEN.VX','BALN.VX','CSGN.VX',
        'LHN.VX',
        'BAER.VX', ## Julius Baer Group
        'NESN.VX', ## Nestle
        'NOVN.VX', ## Novartis
        'CFR.VX','ROG.VX','UHR.VX','SLHN.VX',
        'SCMN.VX','SYNN.VX',
        'UBSG.VX','ZURN.VX',
        ]
    with open('tickers/switzerland.tsv') as f:
        switzerland = [line.strip() for line in f]
    ## Poland/Warzawa WIG20
    ## http://en.wikipedia.org/wiki/WIG20
    wig20 = [
        'ACP.PW',
        'PEO.PW',
        'BZW.PW',
        'CEZ.PW',
        'CPS.PW',
        'GTN.PW',
        'GTC.PW',
        'LTS.PW',
        'KGH.PW',
        'PBG.PW',
        'PGE.PW',
        'PKN.PW',
        'PKO.PW',
        'PXM.PW',
        'PGN.PW',
        'PZU.PW',
        'TPE.PW',
        ]
    ## Turkey/Istanbul ISE30
    ise30 = [
        'AKBNK.IS','AKGRT.IS','AEFES.IS','BAGFS.IS',
        'DOHOL.IS','ENKAI.IS','EREGL.IS','GARAN.IS',
        'ISCTR.IS','ISGYO.IS','KRDMD.IS','KCHOL.IS','PETKM.IS',
        'SAHOL.IS','SKBNK.IS','SISE.IS','HALKB.IS',
        'TSKB.IS','TAVHL.IS','TKFEN.IS','TOASO.IS','TCELL.IS',
        'TUPRS.IS','THYAO.IS','TTKOM.IS','VAKBN.IS','YKBNK.IS',
        ]

    ## Israel / Tel Aviv / TA100
    ## https://en.wikipedia.org/wiki/TA-100_Index
    ta100 = []
#    with open('tickers/ta100.txt') as f:
#        for line in f:
#            ta100.append(line.rstrip().split('\t')[1]+'.TLV')

    ## Russia/Moscow MICEX10 RTSI50
    micex10 = [
        'VTBR.ME',
        'GAZP.ME', ## Gazprom
        'LKOH.ME','GMKN.ME','ROSN.ME',
        'RTKM.ME', ## OAO Rostelekom
        ]

    ## Riyadh/Saudi Arabia The Tadawul All-Share Index
    tasi = [
        'SABIC.AB', ## Saudi Basic Industries
        ]
   

    ##rm
    ## Africa
    ##
    ## South Africa / Johannesburg
    jse = []
#    with open('tickers/jse.txt') as f:
#        for line in f:
#            jse.append(line.rstrip().split('\t')[0]+'.JNB')
 
    ##
    ## South America
    ##
    ## Mexico City/Mexico IPC (Indice de Precios y Cotizaciones) MXX
    ## https://en.wikipedia.org/wiki/Indice_de_Precios_y_Cotizaciones
    ipc35 = [
        'AMXL.MX',
        'ASURB.MX',
        'AXTELCPO.MX',
        'GCARSOA1.MX','GEOB.MX',
        'GFINBURO.MX','GFNORTEO.MX','GMEXICOB.MX',
        'KIMBERA.MX',
        'SIMECB.MX','SORIANAB.MX',
        ]
    ## Sao Paulo/Brazil Bovespa66 BVSP
    bovespa = [
        'FIBR3.SA', ## Fibria
        'BTOW3.SA', ## B2W Varejo (www.americanas.com)
        'BVMF3.SA','BBDC4.SA','BRAP4.SA',
        'BBAS3.SA', ## Banco do Brasil S.A.
        'BRKM5.SA','CCRO3.SA',
        'CMIG4.SA',
        'CESP6.SA','CGAS5.SA','CPLE6.SA','CSAN3.SA','CPFE3.SA','CYRE3.SA',
        'DTEX3.SA',
        'ELET6.SA',
        'EMBR3.SA','GFSA3.SA','GGBR4.SA','GOAU4.SA','GOLL4.SA',
        'ITSA4.SA',
        'JBSS3.SA','KLBN4.SA','LIGT3.SA','LAME4.SA','LREN3.SA','NATU3.SA',
##        'BNCA3.SA',
        'BRFS3.SA', ## Brasil Foods S.A. (BRF)
        'PETR4.SA', ## PetroBras
        'RSID3.SA','SBSP3.SA',
        'CSNA3.SA',
        'TCSL3.SA',
        'TRPL4.SA',
        'USIM5.SA',
        'VALE5.SA', ## Vale S.A.
        'ELPL4.SA',
        'BRML3.SA',
        'BRPR3.SA',
        'HGTX3.SA',
        'CIEL3.SA',
        'HYPE3.SA',
        'PCAR4.SA',
        'MRFG3.SA',
        ]

    b2c = [
        'WMT','MSFT','AAPL','FDO','AAP','ROST','TJX',
        'FR:RMS','SE:HMb','NKE',
        ## Beverages
        'KO',
        ## Food & Staples Retailing
        'WMT',
        ## Household Products
        'PG',
        ## Wireless Telecommunication Services
        '0941.HK',
        ]

    tommy = [
        'CIB', 'EPD', 'BRP','GSH','AEM','VMED','HW','YPF','DOW','DRYS',
        ## Picks and shovels
        'JEC','SGR','JOYG','BUCY','TWI','CBI',
        ## Metals with infrastructure problems
        'FCX','BHP','RIO','FSUMF','AAUK',
        ## Oil and gas companies expanding production
        'DVN','UPL','MUR','HES','PBR','SU','ECA','CNQ','FPL','SRE',
        ## not CAC40?
        'FR:DSY',
        ## brands
        'WDFC',
        ## exotic
        'IMP.J', ## Johannesburg
        'UPL.KL', ## Kuala Lumpur
        '2866.HK','601866.SS', ## China Shipping Container Lines
        'ALY.L', ## popped ub because aly.usa was delisted
        'CHNG', ## China... seen on CAPS
        ## seen on CAPS
        'HEAT','PWAV','VIP','CHLN','FRPT','SENEA','SWHC','EBIX','BKEP','ADTN','TTM','NFLX',
        'JMHLY',
        ## TOPIX 100 (not Nikkei225)
        'JP:6594', ## Nidec
        'SHLAF', ## Swiss
        ## AMEX
        'LPH',
        ## high 5 year average ROE
        'NRT','PBT','BPT','MSB','SBR','CRT','TNH','DOM',
        'STRA','HGT','FHCO','BFR',
        'VWO','TCCO','PAYX',
        ## Buffett
        '1211.HK','BYDDF','BYDDY','002594.SS',
        ## http://www.cnbc.com/id/44484806
        'DTV','GRA','DVA','LCAPA','VCI','CCOI','CBB','WSFS','FTWR',
        ]

    ## 13F-HR SEC form
    Buffett = [
        'KO', 'USG', 'COP', 'AXP', 'USB', 'WLP', 'NRG', 'WSC', 'UNH', 'COST',
##        'BAC', ## sold
        'WBC', 'GSK', 'TMK', 'ETN', 'HD', 'KMX', 'KFT', 'NSC', 'WPO', 'PG', 'WFC', 'UNP', 'GCI', 'IR', 'STI', 'WMT', 'GE', 'IRM', 'MTB', 'UPS', 'MCO',
##        'NKE', ## sold
        'JNJ', 'SNY',
##        'CMCSA', ## sold
##        'LOW', ## sold
##            'CDCO.B', ## liquidation
##            'BUD', ## acquisition
        'GB:TSCO','FR:SAN','005490.KS',
        '1211.HK', ## BYD Auto, 10% annual letter, 240 million USD wikipedia
##        'BDX', ## sold
##        'NESN.VX', ## sold
##        'NLC', ## sold
        'CEG','OMC',
        'GD','IBM',
        'DVA','DTV','DG','LEE','LMCB','MA','NOV','PSX','V','VRSK',
        'KHC', 'IBM', 'AXP', 'KO',
        'AAPL',  # 2016Q1
        ]

    ## http://www.interbrand.com/en/best-global-brands/best-global-brands-2008/best-global-brands-2010.aspx
    ## http://www.interbrand.com/en/best-global-brands/2012/Best-Global-Brands-2012-Brand-View.aspx
    brands = [
        'KO', ## Coca-Cola, ...
        'AAPL', ## iPhone, iPad
        'IBM', ## IBM
        'GOOG', ## Google
        'MSFT', ## Windows, Office
        'AMZN', ## 20 Amazon
        'GE',
        'MCD', ## McDonalds
        'INTC',
        '005930.KS', ## 9 Samsung Electronics
        '005935.KS', ## 9 Samsung Electronics
        'TM', ## Lexus, Prius
        'DE:DAI', ## Mercedes-Benz
        'DE:BMW',
        'DIS', ## Pixar, Disney
        'CSCO',
        'HPQ',
        'PG', ## Gillette, Duracell, Pampers
        'FR:MC', ## Louis Vuitton, Moet Hennessy (34% Diageo ownership)
        'ORCL', ## 18 Oracle
        'NOK', ## 19 Nokia
        'HMC', ## 21 Honda
        'PEP', ## 22 Pepsi, Gatorade, Pringles
        'SE:HMb', ## 23 H&M
        'V','MA','AXP',
        'JP:7201', ## Nissan
        'GB:BRBY', ## Burberry

        ## Restaurants
        'SBUX','MCD','BKW',
        'YUM', ## Pizza Hut, KFC, Taco Bell

        ## Apparel
        'GPS','SE:HMb',
        ## Levi's (private)
        'ES:ITX', ## Zara, Bershka, Massimo Dutti

        'UHR.VX', ## Swatch Group

        'IT:ENI', ## Benetton
        
        'T', ## AT&T
        'VZ',
        'CHL',
        'VOD',
        'WMT',

        ## Sporting Goods
        'DE:ADS','NKE','DE:PUM',

        ## FMCG
        'JNJ','K','EL',
        'CPB', ## Campbell's
        'FR:BN', ## Danone
        'AVP',
        'CL', ## Colgate
        'BRK.B', ## Heinz
        'KMB', ## Kleenex
        'KFT',
        'FR:OR', ## L'Oreal, Biotherm, The Body Shop (Nestle owns 25%), Lancome
        'NESN.VX', ## Nescafe
        'DE:BEI', ## NIVEA, Labello
        ## Wrigley/Mars
        'KFT',

        ## Electronics        
        'CAJ',

        ## Apparel
        'ANF',

        ## Financial Services
        'GS',

        'K','HOG','EBAY',
        'MAT', ## Barbie
        ## Transportation
        'UPS','FDX'
        

        'BUD', ## Budweiser, Stella Artois, Beck's, Leffe, Hoegaarden
        'DEO', ## Smirnoff, Guinness, Johnnie Walker, Jose Cuervo, Baileys
        '1913.HK', ## Prada
        'FR:KER', ## Gucci, Yves Saint Laurent, PUMA, Fnac
        'DE:ADS',
        'TRI', ## Thomson, Reuters
        'CAT',
        'LTD', ## Victoria's Secret
        'MMM',
        'ADBE',
        'TIF',
        'SNY',
        'MO','PMI', ## Marlboro
        'XRX', ## Xerox
        ]
        
    shipping = ['TRMD','TK','TGP','TOO','OSG','FRO','TNP','DRYS',]
    railroad = ['CSX','NSC','CNI','CP','UNP','GSH','TRN','KSU','WAB','GWR','ARII','RAIL','GBX','PWX','GRWIF',]
    ## 
    sp500_aristocrats = [
        ## Utilities
        'ED',
        ## Technology
        'T', 'ADP',
        ## Basic Materials
        'SHW', 'NUE', 'PPG', 'ECL', 'APD',
        'BMS',
        ## Healthcare
        'BDX', 'BCR', 'MDT', 'ABT', 'JNJ', 'ABBV',
        'CAH',  # since 1989
        ## Financials
        'CINF', 'BEN', 'TROW', 'HCP',
        'AFL',  # Life & Health Insurance
        ## Industrials
        'MMM',  # since 1959
        'EMR',
        'DOV',  # since 1956
        'ITW', 'GWW', 'PNR', 'SWK', 'CTAS',
        'GD',
        ## Consumer Staples
        'ADM', 'BF.B', 'HRL', 'KMB', 'MKC', 'SYY',
        'PG',  # since 1957
        'CL', 'CLX',
        'KO', 'PEP',
        ## Consumer Discretionary
        'GPC',
        'LOW', 'TGT', 'WBA', 'MCD', 'WMT', 'VFC', 'LEG',
        'DG',  # not really, but they acquired FDO
        ## Energy
        'XOM',
        ## Other
        'CB',  # Acquired by ACE, but CB ticker kept?!
##        'LLY','TEG', 2010 deletions (dividend not increased)
##        'SVU', 2010 deletions (dividend reduction)
##        'CTL', ## 2011 deletion
##        'PBI', ## 2013jul12 deletion
        ]
   
   ## http://seekingalpha.com/article/1563962-5-upcoming-dividend-aristocrats-you-may-be-looking-past
    sp500_aristocrats_future = [
        'IBM','NVO',
        'UTX',  # since 1994
        ## http://www.suredividend.com/new-dividend-aristocrats
        'LLTC', 'PX', 'ROP',
        ## http://dividendvaluebuilder.com/dividend-contenders-list/
        'UMBF', 'WABC', 'WTR', 'BMI', 'CBU', 'CPKF', 'FELE', 'JKHY', 'MGRC',
        'O', 'PBCT', 'SYK', 'WST', 'TRI', 'AOS', 'AROW', 'ATR', 'BANF', 'BRO', 'CAT', 'CFR', 'JW.A', 'JW-A', 'MDP', 'PRE', 'PSBQ', 'SKT', 'UBA', 'ALB', 'ROST',
        'ENB', 'PII', 'CAH', 'CHD', 'CNI', 'CA:CNR', 'TJX', 'MUR', 'YORW',
        ]
    ## Fact Sheet http://www.standardandpoors.com/indices/sp-europe-350-dividend-aristocrats/en/us/?indexId=spsdivear-eurew--p-reu---
    ## 10 years of increases...
    sp350_aristocrats = [
        'ES:ABE', ## Abertis Infraestructuras, S.A.
        'SE:ATCOa', ## Atlas Copco AB-A Shares
        'GB:BARC', ## Barclays
        'GB:CPI', ## Capita Group
        'GB:CRH', ## CRH PLC
        'GB:COB', ## Cobham
        'GB:DMGT', ## Daily Mail & General Trust A
        'GB:ETI', ## Enterprise Inns PLC
        'FR:EI', ## Essilor ## EI on Euronext
        'GB:FGP', ## FirstGroup
        'ES:GAS', ## Gas Natural SDG, S.A.
        'GB:HMSO', ## Hammerson
        'FR:RMS', ## Hermes International
        'ES:IBE', ## Iberdrola
        'BE:KBC', ## KBC Group NV
        'GB:LGEN', ## Legal & General Group
        'GB:EMG', ## Man Group PLC
##        'GB:MSY', ## Misys (acquired by Vista Equity Partners)
        'GB:NG.', ## National Grid PLC
        'NESN.VX',
        'NOVN.VX',
        'NOVOb.CO', ## Novo Nordisk A/S-B Shares
        'ORK.OL', ## Orkla ASA
        'FR:PUB', ## Publicis Groupe
        'ROG.VX', ## Roche Holding AG
        'GB:RBS',
        'FR:SAN', ## Sanofi-Aventis ## SAN on Euronext
        'GB:SSE', ## Scottish & Southern Energy
        'SE:SHBa', ## Svenska Handelsbanken
        'BE:UCB', ## UCB S.A.
        'WPP.L', ## WPP Group

##        'ES:ACS',
##        'ES:SAB',
##        'ES:BKT','BATS.L',
##        'FR:BN', ## Danone
##        'GB:DGE',
##        'GBLB.BR','HMB.ST',
##        'FR:OREP', ## OR on Euronext
##        'EMG.L','NXT.L',
##        'OPAP.AT', ## OPAP
##        'PSON.L','RSA.L',
##        'DE:RWE',
##        'GB:RB.','REE.MC',
##        'REL.L', ## Reed-Elsevier (vs John Wiley)
##        'REP.MC',
##        'RDSa.L','RDSb.L',
##        'GB:SVT','SE:SWMA','GB:TATE','GB:TSCO',
##        'DE:TKA', ## dividend cut in 2010
##        'GB:ULVR','UNc.AS', ## Unilever, UNA on Euronext
##        'DG.PA', ## Vinci, DG on Euronext
##        'VIV.PA','VOD.L',
        ]

    ## Google categories
    basicmaterials = [
        'BHP','BBL','VALE','RIO','MT','AAUK','MON','POT','MOS','DD','PKX','FCX','ABX','DOW','SCCO','AA','SID','SSL','GG','PX',
        'KMB','GG','SID','NEM',
        ]
    ## construction
    capitalgoods = [
        'BA','CAT','HON','LMT','MITSY','GD','ITW','DE','TS','NOC','CX','CRH','HIT','PCP','IR','FLR','CMI','COL',
        'ITT', ## from conglomerates
        ]
    conglomerates = [
        'GE','UTX','ABB','MMM','EMR','PHG','RTN','TYC','TXT','FO','LUK','DOV','Y','WSC','TIN','FSS','NOG','MXM',
        'NHY',
        ]
    consumercyclical = [
        'TM','DAI','HMC','SNE','NKE','PCAR','F','COH','HOG','VFC','MGA','MAS','MAT','WHR','GPC','GT',
        ]
    consumernoncyclical = [
        'PG','KO','PEP','PM','UL','UN','BTI','DEO',
        'KHC','MO','CL','SYT','ADM','K','GIS','AVP',
        'RAI','HNZ','CPB',
        ]
    energy = [
        'XOM','PTR','PBR','EC','CVX','TOT','E','COP','SLB',
        '0386.HK', ## SNP
        'CEO','STO','OXY','ECA','IMO','CA:IMO','DVN','SU','RIG',
        'MRO','VLO','CNQ',
        ]
    financial = []
    healthcare = [
        'JNJ','PFE','GSK','SNY','ABT','MRK','GB:AZN','LLY','MDT','GILD','AMGN','BMY','ACL','NVO','TEVA','BAX','SYK','MCK','CAH','ABC',
        ]
    services = [
        'CHL','T','WMT','VOD','IBM','AMX','VZ','FTE','DCM','MCD','CMCSA','DIS','DCM','CVS','TWX','CHA','HD','TGT',
        ]
    technology = [
        'MSFT','CSCO','GOOG','NOK','AAPL','INTC','HPQ','ORCL','QCOM','CAJ','TSM','DELL','GLW','TXN','EMC','AMAT','DHR',
        'MOT',
        ]
    transportation = [
        'UPS','UNP','FDX','CNI','CSX','NSC','CP','EXPD','CHRW','LUV','RYAAY','LFL','JBHT','GSH','FRO','ZNH','TK','GOL','KSU',
        'AMR','UAUA','DAL',
        ]
    utilities = [
        'EXC','NGG','VE','SO','FPL','D','DUK','FE','ETR','PEG','KEP',
        'WMB','TRP','PPL','AEP','EIX',
        'CEG','SE','SRE',
        ]

    ## NYSE Europe
    Belgium = ['DEG']
    Greece = ['CCH','DAC','DSX','OTE','NBG','TNP']
    Guernsey = ['DOX']
    Hungary = []
    Luxembourg = ['TS','TX']
    Portugal = []
    # NYSE Asia/Pacific
    MarshallIslands = ['DHT','SSW','TGP','TK']
    NewZealand = ['NZT']
    Philippines = ['PHI']
    SouthKorea = ['KB','KEP','KTC','LPL','PKX','SKM','WF',]
    Taiwan = ['ASX','AUO','CHT','TSM','UMC','SPIL']



    d_indexes = {
        ## US indexes
        'DJI30':dji30,'Nasdaq100':nasdaq100,'S&P500':sp500,'Russell1000':russell1000,'Russell2000':russell2000,'Russell3000':russell3000,'DJU15':dju15,'DJT20':djt20,'DJC65':djc65,
        'S&P500 Dividend Aristocrat':sp500_aristocrats,
        'Future S&P500 Dividend Aristocrat?':sp500_aristocrats_future,
        'S&P Europe 350 Dividend Aristocrat':sp350_aristocrats,
        ## Asia indexes
        'Nikkei':nikkei225,'KOSPI':kospi,'SSE50':sse50,'TAIEX':taiex,'HSI45':hsi45,'STI30':sti30,'SENSEX30':sensex30,
        'SET50': set50,
        'KLCI':klci, 'JII30':jii30, 'KSE30':kse30,
        ## Europe
        'NasdaqEurope': europe,
        'FTSE100':ftse100, 'FTSE250':ftse250, 'LSE':LSE,
        'CAC40':cac40,'CACnext20':cacnext20,'DAX30':dax30,'MDAX50':mdax50,'IBEX35':ibex35,'MIB30':mib30,'ASX50':asx50,
        'TSX60':tsx60,'BEL20':bel20,'AEX25':aex25,
        'PSI20': psi20,
        'OMXS30':omxs30,'OMXC20':omxc20,'OBX25':obx25,'OMXN40':omxn40,'OMXH25':omxh25,'OMXI6':omxi6,'OMXB10':omxb10,
        'OMXsmall': omx_small,
        'OMXmedium': omx_medium,
        'ISEQ20':iseq20,'ATX20':atx20,
        'SMI20':smi20, 'Switzerland':switzerland,
        'WIG20':wig20,
        'ISE30':ise30,
        'MICEX10':micex10,
        'TA100': ta100,
        ## Middle East
        'TASI':tasi,
        ## South America
        'IPC35':ipc35,'Bovespa':bovespa,
        ## Africa
        'Johannesburg': jse,
        ## people
        'B2C': b2c,
        'Tommy':tommy, 'Buffett':Buffett, 'Super Brands':brands,
        ## Google Categories
        'basicmaterials':basicmaterials,'conglomerates':conglomerates,
        'capitalgoods':capitalgoods,'consumercyclical':consumercyclical,'consumernoncyclical':consumernoncyclical,
        'energy':energy,'financial':financial,'healthcare':healthcare,'services':services,'technology':technology,'transportation':transportation,'utilities':utilities,
        ## other
        'IXIC':IXIC,'NYSE':NYSE,
        'shipping':shipping,'railroad':railroad,
        ## Countries
        'Belgium':Belgium,'Greece':Greece,'Guernsey':Guernsey,'Hungary':Hungary,'Luxembourg':Luxembourg,'Portugal':Portugal,
        'Marshall Islands':MarshallIslands,'New Zealand':NewZealand,'Philippines':Philippines,'South Korea':SouthKorea,'Taiwan':Taiwan,
        ## All Nasdaq/NYSE companies
        'US': nasdaq_nyse_amex,
        }

    l_all = []
    for index in d_indexes:
        for ticker in d_indexes[index]:
            if ticker not in l_all:
                if '.OB' in ticker:
                    continue
                if ticker[-2:] in [
##                    '.L',
                    '-S','-C','-V'
                    ]:
                    ticker = ticker[:-2]
                if '-' in ticker and ticker[:3] != 'SE:':
                    ticker = ticker.replace('-','.')
##                    if '.' in ticker and not ('.A' in ticker or '.B' in ticker) and not ':' in ticker:
##                        print ticker
##                        stop
                if ticker not in l_all:
                    l_all.append(ticker)

    ## for avoiding duplicates, and for getting 10 year balance sheet data
    ## keys = US ticker
    ## values = MSN/Yahoo/Reuters tickers
    d_ADR = {
        ## North America
        ## exclude
        'CNI':['CA:CNR'],
        'TCK':['CA:TCKb'],
        'THI':['CA:THI'],
        'IMO':['CA:IMO'], ## imperial oil

        ## 10y
#        'FMX':['FEMSAUBD.MX'],
        'AMX':['AMXL.MX',],
        'GMBXF':['GMEXICOB.MX',],

        ##
        ## South America
        ##
        
        ## 10y
        'VALE':['VALE5.SA'],
        'TNE':['TNLP4.SA'],
        'CIG':['CMIG4.SA'],
        'ELP':['CPLE6.SA'],
        'ABVC':['AMBV4.SA'],
        'SBS':['SBSP3.SA'],
        ##'TMX':['TELMEXL.MX',] ## no 10y data

        ## Asia
        ## exclude
        'CAJ':['JP:7751'], ## Canon
        'TKPHF':['JP:4502'], ## Takeda
        'TM':['JP:7203'],
        'KAKKF':['JP:9107'], ## Kawasaki Kisen Kaisha Ltd
        'FANUF':['JP:6954'], ## Fanuc
        'JPSWF':['JP:5631'],
        'DCM':['JP:9437'], ## NTT DoCoMo

        ## 10y
        'INFY':['INFY.BO'],
        'CHL':['0941.HK'],
        'SNP':['0386.HK','600028.SS',],
        'PTR':['0857.HK','601857.SS',],
        'CEO':['0883.HK'],
        'PKX':['005490.KS'],
        'SKM':['017670.KS'],
        'NJ':['JP:6594'],
        'TLK':['TLKM.IJ'],
        'IIT':['ISAT.IJ'],
        'UNVR.IJ':['UNLRF'],
        '1913.HK':['PRDSF'],

        ## Europe
        ## exclude
        'GSK':['GB:GSK'],
        'AAUK':['GB:AAL'],
        'BRGYY':['GB:BG'],
        'ABB':['SE:ABB','ABBN.VX',],
        'AZN':['GB:AZN','SE:AZN','AZN.L','AZN.ST',],
        'ALFVF':['SE:ALFA'],
        'E':['IT:ENI'],
        'MT':['NL:MT'], ## arcelor mittal (luxembourg hq, paris listed)
        'UN':['NL:UNA'],
        'UL':['NL:UNA'],
##        'DASTY':['FR:DSY'],
        'FMS':['DE:FME',], ## Fresenius Medical Care

        ## 10y
##        'NVO':['NOVOb.CO'],
##        'RHHBY':['ROG.VX'],
##        'TKC':['TCELL.IS'],
##        'SYT':['SYNN.VX'],
##        'SCTBF':['SECUb.ST'],
####        'SCTBF':'SE:SECU-B',
##        'ROS':['RTKM.ME'],
##        'ABB':['ABBN.VX','SE:ABB',],
##        'EONGY':['DE:EOAN'],
##        'ACGY':['ACY.OL'],
##        'LUKOY':['LKOH.ME',],
##        'SCMWY':['SCMN.VX',],
##        'NVS':['NOVN.VX',], ## Novartis
##        'DSDVF':['DSDVF'],
##        'UHR.VX':['SWGAY'],

        ## South Africa
        ## 10y
        'IMPUY':['IMP.J'], ## Impala Platinum Holdings Ltd

        ## Australia
        ## exclude
        'BHP':['AU:BHP'], ## BHP Billiton
        'WOPEF':['AU:WPL'],
        'RIO':['AU:RIO','GB:RIO',],
        }

    ## replace ADR       
    for ticker in list(l_all):
##        continue  # tmp!!!
        if ticker in d_ADR.keys():
            print('remove ADR', ticker)
            l_all.remove(ticker)
            if d_ADR[ticker] not in l_all: ## use sets instead...
                l_all += d_ADR[ticker]
    l_all = list(set(l_all))

    l_all.sort()

    ## add ADR to indexes
    for ADR in d_ADR:
        ticker = d_ADR[ADR]
        for index in d_indexes.keys():
            if ADR in d_indexes[index]:
                if ticker not in d_indexes[index]:
                    d_indexes[index] += [ticker]

    return l_all, d_indexes, d_ADR

##        industries_majmin = {
##            ##
##            ## Reuters
##            ##
##            'Financial':[],
##            'Consumer':[],
####            ##
####            ## MSN
####            ##
####            'Aerospace and Defense': ['Aerospace/Defense - Major Diversified','Aerospace/Defense Products & Services',
####                                      'Diversified',],
####            'Automotive': ['Auto Manufacturers - Major','Auto Parts','Trucks & Other Vehicles','Recreational Vehicles',
####                           'Automotive',],
####            'Banking': [
####                'Foreign Money Center Banks','Foreign Regional Banks','Money Center Banks','Regional - Mid-Atlantic Banks','Regional - Midwest Banks','Regional - Northeast Banks','Regional - Pacific Banks','Regional - Southeast Banks','Savings & Loans', 'Regional - Southwest Banks',
######                'Financial',
####                ],
####            'Chemicals': ['Agricultural Chemicals','Chemicals - Major Diversified','Specialty Chemicals','Synthetics',
####                          'Chemicals',],
####            'Computer Software': ['Personal Products','Application Software','Business Software & Services','Technical & System Software','Multimedia & Graphics Software','Security Software & Services','Information & Delivery Services',
####                                  'Electronics',],
####            'Computer Hardware': ['Networking & Communication Devices', 'Data Storage Devices', 'Diversified Computer Systems', 'Computer Peripherals', 'Personal Computer Systems', 'Computer Based Systems'],
####            'Conglomerates': ['Conglomerates'],
####            'Consumer Durables': [
####                'Photographic Equipment & Supplies','Electronic Equipment','Toys & Games', 'Appliances', 'Sporting Goods', 'Recreational Goods, Other', 'Home Furnishings & Fixtures', 'Housewares & Accessories',
####                'Recreation', ## Konica-Minolta
####                ],
####            'Consumer Non-Durables': [
####                'Information Technology Services','Rubber & Plastics','Paper & Paper Products','Cleaning Products','Textile - Apparel Footware & Accessories', 'Personal Services', 'Packaging & Containers','Office Supplies',
####                'Apparel', ## DE.ADS
####                ],
####            'Diversified Services': ['Management Services','Staffing & Outsourcing Services','Education & Training Services','Business Services', 'Rental & Leasing Services','Consumer Services',],
####            'Drugs': ['Drug Manufacturers - Other','Drug Manufacturers - Major','Biotechnology','Drugs - Generic', 'Drug Delivery',
####                      'Drugs',],
####            'Electronics': ['Scientific & Technical Instruments','Semiconductor - Broad Line','Semiconductor - Integrated Circuits','Semiconductor - Specialized','Semiconductor Equipment & Materials','Printed Circuit Boards', 'Diversified Electronics', 'Semiconductor - Memory Chips'],
####            'Energy': ['Independent Oil & Gas','Major Integrated Oil & Gas','Oil & Gas Drilling & Exploration','Oil & Gas Equipment & Services','Oil & Gas Pipelines', 'Oil & Gas Refining & Marketing',
####                       'Oil, Gas, Coal & Related Services',],
####            'Financial Services': ['Asset Management','Investment Brokerage - National','Credit Services', 'Diversified Investments', 'Business Equipment', 'Investment Brokerage - Regional','Closed-End Fund - Foreign'],
####            'Food and Beverage': ['Beverages - Brewers','Food - Major Diversified','Beverages - Soft Drinks','Dairy Products', 'Confectioners', 'Beverages - Wineries & Distillers', 'Meat Products', 'Processed & Packaged Goods','Farm Products',],
####            'Health Services': ['Medical Instruments & Supplies','Health Care Plans','Medical Appliances & Equipment','Specialized Health Services','Medical Laboratories & Research', 'Hospitals', 'Home Health Care', 'Long-Term Care Facilities','Medical Practitioners',],
####            'Insurance': ['Accident & Health Insurance','Property & Casualty Insurance','Surety & Title Insurance','Life Insurance','Insurance Brokers'],
####            'Internet': ['Internet Software & Services', 'Internet Information Providers'],
####            'Leisure': [
####                'General Entertainment', 'Restaurants', 'Resorts & Casinos', 'Lodging', 'Special Eateries','Gaming Activities','Sporting Activities',
####                ],
####            'Manufacturing': ['Industrial Electrical Equipment','Small Tools & Accessories','Farm & Construction Machinery','Diversified Machinery', 'Industrial Equipment & Components','Metal Fabrication','Machine Tools & Accessories',
####                              'Metal Product Manufacturers','Machinery',],
####            'Materials and Construction': ['Cement','General Building Materials','Residential Construction', 'Heavy Construction', 'Lumber, Wood Production', 'Waste Management', 'General Contractors',
####                                           'Construction',],
####            'Media': ['Advertising Agencies', 'CATV Systems','Broadcasting - TV','Broadcasting - Radio','Entertainment - Diversified','Publishing - Newspapers', 'Marketing Services', 'Publishing - Periodicals', 'Publishing - Books',
####                      'Printing & Publishing',],
####            'Metals and Mining': ['Aluminium','Copper','Gold','Steel & Iron','Industrial Metals & Minerals', 'Nonmetallic Mineral Mining',
####                                  'Metal Production',],
####            'Real Estate': ['REIT - Residential','Property Management','REIT - Office','Mortgage Investment', 'REIT - Retail', 'Real Estate Development', 'REIT - Hotel/Motel', 'REIT - Industrial', 'REIT - Diversified',
####                            'Miscellaneous',],
####            'Retail': ['Grocery Stores','Apparel Stores','Home Furnishing Stores','Electronics Stores','Discount, Variety Stores','Drug Stores','Department Stores', 'Home Improvement Stores', 'Textile - Apparel Clothing', 'Catalog & Mail Order Houses', 'Auto Parts Stores',
####                       'Retailers','Textiles',],
####            'Specialty Retail': ['Specialty Retail, Other', 'Jewelry Stores', 'Auto Dealerships','Sporting Goods Stores'],
####            'Telecommunications': ['Telecom Services - Foreign','Communication Equipment','Diversified Communication Services','Wireless Communications','Telecom Services - Domestic','Processing Systems & Products'],
####            'Tobacco': ['Cigarettes', 'Tobacco Products, Other'],
####            'Transportation': ['Railroads','Air Delivery & Freight Services', 'Air Services, Other', 'Regional Airlines', 'Major Airlines','Shipping','Trucking',
####                               'Transportation',],
####            'Utilities': [
####                'Diversified Utilities','Electric Utilities','Foreign Utilities', 'Gas Utilities','Water Utilities',
####                'Utilities','Electrical',
####                ],
####            'Wholesale': ['Drugs Wholesale', 'Wholesale, Other','Computers Wholesale', 'Auto Parts Wholesale', 'Electronics Wholesale', 'Food Wholesale', 'Industrial Equipment Wholesale', 'Medical Equipment Wholesale','Basic Materials Wholesale',],
####            'N/A': ['N/A']
##            }
##        industries_minmaj = {}
##        for ind_maj in industries_majmin:
##            for ind_min in industries_majmin[ind_maj]:
##                industries_minmaj[ind_min] = ind_maj

if __name__=='__main__':
    main()
