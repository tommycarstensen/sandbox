import import_tickers

tickers, d_indexes, d_ADR = import_tickers.main()

d_yahoo2msn = {
    'GB:BG': 'GB:BG.',
    'GB:RB': 'GB:RB.',
    'GB:SN': 'GB:SN.',
    'FR:EF1': 'FR:EI',
##            'BFb': 'BF.B',
    }

## also bw2reuters
## also ms2retuers
## retuers just different...
d_yahoo2reuters = {
    ## Austria
    'MMK.VI': 'MMKV.VI',
    'EBS.VI': 'ERST.VI',
    'EVN.VI': 'EVNV.VI',
    'RHI.VI': 'RHIV.VI',
    'OMV.VI': 'OMVV.VI',
    ## Belgium
    'BELG.BR': 'BCOM.BR',
    'OME.BR': 'OMEP.BR',
    ## Denmark
    'NOVO-B.CO': 'NOVOb.CO',
    ## England
    'AU.L': 'AUTN.L', 'BA.L': 'BAq.L', 'BBY.L': 'BALF.L', 'BT.A.L': 'BT.L', 'SDRC.L': 'SDR.L',
    'EVR.L': 'HK1q.L',
    ## France
    'EF1.PA': 'ESSI.PA', 'VK.PA': 'VLLP.PA', 'DSY.PA': 'DAST.PA', 'ALO.PA': 'ALSO.PA', 'SU.PA': 'SCHN.PA', 'FP.PA': 'TOTF.PA', 'SAN.PA': 'SASY.PA',
    'BN.PA': 'DANO.PA', 'ALU.PA': 'ALUA.PA',
    'DG.PA': 'SGEF.PA',
    ## Greece
    'OPAP.AT': 'OPAR.AT',
    ## Netherlands
    'MT.AS': 'ISPA.AS',
    'DSM.AS': 'DSMN.AS', 'UNA.AS': 'UNc.AS', 'AH.AS': 'AHLN.AS',
    'HEIA.AS': 'HEIN.AS',
    ## Germany
    'HEN3.DE': 'HNKG_p.DE', 'BAS.DE': 'BASF.DE', 'ADS.DE': 'ADSG.DE',
    'TKA.DE': 'TKAG.DE', 'LHA.DE': 'LHAG.DE', 'SDF.DE': 'SDFG.DE',
    'BEI.DE': 'BEIG.DE', 'BAYN.DE': 'BAYG.DE', 'MAN.DE': 'MANG.DE',
    'SAP.DE': 'SAPG.DE', 'ALV.DE': 'ALVG.DE', 'BMW.DE': 'BMWG.DE',
    'CBK.DE': 'CBKG.DE', 'RWE.DE': 'RWEG.DE', 'DAI.DE': 'DAIGn.DE',
    'DB1.DE': 'DB1Gn.DE', 'DBK.DE': 'DBKGn.DE', 'DPW.DE': 'DPWGn.DE',
    'DTE.DE': 'DTEGn.DE', 'EOAN.DE': 'EONGn.DE', 'IFX.DE': 'IFXGn.DE',
    'LIN.DE': 'LING.DE', 'MRC.DE': 'MRCG.DE', 'SIE.DE': 'SIEGn.DE',
    'SZG.DE': 'SZGG.DE', 'VOW.DE': 'VOWG.DE',
    'FME.DE': 'FMEG.DE',
    'FIE.DE': 'FIEG.DE',
    ## Iceland
    'ATLA.IC': 'FOATLA.IC',
    'FOAIR.IC': 'FO-AIR.IC',
    ## Ireland
    'PLS.IR': 'PAP.I', 'DRS.IR': 'DGO.I',
    'PWL.IR': 'PAP.I',
    'KRX.IR': 'KYGa.I',
    'GN5.IR': 'GRF_u.I',
    ## Italy
    'SPM.MI': 'SPMI.MI',
    'TME.MI': 'TLIT.MI',
    ## Norway
    'SCH.OL': 'SBST.OL',
    ## Spain
    'SAB.MC': 'SABE.MC', 'FER.MC': 'FER1.MC',
    'CRI.MC': 'CRIT.MC',
    'EBRO.MC': 'EVA.MC',
    'ENG.MC': 'ENAG.MC',
    'GRF.MC': 'GRLS.MC',
    'SYV.MC': 'SVO.MC', ## SACYR VALLEHERMOSO
    ## Sweden
    'SKF-B.ST': 'SKFB.ST',
    ## Mexico
    'FEMSAUBD.MX': 'FMSAUBD.MX',
    'PE&OLES*.MX': 'PENOLES.MX',
    ## Kuala Lumpur
    '2089.KL': 'UTPS.KL',
    ## Singapore
    'SCI.SI': 'SCIL.SI',
    'C07.SI': 'JCYC.SI',
    'JCNC.SI': 'JCYC.SI',
    'COS.SI': 'COSC.SI',
    'CIT.SI': 'CTDM.SI',
    'SIE.SI': 'SIAE.SI',
    'DBS.SI': 'DBSM.SI',
    'JM.SI': 'JARD.SI',
    ## India
    'BHARTIARTL.BO': 'BRTI.BO',
    'BHARTI.BO': 'BRTI.BO',
    'RCOM.BO': 'RLCM.BO',
    'HNDL.BO': 'HALC.BO',
    'MSIL.BO': 'MRTI.BO',
    'GRASIM.BO': 'GRAS.BO',
    'HUVR.BO': 'HLL.BO',
    'SCS.BO': 'SATY.BO',
    'SUNP.BO': 'SUN.BO',
    ## South Africa
    'IMP.J': 'IMPJ.J',
    }

d_bw2yahoo = {
    'INFO.BO': 'INFY.BO',
    'MSIL.BO': 'MARUTI.BO',
    'BHARTI.BO': 'BRTI.BO',
    'JCNC.SI': 'C07.SI',
    'SJR/B.TO': 'SJRb.TO',
    'WPRO.BO': 'WIPRO.BO',
    'CA:SJR/B': 'SJRb.TO', ## ???
    'UPL.KL': '2089.KL',
    'PAP.IR': 'PLS.IR', 'PWL.IR': 'PLS.IR',
    }

d_bw2reuters = {
    'BHARTI.BO': 'BRTI.BO',
    'EI.PA': 'ESSI.PA',
    'NESN.VX': 'NESN.S',
    }

d_msn2yahoo = {
    'GB': 'L', 'FR': 'PA', 'DE': 'DE', 'AU': 'AX', 'ES': 'MC', 'JP': 'T', 'IT': 'MI', 'SE': 'ST', 'BE': 'BR', 'CA': 'TO', 'NL': 'AS',
    }

d_yahoo2bw = {
    'VX': 'S',  ## Switzerland
#    'VX': 'VX', ## Switzerland
    'PA': 'PA', ## Pakistan
    'AB': 'AB', ## Saudi Arabia
    'SS': 'CH', ## China (before Stockholm)
    'BO': 'IN', ## India
    'VL': 'LH', ## Lithuania / Vilnius
    'AX': 'AU', ## Australia
    'SA': 'BZ', ## Brazil / Sao Paolo
    'MX': 'MM', ## Mexico
    'SI': 'SP', ## Singapore
    'TO': 'CN', ## Canada
    'HK': 'HK', ## Hong Kong
    'DE': 'GR', ## Germany
    'L': 'LN', ## London
    'T': 'JP', ## Japan / Tokyo

    'QD': 'QD', ## Qatar / Doha
    'IJ': 'IJ', ## Indonesia / Jakarta
    'KL': 'MK', ## Malaysia / Kuala Lumpur
    'KS': 'KS', ## (South) Korea / Seoul
    'AT': 'GA', ## Greece / Athens

    'ME': 'RM', ## Russia / Moscow

    'TLV': 'TLV', ## Israel / Tel Aviv  ## actually ft.com suffix...
    'SET': 'SET', ## actually ft.com suffix...
    'JNB': 'JNB', ## actually ft.com suffix...

    'TL': 'ET', ## Estonia / Tallinn
    'VI': 'AV', ## Austria / Vienna
    'IC': 'IR', ## Iceland / Reykjavik
    'TW': 'TT', ## Taiwan / Taipei
    'MI': 'IM', ## Italy / Milano
    'HE': 'FH', ## Finland / Helsinki
    'PA': 'FP', ## France / Paris
    'ST': 'SS', ## Sweden / Stockholm
    'IR': 'ID', ## Ireland / Dublin
    'OL': 'NO', ## Norway / Oslo
    'IS': 'TI', ## Turkey / Istanbul
    'MC': 'SM', ## Spain / Madrid
    'AS': 'NA', ## Netherlands / Amsterdam
    'CO': 'DC', ## Denmark / Copenhagen
    'BR': 'BB', ## Belgium / Brussels
    'J': 'SJ', ## South Africa / Johannesburg

    'RG': 'RG', ## Latvia / Riga

    ## Poland not on Yahoo...
    'PW': 'PW', ## Poland / Warzaw

    ## Reuters...
    'I': 'ID', ## Ireland / Dublin
    }

## http://www.iso15022.org/MIC/ISO10383_MIC.pdf
d_msnprefix2morningstar = {
    'DE': 'XETR',
    'GB': 'XLON',
    'FR': 'XPAR',
    'ES': 'XMCE',  # MERCATO CONTINUO ESPANOL
    'CA': 'XTSE',  # Toronto Stock Exchange
    'JP': 'XTKS',  # Tokyo Stock Exchange
    'BE': 'XBRU',  # Brussels, Belgium
    'SE': 'XSTO',  # Stockholm, Sweden
    'AU': 'XASX',  # Australia
    }
d_yahoosuffix2morningstar = {
    'HK': 'XHKG',  # Hong Kong
#    'HK': 'SEHK',  # Hong Kong
    'VX': 'XSWX',  # Switzerland
    'SI': 'XSES',  # Singapore
    'SA': 'XBSP',  # Sao Paolo, Brazil, Bovespa...
    'BO': 'XNSE',  # National Stock Exchange of India
    'IS': 'XIST',  # Istanbul
    'L': 'XLON',  # London
    'TO': 'XTSE',  # Toronto
    'TW': 'XTAI',  # Taipei, Taiwan
    'OL': 'XOSL',  # Oslo
    'HE': 'XHEL',  # Helsinki
    'TL': 'XTAL',  # Tallinn, Estonia
    'IC': 'XICE',  # Iceland
    'CO': 'XCSE',  # Copenhagen Stock Exchange
    'IR': 'XDUB',  # Dublin, Ireland
    'MX': 'XMEX',  # Mexico
    'VL': 'XLIT',  # Vilnius, Lithuania
##    'HK': 'SEHKXHKG',  # Hong Kong
    'LS': 'XLIS',  # Lisbon, Portugal
    'VI': 'WBAG',  # Vienna, Austria
    }


def unknown2yahoo(ticker):

    ticker = msn2yahoo(ticker)

    ticker = reuters2yahoo(ticker)

    ticker = bw2yahoo(ticker)

    ticker = ticker.replace('-SEK', '')
    ticker = ticker.replace('-DKK', '')

    if ticker.endswith('.SS'):
        ticker = ticker[:-3]+'.SZ'
        
    return ticker


def msn2morningstar(ticker):

    # http://financials.morningstar.com/financials/getFinancePart.html?&callback=?&t=COLO%20B&region=dnk

    if '.CO' in ticker:
        ticker = ticker[:-3].replace('-', '').replace('b', '%20b')
    if ticker.startswith('SE:'):
        ticker = ticker.replace('a', '%20a')  # e.g. SE:ATCOa -> SE:ATCO%20a

    if ': ' in ticker:
        index = ticker.index(': ')
        try:
            ticker = d_msnprefix2morningstar[ticker[:index]]+ticker[index:]
        except:
            pass
    elif '.' in ticker and not ticker.endswith('.A') and not ticker.endswith('.B'):
        index = ticker.rindex('.')
        if ticker.endswith('.HK'):
            ticker = d_yahoosuffix2morningstar[ticker[index+1:]]+":0"+ticker[:index]
        else:
            ticker = d_yahoosuffix2morningstar[ticker[index+1:]]+":"+ticker[:index]

    if ticker.endswith('.L'):
        ticker = 'XLON:'+ticker[:-2]
        stop1

    if ticker.endswith('.MX'):
        ticker = 'XMEX:'+ticker[:-3]
        if ticker.endswith('.B'):
            ticker = ticker[:-1]+'%s20B'  # http://financials.morningstar.com/financials/getFinancePart.html?&callback=?&t=XMEX:GAP%20B

    return ticker


def bw2yahoo(ticker):

    if ticker in list(d_bw2yahoo.keys()):
        ticker = d_bw2yahoo[ticker]

    return ticker


def reuters2yahoo(ticker_reuters):

    ticker_yahoo = ticker_reuters
    for ticker in list(d_yahoo2reuters.keys()):
        if d_yahoo2reuters[ticker] == ticker_reuters:
            ticker_yahoo = ticker
            break

    return ticker_yahoo


def unknown2reuters(ticker):

    ## businessweek
    ticker = ticker.replace('/A', 'a')
    ticker = ticker.replace('/B', 'b')
    ticker = ticker.replace('*', '') ## e.g. ELEKTRA*.MX

    ticker = msn2reuters(ticker)

    ticker = yahoo2reuters(ticker)

    if ticker in list(d_bw2yahoo.keys()):
        ticker = d_bw2yahoo[ticker]

    if ticker in list(d_bw2reuters.keys()):
        ticker = d_bw2reuters[ticker]

    if ticker.endswith('.ME'):
        ticker.replace('.ME', '.MEP')

    return ticker


def unknown2businessweek(ticker):

    if ticker[-2:] in ['.A', '.B',]:
        ticker = ticker[:-2]+ticker[-2:].replace('.', '/') ## e.g. BT/A:LN
    if ticker[-1] in ['a', 'b',] and 'SE:' not in ticker:
        ticker = ticker[:-1]+'/'+ticker[-1:] ## e.g. JW/A, TCK/A:CN (not TEL2b:SS)

    ## Reuters
    if ticker[-2:] in ['.O', '.N',]:
        ticker = ticker[:-2]

    if ticker[-1] == '.': ## e.g. RB/:LN
        ticker = ticker[:-1]+'/'

    ticker = msn2yahoo(ticker)

    if '.HK' in ticker:
        while ticker[0] == '0':
            ticker = ticker[1:]
    if '-' in ticker:
        ticker = ticker.replace('-', '') ## e.g. SKF-B.ST > SKFB:SS (not BT/A:LN)
    if '.' in ticker:
        ticker = ticker.replace('.', ': ') ## e.g. 0386.HK > 386:HK

    if ': ' in ticker:
        index = ticker.rindex(': ')
        print(ticker)
        ticker = ticker[:index]+': '+d_yahoo2bw[ticker[index+1:]]

##    ticker = ticker.replace(':O', '') ## Reuters can't differentiate if not .O/.N suffix (after Oslo)

    return ticker


def msn2reuters(ticker):

    if ticker[-1] == '.':
        ticker = ticker[:-1] ## e.g. GB:RB. > RB.L

    ## ticker, msn2reuters
    if ': ' in ticker:
        index = ticker.index(': ')
        ticker_reuters = ticker[index+1:]+'.'+d_msn2yahoo[ticker[:index]]
    else:
        if ticker in d_indexes['IXIC'] and '.' not in ticker:
            ticker_reuters = ticker+'.O'
        elif ticker in d_indexes['NYSE'] and '.' not in ticker:
            ticker_reuters = ticker+'.N'
        else:
            ticker_reuters = ticker
    ticker_reuters = ticker_reuters.replace('.a', 'a')
    ticker_reuters = ticker_reuters.replace('.b', 'b')

    return ticker_reuters


def yahoo2reuters(ticker):

    ticker_reuters = ticker
    if ticker_reuters in list(d_yahoo2reuters.keys()):
        ticker_reuters = d_yahoo2reuters[ticker_reuters]
    if '.IR' in ticker_reuters:
        ticker_reuters = ticker_reuters[:-1]

    return ticker_reuters


def msn2yahoo(ticker):

    if ': ' in ticker:
        index = ticker.index(': ')
        ticker = ticker[index+1:]+'.'+d_msn2yahoo[ticker[:index]]
    ## Coloplast and Novo Nordisk and DSV
    if 'b.CO' in ticker:
        ticker = ticker.replace('b.CO', '-B.CO')

    return ticker


def yahoo2msn(ticker):

    if ticker in list(d_yahoo2msn.keys()):
        ticker_msn = d_yahoo2msn[ticker]
    else:
        ticker_msn = ticker
        if (
            'SE:' in ticker
            or
            '.ST' in ticker
            ):
            ticker_msn = ticker_msn.replace('a', '-A')
            ticker_msn = ticker_msn.replace('b', '-B') ## SE:TEL2-B
        else:
            ticker_msn = ticker_msn.replace('a', '.A') ## BFa
            ticker_msn = ticker_msn.replace('b', '.B') ## HUBb,BFb

    return ticker_msn


def unknown2ADR(ticker):

    return ADR


def unknown2advfn(ticker):

    ticker = ticker.replace('GB:', 'LSE:')  # London
    ticker = ticker.replace('DE:', 'FWB:')  # Frankfurt
    ticker = ticker.replace('JP:', 'TSE:')  # Tokyo SE
    if ticker.endswith('.CO'):
        ticker = 'PLUS:'+ticker[:-3]
    if ticker.endswith('.L'):
        ticker = 'LSE:'+ticker[:-2]

    return ticker


def unknown2ft(ticker):

    if ticker.startswith('SE:') and ticker.endswith('a'):
        ticker = ticker[:-1]+'+A'  # Stockholm, Sweden
    if ticker.startswith('SE:') and ticker.endswith('b'):
        ticker = ticker[:-1]+'+B'  # Stockholm, Sweden
    if ticker.endswith('a.CO'):
        ticker = ticker[:-4]+'+A.CO'
    if ticker.endswith('b.CO'):
        ticker = ticker[:-4]+'+B.CO'

    if ticker == 'SE:NDA':
        ticker = 'SE:NDA+SEK'
    if ticker == 'NDA.CO':
        ticker = 'NDA+DKK.CO'

    if ticker.endswith('-B.TO'):
        ticker = ticker[:-5]+'.B:TOR'  # e.g. BBD-B -> BBD.B
    if ticker.endswith('b') and ticker.startswith('CA:'):
        ticker = ticker[3:-1]+'.B:TOR'  # e.g. CA:RCIb -> RCI.B:TOR
    if ticker.endswith('.B'):
        if ticker.startswith('CA:'):
            pass  # e.g. CA:ATD.B -> CA:ATD.B
        else:
##            ticker = ticker[:-2]+'B'  # e.g. VIA.B -> VIAB
            ticker = ticker  # e.g. BF.B -> BF.B

    ticker = ticker.replace('-', '+')

    if ticker.startswith('DE:'):
        if ticker[3:] in ('HEN3', 'VOW3'):
            ticker = ticker[3:]+':GER'
        elif not ticker[3:] in (
            'ARL', 'CBK', 'BEI', 'BMW', 'BYW6', 'CON', 'D9C',
            'ALV', 'AIR', 'FIE', 'FME', 'FPE', 'FRA', 'FRE', 'G1A', 'GBF',
            'GFJ', 'GIL', 'GXI', 'HDD', 'HEI', 'HHFA', 'HOT', 'IVG', 'HLA',
            'KRN', 'LHA', 'LIN', 'LXS', 'MAN', 'MEO', 'MRK', 'NDA', 'PRA',
            'PUM', 'RAA', 'RHK', 'RHM', 'RWE', 'SAP', 'SGL', 'SKYD', 'SY',
            'SZG', 'SZU', 'TKA', 'VOS', 'WCH', 'WIN',
            ):
            ticker = ticker[3:]+'X.N:GER'  # Xetra, Germany
        else:
            ticker = ticker[3:]+'X:GER'
    if ticker.startswith('GB:'):
        ticker = ticker[3:]+':LSE'  # London Stock Exchange
    if ticker.endswith('.L'):
        ticker = ticker[:-2]+':LSE'  # London Stock Exchange
    if ticker.endswith('.CO'):
        ticker = ticker[:-3]+':CPH'  # Copenhagen
    if ticker.endswith('.SA'):
        ticker = ticker[:-3]+':SAO'  # Sao Paolo
    if ticker.endswith('.HK'):
        ticker = ticker[:-3]+':HKG'  # Hong Kong
    if ticker.endswith('.VX'):
        ticker = ticker[:-3]+':SWX'  # Swiss Exchange, Switzerland
    if ticker.endswith('.VX'):
        ticker = ticker[:-3]+':VTX'  # VIRT-X, Switzerland
    if ticker.startswith('JP:'):
        ticker = ticker[3:]+':TYO'  # Tokyo
    if ticker.endswith('.SS'):
        ticker = ticker[:-3]+':SHH'
    if ticker.startswith('ES:'):
        ticker = ticker[3:]+':MCE'  # Mercados Espanoles?
    if ticker.endswith('.OL'):
        ticker = ticker[:-3]+':OSL'  # Oslo
    if ticker.startswith('CA:'):
        ticker = ticker[3:]+':TOR'  # Toronto, Canada
    if ticker.startswith('FR:'):
        ticker = ticker[3:]+':PAR'  # Paris
    if ticker.endswith('.PA'):
        ticker = ticker[:-3]+':KAR'  # Karachi, Pakistan
    if ticker.endswith('.VI'):
        ticker = ticker[:-3]+':VIE'  # Vienna, Austria
    if ticker.endswith('.TLV'):
        ticker = ticker[:-4]+':TLV'  # Tel Aviv, Israel
##    if ticker.endswith('.KS'):
##        ticker = ticker[:-4]+':xxx'  # 
    if ticker.endswith('.SET'):
        ticker = ticker[:-4]+':SET'  # 
    if ticker.startswith('SE:'):
        ticker = ticker[3:]+':STO'  # Stockholm, Sweden
    if ticker.endswith('.ST'):
        ticker = ticker[:-3]+':STO'  # Stockholm, Sweden
    if ticker.endswith('.JNB'):
        ticker = ticker[:-4]+':JNB'  # Johannesburg, South Africa
    if ticker.endswith('.IC'):
        ticker = ticker[:-3]+':ICX'  # Iceland
    if ticker.endswith('.HE'):
        ticker = ticker[:-3]+':HEX'  # Helsinki, Finland
    if ticker.endswith('.SI'):
        ticker = ticker[:-3]+':SES'  # Singapore
    if ticker.endswith('.TL'):
        ticker = ticker[:-3]+':TLX'  # Tallinn, Estonia
    if ticker.endswith('.IS'):
        ticker = ticker[:-3]+':IST'  # Istanbul, Turkey
    if ticker.startswith('AU:'):
        ticker = ticker[3:]+':ASX'  # Australia
    if ticker.startswith('BE:'):
        ticker = ticker[3:]+':BRU'  # Brussels, Belgium
    if ticker.endswith('.IR'):
        ticker = ticker[:-3]+':ISE'  # Dublin, Ireland Stock Exchange
    if ticker.endswith('.MX'):
        ticker = ticker[:-3]+':MEX'  # Mexico
    if ticker.endswith('.ME'):
        ticker = ticker[:-3]+':MCX'  # Russia / Moscov Exchange
    if ticker.endswith('.PW'):
        ticker = ticker[:-3]+':WSE'  # Poland / Warzawa Stock Exchange
    if ticker.endswith('.BO'):
        ticker = ticker[:-3]+':NSI'  # Bombay / National Stock Exchange of India
    if ticker.endswith('.VL'):
        ticker = ticker[:-3]+':VLX'  # Vilnius, Lithuania
    if ticker.endswith('.KL'):
        ticker = ticker[:-3]+':KLS'  # Malaysia, Kuala Lumpur Stock Exchange
    if ticker.endswith('.IJ'):
        ticker = ticker[:-3]+':JKT'  # Jakarta, Indonesia
    if ticker.endswith('.J'):
        ticker = ticker[:-2]+':JNB'  # Johannesburg, South Africa
    if ticker.endswith('.NL'):
        ticker = ticker[:-3]+':AEX'  # Amsterdam, Netherlands
    if ticker.startswith('NL:'):
        ticker = ticker[3:]+':AEX'  # Amsterdam, Netherlands
    if ticker.startswith('IT:'):
        ticker = ticker[3:]+':MIL'  # Milano, Italy
    if ticker.endswith('.TW'):
        ticker = ticker[:-3]+':TAI'  # Taipei, Taiwan
    if ticker.endswith('.LS'):
        ticker = ticker[:-3]+':LIS'  # Lissabon, Portugal

    while ticker.startswith('0'):
        ticker = ticker[1:]

    return ticker


d_msn2currency = {
    'GB': 'British Pounds',
    'AU': 'Australian Dollars',
    'ES': 'Euro', 'IT': 'Euro', 'BE': 'Euro', 'NL': 'Euro', 'FR': 'Euro', 'DE': 'Euro',
    'SE': 'Swedish Krona',
    'JP': 'Yen',
    'CA': 'Canadian Dollars',
    }

if __name__ == '__main__':
    print(unknown2yahoo('PWL.IR'))
