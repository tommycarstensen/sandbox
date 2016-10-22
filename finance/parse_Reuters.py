## import executive payment from SEC filing instead?

import os, time
import screener

class Reuters:


    def parse_overview(self,ticker,url):

        d_factors = {'Mil.':1000000,}

        name = ''
        mc = ''
        currencyCode = 'USD'
        price = ''
        sector = ''
        industry = ''
        statementNA = False
        beta = ''

        print(url)
        lines = screener.finance().read_url(url, ticker)

        for i1 in range(len(lines)):

            if '<div id="sectionTitle">' in lines[i1]:
                index1 = lines[i1+1].index('<h1>')+4
                index2 = lines[i1+1].index('</h1>')
                s = lines[i1+1][index1:index2].strip()
                name = s

            ## price, currency
            if '<div class="sectionQuoteDetail">' in lines[i1]:
                for i2 in range(i1,len(lines)):
                    if '<span style="font-size: 23px;">' in lines[i2]:
                        index1 = 0
                        index2 = lines[i2+1].index('</span>')
                        s = lines[i2+1][index1:index2].replace(',','')
##                                print lines[i2+1]
                        if s in ('--', '\t\t\t\t--'):
                            statementNA = True
                        else:
                            price = float(s)
                            index1 = index2+lines[i2+1][index2:].index('<span>')+6
                            index2 = index1+lines[i2+1][index1:].index('</span>')
                            s = lines[i2+1][index1:index2]
                            currencyCode = s.upper() ## upper if GBp

                        break

            ## sector, industry
            if '<div id="sectionHeaderTopics"><div id="headerTopics">' in lines[i1]:
                index1 = lines[i1+5].index('/sectors')+8+1
                index2 = index1+lines[i1+5][index1:].index('"')
                sector = lines[i1+5][index1:index2]
                index1 = lines[i1+5].index('/sectors/industries/')+len('/sectors/industries/')
                index1 += lines[i1+5][index1:].index('>')+1
                index2 = index1+lines[i1+5][index1:].index('<')
                industry = lines[i1+5][index1:index2]

            ## beta
            if '<td>Beta:</td>' in lines[i1]:
                index1 = lines[i1+1].index('<strong>')+8
                index2 = lines[i1+1].index('</strong>')
                s = lines[i1+1][index1:index2]
                beta = s

            ## market cap
            if '<td>Market Cap' in lines[i1]:
                factor = 'Mil.'
                if not 'Mil.' in lines[i1]:
                    print(lines[i1])
                    stop
                index1 = lines[i1+1].index('<strong>')+8
                index2 = lines[i1+1].index('</strong>')
                s = lines[i1+1][index1:index2]

                s = s.replace('&#8361;','') ## KRW
                s = s.replace('&#8364;','') ## EUR
                s = s.replace('&#72;&#75;&#36;','') ## HKD
                s = s.replace('Â¥','') ## CNY
                s = s.replace('Â','') ## CNY
                s = s.replace('¥','') ## CNY
                s = s.replace('&#165;','') ## JPY
                s = s.replace('&#67;&#72;&#70;','') ## CHF
                s = s.replace('&#163;','') ## GBP
                s = s.replace('Rs','') ## INR
                s = s.replace('&#107;&#114;.','') ## DKK
                s = s.replace('&#107;&#114;','') ## NOK
                s = s.replace('TL','') ## Turkish Lira
                s = s.replace('&#82;','') ## Brazil
                s = s.replace('&#78;&#84;&#36;','') ## TWD
##                        s = s.replace('&#78;&#84;&#36;','') ## SGD
                s = s.replace('&#1088;&#1091;&#1073;','') ## Russia
                s = s.replace('&#76;&#116;','') ## Lithuania
                s = s.replace('&#36;','') ## USD (dollar symbol)
                s = s.replace('&#77;','') ## MYR (Malaysian ringgit - MR)
                s = s.replace('&#3647;','') ## THB
                s = s.replace('&#3647;','') ## IDR Indonesian Rupiah
                s = s.replace('&#8360;','') ## PKR Pakistani ...
                s = s.replace('&#8362;','')

                s = s.replace(',','')
                if s == '--':
                    print('mc', s)
                else:
                    mc = float(s)*d_factors[factor]

##                    ## sector
##                    if '<a href="/finance/industries/allIndustries">' in lines[i1]:
##                        index1 = lines[i1].index('<a href="/finance/industries/allIndustries">')+len('<a href="/finance/industries/allIndustries">')
##                        index2 = index1+lines[i1][index1:].index('</a>')
##                        sector = lines[i1][index1:index2].strip()
##
##                    ## industry
##                    if '<strong>industry:</strong>' in lines[i1]:
##                        index2 = lines[i1].rindex('<')
##                        index1 = lines[i1][:index2].rindex('>')+1
##                        industry = lines[i1][index1:index2]
##                        if industry == 'N/A':
##                            print('industry', industry)
##                            stop_temp

##                    if '<label>Mkt Cap.</label>' in lines[i1]:
##                        index1 = lines[i1+1].index('<span>')+6
##                        index2 = lines[i1+1].index('</span>')
##                        s = lines[i1+1][index1:index2]
##                        while ';' in s:
##                            index1 = s.index('&')
##                            index2 = s.index(';')
##                            if s[index2+1] == '.':
##                                s = s[:index1]+s[index2+2:]
##                            else:
##                                s = s[:index1]+s[index2+1:]
##                        s = s.replace(',','').replace('Â¥','').replace('Ã‚','')
##                        d_factors = {'M':1000000,}
##                        if s == '--M':
##                            statementNA = True
##                        elif s[-2:] == 'pM':
##                            mc = float(s[:-2])*d_factors[s[-1]]
##                        elif s[:2] == 'Rs':
##                            mc = float(s[2:-1])*d_factors[s[-1]]
##                        else:
##                            mc = float(s[:-1])*d_factors[s[-1]]

##                if mc == '' and statementNA == False:
##                    retry
##                else:
##                    break ## break loop of trys

        if name == '':
            statementNA = True

##        if price != '' and mc != '':
##            print(price, mc)
##            stop_loop

        if mc == '' or price == '' or mc == '--':
            statementNA = True

        return (
            name, currencyCode, price, sector, industry, statementNA,
            mc, beta,
            )
    

    def parse_statement(self,url,dic,statement,):

        print(url)

        d_periods = {
            '12 Weeks':'13 Weeks','14 Weeks':'13 Weeks',#'16 Weeks':'13 Weeks',
            '25 Weeks':'26 Weeks','27 Weeks':'26 Weeks',#'24 Weeks':'26 Weeks','28 Weeks':'26 Weeks','29 Weeks':'26 Weeks',
            '38 Weeks':'39 Weeks','40 Weeks':'39 Weeks',#'35 Weeks':'39 Weeks','36 Weeks':'39 Weeks',
            '51 Weeks':'52 Weeks','53 Weeks':'52 Weeks',#'48 Weeks':'52 Weeks',
##                                        '11 Months':'12 Months',
            }

        dic_out = {
            'period':[],
            'date':[],
            }

        lines = screener.finance().read_url(url, '0')

        statementNA = False
        bool_no_financials = False

        for i1 in range(len(lines)):

            if 'No Financials Data Available</div>' in lines[i1]:
                print('No Financials Data Available</div>')
                statementNA = True
                bool_no_financials = True
                stop5
                break

            if '<table class="dataTable financials" cellspacing="1" cellpadding="0" width="100%">' in lines[i1]:

                for i2 in range(i1+1,len(lines)):

                    if '<span class="units">' in lines[i2]:

                        l = lines[i2+1].split()
                        s = l[1]
                        d_factors = {'Millions':1000000.,'Thousands':1000.}
                        factor = d_factors[s.strip()]

                        index1 = str(lines[i2+1]).index(' of')+3
                        index2 = str(lines[i2+1]).index('<')
                        currency = str(lines[i2+1])[index1:index2].strip()

                    elif '<span class="period">' in lines[i2]:
                        index1 = 0
                        index2 = str(lines[i2-1]).index('<')
                        s = str(lines[i2-1])[index1:index2].strip()
                        if int(s[-2:]) <= 6: ## e.g. FXJ.AX, BBY, BKS
                            s = '%5s%02i%3s' %(s[:5],int(s[5:-3])-1,s[-3:],)
                        dic_out['date'] += [s[:-3]]

                        if statement != 'balance':
                            for i3 in range(i2+1,len(lines)):
                                if '</span>' in lines[i3]:
                                    index1 = 0
                                    index2 = lines[i3].index('</span>')
                                    s = lines[i3][index1:index2].strip().replace('&#160;',' ')
                                    if s in list(d_periods.keys()):
                                        s = d_periods[s]
                                    if s not in [
                                        '3 Months','6 Months','9 Months','12 Months',
                                        '13 Weeks','26 Weeks','39 Weeks','52 Weeks',
##                                        '27 Weeks','53 Weeks','25 Weeks',
                                        ]:
                                        if statement == 'income':
                                            statementNA = True
                                            break
                                        else:
                                            print(s)
                                            print((lines[i3]))
                                            print((lines[i3]))
                                            print((lines[i3][index1:index2]))
                                            stop
                                    l = s.split()
                                    dic_out['period'] += [[int(l[0]),l[1],]]
                                    break


                    elif '<tr ' in lines[i2]:

                        col1 = True
                        for i3 in range(i2+1,len(lines)):
                            if '<td ' in lines[i3]:
                                index1 = lines[i3].index('>')+1
                                index2 = lines[i3].rindex('<')
                                s = lines[i3][index1:index2].replace('  ',' ')
                                if col1 == True:
                                    key = s
                                    dic_out[key] = []
                                    col1 = False
                                else:
                                    s = s.replace(',','').replace('(','').replace(')','')
                                    if s == '--':
                                        s = 0
                                    value = factor*float(s)
                                    if 'minus' in lines[i3]:
                                        value *= -1
                                    dic_out[key] += [value]
                            if lines[i3].strip() == '</tr>':
                                break
                            if '<th>' in lines[i3].strip():
                                break

                break ## break loop over lines

        if lines == []:
            if statementNA == True:
                currency = 'N/A'
            else:
                stop_loop
            if bool_no_financials == False:
                stop

        return dic_out, statementNA, currency


    def parse_CEOcompensation(self, ticker):

        url = 'http://www.reuters.com/finance/stocks/companyOfficers?symbol=%s&viewId=comp' %(ticker)

        lines = screener.finance().read_url(url, ticker)

        compensation = 0
        i2 = 0
        l_urls = []
        for i1 in range(len(lines)):
            if i1 < i2:
                continue
            if (
                '<div class="moduleHeader"><h3>Basic Compensation</h3></div>' in lines[i1]
                or
                '<div class="moduleHeader"><h3>Options Compensation</h3></div>' in lines[i1]
                ):

                basic = False
                if '<div class="moduleHeader"><h3>Basic Compensation</h3></div>' in lines[i1]:
                    basic = True
                    i_add = 2

                options = False
                if '<div class="moduleHeader"><h3>Options Compensation</h3></div>' in lines[i1]:
                    options = True
                    i_add = 3

                row = 0
                for i2 in range(i1+1,len(lines)):
                    if '<tr' in lines[i2]:
                        row += 1
                        if row >= 2:
                            index2 = str(lines[i2+i_add]).index('</td>')
                            index1 = str(lines[i2+i_add])[:index2].index('>')+1
                            s = str(lines[i2+i_add])[index1:index2].replace(',','')
                            if s != '--':
                                compensation += float(s)
                            ## if basic compensation table
                            if i_add == 2:
                                index1 = lines[i2+i_add-1].index('<a href="')+len('<a href="')
                                index2 = index1+lines[i2+i_add-1][index1:].index('"')
                                s = lines[i2+i_add-1][index1:index2]
                                url_executive = 'http://www.reuters.com%s' %(s)
                                l_urls += [url_executive]
                    if '</table>' in lines[i2]:
                        break
                if options == True:
                    break

        ##
        ##
        ##
        for i_url_executive in range(len(l_urls)):
            url_executive = l_urls[i_url_executive]
            bool_break = False
            lines = screener.finance().read_url(url_executive, '0')

            for i in range(len(lines)):

                if 'Fiscal Year Total, ' in lines[i]:

                    ## currency
                    index1 = lines[i].index('Fiscal Year Total, ')+len('Fiscal Year Total, ')
                    index2 = index1+lines[i][index1:].index('<')
                    currency = lines[i][index1:index2]
                    break

            if currency != '':
                break
                

##            if (
##                'Chief Executive Officer' in lines[i]
##                or
##                '>President, Director<' in lines[i]
##                or
##                'Chief Exec Officer' in lines[i]
##                or
##                '>Chairman of the Board, President<' in lines[i] ## e.g. NATI
##                or
##                '>President, Representative Director<' in lines[i] ## e.g. 6902.T
##                or
##                '>Representative Executive President, Director<' in lines[i] ## e.g. 4902.T
##                or
##                '>Chairman of the Board, Representative Director<' in lines[i] ## e.g. 7205.T
##                or
##                '>Group Managing Director, Executive Director<' in lines[i] ## 0013.HK
##                or
##                '>Chairman of the Board, Managing Director<' in lines[i] ## 0012.HK
##                or
##                '>Chairman of the Board, Chairman of a Subsidiary, Representative Director<' in lines[i] ## e.g. 8035.T
##                or
##                '>General Manager<' in lines[i] ## e.g. TKC
##                or
##                '>Managing Director (CEO), Chairman of the Executive Committee, Director<' in lines[i] ## e.g. TOTF.PA
##                or
##                '>Managing Director, Executive Director<' in lines[i] ## 0006.HK
##                or
##                '>Deputy Chairman of the Board, Managing Director<' in lines[i] ## 0001.HK
##                or
##                '>Chairman of the Executive Committee, Director<' in lines[i] ## SOLB.BR
##                ):
##                index1 = lines[i-3].index('<a href="')+len('<a href="')
##                index2 = index1+lines[i-3][index1:].index('"')
##                url = 'http://www.reuters.com'+lines[i-3][index1:index2]
##                break
##
##        if i == len(lines)-1:
##
##            print url
##            if ticker_reuters not in ['NL:MT','YZC','FR:FP','0003.HK',]:
##                stop
##            compensation = 0
##
##        else:
##
##            for i in range(10):
##                try:
##                    urllines = urllib2.urlopen(url)
##                    lines = urllines.readlines()
##                    break
##                except:
##                    continue
##            if i == 9:
##                print url
##                stop
##
##
##            for i in range(len(lines)):
##
##                if 'Fiscal Year Total, ' in lines[i]:
##
##                    ## currency
##                    index1 = lines[i].index('Fiscal Year Total, ')+len('Fiscal Year Total, ')
##                    index2 = index1+lines[i][index1:].index('<')
##                    currency = lines[i][index1:index2]
##                    if currency == '':
##                        rate = 0.
##                    else:
##                        rate = d_currency[currency]
##
##                    ## compensation
##                    index1 = lines[i+6].index('>')+1
##                    index2 = lines[i+6].rindex('<')
##                    s = lines[i+6][index1:index2].replace(',','')
##                    if s == '--':
##                        compensation = 0.
##                    else:
##                        compensation = float(s)/rate
##
##                    break
##
##                if i == len(lines)-1:
##                    compensation = 0.
##                    print url

        return compensation, currency

if __name__=='__main__':
    instance_reuters = Reuters()

##    t = instance_reuters.parse_CEOcompensation('AZN')
##    print(t)
##    stop

    t = instance_reuters.parse_overview('IBM')
    print(t)
    stop_success

    d_msn2yahoo = {'GB':'L','FR':'PA','DE':'DE','AU':'AX','ES':'MC','JP':'T','IT':'MI','SE':'ST','BE':'BR','CA':'TO','NL':'AS',}
    d_currency = {}
    ticker_reuters = '0016.HK'

    url = 'http://in.reuters.com/money/quotes/quote?symbol=XOM'
    dic = {}
    statement = 'Balance'
    t = instance_reuters.parse_quote('XOM',)
    print(t)

