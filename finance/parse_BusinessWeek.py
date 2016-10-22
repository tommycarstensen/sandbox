## import executive payment from SEC filing instead?

import os
import time
import screener

class BusinessWeek:


    def parse_statement(self,url,):

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
        currency = None

        lines = screener.finance().read_url(url, '0')

        statementNA = False
        statement_pending = False

        for i1 in range(len(lines)):

            if '>FINANCIALS SECTOR<' in lines[i1]:
                statementNA = True
                break

            if 'There are no Income Statements available at this time for' in lines[i1]:
                print('There are no Income Statements available at this time')
                statement_pending = True
                break

            if ' public company results.' in lines[i1]:
                print('Multiple or zero results')
                print(lines[i1])
                statementNA = True
                break
            
##                    if '<div id="resultCaption">Your search for <strong>' in lines[i1]:
##                        statementNA = True
##                        break

            ## Empty Table (nothing between end of table headers and end of table)
            if '4-Year<br />Trend</td></tr></table>' in lines[i1]:
                print('Blank Statement')
                statementNA = True
                break

##                    if not '<div style="padding-bottom:12px;">' in lines[i1]:
##                        continue
            if not '<div class="financialsSelectContainer">' in lines[i1]:
                continue

            s = lines[i1]

            index1 = s.index('Currency in<br />')+len('Currency in<br />')
            index2 = index1+s[index1:].index('</td>')
            index = index1+s[index1:index2].index(' ')+1
            d_factors = {'Millions':1000000.,'Thousands':1000.}
            factor = d_factors[s[index1:index-1]]
            index = index+s[index:index2].index(' ')+1
            currency = s[index:index2]

            index2 = index1+s[index1:].index('</tr>')
            index1 += s[index1:index2].index('As of:')
            s_year = s[index1:index2]

##                    if (
##                        'Press<br />Release' in s_year
##                        and
##                        ## often Japanese financial statements are press releases in BusinessWeek for a long long time it seems...
##                        url[-3:] != ':JP'
##                        ):
##                        print 'Press Release'
##                        statement_pending = True
##                        break

            index = 0
            while '<br />' in s_year[index:]:
                index += s_year[index:].index('<br />')+1
                ## company has been around for less than 4 years
                if s_year[index:index+10] == 'br />--</t':
                    print('Less than 4 years of data')
                    statement_pending = True
                    break
                if s_year[index:index+10] not in [
                    'br /><span',
                    'br />Trend', ## 4-Year Trend
                    'br />Relea', ## Press Release
                    'br />--<br',
                    ]:
                    year = int(s_year[index+5:index+9])
                    dic_out['date'] += [year]

            td_label = '<td class="statementLabel cell '
            td_label = '<tr>'
            index = 0
##                    while td_label1 in s[index:] or td_label2 in s[index:]:
            while td_label in s[index:]:
                index += s[index:].index(td_label)+1
                index2 = index+s[index:].index('</td>')
                index1 = index+s[index:index2].rindex('>')+1
                key = s[index1:index2]
                if ' ' in key:
                    if key[:key.index(' ')] in list(d_factors.keys()):
                        index += 1
                        continue
                dic_out[key] = []
                s_row = s[index2:index2+s[index2:].index('</tr>')]
                index_row = 0
##                        p = r'<span class="quoteData">$78.08 <span'
                while '<td' in s_row[index_row:]:
                    index_row += s_row[index_row:].index('<td')+1
                    index_row1 = index_row+len('<td')-1
                    index_row1 += s_row[index_row1:].index('>')+1
                    index_row2 = index_row1+s_row[index_row1:].index('</td>')
                    s_value = s_row[index_row1:index_row2]
                    if s_value[:10] == '<img src="':
                        break
                    s_value = s_value.replace(',','')
                    if s_value == '--':
                        s_value = 0
                    if s_value == '&nbsp;':
                        continue
##                            print(s_row)
                    value = factor*float(s_value)
                    dic_out[key] += [value]

        if len(list(dic_out.keys())) == 2:
            statementNA = True

        if currency == None and statementNA == False:
            statementNA = True
##            stop_ticker_wrong_or_not_existing

        return dic_out, statementNA, statement_pending, currency


    def parse_overview(self,ticker,):

        d_factors = {
            'K':1000,
            'M':1000000,
            'B':1000000000,
            'T':1000000000000,
            }

        name = ''
        mc = ''
        currencyCode = ''
        price = ''
        sector = ''
        industry = ''
        statementNA = False

        url = 'http://investing.businessweek.com/research/stocks/snapshot/snapshot.asp?ticker=%s' %(ticker)

        lines = screener.finance().read_url(url, ticker)

        for i1 in range(len(lines)):

            ## name
##                    if '<h1 id="companyTitle"' in lines[i1]:
##            if '<h2 class="pageHeader">' in lines[i1]:
            if '<span itemprop="name">' in lines[i1]:
##                index1 = lines[i1].index('<h2 class="pageHeader">')+len('<h2 class="pageHeader">')
                index1 = lines[i1].index('<span itemprop="name">')+len('<span itemprop="name">')
                index2 = index1+lines[i1][index1:].index('<')
                print(111,lines[i1])
                s = lines[i1][index1:index2].strip()
                print(222,s)
                if '(' in s:
                    s = s[:s.rindex('(')]
                name = s.upper()

            ## price, currency
            if (
                '<div class="dataPoint"><span class="quoteHeading">LAST</span> <span class="quoteData">' in lines[i1]
                or
                '<div class="dataPoint"><span class="quoteHeading">Last</span> <span class="quoteData">' in lines[i1]
                ):
                index1 = lines[i1].index('<span class="quoteData">')+len('<span class="quoteData">')
                index2 = index1+lines[i1][index1:].index('<')
                s = lines[i1][index1:index2].strip()
                for x in 'abcdefghijklmnopqrstuvwxyz&;$':
                    s = s.replace(x,'')
                    s = s.replace(x.upper(),'')
                s = s.replace(',','')
                if s[0] == '.':
                    s = s[1:] ## e.g. SFr.
                if s == '--':
                    statementNA = True
                    break
                price = float(s)
                if '<span class="xSmGreyTxt">' in lines[i1]:
                    index2 = index1+lines[i1][index1:].index('</span>')
                    index1 = lines[i1][:index2].rindex('>')+1
                    s = lines[i1][index1:index2].strip()
                    currencyCode = s.upper()
                else:
                    print('@@@', lines[i1])
                    currencyCode = None

            ## sector, industry
            if '<meta name="sector"' in lines[i1]:
                index2 = lines[i1].index('<meta name="sector"')

                index1 = index2+lines[i1][index2:].index('content="')+len(('content="'))
                index2 = index1+lines[i1][index1:].index('"')
                sector = lines[i1][index1:index2]

                index1 = index2+lines[i1][index2:].index('content="')+len(('content="'))
                index2 = index1+lines[i1][index1:].index('"')
                industry = lines[i1][index1:index2]

            ## market cap
            if '>MARKET CAP<' in lines[i1] or '>Market Cap<' in lines[i1]:

                index1 = lines[i1+1].index('>')+1
                index2 = lines[i1+1].rindex('<')
                factor = lines[i1+1][index2-1]
                s = lines[i1+1][index1:index2-1]

##                        index1 = lines[i1].index('MARKET CAP</div><div class="quoteData">')+len('MARKET CAP</div><div class="quoteData">')
##                        index2 = index1+lines[i1][index1:].index('<')
##                        factor = lines[i1][index2-1]
##                        s = lines[i1][index1:index2-1]

                if s == '-':
                    mc = None
                    statementNA = True
                else:
                    mc = float(s)*d_factors[factor]

        if name == '':
            print(len(lines))
            retry_name

        if lines == [] and price != '' and mc != '':
            stop_loop

        if mc == '' or price == '':
            statementNA = True

        return (
            name, currencyCode, price,
            sector, industry,
            mc,
            statementNA,
            )

if __name__=='__main__':
    instance_businessweek = BusinessWeek()

    t = instance_businessweek.parse_overview('BCR')
    print(t)
    stop

    url = 'http://investing.businessweek.com/businessweek/research/stocks/financials/financials.asp?ticker=COLR:BB&dataset=incomeStatement&period=A&currency=native'
    url = 'http://investing.businessweek.com/research/stocks/financials/financials.asp?ticker=012330:KS&dataset=incomeStatement&period=A&currency=native'
    url = 'http://investing.businessweek.com/research/stocks/financials/financials.asp?ticker=AXA&dataset=incomeStatement&period=A&currency=native'
    url = 'http://investing.businessweek.com/research/stocks/financials/financials.asp?ticker=012330:KS&dataset=incomeStatement&period=A&currency=native'
    t = instance_businessweek.parse_statement(url)
    print(t)
    stop

    t = instance_businessweek.parse_overview('BCR')
    print(t)
    stop

##    t = instance_reuters.parse_CEOcompensation('COP')
##    print t
##    stop

    ticker_businessweek = 'CBKS'
    stmt = 'incomeStatement'
    stmt = 'balanceSheet'
    stmt = 'cashFlow'
    url = 'http://investing.businessweek.com/businessweek/research/stocks/financials/financials.asp?ticker=%s&dataset=%s&period=A&currency=native' %(ticker_businessweek,stmt,)
    d,statement_error,statement_pending,currency = BusinessWeek().parse_statement(url,)
    print(d)
    print(statement_error)
    print(statement_pending)
    stop_end_succes

    d_msn2yahoo = {'GB':'L','FR':'PA','DE':'DE','AU':'AX','ES':'MC','JP':'T','IT':'MI','SE':'ST','BE':'BR','CA':'TO','NL':'AS',}
    d_currency = {}
    ticker_reuters = '0016.HK'

    url = 'http://in.reuters.com/money/quotes/quote?symbol=XOM'
    dic = {}
    statement = 'Balance'
    t = instance_reuters.parse_quote('XOM',)
    print(t)

