## rows not sorted by year...
## http://investing.money.msn.com/investments/financial-statements?symbol=aap

## mixed Mil and Bil
## http://investing.money.msn.com/investments/financial-statements?symbol=AU:CSL

import math, ticker_conversion, os, time
import screener

class MSN:

    def parse_financial_10_year_summary(self,url):

        '''return values in millions'''

        dic_10year_summary = {
            '&nbsp;':[],
##                'Sales':[],
##                'EBIT':[],
##                'Current Assets':[],
##                'Current Liabilities':[],
##                'Shares Outstanding':[],
            'SALES':[],
            'EBIT':[],
            'CURRENT ASSETS':[],
            'CURRENT LIABILITIES':[],
            'SHARES OUTSTANDING':[],
            }
        dic_10year_summary['DATE'] = []

        ## http://investing.money.msn.com/investments/financial-statements?symbol=aap not sorted by year...

        d_factors = {'Mil':1.,'Bil':1000.}

        lines = screener.finance().read_url(url, '0')

        for i1 in range(len(lines)):

##            if 'Balance Sheet - 10 Year Summary' in lines[i1] or 'Income Statement - 10 Year Summary' in lines[i1]:
            if 'INCOME STATEMENT: 10-YEAR SUMMARY' in lines[i1] or 'BALANCE SHEET: 10 YEAR SUMMARY' in lines[i1]:

                d_cols = {}
                col = 0
                bool_tr = False
                bool_td = False
                for i2 in range(i1+1,len(lines)):

                    if '</table>' in lines[i2]:
                        break
                        bool_tr = bool_td = False
                        d_cols = {}
                        continue
                    if '</tr>' in lines[i2]:
                        col = 0
                        bool_tr = False
                        continue
                    elif '<tr' in lines[i2]:
                        bool_tr = True
                        continue
                    elif '<td' in lines[i2] or '<th' in lines[i2]:
                        bool_td = True
                        continue
                    elif '</td>' in lines[i2] or '</th>' in lines[i2]:
                        col += 1
                        bool_td = False
                        continue

                    if (bool_tr == False or bool_td == False):
                        continue

##                    if '</td>' in lines[i2] and col+1 not in d_cols.keys():
                    ##
                    ## header
                    ##
                    if '<span' in lines[i2-1] and col not in list(d_cols.keys()):
##                        col += 1
##                        index2 = lines[i2].index('</td>')
##                        index1 = lines[i2][:index2].rindex('>')+1
##                        k = lines[i2][index1:index2]
                        k = lines[i2].strip().replace('<br />',' ').replace('<br/>',' ')
                        if k in list(dic_10year_summary.keys()):
                            d_cols[col] = k
                    ##
                    ## data
                    ##
                    elif '<span' in lines[i2-1] and col in list(d_cols.keys()):
##                        col += 1
##                        index2 = lines[i2].index('</td>')
##                        index1 = lines[i2][:index2].rindex('>')+1
##                        s = lines[i2][index1:index2]
                        s = lines[i2].strip()
                        if d_cols[col] == 'SHARES OUTSTANDING':
                            if 'Mil' in s or 'Bil' in s:
                                f = float(s[:-len(' Xil')])
                                f *= d_factors[s[-3:]]
                            elif s == '0.00':
                                f = 0.
                            elif s == 'NA':
                                f = 0.
                            ## thousands
                            elif ',' in s and '.' in s:
                                f = float(s.replace(',',''))
                            else:
                                print(url)
                                print(d_cols[col])
                                print(s)
                                stop
                        elif d_cols[col] == 'DATE':
                            ## don't add dates again and assume same sequence in both (income/balance) tables
                            if 'BALANCE SHEET: 10 YEAR SUMMARY' in lines[i1]:
                                continue
                            f = '%s/%s' %(s[-2:],s[:2],)
                        else:
##                            s = s.replace('Bil','').replace('Mil','')
                            ## year
                            if d_cols[col] == '&nbsp;':
                                f = int(s[-2:])
                                if f > 90:
                                    stop
                                    f += 1900
                                else:
                                    f += 2000
                            elif s == 'NA':
                                f = 0.
                            ## float
                            else:
                                if 'Mil' in s or 'Bil' in s:
                                    f = float(s.replace(',','')[:-len(' Xil')])
                                    f *= d_factors[s[-3:]] ## make sure not a ratio with a number from financial stmt later on...
                                else:
                                    f = float(s.replace(',',''))
                        dic_10year_summary[d_cols[col]] += [f]

                ## break when second table is reached
                if 'BALANCE SHEET: 10 YEAR SUMMARY' in lines[i1]:
                    break

        ##
        ## MSN sorts dates and data by month/year; sort by year/date instead
        ##
        l_dates = list(dic_10year_summary['DATE'])
        l_dates.sort()
        l_dates.reverse()
        l_indexes = [dic_10year_summary['DATE'].index(date) for date in l_dates]
        for k in list(dic_10year_summary.keys()):
            if len(dic_10year_summary[k]) == 0:
                continue
            l = [dic_10year_summary[k][index] for index in l_indexes]
            dic_10year_summary[k] = l

        return dic_10year_summary


    def key_ratios_10_year_summary(self,ticker,dic_10year_summary):

        url = 'http://moneycentral.msn.com/investor/invsub/results/compare.asp?Page=TenYearSummary&symbol=%s' %(ticker)
        url = 'http://investing.money.msn.com/investments/key-ratios?symbol=%s&page=TenYearSummary' %(ticker)

        lines = screener.finance().read_url(url, ticker)

        ##
        ##
        ##
        d_cols = {}
        bool_init = 0
        bool_tr = False
        bool_td = False
        row = 0
        col = 0
        for i1 in range(len(lines)):

            if (
                ' is not available.</p>' in lines[i1]
                or
                '<span>SEARCH RESULTS</span>' in lines[i1]
                ):
                break

            elif (
                ' AVG P/E<' in lines[i1+8]
                or
                ('<table' in lines[i1] and bool_init == 1)
                ):
                bool_init += 1
                continue

            elif '</table' in lines[i1] and bool_init == 2:
                break

            elif bool_init == 0:
                continue

            elif '<tr' in lines[i1]:
                bool_tr = True
                row += 1
                continue

            elif '</tr' in lines[i1]:
                bool_tr = False
                col = 0
                continue

            elif bool_tr == False:
                continue

            elif '<td' in lines[i1] or '<th' in lines[i1]:
                bool_td = True
                col += 1
                continue

            elif '</td' in lines[i1] or '</th' in lines[i1]:
                bool_td = False
                continue

            elif bool_td == False:
                continue

            elif '<span' in lines[i1]:
                continue

            elif lines[i1].strip() == '</span>':
                continue

            if '</span>' in lines[i1]:
                s = lines[i1].strip()[:-7].replace('<br/>',' ')
                d_cols[col] = s
                dic_10year_summary[s] = []
            else:
                if col == 1:
                    continue
                s = lines[i1].strip()
                s = s.replace(',','') ## thousand separator
                if s == 'NA':
                    continue
                if col != 1:
                    s = float(s)
                dic_10year_summary[d_cols[col]] += [s]


##
##        d_cols = {}
##        for i1 in range(len(lines)):
##
##            if not ' AVG P/E<' in lines[i1]:
##                continue
##
##            i2 = 0
##            index1_tr = 0
##            while '<tr>' in lines[i1][i2:][index1_tr:]:
##                print 'c'
##                index1_tr = index1_tr+lines[i1][index1_tr:].index('<tr>')+len('<tr>')
##                index2_tr = index1_tr+lines[i1][index1_tr:].index('</tr>')
##                col = 0
##                index1_tx = 0
##                while '</t' in lines[i1][index1_tr:index2_tr][index1_tx:]:
##                    index2_tx = index1_tx+lines[i1][index1_tr:index2_tr][index1_tx:].index('</t')
##                    index1_tx = index1_tx+lines[i1][index1_tr:index2_tr][index1_tx:index2_tx].index('<t')
##                    s = lines[i1][index1_tr:index2_tr][index1_tx:index2_tx]
##                    while s[0] == '<':
##                        s = s[s.index('>')+1:]
##                    if lines[i1][index1_tr:index2_tr][:4] == '<th>':
##                        d_cols[col] = s
##                        dic_10year_summary[s] = []
##                    elif lines[i1][index1_tr:index2_tr][:4] == '<td>':
##                        if col != 0:
##                            if s != 'NA':
##                                ## characters (currency symbols) in front
##                                while len(s) > 0 and s[0] not in '-123456789':
##                                    s = s[1:]
####                                if s[1] not in '0.-123456789':
####                                    s = s[0]+s[2:]
##                                ## thousand seperator
##                                s = s.replace(',','')
##                                ## characters behind
##                                while len(s) > 0 and s[-1] not in '0123456789':
##                                    s = s[:-1]
##                                if len(s) > 0:
##                                    s = float(s)
##                            dic_10year_summary[d_cols[col]] += [s]
##                    else:
##                        print lines[i1][index1_tr:index2_tr][:4]
##                        stop
##                    index1_tx = index2_tx+1
##                    col += 1
##                index1_tr += 1
##
##            break

        return dic_10year_summary


    def parse_company_report(self,ticker,rate,):

        url = 'http://moneycentral.msn.com/companyreport?Symbol='+ticker
        lines = screener.finance().read_url(url, ticker)

        d = {
            'ma50':'N/A',
            'ma200':'N/A',
            'relative strength':'N/A',
            }

        for i in range(len(lines)):

            line = lines[i]

            if 'Exchange : ' in lines[i]:
                index = lines[i].index('Exchange : ')
                index2 = index+lines[i][index:].index('</b>')
                index1 = lines[i][:index2].rindex('>')+1
                exchange = lines[i][index1:index2]
    ##                if exchange == 'OTC BB':
    ##                    break

            if 'Last Price' in lines[i] and '<meta ' not in lines[i]:
                index1 = lines[i+1].index('<td>')+4
                index2 = lines[i+1].index('</td>')
                price = float(lines[i+1][index1:index2].replace(',',''))
                d['price'] = price

            if '50 Day Moving Average' in lines[i]:
                index1 = lines[i+1].index('<td>')+4
                index2 = lines[i+1].index('</td>')
                s = lines[i+1][index1:index2]
                if s == 'NA':
                    ma50 = 'N/A'
                else:
                    ma50 = float(s.replace(',',''))
                d['ma50'] = ma50

            if '200 Day Moving Average' in lines[i]:
                index1 = lines[i+1].index('<td>')+4
                index2 = lines[i+1].index('</td>')
                s = lines[i+1][index1:index2]
                if s == 'NA':
                    ma200 = 'N/A'
                else:
                    ma200 = float(s.replace(',',''))
                d['ma200'] = ma200

            if ': Company Report</' in line:
                index2 = line.index(': Company Report</')
                index1 = line[:index2].rindex('>')+1
                name = line[index1:index2]
                d['name'] = name

            if 'Volatility (beta)' in lines[i]:
                index2 = lines[i+1].index('</')
                index1 = lines[i+1][:index2].rindex('>')+1
                s = lines[i+1][index1:index2]
                if s == 'NA':
                    beta = 'N/A'
                else:
                    beta = float(s.replace(',',''))
                d['beta'] = beta

            if '<td>Sales</td>' in lines[i]:
                index2 = lines[i+2].index('</td>')
                index1 = lines[i+2][:index2].index('>')+1
                s = lines[i+2][index1:index2]
                if s == 'NA':
                    print('sales 5y N/A (maybe because negative)')
                    growth_sales_5y = 'N/A'
                elif '-' in s or '<span ' in s:
                    growth_sales_5y = 0
                else:
                    growth_sales_5y = float(s[:-1])
                d['growth_sales_5y'] = growth_sales_5y

            if '<td>Income</td>' in lines[i]:
                index2 = lines[i+2].index('</td>')
                index1 = lines[i+2][:index2].index('>')+1
                s = lines[i+2][index1:index2]
                if s == 'NA':
                    print('income 5y N/A (maybe because negative)')
                    growth_income_5y = 'N/A'
                elif '-' in s:
                    growth_income_5y = 0
                else:
                    growth_income_5y = float(s[:-1])
                d['growth_income_5y'] = growth_income_5y

            if 'Market Capitalization' in lines[i] and '<meta ' not in lines[i]:
                index = lines[i].index('Market Capitalization')
                index1 = index+lines[i][index:].index('<td>')+len('<td>')
                index2 = index1+lines[i][index1:].index('</td>')
                if lines[i][index1:index1+2] == 'NA':
                    stop_mc
                else:
                    factor = lines[i][index2-3:index2]
                    if factor not in ['Bil','Mil']:
                        factor = .001
                        mc = factor*float(lines[i][index1:index2])/rate
                    else:
                        if factor == 'Bil':
                            factor = 1000000000.
                        elif factor == 'Mil':
                            factor = 1000000.
                        mc = factor*float(lines[i][index1:index2-4])/rate
                d['mc'] = mc

            if '<td>Dividend Yield</td>' in lines[i]:
                index2 = lines[i+1].index('</td>')
                index1 = lines[i+1][:index2].index('>')+1
                s = lines[i+1][index1:index2]
                if s == 'NA':
                    div_yield = 'N/A'
                else:
                    div_yield = float(s[:-1].replace(',',''))/100
                d['div_yield'] = div_yield

            if '<td>Debt/Equity Ratio</td>' in lines[i]:
                index = lines[i].index('Debt/Equity Ratio')
                index1 = index+lines[i][index:].index('<td>')+len('<td>')
                index2 = index1+lines[i][index1:].index('</td>')
                s = lines[i][index1:index2]
                if s == 'NA':
                    debt_equity_ratio = 'N/A'
                else:
                    debt_equity_ratio = float(lines[i][index1:index2])
                d['debt_equity_ratio'] = debt_equity_ratio

            if '<td>Last 12 Months</td>' in lines[i]:
                s = lines[i+2]
                while '<' in s:
                    index1 = s.index('>')+1
                    index2 = s.rindex('<')
                relative_strength = int(s.replace('%%',''))
                d['relative strength'] = relative_strength

        return d


    def parse_statement_yearly(self,url,dic,statement,d_currency):

        dic_out = {}

        lines = screener.finance().read_url(url, '0')
        statementNA = False

        for i in range(len(lines)):

            line = lines[i]

            ## statement not available
            if 'The financial statement for this symbol, is currently not available.' in line:
                statementNA = True
                s_business = 'Industry'
                return dic, statementNA, s_business
            if 'Statement information for this ticker symbol is not available at this time.' in line:
                statementNA = True
                s_business = 'Industry'
                return dic, statementNA, s_business

            ## business type
            if '                        <p><b>Business Type:</b> <span id="lblBusinessType">' in line:
                i2 = line.index('</span></p>')
                i1 = line[:i2].rindex('>')+1
                s_business = line[i1:i2]
                if s_business in ['Bank','Insurance',]:
                    return dic_out, statementNA, s_business

            ## currency
            if 'Financial data in' in line:
                index1 = line.index('>')+1
                index2 = index1+line[index1:].index('<')
                s = line[index1:index2]
                rate = d_currency[s]
                dic_out['rate'] = rate

            ## multiple
            if 'Values in ' in line and ' (Except for per share items)' in line:
                index1 = line.index('Values in')+len('Values in')
                index2 = line.index('(Except for per share items)')
                s = line[index1:index2]
                d_factors = {'Millions':1000000.,'Thousands':1000.}
                factor = d_factors[s.strip()]
                dic_out['factor'] = factor

            ## filing dates
            if 'class="ftable"' in line:
                periods = []
                index = line.index('<tr class="r1">')
                index += line[index:].index('</td>')+1
                for quarter in range(5):
                    index2 = index+line[index:].index('</td>')
                    index1 = line[:index2].rindex('>')+1
                    index = index2+1
                    s = line[index1:index2]
                    if s == '': ## e.g. GB:SSE
                        periods += [None]
                    else:
                        periods += [[int(s[:4]),int(s[-1:])]]

            ## statement sources
            if '>Stmt Source<' in line:
                column1,stmt_source = self.parse_stmt_source(line)
                dic_out['column1'] = column1
                dic_out['periods'] = periods
                dic_out['source'] = stmt_source

                if stmt_source == 'PRESS':
                    stop
                    statementNA = True
                    break ## break line loop

            for key in dic:
                if '>%s<' %key in line:
                    ## key already in another row?
##                    if dic_out[key] not in ['N/A',None]:
##                        print key, dic[key]
##                        print stmt_source
##                        stop
                    if key in ['Total Common Shares Outstanding','Total Preferred Shares Outstanding']:
                        dic_out[key] = self.parse_statement_multiple('>%s<' %(key), line, 1., factor)
                    else:
                        dic_out[key] = self.parse_statement_multiple('>%s<' %(key), line, rate, factor)
                    if dic_out[key][column1] in ['N/A',None]:
                        print(key)
                        stop
                        break

        return dic_out, statementNA, s_business


    def parse_statement_quarterly(self,url,dic,statement, d_currency):

        statementNA = False

        lines = screener.finance().read_url(url, '0')

        ## reset dictionary
        for key in list(dic.keys()):
            dic[key] = None

        for line in lines:
            if 'The financial statement for this symbol, is currently not available.' in line:
##                stop1
                statementNA = True
                s_business = 'Industry'
                return dic, statementNA, s_business
            if 'Statement information for this ticker symbol is not available at this time.' in line:
##                stop2
                statementNA = True
                s_business = 'Industry'
                return dic, statementNA, s_business

        for i in range(len(lines)):

            line = lines[i]

            if '<span id="lblErrorMessage"><br><br>Statement information for this ticker symbol is not available at this time.<br><br><br></span>' in line:
                statementNA = True
                break

            ## business type
            if '                        <p><b>Business Type:</b> <span id="lblBusinessType">' in line:
                i2 = line.index('</span></p>')
                i1 = line[:i2].rindex('>')+1
                s_business = line[i1:i2]
                if s_business in ['Bank','Insurance',]:
                    return dic, statementNA, s_business

            ## currency
            if 'Financial data in' in line:
                index1 = line.index('>')+1
                index2 = index1+line[index1:].index('<')
                s = line[index1:index2]
                rate = d_currency[s]
                dic['rate'] = rate

            ## multiple
            if 'Values in ' in line and ' (Except for per share items)' in line:
                index1 = line.index('Values in')+len('Values in')
                index2 = line.index('(Except for per share items)')
                s = line[index1:index2]
                d_factors = {'Millions':1000000.,'Thousands':1000.}
                factor = d_factors[s.strip()]
                dic['factor'] = factor

            ## filing dates
            if 'class="ftable"' in line:
                periods = []
                index = line.index('<tr class="r1">')
                index += line[index:].index('</td>')+1
                for quarter in range(5):
                    index2 = index+line[index:].index('</td>')
                    index1 = line[:index2].rindex('>')+1
                    index = index2+1
                    s = line[index1:index2]
                    if s == '': ## e.g. GB:SSE
                        periods += [None]
                    else:
                        periods += [[int(s[:4]),int(s[-1:])]]

            ## period lengths
            if statement == 'cashflow' and '>Period Length<' in line:
                period_lengths = []
                index = line.index('>Period Length<')
                index += line[index:].index('</td>')+1
                for quarter in range(5):
                    index2 = index+line[index:].index('</td>')
                    index1 = line[:index2].rindex('>')+1
                    index = index2+1
                    s = line[index1:index2]
                    if s != s.strip():
                        print(s)
                        notexpected
                    period_lengths += [s]

                cf_columns_y,cf_columns_qoq = self.columns_of_interim_statement(periods,period_lengths)
                dic['cf_columns_y'] = cf_columns_y
                dic['cf_columns_qoq'] = cf_columns_qoq

            ## statement sources
            if '>Stmt Source<' in line:
                column1,stmt_source = self.parse_stmt_source(line)
                dic['column1'] = column1
                dic['periods'] = periods
                dic['source'] = stmt_source

                if stmt_source == 'PRESS':
                    stop
                    statementNA = True
                    break ## break line loop

            for key in dic:
                if '>%s<' %key in line:
                    ## key already in another row?
                    if dic[key] not in ['N/A',None]:
                        print(key, dic[key])
                        print(stmt_source)
                        stop
                    if key in ['Total Common Shares Outstanding','Total Preferred Shares Outstanding']:
                        dic[key] = self.parse_statement_multiple('>%s<' %(key), line, 1., factor)
                    else:
                        dic[key] = self.parse_statement_multiple('>%s<' %(key), line, rate, factor)
                    if dic[key][column1] in ['N/A',None]:
                        print(key)
                        stop
                        break

        return dic, statementNA, s_business


    def parse_financial_highlights(self,ticker):

        url = 'http://moneycentral.msn.com/investor/invsub/results/hilite.asp?Symbol=%s' %(ticker)
        lines = screener.finance().read_url(url, ticker)

        for i in range(len(lines)):

            line = lines[i]

            if '<td>Payout Ratio</td>' in line:
                index1 = lines[i+1].index('<td>')+len('<td>')
                index2 = lines[i+1].index('</td>')
                s = lines[i+1][index1:index2]
                if s == 'NA':
                    payout_ratio = 'N/A'
                else:
                    payout_ratio = float(s.replace('%',''))/100.

        return payout_ratio


    def parse_stmt_source(self,line):

        index = line.index('>Stmt Source<')
        index += line[index:].index('</td')+1
        for quarter in range(5):
            index2 = index+line[index:].index('</td')
            index1 = line[:index2].rindex('>')+1
            index = index2+1
            source = line[index1:index2]
            if source in [
                'PRESS', ## e.g. accounts payable on balance sheet missing
                'N/A',
                '', ## e.g. RCL for which only 3 half year periods are listed
                ]:
                continue
            elif source in [
                ## http://en.wikipedia.org/wiki/SEC_filings
                ## http://edgar.sec.gov/info/edgar/forms/edgform.pdf
                '10-K','10-Q', ## yearly, quarterly
                '10-K/A','10-Q/A', ## Amendment
                '10KSB','10QSB', ## Small Business
                '20-F','20-F/A', ## ADR level II and III 10-Q/10-K equivalent
                '6-K', ## Report of foreign issuer [Rules 13a-16 and 15d-16]
                '8-K', ## extraordinary event / current report ## total assets = 0 (e.g. COP)
                'Yuho','Yuho/A', ## Japan annual & interim (12,6 month report)
                'Tanshin','Tanshin/A', ## Japan quarterly report
                'ARS','ARS/A', ## Annual Report to Security Holders
    ##                        'N-CSRS','N-CSR','N-CSR/A','N-30D',
                'PROSPECTUS',
                'Interim Report','Interim Report/A' ## not always trustworthy if ADR !!!
                ]:
                break
            elif source in ['8-K/A','6-K/A']:
                print(source)
                stop
            else:
                print(source)
                stop

        if quarter > 2:
            print(quarter)
            stop_balance_sheet

        return quarter, source


    def parse_ownership(self,data):

        print('parsing ownership')

        tickers = list(data.keys())
        tickers.sort()

        for i in range(len(tickers)):

            ticker = tickers[i]
            ticker_msn = ticker_conversion.yahoo2msn(ticker)
            if i % 10 == 0:
                print('\n%s/%s %s' %(i+1, len(tickers), ticker))

            url = 'http://moneycentral.msn.com/ownership?symbol=%s&Holding=5%%25+Ownership' %(ticker_msn)
            url = 'http://investing.money.msn.com/investments/five-percent-ownership?symbol={}'.format(ticker_msn)

            lines = screener.finance().read_url(url, ticker)

            if lines == [''] or lines == []:
                print('no data', url)
                data[ticker]['holders'] = ''
                continue

            data[ticker]['holders'] = []
            for i in range(len(lines)):
                line = lines[i]
                if 'Ownership' and 'Holder Name' in line:
                    if 'No data available' in line:
                        break
                    index = line.index('Holder Name')
                    while '<tr>' in line[index:]:
                        index += line[index:].index('<tr>')
                        for i_td in range(4):
                            index2 = index+line[index:].index('</td>')
                            index1 = line[:index2].rindex('>')+1
                            index = index2+1
                            s = line[index1:index2]
                            if i_td == 0:
                                holder = s
                            elif i_td == 3:
                                percentage = float(s)
                        s = holder
                        if holder in [
                            ## public
                            'GAMCO Investors, Inc.', ## 1977, GBL, Bill Gates
                            'Royce & Associates, LLC', ## 1972, 1899 LM, Legg Mason acquisition 2001
                            'Franklin Advisory Services, LLC', ## 1947, BEN, Franklin Templeton Investments, franklintempleton.com
                            'State Street Global Advisors (US)', ## 1978, 1792, STT
                            ## private
                            'Neuberger Berman, LLC', ## 1939, private / Lehman Brothers
                            'Dimensional Fund Advisors, LP', ## 1981, private (Scholes, Merton)
                            'Capital World Investors', ## 1931
                            'Capital Research Global Investors', ## 1931
                            'Renaissance Technologies Corp.', ## 1982, private
                            ## policy holder owned (mutual)
                            'State Farm Insurance Companies', ## 1922
                            ## LLP
                            'Wellington Management Company, LLP', ## 1928
                            ## public or private? probably private...
                            'Lord, Abbett & Co. LLC', ## 1929
                            'Keeley Asset Management Corp.', ## 1982, keeleyasset.com
                            'Wells Capital Management Inc.', ## Wells Fargo???
                            'First Eagle Global Fund', ## firsteaglefunds.com
                            'First Eagle Investment Management LLC', ## firsteaglefunds.com
    ##                            'Keeley Small Cap Value Fund, Inc.', ## 1982, keeleyasset.com
    ##                            'Heartland Value Fund', ## heartlandfunds.com
    ##                            'Baron Capital Management, Inc.', ## 1982, Ronald S Baron
    ##                            'Artisan Partners Limited Partnership', ## artisanfunds.com

    ##                            'Morgan Stanley Investment Management Inc. (US)',
    ##                            'Goldman Sachs Asset Management (US)',
    ##                            'J.P. Morgan Investment Management Inc. (New York)',
    ##                            'Dodge & Cox', ## dodgeandcox.com 1930
    ##                            'Wells Fargo Advantage Small Cap Value Fund',
    ##                            'Davis Selected Advisers, L.P.',
                            
    ##                            'Ruane, Cunniff & Goldfarb, Inc.', ## Sequoia Fund
    ##                        'US Trust', ## ustrust.com 1853
    ##                        'Perry Capital', ## perrycap.com
                            ]:
                            continue
                        elif (
    ##                            (
    ##                                'Berkshire' not in s
    ##                                and
    ##                                'Walton' not in s
    ##                                )
    ##                            and
                            (
                                s[:len('Vanguard ')] == 'Vanguard ' ## client owned
                                or
                                s[:len('Fidelity ')] == 'Fidelity ' ## 1946, private
                                or
                                s[:len('American Funds ')] == 'American Funds '
                                or
                                s[:len('Columbia ')] == 'Columbia ' ## Columbia Management Group, Ameriprise Financial (AMP) subsidiary
                                or
                                s[:len('BlackRock ')] == 'BlackRock ' ## BLK
                                or
                                s[:len('T. Rowe Price ')] == 'T. Rowe Price ' ## 1937, TROW
                                or
                                s[:len('Ruane, Cunniff & Goldfarb, Inc. ')] == 'Ruane, Cunniff & Goldfarb, Inc. ' ## 1969, owns the Sequoia Fund [SEQUX]
                                
    ##                            'JPMorgan Chase' in s or
    ##                            'Lord Abbett' in s or
    ##    ##                        'Financial' in s or
    ##    ##                        'International' in s or
    ##    ##                        'Mgmt' in s or
    ##    ##                        'Plc' in s or 
    ##    ##                        'Associates' in s or ## 'T Rowe Price Associates', ## troweprice.com 1937

    ##                            'Management' in s or
    ##                            'Partners' in s or
    ##                            'Advisors' in s or
    ##                            'Investors' in s or
    ##                            'Investment' in s or

    ##    ####                        'REPUBLIC' in s.upper() or
    ##    ####                        'KINGDOM' in s.upper() or
    ##    ##                        'Capital' in s or
    ##    ##                        'Group' in s or

    ##                            s[-10:] == ' Companies' or
    ##                            s[-8:] == ' Company' or
    ##                            s[-11:] == ' Management' or
    ##                            s[-8:] == ' Managem' or
    ##                            s[-5:] == ' Bank' or
    ##                            s[-6:] == ' Trust' or ## http://en.wikipedia.org/wiki/Trust_company
    ##                            s[-9:] == ' Partners' or
    ##                            s[-8:] == ' Savings' or ## Applied Industrial Technologies Retirement Savings
    ##                            s[-5:] == ' ESOP' or ## employee stock ownership plan
    ##                            s[-7:] == ' (ESOP)' or ## employee stock ownership plan
    ##                            ## corporations
    ##                            s[-12:] == ' Corporation' or
    ##                            s[-6:] == ' Corp.' or
    ##                            s[-4:] == ' Co.' or
    ##                            s[-3:] == ' Co' or
    ##                                or
    ##                            ## limited
    ##                            s[-4:] == ' Ltd' or ## Limited (commonwealth)
    ##                            s[-5:] == ' Ltd.' or ## Limited (commonwealth)
    ##                                s[-4:] == ' LLC'
    ##                                or
    ##                                s[-7:] == ' L.L.C.'
    ##                                or
    ##                                s[-4:] == ' LLP'
    ##                                or
    ##                                s[-3:] == ' LP'
    ##                                or
    ##                                s[-5:] == ' L.P.'
    ##                                or
    ##                                s[-len(' Limited Partnership'):] == ' Limited Partnership'
    ##                                )
                                )
                            ):
                            ## check that a private person is not being excluded
                            if (
                                ## a name?
                                ('(' in s or ')' in s) and
                                ## not a name!
                                '(Grove Creek)' not in s and
                                '(US)' not in s and
                                '(UK)' not in s and
                                '(Americas)' not in s and
                                '401(k)' not in s and
                                '(New York)' not in s and
                                '(Switzerland)' not in s and
                                '(Singapore)' not in s and
                                '(International)' not in s
                                and
                                '(ESOP)' not in s
                                ):
                                print(s)
                                stop
                            continue
                        else:
                            data[ticker]['holders'] += ['%s (%.1f)' %(holder,percentage,)]
                            if (
                                '(' not in s and ')' not in s
                                and 'Berkshire' not in s
                                and 'Walton' not in s
                                ):
                                fd = open('investors_notexpected.txt','a')
                                fd.write('%s\t%s\t%s\n' %(holder.split()[-1], ticker, holder))
                                fd.close()

                    ## break loop over lines
                    break
                            
            data[ticker]['holders'] = ', '.join(data[ticker]['holders'])

        return data


    def parse_insidertrading(self,data, tickers):

        print('parsing insider trading')

        for i in range(len(tickers)):

            ticker = tickers[i]
            print(('\n%s/%s' %(i+1, len(tickers)), ticker))

            url = 'http://moneycentral.msn.com/investor/invsub/insider/trans.asp?Symbol=%s' %(ticker)
            lines = screener.finance().read_url(url, ticker)
            for i in range(len(lines)):
                line = lines[i]
                if 'Recent Insider Trading Activity' in line:
                    print(line)
                    index = 0
                    for j in range(3):
                        index += line[index:].index('<tr')
                    for k in range(7):
                        index += line[index:].index('<td')

        return data


    def parse_statement_multiple(self, keyword, line, rate, factor, l_columns_qoq = None, l_columns_qoq_1y = None):

        l_values = []
        index = line.index(keyword)
        index += line[index:].index('<td')
        for column in range(5):
            index2 = index+line[index:].index('</td>')
            index1 = line[:index2].rindex('>')+1
            index = index2+1
            s = line[index1:index2].replace(',','')
            if s in ['','N/A']:
                l_values += [0]
            else:
                l_values += [factor*float(s)/rate]

        return l_values


    def columns_of_interim_statement(self, quarters, period_lengths):

        ## set variables
        columns = []
        column = -1
        months = 0
        weeks = 0
        columns_skip = []
        ## loop over line
        while months != 12 and weeks != 52 and column != 4:
            column += 1
            s = period_lengths[column]
            if s == '' and column in [1,3]: ## e.g. GB:SSE
                columns += [0]
                continue
            elif s == '':
                stopandaddcommentticker
                print(column)
                notexpected
            elif s == 'N/A': ## e.g. PFE
                columns += [0]
                continue
            else:
                period = s.split()[1]
                length = int(s.split()[0])
                ## adjust for 53rd week of Q4
                if period == 'Weeks' and length in [14,53] and quarters[column][1] == 4:
                    length = 13
                ## ignore column
                if column in columns_skip:
                    columns += [0]
                    continue
                ## add column
                elif months < 12 and weeks < 53:
                    columns += [1]
                    if period == 'Months':
                        months += int(length)
                    elif period == 'Weeks':
                        weeks += int(length)
                    else:
                        notexpected
                    ## skip column i+1?
                    if (
                        (period == 'Months' and (
                            (quarters[column][1] == 2 and length == 6) or
                            (quarters[column][1] == 3 and length == 9) or
                            (quarters[column][1] == 4 and length in [6,12])
                            )
                         ) or
                        (period == 'Weeks' and (
                            (quarters[column][1] == 2 and length in [26,27]) or ## e.g. UNF
                            (quarters[column][1] == 3 and length == 39) or
                            (quarters[column][1] == 4 and length in [26,52])
                            )
                         )
                        ):
                        if period == 'Months':
                            quarter_min = quarters[column][1]-length/3
                        elif period == 'Weeks':
                            quarter_min = quarters[column][1]-length/13
                        for q in range(quarters[column][1]-1,quarter_min,-1):
                            quarter = [quarters[column][0],q]
                            if quarter in quarters:
                                column_skip = quarters.index(quarter)
                                if period_lengths[column_skip] == 'N/A':
                                    columns_skip += [column_skip]
                                elif (
                                    months-int(period_lengths[column_skip].split()[0]) < 12 and
                                    weeks-int(period_lengths[column_skip].split()[0]) < 52
                                    ):
                                    columns_skip += [column_skip]
                ## subtract column
                elif months > 12 or weeks > 53:
                    ## subtract column
                    if ((period == 'Months' and int(length) <= months) or (period == 'Weeks' and int(length) <= weeks)):
                        columns += [-1]
                        if period == 'Months':
                            months -= int(length)
                        elif period == 'Weeks':
                            weeks -= int(length)
                    else:
                        notexpected
                ## ignore column
                elif s == '12 Weeks' and weeks in [52,53]:
                    columns += [0]
                    stop_find_example
                else:
                    print(s)
                    print((weeks, months))
                    print(columns)
                    print((column, columns_skip))
                    notexpected

        print((quarters, period_lengths))
        columns_quarter = []
        for i in range(len(period_lengths)):
            s_1 = period_lengths[i]
            if s_1 in ['N/A','']: ## e.g. TM,AU:SGP
                continue
            period_1 = s_1.split()[1]
            length_1 = int(s_1.split()[0])
            if len(columns_quarter) == 2:
                break
            for j in range(i,min(i+3,len(period_lengths))):
                s_2 = period_lengths[j]
                if s_2 in ['N/A','']: ## e.g. TM,AU:SGP
                    continue
                period_2 = s_2.split()[1]
                length_2 = int(s_2.split()[0])
                if period_1 == 'Months' and period_2 == 'Months':
                    if i == j: ## e.g. JP:9104
                        if i < 3 and j < 3:
                            continue
                        if i == 3 and j == 3 and length_1 == 6 and length_2 == 6: ## e.g. JP:9104
                            columns_quarter += [[-i]]
                        else:
                            print((i,j, length_1, length_2))
                            print((s_1, s_2))
                            stop
                    else:
                        if length_1-length_2 == 3: ## e.g. JNJ
                            columns_quarter += [[i,-j]]
                        elif length_1 == 3 and length_2 == 12: ## e.g. ACN (1st quarter)
                            columns_quarter += [[i]]
                        elif length_1 == 6 and length_2 == 12: ## e.g. TM (1st half year)
                            columns_quarter += [[i]]
                        elif length_1 == 12 and length_2 == 6: ## e.g. TM (2nd half year)
                            columns_quarter += [[i,-j]]
                        elif length_1 == 12 and length_2 == 12: ## e.g. PKX
                            columns_quarter += [[i,-j]]
                        else:
                            print((length_1, length_2))
                            print((s_1,s_2))
                            print((quarters, period_lengths))
                            stop
                elif period_1 == 'Weeks' and period_2 == 'Weeks': ## e.g. HD
                    if i == j:
                        if i < 3 and j < 3:
                            continue
                        else:
                            stop
                    else:
                        if length_1-length_2 in [12,13,14]: ## e.g. HD,MSM,COST
                            columns_quarter += [[i,-j]]
                        elif length_1 in [52,53,] and length_2 in [36,40,]: ## e.g. PEP,AAP,UNF
                            columns_quarter += [[i,-j]]
                        elif length_1 in [36,40,] and length_2 in [24,28,]: ## e.g. PEP,AAP
                            columns_quarter += [[i,-j]]
                        elif length_1 in [12,13,14,16,] and length_2 in [52,53]: ## e.g. OXM,CHKE,WFMI,STX (1st quarter)
                            columns_quarter += [[i]]
                        elif length_1 in [24,25,26,27,] and length_2 in [52,53,]: ## e.g. GB:TSCO,GB:MRW (1st half year)
                            columns_quarter += [[i]]
                        elif length_1 in [52,53,] and length_2 in [25,26,27,]: ## e.g. GB:TSCO,GB:MRW,AU:WOW (1st half year)
                            columns_quarter += [[i,-j]]
                        else:
                            print(quarters)
                            print(period_lengths)
                            print(('***', length_1, length_2))
                            stop
                else:
                    print((length_1, length_2))
                    print((period_1, period_2))
                    print((quarters, period_lengths))
                    stop
                break

        if len(columns_quarter) != 2:
            print(columns_quarter)
            stop

        return columns, columns_quarter


if __name__ == '__main__':

    instance_MSN = MSN()

    d = {}
    url = 'http://moneycentral.msn.com/investor/invsub/results/statemnt.aspx?Symbol=SNP&lstStatement=10YearSummary&stmtView=Ann'
    url = 'http://investing.money.msn.com/investments/financial-statements?symbol=PETM'
    dic_10year = {
        '&nbsp;':[],
##        'Sales':[],
        'SALES':[],
        'EBIT':[],
        'Current Assets':[],
        'Current Liabilities':[],
        'Shares Outstanding':[],
        'CURRENT ASSETS':[],
        'CURRENT LIABILITIES':[],
        'SHARES OUTSTANDING':[],
        }
    d = instance_MSN.parse_financial_10_year_summary(url,dic_10year,)
    print(d)
    stop

    d = instance_MSN.key_ratios_10_year_summary('MSFT',{},)
    print(d)
    stop

    d = {'ABT':{}}
    d = instance_MSN.parse_ownership(d,)
    print(d)
    stop_end

    import screener
    instance_finance = screener.finance()
    (
        tickers, months, time,
        d_indexes, d_msn2yahoo, d_msn2currency, d_ADR, d_yahoo2reuters,
        ) = instance_finance.init()
    d = {}
    for ticker in tickers:
        d[ticker] = {}
    d = instance_MSN.parse_ownership(d,)

    d2 = {}
    for ticker in list(d.keys()):
        holders = d[ticker]['holders']
        for holder in holders:
            if not holder in list(d2.keys()):
                d2[holder] = 0
            d2[holder] += 1

    for holder in list(d2.keys()):
        if d2[holder] > 50:
            print(holder, d2[holder])
