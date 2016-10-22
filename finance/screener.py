#!/bin/python3

## built-ins
import os, time, urllib
import import_tickers, parse_Yahoo, FAcandidates, parse_MSN, parse_Reuters

class finance:

    def main(self):

        ## initialize
        (
            tickers, d_months, time,
            d_indexes, d_ADR,
            ) = self.init()

        ## get currency rates
        d_currency_msn, d_currency_reuters = parse_Yahoo.Yahoo().parse_currencies()

        ## get previous FA candidates
        fd = open('l_FAs.txt','r')
        s = fd.read()
        fd.close()
        l_FAs = eval(s)

        ## make sure FAs are being checked
        for ticker in l_FAs:
            tickers.append(ticker)

        ## exclude duplicates and sort
        tickers = sorted(list(set(tickers)))

        d_portfolio = {
            'MSFT': 46.07, 'MMM': 164.48,
##            'DE:HEN3':0,'PG':73.63,
            'FLO': 0,
            'NOVOb.CO': 0,
            'HRL': 0,
            'GB:ULVR': 0, 'UL': 0,
            'RMCF': 0,
            'BVIC.L': 0,
            }
        for ticker in d_portfolio.keys():
            if ticker not in tickers:
                tickers.insert(0,ticker)
            if ticker not in l_FAs:
                l_FAs += [ticker]

##        tickers = list(set(tickers)-set(l_FAs))
##        tickers.sort()

##        tickers = 'GB:PFD'.split(',')

##        ## jump
##        i1 = 2000
##        i2 = 3000
##        print '***', i1, i2
##        tickers = tickers[i1:i2] ## 0-3500 clear, 4000-5000 clear

##        tickers = l_FAs[:5]
##        tickers.sort()

##        for ticker in list(tickers):
##            if ticker[:3] not in ['DE:',]:
##                tickers.remove(ticker)

####        for ticker in d_indexes['SSE50']:
####            print ticker
####            tickers.remove(ticker)

##        tickers = ['BDX',]

##        tickers = tickers[:int(len(tickers)/2)]
##        tickers = tickers[int(len(tickers)/2):]

##        tickers = tickers[0*int(len(tickers)/4):1*int(len(tickers)/4)]
##        tickers = tickers[1*int(len(tickers)/4):2*int(len(tickers)/4)]
##        tickers = tickers[2*int(len(tickers)/4):3*int(len(tickers)/4)]
##        tickers = tickers[3*int(len(tickers)/4):]

        if len(tickers) > 5:
            for fn in [
                'multiples.txt','statement_fundamentals.txt','redflag.txt',
                'statement_pending.txt','statement_error.txt',
                ]:
                fd = open(fn,'r')
                l = fd.readlines()
                fd.close()
                index = 0
    ##            if fn == 'fundamentals.txt':
    ##                l = '0ABCDEFGHIJKLMNOPQRSTUVWXYZA'
    ##                index = s.index(',%s' %(l[l.index(s[0])+1]))
    ##                s = s[index+1:]
    ##            else:
                if fn == 'statement_fundamentals.txt':
                    n = 20 ## long list
                elif fn == 'statement_pending.txt':
                    n = 10 ## keep checking
                else:
                    n = 5
                l = l[n:]
                fd = open(fn,'w')
                fd.writelines(l)
                fd.close()

        ##
        ## identify FA candidates
        ##
        (
            l_FAcandidates, FAdata, FAyahoo,
            l_statement_pending, l_statement_missing,
            l_multiples,
            ) = FAcandidates.FA().find_candidates_FA(
                d_portfolio, l_FAs,
                tickers,
                d_indexes,
                d_currency_msn, d_currency_reuters, 
                d_ADR,
                time,

                ##
                ## multiples
                ##
                ## http://www.multpl.com/
                PE_max = 20.,
                ## http://moneycentral.msn.com/investor/market/treasuries.aspx
#                div_yield_min = 0.001,
                div_yield_min = 0,
                PCF_max = 15.,
                ##
                ## fundamentals
                ##
                NPM_min = 0.01,
                ## http://www.investopedia.com/ask/answers/070914/are-companies-negative-return-equity-roe-always-bad-investment.asp
                ROE_min = 0.15, ## Buffettology (Mary Buffett) 15% minimum (or maybe S&P500 5year average http://moneycentral.msn.com/investor/invsub/results/compare.asp?Page=InvestmentReturns&Symbol=US%3aPG)
                ROA_min = 0.05, ## if high ROE but low ROA then likely a lot of debt, because ROA=inc/(ass-lia)
                ## debt is also kept in check by requiring low enterprise value
##                de_max = .85, ## BNI 2007 0.75 ETN 2008 0.75 BNI 2008Q4 0.82
                DE_max = 3.,  # dont care about debt as long as ROA good
                mc_min = 0.0,
##                payout_ratio_max = .99, ## KO FY2008 58% (does it really matter if equity growing and high ROE?) Nestle 62%  # payout relative to FCF!!! not cincome!!!
                ## High yield low payout ratio stocks outperformed high yield high payout ratio stocks by 8.2% per year from 1990 to 2006. So did low yield low payout and no yield.
                ## http://www.suredividend.com/8rules/
                payout_ratio_max = 0.8,
                payout_ratio_cash_flow_free_max = float('inf'),
                ## more multiples
                EVFCF_max = 15., ## ETN 2008Q2-Q3 EV/FCF=15
                EVOCF_max = 14., ## my own number...
                ## historic fundamentals
##                r2_min = 0.9**2,
                r2_min_sales = 0.95**2,
#                r2_min_equity = 0.9**2,
                r2_min_equity = 0,
            )
        
        print('parse ownership/holders of {}'.format(FAdata.keys()))
        FAdata = parse_MSN.MSN().parse_ownership(FAdata)
##        FAdata = parse_Yahoo.Yahoo().parse_ownership(FAdata)
        ## parse earnings dates
        print('parse earnings dates')
        FAdata = parse_Yahoo.Yahoo().parse_earnings_date(
            FAdata,
            )

    ##    FAdata = parse_insidertrading(FAdata, ticker)
        print('parse insider trading')
#        FAdata = parse_nasdaq.Nasdaq().parse_insider(FAdata)
        ## http://www.nasdaq.com/symbol/ko/insider-trades
        ## <th>Net Activity</th>
	## <td class="center">(625,471)</td>
	## <td class="center">(4,615,358)</td>

        import TAcandidates
        l_tickers = l_statement_pending+l_statement_missing
        print('parse historic dividends')
        l_TAcandidates, TAdata = TAcandidates.TA().find_candidates_TA(
            l_FAcandidates, time, d_months, l_tickers, d_portfolio,
            d_ADR,
            )
    ##    plot_charts(TAcandidates, TAdata, time)

        if len(FAdata.keys()) > 10:
            ## write data (FA and TA)
            self.write_data(
                FAdata, TAdata,
                'FA',
                d_portfolio, d_indexes, time, FAyahoo,
                l_statement_pending, l_FAs, d_ADR, l_statement_missing,
                l_multiples,
                d_months,
                )

            self.upload_data()
        else:
            print('FAdata', FAdata)

        return


    def upload_data(self):

        import ftplib
        print('uploading data')

        with open('password.txt') as fd:
            password = fd.read().strip()
        ftp = ftplib.FTP('ftp.proteinkemi.dk')
        ftp.login('proteinkemi.dk',password)
        ftp.cwd('finance')
        ftp.storlines('STOR ' + 'equities.htm', open('equities.htm','rb'))
        ftp.close()

        return


    def write_data(
        self,
        data, TAdata,
        prefix,
        d_portfolio, d_indexes, time, FAyahoo,
        l_statement_pending, l_FAs, d_ADR, l_statement_missing, l_multiples,
        d_months,
        ):

        year = time[0] ; month = time[1] ; day = time[2]

        ## write candidates sorted by industry to file
        lines = []
        lines += [
            '<html>\n',
            '<head>\n',
            '<title></title>\n',
            ]
        ## Add header stuff.
        with open('head.txt') as f:
            s = f.read()
        lines += [s]
        ## Add javascript to make table sortable.
        with open('sortable.js') as f:
            s = f.read()
        lines += [s]
        lines += [
            '</head>\n',
            '<body>\n',
            ]

        ## table header
        lines += [
            '<br><br>\n\n',
##            '<table border="1" class="tableWithFloatingHeader">\n',
            '<table border="1">\n',
            '\n',

##            ## First table header
##            '<tr>\n',
##            '<td colspan="2"><kbd></kbd></td>\n',
##
##            '<td colspan="5"><kbd>Valuation/Multiples</kbd></td>\n',
##
##            '<td colspan="7"><kbd>Fundamentals</kbd></td>\n',
##
##            '<td colspan="1"><kbd>Past</kbd></td>\n',
##            '<td colspan="1"><kbd>Future</kbd></td>\n',
##            '<td colspan="3"><kbd></kbd></td>\n',
##            '<td colspan="6"><kbd>Technical</kbd></td>\n',
##            '<td colspan="2"><kbd></kbd></td>\n',
##
##            '</tr>\n',


##            ## Second table header.
##            '<tr>\n',
##            '<td colspan="2"><kbd></kbd></td>\n',
##
##            '<td colspan="5"><kbd>Market Ratios</kbd></td>\n', ## e.g. P/E, dividend yield (market cap as numerator or denominator)
##
##            '<td colspan="1"><kbd></kbd></td>\n',
##            '<td colspan="1"><kbd>Debt Ratios</kbd></td>\n', ## e.g. debt equity ratio (balance sheet)
##            '<td colspan="3"><kbd>Profitability Ratios</kbd></td>\n', ## e.g. NPM, ROE (income stmt (earnings) / balance sheet)
##            '<td colspan="1"><kbd>Activity Ratios</kbd></td>\n', ## e.g. Asset Turnover (income stmt (revenue) / balance sheet)
##            '<td colspan="1"><kbd></kbd></td>\n',
##            
##            '<td colspan="13"><kbd></kbd></td>\n',
##            '</tr>\n',


            ## Begin table head.
            '<thead>\n',

##            '<tr>\n',
            '<tr style="background-color:#C0C0C0">\n',
##            '<tr style="position: fixed; background-color: grey;">\n'
            '<td><kbd>Ticker</kbd></td>\n',
            '<td><kbd>Market Cap. (Bil USD)</kbd></td>\n',
            '<td><kbd>P/E</kbd></td>\n',
            '<td><kbd>P/OCF</kbd></td>\n',
            '<td><kbd>EV/OCF</kbd></td>\n',
            '<td><kbd>EV/FCF</kbd></td>\n',
            '<td><kbd>Div. Yield (%)</kbd></td>\n',
            '<td><kbd>Payout Ratio (%)</kbd></td>\n',
            '<td><kbd>Debt / Equity Ratio</kbd></td>\n',
            ## Buffett likes a high E/NTA https://youtu.be/fJGkUSJBrXo?t=4m15s
            '<td><kbd>E/NTA (%)</kbd></td>\n',
##            '<td><kbd>ROIC (%)</kbd></td>\n',
            '<td><kbd>ROE (%)</kbd></td>\n',
            '<td><kbd>NPM (%)</kbd></td>\n',
            '<td><kbd>ATO (%)</kbd></td>\n',
            '<td><kbd>A/E</kbd></td>\n',
            '<td><kbd>r**2</kbd></td>\n',
            '<td><kbd>Earnings</kbd></td>\n',
            '<td><kbd>Industry</kbd></td>\n',
            '<td><kbd>Tags</kbd></td>\n',
            '<td><kbd>5% Holders (excl. common mutual funds; e.g. Fidelity, Vanguard, T. Rowe)</kbd></td>\n',

            '<td><kbd>52w chng %</kbd></td>\n',
            '<td><kbd>10y chng %</kbd></td>\n',
            '<td><kbd>52w high %</kbd></td>\n',
            '<td><kbd>52w low %</kbd></td>\n',
            '<td><kbd>MA50</kbd></td>\n',
            '<td><kbd>MA200</kbd></td>\n',

            '<td><kbd>Exec. comp. (Mil USD)</kbd></td>\n',
            '<td><kbd>Most recent ann. stmt</kbd></td>\n',
            '</tr>\n'

            ## End table head.
            '</thead>\n',

            ## Begin table body.
            '<tbody>\n',
            ]

        d_industries = {}
        for ticker in data.keys():
            sector = data[ticker]['sector']
            industry = data[ticker]['industry']
            if not sector in d_industries.keys():
                d_industries[sector] = {}
            if not industry in d_industries[sector].keys():
                d_industries[sector][industry] = []
            d_industries[sector][industry] += [ticker]
        l_sectors = list(sorted(d_industries.keys()))
        for i_sector in range(len(l_sectors)):
            sector = l_sectors[i_sector]
            print('sector', sector)
            l_industries = list(sorted(d_industries[sector].keys()))
            bool_border = True
            for i_industry in range(len(l_industries)):
                industry = l_industries[i_industry]
                tickers = list(sorted(d_industries[sector][industry]))
                for i_ticker in range(len(tickers)):
                    ticker_msn = ticker = tickers[i_ticker]

                    if ticker in l_statement_missing:
##                        print 'stmt missing', ticker
                        continue
                    if ticker in l_statement_pending:
##                        print 'stmt pending', ticker
                        continue

                    if bool_border == True:
                        s_style = 'style="border-top: 2px solid black;"'
                        bool_border = False
                    else:
                        s_style = ''

                    indexes = ''
                    for index in d_indexes:
                        if index in ['NYSE','IXIC',]:
                            continue

                        if ticker in d_indexes[index]:
                            indexes += '%s, ' %(index)
                    indexes = indexes[:-2]
                    ## not NYSE/IXIC if other US index
                    for USindex in ['S&P500','DJI30','DJC65','Nasdaq100','Russell1000',]:
                        if USindex in indexes:
                            indexes = indexes.replace('NYSE, ','')
                            indexes = indexes.replace(', NYSE','')
                            indexes = indexes.replace('IXIC, ','')
                            indexes = indexes.replace(', IXIC','')
                    ## not DJC65 if DJI30 or DJT20 or DJU15
                    for DJindex in ['DJI30','DJT20','DJU15']:
                        if DJindex in indexes:
                            indexes = indexes.replace('DJC65, ','')
                            indexes = indexes.replace(', DJC65','')
                    ## not Russell1000 if S&P500
                    for DJindex in ['S&P500']:
                        if DJindex in indexes:
                            indexes = indexes.replace('Russell1000, ','')
                            indexes = indexes.replace(', Russell1000','')

##                    s_ticker = ticker
##                    if ticker in d_ADR.keys():
##                        print ticker
##                        break

##                        for ADR in d_ADR.keys():
##                            if d_ADR[ADR] == ticker:
##                                s_ticker += ' (%s)' %(ADR)
##                                break

##                    bool_ADR = False
##                    for ADR in d_ADR.keys():
##                        if ticker in d_ADR[ADR]:
##                            print ticker, ADR, d_ADR[ADR]
####                            ticker = ADR
##                            bool_ADR = True
##                            break
##                    if bool_ADR == True:
##                        break

                    ## table row
                    lines += ['\n<tr>\n']
                    ## ticker,name
##                    lines += ['<td align="left"><kbd><a href="http://moneycentral.msn.com/companyreport?Symbol=%s">%5s</a>,%s</kbd></td>' %(ticker, ticker, data[ticker]['name'])]
                    if ticker in l_multiples:
                        color = '8080ff'
                    else:
                        color = 'ffffff'
                    lines += [
                        '<td align="left" bgcolor="#%s"><kbd><a href="http://www.reuters.com/finance/stocks/overview?symbol=%s">%s</a>,%s</kbd></td>' %(
                            color,data[ticker]['ticker_reuters'], ticker, data[ticker]['name'],
                            )]
                    lines += [
                        ## mc
                        '<td align="right"><kbd>%05.1f</kbd></td>' %(
                            round(data[ticker]['Market Cap (Mil)']/1000000000.,1),
                            ),
                        ## P/E
                        '<td><kbd><a href="http://moneycentral.msn.com/investor/invsub/results/compare.asp?Page=PriceRatios&Symbol={}">{:03.1f}</a></kbd></td>'.format(
                            ticker_msn, data[ticker]['P/E'],
                            ),
                        ]
                    ## P/CF
                    lines += ['<td align="right"><kbd><a href="http://moneycentral.msn.com/investor/invsub/results/compare.asp?Page=PriceRatios&Symbol=%s">%4.1f</a></kbd></td>' %(ticker, data[ticker]['P/CF'],)]
                    ## EV/FCF, div yield, payout ratio
                    lines += [
                        '<td align="right"><kbd>{:04.1f}</kbd></td>'.format(data[ticker]['EV/OCF']),
                        '<td align="right"><kbd>%05.1f</kbd></td>' %(data[ticker]['EV/FCF']),
                        '<td align="right"><kbd>%03.1f</kbd></td>' %(data[ticker]['Div. Yield']),
                        '<td align="right"><kbd>%04.1f</kbd></td>' %(data[ticker]['Payout Ratio']),
                        ]
                    ## debt/equity
                    lines += ['<td><kbd><a href="http://moneycentral.msn.com/investor/invsub/results/compare.asp?Page=FinancialCondition&Symbol=%s">%3.1f</a></kbd></td>' %(ticker, data[ticker]['D/Eq'],)]
                    ## fundamentals
                    lines += [
                        ## E/NTA
                        '<td align="right"><kbd>%05.1f</kbd></td>' %(data[ticker]['E/NTA'],),
##                        ## ROIC
##                        '<td align="right"><kbd>%4.1f</kbd></td>' %(data[ticker]['ROIC'],),
                        ## ROE
                        '<td align="right"><kbd><a href="http://moneycentral.msn.com/investor/invsub/results/compare.asp?Page=InvestmentReturns&Symbol=%s">%05.1f</a></kbd></td>' %(ticker, data[ticker]['ROE']),
                        ## NPM
                        '<td align="right"><kbd><a href="http://moneycentral.msn.com/investor/invsub/results/compare.asp?Page=ProfitMargins&Symbol=%s">%04.1f</a></kbd></td>' %(ticker, data[ticker]['NPM']),
                        ## Asset Turnover
                        '<td align="right"><kbd>%05.1f</kbd></td>\n' %(data[ticker]['ATO']),
                        ## Equity Multiplier
                        '<td align="right"><kbd>%04.1f</kbd></td>\n' %(data[ticker]['Equity Multiplier']),
                        ## r**2
                        '<td align="right"><kbd><a href="http://moneycentral.msn.com/investor/invsub/results/statemnt.aspx?Symbol=%s&lstStatement=10YearSummary&stmtView=Ann">%04.2f</a></kbd></td>' %(ticker,data[ticker]['r_sales_log**2']),
                        ]
                    ## earnings date
                    color = 'ffffff'
                    if data[ticker]['Date'] != 'N/A':
                        if d_months[data[ticker]['Date'][4:7]] == month:
                            color = '8080ff'
                    lines += [
                        ## earnings date
                        '<td bgcolor="#%s"><kbd>%11s</kbd></td>' %(color,data[ticker]['Date']),
                        ## industrty
                        '<td %s><kbd>%s</kbd></td>' %(s_style,data[ticker]['industry']),
                        ## indexes
                        '<td><kbd>%s</kbd></td>\n' %(indexes),
                        ]
                    ## 5% holder
                    lines += [
                        '<td><kbd>%s</kbd></td>\n' %(data[ticker]['holders']),
                        ]
                    if ticker not in TAdata.keys() or 'daily' not in TAdata[ticker].keys() or 'change_52w' not in TAdata[ticker]['daily'].keys(): ## Japan not on Yahoo
                        lines += [
                            '<td align="right"><kbd><font color="black">N/A</font></kbd></td>\n',
                            '<td align="right"><kbd><font color="black">N/A</font></kbd></td>\n',
                            '<td align="right"><kbd><font color="black">N/A</font></kbd></td>\n',
                            ]
                    else:
                        lines += ['<td align="right"><kbd><font color="black">%04.1f</font></kbd></td>\n' %(TAdata[ticker]['daily']['change_52w']),]
                        if 'change_10y' in TAdata[ticker]['daily'].keys():
                            lines += ['<td align="right"><kbd><font color="black">%05.1f</font></kbd></td>\n' %(TAdata[ticker]['daily']['change_10y']),]
                        else:
                            lines += ['<td align="right"><kbd><font color="black">N/A</font></kbd></td>\n',]
                        lines += ['<td align="right"><kbd><font color="black">%04.1f</font></kbd></td>\n' %(TAdata[ticker]['daily']['below_52w_max']),]
                        lines += ['<td align="right"><kbd><font color="black">%04.1f</font></kbd></td>\n' %(TAdata[ticker]['daily']['above_52w']),]
                    ## MA50
                    if not ticker in TAdata.keys() or 'daily' not in TAdata[ticker].keys() :
                        lines += [
                            '<td align="right"><kbd><font color="black">%s</font></kbd></td>\n' %('N/A'),
                            ]
                    elif TAdata[ticker]['daily']['MA50'] < data[ticker]['price']:
                        lines += [
                            '<td align="right"><kbd><font color="green">%.2f</font></kbd></td>\n' %(TAdata[ticker]['daily']['MA50']),
                            ]
                    else:
                        lines += [
                            '<td align="right"><kbd><font color="red">%.2f</font></kbd></td>\n' %(TAdata[ticker]['daily']['MA50']),
                            ]
                    ## MA200
                    if not ticker in TAdata.keys() or 'daily' not in TAdata[ticker].keys():
                        lines += [
                            '<td align="right"><kbd><font color="black">%s</font></kbd></td>\n' %('N/A'),
                            ]
                    elif TAdata[ticker]['daily']['MA200'] < data[ticker]['price']:
                        lines += [
                            '<td align="right"><kbd><font color="green">%.2f</font></kbd></td>\n' %(TAdata[ticker]['daily']['MA200']),
                            ]
                    else:
                        lines += [
                            '<td align="right"><kbd><font color="red">%.2f</font></kbd></td>\n' %(TAdata[ticker]['daily']['MA200']),
                            ]
                    ## executive compensation
                    lines += ['<td><kbd><a href="%s"></a>%s</kbd></td>\n' %(ticker,data[ticker]['compensation'])]
                    ## most recent statement
                    lines += ['<td><kbd>{}</kbd></td>\n'.format(
                        int(data[ticker]['min_year']/1000000))]
                    lines += ['</tr>\n']

        lines += ['</tbody>']

        lines += ['\n\n</table>\n\n']

        ## SEC filings
        s = ''
        for ticker in l_FAs:
            if ticker in l_statement_pending:
                s += '<a href="http://www.reuters.com/finance/stocks/overview?symbol=%s">%s</a>,' %(
                    ticker, ticker,
                    )
        if len(s) > 0:
            lines += [
                '\n<br>Waiting for SEC filings (10-K or similar) from the following companies:\n<br>\n %s' %(s)
                ]

        ## 10 year data
        lines += [
            '\n<br>10 year reports not available for the following non US companies:\n<br>\n'
            ]
        for ticker in l_FAs:
            if ticker in l_statement_missing:
                lines += [
                    '<a href="http://moneycentral.msn.com/investor/invsub/results/statemnt.aspx?Symbol=%s&lstStatement=10YearSummary&stmtView=Ann">%s</a>,' %(
                        ticker, ticker,
                        )
                    ]

        ## disclosure
        lines += [
            '\n<br>Banks and insurance companies are not included on this list, because their financial statements are not transparent.\n',
            ]

        ## Yahoo charts for all picks
        lines += [
            '<br>\n',
            '\n<a href="%s">Yahoo</a>\n' %(FAyahoo),
            ]

        lines += ['\n\n</body>\n</html>']

        fd = open('htm/%s%s%s%s.htm' %(prefix, year, str(month).zfill(2), str(day).zfill(2)), 'w')
        fd.writelines(lines)
        fd.close()

        fd = open('equities.htm', 'w')
        fd.writelines(lines)
        fd.close()

        return


    def delete_unwanted(self):

        cwd = os.getcwd()
        l_fn = os.listdir(cwd)
        for fn in l_fn:
            if fn[-4:] == '.pyc':
                os.remove(fn)

        localtime = time.localtime()
        mktime = time.mktime(localtime)
        dn = os.path.join(cwd,'urls',)
        for dn2 in os.listdir(dn):
            folder_sub = os.path.join(dn, dn2)
            if not os.path.isdir(folder_sub):
                continue
            print(dn2)
            l_fn = os.listdir(folder_sub)
            print('%i files in url dir' %(len(l_fn)))
            if len(l_fn) > 1000:
                for fn in l_fn:
    #                if 'ELEKTRA' in fn: ## tmp, cant delete
    #                    continue
                    if '*' in fn: ## tmp, cant delete
                        continue
                    ## modification time (seconds)
                    mtime = os.path.getmtime(os.path.join(folder_sub ,fn,))
                    ## remove file if older than 24 hours
                    if (mktime-mtime) > 60*60*24*1:
                        os.remove(os.path.join(folder_sub, fn))
                    else:
                        print(dn, fn, int((mktime-mtime)/(60*60*24)), 'days')

        return


    def init(self):

        ## delete unwanted downloaded urls taking up space
        self.delete_unwanted()

        ## get localtime
        localtime = time.localtime()

        all, d_indexes, d_ADR = import_tickers.main()

        d_months = {'Jan': 1, 'Feb': 2, 'Mar': 3, 'Apr': 4, 'May': 5, 'Jun': 6, 'Jul': 7, 'Aug': 8, 'Sep': 9, 'Oct': 10, 'Nov': 11, 'Dec': 12}

        return (
            all, d_months, localtime, d_indexes, d_ADR,
            )


    def read_url(self, url, ticker):

        if ticker[0] != '^':
            folder = ticker[0]
        else:
            folder = '0'
        fp = 'urls/%s/%s' %(
            folder, url.replace(':','').replace('/','').replace('.','').replace('?','').replace('*',''))

        for i in range(10):

            try:
                if (
                    os.path.isfile(fp)
                    and
                    time.localtime(os.stat(fp).st_mtime).tm_mday == time.localtime().tm_mday
                    and
                    time.localtime(os.stat(fp).st_mtime).tm_mon == time.localtime().tm_mon
                    and
                    time.localtime(os.stat(fp).st_mtime).tm_year == time.localtime().tm_year
                    ):
                    fd = open(fp,'r')
                    s = fd.read().rstrip()
                    fd.close()
                else:
                    print(url)
#                    s = urllib.request.urlopen(url).read().decode().rstrip()
#                    urllib.urlretrieve(url,fp)
#                    fd = open(fp,'w')
#                    fd.write(s)
#                    fd.close()

#                    urllib.urlretrieve(url,fp)
#                    with open(fp) as f:
#                        s = rd.read().rstrip()

                    response = urllib.request.urlopen(url)
                    s = response.read().decode().rstrip()
                    with open(fp,'w') as f:
                        f.write(s)

            except:
                s = ''
                pass

        lines = s.split('\n')

        return lines


if __name__=='__main__':
    instance_finance = finance()
    instance_finance.main()

## http://www.sec.gov/cgi-bin/browse-edgar?type=%s&dateb=&owner=include&count=40&action=getcompany&CIK=%s %(form_prefix='10-',ticker)

##http://www.fool.com/shop/newsletters/08/index.htm?source=iiiedilnk9251094
##6 Secrets of Dividend Investing:
##How You Can Earn Great Returns with Less Risk
##Finding the best dividend stocks takes some legwork and careful analysis. But here's how you can find the best long-term winners:
##1 Avoid the Highest Dividend Stocks - You can't pick stocks by dividend yield alone. Above-normal dividends are often a red flag for a company in distress. Studies have consistently shown that you will earn higher long-term returns by avoiding risky stocks with overly high dividends. 
##2 Beware the "Dividend Time Bombs" - Not all dividends are created equal. Even if a company has a generous dividend, it must be able to maintain it. A "doomed-to-be-cut" dividend can be worse than no dividend at all. Once a dividend is cut, it's likely to make the share price fall also. 
##3 Cash Is King - Free cash flow (FCF) is the true health of the business. Find the companies that generate tons of it. Even in the worst of times, those flush with greenbacks have options. Firms with cash can buy back their shares to raise stock prices, make their debt payments, increase dividends, and buy other profitable businesses. That's why cash flow is the single most important factor that determines value in the marketplace. 
##4 Don't Focus on Income without Growth - Only growing businesses are truly healthy. So cash flow needs to be strong enough to not only pay a healthy dividend but also generate enough cash to grow and stay strong strategically. 
##5 Don't Forget Value - An investment's total yield depends on both the dividend amount and the stock price. Stocks of companies making real products and real profits often don't make the headlines. So dividend stocks can also be a great source of hidden value. Finding value by focusing on dividends first can help you avoid catching the "falling knives" that trap some value investors. 
##6 Have a Longer-Term Focus - Many brokerage houses make investment recommendations based on a very short-term view of the world - often a maximum 12-month timeframe. Individual investors should have at least a three- to five-year view when considering investments. More time helps you fully realize the true power of compounding dividends.

## http://www.fool.com/investing/value/2005/04/07/profit-from-pessimism.aspx
## Richard Gibbons
## Generally, I like companies with enterprise value-to-free cash flow ratios of less than 15. I get extremely interested when it's less than 10.

## http://www.fool.com/investing/general/2005/11/29/panning-for-gold.aspx
## Rich Smith
## In a perfect world (meaning a recession, when almost all businesses are priced to sell), I prefer to invest in companies that sell for an EV/FCF of less than 10. This is not a hard-and-fast rule -- sometimes, you'll do quite well buying a fast-growing company at an EV/FCF of as high as 20. But be warned: The higher the EV/FCF, the greater the risk.

## http://www.fool.com/investing/small-cap/2005/01/05/7-steps-to-finding-hidden-gems.aspx
## Rich Smith
## 'And how do I define "bargain?" An EV/FCF ratio of 10 or less gets my attention real quick. Anything pricier than that, and I need to take a good hard look at the company's growth rate and EV/FCF/G ratio.'
