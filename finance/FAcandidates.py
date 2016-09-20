#!/bin/env python3

# Tommy Carstensen, 2005-2016

import math
import parse_morningstar
import parse_MSN
import parse_Reuters
import parse_BusinessWeek
import statistics
import ticker_conversion
import screener
import requests
import bs4
import parse_Yahoo
import parse_advfn
import parse_ft
import sys


source = 'ADVFN'
source = 'Yahoo'


class FA:

    def parse_gurufocus(self, url):

        response = requests.get(url)
        soup = bs4.BeautifulSoup(response.text)

        table = soup.find('table', attrs={'id': 'Rf'})
        if not table:
            return None
        table_body = table.find('tbody')
        rows = table_body.find_all('tr')
        d = {}
        for row in rows:
            cols = row.find_all('td')
            if not cols:
                continue
            if not cols[0].find('a'):
                continue
            if not cols[0].find('a').contents:
                continue
            k = cols[0].find('a').contents[0]
            v = []
            for col in cols[1:]:
                v.append(col.find('div').contents[0])
            d[k] = v

        return d

    def find_candidates_FA(
        self, d_portfolio, l_FAs, tickers,
        d_indexes,
        d_currency_name, d_currency_symbol,
        d_ADR,
        time,
        PE_max, div_yield_min,
        payout_ratio_max, payout_ratio_cash_flow_free_max,
        PCF_max, NPM_min, ROE_min, ROA_min, DE_max, mc_min,
        EVFCF_max, EVOCF_max,
        r2_min_sales, r2_min_equity,
        ):

        FAcandidates = []
        d_FAdata = {}
        skip = []
        d_prices = {}
        l_multiples = []

        d_income_year = {}
        d_balance_year = {}
        d_cash_year = {}
        d_statement_year = {}

        #
        # Reuters
        #
        d_income_in = {
            'Total Revenue': None,  # sum
            'Gross Profit': None,  # sum
                'Unusual Expense (Income)': None,
            'Operating Income': None,  # sum
                'Gain (Loss) on Sale of Assets': None,
            'NET INTEREST EXPENSE': None,  # subtract
            'Net Income Before Extra. Items': None,
            'Net Income': None,  # sum
            'Total Adjustments to Net Income': None,
            'Normalized Income Available to Common': None,
            }
        d_balance_in = {
                    'Cash and Short Term Investments': None,
                        'Cash & Equivalents': None,
                        'Short Term Investments': None,
                    'Total Receivables, Net': None,
                        'Accounts Receivable - Trade, Net': None,
                        'Notes Receivable - Short Term': None,
                    'Total Inventory': None,
                    'Prepaid Expenses': None,
                    'Other Current Assets, Total': None,
                'Total Current Assets': None,
                    'Property/Plant/Equipment, Total - Net': None,
                    'Goodwill': None,
                    'Goodwill, Net': None,
                    'Intangibles, Net': None,
                    'Long Term Investments': None,
                    'Note Receivable - Long Term': None,
                    'Other Long Term Assets, Total': None,
                    'Other Assets, Total': None,
            'Total Assets': None,
            # current liabilities
                    'Accounts Payable': None,
                    'Payable/Acrued': None,
                    'Acrued Expenses': None,
                    'Notes Payable/Short Term Debt': None,
                    'Total Short Term Borrowings': None,
                    'Current Port. of LT Debt/Capital Leases': None,
                    'Other Current Liabilities, Total': None,
                'Total Current Liabilities': None,
                    'Total Long Term Debt': None,
                        'Long Term Debt': None,
                    'Deferred Income Tax': None,
                    'Minority Interest': None,
                    'Other Liabilities, Total': None,
            'Total Liabilities': None,
                    'Redeemable Preferred Stock': None,
                    'Preferred Stock - Non Redeemable, Net': None,
                    'Common Stock': None,
                    'Additional Paid-In Capital': None,
                    'Retained Earnings (Accumulated Deficit)': None,
                    'Treasury Stock - Common': None,
                    'Unrealized Gain (Loss)': None,
                    'Other Equity, Total': None,
            'Total Equity': None,
            'Total Common Shares Outstanding': None,
            'Total Preferred Shares Outstanding': None,
            }
        d_cash_in = {
            'Net Income/Starting Line': None,
            'Depreciation/Depletion': None,
            'Amortization': None,
            'Changes in Working Capital': None,
            'Capital Expenditures': None,
            'Total Cash Dividends Paid': None,
            'Cash from Operating Activities': None,
            }
        d_businessweek2reuters = {
            # income statement
            'TOTAL REVENUES': 'Total Revenue',
            'OPERATING INCOME': 'Operating Income',
            'NET INCOME': 'Net Income',
            'NET INCOME TO COMMON EXCLUDING EXTRA ITEMS':
            'Net Income Before Extra. Items',
            'Cost of Goods Sold': 'Cost of Revenue, Total',
            # balance sheet
            'TOTAL EQUITY': 'Total Equity',
            'TOTAL ASSETS': 'Total Assets',
            'TOTAL CASH AND SHORT TERM INVESTMENTS': 'Cash and Short Term Investments',
            'Long-Term Debt': 'Total Long Term Debt',
            # cash flow
            'CASH FROM OPERATIONS': 'Cash from Operating Activities',
            'TOTAL DIVIDEND PAID': 'Total Cash Dividends Paid',
            'Capital Expenditure': 'Capital Expenditures',
            }

        # skip
        if len(tickers) > 5:
            for path in (
                'statement_fundamentals.txt', 'multiples.txt',
                'statement_error.txt', 'redflag.txt',
                'index_financials.txt',
                ):
                with open(path) as f:
                    skip += f.read().split('\n')
            with open('statement_pending.txt') as f:
                l_statement_pending = fd.read().split('\n')
        else:
            l_statement_pending = []
        l_statement_missing = []

        for j in range(len(tickers)):

            FA = True

            ticker = tickers[j]
            if ticker in skip+l_statement_pending+l_statement_missing:
                continue
            if ticker == '\n':
                continue
##            if ticker != 'GB:RB.': continue
##            if ticker != 'AAPL': continue

            if len(sys.argv) == 2:
                if ticker[:2] < sys.argv[-1]:
                    continue  # tmp!!!

            print('\n', j+1, '\\', len(tickers), ticker)

            ticker_reuters = ticker_conversion.unknown2reuters(ticker)
#            ticker_businessweek = ticker_conversion.unknown2businessweek(ticker)
            ticker_msn = ticker_conversion.yahoo2msn(ticker)
            ticker_advfn = ticker_conversion.unknown2advfn(ticker)
            ticker_yahoo = ticker_conversion.msn2yahoo(ticker)
            ticker_ft = ticker_conversion.unknown2ft(ticker)

            if any([
                ticker_yahoo.endswith('.KL') and ticker_yahoo[0].isdigit,  # need to translate tickers from numerical (Yahoo) to alphabetical (FT/Yahoo)
                ticker_yahoo.endswith('.MC'),  # Madrid not on ft.com
                '.KS' in ticker,  # Korea not on ft.com
                ]):
                continue

            if any([
#                ticker_yahoo.endswith('.SZ'),
#                ticker_yahoo.endswith('.SS'),
#                ticker_yahoo.endswith('.AB'),
#                ticker_yahoo.endswith('.IC'),
#                ticker_yahoo.endswith('.VI'),  # Viennea/Austria
#                ticker_yahoo.endswith('.HK'),  # Hong Kong covered by FT, but I can't identify ISO code for Hong Kong exchange. Not SEHK, XHKG, etc.
#                ticker_yahoo.endswith('.RG'),  # Latvia/Riga not covered by FT
#                ticker_yahoo.endswith('.MC'),  # Madrid Comercial, dissapeared from ft...
                ]):
                # It's not possible to scrape Yahoo any longer!
                source = 'Yahoo'
                index0 = 0  # new to old
                index1 = 1
#                ## Statements from Mexico not in Yahoo Finance.
#                if '.MX' in ticker_yahoo:
#                    continue
#                ## Russia not in Reuters...
#                if '.ME' in ticker_yahoo:
#                    continue
#                ## Lots of countries without financial statements in Yahoo... :(
#                if '.TL' in ticker_yahoo:
#                    continue
#                if '.HK' in ticker_yahoo:
#                    continue
#                if '.TO' in ticker_yahoo:
#                    continue
#                if '.SA' in ticker_yahoo:
#                    continue
#                if '.MC' in ticker_yahoo:
#                    continue
            elif True or ':' in ticker or '.' in ticker:
                source = 'FT'
                index0 = 0  # new to old
                index1 = 1
            else:
                source = 'ADVFN'
                index0 = -1  # old to new
                index1 = -2
                index0 = 0  # new to old
                index1 = 1

            set_indexes = set()
            for index in d_indexes:
                if ticker in d_indexes[index]:
                    set_indexes |= set([index])

            #
            # parse income statement, annual
            #

##            ## Reuters only reports one year without login. I need 2 for YoY calculations.
##            url = 'http://www.reuters.com/finance/stocks/incomeStatement?stmtType=INC&perType=ANN&symbol=%s' %(ticker_reuters)
##            (
##                d_income_year[ticker], statement_error, currency_stmt,
##                ) = parse_Reuters.Reuters().parse_statement(url, d_income_in, 'income',)

##            ## BusinessWeek was acquired by Bloomberg in 2010
##            url = 'http://investing.businessweek.com/businessweek/research/stocks/financials/financials.asp?ticker=%s&dataset=incomeStatement&period=A&currency=native' %(ticker_businessweek)
##            url = 'http://investing.businessweek.com/research/stocks/financials/financials.asp?ticker=%s&dataset=incomeStatement&period=A&currency=native' %(ticker_businessweek)
##            (
##                d_income_year[ticker], statement_error, statement_pending,
##                currency_stmt,
##                ) = parse_BusinessWeek.BusinessWeek().parse_statement(url,)

            if source == 'Yahoo':
                ## Yahoo only does 3 years. But easy to parse and most markets. Do Morningstar or ADVFN for 10y.
                ## Dividends paid missing for non-US companies.
                ## Missing value and zero value identical; i.e. a dash
                url = url_yahoo = 'http://finance.yahoo.com/q/is?s={}+Income+Statement&annual'.format(ticker_yahoo)
                try:
                    d_income_year[ticker], statement_error, currency_stmt = parse_Yahoo.Yahoo().parseIncomeStatement(url_yahoo)
                except ValueError:
                    statement_error = True
                if statement_error is False:
                    if 'Net Income' not in d_income_year[ticker]:
                        statement_error = True
                    elif '-' in d_income_year[ticker]['Net Income'][:2]:
                        print('Net Income', d_income_year[ticker]['Net Income'])
                        statement_error = True
                    elif '-' in d_income_year[ticker]['Operating Income or Loss'][:2]:
                        statement_error = True
                    elif '-' in d_income_year[ticker]['Total Revenue'][:2]:
                        statement_error = True

            if source == 'ADVFN':
                # ADVFN does annual data back to 1993, but only for US markets etc.
    ##            url_advfn = http://uk.advfn.com/p.php?pid=financials&btn=start_date&mode=annual_reports&symbol=NYSE%3AIBM&start_date=18
                if 'LSE:' not in ticker_advfn and ':' not in ticker_advfn:
                    urls = (  # : is %3A
                        'http://www.advfn.com/p.php?pid=financials&btn=start_date&mode=annual_reports&symbol=NYSE:{}'.format(ticker_advfn),
                        'http://www.advfn.com/p.php?pid=financials&btn=start_date&mode=annual_reports&symbol=NASDAQ:{}'.format(ticker_advfn),
                        )
                else:
                    urls = ('http://www.advfn.com/p.php?pid=financials&btn=start_date&mode=annual_reports&symbol={}'.format(ticker_advfn),)
                for url in urls:
                    url_advfn = url
                    (
                        d_income_year[ticker], d_balance_year[ticker], d_cash_year[ticker],
                        statement_error,
                        ## currency_name
                        currency_stmt, d_indicators,
                        ) = parse_advfn.ADVFN().parse_statements(url_advfn)
                    ## break loop if correct stock exchange
                    if statement_error is False:
                        break

            if source == 'FT':
##                if '.KS' in ticker:
##                    continue
##                if '.TLV' in ticker:
##                    continue
##                if '.SET' in ticker:
##                    continue
                if ':' not in ticker_ft and any([
                    not '.' in ticker_ft,
                    ticker_ft.endswith('.A'),
                    ticker_ft.endswith('.B'),
                    ]):
                    for suffix_ft in (
                        ':NYQ',  # NYSE n=3220
                        ':NSQ',  # Nasdaq Global Select n=1612
                        ':NAQ',  # Nasdaq Capital Market n=798
                        ':NMQ',  # Nasdaq Global Market n=732
                        ':ASQ',  # NYSE MKT LLC
                        ):
                        url = url_ft = 'http://markets.ft.com/research/Markets/Tearsheets/Financials?s={}&subview=IncomeStatement'.format(ticker_ft+suffix_ft)
                        url = url_ft = 'http://markets.ft.com/data/equities/tearsheet/financials?s={}&subview=IncomeStatement&periodType=a'.format(ticker_ft+suffix_ft)
                        parse_ft.FT().parse_stmt(url_ft, ticker_ft)
                        try:
                            d_income_year[ticker], statement_error, currency_stmt = parse_ft.FT().parse_stmt(url_ft, ticker_ft)
                            ticker_ft += suffix_ft
                            break
                        except:
                            continue
                    else:
                        statement_error = True
##                        url = url_ft = 'http://markets.ft.com/research/Markets/Tearsheets/Financials?s={}&subview=IncomeStatement'.format(ticker_ft)
##                        d_income_year[ticker], statement_error, currency_stmt = parse_ft.FT().parse_stmt(url_ft)
                elif ticker.endswith('.VX'):
                    for suffix_ft in (':SWX', ':VTX', ':BRN'):
                        url = url_ft = 'http://markets.ft.com/data/equities/tearsheet/financials?s={}&subview=IncomeStatement&periodType=a'.format(ticker[:-3]+suffix_ft)
                        print(suffix_ft, ticker_ft, ticker, url)
                        parse_ft.FT().parse_stmt(url_ft, ticker_ft)
                        try:
                            d_income_year[ticker], statement_error, currency_stmt = parse_ft.FT().parse_stmt(url_ft, ticker_ft)
                            ticker_ft = ticker[:-3]+suffix_ft
                            print('xxxxxxxxxx', ticker_ft)
                            break
                        except:
                            continue
                    else:
                        parse_ft.FT().parse_stmt(url_ft, ticker_ft)
                    parse_ft.FT().parse_stmt(url_ft, ticker_ft)
                else:
                    url = url_ft = 'http://markets.ft.com/research/Markets/Tearsheets/Financials?s={}&subview=IncomeStatement'.format(ticker_ft)
                    url = url_ft = 'http://markets.ft.com/data/equities/tearsheet/financials?s={}&subview=IncomeStatement'.format(ticker_ft)
                    d_income_year[ticker], statement_error, currency_stmt = parse_ft.FT().parse_stmt(url_ft, ticker_ft)

                try:
                    d_income_year[ticker]['Cost of Revenue, Total'] = d_income_year[ticker]['Cost of revenue total']
                    d_income_year[ticker]['Total Revenue'] = d_income_year[ticker]['Total revenue']
                    d_income_year[ticker]['Operating Income'] = d_income_year[ticker]['Operating income']
                    d_income_year[ticker]['Net Income Before Extra. Items'] = d_income_year[ticker]['Net income before extra. Items']
                except KeyError:
                    pass

##            ## MorningStar reports 5 years free and 10 years paid
##            ## but impossible to parse?
##            ticker_morningstar = ticker_conversion.msn2morningstar(ticker)
####            url = 'http://financials.morningstar.com/income-statement/is.html?t=XLON:ADN&region=GBR'
##            url_morningstar = 'http://financials.morningstar.com/income-statement/is.html?t={}'.format(ticket_morningstar)

##            ## key conversion
##            for key in list(d_income_year[ticker].keys()):
##                d_income_year[ticker][key].reverse()
##                if key in list(d_businessweek2reuters.keys()):
##                    d_income_year[ticker][d_businessweek2reuters[key]] = d_income_year[ticker][key]
####                    del d_income_year[ticker][key]

            #
            # check income statement before proceeding with parsing other statements
            #
            if statement_error is True:
                self.stmtERROR(ticker, l_FAs=l_FAs, url=url)
                print('income statement not available')
                continue
##            if statement_pending is True:
##                print(ticker, 'income statement pending')
##                fd = open('statement_pending.txt', 'a')
##                fd.write('%s,' %(ticker))
##                fd.close()
##                l_statement_pending += [ticker]
##                continue

            if source == 'Yahoo':
                url_yahoo = 'http://finance.yahoo.com/q/bs?s={}+Balance+Sheet&annual'.format(ticker_yahoo)
                d_balance_year[ticker], statement_error, currency_stmt = parse_Yahoo.Yahoo().parseIncomeStatement(url_yahoo)
                url_yahoo = 'http://finance.yahoo.com/q/cf?s={}+Cash+Flow&annual'.format(ticker_yahoo)
                d_cash_year[ticker], statement_error, currency_stmt = parse_Yahoo.Yahoo().parseIncomeStatement(url_yahoo)
                if d_income_year[ticker]['Cost of Revenue'][index0] in ('-', '--',):
                    self.stmtERROR(ticker, l_FAs=l_FAs)
                    print('continue', ticker, 'income, Cost of Revenue')
                    continue
                if d_cash_year[ticker]['Capital Expenditures'][index0] in ('-', '--',):
                    self.stmtERROR(ticker, l_FAs=l_FAs)
                    print('continue', ticker, 'cash, Capital Expenditures')
                    continue

            if source == 'FT':
                url = url_ft = 'http://markets.ft.com/research/Markets/Tearsheets/Financials?s={}&subview=BalanceSheet'.format(ticker_ft)
                url = url_ft = 'http://markets.ft.com/data/equities/tearsheet/financials?s={}&subview=BalanceSheet'.format(ticker_ft)
                try:
                    d_balance_year[ticker], statement_error, currency_stmt = parse_ft.FT().parse_stmt(url_ft, ticker_ft)
                except:
                    ## Probably a fund or delisted company...
                    self.stmtERROR(ticker, l_FAs=l_FAs)
                    continue
                url = url_ft = 'http://markets.ft.com/research/Markets/Tearsheets/Financials?s={}&subview=CashFlow'.format(ticker_ft)
                url = url_ft = 'http://markets.ft.com/data/equities/tearsheet/financials?s={}&subview=CashFlow'.format(ticker_ft)
                try:
                    d_cash_year[ticker], statement_error, currency_stmt = parse_ft.FT().parse_stmt(url_ft, ticker_ft)
                except:
                    ## Probably a fund or delisted company...
                    self.stmtERROR(ticker, l_FAs=l_FAs)
                    continue

                d_balance_year[ticker]['Total Equity'] = d_balance_year[ticker]['Total equity']
                d_balance_year[ticker]['Total Assets'] = d_balance_year[ticker]['Total assets']
                d_balance_year[ticker]['Total Long Term Debt'] = d_balance_year[ticker]['Total long term debt']
                d_cash_year[ticker]['Total Cash Dividends Paid'] = d_cash_year[ticker]['Total cash dividends paid']
                d_cash_year[ticker]['Cash from Operating Activities'] = d_cash_year[ticker]['Total cash from operations']
                d_cash_year[ticker]['Capital Expenditures'] = d_cash_year[ticker]['Capital expenditures']

####            ##
####            ## parse balance sheet, annual
####            ## 
##            url = 'http://www.reuters.com/finance/stocks/incomeStatement?stmtType=BAL&perType=ANN&symbol=%s' %(ticker_reuters)
##            d_balance_year[ticker],statement_error,currency_stmt = parse_Reuters.Reuters().parse_statement(url, d_balance_in, 'balance',)
####            url = 'http://investing.businessweek.com/businessweek/research/stocks/financials/financials.asp?ticker=%s&dataset=balanceSheet&period=A&currency=native' %(ticker_businessweek)
####            d_balance_year[ticker],statement_error,statement_pending,currency_stmt = parse_BusinessWeek.BusinessWeek().parse_statement(url,)
####            for key in list(d_balance_year[ticker].keys()):
####                d_balance_year[ticker][key].reverse()
####                if key in list(d_businessweek2reuters.keys()):
####                    d_balance_year[ticker][d_businessweek2reuters[key]] = d_balance_year[ticker][key]
####                    del d_balance_year[ticker][key]
####
####            ##
####            ## parse cash flow, annual
####            ## 
##            url = 'http://www.reuters.com/finance/stocks/incomeStatement?stmtType=CAS&perType=ANN&symbol=%s' %(ticker_reuters)
##            d_cash_year[ticker],statement_error,currency_stmt = parse_Reuters.Reuters().parse_statement(url, d_cash_in, 'cashflow',)
####            url = 'http://investing.businessweek.com/businessweek/research/stocks/financials/financials.asp?ticker=%s&dataset=cashFlow&period=A&currency=native' %(ticker_businessweek)
####            d_cash_year[ticker],statement_error,statement_pending,currency_stmt = parse_BusinessWeek.BusinessWeek().parse_statement(url,)
####            for key in list(d_cash_year[ticker].keys()):
####                d_cash_year[ticker][key].reverse()
####                if key in list(d_businessweek2reuters.keys()):
####                    d_cash_year[ticker][d_businessweek2reuters[key]] = d_cash_year[ticker][key]
####                    del d_cash_year[ticker][key]

            if 'total equity' in list(d_balance_year[ticker].keys()):
                d_balance_year[ticker]['Total Equity'] = d_balance_year[ticker]['total equity']
            if 'total assets' in list(d_balance_year[ticker].keys()):
                d_balance_year[ticker]['Total Assets'] = d_balance_year[ticker]['total assets']

            if 'cash from operations' in list(d_cash_year[ticker].keys()):
                d_cash_year[ticker]['Cash from Operating Activities'] = d_cash_year[ticker]['cash from operations']
            if 'Net Cash From Total Operating Activities'.lower() in list(d_cash_year[ticker].keys()):
                d_cash_year[ticker]['Cash from Operating Activities'] = d_cash_year[ticker]['Net Cash From Total Operating Activities'.lower()]
            if 'total dividend paid' in list(d_cash_year[ticker].keys()):
                d_cash_year[ticker]['Total Cash Dividends Paid'] = d_cash_year[ticker]['total dividend paid']
            if 'cash dividends paid' in list(d_cash_year[ticker].keys()):
                d_cash_year[ticker]['Total Cash Dividends Paid'] = d_cash_year[ticker]['cash dividends paid']
            if 'Payment Of Cash Dividends'.lower() in list(d_cash_year[ticker].keys()):
                d_cash_year[ticker]['Total Cash Dividends Paid'] = d_cash_year[ticker]['Payment Of Cash Dividends'.lower()]
            if 'capital expenditure' in list(d_cash_year[ticker].keys()):
                d_cash_year[ticker]['Capital Expenditures'] = d_cash_year[ticker]['capital expenditure']

            if 'total revenues' in list(d_income_year[ticker].keys()):
                d_income_year[ticker]['Total Revenue'] = d_income_year[ticker]['total revenues']
            if 'total revenue' in list(d_income_year[ticker].keys()):
                d_income_year[ticker]['Total Revenue'] = d_income_year[ticker]['total revenue']
            if 'operating income' in list(d_income_year[ticker].keys()):
                d_income_year[ticker]['Operating Income'] = d_income_year[ticker]['operating income']
            if 'net income to common excluding extra items' in list(d_income_year[ticker].keys()):
                d_income_year[ticker]['Net Income Before Extra. Items'] = d_income_year[ticker]['net income to common excluding extra items']
                d_income_year[ticker]['Net Income To Common Excluding Extra Items'] = d_income_year[ticker]['net income to common excluding extra items']
            if 'Total Net Income'.lower() in list(d_income_year[ticker].keys()):
                d_income_year[ticker]['Net Income Before Extra. Items'] = d_income_year[ticker]['Total Net Income'.lower()]
##            ## add to total net income!!!
##            if 'Extraordinary Income/Losses'.lower() in list(d_income_year[ticker].keys()):
##                if d_income_year[ticker]['Extraordinary Income/Losses'.lower()][0] != 0:
##                    print(ticker, d_income_year[ticker]['Extraordinary Income/Losses'.lower()])
##                    stop
            if 'cost of goods sold' in list(d_income_year[ticker].keys()):
                d_income_year[ticker]['Cost of Revenue, Total'] = d_income_year[ticker]['cost of goods sold']
            ## Conversion from advfn.com
            if 'Cost Of Sales'.lower() in d_income_year[ticker].keys():
                d_income_year[ticker]['Cost of Revenue, Total'] = d_income_year[ticker]['Cost Of Sales'.lower()]
            ## Conversion from advfn.com
            try:
                d_cash_year[ticker]['Purchase Of Property, Plant & Equipment'.lower()] = d_cash_year[ticker]['Purchase Of Property & Equipment'.lower()]
            except KeyError:
                pass
            if source == 'ADVFN':
                print(d_cash_year[ticker]['Purchase Of Property, Plant & Equipment'.lower()])
                print(d_cash_year[ticker]['Acquisitions'.lower()])
                try:
                    d_cash_year[ticker]['Capital Expenditures'] = [sum(x) for x in zip(
                        d_cash_year[ticker]['Purchase Of Property, Plant & Equipment'.lower()],
                        d_cash_year[ticker]['Acquisitions'.lower()],
                        )]
                except KeyError:
                    pass

            if source == 'Yahoo':
                d_income_year[ticker]['Operating Income'] = d_income_year[ticker]['Operating Income or Loss']
                d_income_year[ticker]['Cost of Revenue, Total'] = d_income_year[ticker]['Cost of Revenue']
                d_cash_year[ticker]['Cash from Operating Activities'] = d_cash_year[ticker]['Total Cash Flow From Operating Activities']
                if d_balance_year[ticker]['Short Term Investments'] == ['-']*len(d_balance_year[ticker]['Short Term Investments']):
                    d_balance_year[ticker]['Short Term Investments'] = [0]*len(d_balance_year[ticker]['Short Term Investments'])
                d_cash_year[ticker]['Total Cash Dividends Paid'] = d_cash_year[ticker]['Dividends Paid']
                d_cash_year[ticker]['date'] = d_cash_year[ticker]['Period Ending']

                d_income_year[ticker]['Net Income Before Extra. Items'] = []
                for i in range(len(d_income_year[ticker]['Extraordinary Items'])):
                    if d_income_year[ticker]['Extraordinary Items'][i] not in ('-', '--',):
                        d_income_year[ticker]['Net Income Before Extra. Items'].append(
                            d_income_year[ticker]['Net Income'][i]-d_income_year[ticker]['Extraordinary Items'][i])
                    else:
                        d_income_year[ticker]['Net Income Before Extra. Items'].append(d_income_year[ticker]['Net Income'][i])

                if '-' not in d_balance_year[ticker]['Total Stockholder Equity']:
                    d_balance_year[ticker]['Total Equity'] = d_balance_year[ticker]['Total Stockholder Equity']
                else:
                    print(d_balance_year[ticker]['Total Assets'])
                    d_balance_year[ticker]['Total Equity'] = []
                    for i in range(len(d_balance_year[ticker]['Total Assets'])):
                        if d_balance_year[ticker]['Total Assets'][i] in ('-', '--',):
                            d_balance_year[ticker]['Total Equity'].append('-')
                        else:
                            d_balance_year[ticker]['Total Equity'].append(
                                d_balance_year[ticker]['Total Assets'][i]-d_balance_year[ticker]['Total Liabilities'][i])

            ##
            ## Parse dividend history, because missing from Yahoo CF stmt for non US-eq
            ## However also parsed from morningstar!
            ##
            url_dividend = 'http://real-chart.finance.yahoo.com/table.csv?s={}&g=v&ignore=.csv'.format(ticker_yahoo)
            d_dividends = parse_Yahoo.Yahoo().parse_historical_prices(url_dividend)
            l = []
            if source == 'Yahoo':
                for i, year in enumerate(d_cash_year[ticker]['date']):
                    try:
                        dividend = d_dividends[year]
                    except KeyError:
                        dividend = 0
                    l.append(dividend)
            if d_cash_year[ticker]['Total Cash Dividends Paid'] == ['-']*len(d_cash_year[ticker]['Total Cash Dividends Paid']):
                d_cash_year[ticker]['Total Cash Dividends Paid'] = l

            ##
            ## parse overview/snapshot
            ##
##            (
##                name,
##                currency_overview_symbol, price,
##                sector, industry,  statementNA,
##                mc_overview, beta,
##                ) = parse_BusinessWeek.BusinessWeek().parse_overview(
##                    ticker_businessweek,
##                    )
##            for suffix in ('','L'):
##                ticker_reuters = ticker_reuters.replace('.L','{}.L'.format(suffix))
            ## parse from FT instead... to avoid ticker conversion...
            ## http://markets.ft.com/research/Markets/Tearsheets/Summary?s=SYX:GER
            mc_overview = ''
            for suffix in ('', '.O', '.N'):
                url_reuters = url = 'http://www.reuters.com/finance/stocks/overview?symbol=%s' %(ticker_reuters+suffix)
                try:
                    (
                        name,
                        currency_overview_symbol, price,
                        sector, industry,  statementNA,
                        mc_overview, beta,
                        ) = parse_Reuters.Reuters().parse_overview(
                            ticker_reuters+suffix, url,
                            )
                    if mc_overview != '':
                        break
                except:
                    continue
            if mc_overview == '':
                statement_error = True

            if statement_error is True:
                self.stmtERROR(ticker, l_FAs=l_FAs, url='reuters')
                print('market cap not available')
                continue

            if any([
                sector in ['Bank', 'Insurance', 'Financials', 'financials'],
                'Cost of Revenue, Total' not in d_income_year[ticker].keys(),
                ]):
                if industry not in [
                    'Banks',
                    'REIT - Residential & Commercial',
                    'REITs - Residential',
                    'REITs - Commercial',  # BDN
                    'Real Estate Operations',
                    'Real Estate Development / Operations',
                    'Consumer Financial Services',
                    'Financial Services - Diversified',
                    'Financials - Specialty',
                    'Investment Trusts',
                    'Investment Banking / Brokerage Services',  # BGCP
                    'Investment Services',
                    'Investment Services - Specialty',
                    'Investment Management / Fund Operators',  # GLG
                    'Reinsurance',
                    'Insurance - Life & Health',
                    'Insurance - Life / Health',  # CRVL
                    'Insurance - Property & Casualty',
                    'Insurance - Property / Casualty',  # RSA.L
                    'Insurance - Multiline',
                    'Commercial Services & Supplies',  # SFI, Industrials
                    'Financial & Commodity Market Operators',  # FIS
                    'Corporate Financial Services',  # GMT
                    'Diversified Investment Services',  # ICE
                    'Real Estate Development & Operations',  # IRS
                    'Multiline Insurance & Brokers',  # MMC
                    'Investment Management & Fund Operators',  # PBY
                    'Specialized REITs',  # RLJ
                    'Commercial REITs',  # ROIC
                    'Closed End Funds',  # SLRC
                    'Property & Casualty Insurance',  # GBLI
                    'Investment Banking & Brokerage Services',  # GCAP
                    'Life & Health Insurance',  # SLF
                    'Holding Companies',  # 0001.HK
                    'Residential REITs',  # AMH
                    'Consumer Lending',  # HAWK
                    'Diversified REITs',  # AU:MGR
                    'Real Estate Services',
                    'undefined',
                    ]:
                    print(industry)
                    print(sector)
                    if sector in [
                        'Bank', 'Insurance', 'Financials', 'financials']:
                        stop_tmp
                fd = open('index_financials.txt', 'a')
                fd.write('%s\n' %(ticker))
                fd.close()
                print('financial company', industry, ticker)
                continue

            ##
            ## check statements for unwanted zero or negative values
            ##
            bool_skip = self.check_if_zero(
                ticker,
                d_income_year,
                d_balance_year,
                d_cash_year,
                div_yield_min,
                index0, index1,
                )
            if bool_skip is True:
                self.stmtERROR(ticker, l_FAs=l_FAs, url=url, suffix='fundamentals')
                print('skip skip something zero', ticker)
                continue

            ##
            ##
            ##
            ## Pence to pounds...
            if source in ('Yahoo', 'FT',):
                rate_stmt = d_currency_symbol[currency_stmt.upper()]
            else:
                rate_stmt = d_currency_name[currency_stmt]

            if currency_overview_symbol is not None:
                rate_overview = d_currency_symbol[currency_overview_symbol]
            else:
                ## assume currency is that given on statements
                rate_overview = rate_stmt

            print(type(mc_overview), mc_overview, rate_overview)
            mc_USD = mc_overview/rate_overview ## e.g. AZN (SEK or GBP) > AZN (USD)
            mc_stmt = mc_USD * rate_stmt ## e.g. AZN (USD) > AZN.L (GBP)

##            if round(mc_overview,0) != round(mc_stmt,0):
##                print ticker
##                print 'mc overview', mc_overview
##                print 'mc stmt', mc_stmt
##                print 'mc USD', mc_USD
##                print currency_stmt
##                print currency_overview_symbol
##                stop

            ##
            ## parse CEO compensation
            ##
            try:
                compensation, currency_compensation = parse_Reuters.Reuters().parse_CEOcompensation(
                    ticker_reuters,
                    )
                compensation /= d_currency_symbol[currency_compensation]
            except:
                print('compensation N/A', ticker)
                compensation = 'N/A'

            ##
            ## set first column of statements
            ##
            col1iy = 0
            col1by = 0
            col1cy = 0

            d_fundamentals, d_multiples, d_ratios = self.calc_fundamentals_and_multiples(
                ticker, d_income_year, d_balance_year, d_cash_year,
                mc_stmt,
                index0, index1,
                )

            if d_multiples['EV/FCF'] == 0:
                print('EV/FCF=0, so FCF probably negative...')
                self.stmtERROR(ticker, l_FAs=l_FAs, suffix='fundamentals')
                continue

            ## filtering, balance sheet
            if d_fundamentals['D/Eq'] > DE_max:
                if ticker not in l_FAs:
                    fd = open('statement_fundamentals.txt', 'a')
                    fd.write('%s\n' %(ticker))
                    fd.close()
                else:
                    fd = open('FAissues.txt', 'a')
                    fd.write('%s debt/equity %s\n' %(
                        ticker, d_fundamentals['D/Eq']))
                    fd.close()
                print('Debt/Eq.', d_fundamentals['D/Eq'])
                if ticker not in d_portfolio.keys():
                    print('2015-05-16a', ticker)
                    continue

            ## http://www.fool.com/investing/dividends-income/2005/12/29/foolish-fundamentals-free-cash-flow.aspx
            ## "Just make sure that the company pays out a portion of its free cash flow as a dividend but doesn't pay out more in dividends than it generates in free cash flow."
            ## http://msn.fool.com/investing/dividends-income/2007/11/02/make-millions-with-7-stocks.aspx
            ## "payouts less than 80% of free cash flow (more cushion means a company's better able to pay its dividend consistently)"
            if d_fundamentals['FCF'] == 0:
                payout_ratio_cash_flow_free = float('inf')
            else:
                payout_ratio_cash_flow_free = -d_fundamentals['cashflow_dividends']/d_fundamentals['FCF']
            if payout_ratio_cash_flow_free > payout_ratio_cash_flow_free_max:
##            if payout_ratio_cash_flow_free > 1.1 and ticker not in d_portfolio.keys():
                print('payout_ratio_cash_flow_free', payout_ratio_cash_flow_free)
                print('dividends', -d_fundamentals['cashflow_dividends']/1000000.)
                print('fcf', d_fundamentals['FCF']/1000000.)
                s = '%s payout ratio %.2f > %.2f too high relative to fcf (could be due to acquisation! see also Non-Cash items on cash flow stmt)\n' %(
                    ticker, payout_ratio_cash_flow_free, payout_ratio_cash_flow_free_max)
                if ticker not in l_FAs:
                    fd = open('statement_fundamentals.txt', 'a')
                    fd.write('%s\n' %(ticker))
                    fd.close()
                else:
                    fd = open('FAissues.txt', 'a')
                    fd.write(s)
                    fd.close()
                    pass

##            ##
##            ## parse financial 10 year summary
##            ##

            ## MSN Money stopped working.
##            url = 'http://moneycentral.msn.com/investor/invsub/results/statemnt.aspx?Symbol=%s&lstStatement=10YearSummary&stmtView=Ann' %(ticker_msn)
##            url = 'http://investing.money.msn.com/investments/financial-statements?symbol=%s' %(ticker_msn)
##            for k,v in list(d_ADR.items()):
##                if ticker in v:
##                    url = 'http://moneycentral.msn.com/investor/invsub/results/statemnt.aspx?Symbol=%s&lstStatement=10YearSummary&stmtView=Ann' %(k)
##                    url = 'http://investing.money.msn.com/investments/financial-statements?symbol=%s' %(k)
##                    ticker_msn = k
##                    break
##            ## not multiplied by appropriate factors (shares outstanding and assets/liabilities)
##            dic_10year = parse_MSN.MSN().parse_financial_10_year_summary(url,)

            ## Morningstar reports earnings 10 years back, but per share ratios for 5 years.
            ticker_morningstar = ticker_conversion.msn2morningstar(ticker)
            url = url_morningstar = 'http://financials.morningstar.com/financials/getFinancePart.html?&callback=?&t=%s' %(ticker_morningstar)
            try:
                dic_10year = parse_morningstar.morningstar().parseKeyRatios(
                    url_morningstar, ticker_morningstar)
                assert len(dic_10year.keys()) > 0
            except:
                fd = open('statement_pending.txt', 'a')
                fd.write('%s\n' %(ticker))
                fd.close()
                l_statement_pending += [ticker]
                print(url_morningstar)
                dic_10year = parse_morningstar.morningstar().parseKeyRatios(
                    url_morningstar, ticker_morningstar)
                print('morningstar', dic_10year)
                if ':' in ticker or '.' in ticker and (ticker == ticker_morningstar or ticker_yahoo == ticker_morningstar):
                    print(ticker, ticker_morningstar)
                    print(url_morningstar)
                continue

####            url_gurufocus = 'http://www.gurufocus.com/financials/{}'.format(ticker)
####            d_gf = self.parse_gurufocus(url_gurufocus)
####            if not d_gf: continue # tmp!
##
##            l_equity_log = [math.log(x) for x in d_gf['Total Equity']]
##            l_sales_log = [math.log(x) for x in d_gf['Revenue']]
##            r_equity_log = statistics.correlation(list(range(len(l_equity_log),1-1,-1)),l_equity_log)
##            r_sales_log = statistics.correlation(list(range(len(l_sales_log),1-1,-1)),l_sales_log)
##            
##            print(r_equity_log)
##            print(r_sales_log)
##            stop

            bool_continue = False
            for k in ('Book Value Per Share', 'Earnings Per Share',):
                if '-' in dic_10year['Book Value Per Share']:
                    print(k, dic_10year[k])
                    self.stmtPEND(ticker)
                    bool_continue = True
                    break
                if min([x for x in dic_10year[k] if x != '-']) <= 0:
                    print('negative', k, dic_10year[k])
                    self.stmtERROR(
                        ticker, l_FAs=l_FAs, url=url, suffix='fundamentals')
            if bool_continue:
                continue

            if len(dic_10year['SHARES OUTSTANDING']) == 0:
                r_sales_log = 0
                if ticker in l_FAs:
                    l_statement_missing += [ticker]
                    print('MSN 10y statement missing, zero shares')
                    r_equity_log = 0
                else:
                    print('MSN 10y statement missing, zero shares')
                    r_equity_log = 0

            elif len(dic_10year['SHARES OUTSTANDING']) > 0:
                if any([
                    len(dic_10year['SHARES OUTSTANDING']) < 9,
                    ## Could be a non-parsed dash...
                    len(dic_10year['Book Value Per Share']) < 9,
                    ]):
                    if ticker not in l_FAs:
                        if len(dic_10year['SHARES OUTSTANDING']) == 9:
                            fd = open('statement_pending.txt', 'a')
                            fd.write('%s\n' %(ticker))
                            fd.close()
                            l_statement_pending += [ticker]
                        else:
                            fd = open('statement_fundamentals.txt', 'a')
                            fd.write('%s\n' %(ticker))
                            fd.close()
                        print('company has only been around for', len(dic_10year['SHARES OUTSTANDING']), 'years')
                        print(dic_10year['SHARES OUTSTANDING'])
##                        print(dic_10year['CURRENT ASSETS'])
                        continue
                    else:
##                        while len(dic_10year['SHARES OUTSTANDING']) < 10:
##                            dic_10year['SHARES OUTSTANDING'].append(0)
                        l_statement_pending += [ticker]
                        print('len', 'Book Value Per Share', len(dic_10year['Book Value Per Share']))
                        print(
                            'SHARES OUTSTANDING',
                            len(dic_10year['SHARES OUTSTANDING']), dic_10year['SHARES OUTSTANDING'])
                        print(url)
                        print('2016-03-29a', ticker)
##                        if ticker in ('ULVR.L', 'GB:ULVR'):
##                            stop5
                        continue
                l_equity = []
                l_sales = []
                l_NetIncome = []
                l_EPS = []
##                l_EBIT = []
                print(dic_10year.keys())
                print(len(dic_10year['Book Value Per Share']), dic_10year['Book Value Per Share'])
                dic_10year['RETURN ON EQUITY (%)'] = []
                dic_10year['NET PROFIT MARGIN (%)'] = []
                ## Do a loop to calculate equity as assets-liabilities and equity per share
                for i in range(max(9, len(dic_10year['SHARES OUTSTANDING']))):
##                    assets = dic_10year['CURRENT ASSETS'][i]
##                    liabilities = dic_10year['CURRENT LIABILITIES'][i]
##                    equity = assets-liabilities
                    equity = dic_10year['Book Value Per Share'][i]
                    shares = dic_10year['SHARES OUTSTANDING'][i]
##                    if assets == 0 or shares == 0:
##                        break
                    if shares == 0:
                        break
##                    l_equity += [equity/shares]
                    l_equity += [equity]
                    sales = dic_10year['SALES'][i]
##                    EBIT = dic_10year['EBIT'][i]
                    l_sales += [sales/shares]
                    l_NetIncome += [dic_10year['Net Income'][i]]
                    l_EPS += [dic_10year['Earnings Per Share'][i]]
##                    l_EBIT += [EBIT/shares]
                    if dic_10year['Book Value Per Share'][i] != '-' and dic_10year['Earnings Per Share'][i] != '-':
                        roe = dic_10year['Earnings Per Share'][i]/dic_10year['Book Value Per Share'][i]
                    else:
                        roe = '-'
                    dic_10year['RETURN ON EQUITY (%)'].append(roe)
                    if dic_10year['Net Income'][i] != '-':
                        dic_10year['NET PROFIT MARGIN (%)'].append(dic_10year['Net Income'][i]/dic_10year['Revenue'][i])

                ## add cash that was used for *recent* share repurchase to shareholders' equity (FDO 2011)
                ## assume same currency at msn and bw...
                if len(l_equity) > 0 and 'Repurchase of Common Stock' in d_cash_year[ticker].keys():
                    l_equity[0] += -(
                        d_cash_year[ticker]['Repurchase of Common Stock'][index0]
                        ) / (1000000*dic_10year['SHARES OUTSTANDING'][index0])

                    bw_equity = d_balance_year[ticker]['Total Equity'][index0]
                    msn_equity = 1000000*l_equity[0]*dic_10year['SHARES OUTSTANDING'][index0]
                    if any([
                        (bw_equity-msn_equity)/bw_equity > 0.1,
                        (bw_equity-msn_equity)/msn_equity > 0.1,
                        ]):
                        print('equity', l_equity)
                        print('bw', bw_equity)
                        print('msn', msn_equity)
                        print((bw_equity-msn_equity)/bw_equity)
                        print((bw_equity-msn_equity)/msn_equity)
                        print(dic_10year['SHARES OUTSTANDING'][index0])
##                        stop

                if min(dic_10year['SHARES OUTSTANDING']) == 0:
                    if i != 0:
                        print('zero shares, assets and/or liabilities. listed less than 10 years?', assets, shares)
                    else:
                        print('zero shares, assets and/or liabilities. numbers not updated or recent listing.')
                    if ticker not in l_FAs:
                        fd = open('statement_pending.txt', 'a')
                        fd.write('%s\n' %(ticker))
                        fd.close()
                    l_statement_pending += [ticker]
                    print('2016-03-29b', ticker)
                    continue
                elif any([
                    min((x for x in l_equity if x != '-')) < 0,
                    min((x for x in dic_10year['Operating Income'] if x != '-')) < 0,
                    min((x for x in dic_10year['Operating Cash Flow'] if x != '-')) < 0,
                    ]):
                    print('Operating Income', dic_10year['Operating Income'])
                    print('Operating Cash FLow', dic_10year['Operating Cash Flow'])
                    print('equity', l_equity)
                    print('negative equity or earnings or OCF within the past 10 years', l_equity)
                    ## It's OK for a startup such as Illumina and PacBio
                    ## to have had negative operating income and negative equity!
                    if len(tickers) > 5:
                        print('2016-03-29c', ticker)
                        if ticker in d_portfolio.keys():
                            pass
                        else:
                            if ticker in l_FAs: ## tmp!!!
                                fd = open('FAissues.txt', 'a')
                                fd.write('%s equity %s operating income %s\n' %(
                                    ticker, str(l_equity),
                                    str(dic_10year['Operating Income'])))
                                fd.close()
                            else:
                                fd = open('statement_fundamentals.txt', 'a')
                                fd.write('%s\n' %(ticker))
                                fd.close()
                            continue
                elif min(l_sales) == 0:
                    print('stmt pending, sales', l_sales)
                    l_statement_pending += [ticker]
                    fd = open('statement_pending.txt', 'a')
                    fd.write('%s\n' %(ticker))
                    fd.close()
                    print('2016-03-29d', ticker)
                    continue
                else:
                    ## the correlation is independent of the base of the log
                    l_equity_log = [math.log(l_equity[year]) for year in range(len(l_equity)) if l_equity[year] != '-']
                    l_sales_log = [math.log(l_sales[year]) for year in range(len(l_sales))]

                    if '-' in l_equity_log[:5] or len(l_equity_log) < 5:
                        continue

                    try:
                        l_income_log = [math.log(l_NetIncome[year]) for year in range(len(l_NetIncome))]
                        l_eps_log = [math.log(l_EPS[year]) for year in range(len(l_EPS))]
                        r_income_log = statistics.correlation(list(range(len(l_income_log),1-1,-1)),l_income_log)
                        r_eps_log = statistics.correlation(list(range(len(l_eps_log),1-1,-1)),l_eps_log)
                    except ValueError:
                        r_income_log = 1
                        r_eps_log = 1

                    r_equity_log = statistics.correlation(list(range(len(l_equity_log),1-1,-1)),l_equity_log)
                    r_sales_log = statistics.correlation(list(range(len(l_sales_log),1-1,-1)),l_sales_log)

##                    l_equity_sorted = [x for x in l_equity if x != '-']
##                    l_equity_sorted.sort()
##                    l_equity_sorted.reverse()
                    l_sales_sorted = list(l_sales)
                    l_sales_sorted.sort()
                    l_sales_sorted.reverse()

##                    ## do exponential regression instead of linear regression...
##                    n = len(l)
##                    l2 = [ [ year+1, l[10-1-year], ] for year in range(10) ]
##                    sumXX, sumXY, sumX, sumY = 0.0, 0.0, 0.0, 0.0
##                    print l2
##                    for x, y in l2:
##                        sumXY += x * y
##                        sumXX += x * x
##                        sumX += x
##                        sumY += y
##                    a = (n * sumXY - sumX * sumY) / (n * sumXX - sumX * sumX)
##                    b = (sumY - a * sumX) / n
##                    print a,b
##                    stop

                    if any([
##                        ## Buffett bought IBM despite r_equity_log**2 being 0.18
                        all([
                            r_equity_log**2 < r2_min_equity,
                            ## if steady equity growth then OK it was lower one year
                            list(reversed(sorted([
                                x for x in l_equity if x != '-']))) != l_equity,
                            ]),
                        all([
                            r_sales_log**2 < r2_min_sales,
                            r_income_log**2 < 0.9,
                            r_eps_log**2 < 0.9,
                            l_sales_sorted != l_sales,
                            ]),
                        ]):
                        print('r_equity, r_sales, net_inc, eps', r_equity_log, r_sales_log, r_income_log, r_eps_log)
                        if ticker not in l_FAs:
                            fd = open('statement_fundamentals.txt', 'a')
                            fd.write('%s\n' %(ticker))
                            fd.close()
                        else:
                            fd = open('FAissues.txt', 'a')
                            fd.write('%s r**2 (10y) equity %s sales %s\n' %(
                                ticker,round(r_equity_log**2,2),round(r_sales_log**2,2),))
                            fd.close()
                        if ticker not in d_portfolio.keys():
                            print('r_equity_log**2', r_equity_log**2)
                            print('r_sales_log**2', r_sales_log**2)
                            print('l_sales', '\n'.join(str(_) for _ in l_sales))
                            print('l_equity', '\n'.join(str(_) for _ in l_equity))
                            continue

                    ##
                    ## 5 year correlations
                    ##
                    r5_equity = statistics.correlation(list(range(5,1-1,-1)),l_equity_log[:5])
                    r5_sales = statistics.correlation(list(range(5,1-1,-1)),l_sales_log[:5])
                    r5_eps = statistics.correlation(list(range(5,1-1,-1)),l_eps_log[:5])
##                    if r5_equity < .5 or r5_sales < .5:
##                    if r5_sales < .5:
                    if r5_equity < .5 and r5_sales < .5 and r5_eps < .5:
                        print('r_equity5, r_sales5', r5_equity, r5_sales)
                        if ticker not in l_FAs:
                            fd = open('statement_fundamentals.txt', 'a')
                            fd.write('%s\n' %(ticker))
                            fd.close()
                        else:
                            fd = open('FAissues.txt', 'a')
                            fd.write('%s r5 %s %s\n' %(ticker,r5_equity,r5_sales,))
                            fd.close()
                        if ticker not in d_portfolio.keys():
                            print('r5', r5_equity, r5_sales)
                            continue

                    a10_sales,b10_sales,RR10_sales = self.linreg(list(range(max(9,len(l_sales_log)),0,-1)), l_sales_log[:10],)
                    a5_sales,b5_sales,RR5_sales = self.linreg(list(range(5,0,-1)), l_sales_log[:5],)
                    ## convert a in b*e^(a*x) to annual percentage
                    a10_sales = round(100*math.exp(a10_sales),1)-100
                    a5_sales = round(100*math.exp(a5_sales),1)-100
                    if a10_sales < 5 or a5_sales < 3.5:
                        print('sales slope less than 5% (10y) or less than 3.5% (5y)', a10_sales, a5_sales, l_sales)
                        if ticker in l_FAs:
                            fd = open('FAissues.txt', 'a')
                            fd.write(
                                '%s 10ysales %.1f 5ysales %.1f\n' %(
                                    ticker,a10_sales,a5_sales,))
                            fd.close()
                        else:
                            fd = open('statement_fundamentals.txt', 'a')
                            fd.write('%s\n' %(ticker))
                            fd.close()
                            continue
                        
            else:
                print('MSN 10y statement missing', len(dic_10year['SHARES OUTSTANDING']))
                r_equity_log = 0
                r_sales_log = 0

##            try:
##                dic_10year = parse_MSN.MSN().key_ratios_10_year_summary(ticker_msn, dic_10year,)
##            except:
##                l_statement_missing += [ticker]
##                continue

            if 'RETURN ON EQUITY (%)' not in dic_10year.keys():

                if ticker in l_FAs:
                    l_statement_missing += [ticker]
                    print('MSN 10y statement missing, ROE missing', ticker_msn)
                else:
                    print('MSN 10y statement missing, ROE missing', ticker_msn)

            else:

##                for k in ['Return on Equity (%)','Net Profit Margin (%)']:
                for k in ['RETURN ON EQUITY (%)', 'NET PROFIT MARGIN (%)']:
                    l = list(reversed(dic_10year[k]))  # oldest at top
                    l = dic_10year[k]  # new to old (morningstar)

                    bool_decrease = False

                    if len(l) < 10:
                        if ticker not in l_FAs:
                            fd = open('statement_pending.txt', 'a')
                            fd.write('%s\n' %(ticker))
                            fd.close()
                        l_statement_pending += [ticker]
                        break

                    if min(l) < 0 or 'NA' in l:
                        bool_decrease is True
                        break

                    ## decreasing every year over a 5 year period?
                    bool_decrease = True
                    for i in range(5):
                        if l[i] > l[i+1]:
                            bool_decrease = False
                            break
                    if bool_decrease is True:
                        break

                    ## decreasing on average more than a treshold over a 5 and 10 year period?
                    a10,b10,RR10 = self.linreg(list(range(10,0,-1)), l[:10],)
                    a5,b5,RR5 = self.linreg(list(range(5,0,-1)), l[:5],)
                    ## you don't want a big recent drop
                    ## you don't want a big steady drop over many years
##                    if a5 < -b5/20. or a10 < -b10/20.:
                    if (a5 < -b5/20. and a10 < 0) or (a10 < -b10/30. and a5 < 0):
                        print(k, l, 'slope', '5y', a5, '10y', a10)
                        if ticker not in l_FAs:
                            bool_decrease = True
                            break

                if bool_decrease is True:
                    print('decrease', ticker, k, l)
                    print(
                        'BVPS', dic_10year['Book Value Per Share'],
                        '\n',
                        'EPS', dic_10year['Earnings Per Share'],
                        )
                    if ticker not in l_FAs:
                        fd = open('statement_fundamentals.txt', 'a')
                        fd.write('%s\n' %(ticker))
                        fd.close()
                    else:
                        fd = open('FAissues.txt', 'a')
                        fd.write('%s %s 5yslope %.2f 10yslope %.2f, b5=%s, b10=%s\n' %(
                            ticker, k, a5, a10, b5, b10))
                        fd.close()
                    continue

            ##
            ## statement pending
            ##
            if source in ('Yahoo', 'FT',):
                if max(d_cash_year[ticker]['date']) <= 2012:
                    print('check if this company is still reporting')
                    continue
            else:
                if 'date' not in d_cash_year[ticker].keys():
                    print(d_cash_year[ticker])
                    if source == 'ADVFN':
                        ## not entirely correct...
                        d_cash_year[ticker]['date'] = d_indicators['date preliminary data loaded']
                else:
                    if dic_10year['&nbsp;'] == []:
        ##                min_year = int(max(d_cash_year[ticker]['date'])[:4])
                        min_year = int(max(d_cash_year[ticker]['date']))
                    else:
        ##                min_year = min(int(max(d_cash_year[ticker]['date'])[:4]),max(dic_10year_balance['&nbsp;']))
                        min_year = min(int(max(d_cash_year[ticker]['date'])),max(dic_10year['&nbsp;']))

                    if min_year <= 2009:
                        print('check if this company is still reporting')
                        continue

                    ## this check only applicable if not parsed from advfn
                    if d_income_year[ticker]['date'] != d_balance_year[ticker]['date']:
                        if ticker not in l_FAs:
                            fd = open('statement_pending.txt', 'a')
                            fd.write('%s\n' %(ticker))
                            fd.close()
                        l_statement_pending += [ticker]
                        print('waiting for income statement or balance sheet')
                        continue
                    if d_cash_year[ticker]['date'][:4] != d_income_year[ticker]['date'][:4]:
                        if ticker not in l_FAs:
                            fd = open('statement_pending.txt', 'a')
                            fd.write('%s\n' %(ticker))
                            fd.close()
                        l_statement_pending += [ticker]
                        print('waiting for cashflow statement')
                        continue

            ##
            ## parse financial highlights
            ##
            payout_ratio = abs(d_fundamentals['cashflow_dividends'])/d_income_year[ticker]['Net Income Before Extra. Items'][index0]

            if payout_ratio == 'N/A':
                print('payout N/A')
                continue

            if payout_ratio > payout_ratio_max:
                if len(tickers) > 5:
                    print('divi', abs(d_fundamentals['cashflow_dividends']))
                    print('income', d_income_year[ticker]['Net Income Before Extra. Items'][index0])
                    if ticker not in l_FAs:
                        fd = open('statement_fundamentals.txt', 'a')
                        fd.write('%s\n' %(ticker))
                        fd.close()
                    else:
                        fd = open('FAissues.txt', 'a')
                        fd.write('%s payot (%s) > payout max\n' %(ticker,payout_ratio))
                        fd.close()
                    print('payout (%.2f) > payout max (%s)' %(payout_ratio, payout_ratio_max,))
                    if ticker not in d_portfolio.keys():
                        continue

            ## zero or negative free/operating cash flow
            if d_fundamentals['OCF'] <= 0:
                if ticker not in l_FAs:
                    fd = open('statement_fundamentals.txt', 'a')
                    fd.write('%s\n' %(ticker))
                    fd.close()
                else:
                    fd = open('FAissues.txt', 'a')
                    fd.write('%s negative fcf\n' %(ticker))
                    fd.close()
                print('zero or negatgive ocf')
                if ticker not in d_portfolio.keys():
                    continue

            x1 = d_cash_year[ticker]['Cash from Operating Activities'][-1]
            x2 = abs(d_cash_year[ticker]['Capital Expenditures'][-1])
            if x1 < x2:
                fd = open('statement_fundamentals.txt', 'a')
                fd.write('%s\n' %(ticker))
                fd.close()
                print('operating cash flow doesnt cover CapEx')
                continue

            if div_yield_min > 0 and not 'Total Cash Dividends Paid' in d_cash_year[ticker].keys():
                fd = open('statement_fundamentals.txt', 'a')
                fd.write('%s\n' %(ticker))
                fd.close()
                print('zero dividend')
                continue

##            ## "ROE is presumably irrelevant if the earnings are not reinvested."
##            ## "The sustainable growth model shows us that when firms pay dividends, earnings growth lowers. If the dividend payout is 20%, the growth expected will be only 80% of the ROE rate."
##            ## http://en.wikipedia.org/wiki/Return_on_equity
##            ## measure growth in sales and earnings instead... large reinvestments might not be necessary if net tangible assets is a low number (E/NTA high)
##            if d_fundamentals['roe_curr']*(1-payout_ratio) < .1:
##                print 'high ROE but money not reinvested but paid out'
##                print 'roe', d_fundamentals['roe_curr'], 'payout', payout_ratio, 'roe adj', d_fundamentals['roe_curr']*(1-payout_ratio)
##                if len(tickers) > 5:
##                    if ticker not in l_FAs:
##                        fd = open('statement_fundamentals.txt', 'a')
##                        fd.write('%s,' %(ticker))
##                        fd.close()
##                    else:
##                        fd = open('FAissues.txt','a')
##                        fd.write(
##                            '%s high ROE but money not reinvested but paid out, %s, income curr=%s prev=%s, equity curr=%s prev=%s\n' %(
##                                ticker, d_fundamentals['roe_curr'],
##                                d_income_year[ticker][s_income][0]/1000000.,
##                                d_income_year[ticker][s_income][1]/1000000.,
##                                d_balance_year[ticker]['Total Equity'][0]/1000000.,
##                                d_balance_year[ticker]['Total Equity'][1]/1000000.,
##                                )
##                            )
##                        fd.close()
##                    if ticker not in d_portfolio.keys():
##                        continue

            ##
            ## fundamentals
            ##
            FA = True
            d_fundamentals_limits = {
                'ROA':[d_fundamentals['ROA'],ROA_min],
                'ROE':[d_fundamentals['ROE'],ROE_min],
                'NPM':[d_fundamentals['NPM'],NPM_min],
                }
            for key in d_fundamentals_limits.keys():
                if d_fundamentals_limits[key][index0] == 'N/A':
                    stop
                if d_fundamentals_limits[key][0] < d_fundamentals_limits[key][1]:
                    print(key, d_fundamentals[key])
                    if ticker not in l_FAs:
                        fd = open('statement_fundamentals.txt', 'a')
                        fd.write('%s\n' %(ticker))
                        fd.close()
                    else:
                        fd = open('FAissues.txt', 'a')
                        fd.write('%s %s %s\n' %(
                            ticker, key, round(d_fundamentals_limits[key][0],2),
                            ))
                        fd.close()
                    FA = False
                    print('fundamental', key)
                    break
            if FA is False:
                continue

            ##
            ## multi year fundamentals
            ##
            bool_FA = True
            for year in range(len(d_ratios['profitability']['ROE'])):
                if d_ratios['profitability']['ROE'][year] < ROE_min:
                    print('year', year, 'ROE', d_ratios['profitability']['ROE'][year], '< ROE_min', ROE_min)
                    if ticker not in l_FAs:
                        fd = open('statement_fundamentals.txt', 'a')
                        fd.write('%s\n' %(ticker))
                        fd.close()
                    else:
                        fd = open('FAissues.txt', 'a')
                        fd.write('%s ROE %s year %s %s\n' %(ticker, 'ROE', year, d_ratios['profitability']['ROE'][year],))
                        fd.close()
                    bool_FA = False
                    break
            if bool_FA is False:
                continue

            ##
            ## multiples
            ##
            FA = True
            d_multiples_limits = {
                'P/E':[d_multiples['P/E'],PE_max],
                'P/CF':[d_multiples['P/CF'],PCF_max],
                'EV/FCF':[d_multiples['EV/FCF'],EVFCF_max],
                'EV/OCF':[d_multiples['EV/OCF'],EVOCF_max],
                }
            for key in d_multiples_limits.keys():
                if d_multiples_limits[key][0] > d_multiples_limits[key][1]:
                    l_multiples += [ticker]
                    print(key, d_multiples_limits[key][0])
                    print('multiple', key)
##                    if ticker not in d_portfolio.keys():
##                        if ticker not in l_FAs:
##                            fd = open('multiples.txt', 'a')
##                            fd.write('%s,' %(ticker))
##                            fd.close()
##                            FA = False
##                        else:
##                            fd = open('FAissues.txt','a')
##                            fd.write('%s %s multiple %s\n' %(ticker,key, d_multiples[key],))
##                            fd.close()
##                        break
            if FA is False and ticker not in d_portfolio.keys():
                continue

            ##
            ## multiples
            ##
            if d_multiples['div_yield'] < div_yield_min:
                if ticker not in d_portfolio.keys():
                    if ticker not in l_FAs:
                        fd = open('multiples.txt', 'a')
                        fd.write('%s\n' %(ticker))
                        fd.close()
                        print(
                            'div yield < minimum div yield',
                            d_multiples['div_yield'], '<', div_yield_min)
                        continue
                    else:
                        fd = open('FAissues.txt', 'a')
                        fd.write('%s %s div yield\n' %(
                            ticker, d_multiples['div_yield']))
                        fd.close()

            if all([
                'Total Cash Dividends Paid' in d_cash_year[ticker].keys(),
                abs(d_cash_year[ticker]['Total Cash Dividends Paid'][col1cy]) > 0,
                compensation != 'N/A',
                ]):
                x2 = 1.2*abs(d_cash_year[ticker]['Total Cash Dividends Paid'][col1cy])
                if compensation > x2:
                    print('compensaiton', compensation)
                    print('dividend   ', d_cash_year[ticker]['Total Cash Dividends Paid'][col1cy])
                    print('%s executive compensation is %s%% of annual dividends' %(
                        ticker,
                        100*compensation/abs(d_cash_year[ticker]['Total Cash Dividends Paid'][col1cy]),
                        ))
                    if ticker not in l_FAs:
                        fd = open('statement_fundamentals.txt', 'a')
                        fd.write('%s\n' %(ticker))
                        fd.close()
                    else:
                        fd = open('FAissues.txt', 'a')
                        fd.write('%s executive compensation %s\n' %(
                            ticker,
                            100*compensation/abs(d_cash_year[ticker]['Total Cash Dividends Paid'][col1cy]),
                            ))
                        fd.close()
                    continue

##            ## market capitalization calculation (problem when 1ADR!=1share)
##            mc2 = price*(d_balance_quarter[ticker]['Total Common Shares Outstanding'][col1b]+d_balance_quarter[ticker]['Total Preferred Shares Outstanding'][col1b])
##            if price != 0 and abs(mc-mc2)/mc > 0.07:
##                print price
##                print mc
##                print mc2
##                print abs(mc-mc2)/mc
##                print 'com', d_balance_quarter[ticker]['Total Common Shares Outstanding'][col1b]
##                print 'prf', d_balance_quarter[ticker]['Total Preferred Shares Outstanding'][col1b]
##                print 'rate', rate
##                stop_mc_diff

            ##
            ## check dividend growth (not decline)
            ##

            ## can't take into account extraordinary payments
            ## http://finance.yahoo.com/q/hp?s=COLO-B.CO&a=00&b=1&c=2003&d=07&e=31&f=2011&g=v
            ## Yahoo doesn't take into account splits for foreign stocks
            ## http://finance.yahoo.com/q/hp?s=HM-B.ST&a=00&b=1&c=2003&d=07&e=31&f=2011&g=v

            bool_dividend_not_decreasing = True

            ticker_yahoo = ticker_conversion.msn2yahoo(ticker)
            if '.' not in ticker_yahoo:
                bool_dividend_not_decreasing = self.check_dividend_historic(ticker_yahoo, d_ADR, time,)

            ## cant be used if share buyback
##            for year in range(len(d_cash_year[ticker]['Total Cash Dividends Paid'])-1):
##                if abs(d_cash_year[ticker]['Total Cash Dividends Paid'][year]) <= abs(0.99*d_cash_year[ticker]['Total Cash Dividends Paid'][year+1]):
##                    bool_dividend_not_decreasing = False
##                    break

            if bool_dividend_not_decreasing is False and ticker not in l_FAs:
                print(ticker, 'dividend not increasing')
                continue

            if d_multiples['div_yield'] == 0 and d_fundamentals['D/Eq'] >= 0.1:
                print(ticker, 'stop_zero_dividend_but_debt')
                continue
##                stop_zero_dividend_but_debt

################################################################################


            print('*****FA CANDIDATE*****')
            if ticker not in l_FAs:
                fd = open('FAcandidates.txt', 'a')
                fd.write('%s\n' %(ticker))
                fd.close()
            FAcandidates.append(ticker)

##            d_msn = parse_MSN.MSN().parse_company_report(ticker,rate,)
            d_msn = {
                'ma50': 'N/A', 'ma200': 'N/A', 'relative strength': 'N/A',}

            if compensation != 'N/A':
                compensation = round(compensation/1000000.,1)

            d_FAdata[ticker] = {
                'Market Cap (Mil)': mc_USD, 'name': name,
                'P/E': d_multiples['P/E'],
                'Div. Yield': d_multiples['div_yield']*100,
                'NPM': d_fundamentals['NPM']*100,
                'ROA': d_fundamentals['ROA']*100,
                'ROE': d_fundamentals['ROE']*100,
                'E/NTA': d_fundamentals['E/NTA']*100,
                'grow5y_sales': 0, 'grow5y_income': 0,
                'MA200': d_msn['ma200'], 'MA50': d_msn['ma50'],
                'price': price,
                'Payout Ratio': payout_ratio*100,
                'P/CF': d_multiples['P/CF'],
                'EV/FCF': d_multiples['EV/FCF'],
                'EV/OCF': d_multiples['EV/OCF'],
                'ATO': d_fundamentals['ATO'],
                'Equity Multiplier': d_fundamentals['equity_multiplier'],
                'compensation': compensation,
                'sector': sector, 'industry': industry,
                'r_sales_log**2': round(r_sales_log**2,2),
                'min_year': max(d_cash_year[ticker]['date']),
                'ticker_reuters': ticker_reuters,
                'relative strength': d_msn['relative strength'],
                }

            for key in ('NPM', 'ROE', 'ROA', 'ATO'):
                d_FAdata[ticker][key] = d_fundamentals[key]*100

            d_FAdata[ticker]['Equity Multiplier'] = d_fundamentals['equity_multiplier']

            d_FAdata[ticker]['D/Eq'] = d_fundamentals['D/Eq']
            
    ################################################################################

        print('not candidates anymore', set(l_FAs)-set(FAcandidates)-set(l_statement_pending)-set(l_statement_missing))
        print('new candidates', set(FAcandidates)-set(l_FAs))
        yahoo = 'http://finance.yahoo.com/q/cq?d=v1&s='
        for ticker in FAcandidates:
            ticker_yahoo = ticker_conversion.msn2yahoo(ticker)
            yahoo += '%s+' %(ticker_yahoo)

        print()
        for ticker in d_portfolio.keys():
            if not ticker in l_statement_pending+FAcandidates:
                print('SELL!!!', ticker)

        return FAcandidates, d_FAdata, yahoo, l_statement_pending, l_statement_missing, l_multiples


    def linreg(self, X, Y):
        """
        Summary
            Linear regression of y = ax + b
        Usage
            real, real, real = linreg(list, list)
        Returns coefficients to the regression line "y=ax+b" from x[] and y[], and R^2 Value
        """
        if len(X) != len(Y):
            raise ValueError('unequal length')
        N = len(X)
        Sx = Sy = Sxx = Syy = Sxy = 0.0
##        for x, y in map(None, X, Y):
        for x, y in zip(X, Y):
            Sx = Sx + x
            Sy = Sy + y
            Sxx = Sxx + x*x
            Syy = Syy + y*y
            Sxy = Sxy + x*y
        det = Sxx * N - Sx * Sx
        a, b = (Sxy * N - Sy * Sx)/det, (Sxx * Sy - Sx * Sxy)/det
        meanerror = residual = 0.0
##        for x, y in map(None, X, Y):
        for x, y in zip(X, Y):
            meanerror = meanerror + (y - Sy/N)**2
            residual = residual + (y - a * x - b)**2
        RR = 1 - residual/meanerror
        ss = residual / (N-2)
        Var_a, Var_b = ss * N / det, ss * Sxx / det
        return a, b, RR


    def stmtERROR(self, ticker, l_FAs=[], url='', suffix='error'):

        if ticker not in l_FAs:
            fd = open('statement_{}.txt'.format(suffix), 'a')
            fd.write('%s\n' %(ticker))
            fd.close()
        else:
            fd = open('FAissues.txt', 'a')
            fd.write('%s %s\n' %(ticker, url,))
            fd.close()

        return


    def stmtPEND(self, ticker):

        print('pending', ticker)
        fd = open('statement_pending.txt', 'a')
        fd.write('{}\n'.format(ticker))
        fd.close()

        return


    def check_dividend_historic(self, ticker, d_ADR, time,):

        year2 = time[0]
        month2 = time[1]
        day2 = time[2]
        year1 = year2-30
        month1 = month2
        day1 = day2

        url = 'http://ichart.finance.yahoo.com/table.csv?s=%s&a=%s&b=%s&c=%s&d=%s&e=%s&f=%s&g=v&ignore=.csv' %(
                ticker, month1-1, day1, year1, month2-1, day2, year2,
                )
        lines_dividend = screener.finance().read_url(url, ticker)
        if not lines_dividend:
            url = 'http://ichart.finance.yahoo.com/table.csv?s=%s&a=%s&b=%s&c=%s&d=%s&e=%s&f=%s&g=v&ignore=.csv' %(
                ticker_ADR, month1-1, day1, year1, month2-1, day2, year2,
                )
            lines_dividend = screener.finance().read_url(url, ticker)

        if not lines_dividend:
            bool_dividend_not_decreasing = True
            return bool_dividend_not_decreasing

        bool_dividend_not_decreasing = True

        if len(lines_dividend) == 1:
            return bool_dividend_not_decreasing

        d_dividends = {}
        l_months = []
        for i_line in range(1,len(lines_dividend)):

            line = str(lines_dividend[i_line])
            if 'Split' in line:
                continue

            try:
                year = int(line[:4])
                month = int(line[5:7])
                day = int(line[8:10])
            except:
                print(lines_dividend)
                print(url)
                stop
            dividend = float(line.strip().split(',')[1])
            if month == 12 and 1 in l_months:
                print(line)
                print(l_months)
                print(d_dividends)
                print(line)
                print(year, month)
                year += 1
                month = 1
            elif month == 1 and 12 in l_months:
                year -= 1
            elif month not in l_months and month-1 in l_months and day in range(1,11+1,):  # e.g. DGX, GD
                month -= 1
            elif month not in l_months and month+1 in l_months and day in range(27,31+1):  # e.g. CPSI, GD
                month += 1
            else:
                l_months += [month]

            if not year in d_dividends.keys():
                d_dividends[year] = {}

            if not month in d_dividends[year].keys():
                d_dividends[year][month] = 0
            ## extra dividend or Yahoo error
            elif dividend == d_dividends[year][month]: ## e.g. CHRW error
                continue
            ## extra dividend (prev, curr = next)
            elif dividend == float(lines_dividend[i_line+1].strip().split(',')[1]): ## e.g. 2002 HNZ extra dividend
                d_dividends[year][month] = float(line.strip().split(',')[1])
                continue
            ## extra dividend (curr, prev = next) ## e.g. 2010 RAVN extra dividend
            elif float(lines_dividend[i_line-1].strip().split(',')[1]) == float(lines_dividend[i_line+1].strip().split(',')[1]):
                d_dividends[year][month] = float(line.strip().split(',')[1])
                continue
            ## who cares...
            elif year < year2-10:
                continue
            else:
                d_dividends[year][month] = float(line.strip().split(',')[1])
##            elif not '.' in ticker and not ':' in ticker:
##                print(lines_dividend)
##                print(line)
##                print(year, month)
##                print(ticker)
##                print(d_dividends[year][month])
##                print(d_dividends[year])
            d_dividends[year][month] += float(line.strip().split(',')[1])

        if d_dividends == {}:
            print(lines_dividend)
            stop

        for year in range(
            max(
                ## go back 10 years
                max(d_dividends.keys())-10,
                min(d_dividends.keys())+1,
                ),
            max(d_dividends.keys())-1,
            ):
            if year+1 not in d_dividends.keys():
                bool_dividend_not_decreasing = False
                print('year+1', year+1)
                break
            if year not in d_dividends.keys():
                bool_dividend_not_decreasing = False
                print('year', year)
                break
            if max(d_dividends.keys())-1 not in d_dividends.keys():
                bool_dividend_not_decreasing = False
                print('max year', max(d_dividends.keys())-1)
                break
            if all([
                ## 1) no change from annual to quarterly dividend (e.g. SYK 2008-2009)
                ## 2) no seperate extraordinary dividend
                ## 3) no Yahoo forgetting about a dividend (e.g. HNZ Mar2004)
                len(d_dividends[year+1].keys()) == len(d_dividends[year].keys()),
                ## dividend not lowered on an annual basis
                any([
                    all([
                        len(d_dividends.keys()) > 15,
                        sum(d_dividends[year+1].values()) <= sum(d_dividends[year].values()),
                        ]),
                    all([
                        len(d_dividends.keys()) <= 15,
                        sum(d_dividends[year+1].values()) < sum(d_dividends[year].values()),
                        ]),
                    ]),
                ]):
                bool_dividend_not_decreasing = False
                print(year, year+1, d_dividends[year], d_dividends[year+1])
                break

        return bool_dividend_not_decreasing


    def check_if_zero(
        self,
        ticker,
        d_income_year,
        d_balance_year,
        d_cash_year,
        div_yield_min,
        index0, index1,
        ):

        bool_skip = False

        ## check that not negative (or zero)
        d_keys = {
            'income':{
                'dic': d_income_year,
                'keys': [
                    'Total Revenue',
                    'Operating Income',
                    'Net Income Before Extra. Items',
                    ],
                },
            'balance': {
                'dic': d_balance_year,
                'keys': [
                    'Total Equity',
                    'Total Assets',
                    ],
                },
            'cashflow': {
                'dic': d_cash_year,
                'keys': [
                    'Cash from Operating Activities',
                    ],
                }
            }

        for stmt in d_keys.keys():
            dic = d_keys[stmt]['dic']
            l_keys = d_keys[stmt]['keys']
            for key in l_keys:
                if not key in dic[ticker].keys():
                    print(key, 'not in', stmt, 'stmt')
                    bool_skip = True
                    return bool_skip
                print(ticker, key, dic[ticker][key][index0], type(dic[ticker][key][index0]))
                
                if any([
                    dic[ticker][key][index0] in ('-', '--',),
                    dic[ticker][key][index0] <= 0,
                    ]):
                    print(ticker, key, dic[ticker][key][index0], dic[ticker][key][index1])
                    bool_skip = True
                    return bool_skip

        ## check that present
        d_keys = {
            'income':{
                'dic': d_income_year,
                'keys': [
                    ],
                },
            'balance': {
                'dic': d_balance_year,
                'keys': [
                    ],
                },
            'cashflow': {
                'dic': d_cash_year,
                'keys': [
                    ],
                },
            }

        if div_yield_min > 0:
            d_keys['cashflow']['keys'] += ['Total Cash Dividends Paid']

        for stmt in list(d_keys.keys()):
            dic = d_keys[stmt]['dic']
            l_keys = d_keys[stmt]['keys']
            for key in l_keys:
                if not key in list(dic[ticker].keys()):
                    print(key, 'not in', stmt)
                    bool_skip = True
                    return bool_skip
                if dic[ticker][key][index0] == 0:
                    print(ticker, key, dic[ticker][key][index0], dic[ticker][key][index1])
                    if key != 'Total Cash Dividends Paid':
                        bool_skip = True
                        return bool_skip

        return bool_skip


    def calc_fundamentals_and_multiples(
        self,
        ticker, d_income_year, d_balance_year, d_cash_year,
        mc_local_currency,
        index0, index1,
        ):

        d_ratios = {
            'market': {}, 'debt': {}, 'profitability': {}, 'market': {},}

        d_fundamentals = {}
        d_multiples = {}

        ##
        ## income statement
        ##

        minuend = d_income_year[ticker]['Total Revenue'][index0]
        subtrahend = d_income_year[ticker]['Cost of Revenue, Total'][index0]
        gross_profit = minuend - subtrahend
        
        gpm_curr = (
            gross_profit / d_income_year[ticker]['Total Revenue'][index0]
            )
        opm_curr = (
            d_income_year[ticker]['Operating Income'][index0]
            /
            d_income_year[ticker]['Total Revenue'][index0]
            )
        npm_curr = (
            d_income_year[ticker]['Net Income Before Extra. Items'][index0]
            /
            d_income_year[ticker]['Total Revenue'][index0]
            )
        
        pe = mc_local_currency/(
            d_income_year[ticker]['Net Income Before Extra. Items'][index0]
            )

        gpm_prev = (
            (
                d_income_year[ticker]['Total Revenue'][index1]
                -
                d_income_year[ticker]['Cost of Revenue, Total'][index1]
                )
            /
            d_income_year[ticker]['Total Revenue'][index1]
            )

        #
        # balance sheet
        #

        if 'Total Long Term Debt' in list(d_balance_year[ticker].keys()):  # not on BW balance sheet if zero
            de = (
                d_balance_year[ticker]['Total Long Term Debt'][index0]
                /
                d_balance_year[ticker]['Total Equity'][index0]
                )
            dividend = d_balance_year[ticker]['Total Long Term Debt'][index1]
            divisor = d_balance_year[ticker]['Total Equity'][index1]
            de_prev = dividend / divisor
        else:
            de = 0
            de_prev = 0

        ##
        ## calculate enterprise value (from balance sheet)
        ## +minority interest
        ## +preferred equity
        ev = mc_local_currency
        if 'Total Long Term Debt' in list(d_balance_year[ticker].keys()):  # not on BW balance sheet if zero
            ev += d_balance_year[ticker]['Total Long Term Debt'][index0]
        if 'Cash and Short Term Investments' in list(d_balance_year[ticker].keys()):
            ev -= d_balance_year[ticker]['Cash and Short Term Investments'][index0]
        else:
            if 'Short Term Investments' in list(d_balance_year[ticker].keys()):
                if source == 'Yahoo' and d_balance_year[ticker]['Short Term Investments'][index0] in ('-', '--',):
                    d_balance_year[ticker]['Short Term Investments'][index0] = 0
                ev -= d_balance_year[ticker]['Short Term Investments'][index0]
            if 'Cash' in list(d_balance_year[ticker].keys()):
                ev -= d_balance_year[ticker]['Cash'][index0]
                if 'Cash & Equivalents' in list(d_balance_year[ticker].keys()):
##                        if d_balance_year[ticker]['Cash'][0] == 0:
##                            ev -= d_balance_year[ticker]['Cash & Equivalents'][0]
##                        else:
##                            print d_balance_year[ticker]['Cash'][0]
                        print(d_balance_year[ticker]['Cash & Equivalents'][index0])
                        stopstopstop1
            elif 'Cash & Equivalents' in list(d_balance_year[ticker].keys()):
                -d_balance_year[ticker]['Cash & Equivalents'][index0]
                if 'Short Term Investments' in list(d_balance_year[ticker].keys()):
                    print(d_balance_year[ticker]['Short Term Investments'][index0])
                    stopstopstop2
                if 'Cash' in list(d_balance_year[ticker].keys()):
                    print(d_balance_year[ticker]['Cash'][index0])
                    stopstopstop3

        ## net tangible assets (E/NTA is same as ROE, but excluding intangibles)
        ## E/NTA tells how large an investment is required to grow earnings...
        NTA = d_balance_year[ticker]['Total Equity'][index0]
        if 'Other Intangibles' in list(d_balance_year[ticker].keys()):
            NTA -= d_balance_year[ticker]['Other Intangibles'][index0]
        if 'Goodwill' in d_balance_year[ticker].keys():
            if d_balance_year[ticker]['Goodwill'][index0] != '-':
                NTA -= d_balance_year[ticker]['Goodwill'][index0]

        ##
        ## income statement and balance sheet
        ##

        roe_curr = d_income_year[ticker]['Net Income Before Extra. Items'][index0] / d_balance_year[ticker]['Total Equity'][index0]
        roe_prev = d_income_year[ticker]['Net Income Before Extra. Items'][index1] / d_balance_year[ticker]['Total Equity'][index1]
        roa_curr = d_income_year[ticker]['Net Income Before Extra. Items'][index0] / d_balance_year[ticker]['Total Assets'][index0]
        ENTA = d_income_year[ticker]['Net Income Before Extra. Items'][index0] / NTA
        d_ratios['profitability']['ROE'] = []
        for year in range(len(d_balance_year[ticker]['Total Equity'])):
            if d_balance_year[ticker]['Total Equity'][year] == 0:
                break
            d_ratios['profitability']['ROE'] += [d_income_year[ticker]['Net Income Before Extra. Items'][year] / d_balance_year[ticker]['Total Equity'][year]]

        ## calculate asset turnover (second term of dupont ROE formula)
        ATO = (
            d_income_year[ticker]['Total Revenue'][index0]
            /
            d_balance_year[ticker]['Total Assets'][index0]
            )

##            ## calculate third term of dupont ROE formula
##            ## leverage ratio = assets/equity = equity+liabilities/e (=) e+d/e = d/e + 1
        equity_multiplier = d_balance_year[ticker]['Total Assets'][index0]/d_balance_year[ticker]['Total Equity'][index0]

##        ROIC = NOPAT/IC ## wrong!!!

        ##
        ## cash flow
        ##

        fcf = (
            d_cash_year[ticker]['Cash from Operating Activities'][index0]
            -
            abs(d_cash_year[ticker]['Capital Expenditures'][index0])
            )
        print('Cash from Operating Activities', d_cash_year[ticker]['Cash from Operating Activities'][index0])
        print('Capital Expenditures', d_cash_year[ticker]['Capital Expenditures'][index0])

        ocf = d_cash_year[ticker]['Cash from Operating Activities'][index0]
        if 'Total Cash Dividends Paid' in d_cash_year[ticker].keys():
            cashflow_dividends = d_cash_year[ticker]['Total Cash Dividends Paid'][index0]
        else:
            cashflow_dividends = 0
        if cashflow_dividends in ('-', '--',):
            cashflow_dividends = 0

        #
        # cashflow market ratios
        #
        pcf = mc_local_currency/d_cash_year[ticker]['Cash from Operating Activities'][index0]

        if ev > 0 and fcf > 0:
            EVFCF = ev / fcf
            EVOCF = ev / ocf
        else:
            if 'Total Long Term Debt' in d_balance_year[ticker].keys(): ## not on BW balance sheet if zero
                print('long term debt', d_balance_year[ticker]['Total Long Term Debt'][index0])
            try:
                print(
                    'cash and short term investments',
                     d_balance_year[ticker][
                         'Cash and Short Term Investments'][index0],
                     )
            except KeyError:
                pass
            EVFCF = 0.
            EVOCF = 0.

        print('cashflow_dividends', cashflow_dividends)
        print('mc_local_currency', mc_local_currency)
        div_yield = abs(cashflow_dividends)/mc_local_currency

        #
        # append to dictionary
        #
        
        # ratios
        d_fundamentals['D/Eq'] = de
        d_fundamentals['ROE'] = roe_curr
        d_fundamentals['ROA'] = roa_curr
        d_fundamentals['NPM'] = npm_curr
        d_fundamentals['ATO'] = ATO
        d_fundamentals['equity_multiplier'] = equity_multiplier
        d_fundamentals['GPM'] = gpm_curr
        d_fundamentals['E/NTA'] = ENTA

        # sums
        d_fundamentals['FCF'] = fcf
        d_fundamentals['OCF'] = ocf

        # numbers
        d_fundamentals['cashflow_dividends'] = cashflow_dividends

        d_multiples['div_yield'] = div_yield
        d_multiples['EV/FCF'] = EVFCF
        d_multiples['EV/OCF'] = EVOCF
        d_multiples['P/E'] = pe
        d_multiples['P/CF'] = pcf

        return d_fundamentals, d_multiples, d_ratios
