import time
import os
import screener
import re


class FT:


    def parse_stmt(self, url, ticker):

        print(url)

        d_factors = {'millions':1000000.,'Bil':1000000000.}

        d = {}

        statement_error = False

        lines = screener.finance().read_url(url, ticker)

        ## ticker wrong or simply no data on URL
        if not lines:
            return None, None, None, True, None, None

#<div class="currencyDisclaimer contains"><span class="fleft">In millions of EUR<
        for i, line in enumerate(lines):
#            if 'currencyDisclaimer' in line:
#                break
            if 'mod-main-content' in line:
##            if '<table class="mod-ui-table">' in line:
                break
        line = lines[i+1]
        match = re.findall(r'>In (.*?) of (.*?)<', line)
        ## No financial data available.
        if not match:
            return None, None, None, True, None, None
        ## Millions or Billions or Thousands?
        factor = d_factors[match[0][0]]
##        matches = re.findall(
##            r'''<font face='arial' size='2'>(.*?)</font>''', line)
        currency_string = match[0][1]

##        regex = r'<table data-ajax-content="true">(.*?)</table>'
        regex = r'<table class="mod-ui-table">(.*?)</table>'
        pattern = re.compile(regex)
        match = re.search(pattern, line)
        s = match.group(1)
        regex = r'<tr class="(odd|even|Bold even|Bold odd)">(.*?)</tr>'
        regex = r'<tr(.*?)>(.*?)</tr>'
        pattern = re.compile(regex)
        match = re.findall(pattern, s)
        for m in match:
            l = re.findall(r'<t[dh].*?>(.*?)</t[dh]>', m[1])
            if not l:
                continue
            print(ticker, l)
            if 'Fiscal data as of' in l[0]:
                l[0] = 'date'
            d[l[0]] = l[1:]
            for i, x in enumerate(d[l[0]]):
                if x.startswith('(') and x.endswith(')'):
                    d[l[0]][i] = factor*-float(x[1:-1].replace(',',''))
                elif x == '--':
                    d[l[0]][i] = 0
                else:
                    try:
                        d[l[0]][i] = factor*float(x.replace(',',''))
                    except:
                        d[l[0]][i] = x

        return d, statement_error, currency_string


if __name__ == '__main__':
    x = FT()
    x.parse_stmt('http://markets.ft.com/data/equities/tearsheet/financials?s=AAP:NYQ&subview=IncomeStatement&periodType=a', 'AAP',)
