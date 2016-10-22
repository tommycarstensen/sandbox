## rows not sorted by year...
## http://investing.money.msn.com/investments/financial-statements?symbol=aap

## mixed Mil and Bil
## http://investing.money.msn.com/investments/financial-statements?symbol=AU:CSL

import urllib.request, urllib.error, urllib.parse, math, ticker_conversion, os, time
import re
import screener

class morningstar:

    def parseKeyRatios(self,url, ticker):

        print(url)

        dic_10year = {}

        d_factors = {'Mil':1.,'Bil':1000.}

##        url = 'http://financials.morningstar.com/financials/getFinancePart.html?&callback=?&t=AAP'
##        url = 'http://financials.morningstar.com/financials/getFinancePart.html?&callback=?&t=ABT'

        lines = screener.finance().read_url(url, ticker)
        s_html = '\n'.join(lines)

        ## ticker wrong or simply no data on morningstar.com
        if s_html == '' or s_html == '?({"componentData":null})':
            return dic_10year

##        print(s_html)
        for k in (
            'Revenue', 'Operating Income', 'Net Income', 'Earnings Per Share',
            'Dividends', 'Shares', 'Book Value Per Share',
            'Operating Cash Flow', 'Cap Spending', 'Free Cash Flow',
            'Free Cash Flow Per Share', 'Working Capital',):
##            print(k)
            ## \w\w\w is the currency code
            p = r'>{}&nbsp;<span>\w\w\w.*?<\\/tr\>'.format(k)
##            print(p)
            s_tr = re.search(p, s_html).group()
            try:
                factor = d_factors[re.search(r'[MB]il', s_tr).group()]
            except AttributeError:
                factor = 1
            l = []
            for match in re.finditer(r'>([-\d,.]+|&mdash;)<',s_tr):
                s = match.group(1)
                if s == '&mdash;':
                    l.append('-')
                else:
                    l.append(factor * float(s.replace(',','')))
            ## Do not include TTM.
            if len(l) == 11:
                l = l[:-1]
            else:
                print(l)
                print(s_tr)
                stop
            l = list(reversed(l))
            assert len(l) == 10
            if l.count('-') == 1:
                i = l.index('-')
                ## if earliest year, then same as year after
                if i == 9:
                    l[9] = l[8]
                ## otherwise average of neigbouring years
                else:
                    l[i] = (l[i-1] + l[i+1])/2
            if k == 'Revenue' and '-' in l:
                print(ticker, k, l)
                return {}
##            l = list(
##                factor*float(x.replace(',','')) for x in reversed(
##                    re.findall(r'>([-\d,.]+)<',substr)))[1:]
##                    re.findall(r'>([-\d,.]+|&mdash;)<',substr)))[1:]
            print(k, len(l), l)
            ## often initial bvps missing, so just assume 10% lower last year
            ## fudge factor big time!
            if len(l) == 9 and k in (
                'Book Value Per Share',
                '',
                ):
                l.append(0.9*l[-1])
            if len(l) < 10 and k not in (
                ## Some of them not used...
                'Free Cash Flow Per Share',
                'Working Capital',
                'Dividends',
                'Book Value Per Share',
                ):
                print(k, len(l), l)
                return {}
            dic_10year[k] = l

##        print('string\n',s[s.index('Shares')-2*80:s.index('Shares')+10*80],'\nstring',)
##        print(s[s.index('2003')-2*80:s.index('Shares')+10*80])
        s_tr = re.search(r'>Shares&nbsp;<span>[MB]il.*?<\\/tr\>', s_html).group()
        factor = d_factors[re.search(r'[MB]il', s_tr).group()]
        l_shares = list(
            factor*float(x.replace(',','')) for x in reversed(
                re.findall(r'>([\d,]+)<', s_tr)))[1:]
##        print(l_shares)

        p = r'<th scope=\\"col\\" align=\\"right\\" id=\\"Y\d+\\">(\d\d\d\d-\d\d)<\\/th>'
        l_dates = list(reversed(re.findall(p,s)))
##        print(l_dates)

        dic_10year['DATE'] = l_dates
        dic_10year['SHARES OUTSTANDING'] = l_shares
        dic_10year['SALES'] = dic_10year['Revenue']

        return dic_10year

if __name__ == '__main__':
    x = morningstar()
    x.parseKeyRatios('http://financials.morningstar.com/financials/getFinancePart.html?&callback=?&t=EBAY', '0')
    exit()
    x.parseKeyRatios('http://financials.morningstar.com/financials/getFinancePart.html?&callback=?&t=XTKS:8267', '0')
    exit()
    x.parseKeyRatios('http://financials.morningstar.com/financials/getFinancePart.html?&callback=?&t=XTKS:5333', '0')
    exit()
    x.parseKeyRatios('http://financials.morningstar.com/financials/getFinancePart.html?&callback=?&t=XBRU:ONTEX', '0')
