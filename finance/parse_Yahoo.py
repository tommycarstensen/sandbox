import time
import os
import screener
import re
import html
from html.parser import HTMLParser
from bs4 import BeautifulSoup


class MyHTMLParser(HTMLParser):
    def handle_starttag(self, tag, attrs):
        print("Encountered a start tag:", tag, attrs)
    def handle_endtag(self, tag):
        print("Encountered an end tag :", tag)
    def handle_data(self, data):
        print("Encountered some data  :", data)


class Yahoo:


    def parse_ownership(self, d):

        for ticker in d.keys():
            url = 'https://finance.yahoo.com/q/mh?s={}+Major+Holders'.format(ticker)
            lines = screener.finance().read_url(url, '0')
            for line in lines:
##                print(line)
                r = r'<table.*?><tr><th.*?>Holder</th><th.*?>Shares</th>.*?(<tr>.*?</tr>)</table>'
                p = re.compile(r)
                matches = re.finditer(p, line)
                for match in matches:
                    soup = BeautifulSoup(match.group(0))
                    letters = soup.find_all("td", class_="yfnc_tabledata1")
                    print(letters)
                    for element in letters:
                        print(element)

        return d


    def parseIncomeStatement(self, url):

        print(url)

        d_factors = {'Mil':1.,'Bil':1000.}

        lines = screener.finance().read_url(url, '0')
        s_html = s = '\n'.join(lines)
##        print(lines)
##        stop

        ## ticker wrong or simply no data on morningstar.com
        if s == '':
            return {}

##        ## DOTALL = Make the '.' special character match any character at all, including a newline; without this flag, '.' will match anything except a newline.
##        p = re.compile(regex, re.DOTALL)
##        m = re.search(p, s)
##        print(m.group())
##        from xml.etree import ElementTree as ET
##        table = ET.XML(m.group())
##        rows = iter(table)
##        headers = [col.text for col in next(rows)]
##        for row in rows:
##            values = [col.text for col in row]
##            print(values)
##        print(headers)
##        stop
##
##        parser = MyHTMLParser()
##        parser.feed(s)
##        stop

        r = r'<TABLE class="yfnc_tabledata1".*?><TR><TD>(<TABLE.*?></TABLE>)</TD></TR></TABLE>'
        p = re.compile(r, re.DOTALL)
        m = re.search(p, s_html)
        if not m:
            return None, True, None
        s_table = m.group(1)

##        r = r'Period Ending.*?<b>(.*?)</b></TD></TR>'
##        p = re.compile(r, re.DOTALL)
##        matches = re.finditer(p, s_table)
##        for match in matches:
##            print(match.group(0))

        ## regex for tr
##        regex1 = r'<tr>(<td.*?>.*?</td>)</tr>'
##        regex1 = r'<TABLE class="yfnc_tabledata1".*?(<tr><td.*?>.*?</td></tr>)</TABLE>'
##        regex1 = r'<TABLE class="yfnc_tabledata1".*?></TR>(<t[r]>.*?</t[r]>)</TABLE>'
        r = r'<[tT][rR]>(.*?)</[tT][rR]>'
        ## DOTALL = Make the '.' special character match any character at all, including a newline; without this flag, '.' will match anything except a newline.
        p = re.compile(r, re.DOTALL)
        matches = re.finditer(p, s_table)
        d = {}
        for match in matches:
##            print('match', match.group())
            if 'style="display:block' in match.group(0):
                continue
            s_tr = match.group(0)
            s_tr = s_tr.replace('<strong>','').replace('</strong>','')
            # Replace HTML tags with an empty string.
##            regex2 = r'<td.*?>.*?([-\d,.]+).*?</td>'
##            matches2 = re.findall(r'<td.*?>(.*?)</td>', s_tr, re.DOTALL)
            matches2 = re.findall(
                r'<[tT][dD].*?>(.*?)</[tT][dD]>', s_tr, re.DOTALL)
            if len(matches2) <= 3:
                continue
##            print(matches2)
            if len(matches2) >= 5 and 'spacer' in match.group(0):
                matches2 = matches2[1:]
            k = re.sub('<.*?>', '', matches2[0].strip())
##            print(matches2)
            l = []
            for i in range(1,4):
                s = matches2[i]
                s = re.sub('<.*?>', '', s)
                s = s.replace(',','').replace('&nbsp;','').strip()
                if k == 'Period Ending':
                    l.append(int(s[-4:]))
                elif s == '-':
                    l.append(s)
                else:
                    ## 1000 if all numbers in thousands
                    assert 'All numbers in thousands' in s_html
                    factor = 1000
                    if s[0] == '(' and s[-1] == ')':
                        l.append(-factor*float(s[1:-1]))
                    else:
                        l.append(factor*float(s))
            d[k] = l

        m = re.search('>Currency in (\w*)', s_html)
        currency = m.group(1)

        return d, False, currency


    def parse_historical_prices(self, url):

        print(url)

        ## http://real-chart.finance.yahoo.com/table.csv?s=AVV.L&g=v&ignore=.csv

        ## key = year, value = sum of dividends
        d = {}

        lines = screener.finance().read_url(url, '0')
        for line in lines[1:]:
            date, dividend = line.split(',')
            year = int(date[:4])
            dividend = float(dividend)
            try:
                d[year] += dividend
            except KeyError:
                d[year] = dividend

        return d


    def parse_currencies(self):

        d_currency_reuters = {}

        ## ISO 4217 codes
        l_currency_msn = [
            ## most traded currencies
            'Euro','EUR',
            'Japanese Yen','JPY',
            'British Pounds','GBP',
            'Swiss Francs','CHF',
            'Australian Dollars','AUD',
            ## Asia
            'Chinese Renminbi','CNY', ## China
            'Taiwanese Dollars','TWD',
            'Hong Kong Dollars','HKD',
            'Philippine Pesos','PHP',
            'Singapore Dollars','SGD',
            'Indonesian Rupiah','IDR', ## Indonesia
            'Indian Rupee','INR',
            'Thai Bahts','THB',
            'South Korean Won','KRW',
            'Malaysian Ringgit','MYR', ## Malaysia
            ## Europe
            'Swedish Krona','SEK',
            'Norwegian Krone','NOK',
            'Danish Krone','DKK',
##            'Iceland Krona','ISK',
            'Icelandic Kronas','ISK',
            'Hungarian Forint','HUF',
            'Czech Korunas','CZK',
##            'Estonian Kroon','EEK', ## Euro 2011-
            'Lithuanian Lita','LTL',  # Euro 2015-
            'Russian Rouble','RUB', ## msn
            'Turkish Lira','TRY',
            'Polish Zlotys','PLN',
            ## North America
            'Canadian Dollars','CAD',
            'Mexican Pesos','MXN',
            ## Oceania
            'N.Z. Dollars','NZD',
            ## South America
            'Colombian Peso','COP',
            'Argentine Peso','ARS',
            'Chilean Peso','CLP',
            'Brazilian Real','BRL',
            'Peruvian Nuevo Sol','PEN',
##            'New Sol','PEN', 'Bolivar','VEB',
            ## Middle East
            'Israeli Shekel','ILS', ## Israel
            'Kuwait Dinars','KWD',
            'Saudi Arabian Riyals','SAR',
            'Qatari Rials','QAR',
            ## Africa
            'South African Rand','ZAR', ## South Africa (ZAR...)
            'Nigerian Naira','NGN',
            ]

        d_currency_msn = {}

        for i in range(0,len(l_currency_msn),2):

            name = l_currency_msn[i]
            symbol = l_currency_msn[i+1]

            url = 'http://download.finance.yahoo.com/d/quotes.csv?s=USD%s=X&f=sl1d1t1c1ohgv&e=.csv' %(symbol,)
##            print url

            lines = screener.finance().read_url(url, '0')
            
            line = lines[0]

            rate = float(str(line).split(',')[1])
            if rate == 0:
                print(symbol)
                print(rate)
                print(url)
                print(line)
                stop_zero_rate
            d_currency_reuters[symbol] = rate
            d_currency_msn[name] = rate
            print('%s%s %3s %7.2f' %(name, (24-len(name))*'-', symbol, rate))

        d_currency_msn['U.S. Dollars'] = 1.
        d_currency_reuters['USD'] = 1.
        d_currency_msn['US Dollars'] = 1.
##        d_currency_msn['GBX'] = d_currency_msn['GBP']*100.

        d_currency_reuters['TRL'] = d_currency_reuters['TRY']
        d_currency_reuters['ZAX'] = d_currency_reuters['ZAR']
        d_currency_reuters['GBX'] = d_currency_reuters['GBP']*100.

        for currency in list(d_currency_msn.keys()):
            d_currency_msn['%ss' %(currency)] = d_currency_msn[currency]
        d_currency_msn['Chinese Renminbi (Yuan)s'] = d_currency_msn['Chinese Renminbi']
        d_currency_msn['Taiwan Dollars'] = d_currency_msn['Taiwanese Dollars']
##        d_currency_msn['Lithuanian Litass'] = d_currency_msn['Lithuanian Lita']
        d_currency_msn['Turkish New Liras'] = d_currency_msn['Turkish Lira']
        d_currency_msn['Philippines Pesos'] = d_currency_msn['Philippine Pesos']
        d_currency_msn['New Zealand Dollars'] = d_currency_msn['N.Z. Dollars']
##        d_currency_msn['Qatari Rial'] = d_currency_msn['Qatari Rial']
        d_currency_msn['Won'] = d_currency_msn['South Korean Won']

        return d_currency_msn, d_currency_reuters


    def parse_earnings_date(self,data,):

        print('parsing earnings dates')

        tickers = list(data.keys())
        tickers.sort()

        for i in range(len(tickers)):

            ticker = tickers[i]
            if i % 10 == 0:
                print('%s/%s %s' %(i+1, len(tickers), ticker))

            url = 'http://biz.yahoo.com/research/earncal/%s/%s.html' %(ticker[0], ticker.lower())
            lines = screener.finance().read_url(url, ticker)

            for i in range(len(lines)):
                line = lines[i]
                if 'Earnings Calendar for' in line:
                    index1 = 0
                    index2 = lines[i+1].index('<')
                    results = lines[i+1][index1:index2]
                    year = int(results[-4:])
                    month = results[:3]
                    day = int(results[-8:-6])
                    results = '%s%s%s' %(year, month, str(day).zfill(2))

            if lines in [[],[''],]:
                data[ticker]['Date'] = '%s' %('N/A')
            else:
                data[ticker]['Date'] = '%s' %(results)

        return data

if __name__ == '__main__':
    x = Yahoo()
    x.parse_ownership({'MSFT': 9})
