import time
import os
import screener
import re


class ADVFN:


    def parse_statements(self, url):

        print(url)

        d_factors = {'Millions':1000000.,'Bil':1000000000.}

        lines = screener.finance().read_url(url, '0')

        ## ticker wrong or simply no data on URL
        if not lines:
            return None, None, None, True, None, None

        for line in lines:
            if 'All amounts in' in line:
                break
        match = re.search(r'All amounts in (.*?) of', line)
        ## No financial data available.
        if not match:
            return None, None, None, True, None, None
        ## Millions or Billions or Thousands?
        factor = d_factors[match.group(1)]
        matches = re.findall(
            r'''<font face='arial' size='2'>(.*?)</font>''', line)
        currency_string = matches[1].strip()

        regex1 = r'<TR.*?(<td.*?</td>)</tr>'
        pattern1 = re.compile(regex1, re.DOTALL)
        regex2 = r'<td.*?>(.*?)</td>'
        pattern2 = re.compile(regex2, re.DOTALL)
        bool_init = False
        d_income = {}
        d_balance = {}
        d_cash = {}
        d = d_indicators = {}
        for line in lines:
            if 'INDICATORS' in line:
                bool_init = True
            if bool_init == False:
                continue
            if 'INCOME STATEMENT' in line:
                d = d_income
            if 'CASH-FLOW STATEMENT' in line or 'CASH FLOW STATEMENT' in line:
                d = d_cash
            if 'BALANCE SHEET' in line:
                d = d_balance
##            if '</table>' in line:
##                break
            if 'RATIOS CALCULATIONS' in line:
                break
            match1 = re.match(pattern1, line)
            if not match1:
                continue
            match2 = re.findall(pattern2, match1.group(0))
            try:
                d[match2[0]] = [factor*float(_.replace(',','')) for _ in match2[1:]]
            except ValueError:
                d[match2[0]] = match2[1:]

        statement_error = False
        print(d_income.keys())
##        ## less than 5 years of data
##        if len(d_income['total net income']) < 5:
##            statement_error = True
        ## less than 5 years of data
        if d_income['total net income'][-1] == '':
            statement_error = True

        for d in (d_income, d_balance, d_cash):
            for k in d.keys():
                d[k] = list(reversed(d[k]))  # reverse from old to new to new to old

        return d_income, d_balance, d_cash, statement_error, currency_string, d_indicators
