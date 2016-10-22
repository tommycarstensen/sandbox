## only want this script to be dependent on a set of tickers

import ticker_conversion, os, time
import screener

class TA:

    def find_candidates_TA(
        self,
        l_tickers, l_time, months, l_statementNA, d_portfolio,
        d_ADR,
        ):

        print('finding TA candidates')

        year2 = l_time[0] ; month2 = l_time[1] ; day2 = l_time[2]
        year1 = year2-11 ; month1 = month2  ; day1 = day2

        TAcandidates = []
        TAdata = {}

        l_supports = []
        l_breakouts = []

        l_MA50_increasing = []
        l_MA50_bounce = []
        l_52w_low = []
        l_down10percent_morethansp500 = []

        ##
        ## S&P 500
        ##
        ticker = ticker_yahoo = '^GSPC' ## S&P500
        TAdata[ticker] = {}
        period = 'daily'
        TAdata[ticker][period] = {}
        ## read url
        url = 'http://ichart.finance.yahoo.com/table.csv?s=%s&d=%s&e=%s&f=%s&g=d&a=%s&b=%s&c=%s&ignore=.csv' %(
            ticker_yahoo, month2-1, day2, year2, month1-1, day1, year1,
            ) ## g=d signifies daily_weekly_monthly graph
        linesd = screener.finance().read_url(url, ticker)
        data = TAdata[ticker][period]['raw'] = linesd[1:]
        ## parse lines
        TAdata[ticker][period]['price'] = {
            'date':[], 'open':[], 'high':[], 'low':[], 'close':[], 'volume':[],'adjclose':[],
            }
        TAdata = self.data_conversion(ticker,period,data,TAdata,)
        ## price today
        price_today = TAdata[ticker][period]['price']['adjclose'][-1]
        date_today = TAdata[ticker][period]['price']['date'][-1]
        ## price 52w
        date_52w = '%4s%s' %(int(date_today[:4])-1,date_today[4:])
        price_52w = None
        for i in range(2,len(TAdata[ticker][period]['price']['date'])):
            if TAdata[ticker][period]['price']['date'][-i] <= date_52w:
                price_52w = TAdata[ticker][period]['price']['adjclose'][-i]
                break
        ## change 52w
        sp500_52w_change = (price_today-price_52w)/price_52w

        for ticker in l_tickers:

            ticker_FA = ticker

##            if ticker[-2:] == '.I':
##                continue
##            if ticker[-3:] in [
##                '.HE','.VX','.IS','.BR','.MM',
##                '.MX','.SA',
##                '.HK','.BO',
##                ]:
##                continue

            if ticker[-3:] in [
##                '.IC', ## Iceland not on Yahoo
##                '.SI', ## Singapore not on Yahoo
                '.BO', ## India not on Yahoo
##                '.ME', ## Russia not on Yahoo
                ]:
                continue
##            if ticker == 'SUN.BO':
##                continue
##            if ticker == 'WIPR.BO':
##                continue
            if ticker == 'INGC.BO':
                continue
            if ticker == 'HUVR.BO':
                continue

            if '.' in ticker and ticker[-2:] in ['.A','.B',] and ticker[-2:] not in ['.O']:
                index = ticker.index('.')
                ticker = ticker[:index]+'-'+ticker[index+1:]
            ticker = ticker.replace('.a','-a')
            ticker = ticker.replace('.b','-b')
            ticker = ticker.replace('b','-B') ## HUBb, NOVO-B.CO
            ticker = ticker.replace('a','-A') ## BFa
            if ':' in ticker:
                index = ticker.index(':')
##                if ticker[:index] == 'JP': ## Japan not on Yahoo
##                    continue
                if ticker[:index] == 'CA' and '.' in ticker:
                    ticker.replace('.','-')
                    stop
##                if ticker[:index] == 'SE' and '-' in ticker:
##                    ticker = ticker.replace('-','')
                ticker = ticker_conversion.unknown2yahoo(ticker)
            elif '.' in ticker:
                ticker = ticker_conversion.unknown2yahoo(ticker)
            ticker = ticker.replace('..','.') ## RB..L

            ticker_yahoo = ticker

            ticker = ticker_FA

##            if ticker in d_yahoo2reuters:

            ##
            ## parse historical data
            ##

            TAdata[ticker] = {}

            ## daily
            url = 'http://ichart.finance.yahoo.com/table.csv?s=%s&d=%s&e=%s&f=%s&g=d&a=%s&b=%s&c=%s&ignore=.csv' %(
                ticker_yahoo, month2-1, day2, year2, month1-1, day1, year1,
                ) ## g=d signifies daily_weekly_monthly graph

            linesd = screener.finance().read_url(url, ticker)
            fp = 'urls/%s' %(url.replace(':','').replace('/','').replace('.','').replace('?',''))

            ## no data
            if linesd == ['']: continue

##            ## weekly
##            url = 'http://ichart.finance.yahoo.com/table.csv?s=%s&a=%s&b=%s&c=%s&d=%s&e=%s&f=%s&g=w&ignore=.csv' %(
##                ticker_yahoo, month1-1, day1, year1, month2-1, day2, year2,
##                ) ## g=w signifies daily_weekly_monthly graph
##            for x in range(10):
##                try:
##                    urllines = urllib2.urlopen(url)
##                    linesw = urllines.readlines()
##                    break
##                except:
##                    print x, url
##                    continue
##            if x == 9:
##                continue
##
##            ## monthly
##            url = 'http://ichart.finance.yahoo.com/table.csv?s=%s&a=%s&b=%s&c=%s&d=%s&e=%s&f=%s&g=w&ignore=.csv' %(
##                ticker_yahoo, month1-1, day1, year1, month2-1, day2, year2,
##                ) ## g=w signifies daily_weekly_monthly graph
##            for x in range(10):
##                try:
##                    urllines = urllib2.urlopen(url)
##                    linesm = urllines.readlines()
##                    break
##                except:
##                    print x, url
##                    continue
##            if x == 9:
##                continue

            TAdata[ticker]['daily'] = {
                'raw':linesd[1:],
                }

            ## find TA candidates
            periods = list(TAdata[ticker].keys())
            TAcandidate = True
            for period in [
                'daily',
##                'weekly','monthly',
                ]:

                TAdata[ticker][period]['price'] = {
                    'date':[], 'open':[], 'high':[], 'low':[], 'close':[], 'volume':[],'adjclose':[],
                    }

                data = TAdata[ticker][period]['raw']
                n = len(data)

                TAdata = self.data_conversion(ticker,period,data,TAdata,)

                ## calculate MA
                if period == 'daily':
                    TAdata, MA50, MA200, l_MA50_increasing, l_MA50_bounce = self.MA(
                        ticker,period,TAdata,l_MA50_increasing,l_MA50_bounce,
                        )

                    TAdata[ticker][period]['MA50'] = MA50
                    TAdata[ticker][period]['MA200'] = MA200

##                    print ticker, 'ma50', MA50, 'ma200', MA200
                    price_today = TAdata[ticker][period]['price']['adjclose'][-1]
                    date_today = TAdata[ticker][period]['price']['date'][-1]

                    date_52w = '%4s%s' %(int(date_today[:4])-1,date_today[4:])
                    price_52w = None
                    for i in range(2,len(TAdata[ticker][period]['price']['date'])):
                        if TAdata[ticker][period]['price']['date'][-i] <= date_52w:
                            price_52w = TAdata[ticker][period]['price']['adjclose'][-i]
                            break

                    if price_52w == None:
                        continue

                    date_10y = '%4s%s' %(int(date_today[:4])-10,date_today[4:])
                    price_10y = None
                    for i in range(2,len(TAdata[ticker][period]['price']['date'])):
                        if TAdata[ticker][period]['price']['date'][-i] <= date_10y:
                            price_10y = TAdata[ticker][period]['price']['adjclose'][-i]
                            break
                        
                    l_prices_52w = []
                    for i in range(2,len(TAdata[ticker][period]['price']['date'])):
                        if TAdata[ticker][period]['price']['date'][-i] >= date_52w:
                            l_prices_52w += [TAdata[ticker][period]['price']['adjclose'][-i]]
                            continue

                    if price_52w:
                        change_52w = (
                            price_today
                            -
                            price_52w
                            ) / price_52w
                        TAdata[ticker][period]['change_52w'] = round(100*change_52w,0)
                    else:
                        change_52w = None

                    if price_10y:
                        change_10y = (price_today-price_10y)/price_10y
                        TAdata[ticker][period]['change_10y'] = round(100*change_10y,0)
                    else:
                        change_10y = None

                    price_52w_min = min(l_prices_52w)
                    price_52w_max = max(l_prices_52w)

                    above_52w = (price_today-price_52w_min)/price_52w_min
                    below_52w_max = (price_today-price_52w_max)/price_52w_max

                    TAdata[ticker][period]['above_52w'] = round(100*above_52w,0)
                    TAdata[ticker][period]['below_52w_max'] = round(100*below_52w_max,0)

                    if price_today < 1.05*price_52w_min:
                        l_52w_low += [ticker]

                    ## dropped more than 10% relative to market
                    if (price_today-price_52w)/price_52w < sp500_52w_change-0.1:
                        l_down10percent_morethansp500 += [ticker]

                ## find support and resistance
                ## conflicts if support or resistance while paying out dividend...
                if period == 'daily':
                    l_supports, l_breakouts = self.support_and_resistance(
                        ticker,TAdata,
                        l_supports,l_breakouts,
                        )

                ## find gap support/resistance
                if period == 'daily':
                    l_gaps = self.gaps(ticker,data,)

                ## calculate RSI
                if period == 'daily':
                    TAdata = self.RSI(ticker,period,TAdata)

                ## calculate MFI
                if period == 'daily':
                    TAdata = self.MFI(ticker,period,TAdata)

                ## calculate MACD
                TAdata = self.MACD(ticker,period,TAdata,)

                ## evaluate MACD (bullish)
                if period != 'monthly' and not (TAdata[ticker][period]['MACD']['DIV'][-1] > TAdata[ticker][period]['MACD']['DIV'][-2] and TAdata[ticker][period]['MACD']['DIV'][-2] < 0):
                    TAcandidate = False
                elif period == 'monthly' and not TAdata[ticker][period]['MACD']['DIV'][-2] < 0:
                    TAcandidate = False

                ## end of loop over periods


##            ## evaluate MACD (bullish)
##            if (
####                TAdata[ticker]['daily']['MACD']['DIV'][-1] > TAdata[ticker]['daily']['MACD']['DIV'][-2]
####                and
####                TAdata[ticker]['daily']['MACD']['DIV'][-2] < 0
####                and
##                TAdata[ticker]['weekly']['MACD']['DIV'][-1] > TAdata[ticker]['weekly']['MACD']['DIV'][-2]
##                and
##                TAdata[ticker]['weekly']['MACD']['DIV'][-2] < 0
##                and
##                TAdata[ticker]['monthly']['MACD']['DIV'][-2] < 0
##                ):
##                    TAcandidates.append(ticker)
##                    print 'TAcandidate!!!'
##
##            ## evaluate MACD (bearish)
##            if ticker in d_portfolio.keys() and ticker not in l_statementNA:
##                if (
##                    TAdata[ticker]['daily']['MACD']['DIV'][-1] > 0 and
##                    TAdata[ticker]['weekly']['MACD']['DIV'][-1] > 0 and
##                    TAdata[ticker]['monthly']['MACD']['DIV'][-1] > 0
##                    ):
##                    print 'SELL %s !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' %(ticker)

        d_TA = {
            'MACD':TAcandidates,
            'bouncing at support level':l_supports,
            'breaking resistance level':l_breakouts,
            'MA50 increasing':l_MA50_increasing,
            'MA50 bounce':l_MA50_bounce,
            }
        for s_TA in list(d_TA.keys()):
            l_TA = d_TA[s_TA]
            yahoo = 'http://finance.yahoo.com/q/cq?d=v1&s='
            for ticker in l_TA:
                yahoo += '%s+' %(ticker)
            print('\n')
            print(s_TA)
            print((yahoo[:-1]))

        print('\n')
        print(('l_52w_low', l_52w_low))
        print(('l_down10percent_morethansp500', l_down10percent_morethansp500))
        print('\n')

        s_CAPS = ''
        for ticker in l_TA:
            if ticker in list(d_ADR.values()):
                for ADR,v in list(d_ADR.items()):
                    if v == ticker:
                        break
                ticker_US = ADR
            else:
                ticker_US = ticker
            s_CAPS += '%s,' %(ticker_US)
        print((s_CAPS[:-1]))

        fd = open('TAcandidates.txt', 'w')
        fd.write('%s\n%s' %(l_tickers, TAcandidates))
        fd.close()

##        matrix = self.covar_matrix(d_pca)
##        eigenvalues,eigenvectors = self.diagonalization(matrix)
##        print 'matrix', matrix
##        print 'eval', eigenvalues
##        print 'evec', eigenvectors[0]
##        print 'tickers', l_tickers
##        for i in range(len(l_tickers)):
##            print '%s \t %s' %(l_tickers[i],eigenvectors[0][i],)

        return TAcandidates, TAdata


    def data_conversion(
        self,ticker,period,data,TAdata,
        ):

        ## loop forwards in time
        for j in range(len(data)-1,-1,-1):

            ## convert data string to data list
            data[j] = data[j].split(',')[:] # Date,Open,High,Low,Close,Volume,Adj. Close

            ## convert to floats
            data[j][1] = float(data[j][1]) ## open
            data[j][2] = float(data[j][2]) ## high
            data[j][3] = float(data[j][3]) ## low
            data[j][4] = float(data[j][4]) ## close
            data[j][6] = float(data[j][6]) ## adjusted close

            ## convert date format from dd-mmm-yy to dd/mm/yy
            date = data[j][0]
            index = date.index('-')+1
            month = date[index:index+3]
##                    if month in months:
##            date = date.replace(month,str(months[month]).zfill(2)).replace('-','/').replace(',',' ')
            date = date.replace('-','/').replace(',',' ')
            data[j][0] = date

            TAdata[ticker][period]['price']['date'] += [date]
            ## write data to lists to make it possible to sum easily
            TAdata[ticker][period]['price']['open'] += [data[j][1]]
            TAdata[ticker][period]['price']['high'] += [data[j][2]]
            TAdata[ticker][period]['price']['low'] += [data[j][3]]
            TAdata[ticker][period]['price']['close'] += [data[j][4]]
            TAdata[ticker][period]['price']['volume'] += [data[j][5]]
            TAdata[ticker][period]['price']['adjclose'] += [data[j][6]]

        return TAdata


    def MA(
        self,ticker,period,TAdata,
        l_MA50_increasing,
        l_MA50_bounce,
        ):

        '''adjusted close is used'''
        
        TAdata[ticker][period]['MA50'] = []
        TAdata[ticker][period]['MA200'] = []
        for j in range(len(TAdata[ticker][period]['price']['close'])):
            if j >= 50-1:
                MA50 = sum(TAdata[ticker][period]['price']['adjclose'][j-(50-1):j+1])/50.
            else:
                MA50 = 0.
            if j >= 200-1:
                MA200 = sum(TAdata[ticker][period]['price']['adjclose'][j-(200-1):j+1])/200.
            else:
                MA200 = 0.
            TAdata[ticker][period]['MA50'] += [MA50]
            TAdata[ticker][period]['MA200'] += [MA200]

        if TAdata[ticker][period]['MA50'][-1] > TAdata[ticker][period]['MA50'][-2]:
            l_MA50_increasing += [ticker]

        if (
            ## MA50 below price
            TAdata[ticker][period]['MA50'][-1] < TAdata[ticker][period]['price']['adjclose'][-1]
            and
            ## MA50 and price within 2% of each other
            (TAdata[ticker][period]['MA50'][-1]-TAdata[ticker][period]['price']['adjclose'][-1])/TAdata[ticker][period]['price']['adjclose'][-1] < 0.02
            ):
            l_MA50_bounce += [ticker]

        return TAdata, MA50, MA200, l_MA50_increasing, l_MA50_bounce


    def gaps(self,ticker,data,):

        levels = []
        for j in range(len(data)):
            if data[j][0] == data[0][0][:5]+str(int(data[0][0][5:7])-6).zfill(2)+data[0][0][7:]:
                break
##        for k in range(j-1,-1,-1):
##                        if float(data[k][1]) == 87.31:
##                            print data[k][0]

        return levels


    def support_and_resistance(
        self,ticker,TAdata,
        l_supports,l_breakouts,
        ):

        l_lows = TAdata[ticker]['daily']['price']['low']
        l_highs = TAdata[ticker]['daily']['price']['high']
        l_close = TAdata[ticker]['daily']['price']['close']
        l_adjclose = TAdata[ticker]['daily']['price']['adjclose']
    
        if len(l_adjclose) > 200:
            scale = 0.05*l_adjclose[-1]
##                    if l_adjclose[-1] > 5 and l_adjclose[-1] < 20:
##                        scale = .5
##                    elif l_adjclose[-1] > 20 and l_adjclose[-1] < 100:
##                        scale = 1.
##                    elif l_adjclose[-1] > 100 and l_adjclose[-1] < 200:
##                        scale = 2.
##                    elif l_adjclose[-1] > 200 and l_adjclose[-1] < 500:
##                        scale = 4.
##                    elif l_adjclose[-1] > 500 and l_adjclose[-1] < 1000:
##                        scale = 5.
##                    else:
##                        print l_adjclose[-1]
##                        print 'http://stockcharts.com/help/doku.php?id=support:understanding_pnf_charts'
##                        stop_scale
            minimum = min(l_lows[-200:])
            maximum = max(l_highs[-200:])
            ## no min/max if within first five days
            index = l_lows[-200:].index(minimum)
            if index in range(5):
                minimum = 'N/A'
##                    print minimum, TAdata[ticker][period]['price']['date'][-200+index]
            index = l_highs[-200:].index(maximum)
            if index in range(5):
                maximum = 'N/A'
##                    print maximum, TAdata[ticker][period]['price']['date'][-200+index]
            range_max_low = l_highs[-199]+l_adjclose[-199]-l_close[-199]
            range_max_high = l_highs[-199]+l_adjclose[-199]-l_close[-199]
            range_min_low = l_lows[-199]+l_adjclose[-199]-l_close[-199]
            range_min_high = l_lows[-199]+l_adjclose[-199]-l_close[-199]
##                    print TAdata[ticker][period]['price']['date'][-199]
            for j in range(-200+1+5,-6):
                diff = l_adjclose[j]-l_close[j]
                adjhigh = l_highs[j]+diff
                adjlow = l_lows[j]+diff
##                        print j,TAdata[ticker][period]['price']['date'][j], '%.2f %.2f %.2f' %(l_highs[j]+diff,l_lows[j]+diff,l_adjclose[j])
                if (
                    l_adjclose[j] < max(l_adjclose[j-5:j+6]) and
                    l_adjclose[j] > min(l_adjclose[j-5:j+6])
                    ):
                    continue
                if (
                    l_highs[j] < range_max_low and
                    l_lows[j] > range_min_high
                    ):
                    continue
                if (
                    l_highs[j] > range_max_low and
                    (
                        l_adjclose[j] >= max(l_adjclose[j-5:j+6]) or
                        l_highs[j] >= max(l_highs[j-5:j+6])
                        )
                    ):
                    if l_adjclose[j] > range_max_low+scale/2:
                        range_max_low = l_adjclose[j]-scale/2
                        range_max_high = max(adjhigh,l_adjclose[j]+scale/2)
##                                print '***HI***', TAdata[ticker][period]['price']['date'][j], l_adjclose[j], range_max_low, range_max_high, '%.2f %.2f %.2f' %(l_highs[j]+diff,l_lows[j]+diff,l_adjclose[j]), max(l_highs[j-5:j+6])+diff, max(l_adjclose[j-5:j+6])
                if (
                    l_lows[j] < range_min_high and
                    (
                        l_adjclose[j] <= min(l_adjclose[j-5:j+6]) or
                        l_lows[j] <= min(l_lows[j-5:j+6])
                        )
                    ):
                    if l_adjclose[j] < range_min_high+scale/2:
                        range_min_high = l_adjclose[j]+scale/2
                        range_min_low = min(adjlow,l_adjclose[j]-scale/2)
##                                print '***LO***', TAdata[ticker][period]['price']['date'][j],
##                                print l_adjclose[j],
##                                print range_min_low, range_min_high,
##                                print '%.2f %.2f' %(
##                                    l_highs[j]+diff,l_lows[j]+diff,
##                                    ),
##                                print min(l_lows[j-5:j+6])+diff, min(l_adjclose[j-5:j+6])

##                        print range_max_low,range_max_high
##                        print range_min_low,range_min_high
##                        print l_adjclose[-1],minimum,maximum
##                    if l_adjclose[-1] < range_min_high and l_adjclose[-1] > range_min_low:
##                    if l_adjclose[-1] < range_min_high and l_lows[-1] > minimum:
            if l_adjclose[-1] < range_min_high and (
                l_adjclose[-1] > l_adjclose[-2] and
                (
                    l_highs[-1] > l_highs[-2] or
                    l_lows[-1] > l_lows[-2]
                    )
                ):
                l_supports += [ticker]
            if l_adjclose[-1] > range_max_high:
                l_breakouts += [ticker]

        return l_supports, l_breakouts


    def MACD(self,ticker,period,TAdata,):
        
        ## append initial TA values to historical data
        TAdata[ticker][period]['MACD'] = {
            'EMAshort':[], 'EMAlong':[], 'MACD':[], 'EMA':[], 'DIV':[],
            }

        close = float(TAdata[ticker][period]['price']['close'][0])
        TAdata[ticker][period]['MACD']['EMAshort'] += [close]
        TAdata[ticker][period]['MACD']['EMAlong'] += [close]
        TAdata[ticker][period]['MACD']['MACD'] += [0.]
        TAdata[ticker][period]['MACD']['EMA'] += [0.]
        TAdata[ticker][period]['MACD']['DIV'] += [0.]

        ## MACD
        for j in range(1,len(TAdata[ticker][period]['price']['close'])):

            close = float(TAdata[ticker][period]['price']['close'][j])

            prevEMAshort = TAdata[ticker][period]['MACD']['EMAshort'][j-1]
            EMAshort = (2./(1.+12.))*(close-prevEMAshort)+prevEMAshort
            prevEMAlong = TAdata[ticker][period]['MACD']['EMAlong'][j-1]
            EMAlong = (2./(1.+26.))*(close-prevEMAlong)+prevEMAlong

            MACD = (EMAshort-EMAlong)
            prevEMA = TAdata[ticker][period]['MACD']['EMA'][j-1]
            EMA = (2./(1.+9.))*(MACD-prevEMA)+prevEMA
            DIV = MACD - EMA

            TAdata[ticker][period]['MACD']['EMAshort'] += [EMAshort]
            TAdata[ticker][period]['MACD']['EMAlong'] += [EMAlong]
            TAdata[ticker][period]['MACD']['MACD'] += [MACD]
            TAdata[ticker][period]['MACD']['EMA'] += [EMA]
            TAdata[ticker][period]['MACD']['DIV'] += [DIV]

        return TAdata


    def RSI(self,ticker,period,TAdata):

        TAdata[ticker]['daily']['RSI'] = {'gain':[], 'loss':[], 'RSI':[]}

        RSIperiod = 14
        
        for j in range(RSIperiod):
            TAdata[ticker][period]['RSI']['gain'] += ['N/A']
            TAdata[ticker][period]['RSI']['loss'] += ['N/A']
            TAdata[ticker][period]['RSI']['RSI'] += ['N/A']

        gain = 0.
        loss = 0.
        for j in range(1,RSIperiod+1):
            close = float(TAdata[ticker][period]['price']['close'][j])
            close_prev = float(TAdata[ticker][period]['price']['close'][j-1])
            change = close - close_prev
            if change > 0:
                gain += change
            else:
                loss += change
        if gain-loss == 0:
            RSI = 0.5
        else:
            RSI = gain/(gain-loss)
        TAdata[ticker][period]['RSI']['gain'] += [gain]
        TAdata[ticker][period]['RSI']['loss'] += [loss]
        TAdata[ticker][period]['RSI']['RSI'] += [RSI]

        for j in range(RSIperiod+1,len(TAdata[ticker][period]['price']['close'])):
            close = float(TAdata[ticker][period]['price']['close'][j])
            close_prev = float(TAdata[ticker][period]['price']['close'][j-1])
            change = close - close_prev
            if change > 0:
                gain = (gain*(RSIperiod-1)+change)/RSIperiod
            else:
                loss = (loss*(RSIperiod-1)+change)/RSIperiod
            if gain-loss == 0:
                RSI = 0.5
            else:
                RSI = 100*gain/(gain-loss)
            TAdata[ticker][period]['RSI']['gain'] += [gain]
            TAdata[ticker][period]['RSI']['loss'] += [loss]
            TAdata[ticker][period]['RSI']['RSI'] += [RSI]

        return TAdata


    def MFI(self,ticker,period,TAdata,):

        TAdata[ticker]['daily']['MFI'] = {'moneyflow':[], 'moneyflow_index':[]}
        
        MFIperiod = 14
        for j in range(MFIperiod):
            high = float(TAdata[ticker][period]['price']['high'][j])
            low = float(TAdata[ticker][period]['price']['low'][j])
            close = float(TAdata[ticker][period]['price']['close'][j])
            volume = float(TAdata[ticker][period]['price']['volume'][j])
            TAdata[ticker][period]['MFI']['moneyflow'] += [volume*(high+low+close)/3.]
        for j in range(MFIperiod,len(TAdata[ticker][period]['price']['close'])):
            high = float(TAdata[ticker][period]['price']['high'][j])
            low = float(TAdata[ticker][period]['price']['low'][j])
            close = float(TAdata[ticker][period]['price']['close'][j])
            volume = float(TAdata[ticker][period]['price']['volume'][j])
            TAdata[ticker][period]['MFI']['moneyflow'] += [volume*(high+low+close)/3.]
            moneyflow_positive = 0
            moneyflow_negative = 0
            for k in range(MFIperiod):
                moneyflow = TAdata[ticker][period]['MFI']['moneyflow'][(j-1)-k]
                moneyflow_prev = TAdata[ticker][period]['MFI']['moneyflow'][(j-1)-k-1]
                if moneyflow > moneyflow_prev:
                    moneyflow_positive += moneyflow
                else:
                    moneyflow_negative += moneyflow
            if moneyflow_positive+moneyflow_negative == 0:
                moneyflow_index = 0
            else:
                moneyflow_index = 100.*moneyflow_positive/(moneyflow_positive+moneyflow_negative)
            TAdata[ticker][period]['MFI']['moneyflow_index'] += [moneyflow_index]

        return TAdata
    

    def diagonalization(self,matrix):

        print(('diagonalizing %sx%s matrix' %(len(matrix),len(matrix))))

        import numpy

        eigenvalues,eigenvectors = numpy.linalg.eig(matrix)
        ## transpose
        eigenvectors = numpy.transpose(eigenvectors)
        ## organize eigenvalues and eigenvectors in list
        eigen_list = list(zip(eigenvalues, eigenvectors))
        ## sort list
        eigen_list.sort()
        ## reverse list
        eigen_list.reverse()
        ## parse sorted eigenvalues and eigenvectors
        eigenvalues = [eigen_list[eigen][0] for eigen in range(len(eigen_list))]
        eigenvectors = [eigen_list[eigen][1] for eigen in range(len(eigen_list))]

        return eigenvalues, eigenvectors


    def covar_matrix(self,d_pca):

        import numpy

        for ticker in list(d_pca.keys()):
            print((ticker, len(d_pca[ticker])))
            for i in range(len(d_pca[ticker])):
                d_pca[ticker][i] /= sum(d_pca[ticker])

        l_tickers = list(d_pca.keys())
        N = len(l_tickers)
        matrix = numpy.zeros((N,N))
        for i in range(N):
            l1 = d_pca[l_tickers[i]]
            n = len(l1)
            print((l_tickers[i]))
            average1 = sum(l1)/n
            for j in range(i,N):
                l2 = d_pca[l_tickers[j]]
                if n != len(l2):
                    raise 'data lists are of different length'
                average2 = sum(l2)/n

                SS = 0
                for k in range(n):
                    SS += (l1[k]-average1)*(l2[k]-average2)
                covar = SS/(n-1)
                matrix[i][j] = covar
                matrix[j][i] = covar

        eigenvalues, eigenvectors = numpy.linalg.eig(matrix)
       
        return matrix


    def plot_charts(self,TAcandidates, data, l_time):

        import os, math

        year2 = l_time[0] ; month2 = l_time[1] ; day2 = l_time[2]
        year1 = year2-1 ; month1 = month2  ; day1 = day2

        for ticker in TAcandidates:

            lines = []
            xtics = ''
            for j in range(len(data[ticker]['daily']['price']['date'])):
                date = data[ticker]['daily']['price']['date'][j]
                if date == '%s/%s/%s' %(day1, str(month1).zfill(2), str(year1)[-2:]):
                    break

            for i in range(j):
                date = data[ticker]['daily']['price']['date'][j]
                op = float(data[ticker]['daily']['price']['open'][j]); cl = float(data[ticker]['daily']['price']['close'][j])
                hi = float(data[ticker]['daily']['price']['high'][j]); lo = float(data[ticker]['daily']['price']['low'][j])
                line = '%4i %s %s %s %s %s' %(j-i, date, op, hi, lo, cl)
                MACD = data[ticker]['daily']['MACD']['MACD'][i]; EMA = data[ticker]['daily']['MACD']['EMA'][i]; diff = data[ticker]['daily']['MACD']['DIV'][i]
                line += ' %6.3f %6.3f %6.3f\n' %(MACD, EMA, diff)
                lines.append(line)
                if i-int(i/28)*28 == 0:
                    xtics += '"%s" %s, ' %(date, j-i)
                if date == '%s/%s/%s' %(day1, str(month1).zfill(2), str(year1)[-2:]):
                    break

            fd = open('finance.dat', 'w')
            fd.writelines(lines)
            fd.close()

            lines = [
                'set terminal postscript eps enhanced color "Helvetica" 48\n',
                'set output "%s.ps"\n' %(ticker),
    #            'set boxwidth 0.5\n'
                'set style fill solid 1.0\n'
                'set size 4,4\n',
                'set xtics(%s)\n' %(xtics[:-2]),
    #            'set style data fsteps\n',
    #            'set timefmt "%d/%m/%y"\n',
    #            'set xdata time\n',
    #            'set xrange [ "%s/%s/%s":"%s/%s/%s" ]\n' %(day1,int(month1)-3,str(int(str(year2)[-2:])-0).zfill(2),day2,month2,str(year2)[-2:]),
    #            'set format x "%d/%m/%y"\n',
                'set grid\n',
                'plot "finance.dat" using 1:3:5:4:6 with candlesticks title "%s"' %(ticker), # date, open, low, high, close
                ',"finance.dat" using 1:9 with impulses title ""', # date, diff
                ',"finance.dat" using 1:7 with line title "MACD"', # date, macd
                ',"finance.dat" using 1:8 with line title "EMA9"\n', # date, ema9
                ]

            fd = open('gnuplot.script', 'w')
            fd.writelines(lines)
            fd.close()

            os.system('gnuplot gnuplot.script')
            os.system('convert %s.ps %s.png' %(ticker, ticker))
            os.remove('%s.ps' %(ticker))

if __name__=='__main__':
    instance_TA = TA()

    s = 'AKS,SWN,A,AA,AAPL,ABC,ABK,ABT,ACE,ACS,ADBE,ADI,ADM,ADSK,AEE,AEP,AES,AET,AFL,AGN,AIG,AIV,AIZ,AKAM,ALL,ALTR,AMAT,AMD,AMGN,AMP,AMT,AMZN,AN,ANF,AOC,APA,APC,APD,APOL,ATI,AVB,AVP,AVY,AXP,AYE,AZO,BA,BAC,BAX,BBBY,BBT,BBY,BC,BCR,BDK,BDX,BEN,BHI,BIG,BIIB,BJS,BK,BLL,BMC,BMS,BMY,BNI,BRCM,BSX,BTU,BXP,C,CA,CAG,CAH,CAM,CAT,CB,CBG,CBS,CCE,CCL,CEG,CELG,CHK,CHRW,CI,CIEN,CINF,CL,CLX,CMA,CMCSA,CME,CMI,CMS,CNP,CNX,COF,COH,COL,COP,COST,CPB,CPWR,CSC,CSCO,CSX,CTAS,CTL,CTSH,CTXS,CVG,CVH,CVS,CVX,D,DD,NDAQ,DE,DELL,DF,DFS,DGX,DHI,DHR,DIS,DOV,DOW,DRI,DTE,DTV,DUK,DVN,DYN,EBAY,ECL,ED,EFX,EIX,EK,EL,EMC,EMN,EMR,EOG,EP,EQR,ERTS,ESRX,ESV,ETFC,ETN,ETR,EXC,EXPD,EXPE,F,FCX,FDO,FDX,FE,FHN,FII,FIS,FISV,FITB,FLR,FO,FPL,FRX,GAS,GCI,GD,GE,GENZ,GILD,GIS,GLW,GME,GNW,GOOG,GPC,GPS,GR,GS,GT,GWW,HAL,HAR,HAS,HBAN,HCBK,HCP,HD,HES,HIG,HNZ,HOG,HON,HOT,HPQ,HRB,HSP,HST,HSY,HUM,IBM,ICE,IFF,IGT,INTC,INTU,IP,IPG,ITT,ITW,JAVA,JBL,JCI,JCP,JDSU,JEC,JNJ,JNPR,JNS,JPM,JWN,K,KBH,KEY,KFT,KG,KIM,KLAC,KMB,KO,KR,KSS,LEG,LEN,LH,LLL,LLTC,LLY,LM,LMT,LNC,LOW,LSI,LTD,L,LUK,LUV,LXK,M,MAR,MAS,MAT,MBI,MCD,MCHP,MCK,MCO,MDP,MDT,MET,MHP,MHS,MI,MIL,MKC,MMC,MMM,MO,MOLX,MON,MOT,MRK,MRO,MS,MSFT,MTB,WEC,MU,MUR,MWV,MYL,NBL,NBR,NEM,NI,NKE,NOC,NOV,NOVL,NSC,NSM,NTAP,NTRS,NUE,NVDA,NVLS,NWL,NWSa,NYT,NYX,ODP,OMC,OMX,ORCL,OXY,PAYX,PBG,PBI,PCAR,PCG,PCL,PCP,PDCO,PEG,PEP,PFE,PFG,PG,PGN,PGR,PH,PHM,PKI,PLD,PLL,PM,PNC,PNW,POM,PPG,PPL,PRU,PSA,PTV,PX,Q,QCOM,QLGC,R,RAI,RDC,RF,RHI,RL,ROK,RRC,RRD,RSH,RTN,RX,S,PXD,SBUX,SCHW,SE,SEE,SGP,SHLD,SHW,SIAL,SII,SLB,SLE,SLM,SNA,SNDK,SO,SPG,SPLS,SRE,SSP,STI,STJ,STR,STT,STZ,SUN,SVU,SWK,SWY,SYK,SYMC,SYY,T,TAP,TDC,TE,TEG,TER,SJM,TGT,THC,TIE,TIF,TJX,TLAB,TMK,TMO,TROW,TRV,TSN,TSO,TSS,TWX,TXN,TXT,UNH,UNM,UNP,UPS,USB,UTX,VAR,VFC,VIAb,VLO,VMC,VNO,VRSN,VZ,WAG,WAT,APH,WFC,WFMI,WFR,WHR,WIN,WLP,FLS,WMB,WMI,WMT,WPI,WPO,WU,DPS,WY,WYN,X,XEL,XL,XLNX,XOM,XRX,XTO,YHOO,YUM,ZION,ZMH,DVA,FTR,CRM,FAST,ADP,PBCT,XRAY,WYNN,CEPH,MWW,BFb,BFa,SRCL,LIFE,DNB,RSG,EQT,MFE,OI,FLIR,SCG,IRM,HCN,DO,VTR,HRL,NU,ORLY,TWC,DNR,FTI,DV,PCS,PWR,WDC,RHT,FMC,CFN,ARG'
    l_tickers = s.split(',')
    l_tickers = ['AA','INTC','MSFT','AIG','AXP','BA','C','CAT','DD','DIS','F','GE','GM','HD','HPQ','IBM','JNJ','JPM','KO','MCD','MMM','MRK','PFE','PG','UTX','VZ','WMT','XOM','BAC','CVX',]
    l_tickers = ['TLKM.IJ',]
    d_ADR = {'TLK':['TLKM.IJ']}
    l_tickers = ['NEU',]

    import time
    l_time = time.localtime()
    months = {'Jan': 1, 'Feb': 2, 'Mar': 3, 'Apr': 4, 'May': 5, 'Jun': 6, 'Jul': 7, 'Aug': 8, 'Sep': 9, 'Oct': 10, 'Nov': 11, 'Dec': 12}
    d_msn2yahoo = {'GB':'L','FR':'PA','DE':'DE','AU':'AX','ES':'MC','JP':'','IT':'MI','SE':'ST','BE':'BR','CA':'TO','NL':'AS',}
    l_statementNA = []
    d_portfolio = {}

    instance_TA.find_candidates_TA(
        l_tickers, l_time, months, l_statementNA, d_portfolio,
        d_ADR,
        )
