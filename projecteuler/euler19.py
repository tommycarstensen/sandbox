sundays = 0
for year in range(1900,2001):
    if not year%4:
        if not year%100:
            if not year%400:
                days = 366
            else:
                days = 365
        else:
            days = 366
    else:
        days = 365
    print(year,days)
##    31 - 365%7 366%7
##    stop
