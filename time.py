import time, datetime

time_str = datetime.date(2012,1,1,)

time_tuple = time.gmtime()

time_seconds = time.mktime(time_tuple)

time_seconds = os.path.getmtime(fn)

week_number = int(time.strftime('%U',time_tuple,))

time_tuple = time.strptime("%s %s %s" %(year,month,day,), "%Y %m %d")

time_tuple = time.localtime(time_seconds)

time_str = time_str-datetime.timedelta(days=int(42))

'''
time() -- return current time in seconds since the Epoch as a float
clock() -- return CPU time since process start as a float
sleep() -- delay for a number of seconds given as a float
gmtime() -- convert seconds since Epoch to UTC tuple
localtime() -- convert seconds since Epoch to local time tuple
asctime() -- convert time tuple to string
ctime() -- convert time in seconds to string
mktime() -- convert local time tuple to seconds since Epoch
strftime() -- convert time tuple to string according to format specification
strptime() -- parse string to time tuple according to format specification
tzset() -- change the local timezone
'''
