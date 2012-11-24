qid = 760000
qid = 781352
qid = 1559561 ## 9. maj 2010, time1
qid = 1559562 ## 9. maj 2010, time2
qid = 1338348 ## Pagten (2009) - Uden Hinanden - Maria Lucia

##
## 2) mms stream in old net player
##
url = 'mms://wms.dr.dk/nas01/auto/cms/Resources/dr.dk/P3/2010/12/f4221f06-53ff-4e96-a7b6-afce8c20c61f/Det_Elektriske_Barometer_201012052103_2300.wma?streamingsessionid=a08496b2-498e-6040-a3be-3f395abf%26playersessionid=1fdfd246-b0b0-eecd-a5c0-8275af7d'
url = 'mms://wms.dr.dk/nas01/auto/cms/Resources/dr.dk/NETRADIO/2010/05/b2b439b3-235b-4779-b8ce-2ee78191278d/Det%20Elektriske%20Barometer.wma?streamingsessionid=a08496b2-498e-6040-a3be-3f395abf%26playersessionid=1fdfd246-b0b0-eecd-a5c0-8275af7d'
if 'Det_Elektriske_Barometer' in url:
    index1 = url.index('Det_Elektriske_Barometer')
    index2 = url.index('.wma')
    fn = url[index1:index2-9]
else:
    fn = 'DEB'
s = 'mplayer "%s" -dumpfile /home/tc/Podcasts/%s.mp3' %(url,fn,)
print s

##
## 1) playlist (qid) in old net player
##
s = 'mplayer -playlist "http://www.dr.dk/Forms/Published/PlaylistGen.aspx?bitrate=low&qid=%s&location=udland&uri=http://dr.dk/Forms/Published/PlaylistGen.aspx" -dumpstream -dumpfile %s.mp3' %(qid,qid,)
print
print s
print
s = 'mplayer -playlist "http://www.dr.dk/Forms/Published/PlaylistGen.aspx?qid=%s" -dumpstream -dumpfile %s.mp3' %(qid,qid,)
print
print s
print

##
## 3) CooJah - new Flash Player (2011May)
##
playpath_fn = 'Det_Elektriske_Barometer_201105082103_2300_high'
playpath_dn = 'mp3:cms/Resources/dr.dk/P3/2011/05/19779d9c-a3cf-48c8-a588-a4f8811b2798'
tcUrl = 'rtmp://vod.dr.dk/cms'
pageUrl = 'http://www.dr.dk/LiveNetRadio/popup.html#id/593'
url = 'rtmp://vod.dr.dk/cms/mp3:cms/Resources/dr.dk/P3/2011/05/a909a85f-210b-4187-ba29-8fc824a1a007/Det_Elektriske_Barometer_201105152103_2300_high'
tcUrl = url[:url.index('/mp3')]
playpath_dn = url[url.index('mp3'):url.rindex('/')]
playpath_fn = url[url.rindex('/')+1:]
print 'rtmpdump',
print '-r "%s/%s/%s"' %(tcUrl,playpath_dn,playpath_fn),
##print '-a "cms"',
##print '-f "WIN 10,3,181,14"',
##print '-p "%s"' %(pageUrl),
##print '-s "http://www.dr.dk/LiveNetRadio/NetRadio.swf?ver=2.0.10"',
##print '-t "%s"' %(tcUrl),
print '-y "%s/%s"' %(playpath_dn,playpath_fn),
print '-o "%s.flv"' %(playpath_fn[:33])
