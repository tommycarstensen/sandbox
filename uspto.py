import urllib2
## ICL C07D
lines_url = urllib2.urlopen('http://patft.uspto.gov/netacgi/nph-Parser?Sect1=PTO2&Sect2=HITOFF&u=%2Fnetahtml%2FPTO%2Fsearch-adv.htm&r=0&p=1&f=S&l=50&Query=ICL%2FC07D&d=PTXT')
lines = lines_url.readlines()
## CCL 514/270
url = 'http://patft.uspto.gov/netacgi/nph-Parser?Sect1=PTO2&Sect2=HITOFF&u=%2Fnetahtml%2FPTO%2Fsearch-adv.htm&r=0&p=1&f=S&l=50&Query=CCL%2F514%2F270&d=PTXT'
## AN Pfizer
url = 'http://patft.uspto.gov/netacgi/nph-Parser?Sect1=PTO2&Sect2=HITOFF&u=%2Fnetahtml%2FPTO%2Fsearch-adv.htm&r=0&p=1&f=S&l=50&Query=AN%2FPfizer&d=PTXT'
url = 'http://patft.uspto.gov/netacgi/nph-Parser?Sect1=PTO2&Sect2=HITOFF&u=%2Fnetahtml%2FPTO%2Fsearch-adv.htm&r=0&f=S&l=50&d=PTXT&OS=AN%2FPfizer&RS=AN%2FPfizer&Query=AN%2FPfizer&TD=3760&Srch1=Pfizer.ASNM.&NextList2=Next+50+Hits'
## AN/Pfizer and Sildenafil and CCL/514/252.05
url = 'http://patft.uspto.gov/netacgi/nph-Parser?Sect1=PTO2&Sect2=HITOFF&u=%2Fnetahtml%2FPTO%2Fsearch-adv.htm&r=0&p=1&f=S&l=50&Query=AN%2FPfizer+and+Sildenafil+and+CCL%2F514%2F252.05&d=PTXT'
## PN/5250534
url = 'http://patft.uspto.gov/netacgi/nph-Parser?Sect1=PTO2&Sect2=HITOFF&u=%2Fnetahtml%2FPTO%2Fsearch-adv.htm&r=1&p=1&f=G&l=50&d=PTXT&S1=5250534.PN.&OS=PN/5250534&RS=PN/5250534'
print lines
