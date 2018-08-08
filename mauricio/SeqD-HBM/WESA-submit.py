#/usr/bin/env python

import sys
import MultipartPostHandler, urllib2, cookielib

cookies = cookielib.CookieJar()
opener = urllib2.build_opener(urllib2.HTTPCookieProcessor(cookies),
                                MultipartPostHandler.MultipartPostHandler)
lines=open(sys.argv[3]).readlines()
seq=''
for line in lines:
	seq=seq+line[0:]

params = { "submitter" : sys.argv[1],\
	   "emailAddr" : sys.argv[2],\
           "userInput" : seq,\
	  }
print opener.open("http://pipe.sc.fsu.edu/cgi-bin/wesa/run", params).read()

