#!/usr/bin/python
import sys
import hashlib

str="".join(sys.argv[1:])

m=hashlib.md5()
m.update(str)
hex=m.hexdigest()
print int(hex,16)%(2**30)
