#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 15:17:41 2018

@author: stangen
"""

import sys
import json
data = sys.argv[1]
print(data)
data=json.loads(data)
#
#data = '{"T2M":1, "ALT":1}'
print(data)
print(type(data))
print(type(data['T2M']))
#print(json.loads(data))
#print(type(json.loads(data)))


#print(type(data))


