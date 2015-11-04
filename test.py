'''
Created on 05.06.2014

@author: david
'''

#!/usr/bin/env python
# a bar plot with errorbars


from Helper import Helper



valueMatrix=[[4,8,4,3],[4,8,4,3],[4,8,4,3],[4,8,4,3]]
fileName="/home/david/Desktop/dink.png"
barNamesTuple=('a','b','c','d')
legendTuple=('w','x','w','z')


Helper.createBarplot(valueMatrix, fileName, barNamesTuple, legendTuple)