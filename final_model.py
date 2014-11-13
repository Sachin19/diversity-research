import matplotlib.pyplot as plt
import numpy as ny
import math
import os
from math import sqrt
import time
import operator
import random
import sys
from collections import deque

def pearsoncor(X, Y, code = 0):
    """
    Computes pearson's correlation coefficient.
    code
       0 - using deviations from means.
       1 - common formula.
    """
    n = len(X)
    sx = ssd(X)
    sy = ssd(Y)
    xbar = float(sum(X)) / n
    ybar = float(sum(Y)) / n
    if code == 0:
       return sum([(x - xbar) * (y - ybar) for x, y in zip (X,Y)])/(sx * sy*(n-1.0))
    else:
       numer = sum([x*y for x,y in zip(X,Y)]) - n*(xbar * ybar)
       denom = sqrt((sum([x*x for x in X]) - n* xbar**2)*(sum([y*y for y in Y]) -n* ybar**2))
       return (numer /denom)

def svar(X):
	n = len(X)
	if n <= 1:
	   raise ValueError, "sd(): n must be greater than 1"
	xbar = float(sum(X)) /n
	return (sum([(x-xbar)**2 for x in X])/(n-1))

def ssd(X):
	return sqrt(svar(X))
	
def fieldhieararchy():
	readfile = open('data/field_hierarchy.txt','r')
	d = {}
	for line in readfile:
		line = line.strip()
		l = line.split()
		d[l[0]] = l[1]
	readfile.close()
	return d
	
def getcitesdict():
	files = os.listdir('data/domains_with_cites')
	d = {}
	k = 0
	for file in files:
		readfile = open('data/domains_with_cites/'+file,'r')
		for line in readfile:
			l = line.split()
			cites = l[-1]
			l = l[:-1]
			if len(l) < 2:
				k = k + 1
			else:
				if l[-1][0] == '(': l = l[:-1]
				line = ' '.join(l)
			line = line.strip()
			line = line.lower()
			d[line] = cites
		readfile.close()
	return d

def authres(level=1):
	readfile = open('data/dataset_alltagged1','r')
	thispaper = []
	dyear = {}
	mscites = {}
	bestcoauth = {}
	dauthresl ={}
	mscites = {}
	din = {}
	imp = {}
	dcite = getcitesdict()
	dfield = fieldhieararchy()
	fielddict = {}	#fielddict[1234] = 'algo'
	dauths = {}		#dauths[1234] = ['jack','richard'...]
	authglob = {}
	firstauth = {}
	for line in readfile:
		line = line.strip()
		if line == '':
			auths,fields,cites = [],[],[]
			for a in thispaper:
				if a[1] == 't':
					year = int(a[2:])
				elif a[1] == 'i':
					index = a[6:]
				elif a[1] == '%':
					cites.append(a[2:])
				elif a[1] == '@':
					auths = a[2:].split(',')
				elif a[1] == 'f':
					if level == 2:
						fields.append(dfield[a[2:]])
					else: fields.append(a[2:])
				elif a[1] == '*':
					title = a[2:]
					if title[-1] == '.':
						title = title[:-1]
					title = title.lower()
			if len(fields) == 1: 
				fielddict[index] = fields[0]
				dauths[index] = auths
				dyear[index] = year
				for auth in auths:
					if auth != '':
						auth = auth.replace(' ','_')
						if auth in authglob:
							authglob[auth].append(index)
						else: authglob[auth] = [index]
				if title in dcite:
					mscites[index] = int(dcite[title])
				else:
					mscites[index] = 'not_in_ms'
				for cite in cites:
					if cite not in din:
						din[cite] = [index]
					else: din[cite].append(index)
				for auth in auths:
					if auth != '':
						auth = auth.replace(' ','_')
						if auth not in dauthresl:
							dauthresl[auth] = []
						for f in fields:
							dauthresl[auth].append((f,index,year,1))
			thispaper = []
			continue
		thispaper.append(line)	
	for i in dyear:
		sum = 0
		if i in din:
			for y in din[i]:
				sum += 1
		imp[i] = sum
	
	readfile.close()
	return (dyear,fielddict,dauths,imp,mscites,authglob,dauthresl)

def binspapers(n,level=1):
	(dyear,fielddict,dauths,imp,mscites,authglob,dauthresl) = authres(level)
	bins = []
	for i in range(n):
		bins.append([])
	temp = []
	for i in dyear:
		temp.append((i,dyear[i]))
	temp = sorted(temp,key=byyear)
	k = len(temp)/n
	k1 = len(temp) % n
	for i in range(k+k1):
		bins[0].append(temp[i][0])
	j = k+k1
	for i in range(1,n):
		for m in range(k):
			bins[i].append(temp[j][0])
			j += 1
	return (dyear,fielddict,dauths,imp,mscites,authglob,bins,dauthresl)

def byyear(t):
	return int(t[1])

def getfields(level=1):
	d = fieldhieararchy()
	if level == 2:
		return list(set(d.values()))
	return d.keys()

def normalizepref(ovals):
	sum = 0
	for i in ovals:
		sum = sum + i
	for i in range(len(ovals)):
		ovals[i] = float(ovals[i])/sum
	cumovals = ny.cumsum(ovals)
	return cumovals
    
def prefattach(circle):
	r = random.random()
	for i in range(len(circle)):
		if i == 0 and r < circle[0]:
			return i
		elif i > 0 and circle[i-1] <= r and r < circle[i]:
			if i+1 <= len(circle)-1 and circle[i] == circle[i+1]:
				return random.choice([i,i+1])
			else:  return i

def getpreffield(d):
	x,y = [],[]
	for i in d:
		x.append(i)
		y.append(d[i])
	y = normalizepref(y)
	i = prefattach(y)
	return x[i]

def updateauthdict(mauth,f,auths):
	for auth in auths:
		if auth in mauth:
			d = mauth[auth]
			if f in d:
				d[f] += 1
			else: d[f] = 1
			mauth[auth] = d
		else:
			d = {}
			d[f] = 1
			mauth[auth] = d
	return mauth

def updatetotcitedict(mauthc,ncites,f,auths):
	for auth in auths:
		if auth in mauthc:
			mauthc[auth] += ncites
		else:
			mauthc[auth] = ncites
	return mauthc

def basic4(n,k0,rand,level=1):			#weightage to authors who got more total cites
	(dyear,fielddict,dauths,imp,mscites,authglob,bins,dauthresl) = binspapers(n,level)
	writefile = open('data/authresearch_1','w')
	allf = getfields(level)
	mfield = {}
	mauth = {}
	mauthc = {}
	count,count1 = 0,0
	mc,mc1 = 0,0
	for i in range(k0):
		for paper in bins[i]:
			f = fielddict[paper]
			mfield[paper] = f
			if mscites[paper] == 'not_in_ms':	ncites = imp[paper]
			else: ncites = max(mscites[paper],imp[paper])
			mauth = updateauthdict(mauth,f,dauths[paper])
			mauthc = updatetotcitedict(mauthc,ncites,f,dauths[paper])
			count1 += 1
	random.seed()
	for i in range(k0,n):
		for paper in bins[i]:
			r = random.random()
			if r >= rand:
				d,d2 = {},{}
				den,c1 = 0,-1
				maxauth = ''
				for auth in dauths[paper]:
					if auth in mauthc:
						if mauthc[auth] > c1:
							maxauth = auth
							c1 = mauthc[auth]
						den += mauthc[auth]
				if den > 0:
					for auth in dauths[paper]:
						if auth in mauth:
							d1 = mauth[auth]
							num = mauthc[auth]
							for j in d1:
								if j in d: 
									d[j] += d1[j]*(float(num)/den)
									d2[j] += d1[j]
								else:
									d[j] = d1[j]*(float(num)/den)
									d2[j] = d1[j]
				if len(d) > 0:
					f = getpreffield(d)
					realf = fielddict[paper]
					dex = mauth[maxauth]
					maxf = max(dex.iteritems(), key=operator.itemgetter(1))[0]
					if maxf == realf: mc += 1
					else: mc1 += 1
					count += 1
			if  r < rand or len(d) == 0:
				r1 = random.randint(0,len(allf)-1)
				f = allf[r1]
			mfield[paper] = f
			mauth = updateauthdict(mauth,f,dauths[paper])
			if mscites[paper] == 'not_in_ms':	ncites = imp[paper]
			else: ncites = mscites[paper]
			mauthc = updatetotcitedict(mauthc,ncites,f,dauths[paper])
	#print str(mc),str(mc1),str(float(mc)/(mc+mc1))
	for auth in authglob:
		writefile.write(auth+' ')
		for index in authglob[auth]:
			writefile.write(mfield[index]+'@'+str(dyear[index])+'@'+str(imp[index])+'@'+str(mscites[index])+' ')
		writefile.write('\n')
	writefile.close()
	return (mscites,dauthresl,imp,mfield,fielddict)

def degreedist(l,x,flag):
	b = 10
	y = []
	buck = float(max(l)-min(l))/b
	i = min(l)
	if flag: x.append(i+buck)
	y.append(0)
	for j in range(1,b):
		if flag: x.append(x[j-1]+buck)
		y.append(0)
	for j in l:
		for k in range(b):
			if k == 0 and j < x[k]:
				y[k] += 1
				break
			if x[k-1] <= j < x[k] :
				y[k] += 1
				break
			if k == b-1 and j == x[k]:
				y[k] += 1
	x1 = [float(x[0])/2]
	for i in range(1,b):
		x1.append(float(x[i]+x[i-1])/2)
	s = sum(y)
	for i in range(len(y)):
		y[i] = float(y[i])/s
	return (x,x1,y)

def plotlines1(l,l1,flag,wl):
	x,x1,y = degreedist(l,[],True)
	x2,xp1,yp1 = degreedist(l1,[],True)
	return (xp1,yp1,pearsoncor(y,yp1))
	
def bothdiv2(n,k0,rand,wl,mw,enttype,level=1):
	w = 2 #for pearsoncorr
	(mscites,dauthresl,imp,mfield,fielddict) = basic4(n,k0,rand,level)  #for basic2 and basic5 5 args
	if enttype == 'windowentropy':
		dres1 = windowentropy(wl,1,level)
		dres2 = windowentropy(wl,0,level)
	else:
		dres1 = plaindiv(w,1,level)
		dres2 = plaindiv(w,0,level)
	
	(xp1,yp1,corv1) = plotlines1(dres1.values(),dres2.values(),1,wl)
	#print str(xp1)
	#print str(yp1)
	return (xp1,yp1,corv1)

def bothdiv(enttype):
	nl = [100]
	randl = [0]
	wl = [5]
	mw = [5]
	level = 1
	its = 1
	le = -1
	x,y = [],[]
	cov = 0
	
	for j in range(its):
		print str(j)
		xp1,yp1,corv1 = bothdiv2(100,10,0,5,5,enttype,1)
		print str(xp1)
		if len(x) == 0:
			for i in range(len(xp1)):
				x.append(xp1[i])
				y.append(yp1[i])
			cov = corv1
			le = len(x)
		else:
			for i in range(le):
				x[i] = float((j*x[i])+xp1[i])/(j+1)
				y[i] = float((j*y[i])+yp1[i])/(j+1)
			cov = float((j*cov)+corv1)/(j+1)
		
	plotlines2(x,y,enttype)
	print 'Covariance: ' + str(cov)
	return

def plot2text(fname,x,y):
	writefile = open('output_'+fname,'a')
	for i in range(len(x)):
		writefile.write(str(x[i])+' ')
	writefile.write('\n')
	for i in range(len(y)):
		writefile.write(str(y[i])+' ')
	writefile.write('\n\n')
	writefile.close()
	return

def plotlines2(xm,ym,enttype):
	if enttype == 'windowentropy':
		dres1 = windowentropy(5,1,1)	#for window entropy
	else:
		dres1 = plaindiv(2,1,1)		#for plain entropy
	
	x,x1,y = degreedist(dres1.values(),[],True)
	
	plt.figure()
	plt.plot(x1,y,'g-',xm,ym,'r-')
	plot2text(enttype+'.txt',x1,y)
	plot2text(enttype+'.txt',xm,ym)
	plt.legend(['real-data','model-data'])
	plt.ylabel('Fraction of Authors')
	plt.xlabel(enttype)
	plt.show()
	return

def caldiversity(d,flag = 0):
	sum = 0
	for i in d.values():
		if flag == 0: sum += len(i)
		else:
			if i > 0: sum += i
	res = 0
	for i in d.values():
		if flag == 0:f = float(len(i))/sum
		elif i > 0:
			f = float(i)/sum
			res += f*math.log(f)
	if res != 0: res = -1*res
	return res

def getresdict(flag,level=1):
	if flag == 1:
		if level == 1: readfile = open('data/authresearch_data','r')
	else:
		readfile = open('data/authresearch_1','r')
		
	d = {}
	for line in readfile:
		line = line.strip()
		line = line.split()
		if line[0] != '':
			d[line[0]] = []
			for s in line[1:]:
				ls = s.split('@')
				d[line[0]].append(tuple(ls))
		ls = d[line[0]]
		ls = sorted(ls,key=byyear)
		d[line[0]] = ls
	readfile.close()
	return d

def plaindiv(w,flag,level=1):
	d = getresdict(flag,level)
	res,y = [],[]
	gres = {}
	for auth in d:
		dres = {}
		prev = 'x'
		sum,k = 0,0
		if len(d[auth]) < w: continue
		for a,b,c,x in d[auth]:
			if a in dres: dres[a] += 1
			else: dres[a] = 1 
		res.append(caldiversity(dres,1))
		gres[auth] = caldiversity(dres,1)
	return gres

def windowentropy(w,flag,level=1):
	d = getresdict(flag,level)
	resl = []
	dres = {}
	for auth in d:
		res,k = 0,0
		if len(d[auth]) < w: continue
		temp = {}
		for a,b,c,x in d[auth][:w]:
			if a in temp: temp[a] += 1
			else: temp[a] = 1
		res = caldiversity(temp,1)
		k += 1
		l,r = 0,w-1
		while r < len(d[auth]):
			lv,rv = d[auth][l][0],d[auth][r][0]
			temp[lv] -= 1
			if rv in temp:
				temp[rv] += 1
			else: temp[rv] = 1
			res += caldiversity(temp,1)
			k += 1
			l += 1
			r += 1
		res = float(res)/k
		resl.append(res)
		dres[auth] = res
	#drawhist(resl)
	return dres

def main(enttype):
	if (enttype == 'windowentropy') or (enttype == 'plainentropy'):
		bothdiv(enttype)
		return
	print 'Invalid argument. Valid arg are windowentropy or plainentropy'

if __name__ == '__main__':
	if	len(sys.argv) > 1:
		main(sys.argv[1])
	else: print 'Enter an argument. Valid arg are windowentropy or plainentropy'