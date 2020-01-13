# -*- coding: utf-8 -*-
"""
Created on Sun Dec 29 16:37:35 2019

This script parses aliqueit logs in order to obtain a scatter plot of the 
time spent factoring as a function of the base-10 logarithm. Two log/linear fits
are done to determine the crossover point between SIQS/GNFS.

contact: patrick.f.poitras [at] gmail.com

Requirements:
Python 3, with numpy, scipy and matplotlib.
Anaconda with Python 3 has all those included and more. I recommend using that.

(This code might also work in Python 2, but I haven't tested it. You might have to
modify it a bit.)

Licensing: Do what you want with it, but if you modify the code, at least mention
that it's based on or inspired by this work. (CC-BY 4.0)
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize

#For saving the figure
SAVE_FIGURE = False
FILENAME = "graph_zoomed.png"
DPI = 600

#For zooming on the crossover point
ZOOM = False
ZOOM_X_LIMITS = (94,99)
ZOOM_Y_LIMITS = (1000,3600)

#Graph parameters
XTICK_POSITION = [80,85,90,95,100,105,110,115,120,125,130,135,140]
XTICK_LABELS = [80,None,90,None,100,None,110,None,120,None,130,None,140]
YTICK_LABELS = ["1m","2m",None,None,"5m","10m","20m","30m",None,None,"1h","2h","3h",None,"5h","10h","1d","2d",None,"4d"]
YTICK_POSITION = 60*np.array([1,2,3,4,5,10,20,30,40,50,60,120,180,4*60,5*60,10*60,24*60,48*60,24*60*3,24*60*4])

#Functions for the curve fits
def fitfunc_gnfs(x,a,b,c,d):
    return a*np.exp((x-c)*b) + d

def fitfunc_siqs(x,a,b,c,d,e,f,g,h):
    return a*np.exp((x-c)*b) + e*np.log10(np.abs(x*f)) + d 

#These are the original parameter guesses for the curve fit, it must have
#the same number of elements as the number of parameters (not counting x) in the
# fit function above
p0_gnfs = (1.38,0.03,0,0)
p0_siqs = (1.38,0.04,0,0,1,1,60,1)

#################Code start#######################
f = open("aliqueit.log")
L = []
for line in f:
    L.append(line) 
f.close()

events = []
for i in L:
    if "Starting factorization of" in i:
        events.append(("Start",i))
    elif "prp" in i and not "***" in i:
        events.append(("prp", i))
    elif " pm1: starting B1 =" in i:
        events.append(("ECM-start",i))
    elif "final ECM pretested depth" in i:
        events.append(("ECM-end",i))
    elif "nfs: commencing nfs" in i:
        events.append(("NFS-start",i))
    elif "NFS elapsed time" in i:
        events.append(("NFS-end",i))
    elif "Total factoring time" in i:
        events.append(("factoring-end",i))
    elif "scheduler: switching to sieve method" in i:
        events.append(("sieve-start",i))
    elif "starting SIQS on" in i:
        events.append(("SIQS-start",i))
    elif "SIQS elapsed time" in i:
        events.append(("SIQS-end",i))
     
u = []
v = []
for i in events:
    if i[0] == "NFS-start" or i[0] == "NFS-end":
        u.append(i)
    elif i[0] == "SIQS-start" or i[0] == "SIQS-end":
        v.append(i)

start = None
timestamps=[]
errors = 0
for i in u:
    if i[0] == "NFS-start":
        start = i
    elif i[0] == "NFS-end" and start != None:
        #length parsing
        startlen = start[1].split("commencing nfs on")[1].split(":")[1].replace("\n","").replace(" ","")
        startlenf = float(startlen)
        logstartlenf = np.log10(startlenf)
        #time parsing
        time = float(i[1].split("= ")[1].replace(" seconds.\n",""))
        logtime = np.log10(time)
        #sending to timestamps
        timestamps.append((logstartlenf,time))
    else:
        errors += 1
print("Errors: ",errors)

    
siqs_ts = []
start  = None
for i in v:
    if i[0] == "SIQS-start":
        start = i
    elif i[0] == "SIQS-end" and start != None:
        startlenf = float(start[1].split("SIQS on c")[1].split(": ")[1])
        logstartlenf = np.log10(startlenf)
        #time parsing
        time = float(i[1].split("SIQS elapsed time = ")[1].replace(" seconds.\n",""))
        logtime = np.log10(time)
        siqs_ts.append((logstartlenf,time))
    else:
        errors += 1

print("NFS points:", len(timestamps))
print("SIQS points:", len(siqs_ts), " ({0} >= 80)".format(len([i for i in siqs_ts if i[0]>=80])))

#xx = np.linspace(92,1.01*max([i[0] for i in timestamps]),500)
xx = np.linspace(92,140,1000)
h = scipy.optimize.curve_fit(fitfunc_gnfs,np.array([i[0] for i in timestamps]),np.array([i[1] for i in timestamps]),p0=p0_gnfs)

l = scipy.optimize.curve_fit(fitfunc_siqs,np.array([i[0] for i in siqs_ts if i[0] > 60]), np.array([i[1] for i in siqs_ts if i[0] > 60]), p0=p0_siqs,maxfev=10000)
plt.scatter([i[0] for i in timestamps],[i[1] for i in timestamps],color='red',alpha=0.3)
plt.scatter([i[0] for i in siqs_ts], [i[1] for i in siqs_ts], color='blue',alpha=0.3)
plt.plot(xx,fitfunc_gnfs(xx,*h[0]),color='green',linestyle='--')
xx = np.linspace(80,100,500)
plt.plot(xx,fitfunc_siqs(xx,*l[0]),color='magenta',linestyle='--')
plt.grid()
plt.yscale("log",subsy=3600*np.array([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]))
plt.yticks(YTICK_POSITION,YTICK_LABELS)
#plt.yticks(60*np.array([1,2,3,4,5,10,20,30,40,50,60,120,180,4*60,5*60,10*60,24*60]),["1m","2m",None,None,"5m","10m","20m","30m",None,None,"1h","2h","3h",None,"5h","10h","1d"])
plt.xlabel("log10(n)")

if not ZOOM:
    plt.xticks(XTICK_POSITION, XTICK_LABELS)
    plt.xlim(80,None)
    plt.ylim(0,None)
else:
    plt.xlim(*ZOOM_X_LIMITS)
    plt.ylim(*ZOOM_Y_LIMITS)
    
if SAVE_FIGURE:
    plt.savefig(FILENAME,dpi=DPI)
