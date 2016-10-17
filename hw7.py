# -*- coding: utf-8 -*-
"""
Created on Mon Nov 16 13:54:27 2015

@author: itanveer
"""

from collections import defaultdict as ddict
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

def gibbs(prbList):
    x = np.cumsum(prbList)
    y = np.random.rand(1)
    z = np.argmax(x>y)
    return z

def train(eng,fr,prtable,debug=False):
    ec_num = ddict(lambda:1e-6)
    ec_den = ddict(lambda:1e-6)
    # E step
    if debug:
        print 'E Step'
    for eSent,fSent in zip(eng,fr):
        # split the words and insert NULL word
        e = eSent.strip().split(' ')
        e.insert(0,'NULL')
        f = fSent.strip().split(' ')
        l = len(e)
        m = len(f)
        # Gibb's sampling with 5 samples
        K = 5
        # t represents the sample number
        for t in xrange(K):
            for j in xrange(m):
                # sample the jth alignment using Gibb's sampler
                prbList = np.zeros(l)
                for i in xrange(l):
                    # for IBM model-1, alignments are independent. So,
                    # in Eq 7, the translation probabilities t(f1|ea1)t(f2|ea2)... will be
                    # cancelled in the numerator and denominator. so Eq 7 reduces down to:
                    # p(aj|a1, ...,aj-1,aj+1, ... aJ, F, E) = t(fj|eaj)/sum_aj t(fj|eaj)
                    prbList[i]=prtable[f[j],e[i]]
                prbList = prbList/np.sum(prbList)
                aj_sample = gibbs(prbList)
                # adding the strength for a single sample
                #ec_num[f[j],e[aj_sample]]+=(1.0/K)
                #ec_den[e[aj_sample]]+=(1.0/K)
                ec_num[f[j],e[aj_sample]]+=(1.0/K)
                ec_den[e[aj_sample]]+=(1.0/K)
                
    # M Step
    if debug:
        print 'M Step'
    for fj,ei in prtable.keys():
        # calculate p(fj|ei) table from expected counts
        prtable[fj,ei] = ec_num[fj,ei]/ec_den[ei]
        if debug:
            print ei,fj,ec_num[fj,ei],prtable[fj,ei]
    return prtable

def calclikelihood(prtable,eng,fr):
    L = 0.0
    for eSent,fSent in zip(eng,fr):
        # split the words and insert NULL word
        e = eSent.strip().split(' ')
        e.insert(0,'NULL')
        f = fSent.strip().split(' ')
        l = len(e)
        m = len(f)
        for j in xrange(m):
            k = 0.0
            for i in xrange(l):
                k+=prtable[f[j],e[i]]/l
            L+=np.log(k)
    return np.float(L)

def readfile(filename):
    sentList=[]
    with open(filename,'r') as f:
        for aline in f:
            sentList.append(aline.strip())
    return sentList

# Run this module to check the sample results against the given test case
def sampleresults():
    prtable = ddict(lambda: 1e-16)
    for iter in xrange(10):
        prtable = train(['a b','a c'],['A B','A C'],prtable,True)
#        xx=0
#        for akey in prtable.keys():
#            if akey
#            prtable[akey]
        print
        print 'Iteration #',iter
        print '================'
        print 'train loglikelihood',calclikelihood(prtable,['a b','a c'],['A B','A C'])
        print 'test loglikelihood',calclikelihood(prtable,['b c'],['B C'])
        print

# Run this module to apply the algorithm on a small training and test data
def runondata(train_eng,train_fra,test_eng,test_fra,output,iterNum):
    prtable = ddict(lambda: 1e-6)
    # read training data
    trainlist_eng = readfile(train_eng)
    trainlist_fra = readfile(train_fra)
    # read test data
    testlist_eng = readfile(test_eng)
    testlist_fra = readfile(test_fra)    

    plt.figure()    
    # iter
    maxloglikeli=-1*np.inf
    for iter in xrange(iterNum):
        print 'Iteration #',iter
        prtable = train(trainlist_eng,trainlist_fra,prtable)
        trainll = calclikelihood(prtable,trainlist_eng,trainlist_fra)
        testll = calclikelihood(prtable,testlist_eng,testlist_fra)
        print 'train_loglikelihood',trainll,'test_loglikelihood',testll
        
        if testll>maxloglikeli:
            maxloglikeli = testll
            maxptable = prtable.copy()
        
        # visualization
        if iter>0:
            plt.subplot(211)
            plt.scatter(iter,trainll,c='r')
            plt.hold(True)
            plt.xlabel('Iteration')
            plt.ylabel('Log Likelihood (Train)')
            
            plt.subplot(212)
            plt.scatter(iter,testll,c='b')
            plt.hold(True)
            plt.xlabel('Iteration')
            plt.ylabel('Log Likelihood (Test)')
            
            if train_eng=='training_short.eng':
                plt.suptitle('IBM Model 1 Log Likelihood for a small dataset (~100 sentences)')
            else:
                plt.suptitle('IBM Model 1 Log Likelihood for full data')
            plt.draw()
            plt.pause(0.01)
            plt.savefig(output,dpi=300)
            plt.pause(0.01)
    return maxptable

def decode(ptable,E,F):
    a = np.zeros(len(F))
    prb = []
    for indx,f in enumerate(F):
        for e in E:
            prb.append(ptable[f,e])
        a[indx]=np.argmax(np.array(prb))
    return a

def readprtable(filename):
    f = open(filename,'r')
    prtable=ddict(lambda:1e-16)
    for line in f:
        x = line.strip().split(',')
        if x[0]!=' ' and x[1]!=' ' and len(x)==3:
            prtable[x[0].strip(),x[1].strip()]=float(x[2].strip())
    return prtable

if __name__=='__main__':
    # Check the results against given test case
#    sampleresults()
    # Run the IBM model on a small data
#    ptable = runondata('training_short.eng','training_short.fra','test_short.eng','test_short.fra',\
#        'result_small.png',5)
    ptable = runondata('training.eng','training.fra','test.eng','test.fra',\
        'result_full.png',5)
    
    #ptable = readprtable('/Users/itanveer/devel/ibm-Model1-ttable/ttable')

    fe = open('aer/test.e','r')
    ff = open('aer/test.f','r')
    fo = open('aer/test.o','w')

    # Calculate AER
    for idx in xrange(447):
        le = fe.readline().strip().split(' ')[2:-2]
        lf = ff.readline().strip().split(' ')[2:-2]
        a = decode(ptable,le,lf)
        e=ddict(list)
        for x in xrange(len(a)):
            if a[x]!=0:
                e[a[x]].append(x+1)
            #e[a[x]].append(x)
        print>>fo,'begin',idx+1
        for anitem in e.items():
            print>>fo,int(anitem[0]),
            for vals in anitem[1]:
                print>>fo,vals,
            print>>fo,' '
            fo.flush()
        
#    runondata('training.eng','training.fra','test.eng','test.fra','result_Full.png',20)
