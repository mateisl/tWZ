from __future__ import division
#!/usr/bin/env python

import numpy as np

tW = 0.53
sW = 0.47


def calc_c3pl(cpWB, cpD):
    return -cpWB/tW - cpD/(4*tW*tW)

def calc_c1pl(cpWB, cpD):
    return -cpD/4

def calc_cpl(cpWB, cpD):
    return -cpD/2

def calc_cdiffpq(cpWB, cpD):
    return  cpWB/tW + cpD/(4*sW*sW) - cpD/6

def calc_c3pq(cpWB, cpD):
    return -cpWB/tW - cpD/(4*tW*tW)

def calc_cpu(cpWB, cpD):
    return  cpD/3

def calc_cpd(cpWB, cpD):
    return -cpD/6

def calc_cpD(cdiffpq,c3pq):
    return (-cdiffpq - c3pq)/(-1/(4*sW*sW) + 1/(4*tW*tW) + 1/6 )

def calc_cpWB(cdiffpq,c3pq):
    return tW*(-c3pq - 1/(4*tW*tW)*calc_cpD(cdiffpq,c3pq) )

cdiffpq_list = [0   ,  0]
c3pq_list    = [-0.1,  0.1]

for i, cdiffpq in enumerate(cdiffpq_list):
    c3pq = c3pq_list[i]

    cpWB = calc_cpWB(cdiffpq, c3pq)
    cpD = calc_cpD(cdiffpq, c3pq)
    c1pq = cdiffpq + c3pq

    print '----------------------'
    print 'cpWB   =', cpWB
    print 'cpD    =', cpD
    print 'c3pl   =', calc_c3pl(cpWB, cpD)
    print 'c1pl   =', calc_c1pl(cpWB, cpD)
    print 'cpl    =', calc_cpl(cpWB, cpD)
    print 'c(-)pq =', calc_cdiffpq(cpWB, cpD)
    print 'c1pq   =', c1pq
    print 'c3pq   =', calc_c3pq(cpWB, cpD)
    print 'cpu    =', calc_cpu(cpWB, cpD)
    print 'cpd    =', calc_cpd(cpWB, cpD)

for cpD in [-1, 0, 1]:
    for cpWB in [-0.472, 0, 0.472]:
        cdiff = calc_cdiffpq(cpWB, cpD)
        c3pq = calc_c3pq(cpWB, cpD)
        c1pq = cdiff + c3pq
        print 'c1pq   =', c1pq
        print 'c3pq   =', c3pq
        print 'cpWB   =', calc_cpWB(cdiff, c3pq)
        print 'cpD    =', calc_cpD(cdiff, c3pq)
        print '-----'
