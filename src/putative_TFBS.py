import argparse
import sys
import glob
from pwm import *

parser = argparse.ArgumentParser(description="Find the index of the best window for a TFBS.")

parser.add_argument("-t", "--threshold", type=float, help='score threshold for pssm',action="store",dest='t', required=True)
parser.add_argument("-l", "--promotor-length",type=int, help='promotor sequence length',default=1000, action="store", dest='l')
parser.add_argument("-w", "--window-size",type=int, help='sliding window size', default=40, action="store", dest='w')
parser.add_argument("-m", "--pfm",type=str, help='jaspar matrix name',action="store",dest='m', required=True)
parser.add_argument("-s", "--window-threshold", type=float, help='window score threshold', action="store", dest='s', required=True)
parser.add_argument("-p", "--pseudocount",type=float, help='pseudo-count weight', action="store", dest='p')
parser.add_argument("mrna", nargs='*', help='list of mRNA id.')

args = parser.parse_args()

"""
print(args.mrna)
print(args.m)
print(args.p)
print(args.s)
print(args.t)
print(args.l)
print(args.w)
"""

tf_name , wsi = search_luncher(args.mrna,args.m,args.p,args.t,args.l,args.w,args.s)

if len(wsi)== 0:
	print("****WARNING**** : No TFBS found with this parameters, please choose different ones.")

for i in wsi :
	print("Window id : "+str(i)+" | TF : "+tf_name+" | Window pos : [%d,%d] | Window score : %f" %(wsi[i][0],wsi[i][1],wsi[i][2]))
	wi = wsi[i][3]
	for j in wi :
		print(""+str(i)+" "+tf_name+"  %s  %d  %f" %(j,wi[j][0],wi[j][1]))

