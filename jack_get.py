import sys
sys.path.append("/home/trung/_qhe-library")
from subprocess import call, Popen, PIPE
import numpy as np
import FQH_states as FQH
import misc
import time
from itertools import product, combinations_with_replacement

import string
import random

def read_jack(root, debug=False):
	# To adapt to the changes in convention
	# If the root starts with 0, remove the 0 e.g. 001001001 --> 01001001 --> 1001001
	# If the root is in negative Lz sector (convention left most digit = -S), reverse it. e.g 10001001 --> 10010001
	jack = FQH.fqh_state(f"jacks/J_{root}", quiet=True)
	if jack.dim()==0:
		if debug: print(f"{root} not found")
		S = (len(root)-1)/2
		if misc.findLZ(int(root,2), S) > 0:
			root_read = root[::-1]
			reverse = True
			if debug: print(f"reading root {root_read} instead")
		else:
			root_read = root
			reverse = False
			if debug: print(f"reading {root}")
		if root_read[0]=="0":
			if debug: print(f"reading root {root_read[1:]} instead")
			jack = read_jack(root_read[1:])
			jack.basis = ["0"+x for x in jack.basis] # put back the zero
		else:
			jack = FQH.fqh_state(f"jacks/J_{root_read}")
		if reverse:
			jack.basis = [x[::-1] for x in jack.basis]

	if debug:
		if jack.dim()==0:
			print(f"WARNING: reading {root} outputs a null state")
		print("-----")
	return jack

def generate_jack(root, alpha,fname=None,debug=False): 
	# root must be in binary string format
	p = Popen(["./jack"], stdin = PIPE, shell = True)
	No = len(root)
	rootconfig = " ".join(root)
	# Use this to get the coefficients
	if fname==None:
		filename = "none"
	else:
		filename = fname[:4]

	p.communicate(input=bytes(f"{No} -1\n{rootconfig}\n{filename}\n1\n{alpha}\n",encoding = 'utf-8'))
	call(["mv", filename, f"jacks/J_{root}"])
	return


def get_jack_list(root_list, alpha, geom="none", placeholder=None):
	# special means special case where one quasihole is at center
	try:
		root_read = FQH.fqh_state(root_list)
	except FileNotFoundError:
		print("file not found")
		return
	jack_list = []
	if placeholder==None:
		placeholder = "".join(random.choices(string.ascii_uppercase + string.digits + string.ascii_lowercase, k=4))
	for root_bin in root_read.basis:
		jack_read = FQH.fqh_state(f"jacks/J_{root_bin}", quiet=True)
		if jack_read.dim()==0:
			generate_jack(root_bin, alpha,fname=placeholder)
			jack_read = FQH.fqh_state(f"jacks/J_{root_bin}")
		if geom=="disk":
			jack_read.disk_normalize()
		elif geom=="sphere":
			jack_read.sphere_normalize()
		jack_list.append(jack_read)
	return jack_list
"""
def gen_Jack_2qh((k,l), Ne):
	# Assume k<=l
	root = "1"*k+"0"+"1"*(l-k)+"0"
	return
"""

"""

if __name__ == "__main__":
	nqh = int(input("Generate jack list for 1 or 2 quasihole? "))
	Ne = int(input("Input N_e: "))
	ph = str(input("Placeholder name for jacks (max 4 characters): "))
	JL = get_jack_list(Ne,nqh, placeholder = ph)
	print("------")
	print(f"{Ne} electrons and {nqh} quasiholes")
	print(f"Number of jacks = {len(JL)}")
	#st = time.time()
	#a  = two_quasihole_state(4, -2.5, 2.5,JL)
	#print(f"{time.time()-st} seconds")
	#print("Plotting density")
	#a.plot_disk_density("density_4_2.pdf", ref=[-2.5,2.5])
	#print(f"Total: {time.time()-st} seconds")
	#print("------")

"""

if __name__ == "__main__":
	import os
	if not os.path.isdir("jacks"):
		call(["mkdir", "jacks"])
		
	print("Generate a list of Jack polynomials from a given list of root configurations")
	fname = str(input("Input file name: "))

	alpha = str(input("alpha parameter = num/den. Input num den (separated with space)"))
	jacks = get_jack_list(fname, alpha)
