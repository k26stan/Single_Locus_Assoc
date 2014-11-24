#!/usr/bin/python
# Get_Annot.py - This will take in a large annotation file from Cypher and a file containing positions of interest (chr:pos), and spit out annotations for ONLY positions in the position file.
# Usage: python Get_Annot.py /PATH/TO/POSIT.file /PATH/TO/ANNOT.file /PATH/TO/OUTPUT.file
# By Regina Nguyen and Mollee Jain (TSRI - July 2013)
# Modified by Tristan Carland (TSRI - Jan 2014)
# Modified by Kristopher Standish (UCSD/JCVI - March 2014)

#### HAVE TO ORGANIZE THE OUTPUTS THAT CHRIS WANTS ####
# Chromosome (.P file)
# Variant Position (.P file)
# rsID (.P file)
# Test P-Value (.P file)
# Expected P-Value (Calculate w/ R)
# HWE P-Value (.hwe file)
# Allele Frequency in Cases (NA)
# Allele Frequency in Ctrls (NA)
  # Genotype Frequency Instead (.hwe file)

from sys import argv # standard packages to allow imput params
from sys import exit
import gzip
import datetime
# import exit function to stop script if there is a ValueError in input parameters

# Checks for the correct number of input params, displays usage info if needed
try:
    script, var_file, hw_file, adj_file, frq_file, pv_file, hw_out, adj_out, frq_out, pv_out = argv
except ValueError:
    print "\nScript used to create subsets of variants from vcf files that correspond to each of the snps lists."
    exit("cmd_usage$> script, var_file, hw_file, adj_file, frq_file, pv_file, hw_out, adj_out, frq_out, pv_out\n")

## Open Position file
var = open(var_file)
print "Var_File open - %s" % datetime.datetime.now().time().isoformat()

## Open Hardy/Adjusted file
hw = open(hw_file)
adj = open(adj_file)
frq = open(frq_file)
pv = open(pv_file)
print "HW_File, Adj_File, Frq_File, & PV_File open - %s" % datetime.datetime.now().time().isoformat()

## Open files to write to
hw_write=open(hw_out, 'w')
adj_write=open(adj_out, 'w')
frq_write=open(frq_out, 'w')
pv_write=open(pv_out, 'w')

## Skip Header for all files
var_line = var.next()
hw_line = hw.next()
adj_line = adj.next()
frq_line = frq.next()
pv_line = pv.next()

## Write Headers to Outputs
hw_write.write(hw_line)
adj_write.write(adj_line)
frq_write.write(frq_line)
pv_write.write(pv_line)

## Make a dictionary out of the positions you're interested in
dictionary = { }
how_far = 0
for var_line in var:
  how_far += 1
  temp = var_line.strip().split()
  # entry = "%s:%s" % (temp[0], temp[2])
  entry = temp[1]
  # print var_line
  # print entry
  dictionary[entry] = 1
print "Dictionary Made - %s" % datetime.datetime.now().time().isoformat()
var.close()

## Go thru hw file and pull out lines for positions in var_file
print "Going Thru HW_File open - %s" % datetime.datetime.now().time().isoformat()
how_far = 0
for hw_line in hw:
  temp = hw_line.split()[0:8]
  entry = temp[1]
  how_far += 1
  # if how_far%1000000==0:
  #   print entry
  # entry = "%s:%s" % (temp[1][3:], temp[3])
  if entry in dictionary:
    # print entry
    split_line = hw_line.split()
    hw_write.write(hw_line)
print "HW_File done - %s" % datetime.datetime.now().time().isoformat()

## Go thru adj file and pull out lines for positions in var_file
print "Going Thru ADJ_File open - %s" % datetime.datetime.now().time().isoformat()
how_far = 0
for adj_line in adj:
  temp = adj_line.split()[0:8]
  entry = temp[1]
  how_far += 1
  # if how_far%1000000==0:
  #   print entry
  # entry = "%s:%s" % (temp[1][3:], temp[3])
  if entry in dictionary:
    # print entry
    split_line = adj_line.split()
    adj_write.write(adj_line)
print "ADJ_File done - %s" % datetime.datetime.now().time().isoformat()

## Go thru frq file and pull out lines for positions in var_file
print "Going Thru FRQ_File open - %s" % datetime.datetime.now().time().isoformat()
how_far = 0
for frq_line in frq:
  temp = frq_line.split()[0:8]
  entry = temp[1]
  how_far += 1
  # if how_far%1000000==0:
  #   print entry
  # entry = "%s:%s" % (temp[1][3:], temp[3])
  if entry in dictionary:
    # print entry
    split_line = frq_line.split()
    frq_write.write(frq_line)
print "FRQ_File done - %s" % datetime.datetime.now().time().isoformat()

## Go thru pv file and pull out lines for positions in var_file
print "Going Thru PV_File open - %s" % datetime.datetime.now().time().isoformat()
how_far = 0
for pv_line in pv:
  temp = pv_line.split()[0:8]
  entry = temp[1]
  how_far += 1
  # if how_far%1000000==0:
  #   print entry
  # entry = "%s:%s" % (temp[1][3:], temp[3])
  if entry in dictionary:
    # print entry
    split_line = pv_line.split()
    pv_write.write(pv_line)
print "PV_File done - %s" % datetime.datetime.now().time().isoformat()

hw_write.close()
adj_write.close()
frq_write.close()
pv_write.close()
dictionary.clear()
print "All Done...I hope - %s" % datetime.datetime.now().time().isoformat()

print "Time to close everything"

hw.close()
adj.close()
frq.close()
pv.close()