#!/usr/bin/python
# Get_Annot.py - This will take in a large annotation file from Cypher and a file containing positions of interest (chr:pos), and spit out annotations for ONLY positions in the position file.
# Usage: python Get_Annot.py /PATH/TO/POSIT.file /PATH/TO/ANNOT.file /PATH/TO/OUTPUT.file
# By Regina Nguyen and Mollee Jain (TSRI - July 2013)
# Modified by Tristan Carland (TSRI - Jan 2014)
# Modified by Kristopher Standish (UCSD/JCVI - March 2014)

from sys import argv # standard packages to allow imput params
from sys import exit
import gzip
import datetime
# import exit function to stop script if there is a ValueError in input parameters

# Checks for the correct number of input params, displays usage info if needed
try:
    script, var1_file, annot_file, output, = argv
except ValueError:
    print "\nScript used to create subsets of variants from vcf files that correspond to each of the snps lists."
    exit("cmd_usage$> script  var1_file  annot_file  output\n")
# Input params:
# script - this script
# var1_file - file with positions of variants you want annotations for
# annot_file - file with annotations at a bunch of variant positions
# output - outut file containing annotations for positions in var1

# open Position file
if "gz" in var1_file:
  v1 = gzip.open(var1_file)
if "gz" not in var1_file:
  v1 = open(var1_file)
print "V1 open - %s" % datetime.datetime.now().time().isoformat()

# open Annotation file
if "gz" in annot_file:
  annot = gzip.open(annot_file)
if "gz" not in annot_file:
  annot = open(annot_file)
print "Annotation File Open - %s" % datetime.datetime.now().time().isoformat()

# skip all the comment lines up top on the Annotation file
annot_line = annot.next()
while annot_line.startswith("#"):
  annot_line = annot.next()
print "Comments Skipped - %s" % datetime.datetime.now().time().isoformat()

# make a dictionary out of the positions you're interested in
dictionary = { }
v1_line = v1.next()
for v1_line in v1:
  temp = v1_line.strip().split()
  entry = "%s:%s" % (temp[1], temp[2])
  dictionary[v1_line.strip()] = 1
print "Dictionary Made - %s" % datetime.datetime.now().time().isoformat()

# open a file to write to
OUT=open(output, 'w')

# go thru annotation file and pull out lines for positions in var_file1 (aka, in dictionary)
for annot_line in annot:
  temp = annot_line.split()[0:5]
  if temp[4]=="snp":
    entry = "%s:%s" % (temp[1][3:], temp[3])
#    print entry
    if entry in dictionary:
#      print entry
      OUT.write(annot_line)

OUT.close()
dictionary.clear()

print "All Done...I hope - %s" % datetime.datetime.now().time().isoformat()

print "Time to close everything"

v1.close()
v2.close()
