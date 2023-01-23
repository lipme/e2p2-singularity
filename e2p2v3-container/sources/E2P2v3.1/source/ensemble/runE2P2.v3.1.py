#!/usr/bin/python

"""
Version      3.1
Name:        runE2P2.py
Original Date:        20161115
Authors:     Chuan Wang, Lee Chae
             Plant Metabolic Network
             Department of Plant Biology
             Carnegie Institution for Science
             Stanford, CA 94305
Contributor: Bo Xue

Description: Runs the Ensemble Enzyme Prediction Pipeline (E2P2) on a set of input protein sequences,
             outputting enzyme functional annotations in the forms of EC numbers or MetaCyc reaction
             IDs for any predicted enzyme sequences.

Usage:       runE2P2.py -i <input file of sequences> -o <output filename> -r <run directory>

"""

from multiprocessing import Process
from operator import itemgetter
import sys
import os
import shutil
import subprocess
import re
import time
import datetime

# Retrieve E2P2 directory path
script_path = os.path.abspath(__file__)
e2p2_path = os.path.dirname(script_path)


# Set up application path during runtime and import modules.
sys.path.insert(0, os.path.join(e2p2_path, 'source', 'ensemble'))
import prog
import ensemble
import refinepf

# Define classes and functions.
class Prediction:
    def __init__(self, id):
        self.id = id
        self.predictions = {}
        self.weights = {}

class FinalPredictions:
    def __init__(self,id):
        self.id = id
        self.predictions = []

class Weight:
    def __init__(self, name):
        self.name = name
        self.parameter = parameter
        self.classes = {}

class Classifier:
    def __init__(self, id):
        self.id = id
        self.predictions = {}
        self.weights = {}

def run_process(cmd):
    # Makes the actual system call.

    ret = subprocess.call(cmd, stdout=open('/dev/null', 'w'), stderr=subprocess.STDOUT)
    # ret = subprocess.call(cmd, stderr=subprocess.STDOUT)

def create_process(cmd):
    # Launches the process by sending the command arguments to the classify subroutine.
    # Returns a Process object that can then be queried for process status.
    p = Process(target=run_process, args=(cmd,))
    p.start()
    return p

def mkdirp(directory):
    if not os.path.isdir(directory):
        os.mkdir(directory)


# Still Testing
def handle_spaces_in_paths(cmd):
    return [r'"%s"' % c if ' ' in c else r'%s' % c for c in cmd]


# Assemble help message for the program using the get_help object 
# found in the prog module.
name = '''
    runE2P2.py
    '''
description = '''
    Runs the Ensemble Enzyme Prediction Pipeline (E2P2) on a set of input protein sequences,
    outputting enzyme functional annotations in the forms of EC numbers or MetaCyc reaction
    IDs for any predicted enzyme sequences.
    '''
options = '''
    -h --Displays this help message.
    -i --Name of the file containing the input protein sequences.
    -o --Name for the output file. [/tmp/E2P2v3.out]
    -r --Run directory [/tmp]
    -e --evalue cutoff [1e-5]
    -t --Number of threads (CPUs) to use in the BLAST search [1]
    '''
usage = '''
    runE2P2.py -i <input file of sequences> -o <output filename>
    '''
notes = '''
    - Input protein sequences should be in FASTA format.
    - Headers in the FASTA file should begin with the sequence ID followed by a space.
    - Intermediate results files can be found in the run/ directory in its own subdirectory labeled with a
      date and time stamp.
'''
message = prog.get_help(name, description, options, usage, notes)

# Collect command line options using get_options in prog.
flags = 'hi:o:r:e:t:'
args = sys.argv[1:]
options = prog.get_options(args, flags)

# Check for help request.
prog.check_help(options, message)

# Store the command line arguments in useful form.
## Edit: 9/16/16 Every Input Has a Seperate Folder for intermediate files
#if "-i" not in [a[0] for a in options[:]] or "-o" not in [a[0] for a in options[:]] or "-r" not in [a[0] for a in options[:]] :
if "-i" not in [a[0] for a in options[:]] :
    print "Please Specify Input"
    sys.exit()

threads = "1"
evaluecutoff = 1e-5
rundir="/tmp"
filename_output="/tmp/E2P2v3.out"

for a in options[:]:
    if a[0] == "-i":
        filename_input = os.path.abspath(a[1])
    if a[0] == "-o":
        if len(a[1]) < 1 or not os.path.isdir(os.path.dirname(os.path.abspath(a[1]))):
            print "Output Path Invalid: %s." % (a[1])
            sys.exit()
        filename_output = os.path.abspath(a[1])
    if a[0] == "-r":
        rundir = os.path.abspath(a[1])
    if a[0] == "-e":
        evaluecutoff = a[1]
    if a[0] == "-t":
        threads = a[1]
    

# Record date and time.
now = datetime.datetime.now()
timestamp = time.time()
time_stamp = str(timestamp)


# Read in the sequence IDs. Note that we do not need to load the actual sequence
# data into memory.
print "Reading input data."
sequences = {}
try:
    input = open(filename_input, 'r')
    for i in input:
        if i.startswith(">"):
            header = i.rstrip().lstrip(">").split("|")[0]
            id = re.split("\s+", header)[0]
            sequences[id] = 1
    input.close()
except IOError:
    print "Can't find the input file: %s" % (filename_input)
    sys.exit()

# Create the classifier objects.
classifiers = {}
cnames = ["BLAST", "CatFam", "Priam"]
for cn in cnames:
    c = Classifier(cn)
    c.predictions = {}
    c.weights = {}
    classifiers[cn] = c

# Populate the weight data for each classifier.
data = ""
try:
    input = open(os.path.join(e2p2_path, "source", "ensemble", "data", "weights"), 'r')
    for i in input:
        data += i
    input.close()

    # Populate temporary data structures to process weight data for each method.
    temp = data.split(">")
    for t in temp:
        if t:
            entries = t.split("\n")

            # Record weights.
            cname, parameter = entries[0].split("|")[0], entries[0].split("|")[1]
            c = classifiers[cname]
            for e in entries[1:]:
                try:
                    cid, wval = e.split("\t")[0], e.split("\t")[1]
                    if wval == "NA" or wval == "0" or wval == "-1.000":
                        wval = "0.000"
                    c.weights[cid] = wval
                except:
                    continue
except IOError:
    print "Can't find the following weight file: source/ensemble/data/weights."
    sys.exit()

del classifiers["CatFam"]

# Read in the mapping data that relates Enzyme Functional class numbers to EC classes and reaction IDs.
fc_map = {}
try:
    input = open(os.path.join(e2p2_path, "source", "ensemble", "data", "fcmap"), 'r')
    # create a dictionary mapping FC IDs to EC and reaction IDs.
    for i in input:
        if "#" not in i:
            try:
                ef, id = i.split("\t")[0], i.rstrip().split("\t")[1]
                fc_map[ef] = id
            except:
                continue
    input.close()
except IOError:
    print "Can't find the following functional class mapping file: source/ensemble/data/fcmap."
    sys.exit()

## Process the input file with each level-0 classifier. Run the classifiers concurrently as
## separate processes to save time.
print "Running level-0 classification processes."

## Edit: 9/16/16 Every Input Has a Seperate Folder for intermediate files
input_run_folder = os.path.join(rundir, 'run', os.path.basename(filename_input) + '.' + time_stamp)
if not os.path.exists(input_run_folder):
    os.makedirs(input_run_folder)
output_blast = os.path.join(input_run_folder, "blast." + time_stamp)

touch_cmd = handle_spaces_in_paths(['touch', output_blast])
touch_ret = run_process(touch_cmd)

blast_cmd = handle_spaces_in_paths([os.path.join(e2p2_path, 'source', 'blast', 'ncbi-blast-2.2.30+', 'bin', 'blastp'), '-db', os.path.join(e2p2_path, 'source', 'blast', 'db', 'rpsd-3.1.fa'), '-query', filename_input, '-out', output_blast, '-outfmt', '6', '-num_threads', threads])
#print(blast_cmd)
blast = create_process(blast_cmd)

output_priam = os.path.join(input_run_folder, "PRIAM_%s" % (time_stamp), "ANNOTATION", "sequenceECs.txt")
## Edit: 9/16/16 Add Memory Settings for Java
priam_cmd = handle_spaces_in_paths([os.path.join(e2p2_path, 'source', 'java', 'jre1.6.0_30', 'bin', 'java'), '-Xms3072m', '-Xmx3072m', '-jar', os.path.join(e2p2_path, 'source', 'priam', 'PRIAM_search.jar'), '--bd', os.path.join(e2p2_path, 'source', 'blast', 'blast-2.2.26', 'bin'), '-n', time_stamp, '-i', filename_input, '-p', os.path.join(e2p2_path, 'source', 'priam', 'profiles'), '--bh', '-o', input_run_folder, '--np', threads])
priam = create_process(priam_cmd)

# Hold until the last classifier finishes.
blast.is_alive()
priam.is_alive()
while blast.is_alive() or priam.is_alive():
    time.sleep(5)

## Process the output files from each classifer.
predictions = {}
print "Compiling predictions."

# Blast
temp_preds = {}
temp_hits = {}
preds_blast = {}
input = open(output_blast, 'r')
temp = ""
for i in input:
    temp += i
input.close()
data = temp.split("\n")

# Iterate over all queries.
for d in data:
    if d:
        # Get the ID.
        temp = d.split("\t")
        ##  Edit: 9/16/16 Make Sure Python Can Read the Format of e-value
        if temp[-2].startswith('e'):
            temp_float = "1" + temp[-2]
        else:
            temp_float = temp[-2]
        qid, temp_hit, eval = temp[0].split("\s+")[0].split("|")[0], temp[1], float(temp[-2])

        if eval > float(evaluecutoff):
            continue

        # Retrieve or create prediction object for a given query.
        if qid in temp_hits:
            continue
        else:
            temp_hits[qid] = eval

        if qid in temp_preds:
            tp = temp_preds[qid]
        else:
            tp = Prediction(qid)

        # Get the EFs from the hit line.
        hits = []
        ## Edit: 9/16/16 Split Query ID by "|" and whitespace
        #hits_temp = temp_hit.split("|")
        hits_temp = re.split("[|\s]+", temp_hit)
        for h in hits_temp[1:]:
            if "EF" in h:
                hits.append(h)
        for h in hits: # Skip the first field, which contains the hit ID. We only need the predicted EF class.
            if h in tp.predictions: # Record the lowest evalue for that EF class.
                best_eval = tp.predictions[h]
                if eval < best_eval:
                    tp.predictions[h] = eval
            else:
                tp.predictions[h] = eval
        temp_preds[qid] = tp

# Iterate over queries to find and record their top hits.
c = classifiers["BLAST"]
for qid in temp_preds:
    c.predictions[qid] = []
    p = Prediction(qid)
    tp = temp_preds[qid]

    # Find the lowest evalue.
    best_eval = 1.0
    for hit in tp.predictions:
        if tp.predictions[hit] < best_eval:
            best_eval = tp.predictions[hit]

    # Iterate over predictions to find the top hits.
    for hit in tp.predictions:
        if tp.predictions[hit] == best_eval:
            c.predictions[qid].append(hit)
            if hit in p.predictions:
                p.predictions[hit] = tp.predictions[hit]
            else:
                p.predictions[hit] = tp.predictions[hit]
    
    preds_blast[qid] = p            
predictions["BLAST"] = preds_blast

# Priam
#preds_priam = {}
c = classifiers["Priam"]
input = open(output_priam, 'r')
temp = ""
for i in input:
    temp += i
input.close()
data = temp.split(">")

# Iterate over all queries.
for d in data[1:]:

    # Process the hit info for the query.
    temp = d.split("\n")
    hits = []
    
    # Get the query ID.
    #qid = temp[0].split("\s+")[0].split("|")[0]
    
    ## Edit: 9/16/16 Split Query ID by "|" and whitespace
    qid = re.split("[|\s]+", temp[0])[0]
    
    c.predictions[qid] = []
    p = Prediction(qid)
    for t in temp[1:]:
        # Get the hits.
        if t.startswith('EF'):
            hits.append(t.split("\t")[0].rstrip())

    c.predictions[qid] = []
    #p = Prediction(qid)
    for h in hits:
        c.predictions[qid].append(h)
        #p.predictions[h] = 1.0
    #preds_priam[qid] = p
input.close()

# Calculate the ensemble prediction for each query sequence.
print "Computing ensemble predictions."
threshold = float(0.5)
final_predictions = {}
for qid in sequences:

    # Call the ensemble subroutine from the ensemble module.
    fpred = ensemble.perform_max_weight_absolute_threshold(qid, classifiers, threshold)
    final_predictions[qid] = fpred

## Output results files.
print "Preparing results files."
# Assemble run information.
run_data = "# Run date, time:  %s\n\
# Ensemble method used:  %s\n" % (now, "Maximum weight with absolute threshold (0.5)")

# Prepare short version of results.
short_output = "%s" % (run_data)
for seq_id in final_predictions:
    fpred = final_predictions[seq_id]
    pred_output = ""
    for real_class in fpred.predictions:
        real_class_edited = real_class.split(" ")[0]
        pred_output += real_class_edited + "|"
    final_preds = pred_output.rstrip("|")
    short_output += seq_id + "\t" + final_preds + "\n"

# Print short output.
#filename_output = filename_output + "." + voting_scheme
output = open(filename_output, 'w')
output.write(short_output)
output.close()

# Prepare long version of results.
long_output = "%s" % (run_data)
for seq_id in final_predictions:
    fpred = final_predictions[seq_id]
    pred_output = ""
    for real_class in fpred.predictions:
        pred_output += real_class + "|"
    final_preds = pred_output.rstrip("|")
    long_output += ">" + seq_id + "\t" + final_preds + "\n"
    for cname in fpred.classifiers:
        classifier_pred_output = ""
        for real_class in fpred.classifiers[cname]:
            classifier_pred_output += real_class + "|"
        final_classifier_preds = classifier_pred_output.rstrip("|")
        if final_classifier_preds != "":
            long_output += cname + "\t" + final_classifier_preds + "\n"
    long_output += "\n"

# Print long output.
filename_output_full = filename_output + ".long"
output = open(filename_output_full, 'w')
output.write(long_output)
output.close()

# Prepare Pathologic input file.
# Read in input data and store sequence ID and predictions if it's a valid enzyme prediction.
entries = {}

input = open(filename_output, 'r')
for line in input:
    labels = []
    
    # Skip the commented information.
    if '#' in line or line.startswith('\n'):
        continue
    else:
        id = ""
        preds = []
        check = 0
        l = line.rstrip()
        if '\t' in l:
            data = l.split('\t')
            id = data[0].split('|')[0]
            if "EF" in data[1]:
                preds = data[1].split('|')
                check = 1
            if check == 1:
                e = Prediction(id)
                e.labels = []
                for p in preds:
                    if "EF" in p:
                        e.labels.append(p)
                entries[id] = e
                check = 0

# Create the output file, translating EF classes into ECs and reaction IDs.
results = ''
for id in entries:
    e = entries[id]
    results += "ID\t%s\nNAME\t%s\nPRODUCT-TYPE\tP\n" % (id, id)
    for l in e.labels:
        if "EF" in l:
            try:
                translated_reaction = fc_map[l]
                if "RXN" in translated_reaction:
                    results += "METACYC\t%s\n" % (translated_reaction)
                else:
                    results += "EC\t%s\n" % (translated_reaction)
            except:
                print "EF class %s assigned to %s not found.\n" % (translated_reaction, id)
    results += "//\n"

# Print pathologic input file.
filename_pathologic = filename_output + ".pf"
output = open(filename_pathologic, 'w')
output.write(results)
output.close()

# Print pathologic input file translated to reaction ids only.
filename_pathologic_orxn = filename_output + ".orxn.pf"
pf_cmd = ['perl', os.path.join(e2p2_path, 'tools', 'pf-EC-to-official-RXN.pl'), filename_pathologic]
pf = create_process(pf_cmd)

# Hold until the last classifier finishes.
pf.is_alive()
while pf.is_alive():
    time.sleep(1)
refinepf.remove_empty_from_pf(filename_pathologic_orxn)
# Notify user of completion and exit.
print "Operation complete."
print "Main results are in the file: %s" % filename_output
print "Detailed results are in the file: %s" % filename_output_full
print "To build PGDB, use .pf file: %s" % filename_pathologic_orxn
print "Intermediate files are in the directory: %s" % input_run_folder
sys.exit()




