#!/usr/bin/env python

'''
Author: Pelle Scholten
Usage: python %prog -i in.fa.out -o out
Requirement: Python 3.6 or later
License: see LICENSE file
'''

import os,sys,gzip,datetime,collections,argparse,errno
import fuseTE, mergeTE, rcStatm, filter
import subprocess
import re

version='2021/12/03'

'''
Version log:
version='2024/06/28': first release
'''


description='''
    This script collapses TE annotations from RepeatMasker (i.e. genome.fa.out / genome.fa.align file).
    The TE annotations from RepeatMasker frequently overlap each other. This scripts finds the most likely
    TE annotation for every overlapping base and adjusts the annotations accordingly. This is done recursively 
    and efficiently using a sorting tree algorithm. Optionally the annotations are further processed (merging, filtering).
    This script only uses single thread.
    Requirement: Python 3.6 or later
'''

# args
parser=argparse.ArgumentParser(description=description)
parser.add_argument('-i', metavar='input_file', type=str, help='Specify input genome.fa.out file.', required=True)
parser.add_argument('-o', metavar='output_file', type=str, help='Specify output basename.', required=True)
parser.add_argument('-lvl', metavar='int', type=int, help='Optional. Specify percentage of overlap necessary to complete remove a repeat. Default: 80', default=80)
parser.add_argument('-min', metavar='int', type=int, help='Optional. Specify minimum length of basepairs a repeat must have after being cut. Default: 50', default=50)
parser.add_argument('-remove_simple_repeat', help='Optional. Specify if you want to rempve "Simple_repeat" and "Low_complexity".', action='store_true')
parser.add_argument('-mode', help='Optional. Specify if you want to generate out files for TE content extraction (TEcontent) or for consensus alignment (alignment). Default = TEcontent', type=str, default = 'TEcontent')
parser.add_argument("-mergemode",
                    help="Merge mode for defragmentation (if alignment mode is used). Use repeatmasker ID or a threshold, or both. Treshhold can be determined with -gapsize. Choose between 'ID', 'threshold' and 'both' Default = none",
                    default="none", type=str)
parser.add_argument('-remove', help='Optional. Call if short fragments should be removed', action='store_true')
parser.add_argument('-aligninput', help='Optional. Call if .align file input is used', action='store_true')
parser.add_argument('-famlength', help='Optional. Specify family length file for relative removal', default="none", type=str)
parser.add_argument('-gapsize', metavar='int', type=int, help='Optional. Specify gapsize for defragmentation. Default: 150', default=150)
parser.add_argument('-leave_overlap', help='Optional. Specify if you do not want to resolve overlapping repeat annotations.', action='store_true')
parser.add_argument('-v', '--version', action='version', version='Version: %s %s' % (os.path.basename(__file__), version))
parser.add_argument('-testrun', action='store_true', help=argparse.SUPPRESS)
args=parser.parse_args()

# set up
obed='%s.bed' % args.o
ogff='%s.gff' % args.o
olabel='%s.label.gff' % args.o
omerge='%s.merged.gff' % args.o
ofilter='%s.filter.gff' % args.o

if args.famlength != "none":
    lengthfile = args.famlength
    checklength = True
else:
    checklength = False

frequencyfilealignclass ='%s.freqs.alignment.class.tsv' % args.o
frequencyfilealignfamily ='%s.freqs.alignment.family.tsv' % args.o
frequencyfilecontentclass ='%s.freqs.TEcontent.class.tsv' % args.o
frequencyfilecontentfamily ='%s.freqs.TEcontent.family.tsv' % args.o
contentfile='%s.content.tsv' % args.o

gap=0
gapsize= args.gapsize
level = args.lvl / 100
threshold = args.min
annots=('gene', 'transcript', 'exon')
remove = args.remove
mode = args.mode
aligninput = args.aligninput
mergemode = args.mergemode
leave_overlap = args.leave_overlap

if mode == "alignment": 
    alignment = True
    tomerge = True
else:
    alignment = False

if mode != "TEcontent" and mode != "alignment":
    sys.stderr.write("Mode is not recognised\nIf you want to have out files for TE content, please specify \'-mode TEcontent\'\nIf you want to have out files for consensus alignment, please specify \'-mode alignment\'\n")
    sys.exit(1)
if mergemode != "none" and mergemode != "both" and mergemode != "threshold" and mergemode != "ID":
    sys.stderr.write("Mergemode is not recognised. Please use either \"ID\", \"threshold\", \"both\" or \"none\" (default = \"none\"). \n")
    sys.exit(1)

if mergemode == "none":
    tomerge=False
else:
    tomerge=True

if args.remove_simple_repeat is True:
    remove_simple_low_complex={'Simple_repeat', 'Low_complexity'}
else:
    remove_simple_low_complex={}

_date=datetime.datetime.now()

# check input and output
if os.path.exists(args.i) is False:
    print(FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.i))
    exit(1)
outdir=os.path.dirname(args.o)
if '/' in outdir:
    os.makedirs(outdir, exist_ok=True)

log_file = open("errors.log","w")

class IntervalTree:
    def __init__(self):
        self.root=None
        
    def insert(self, start, end, score, info, consensus_info):
        if self.root is None:
            self.root=IntervalNode(start, end, score, info, consensus_info)
        else:
            self.root=self.root.insert(start, end, score, info, consensus_info)
    
    def find_top_score(self, start, end):
        return self.root.find_top_score(start, end)
    
    def reassign_short_fragments(self, fragment, connected):
        return self.root.reassign_short_fragments(fragment, connected)


class IntervalNode:
    def __init__(self, start, end, score, info, consensus_info):
        self.start=start
        self.end=end
        self.score=score
        self.info=info
        self.consensus_info=consensus_info
        self.maxend=self.end
        self.minstart=self.start
        self.left=None
        self.right=None

    def replace_root(self,start,end,score,info, consensus_info):
        # replace root and corresponding branches with new node
        # while deleting old root from tree
        root = IntervalNode(start, end, score, info, consensus_info)
        if self.left:
            root.left = self.left
        if self.right:
            root.right = self.right
        return root

    def insert(self, start, end, score, info, consensus_info):
        root=self

        # Check whether there are overlaps of more than 80 (or other value determined in command line) percent shared sequence
        # If sequence is shared over the threshold than remove the sequence with lower score 
        # If not, keep sequence in tree to be resolved
        if start >= self.start and start <= self.end and end >= self.end:  # start>self.start unnecessary as file is ordered, if unordered file: would require code for start<self.start
            overlap = self.end - start

            # if root (almost) contained in new element and new element is better--> delete and replace root
            # root:              ------------     SW = 500       -->  
            # new element:         ------------      SW = 1000   -->      ------------      SW = 1000   
            if overlap > level * (self.end - self.start) and self.score < score:
                #sys.stderr.write(f"\r{self.info}\n")
                return self.replace_root(start,end,score,info,consensus_info)
            
            # if new element (almost) contained in root and root is better --> delete and replace root
            # root:              ------------     SW = 1000     -->       ------------      SW = 1000
            # new element:         ------------      SW = 500   -->      
            if overlap > level * (end - start) and self.score > score:
                return root
        
            # if an element is (almost) contained by another and score is the same --> keep longer element
            # root:              ------------       SW = 1000     -->       
            # new element:        --------------   SW = 1000     -->   ---------------   SW = 1000
            if (overlap > level * (end - start) or overlap > level * (self.end - self.start)) and self.score == score:
                if end - start >= self.end - self.start:
                    return self.replace_root(start,end,score,info, consensus_info)
                else:
                    return root
                
            # If none of these, keep sequence in tree to be resolved

        # Same as previous section but for chimers
        if start >= self.start and start <= self.end and end <= self.end:

            # If nested is worse or equal --> do not add to tree
            # nesting:      --------------     SW = 1000       -->   only NESTING remains     --------------      SW = 1000 
            # nested:        ---------         SW = 500        -->   
            if score <= self.score:
                return root
            
            # If nested is better and the nesting is for more than 80 (or other) percent contained in nested --> delete nesting
            # nesting:      --------------     SW = 500        -->
            # nested:        ---------         SW = 1000       -->   only NESTED remains         ---------        SW = 1000 
            elif (end - start)/(self.end - self.start) > level:
                root = IntervalNode(start, end, score, info, consensus_info)
                if self.left:
                    root.left = self.left
                if self.right:
                    root.right = self.right
                return root
            
            # If none of these, keep chimer in tree, to be resolved later...

        if start >= self.start:
            # insert to right tree
            if self.right: #if there is already a node the right of the root, repeat with that node
                self.right=self.right.insert(start, end, score, info, consensus_info)
            else: # add node to the right
                self.right=IntervalNode(start, end, score, info, consensus_info)
            # build heap

            # change root if score is higher
            if self.score < self.right.score:
                root=self.rotateleft()
        else: # not in ordered dataset
            print("Dataset non ordered at: " + self.start + " " + self.end + " " + self.score + " " + self.info + " please order dataset (or change the code)", outfile = log_file)
            sys.exit('\n Dataset non ordered, please order dataset (or change the code)\n')
            # insert to left tree
            if self.left:
                self.left=self.left.insert(start, end, score, info)
            else:
                self.left=IntervalNode(start, end, score, info)
            # build heap
            if self.score < self.left.score:
                root=self.rotateright()

        # for finding intervals
        # max end of all nodes in a branch, to know whether a segment is contained in a branch
        if root.right and root.left:
            root.maxend=max(root.end, root.right.maxend, root.left.maxend)
            root.minstart=min(root.start, root.right.minstart, root.left.minstart)
        elif root.right:
            root.maxend=max(root.end, root.right.maxend)
            root.minstart=min(root.start, root.right.minstart)
        elif root.left:
            root.maxend=max(root.end, root.left.maxend)
            root.minstart=min(root.start, root.left.minstart)
        return root
    
    def rotateright(self):
         # new node becomes root
        root=self.left
        # node on right becomes node on left
        # old root becomes node on right
        self.left=self.left.right
        root.right=self
        if self.right and self.left:
            self.maxend=max(self.end, self.right.maxend, self.left.maxend)
            self.minstart=min(self.start, self.right.minstart, self.left.minstart)
        elif self.right:
            self.maxend=min(self.end, self.right.maxend)
            self.minstart=min(self.start, self.right.minstart)
        elif self.left:
            self.maxend=max(self.end, self.left.maxend)
            self.minstart=min(self.start, self.left.minstart)
        return root

    def rotateleft(self):
        # new node becomes root
        root=self.right
        # node on left becomes node on right
        # old root becomes node on left
        self.right=self.right.left
        root.left=self
        if self.right and self.left:
            self.maxend=max(self.end, self.right.maxend, self.left.maxend)
            self.minstart=min(self.start, self.right.minstart, self.left.minstart)
        elif self.right:
            self.maxend=max(self.end, self.right.maxend)
            self.minstart=min(self.start, self.right.minstart)
        elif self.left:
            self.maxend=max(self.end, self.left.maxend)
            self.minstart=min(self.start, self.left.minstart)
        return root
    
    def find_top_score(self, start, end):
        # if segment is contained in current node
        if start < self.end and end > self.start:
            return [(self.score, self.info, self.consensus_info)]
        else:
        # look if segment is contained in node to the left or right (with a lower score)
        # by looking at highest end value (maxend) and lowest start value (minstart) of the branches left and right respectively.
            res=[]
            if self.left and start < self.left.maxend:
                res.extend(self.left.find_top_score(start, end))
            if self.right and end > self.right.minstart:
                res.extend(self.right.find_top_score(start, end))
            return res
    
    def reassign_short_fragments(self, fragment, connected):
        # if same element as the one it is currently assigned to
        if self.info == fragment.info or self.score >= fragment.score:
            res=[]
            if self.left and fragment.start < self.left.maxend:
                res.extend(self.left.reassign_short_fragments(fragment, connected))
            if self.right and fragment.end > self.right.minstart:
                res.extend(self.right.reassign_short_fragments(fragment, connected))   
            return res

        # get list of all elements in the current 'island'
        # to make sure 
        elements = []
        for element in connected:
            elements.append(element.info)

        # if fragment is contained in current node --> reassign fragment
        # but only part of fragment that is contained by the current node
        if fragment.start < self.end and fragment.end > self.start:# and self.info in elements:
                if fragment.end < self.end and fragment.start > self.start: #chimer
                    return [(fragment.start, fragment.end, self.score, self.info, self.consensus_info)]
                elif fragment.start > self.start: # non chimer where fragment more right than current node
                    if fragment.end - self.end > 0:
                        unresolved_fragment =[(self.end, fragment.end, fragment.score, fragment.info, self.consensus_info, "unresolved")]
                        return [(fragment.start, self.end, self.score, self.info, self.consensus_info)], unresolved_fragment
                    else:
                        return [(fragment.start, self.end, self.score, self.info, self.consensus_info)]
                else: # non chimer where fragment more left than current node
                    if fragment.start - self.start > 0:
                        unresolved_fragment = [(fragment.start, self.start, fragment.score, fragment.info, self.consensus_info, "unresolved")]
                        return [(self.start, fragment.end, self.score, self.info, self.consensus_info)], unresolved_fragment
                    else:
                        return [(self.start, fragment.end, self.score, self.info, self.consensus_info)]
        else:
        # look if segment is contained in node to the left or right (with a lower score)
        # by looking at highest end value (maxend) and lowest start value (minstart) of the branches left and right respectively.
            res=[]
            if self.left and fragment.start < self.left.maxend:
                res.extend(self.left.reassign_short_fragments(fragment, connected))
            if self.right and fragment.end > self.right.minstart:
                res.extend(self.right.reassign_short_fragments(fragment, connected))
            return res

def connect(l):
    connected=[]
    prev_e= (gap * -1) - 1
    prev_i=None
    prev_c=None
    prev_r=None
    for _rep in l:
        if _rep.end - _rep.start == 0:
            continue

        # check if elements can and should be merged
        dist= _rep.start - prev_e
        if _rep.info == prev_i and dist <= gap and _rep.score == prev_r: #re-connected the fragmented pieces
            prev_e=_rep.end
            prev_r = max(prev_r, _rep.score)
        else:
            if not prev_i is None: #do nothing if its first element
                connected.append(Rep(prev_s, prev_e, prev_r, prev_i, prev_c))
            prev_s,prev_e,prev_r,prev_i,prev_c=_rep #update which fragment is the previous one
    connected.append(Rep(prev_s, prev_e, prev_r, prev_i, prev_c))

    return connected

def connect_and_reassign(results,tmp,tree):

    # connect adjacent fragments with same info
    connected=connect(results)
    # filter out elements that become too short after fixing overlap
   
    connected_dummy = connected.copy()
    #sort from high score to low to make sure that:
    # 1000      ----------                       ----------                       ----------                               ----------
    #  600         ----------            -->               ---            -->                                and not
    #  300            ----------                              ---                           -------  
    #  400               ----------------                         ---------                        ---------                         ----------------
    connected_dummy.sort(key = lambda element: element[2], reverse=True)

    for fragment in connected_dummy:
        #if fragment has been cut...


        if fragment not in tmp and fragment.end - fragment.start < threshold: # too short

            # try to find new element for fragment
            outcome= tree.reassign_short_fragments(fragment, connected)
            # outcome is a list with fragments
            # there can be multiple hits for a fragment
            # there can also be unresolved fragments in the list, these always directly follow the fragment that are belonging too
            # for example
            # [hit 1, unresolved 1, hit 2, unresolved 2]
            # [hit 1, hit 2, unresolved 2]

            # sort unresolved and resolved items
            result =[]
            unresolved = []
            dummy_outcome = []

            if len(outcome) == 0:
                result = outcome
            else:
                for i in outcome:
                    if len(i) == 5:         # when only 1 fragment in outcome list, formatting isue
                        result.append(i)
                        dummy_outcome.append(i)
                    elif len(i) == 6:        # when only 1 fragment in outcome list , formatting isue          
                        unresolved.append(i)
                        dummy_outcome.append(i)
                    elif len(i[0]) == 5:     # when multiple fragments in outcome list
                        result.append(i[0])
                        dummy_outcome.append(i[0])
                    else:
                        unresolved.append(i[0])
                        dummy_outcome.append(i[0])

            connected.remove(fragment)

            if len(result) == 0: #if fragment does not have a matching element --> do not add back
                continue
            elif len(result) == 1: # if 1 hit
                result=result[0]
                if len(unresolved) > 0:
                    unresolved = unresolved[0]    
            else:   # if multiple hits
                # get result with highest score
                result.sort(key = lambda element: element[2])
                result = result[-1]

                #find location in outcome list
                location_result = dummy_outcome.index(result)
                location_potential_unresolved = location_result+1

                # if there are unresolved fragments in the list, and if there is a value after the location of the result
                # check whether this location is unresolved, if so assign to unresolved
                if len(unresolved) > 0 and len(dummy_outcome) > location_potential_unresolved:
                    if len(dummy_outcome[location_potential_unresolved]) == 5:
                        unresolved = []
                    if len(dummy_outcome[location_potential_unresolved]) == 6:
                        unresolved = dummy_outcome[location_potential_unresolved]
                else:
                    unresolved = []
                  
            # if succesful, add fragment with new info to results & repeat

            connected.append(Rep(*result))

            #check
            if len(unresolved) > 1 and len(unresolved) != 6:
                print("Something went wrong, more than one unresolved fragment found at: " + str(fragment) + " " + str(unresolved), file=sys.stderr) 
                print("Something went wrong, more than one unresolved fragment foundat: " + str(fragment) + " " + str(unresolved), file=log_file)
    
            if len(unresolved) > 0:
                connected.append(Rep(*unresolved[0:5])) # append without unresolved label

            connected.sort(key = lambda element: element[0])
            connected = connect_and_reassign(connected,tmp,tree)
            break

    return connected

def remerge(connected, tmp):
    # remerge fragments that have been split up by a chimer (for aligment only)
    for i in range(len(connected)):
        for j in range(i+1,len(connected)):
            if connected[i] not in tmp and connected[j] not in tmp: # do not deal with untouched chimers, they are fixed later
                if connected[i].info == connected[j].info and connected[i].score == connected[j].score: # check information of elements
                    connected[i] = Rep(connected[i].start, connected[j].end, connected[j].score, connected[j].info, connected[j].consensus_info)
                    connected.remove(connected[j])
                    connected = remerge(connected,tmp)
                    break
    return connected

                
def collapse(tmp, chr, gft_id_n):
    if len(tmp) == 1:
        connected=tmp
    else:
        tree=IntervalTree()
        poss=set()

        for _rep in tmp:
            tree.insert(*_rep) # insert repeat as node into the tree
            poss |= {*_rep[:2]}
          
        # get start and end positions of repeats in overlap island
        poss=sorted(list(poss))
        results=[]

        for n in range(len(poss) - 1):

            # get start and end position + score
            se = poss[n:n + 2]

            #find top score
            #find scores associated to fragment
            result=tree.find_top_score(*se)

            if len(result) == 0: #if fragment does not have a matching element 
                continue
            elif len(result) == 1: # if fragment in 1 node 
                result=result[0]
            else:
                result=sorted(result)[-1]  # if fragment is contained in multiple nodes              
  
            #add outcome of fragment to results
            results.append(Rep(*se, *result))
        
        connected = connect_and_reassign(results,tmp,tree)

        # remerge elements that got fragmented
        if alignment == True:
            connected = remerge(connected, tmp)

    lines = []
    if alignment == False:
        #make bed for TE content
        for _rep in connected:  # start, end, score, info
            strand,rep_name=_rep.info
            gft_id='RM_%s.%s' % (gft_id_n, rep_name)
            l=[chr, str(_rep.start), str(_rep.end), gft_id, str(_rep.score), strand] 
            lines.append('\t'.join(l) +'\n')
    else: 
        #make gff for alignment
        for _rep in connected: 
            
            strand,rep_name=_rep.info
            gft_id='RM_%s.%s' % (gft_id_n, rep_name)
            #output for repeatcraft
            ID = gft_id.split(".")
    
            if _rep in tmp:
                changed = False
            else:
                changed = True

            score = _rep.score

            info = "Tstart="+str(_rep.consensus_info[0])+";Tend="+str(_rep.consensus_info[1])+";Ttotal="+str(_rep.consensus_info[2])+";ID="+ID[1]+";annot_changed="+str(changed)

            l=[chr, "RepeatMasker", ID[2], str(_rep.start + 1), str(_rep.end), str(score), strand, ".",info, str(ID[3])]
            
            lines.append('\t'.join(l) +'\n')

    gft_id_n += 1

    return lines,gft_id_n 

def get_lines_no_resolving(tmp, chr, gft_id_n):
    lines = []
    if alignment == False:
        #make bed for TE content
        for _rep in tmp:  # start, end, score, info
            #construct line for bed
            strand,rep_name=_rep.info
            gft_id='RM_%s.%s' % (gft_id_n, rep_name)
            l=[chr, str(_rep.start), str(_rep.end), gft_id, str(_rep.score), strand] 
            lines.append('\t'.join(l) +'\n')
    else: 
        #make gff for alignment
        for _rep in tmp: 
            
            strand,rep_name=_rep.info
            gft_id='RM_%s.%s' % (gft_id_n, rep_name)
            #output for repeatcraft
            ID = gft_id.split(".")
    
            if _rep in tmp:
                changed = False
            else:
                changed = True

            #construct line for gff
            score = _rep.score
            info = "Tstart="+str(_rep.consensus_info[0])+";Tend="+str(_rep.consensus_info[1])+";Ttotal="+str(_rep.consensus_info[2])+";ID="+ID[1]+";annot_changed="+str(changed)
            l=[chr, "RepeatMasker", ID[2], str(_rep.start + 1), str(_rep.end), str(score), strand, ".",info, str(ID[3])]
            lines.append('\t'.join(l) +'\n')

    gft_id_n += 1

    return lines,gft_id_n 

def parse_line(ls):
    start= int(ls[5]) - 1
    end= int(ls[6])
    score= int(ls[0])
    strand=ls[8]
        
    if aligninput:
        if ls[8] == 'C' and len(ls) == 15:

            repname='%s.%s.%s' % (ls[9], ls[13], ls[14])
            strand = 'C'

            consensus_columns = [ls[10], ls[11], ls[12]]

        else:
            repname='%s.%s.%s' % (ls[8], ls[12], ls[13])
            strand = '+'

            consensus_columns = [ls[9], ls[10], ls[11]]
        
    else:
        repname='%s.%s.%s' % (ls[9], ls[10], ls[14])
        consensus_columns = [ls[11], ls[12], ls[13]]
    
    
    try:
        consensus_columns[0] = int(consensus_columns[0])
    except:
        consensus_columns[2] = int(consensus_columns[2])

    if type(consensus_columns[0]) == int:
        consensus_length = int(consensus_columns[1]) + int(consensus_columns[2].strip("()"))
        consensus_info = [consensus_columns[0], consensus_columns[1], consensus_length]
    else:
         consensus_length = int(consensus_columns[1]) + int(consensus_columns[0].strip("()"))
         consensus_info = [consensus_columns[2], consensus_columns[1], consensus_length]
    info=(strand, repname)
    
    return Rep(start, end, score, info, consensus_info)


def per_chr(reps, gft_id_n):

    if leave_overlap:
        reps=sorted(reps)
        for _rep in reps:
            
            tmp=[_rep]
            lines,gtf_id_n = get_lines_no_resolving(tmp, prev_chr, gft_id_n)
            if alignment:
                outgff.write(''.join(lines))
            else:
                outbed.write(''.join(lines))     
    else:
        print('resolving overlaps on %s...' % prev_chr)
        # prepare
        prev_end= (gap * -1) - 1 #?????????????
        reps=sorted(reps)

        maxend = 0
        
        # loop over 
        for _rep in reps:

            dist= _rep[0] - maxend + 1

            # check overlap with all previous elements, not only the last one
            if _rep[1] > maxend:
                maxend = _rep[1]

            # add sequence to overlap island if overlapping
            if dist <= gap:
                tmp.append(_rep)
            else:
            # if non overlapping, process previous island (can be single annotation)
                if prev_end > 0:#?????????????
                    # process island
                    lines,gft_id_n=collapse(tmp, prev_chr, gft_id_n)
                    if alignment:
                        outgff.write(''.join(lines))
                    else:
                        outbed.write(''.join(lines))     

                tmp=[_rep] # form new island
                
            prev_end=_rep[1] #?????????????

        lines,gft_id_n=collapse(tmp, prev_chr, gft_id_n) # collapse last island

    if alignment:
        outgff.write(''.join(lines))
    else:
        outbed.write(''.join(lines))  
    return gft_id_n


# main
if args.testrun is False:

    if alignment:
        outgff=open(ogff, 'wt')
    else:
        outbed=open(obed, 'wt')
        outbed.write('#date: %s\n' % _date)
        outbed.write('#version: %s %s\n' % (os.path.basename(__file__), version))
        outbed.write('#original file: %s\n' % args.i)
        outbed.write('#format: chr start end name bit_score strand\n')
    
    prev_chr='dummy'
    prev_id='dummy'
    gft_id_n=0
    Rep=collections.namedtuple('Rep', ['start', 'end', 'score', 'info', 'consensus_info'])


    if aligninput:
        clean_infile = "clean_" + os.path.basename(args.i)
        
        #clean align input file so that you only take the informative lines
        with open(args.i, "r") as infile, open(clean_infile, "w") as outfile:
            for line in infile:
                # check if line starts with a number (0-9)
                if re.match(r'^\d', line):
                    outfile.write(line)

        inputfile = clean_infile

    else:
        inputfile = args.i
    
    with open(inputfile) as infile:
        if not aligninput:
            for _ in range(3):
                next(infile)

        #Printing runmodes

        if leave_overlap:
            if tomerge:
                print(f'\nTEatree: not resolving overlaps, but only merging fragmented annotations for \'{args.i}\' for \'{args.mode}\' mode\n')
            else:
                if mode == "TEcontent":
                    print(f'\nTEatree: \'leave_overlap\' specified and \'mode\' = \'TEcontent\' ...\nNo merging for TEcontent mode\nNot resolving overlaps and not merging ... \nNothing to be done\nPlease specify other options\n')
                else:
                    print(f'\nTEatree: \'leave_overlap\' specified and \'mergemode\' = \'none\' ...\nNot resolving overlaps and not merging ... \nNothing to be done\nPlease specify other options\n')
                sys.exit(1)
        else:
            if tomerge:
                print(f'\nTEatree: resolving overlaps for \'{args.i}\' and merging fragmented annotations for \'{args.mode}\' mode\n')
            else:
                print(f'\nTEatree: resolving overlaps for \'{args.i}\' for \'{args.mode}\' mode\n')

        if tomerge == True:
            if mergemode == "both":
                print(f'TEatree will merge fragmented annotations with the same Repeatmasker ID and/or that share their consensus family and are close together relative position(< {args.gapsize} bp apart)\n')
            elif mergemode == "threshold":
                print(f'TEatree will merge fragmented annotations if they share their consensus family and are close together relative position (< {args.gapsize} bp apart)\n')
            elif mergemode == "ID":
                print(f'TEatree will merge fragmentedannotations with the same Repeatmasker ID\n')
        else:
            print("")

        if args.mode == "alignment":
            print(f'\'{args.mode}\' mode best for re-alignment of repeats to their consensus families, use \'TEcontent\' mode if interested in repeat content')
        if args.mode == "TEcontent":
            print(f'\'{args.mode}\' mode best for repeat content estimations, use \'alignment\' mode if interested in re-aligning repeats to their consensus family')
        

        for line in infile:
            ls=line.split()
            # some fa.out contains strange lines without ID --> take previous ID

            if len(ls) == 15 and ls[-1] == '*':
                ls=ls[:13]
                ls.append(prev_id)
                ls.append('*')
            elif len(ls) == 14:
                ls.append(prev_id)
            elif len(ls) <= 13:
                print("strange line found. Please check the format again:")
                print(line, end='')
                exit(1)
            if ls[10] in remove_simple_low_complex: # do not parse simple repeat or low complexity of specified
                continue
            _rep=parse_line(ls)


            # remove small fragments if specfied, only if alignment is not true
            # removal happens at the end if remove == True AND alignment == True
            if remove == True and _rep[1] - _rep[0] < threshold and alignment == False:
                continue

            # add everything of one chromosome to reps
            # once a new chromosome is reached, process the previous chromosome reps
            if not ls[4] == prev_chr:
                if not prev_chr == 'dummy':
                    gft_id_n=per_chr(reps, gft_id_n)
                prev_chr=ls[4]
                reps=[_rep]
            else:
                reps.append(_rep)
            prev_id=ls[-1]
   
    gft_id_n=per_chr(reps, gft_id_n) # process last chromosome 
    

    if alignment:
        outgff.close()
    else:
        outbed.close()
    log_file.close()


    # Check if errors created in error file
    if os.path.getsize("errors.log") != 0:
        print('\nerrors occurred: check the errors.log file\n')
    else:
        os.remove("errors.log")
    
    if not leave_overlap:
        print('Done resolving overlaps!')

    # Process output
    if alignment:
        if tomerge: # MERGE!
            lcnt = fuseTE.truefusete(gffp=ogff, gapsize=gapsize, outfile=olabel, mergemode=mergemode)
            mergeTE.extratruemergete(gffp=olabel, outfile=omerge, remove=remove, threshold=threshold, lcnt = lcnt)
            target_file = omerge
        else:
            target_file = ogff

        if checklength: # filter annotations based on consensus length covered
            print('\nFiltering annotations based on consensus length...\n')
            filter.filterlength(target_file, ofilter, lengthfile)

        # summary stats
        print('\nWriting summary files...\n')

        rcStatm.freqalign(inputfile, ofilter if checklength else target_file, frequencyfilealignclass, frequencyfilealignfamily, aligninput)

    else: # summary stats only + bp content
        print('\nWriting summary files...\n')
        rcStatm.freqcontent(args.i, obed, frequencyfilecontentclass, frequencyfilecontentfamily, aligninput)
        rcStatm.bpcontent(args.i, obed, contentfile, aligninput)

    print(f"TEatree done")


else: # for debug
        Rep=collections.namedtuple('Rep', ['start', 'end', 'score', 'info'])
        tree=IntervalTree()
        tree.insert(100, 200, 150, ('+', 'LTR'))
        tree.insert(110, 300, 20, ('-', 'SINE'))
        tree.insert(120, 290, 100, ('-', 'LINE'))
        tree.insert(280, 310, 60, ('+', 'DNA'))
        poss=[100, 110, 120, 200, 280, 290, 300, 310]
        results=[]
        for n in range(len(poss) - 1):
            se=poss[n:n + 2]
            result=tree.find_top_score(*se)
            if len(result) == 0:
                continue
            elif len(result) == 1:
                result=result[0]
            else:
                result=sorted(result)[-1]
            results.append(Rep(*se, *result))
        for rep in results:
            print(rep)
        print()
        final=connect(results)
        for rep in final:
            print(rep)
