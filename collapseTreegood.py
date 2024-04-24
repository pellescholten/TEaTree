#!/usr/bin/env python

'''
Author: Shohei Kojima @ RIKEN
Usage: python %prog -i in.fa.out -o out
Requirement: Python 3.6 or later
License: see LICENSE file
'''

import os,sys,gzip,datetime,collections,argparse,errno

version='2021/12/03'

'''
Version log:
version='2021/09/20': first release
version='2021/12/03': fixed a small bug for corner cases. Occasionally .fa.out contains strange lines without ID.
'''



description='''
    This script collapses TE annotations from RepeatMasker (i.e. genome.fa.out file).
    The TE annotations from RepeatMasker are frequently overlap each other. If so, when counting
    NGS reads mapping to TEs, reads mapping to two or more TE annotations may not be counted.
    In such case, it is preferable to use non-overlapping TE annotations. This script makes
    non-overlapping TE annotations by removing TEs with lower bit score. Processing of one file
    (e.g. hg38, mm10) requires ~3 min and ~100 MB RAM. This script only uses single thread.
    Requirement: Python 3.6 or later
'''

# args
parser=argparse.ArgumentParser(description=description)
parser.add_argument('-i', metavar='str', type=str, help='Specify input genome.fa.out file.', required=True)
parser.add_argument('-o', metavar='str', type=str, help='Specify output basename. [basename].gtf.gz and [basename].bed.gz will be generated.', required=True)
parser.add_argument('-gap', metavar='int', type=int, help='Optional. Specify gap distance to connect two TEs belonging to the same repeat. Default: 0', default=0)
parser.add_argument('-lvl', metavar='int', type=int, help='Optional. Specify percentage of overlap necessary to complete remove a repeat. Default: 80', default=80)
parser.add_argument('-threshold', metavar='int', type=int, help='Optional. Specify minimum length of basepairs a repeat must have after being cut. Default: 50', default=50)
parser.add_argument('-min', metavar='int', type=int, help='Optional. Minimal length of repeat to output. Default: 1', default=1)
parser.add_argument('-keep_simple_repeat', help='Optional. Specify if you want to keep "Simple_repeat" and "Low_complexity".', action='store_true')
parser.add_argument('-quiet', help='Optional. Specify if you do not want to output processing status.', action='store_true')
parser.add_argument('-v', '--version', action='version', version='Version: %s %s' % (os.path.basename(__file__), version))
parser.add_argument('-testrun', action='store_true', help=argparse.SUPPRESS)
args=parser.parse_args()

# set up
ogtf='%s.gtf.gz' % args.o
obed='%s.bed.gz' % args.o
gap=args.gap
level = args.lvl / 100
threshold = args.threshold
annots=('gene', 'transcript', 'exon')
remove_simple_low_complex={'Simple_repeat', 'Low_complexity'}
if args.keep_simple_repeat is True:
    remove_simple_low_complex={}
_date=datetime.datetime.now()

# check input and output
if os.path.exists(args.i) is False:
    print(FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), args.i))
    exit(1)
outdir=os.path.dirname(args.o)
if '/' in outdir:
    os.makedirs(outdir, exist_ok=True)



class IntervalTree:
    def __init__(self):
        self.root=None
        
    def insert(self, start, end, score, info):
        if self.root is None:
            self.root=IntervalNode(start, end, score, info)
        else:
            self.root=self.root.insert(start, end, score, info)
    
    def find_top_score(self, start, end):
        return self.root.find_top_score(start, end)
    
    def reassign_short_fragments(self, fragment, connected):
        return self.root.reassign_short_fragments(fragment, connected)


class IntervalNode:
    def __init__(self, start, end, score, info):
        self.start=start
        self.end=end
        self.score=score
        self.info=info
        self.maxend=self.end
        self.minstart=self.start
        self.left=None
        self.right=None

    def replace_root(self,start,end,score,info):
        # replace root and corresponding branches with new node
        # while deleting old root from tree
        root = IntervalNode(start, end, score, info)
        if self.left:
            root.left = self.left
        if self.right:
            root.right = self.right
        return root

    def insert(self, start, end, score, info):
        root=self

        # Check whether there are overlaps of more than 80 (or other value determined in command line) percent shared sequence
        # If sequence is shared over the threshold than remove the sequence with lower score 
        # If not, keep sequence in tree to be resolved
        if start > self.start and start < self.end and end > self.end:  # start>self.start unnecessary as file is ordered, if unordered file: would require code for start<self.start
            overlap = self.end - start

            # if root (almost) contained in new element and new element is better--> delete and replace root
            # root:              ------------     SW = 500       -->  
            # new element:         ------------      SW = 1000   -->      ------------      SW = 1000   
            if overlap > level * (self.end - self.start) and self.score < score:
                return self.replace_root(start,end,score,info)
            
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
                    return self.replace_root(start,end,score,info)
                else:
                    return root

        # Sane as previous section but for chimers
        if start >= self.start and start <= self.end and end <= self.end:
            # If nested is worse or equal --> do not add to tree
            # nesting:      --------------     SW = 1000       -->   only NESTING remains     --------------      SW = 1000 
            # nested:        ------------      SW = 500        -->   
            if score <= self.score:
                return root
            
            # If nested is better and the nesting is for more than 80 (or other) percent contained in nested --> delete nesting
            # nesting:      --------------     SW = 500        -->
            # nested:        ------------      SW = 1000       -->   only NESTED remains     ------------      SW = 1000 
            elif (end - start)/(self.end - self.start) > level:
                root = IntervalNode(start, end, score, info)
                if self.left:
                    root.left = self.left
                if self.right:
                    root.right = self.right
                return root
            # If none of these, keep chimer in tree, to be resolved later...
        

        if start >= self.start:
            # insert to right tree
            if self.right: #if there is already a node the right of the root, repeat with that node
                self.right=self.right.insert(start, end, score, info)
            else: # add node to the right
                self.right=IntervalNode(start, end, score, info)
            # build heap

            # change root if score is higher
            if self.score < self.right.score:
                root=self.rotateleft()
        else: # not in ordered dataset
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
            root.minstart=min(root.end, root.right.minstart, root.left.minstart)
        elif root.right:
            root.maxend=max(root.end, root.right.maxend)
            root.minstart=min(root.end, root.right.minstart)
        elif root.left:
            root.maxend=max(root.end, root.left.maxend)
            root.minstart=min(root.end, root.left.minstart)
        return root
    
    def rotateright(self):
        root=self.left
        self.left=self.left.right
        root.right=self
        if self.right and self.left:
            self.maxend=max(self.end, self.right.maxend, self.left.maxend)
            self.minstart=min(self.end, self.right.minstart, self.left.minstart)
        elif self.right:
            self.maxend=min(self.end, self.right.maxend)
            self.minstart=min(self.end, self.right.minstart)
        elif self.left:
            self.maxend=max(self.end, self.left.maxend)
            self.minstart=min(self.end, self.left.minstart)
        return root

    def rotateleft(self):

        # remember: score of right is higher than score of root
        # new node becomes root
        root=self.right

        # node on left becomes node on right
        # old root becomes node on left
        self.right=self.right.left
        root.left=self

        if self.right and self.left:
            self.maxend=max(self.end, self.right.maxend, self.left.maxend)
            self.minstart=min(self.end, self.right.minstart, self.left.minstart)
        elif self.right:
            self.maxend=max(self.end, self.right.maxend)
            self.minstart=min(self.end, self.right.minstart)
        elif self.left:
            self.maxend=max(self.end, self.left.maxend)
            self.minstart=min(self.end, self.left.minstart)
        return root
    
    def find_top_score(self, start, end):
        # if segment is contained in current node
        if start < self.end and end > self.start:
            return [(self.score, self.info)]
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
        if self.info == fragment.info or self.score > fragment.score:
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
                if fragment.end < self.end and fragment.start > self.start: # chimer
                    return [(fragment.start, fragment.end, self.score, self.info)]
                else: # non chimer
                    return [(fragment.start, self.end, self.score, self.info)]
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
    for _rep in l:
        if _rep.end - _rep.start == 0:
            continue
        dist= _rep.start - prev_e
        if _rep.info == prev_i and dist <= gap:
            prev_e=_rep.end
            prev_r = max(prev_r, _rep.score)
        else:
            if not prev_i is None:
                connected.append(Rep(prev_s, prev_e, prev_r, prev_i))
            prev_s,prev_e,prev_r,prev_i=_rep
    connected.append(Rep(prev_s, prev_e, prev_r, prev_i))
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
        if fragment not in tmp and fragment.end - fragment.start < threshold:
            # try to find new element for fragment
            result = tree.reassign_short_fragments(fragment, connected)

            connected.remove(fragment)
            
            if len(result) == 0: #if fragment does not have a matching element --> do not add back
                continue
            elif len(result) == 1: # if fragment in 1 node 
                result=result[0]    
            else:
                result=sorted(result)[-1]  # if fragment is contained in multiple nodes

            # if succesfull, add fragment with new info to results & repeat
            breakpoint()
            connected.append(Rep(*result))
            connected.sort(key = lambda element: element[0])
            connected = connect_and_reassign(connected,tmp,tree)
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
            se=poss[n:n + 2]

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
        

        #if results[0].start == 12248790:
        #    breakpoint()
        connected = connect_and_reassign(results,tmp, tree)

    gtf_lines=[]
    bed_lines=[]
    for _rep in connected:  # start, end, score, info
        if _rep.end - _rep.start < args.min:
            continue
        strand,rep_name=_rep.info
        gft_id='RM_%s.%s' % (gft_id_n, rep_name)
        gene_attr='gene_id "%s"; gene_name "%s"; bit_score "%s";' % (gft_id, gft_id, _rep.score)
        tran_attr='gene_id "%s"; transcript_id "t_%s"; gene_name "%s"; bit_score "%s";' % (gft_id, gft_id, gft_id, _rep.score)
        exon_attr='gene_id "%s"; transcript_id "t_%s"; gene_name "%s"; exon_number 1; exon_id "e_%s"; bit_score "%s";' % (gft_id, gft_id, gft_id, gft_id, _rep.score)
        attrs=(gene_attr, tran_attr, exon_attr)
        for annot,attr in zip(annots, attrs):
            l=[chr, 'RepeatMasker', annot, str(_rep.start + 1), str(_rep.end), '.', strand, '.', attr]
            gtf_lines.append('\t'.join(l) +'\n')
        l=[chr, str(_rep.start), str(_rep.end), gft_id, str(_rep.score), strand]
        bed_lines.append('\t'.join(l) +'\n')
        gft_id_n += 1
    return gtf_lines, bed_lines, gft_id_n


def parse_line(ls):
    start= int(ls[5]) - 1
    end= int(ls[6])
    score= int(ls[0])
    strand=ls[8]
    if strand == 'C':
        strand='-'
    repname='%s.%s.%s' % (ls[9], ls[10], ls[14])
    info=(strand, repname)
    return Rep(start, end, score, info)


def per_chr(reps, gft_id_n):
    # prepare
    if args.quiet is False:
        print('processing %s...' % prev_chr)
    prev_end= (gap * -1) - 1
    reps=sorted(reps)

    maxend = 0
    
    # loop over 
    for _rep in reps:

        dist= _rep[0] - maxend + 1

        # make sure to check overlap with all previous elements, not only the last one (maxend)
        if _rep[1] > maxend:
            maxend = _rep[1]

        # add sequence to overlap island if ... overlapping
        if dist <= gap:
            tmp.append(_rep)
        else:
        # if non overlapping, process previous island
            if prev_end > 0:

                # process island
                gtf_lines,bed_lines,gft_id_n=collapse(tmp, prev_chr, gft_id_n)
                outgtf.write(''.join(gtf_lines))
                outbed.write(''.join(bed_lines))
            tmp=[_rep]
             
        prev_end=_rep[1]

    gtf_lines,bed_lines,gft_id_n=collapse(tmp, prev_chr, gft_id_n)
    if len(bed_lines) >= 1:
        outgtf.write(''.join(gtf_lines))
        outbed.write(''.join(bed_lines))
    return gft_id_n



# main
if args.testrun is False:
    outgtf=gzip.open(ogtf, 'wt')
    outbed=gzip.open(obed, 'wt')
    
    outgtf.write('##format: gtf\n')
    outgtf.write('##date: %s\n' % _date)
    outgtf.write('##version: %s %s\n' % (os.path.basename(__file__), version))
    outgtf.write('##original file: %s\n' % args.i)
    outbed.write('#date: %s\n' % _date)
    outbed.write('#version: %s %s\n' % (os.path.basename(__file__), version))
    outbed.write('#original file: %s\n' % args.i)
    outbed.write('#format: chr start end name bit_score strand\n')
    
    prev_chr='dummy'
    prev_id='dummy'
    gft_id_n=0
    Rep=collections.namedtuple('Rep', ['start', 'end', 'score', 'info'])
    with open(args.i) as infile:
        for _ in range(3):
            next(infile)
        for line in infile:
            ls=line.split()
            # some fa.out contains strange lines without ID
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
            if ls[10] in remove_simple_low_complex:
                continue
            _rep=parse_line(ls)
            # per chr

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
    gft_id_n=per_chr(reps, gft_id_n)
    
    outgtf.close()
    outbed.close()
    print('\n%s repeats were retained.\nThank you for using this script!\n' % gft_id_n)

else:
    # for debug
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
