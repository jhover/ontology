#!/usr/bin/env python
#
#  Various facilities for handling ontology info. 
#  
#  1) Consume .obo file and build walkable tree of objects.
#  output info in dataframe or dict form.     
#
#  2) Consume .gaf or (.gpad + .gpi) and (optionally alternate .obo) files. 
#  create matrix
#  extend with is_a relationships so all parent functions are flagged. 
#  output dataframe
#  save/load CSV
#          
#            F1   F2   f3   f4  ...
#   gene 1   0     1    0    1
#   gene 2   1     0    0    0
#   gene 3   0     1    1    0 
#    ...  
#
#   Use GAF or OBO files from    
#    https://current.geneontology.org/ontology/go.obo
#    
#   And species-specific annotations from https://current.geneontology.org/products/pages/downloads.html
#   E.g. https://current.geneontology.org/annotations/mgi.gaf.gz 
#     
#


import argparse
import collections
from configparser import ConfigParser
import itertools
import logging
import os
import pickle
import pprint
import sys
import traceback

import pandas as pd
import numpy as np

from scipy import sparse
#import pyarrow as pa
#import pyarrow.parquet as pq

GOASPECTMAP= { 'biological_process' : 'bp',
               'cellular_component' : 'cc',
               'molecular_function' : 'mf',
               'external'           : 'ex'
            }

UPASPECTMAP = { 'C': 'cc',
                'F': 'mf',
                'P': 'bp'
              }


class Ontology(dict):
    """
    gomatrix:  goterm x goterm np.ndarray fully propagated relationships. 
    gotermidx: { <str> goterm : <int> index in gomatrix, ... } 
    gotermlist:
    
    NOTE: dict key order is now assumed stable as of Python 3.7. If this is run on 
    earlier version, unstable key order will cause chaos. 
    
    """
    
    instance = None
    
    def __init__(self, gomatrix, gotermidx, altidx):
        self.kname = self.__class__.__name__
        self.lkname = self.kname.lower()
        self.log = logging.getLogger(self.kname)
        
        # a list of np.arrays
        self.data = gomatrix
        
        # a dictionary of goterms -> array_index
        self.gotermidx = gotermidx
        
        # a dictionary of (alternate) goterms -> (real) goterms
        self.altidx = altidx
        # list of go terms
        # **depends on key sort order stability**
        self.gotermlist = list(gotermidx)
        
        
    def __getitem__(self, key):
        #
        val = None
        try:
            val = self.data[self.gotermidx[key]]
        except KeyError:
            #self.log.debug(f"bad primary goterm key {key}, looking up alt. ")
            try:
                realgt = self.altidx[key]
                realidx = self.gotermidx[realgt]
                val = self.data[realidx]
            except KeyError:
                if not math.isnan(key):
                    self.log.error(f"bad altidx key: '{key}' type:{type(key)} should not happen.")
                # otherwise a nan is sometimes expected for unannotated entries. 
        return val

    def __repr__(self):
        #return repr(self.__dict__)
        s = "{"
        for k in self.gotermidx.keys(): 
            s += f"'{k}' : {self.data[self.gotermidx[k]]} "
        s += "}"
        return s 

    def __len__(self):
        return len(self.data)

    def keys(self):
        return self.gotermidx.keys()
    
    def indexof(self, key):
        return self.gotermidx[key]

    def keyof(self, indexid ):
        return self.gotermlist[indexid]











def get_default_config():
    dc = os.path.expanduser('~/git/ontology/etc/ontology.conf')
    cp = ConfigParser()
    cp.read(dc)
    return cp
 
def format_config(cp):
    cdict = {section: dict(cp[section]) for section in cp.sections()}
    s = pprint.pformat(cdict, indent=4)
    return s


def test(config):
    go = GeneOntology(config)
    df = go.get_df()
    print(str(df))
    df = go.get_df()
    
    goidx = go.get_tree_termindex()
    obj = goidx['GO:0009987']
    print(obj)
    
    testterms = os.path.expanduser('~/play/cafa4/gocodes.txt')
    #testterms = os.path.expanduser('~/play/cafa4/smallgocodes.txt')

    golists = {}
    fh = open(testterms, 'r')
    for line in fh:
        gt = line.strip()
        print("looking up term: %s" % gt)
        gtobj = goidx[gt]
        print("goterm: %s -> is_a: %s" % (gt,  gtobj.superclassesstr()) )

        
def get_ontology_object(config, usecache=True):
    if Ontology.instance is None:
        build_ontology(config, usecache)
    return Ontology.instance
    

def load_ontology():
    conffile = os.path.expanduser("~/git/cshlwork/etc/ontology.conf")
    config = ConfigParser() 
    config.read(conffile)
    go = GeneOntology(config)
    #go.get_tree_termindex()
    return go


def make_hotmatrix(obofile, gaffile, outfile, config=None, usecache=False):
    if config is None:
        config = get_default_config()
    build_ontology(config, usecache=usecache)
    df = build_gomatrix(config, usecache, gaffile=gaffile)
    if outfile is not None:    
        if outfile.endswith(".csv"):
            logging.debug(f"Saving CSV matrix to {outfile}")       
            df.to_csv(outfile, index=True, header=True, sep=',')
       
        elif outfile.endswith(".tsv"):
            logging.debug(f"Saving TSV matrix to {outfile}")     
            df.to_csv(outfile, index=True, header=True, sep='\t')
    else:
        logging.debug("no output required. Just returning df.")    
    return df


def build_gomatrix(config, usecache=False, gaffile=None):
    logging.debug(f"Getting GO ontology object... ")
    ontobj = get_ontology_object(config, usecache=True)
    logging.debug("Parsing GAF file...")
    lol = parse_gaf(gaffile, config)
    logging.debug(f"List of lists: {len(lol)}")
    
    logging.debug(f"Building matrix...")
    genebygo, genelist = build_genematrix(lol, ontobj)
    
    logging.debug(f"genebygo sum={genebygo.sum()} should be > num genes.   ")

    logging.debug("converting to sparse matrix.")    
    genebygo = sparse.lil_matrix(genebygo, dtype=bool)
    gomatrix = sparse.lil_matrix(ontobj.data, dtype=bool)
    logging.info("conversion complete. Performing matrix multiplication.")   
    
    genebygo =  genebygo @ gomatrix  
    logging.debug(f"genebygo sum={genebygo.sum()} after matrix multiplication ")    

    logging.debug("Converting back to dense matrix....")
    genebygo = genebygo.todense()    
    logging.info(f"Done. genebygo: t{type(genebygo)} shape {genebygo.shape} dtype {genebygo.dtype} ")
    genebygo = genebygo.astype('uint8')
    logging.info(f"Done. genebygo: dtype {genebygo.dtype} ")               
    df = pd.DataFrame(genebygo, dtype='uint8', index=genelist, columns=ontobj.gotermlist)
    logging.debug("Done.")       
    return df


def build_genematrix(goadata, ontobj):
    '''
    Takes goadata [ [ <gene>, <goterm> ], [ <gene>, <goterm> ] ...]
    produces boolean matrix of genes by goterm vector. 

    '''
    logging.debug(f"In build_genematrix for genelist of {len(goadata)} genes...")
    gotermlist = ontobj.gotermlist  # columnlabels
    genelist = []
    idxdict = {}
    golen = len(gotermlist)
    genedict = {}
    row = 0
    for e in goadata:
        (gene, goterm ) = e
        try:
            gd = genedict[gene]
            gd.append(goterm)
        except KeyError:
            genedict[gene] = [goterm]
            genelist.append(gene)
            idxdict[gene] = row
            row += 1
            
    genelen= len(genedict)  
    logging.info(f"{genelen} genes with {len(goadata)} associated goterms.") 
    logging.debug(f"building {genelen} by {golen} boolean array...")
    bool_array = np.full( (genelen, golen), False, dtype=bool)
    logging.debug(f"Filling in gene x goterm boolean matrix...")
    for gene in genedict.keys():
        for goterm in genedict[gene]:
            bool_array[ idxdict[gene] , ontobj.indexof(goterm)] = True                      
    logging.debug(f"Done. Handled {row} genes. ")
    #return govectors
    return bool_array, genelist



def build_ontology(config, usecache):
    """
    obofile=~/data/go/go.obo
    cachedir = ~/play/cafa4      
    
    from parse_obo():
    { 'GO:2001315': 
         {'is_a': ['GO:0009226', 'GO:0046349', 'GO:2001313'], 
         'part_of': [], 
         'goterm': 'GO:2001315', 
         'goname': 'UDP-4-deoxy-4-formamido-beta-L-arabinopyranose biosynthetic process', 
         'goasp': 'bp', 
       ...  
    }
 
    result:  Numpy boolean matrix of all relations in ontology, ordered by sorted goterm name.  
       
    """
    #global GOMATRIX
    #global GOTERMIDX
    #global ALTIDDICT
    #global GOTERMLIST
    # def __init__(self, gomatrix, gotermidx, altidx):
    
    logging.debug(f"usecache={usecache}")
    cachedir = os.path.expanduser(config.get('ontology','cachedir'))
    ontologyfile = f"{cachedir}/ontology.npy"
    termidxfile = f"{cachedir}/gotermindex.pickle"
    altiddictfile = f"{cachedir}/altiddict.pickle"
    include_partof = config.getboolean('ontology','include_partof')
    
    gomatrix = None
    
    if os.path.exists(ontologyfile) and usecache:
        logging.debug("Cache hit. Using existing matrix...")
        gomatrix = np.load(ontologyfile)
        logging.debug(f"loaded matrix: {matrix_info(gomatrix)} from {ontologyfile}")
        
        f = open(termidxfile, 'rb')
        gotermidx = pickle.load(f)
        f.close()

        f = open(altiddictfile, 'rb')
        altiddict = pickle.load(f)
        f.close()
                
        logging.debug(f"goterm index, e.g. : \n{list(gotermidx)[0]} :  {gotermidx[list(gotermidx)[0]]} ")
    
    else:
        (godict, altiddict) = parse_obo(config)
        
        # get keys from dict
        gotermlist = list(godict)
        logging.debug(f"parsed obo with {len(gotermlist)} entries. ")
        logging.debug(f"example entry:\n{gotermlist[0]}")
        logging.debug("sorting goterms")
        gotermlist.sort()
        logging.debug(f"sorted: e.g. {gotermlist[0:5]} ")
        logging.debug("creating goterm index dict.")
        #
        i = 0
        gotermidx = {}
        for gt in gotermlist:
            gotermidx[gt] = i
            i = i + 1
              
        logging.debug(f"creating zero matrix of dimension {len(gotermlist)}")
        shape = (len(gotermlist), len(gotermlist))
        gomatrix = np.zeros( shape, dtype=bool )
        logging.debug(f"filling in parent matrix for all goterms...")
        for gt in godict.keys():
            for parent in godict[gt]['is_a']:
                    gomatrix[gotermidx[gt]][gotermidx[parent]] = True
        if include_partof:
            logging.debug("Including part_of relationships as is_a")
            for gt in godict.keys():
                for parent in godict[gt]['part_of']:
                        gomatrix[gotermidx[gt]][gotermidx[parent]] = True
        
        #logging.debug(f"initial matrix:\n{print_square(gomatrix, GOTERMLIST)}")
        logging.debug("Calculating sparsity...")
        logging.debug(f"sparsity = { 1.0 - np.count_nonzero(gomatrix) / gomatrix.size }")
        logging.debug("converting to sparse matrix.")
        gomatrix = sparse.lil_matrix(gomatrix, dtype=bool)
        logging.debug(f"converging matrix: {matrix_info(gomatrix)}")
        gomatrix = converge_sparse(gomatrix)
        logging.info(f"got converged matrix:\n{matrix_info(gomatrix)} ")
        logging.debug(f"converged matrix sum={gomatrix.sum()}")
        #logging.debug("Calculating sparsity...")
        #sparsity = 1.0 - np.count_nonzero(gomatrix) / gomatrix.size
        #logging.debug(f"sparsity = { 1.0 - np.count_nonzero(gomatrix) / gomatrix.size }")        
        gomatrix = gomatrix.todense()
        gomatrix = np.asarray(gomatrix, dtype='bool')
            
        logging.debug(f"Caching all values/indexes...")
        logging.debug(f"Saving matrix: {matrix_info(gomatrix)} to {ontologyfile}")
        np.save(ontologyfile, gomatrix)

        logging.debug(f"Saving gotermidx {len(gotermidx)} items to {termidxfile}")
        f = open(termidxfile, 'wb')   
        pickle.dump(gotermidx, f )
        f.close()
        
        logging.debug(f"Saving altiddict {len(altiddict)} items to {altiddictfile}.")
        f = open(altiddictfile, 'wb')
        pickle.dump(altiddict, f)
        f.close()
        
        logging.debug("Done constructing input for Ontology().")

    ontobj = Ontology(gomatrix, gotermidx, altiddict)
    # set global instance
    Ontology.instance = ontobj  
    logging.debug("Done creating Ontology object.")





def get_altiddict(config, usecache):
    if ALTIDDICT is not None:
        return ALTIDDICT
    else:
        build_ontology(config, usecache)    
        return ALTIDDICT


def parse_obo(config, obofile=None):
    """
    creates dict of dicts. key is goterm, contents is dict of 
       
       goterm  ""
       goname  ""
       goasp   ""
       godef   ""
       goisa   [] of goterms
       gohasa  [] of goterms
    
    
    """
    if obofile is None:
        obofile = os.path.expanduser(config.get('ontology','obofile'))
    filehandle = open(obofile)
    godict = {}
    altids = {}
    current = None
    logging.info(f"Parsing file {obofile}")
    try:
        for line in filehandle:
            if line.startswith("[Typedef]"):
                godict[current['goterm']]= current
                break
            elif line.startswith("[Term]"):     
                if current is not None:
                    godict[current['goterm']]= current
                # create new item...
                current = {}
                current['is_a'] = []
                current['part_of'] = []
                
            elif line.startswith("id: "):
                current['goterm'] = line[4:].strip()
                
            elif line.startswith("name: "):
                current['goname'] = line[6:].strip()
            
            elif line.startswith("namespace: "):
                asp = line[11:].strip()
                
                current['goasp'] = GOASPECTMAP[asp]

            # must create a separate mapping from alt_ids that maps to
            # primary, so it can be added to gomatrix properly
            elif line.startswith("alt_id: "):
                #current['alt_id'] = line[8:18].strip()
                #current['goasp'] = GOASPECTMAP[asp]            
                altid = line[8:18].strip()
                altids[altid] = current['goterm'] 
            
            #elif line.startswith("def: "):
            #    current['godef'] = line[5:].strip()

            #elif line.startswith("synonym: "):
            #    current.synonym.append(line[9:].strip())

            elif line.startswith("is_a: "):
                current['is_a'].append(line[6:16].strip())
            
            elif line.startswith("relationship"):
                if "part_of" in line:
                    current['part_of'].append(line[22:32])
                                 
    except Exception as e:
        traceback.print_exc(file=sys.stdout)                
    
    logging.info(f"Parsed file with {len(godict)} terms and {len(altids)} alt terms.")    
    return (godict, altids)
    

def parse_gaf(infile=None, config=None):
    '''
    create list of lists of GOA database. 
    [ <gene>, <goterm> ]

    '''
    if config is None:
        config = get_default_config()
    if infile is None:
        filepath = os.path.expanduser(config.get('ontology','gaffile'))
    else:
        filepath = infile
    
    allentries = []    
    try:
        logging.debug(f" attempting to open '{filepath}'")
        with open(filepath, 'r') as f:
            current = None
            sumreport = 1
            suminterval = 10000
            repthresh = sumreport * suminterval
            while True:
                line = f.readline()
                if line == '':
                    break
                if line.startswith("!"):
                    pass
                else:
                    fields = line.split('\t')
                    fields = fields[:7]
                    geneterm = [ fields[2], fields[4] ]
                    allentries.append(geneterm)

    except FileNotFoundError:
        logging.error(f"No such file {filepath}")   
        
    except Exception as e:
        traceback.print_exc(file=sys.stdout)                
    logging.info(f"Parsed file with {len(allentries)} entries" )
    logging.debug(f"Some entries:  {allentries[10:15]}")
    return allentries


def converge_sparse(matrix):
    logging.debug(f"starting matrix: \n{matrix_info(matrix)}")
    #logging.debug(f"{print_square(matrix.todense(), GOTERMLIST)}")
    oldval = 0
    logging.debug("Summing inbound matrix...")
    newval = matrix.sum()
    logging.debug("Beginning convergence loop.")
    icount = 0
    while oldval != newval:
        #logging.debug(f"Inbound matrix:\n{matrix_info(matrix)}")
        #logging.debug(f"oldval={oldval} newval={newval}")
        oldval = newval
        if not isinstance(matrix,  sparse.lil.lil_matrix): 
            logging.debug(f"{type(matrix)} is not scipy.sparse.lil.lil_matrix, converting... ")
            #matrix = sparse.lil_matrix(matrix, dtype=np.bool)
            matrix = sparse.lil_matrix(matrix, dtype=bool)
        else:
            pass
            #logging.debug("matrix already lil_matrix...")
        #logging.debug("Multiplying...")
        mult = matrix @ matrix
        #logging.debug("Adding back original...")
        matrix = mult + matrix
        #logging.debug("Getting new sum...")
        newval = matrix.sum()
        #logging.debug(f"New matrix {icount}:\n{matrix.todense()}")
        #logging.debug(f"{print_square(matrix.todense(), GOTERMLIST)}")
        icount += 1
    logging.debug(f"Convergence complete. {icount} iterations. matrix:\n{matrix_info(matrix)}")
    return matrix


def matrix_info(matrix):
    """
    Returns string of inexpensive info on large matrix for logging. 
    """
    return f"type: {type(matrix)} shape: {matrix.shape} dtype: {matrix.dtype}"


###########################################

class GOMatrix(object):
    """
    produces matrix of genes/proteins x goterms
    optionally extended by all is_a relationships..
    
    """
    def __init__(self, config, gaffile):
        #self.kname = self.__class__.__name__
        self.log = logging.getLogger(self.__class__.__name__)
        self.config = config
        self.cachedir = os.path.expanduser(config.get('ontology','cachedir'))
        self.cachefile = 'gomatrix.csv'
        self.cachepath = "%s/%s" % (self.cachedir, self.cachefile)
        self.df = None
        self.log.debug("GOMatrix initialized.")


    def get_df(self, cache=True):
        self.df = self._df_from_cache()
        if self.df is None:
            self.execute()
        return self.df
        
    def _df_from_cache(self):
        self.log.debug("Trying to load from cache: %s" % self.cachepath )  
        ret = None
        try:
            self.df = pd.read_csv(self.cachepath, index_col=0)
        except FileNotFoundError:
            logging.info('no cached file found')
        return ret


class GOTerm(object):
    #  Cache of showing all (non-redundant) is_a parents for given GO Term. 
    ISA_LISTCACHE = {}
    
    def __init__(self ):
        self.log = logging.getLogger(self.__class__.__name__)
        self.goterm = None
        self.goname = None
        self.goaspect = None
        self.godef = None
        self.synonym = []
        self.is_a = []
        self.part_of = []
    
    def superclasses(self):
        '''
        Return non-redundant list of all goterm strings that are is_a|part_of parents of this term
        '''
        # Root terms
        nl = None
        
        if self.goterm in GOTerm.ISA_LISTCACHE.keys() :
                self.log.debug("%s Cache hit, using..." % self.goterm)
                nl = GOTerm.ISA_LISTCACHE[self.goterm]
        
        elif len(self.is_a) == 0:
            self.log.debug("%s Cache miss" % self.goterm)
            self.log.debug("%s ROOT TERM!" % self.goterm)
            nl = []
            self.log.debug("%s Storing cache entry= %s" % (self.goterm, nl))
            GOTerm.ISA_LISTCACHE[self.goterm] = nl 

        # Only has one parent, simply construct list
        elif len(self.is_a) == 1:
            # no entry for this isa_list. Must construct...
            self.log.debug("%s Cache miss" % self.goterm)
            item = self.is_a[0]
            itemlist = [item.goterm] + item.superclasses()
            self.log.debug("%s Got itemlist %s for %s" % (self.goterm, itemlist, item.goterm))
            nl =  itemlist
            self.log.debug("%s new list is %s" % (self.goterm, nl))
            # add item term to item's isalist. 

            self.log.debug("%s Storing cache entry= %s" % (self.goterm, nl))
            GOTerm.ISA_LISTCACHE[self.goterm] = nl
        
        # Has multiple paths to root, must use sets to de-duplicate. 
        else:
            self.log.debug("%s Cache miss" % self.goterm)   
            gl = []              
            for item in self.is_a:
                gl.extend( [item.goterm] + item.superclasses())
            
            nl = list(collections.OrderedDict.fromkeys(gl))
            # store de-duped list to cache, preserving order
            self.log.debug("%s Storing cache entry= %s" % (self.goterm, nl))
            GOTerm.ISA_LISTCACHE[self.goterm] = nl

        self.log.debug("%s returning %s" % (self.goterm, nl))
        return nl    
    
    def superclassesstr(self):
        '''
        Return string representation of all is_a|part_of parents of this term. 
        
        '''
        il = self.superclasses()
        s = ' '.join(il)
        return s   
    
    def __repr__(self):
        s = "GOTerm:"
        for atr in ['goterm','goname', 'goaspect']:
            s += " %s=%s" % (atr, self.__getattribute__(atr))
        for isa in self.is_a:
            s += f" isa={isa.goterm}"
        #s += f" {self.is_a}"
        return s
        

class GeneOntology(object):
              
    NSMAP= { 'biological_process' : 'bp',
             'cellular_component' : 'cc',
             'molecular_function' : 'mf',
             'external'           : 'ex'
            }
        
    def __init__(self, config, obofile=None):
        self.config = config
        self.kname = self.__class__.__name__
        self.lkname = self.kname.lower()
        self.log = logging.getLogger(self.kname)
        self.log.debug(f'config={config} obofile={obofile}')
        self.cachedir = os.path.expanduser(self.config.get('global' ,'cachedir'))
        self.cachefile = "%s.csv" % self.lkname
        self.cachepath = "%s/%s" % (self.cachedir, self.cachefile)
        self.obofile = os.path.expanduser(obofile)
        if obofile is None:
            self.obofile = os.path.expanduser(self.config.get('ontology','obofile'))
        
        self.goidx = {}   #  { 'GO:XXXXXX' : GoTermObject, }      
        self.df = None
        self.get_tree_termindex()
        self.log.debug(f"GeneOntology with {len(self.goidx)} terms initialized.")
       
    
    def get_tree_termindex(self):
        '''
        Builds full object tree, 
        provides dictionary with GOTerms as keys, object in tree as value. 
        object is_a/part of relations can be used to traverse tree.
        
        '''
        filehandle = None
        try:
            self.log.debug("opening file %s" % self.obofile)
            filehandle = open(self.obofile, 'r')
            self._parse2tree(filehandle)
            filehandle.close()
                
        except FileNotFoundError:
            self.log.error("No such file %s" % self.obofile)        
        
        finally:
            if filehandle is not None:
                filehandle.close()
        self._add_references()        
        return self.goidx


    def get_term(self, goterm):
        self.log.debug(f"Looking up term {goterm}")
        return self.goidx[goterm]

    def get_parent_list(self, goterm):
        self.log.debug(f"Parents for {goterm}")
        gtobj = self.goidx[goterm]
        pl = gtobj.superclasses()
        # print("goterm: %s -> is_a: %s" % (gt,  gtobj.get_isastr()) )
        return pl

    def _parse2tree(self, filehandle):
        '''
        Create in-memory GO tree with GOTerm index. 
        Two passes:
            1. Do not resolve foreign references, just add them as strings...
            2. Resolve foreign refs, replacing strings w/ object refs. 
                
        goidx = { 'GO:0000001' : <GOTerm Object>,
                  'GO:0000002' : <<GOTerm Object>,
        }
        '''      
        self.goidx = {}
        current = None
        self.log.info("Parsing file %s" % self.obofile)
        try:
            for line in filehandle:
                #self.log.debug("line is %s" % line)
                if line.startswith("[Term]"):     
                    #self.log.debug("found term")
                    if current is not None:
                        self.goidx[current.goterm]= current
                        current = GOTerm()
                    else:
                        #self.log.debug("Creating new GOTerm object...")
                        current = GOTerm()
                    
                elif line.startswith("id: "):
                    current.goterm = line[4:].strip()
                    
                elif line.startswith("name: "):
                    current.goname = line[6:].strip()
                
                elif line.startswith("namespace: "):
                    asp = line[11:].strip()
                    current.goaspect = GeneOntology.NSMAP[asp]
                
                elif line.startswith("def: "):
                    current.godef = line[5:].strip()

                elif line.startswith("synonym: "):
                    current.synonym.append(line[9:].strip())

                elif line.startswith("is_a: "):
                    current.is_a.append(line[6:16].strip())
                    
        except Exception as e:
            traceback.print_exc(file=sys.stdout)                
        
        self.log.info("Parsed file with %d terms" % len(self.goidx) )

        
    def _add_references(self):
        '''
        Go through goidx and add foreign refs...
        '''
        for gt in self.goidx.keys():
            self.log.debug(f'handling term={gt}')
            gto = self.goidx[gt]
            isalist = gto.is_a
            newisa = []
            for igt in isalist:
                self.log.debug(f'handling isalist term={igt}')
                try:
                    newisa.append(self.goidx[igt])
                except KeyError as ke:
                    self.log.debug(f'no key for {igt}')
            
            gto.is_a = newisa    

    def _df_from_cache(self):
        if os.path.exists(self.cachepath):
            self.log.debug("Trying read from cachepath %s" % self.cachepath)
            self.df = pd.read_csv(self.cachepath, index_col=0)
            self.log.debug("Loaded DF from %s" % self.cachepath)
        else:
            self.log.debug("No file at cachepath: %s" % self.cachepath)


    def get_df(self, usecache=False):
        '''
        Create dataframe from local file for further usage...
        Used cached version on disk if available. 
        '''
        if usecache:
            if self.df is not None:
                self.log.debug("Cache hit. Using DataFrame from cache...")
        else:
            self.log.debug("Cache miss or force. Regenerating DF...")
            data = self.get_lol()
            #data = self.get_dict()
            #df = pd.DataFrame.from_dict(data, orient='index', columns=['goterm','goname','goaspect']) 
            df = pd.DataFrame(data, columns=['goterm','goname','goaspect'] )
            df.set_index('goterm')
            df.to_csv(self.cachepath)
            self.df = df
        self.log.debug(str(self.df))
        return self.df

    def get_lol(self):
        '''
        get info in list of lists form. 
        ''' 
        dt = []
        for gto in self.goidx.values():
            dt.append( [ gto.goterm, gto.goname, gto.goaspect ]) 
        return dt        

    def get_dict(self):
        '''
        get info in dictionary form. 
        ''' 
        dt = {}
        for gto in self.goidx.values():
            dt[gto.goterm] = [ gto.goterm, gto.goname, gto.goaspect ] 
        return dt

    def __repr__(self):
        s = "GeneOntology:"
        for gto in self.goidx.values():
            s +=f'{gto}\n'
        #for atr in ['goterm','goname', 'goaspect']:
        #    s += " %s=%s" % (atr, self.__getattribute__(atr))
        return s

    
    @classmethod
    def get_default_df(cls):
        cp = ConfigParser()
        cp.read(os.path.expanduser('~/git/ontology/etc/ontology.conf'))
        go = GeneOntology(cp)
        df = go.get_df()
        return df         

    @classmethod 
    def get_default_go(cls):
        cp = ConfigParser()
        cp.read(os.path.expanduser('~/git/ontology/etc/ontology.conf'))
        go = GeneOntology(cp)
        return go


if __name__ == '__main__':
    FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(filename)s:%(lineno)d %(name)s.%(funcName)s(): %(message)s'
    logging.basicConfig(format=FORMAT)
    
    parser = argparse.ArgumentParser()
      
    parser.add_argument('-d', '--debug', 
                        action="store_true", 
                        dest='debug', 
                        help='debug logging')

    parser.add_argument('-v', '--verbose', 
                        action="store_true", 
                        dest='verbose', 
                        help='verbose logging')

    parser.add_argument('-c', '--config', 
                        action="store", 
                        dest='conffile', 
                        default='~/git/cshlwork/etc/ontology.conf',
                        help='Config file path [~/git/cshlwork/etc/ontology.conf]')

    parser.add_argument('-t', '--test', 
                        action="store_true", 
                        dest='test', 
                        help='run tests')

    parser.add_argument('-g', '--gaffile', 
                        action="store", 
                        dest='gaffile', 
                        default=None,
                        help='Gene annotation file')

    parser.add_argument('-O', '--obofile', 
                        action="store", 
                        dest='obofile', 
                        default=None,
                        help='GO OBO File.')
    
    parser.add_argument('-o','--outfile', 
                    metavar='outfile',
                    required=False,
                    default=None, 
                    type=str, 
                    help='PDF plot out file. "simple_plots.pdf" if not given')  

    parser.add_argument('goterms', 
                        metavar='goterms', 
                        type=str, 
                        help='one or more space-separated goterms GO:XXXXXXXX' ,
                        nargs='*')
                   
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)    

    config = ConfigParser()
    logging.debug(f'config file = {args.conffile}') 
    config.read(os.path.expanduser(args.conffile))
       
    if args.test:
        test(config)
    
    if len(args.goterms) > 0:
        go = GeneOntology(config, obofile=args.obofile)
        logging.debug(go)
        goidx = go.get_tree_termindex()
        for gt in args.goterms:
            gobj = goidx[gt]
            print(f"{gt} -> {gobj.is_a} ")
    else:
        go = GeneOntology(config, obofile=args.obofile)
        goidx = go.get_tree_termindex()
        s = pprint.pformat(goidx, indent=4)
        print(s  )
        
    #if args.gaffile is not None:
    #    gm = GOMatrix(config, args.gaffile)
    #    df = gm.get_df()
    #    print(df)
    