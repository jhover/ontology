#!/usr/bin/env python
#
#      

import argparse
from configparser import ConfigParser
import logging
import os
import sys

gitpath=os.path.expanduser("~/git/ontology")
sys.path.append(gitpath)

from ontology.ontology import *
  

if __name__ == '__main__':
    # defaults
    OBOFILE=os.path.expanduser('~/data/ontology/go.obo') 
    GAFFILE=os.path.expanduser('~/data/ontology/Mus_musculus/mgi.gaf')
    CONFFILE=os.path.expanduser('~/git/ontology/etc/ontology.conf')
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

    parser.add_argument('-b', '--obofile', 
                        action="store", 
                        dest='obofile', 
                        default=OBOFILE,
                        help='Gene ontology OBO file.')
    
    parser.add_argument('-g', '--gaffile', 
                        action="store", 
                        dest='gaffile', 
                        default=GAFFILE,
                        help='GAF file. ')

    parser.add_argument('-C', '--usecache',
                        action='store_true', 
                        dest='nocache',
                        default=False, 
                        help='Use cached information.' )


    parser.add_argument('-o', '--outfile', 
                        action="store", 
                        dest='outfile', 
                        default='genegomatrix.tsv',
                        help='One-hot matrix file. ')

   
    parser.add_argument('-c', '--config', 
                        action="store", 
                        dest='conffile', 
                        default=CONFFILE,
                        help='Config file path [~/etc/gotool.conf]') 
                    
    #parser.add_argument('goterms', 
    #                    metavar='goterms', 
    #                    type=str, 
    #                    help='one or more space-separated goterms GO:XXXXXXXX' ,
    #                    nargs='*'
    #               )
    
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)
    
    cp = ConfigParser()
    cp.read(args.conffile)
    cdict = format_config(cp)
    logging.debug(f'Running with config={args.conffile}:\n{cdict}')

    make_hotmatrix(obofile=args.obofile, gaffile=args.gaffile, outfile=args.outfile, usecache=args.nocache)