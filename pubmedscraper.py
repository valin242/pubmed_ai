from Bio import Entrez
import pandas as pd
import numpy as np


# get list of IDs of studies that match your query
def search(query):
    Entrez.email = 'vdp14@case.edu'
    handle = Entrez.esearch(db='pubmed', sort='relevance', retmax='250000', retmode='xml', term=query)
    results = Entrez.read(handle)
    return results

# search for cardiac tissue engineering papers
q = 'decellularized ECM for cardiac tissue engineering'
studies = search(q)
studiesIdList = studies['IdList']

# Use efetch to get the details of each study
def fetch_details(id_list):
    ids = ','.join(id_list)
    Entrez.email = 'vdp14@case.edu'
    handle = Entrez.efetch(db='pubmed', retmode='xml', id=ids)
    results = Entrez.read(handle)
    return results

# create a pandas df that has the pubmed article info (title, abstract, journal, date, etc.)
# Efetch runs 10,000 studies max, so separate idList in sections
titles = []
abstracts = []
journals = []
pubdate_years = []
pudate_months = []

studies = fetch_details(studiesIdList)
chuncks_size = 10000