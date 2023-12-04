"""
PubMed Language Model for iLabNote

PubMed Scraper: Retrieve papers with a specific search match
@author(s): Valinteshley Pierre

References
Pubmed
License argreements
"""

from Bio import Entrez
import pandas as pd
import numpy as np
import random
import urllib.request, json


# get list of IDs of studies that match your query
def search(query):
    Entrez.email = 'vdp14@case.edu'
    handle = Entrez.esearch(db='pubmed', sort='relevance', retmax='250000', retmode='xml', term=query)
    results = Entrez.read(handle)
    studiesIdList = results['IdList']
    return studiesIdList

# For Example: search for cardiac tissue engineering papers
# q = 'decellularized ECM for cardiac tissue engineering'
# studies = search(q)
# studiesIdList = studies['IdList']
# print(len(studiesIdList))

# If you want to reduce the number of papers for easy computation
# studiesIdList_shortened = random.sample(studiesIdList, 300)


# # Use efetch to get the details of each study (this only fetches the abstracts)
# ****In the future a function that suggests relevant papers based on abstracts****
# def fetch_details(id_list):
#     ids = ','.join(id_list)
#     Entrez.email = 'vdp14@case.edu'
#     handle = Entrez.efetch(db='pubmed', retmode='xml', id=ids)
#     results = Entrez.read(handle)
#     return results

# # create a pandas df that has the pubmed article info (title, abstract, journal, date, etc.)
# # Efetch runs 10,000 studies max, so separate idList in sections
# titles = []
# abstracts = []
# journals = []
# pubdate_years = []
# pudate_months = []

# studies = fetch_details(studiesIdList)
# chuncks_size = 10000

'''
Challenge:
Not all full text articles are available so we will have to write a function that retrieves free full text articles. There will probably be another list that has texts
that were not able to be retrieve that we will then feed PDFs to. The function will break down the PDF into similar format as the .xml files for consistency.
'''

# URL to retrieve publically available articles
# format_json = 'json'
# id_list = studiesIdList
# encoding_uni = 'unicode'


def retrieve_text(id):

    # Dataframe to store retrieved data
    studies_table = pd.DataFrame(columns=['id', 'title', 'body_content'])

    for i in id: # pull the json docs from the list of pubmed IDs
        try:
            url_bioc = 'https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_json' + '/' + i + '/' + 'unicode'
            response = urllib.request.urlopen(url_bioc)
            data = json.loads(response.read())

            # Parse out the body of text from the json file
            body_size = len(data['documents'][0]['passages'])
            body_content = ''
            for i in range(body_size):
                body_content = body_content + data['documents'][0]['passages'][i]['text']

            body_content = body_content.replace('\t', '')

            title = data['documents'][0]['passages'][0]['text']

            pmid = data['documents'][0]['passages'][0]['infons']['article-id_pmid']

            content_list = [pmid, title, body_content]

            # store the data in a dataframe
            studies_table.loc[len(studies_table)] = content_list
        except:
            no_access = []
            no_access.append(i)
            # print('No open access file for ', i)

    return studies_table
    
# raw_table = retrieve_text(format_json, id_list, encoding_uni)
# print('Size of the tables is ', raw_table.size)
        
