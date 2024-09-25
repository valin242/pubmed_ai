"""
PubMed Language Model for iLabNote

PubMed Transform: Transform data for NLP
@author(s): Valinteshley Pierre

References
Pubmed
License argreements
"""
from pubmedscraper import RetrieveFromPub as rfp
import pandas as pd
# import torch
# import torchtext
# from transformers import GPT2Tokenizer, GPT2Model
# from torch.utils.data import Dataset, DataLoader
# from torch.optim import Adam


# Retrieve papers on cardiac tissue engineering to train
# Explanation: We'll create a function to get user input for the search query
# This allows for more flexibility in searching different topics

def get_user_query():
    """
    Prompts the user to enter a search query for PubMed articles.
    
    Returns:
    str: The user's input query
    """
    print("Welcome to the PubMed article search tool!")
    query = input("Please enter the topic you'd like to search for: ")
    return query

# Get the user's search query
user_query = get_user_query()

# Create an instance of RetrieveFromPub with the user's query
pubmed_retriever = rfp(user_query)

# TODO: Implement error handling for invalid queries
# TODO: Add option for user to refine or change their search query

# Search for studies based on the user's query
studies_list = pubmed_retriever.search()

print(f"Found {len(studies_list)} studies related to '{user_query}'")

# TODO: Implement functionality to display a summary of the search results
# TODO: Add option for user to retrieve full text of selected articles

def retrieve_and_process_articles(pubmed_retriever):
    """
    Retrieves and processes articles from PubMed.
    
    Args:
    pubmed_retriever (RetrieveFromPub): An instance of the RetrieveFromPub class.
    
    Returns:
    pandas.DataFrame: A DataFrame containing the retrieved articles.
    """
    print("Retrieving articles...")
    try:
        studies_table, no_access = pubmed_retriever.retrieve_text()
        print(f"Retrieved {len(studies_table)} articles. {len(no_access)} articles were not accessible.")
        
        # Debug: Print the type and shape of studies_table
        print(f"Type of studies_table: {type(studies_table)}")
        if isinstance(studies_table, pd.DataFrame):
            print(f"Shape of studies_table: {studies_table.shape}")
            print("\nFirst few rows of studies_table:")
            print(studies_table.head())
        else:
            print("studies_table is not a pandas DataFrame")
        
        return studies_table
    except Exception as e:
        print(f"An error occurred: {str(e)}")
        return None

def main():
    """
    Main function to run the PubMed article retrieval and processing pipeline.
    """
    user_query = get_user_query()
    pubmed_retriever = rfp(user_query)
    studies_list = pubmed_retriever.search()
    print(f"Found {len(studies_list)} studies related to '{user_query}'")
    
    studies_table = retrieve_and_process_articles(pubmed_retriever)
    
    print("Sample of retrieved articles:")
    print(studies_table.head())
    
    # TODO: Implement further processing of the retrieved articles
    # TODO: Add functionality to save the retrieved articles to a file

if __name__ == "__main__":
    main()

# /todo: Implement error handling for network issues or API failures
# /todo: Add logging to track the retrieval process and any issues encountered

# Fine tuning GPT-2 model or knowledge base embedding
# tokenizer = GPT2Tokenizer.from_pretrained('gpt2')
# model = GPT2Model.from_pretrained('gpt2')
# text = '' # pubmed articles
# encoded_input = tokenizer(text, return_tensors='pt')
# output = model(**encoded_input)


"""
Citations:
@article{radford2019language,
  title={Language Models are Unsupervised Multitask Learners},
  author={Radford, Alec and Wu, Jeff and Child, Rewon and Luan, David and Amodei, Dario and Sutskever, Ilya},
  year={2019}
}

"""