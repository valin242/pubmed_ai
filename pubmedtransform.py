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
import torch
import torchtext
from transformers import GPT2Tokenizer, GPT2Model
from torch.utils.data import Dataset, DataLoader
from torch.optim import Adam


# Retrieve papers on cardiac tissue engineering to train
studies_list = rfp('cardiac tissue engineering')

# Fine tuning GPT-2 model or knowledge base embedding
tokenizer = GPT2Tokenizer.from_pretrained('gpt2')
model = GPT2Model.from_pretrained('gpt2')
text = '' # pubmed articles
encoded_input = tokenizer(text, return_tensors='pt')
output = model(**encoded_input)


"""
Citations:
@article{radford2019language,
  title={Language Models are Unsupervised Multitask Learners},
  author={Radford, Alec and Wu, Jeff and Child, Rewon and Luan, David and Amodei, Dario and Sutskever, Ilya},
  year={2019}
}

"""