# Scientific Text Generation with GPT-2 (iLabNote)

## Overview

This project focuses on training the GPT-2 language model on scientific concepts extracted from PubMed articles. The goal is to generate text relevant to specific scientific domains based on the preexisting knowledge in PubMed.

## Technologies and Tools

This project utilizes a variety of technologies and tools for data processing, model training, and deployment:

- **Python**: Primary programming language for data preprocessing, training, and inference.
- **Pandas**: Used for efficient data manipulation and analysis.
- **PyTorch**: Deep learning framework for fine-tuning the GPT-2 model.
- **Java**: Employed for specific tasks or integrations as needed.
- **HTML/CSS**: Utilized for building any web-based components or interfaces.
- **JavaScript**: Applied for dynamic and interactive elements in web applications.

## Steps to Replicate

### 1. Data Collection

#### a. PubMed API

Utilize the PubMed API to retrieve scientific papers related to your specific concepts. Filter papers based on keywords, authors, journals, etc.

#### b. Data Preprocessing

Extract text content from the retrieved papers and clean the text by removing irrelevant information, such as headers, footers, and non-textual elements.

### 2. Data Formatting

#### a. Tokenization

Tokenize the text into smaller units, such as words or subwords. Ensure consistent tokenization with the GPT-2 pretraining.

#### b. Sequence Length

Divide the text into sequences of a fixed length compatible with GPT-2 processing.

### 3. Model Fine-Tuning

#### a. Choose a GPT-2 Variant

Select the GPT-2 model variant based on computational resources and dataset size.

#### b. Set Up Environment

Create a suitable environment for training using frameworks like TensorFlow or PyTorch.

#### c. Fine-Tuning Script

Adapt the GPT-2 fine-tuning script to your needs.

#### d. Loss Function and Optimization

Define a suitable loss function and optimization algorithm. Adjust hyperparameters like learning rate.

#### e. Train the Model

Train the GPT-2 model on the preprocessed dataset. This can be a time-intensive process and requires robust hardware.

### 4. Evaluation and Iteration

#### a. Validation Set

Set aside a portion of data for validation to monitor the model's performance during training.

#### b. Evaluation Metrics

Define metrics to assess the quality of the generated text.

#### c. Iterative Refinement

Fine-tune the model iteratively based on evaluation results. Adjust hyperparameters or increase the training data if needed.

### 5. Model Deployment (Optional)

#### a. Save the Model

Save the trained model and its configuration.

#### b. Inference

Set up an inference pipeline to generate text based on the fine-tuned GPT-2 model.

## Key Skills

This project involves the application of key skills relevant to software engineering and data science:

- **Natural Language Processing (NLP)**
- **Deep Learning**
- **Data Preprocessing**
- **Model Evaluation and Optimization**
- **Web Development (HTML, CSS, JavaScript)**

## Usage

Once the model is trained and validated, you can use it to generate scientific text relevant to your specified concepts.

## Contributing

Contributions to this project are welcome. If you encounter issues or have suggestions for improvements, feel free to open an issue or submit a pull request.

## License

This project is licensed under the [MIT License](LICENSE).


## Citaions
1. Comeau DC, Wei CH, Islamaj Doğan R, and Lu Z. PMC text mining subset in BioC: about 3 million full text articles and growing, Bioinformatics, btz070, 2019.
2. PMC Open Access Subset [Internet]. Bethesda (MD): National Library of Medicine. 2003 - [cited 2023 12 04]. Available from https://www.ncbi.nlm.nih.gov/pmc/tools/openftlist/.



# LabNote AI (PubMed AI)
Scrape Pubmed articles with the key words 'Cardiac Tissue Engineering',  'Biomaterial', 'dECM' and train a GPT2 model on the paper. The goal is to be able to write a review paper based on information scraped from biomaterials paper on cardiac tissue engineering.