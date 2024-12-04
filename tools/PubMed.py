from Bio import Entrez
from Bio import Medline
import pandas as pd

def fetch_pubmed_data(query, max_results=10):
    # Set your email to use Entrez
    Entrez.email = "your_email@example.com"  # Replace with your email

    # Search PubMed for articles related to the query
    search_handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results, usehistory="y")
    search_results = Entrez.read(search_handle)
    search_handle.close()

    # Get the list of PubMed IDs
    id_list = search_results["IdList"]
    if not id_list:
        return []

    # Fetch the details of the articles
    fetch_handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="text")
    records = Medline.parse(fetch_handle)
    articles = []

    for record in records:
        article_info = {
            "title": record.get("TI", "No title available"),
            "abstract": record.get("AB", "No abstract available"),
            "pub_date": record.get("DP", "No publication date available"),
        }
        articles.append(article_info)

    fetch_handle.close()
    return articles


def get_PubMed_text(smiles_list, csv_name, max_results=5):
    data = {
        'smiles': smiles_list,
        'pubmed_text': []
    }
    for name in smiles_list:
        results = fetch_pubmed_data(name, max_results)
        text = ""
        for idx, article in enumerate(results):
            text = text + f"Article {idx + 1}: Title: {article['title']} Abstract: {article['abstract']} Publication Date: {article['pub_date']} "
        data["pubmed_text"].append(text)
    df = pd.DataFrame(data)
    df.to_csv(csv_name, index=False)


smiles_list = ["Aspirin"]
csv_name = "test.csv"
get_PubMed_text(smiles_list, csv_name)