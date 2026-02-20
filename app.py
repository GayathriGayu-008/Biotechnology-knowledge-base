"""
Biotechnology Search Web Application
Flask backend for searching NCBI and UniProt databases
"""

from flask import Flask, render_template, request, jsonify
import requests
import json
import time
from urllib.parse import quote, urlencode

app = Flask(__name__)

# NCBI E-utilities base URL
NCBI_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
# UniProt API base URL
UNIPROT_BASE_URL = "https://rest.uniprot.org"


def search_ncbi(query, db="protein", max_results=10):
    """
    Search NCBI database using E-utilities
    """
    try:
        # First, use esearch to get IDs
        search_url = f"{NCBI_BASE_URL}/esearch.fcgi"
        params = {
            "db": db,
            "term": query,
            "retmax": max_results,
            "retmode": "json",
            "usehistory": "y"
        }
        
        response = requests.get(search_url, params=params, timeout=30)
        response.raise_for_status()
        search_data = response.json()
        
        id_list = search_data.get("esearchresult", {}).get("idlist", [])
        
        if not id_list:
            return {"status": "success", "results": [], "message": "No results found in NCBI"}
        
        # Now fetch summary for each ID
        summary_url = f"{NCBI_BASE_URL}/esummary.fcgi"
        summary_params = {
            "db": db,
            "id": ",".join(id_list),
            "retmode": "json"
        }
        
        summary_response = requests.get(summary_url, params=summary_params, timeout=30)
        summary_response.raise_for_status()
        summary_data = summary_response.json()
        
        results = []
        for uid, info in summary_data.get("result", {}).items():
            if uid == "uids":
                continue
            results.append({
                "uid": uid,
                "title": info.get("title", "N/A"),
                "accession": info.get("accessionversion", "N/A"),
                "length": info.get("length", "N/A"),
                "create_date": info.get("createdate", "N/A"),
                "update_date": info.get("updatedate", "N/A"),
                "source": "NCBI"
            })
        
        return {"status": "success", "results": results, "total": len(results)}
        
    except requests.exceptions.RequestException as e:
        return {"status": "error", "message": f"NCBI API error: {str(e)}"}
    except Exception as e:
        return {"status": "error", "message": f"Error searching NCBI: {str(e)}"}


def search_uniprot(query, max_results=10):
    """
    Search UniProt database using their REST API
    """
    try:
        # Use UniProt search API
        search_url = f"{UNIPROT_BASE_URL}/search"
        params = {
            "query": query,
            "format": "json",
            "size": max_results,
            "fields": "accession,gene_names,gene_primary,protein_name,organism_name,organism_id,protein_sequence_length,sequence_version"
        }
        
        response = requests.get(search_url, params=params, timeout=30)
        response.raise_for_status()
        
        data = response.json()
        results = []
        
        for entry in data.get("results", []):
            results.append({
                "accession": entry.get("primaryAccession", "N/A"),
                "protein_name": entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value", "N/A"),
                "gene_name": entry.get("genes", [{}])[0].get("geneName", {}).get("value", "N/A") if entry.get("genes") else "N/A",
                "organism": entry.get("organism", {}).get("scientificName", "N/A"),
                "sequence_length": entry.get("sequence", {}).get("length", "N/A"),
                "source": "UniProt"
            })
        
        return {"status": "success", "results": results, "total": len(results)}
        
    except requests.exceptions.RequestException as e:
        return {"status": "error", "message": f"UniProt API error: {str(e)}"}
    except Exception as e:
        return {"status": "error", "message": f"Error searching UniProt: {str(e)}"}


def search_ncbi_nucleotide(query, max_results=10):
    """
    Search NCBI Nucleotide database
    """
    return search_ncbi(query, db="nucleotide", max_results=max_results)


def search_ncbi_genes(query, max_results=10):
    """
    Search NCBI Gene database
    """
    return search_ncbi(query, db="gene", max_results=max_results)


@app.route('/')
def index():
    """Home page with search form"""
    return render_template('index.html')


@app.route('/search', methods=['POST'])
def search():
    """Handle search requests"""
    query = request.form.get('query', '').strip()
    databases = request.form.getlist('databases')
    max_results = int(request.form.get('max_results', 10))
    
    if not query:
        return render_template('index.html', error="Please enter a search term")
    
    if not databases:
        databases = ['ncbi_protein', 'uniprot']
    
    all_results = {}
    
    # Search selected databases
    if 'ncbi_protein' in databases:
        all_results['NCBI Protein'] = search_ncbi(query, db="protein", max_results=max_results)
    
    if 'ncbi_nucleotide' in databases:
        all_results['NCBI Nucleotide'] = search_ncbi_nucleotide(query, max_results=max_results)
    
    if 'ncbi_genes' in databases:
        all_results['NCBI Genes'] = search_ncbi_genes(query, max_results=max_results)
    
    if 'uniprot' in databases:
        all_results['UniProt'] = search_uniprot(query, max_results=max_results)
    
    return render_template('results.html', 
                           query=query, 
                           results=all_results,
                           databases=databases)


@app.route('/api/search')
def api_search():
    """API endpoint for programmatic searches"""
    query = request.args.get('q', '').strip()
    database = request.args.get('db', 'all')
    max_results = int(request.args.get('max', 10))
    
    if not query:
        return jsonify({"status": "error", "message": "No query provided"})
    
    results = {}
    
    if database in ['all', 'ncbi_protein']:
        results['ncbi_protein'] = search_ncbi(query, db="protein", max_results=max_results)
    
    if database in ['all', 'ncbi_nucleotide']:
        results['ncbi_nucleotide'] = search_ncbi_nucleotide(query, max_results=max_results)
    
    if database in ['all', 'ncbi_genes']:
        results['ncbi_genes'] = search_ncbi_genes(query, max_results=max_results)
    
    if database in ['all', 'uniprot']:
        results['uniprot'] = search_uniprot(query, max_results=max_results)
    
    return jsonify(results)


@app.route('/api/health')
def health_check():
    """Health check endpoint"""
    return jsonify({"status": "healthy", "service": "Biotech Search API"})


if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)
