# 🧬 Biotechnology Database Search Application

A Flask-based web application for searching biotechnology databases including NCBI (Protein, Nucleotide, Genes) and UniProt.

## Features

- **Multi-database Search**: Search across multiple biotechnology databases simultaneously
- **NCBI Integration**: Search Protein, Nucleotide, and Gene databases via NCBI E-utilities
- **UniProt Integration**: Search protein sequences and annotations via UniProt REST API
- **User-friendly Interface**: Modern, responsive web interface with gradient design
- **Quick Search Links**: Pre-configured search examples for common biotechnology terms
- **API Endpoint**: Programmatic access via REST API

## Prerequisites

- Python 3.8 or higher
- pip (Python package manager)

## Installation

1. Navigate to the project directory:
   
```
bash
   cd biotech_search_app
   
```

2. Create a virtual environment (optional but recommended):
   
```
bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   
```

3. Install dependencies:
   
```
bash
   pip install -r requirements.txt
   
```

## Running the Application

1. Start the Flask server:
   
```
bash
   python app.py
   
```

2. Open your web browser and navigate to:
   
```
   http://localhost:5000
   
```

## Usage

### Web Interface

1. Enter a search term (gene name, protein name, keyword, etc.)
2. Select which databases to search:
   - NCBI Protein
   - NCBI Nucleotide  
   - NCBI Genes
   - UniProt
3. Choose the maximum number of results to display
4. Click "Search Databases" to view results

### Quick Search Examples

- Insulin
- BRCA1
- Kinase
- Enzyme
- CRISPR
- Hemoglobin

### API Usage

You can also access the search functionality programmatically:

```
bash
# Search all databases
curl "http://localhost:5000/api/search?q=insulin"

# Search specific database
curl "http://localhost:5000/api/search?q=BRCA1&db=uniprot"

# Limit results
curl "http://localhost:5000/api/search?q=kinase&max=5"

# Health check
curl "http://localhost:5000/api/health"
```

## Project Structure

```
biotech_search_app/
├── app.py                 # Main Flask application
├── requirements.txt       # Python dependencies
├── templates/
│   ├── index.html        # Home/search page
│   └── results.html     # Results display page
├── static/
│   └── style.css        # Additional styles
└── README.md            # This file
```

## Supported Search Terms

You can search for:
- Gene names (e.g., BRCA1, TP53, MYC)
- Protein names (e.g., insulin, hemoglobin, kinase)
- Enzyme classifications (e.g., oxidoreductase, transferase)
- Keywords (e.g., CRISPR, antibiotic, enzyme)
- Accession numbers (e.g., P01308, NP_000517)

## API Response Format

```
json
{
  "NCBI Protein": {
    "status": "success",
    "results": [
      {
        "uid": "123456",
        "title": "protein name",
        "accession": "XP_123456",
        "length": 500,
        "source": "NCBI"
      }
    ],
    "total": 1
  },
  "UniProt": {
    "status": "success",
    "results": [
      {
        "accession": "P01308",
        "protein_name": "Insulin",
        "gene_name": "INS",
        "organism": "Homo sapiens",
        "sequence_length": 110,
        "source": "UniProt"
      }
    ],
    "total": 1
  }
}
```

## License

This project is for educational and research purposes.

## Acknowledgments

- [NCBI E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25501/)
- [UniProt REST API](https://www.uniprot.org/help/api)
