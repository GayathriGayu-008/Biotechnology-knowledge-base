"""
Biotechnology Knowledge Base - Data Collection Module
Collects terms from NCBI, UniProt, and other databases
"""

import requests
import json
import time
import sqlite3
from datetime import datetime
from urllib.parse import quote, urlencode

# Database configuration
DATABASE_PATH = "knowledge_base.db"

# API endpoints
NCBI_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
UNIPROT_BASE_URL = "https://rest.uniprot.org"
MESH_URL = "https://id.nlm.nih.gov/mesh/vocab/2024"
GO_URL = "http://purl.obolibrary.org/obo/go/go-basic.obo"

# Category definitions
CATEGORIES = {
    "molecular_biology": ["gene", "protein", "DNA", "RNA", "enzyme", "kinase"],
    "genetics": ["chromosome", "allele", "mutation", "genome", "genotype"],
    "biochemistry": ["metabolite", "pathway", "reaction", "compound"],
    "cell_biology": ["cell", "organelle", "membrane", "cytoplasm"],
    "microbiology": ["bacteria", "virus", "fungus", "pathogen"],
    "medicine": ["disease", "disorder", "syndrome", "treatment", "drug"],
    "biotechnology": ["cloning", "PCR", "CRISPR", "vector", "expression"],
    "bioinformatics": ["sequence", "database", "alignment", "annotation"]
}

class KnowledgeBaseCollector:
    def __init__(self, db_path=DATABASE_PATH):
        self.db_path = db_path
        self.init_database()
        
    def init_database(self):
        """Initialize the knowledge base database"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        # Main terms table
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS terms (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                term TEXT NOT NULL,
                term_type TEXT,
                category TEXT,
                definition TEXT,
                function TEXT,
                importance TEXT,
                source TEXT,
                external_id TEXT,
                external_link TEXT,
                metadata TEXT,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        ''')
        
        # Search index table
        cursor.execute('''
            CREATE VIRTUAL TABLE IF NOT EXISTS terms_fts USING fts5(
                term, definition, function, category,
                content='terms',
                content_rowid='id'
            )
        ''')
        
        # Collection status table
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS collection_status (
                source TEXT PRIMARY KEY,
                last_collection TIMESTAMP,
                terms_collected INTEGER,
                status TEXT
            )
        ''')
        
        # Categories table
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS categories (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                name TEXT UNIQUE NOT NULL,
                description TEXT,
                parent_id INTEGER,
                FOREIGN KEY (parent_id) REFERENCES categories(id)
            )
        ''')
        
        # Relationships table
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS term_relationships (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                term1_id INTEGER,
                term2_id INTEGER,
                relationship_type TEXT,
                FOREIGN KEY (term1_id) REFERENCES terms(id),
                FOREIGN KEY (term2_id) REFERENCES terms(id)
            )
        ''')
        
        conn.commit()
        conn.close()
        print("Database initialized successfully!")
    
    def add_term(self, term, term_type=None, category=None, definition=None, 
                 function=None, importance=None, source=None, external_id=None,
                 external_link=None, metadata=None):
        """Add a single term to the knowledge base"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        try:
            cursor.execute('''
                INSERT OR REPLACE INTO terms 
                (term, term_type, category, definition, function, importance, 
                 source, external_id, external_link, metadata, updated_at)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (term, term_type, category, definition, function, importance,
                  source, external_id, external_link, metadata, datetime.now()))
            
            term_id = cursor.lastrowid
            conn.commit()
            conn.close()
            
            # Update FTS index
            self.update_fts_index(term_id, term, definition, function, category)
            
            return term_id
        except Exception as e:
            conn.close()
            print(f"Error adding term {term}: {e}")
            return None
    
    def update_fts_index(self, term_id, term, definition, function, category):
        """Update the full-text search index"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        try:
            # Delete old FTS entry
            cursor.execute('DELETE FROM terms_fts WHERE rowid = ?', (term_id,))
            
            # Add new FTS entry
            cursor.execute('''
                INSERT INTO terms_fts (rowid, term, definition, function, category)
                VALUES (?, ?, ?, ?, ?)
            ''', (term_id, term, definition or "", function or "", category or ""))
            
            conn.commit()
        except Exception as e:
            print(f"FTS update error: {e}")
        finally:
            conn.close()
    
    def add_category(self, name, description=None, parent_id=None):
        """Add a category"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        cursor.execute('''
            INSERT OR IGNORE INTO categories (name, description, parent_id)
            VALUES (?, ?, ?)
        ''', (name, description, parent_id))
        
        conn.commit()
        conn.close()
    
    def search(self, query, limit=10):
        """Search the knowledge base"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        cursor.execute('''
            SELECT id, term, term_type, category, definition, function, importance, source
            FROM terms
            WHERE term LIKE ? OR definition LIKE ?
            LIMIT ?
        ''', (f'%{query}%', f'%{query}%', limit))
        
        results = cursor.fetchall()
        conn.close()
        
        return [{
            'id': r[0], 'term': r[1], 'term_type': r[2], 'category': r[3],
            'definition': r[4], 'function': r[5], 'importance': r[6], 'source': r[7]
        } for r in results]
    
    def get_total_terms(self):
        """Get total number of terms"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        cursor.execute('SELECT COUNT(*) FROM terms')
        count = cursor.fetchone()[0]
        conn.close()
        return count

    # ============ NCBI Data Collection ============
    
    def collect_ncbi_taxonomy(self, max_taxa=10000):
        """Collect organism data from NCBI Taxonomy"""
        print(f"Collecting NCBI Taxonomy (max: {max_taxa})...")
        
        try:
            # Get list of common organisms
            common_taxids = [
                9606,  # Homo sapiens
                10090, # Mus musculus
                10116, # Rattus norvegicus
                7227,  # Drosophila melanogaster
                559292, # Saccharomyces cerevisiae
                9823,  # Sus scrofa
                9913,  # Bos taurus
                836,   # Xenopus laevis
                7955,  # Danio rerio
                3702,  # Arabidopsis thaliana
                4577,  # Zea mays
                4530,  # Oryza sativa
                4896,  # Schizosaccharomyces pombe
                4932,  # Saccharomyces cerevisiae
                6239,  # Caenorhabditis elegans
                86679, # Listeria monocytogenes
                511145, # Escherichia coli str. K-12
                224308, # Bacillus anthracis
                169963, # Staphylococcus aureus
                83333,  # Mycobacterium tuberculosis
            ]
            
            collected = 0
            for taxid in common_taxids[:max_taxa]:
                try:
                    url = f"{NCBI_BASE_URL}/esummary.fcgi"
                    params = {"db": "taxonomy", "id": taxid, "retmode": "json"}
                    
                    response = requests.get(url, params=params, timeout=30)
                    if response.status_code == 200:
                        data = response.json()
                        result = data.get("result", {})
                        
                        if str(taxid) in result:
                            info = result[str(taxid)]
                            self.add_term(
                                term=info.get("scientificname", ""),
                                term_type="organism",
                                category="life_science",
                                definition=f"Taxonomy ID: {taxid}. {info.get('description', '')}",
                                importance=f"Common model organism, ranked {info.get('rank', 'N/A')}",
                                source="NCBI Taxonomy",
                                external_id=str(taxid),
                                external_link=f"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={taxid}",
                                metadata=json.dumps({"taxid": taxid, "rank": info.get("rank")})
                            )
                            collected += 1
                            print(f"Collected: {info.get('scientificname')}")
                except Exception as e:
                    print(f"Error collecting taxonomy {taxid}: {e}")
                
                time.sleep(0.1)  # Rate limiting
            
            self.update_collection_status("ncbi_taxonomy", collected)
            print(f"NCBI Taxonomy collection complete: {collected} terms")
            return collected
            
        except Exception as e:
            print(f"NCBI Taxonomy collection error: {e}")
            return 0
    
    def collect_ncbi_genes(self, search_terms=None, max_per_term=100):
        """Collect gene data from NCBI Gene database"""
        if search_terms is None:
            search_terms = [
                "kinase", "receptor", "transporter", "enzyme",
                "hormone", "growth factor", "immune", "neural"
            ]
        
        print(f"Collecting NCBI Genes (max {max_per_term} per term)...")
        collected = 0
        
        for search_term in search_terms:
            try:
                # Search for genes
                search_url = f"{NCBI_BASE_URL}/esearch.fcgi"
                params = {
                    "db": "gene",
                    "term": search_term,
                    "retmax": max_per_term,
                    "retmode": "json"
                }
                
                response = requests.get(search_url, params=params, timeout=30)
                if response.status_code == 200:
                    data = response.json()
                    id_list = data.get("esearchresult", {}).get("idlist", [])
                    
                    if id_list:
                        # Get summaries
                        summary_url = f"{NCBI_BASE_URL}/esummary.fcgi"
                        summary_params = {
                            "db": "gene",
                            "id": ",".join(id_list),
                            "retmode": "json"
                        }
                        
                        summary_response = requests.get(summary_url, params=summary_params, timeout=30)
                        if summary_response.status_code == 200:
                            summary_data = summary_response.json()
                            
                            for uid, info in summary_data.get("result", {}).items():
                                if uid == "uids":
                                    continue
                                
                                self.add_term(
                                    term=info.get("name", f"Gene {uid}"),
                                    term_type="gene",
                                    category="genetics",
                                    definition=info.get("description", ""),
                                    importance=f"Official Symbol: {info.get('name', 'N/A')}",
                                    source="NCBI Gene",
                                    external_id=str(uid),
                                    external_link=f"https://www.ncbi.nlm.nih.gov/gene/{uid}",
                                    metadata=json.dumps(info)
                                )
                                collected += 1
                
                time.sleep(0.2)  # Rate limiting
                
            except Exception as e:
                print(f"Error collecting genes for {search_term}: {e}")
        
        self.update_collection_status("ncbi_genes", collected)
        print(f"NCBI Genes collection complete: {collected} terms")
        return collected
    
    def collect_ncbi_proteins(self, search_terms=None, max_per_term=50):
        """Collect protein data from NCBI Protein database"""
        if search_terms is None:
            search_terms = [
                "kinase", "receptor", "enzyme", "antibody",
                "hormone", "transport protein", "structural protein"
            ]
        
        print(f"Collecting NCBI Proteins...")
        collected = 0
        
        for search_term in search_terms:
            try:
                search_url = f"{NCBI_BASE_URL}/esearch.fcgi"
                params = {
                    "db": "protein",
                    "term": search_term,
                    "retmax": max_per_term,
                    "retmode": "json"
                }
                
                response = requests.get(search_url, params=params, timeout=30)
                if response.status_code == 200:
                    data = response.json()
                    id_list = data.get("esearchresult", {}).get("idlist", [])
                    
                    if id_list:
                        summary_url = f"{NCBI_BASE_URL}/esummary.fcgi"
                        summary_params = {
                            "db": "protein",
                            "id": ",".join(id_list[:20]),  # Limit per batch
                            "retmode": "json"
                        }
                        
                        summary_response = requests.get(summary_url, params=summary_params, timeout=30)
                        if summary_response.status_code == 200:
                            summary_data = summary_response.json()
                            
                            for uid, info in summary_data.get("result", {}).items():
                                if uid == "uids":
                                    continue
                                
                                self.add_term(
                                    term=info.get("title", f"Protein {uid}")[:200],
                                    term_type="protein",
                                    category="molecular_biology",
                                    definition=f"Accession: {info.get('accessionversion', 'N/A')}. Length: {info.get('length', 'N/A')} aa",
                                    function="Protein with enzymatic or binding function",
                                    importance=f"Title: {info.get('title', 'N/A')[:100]}",
                                    source="NCBI Protein",
                                    external_id=str(uid),
                                    external_link=f"https://www.ncbi.nlm.nih.gov/protein/{uid}",
                                    metadata=json.dumps({"length": info.get("length")})
                                )
                                collected += 1
                
                time.sleep(0.2)
                
            except Exception as e:
                print(f"Error collecting proteins for {search_term}: {e}")
        
        self.update_collection_status("ncbi_proteins", collected)
        print(f"NCBI Proteins collection complete: {collected} terms")
        return collected

    # ============ UniProt Data Collection ============
    
    def collect_uniprot_proteins(self, organism=None, max_results=5000):
        """Collect protein data from UniProt"""
        print(f"Collecting UniProt proteins (max: {max_results})...")
        collected = 0
        
        query = "reviewed:true"
        if organism:
            query += f" organism:{organism}"
        
        try:
            for start in range(0, min(max_results, 500), 50):
                url = f"{UNIPROT_BASE_URL}/search"
                params = {
                    "query": query,
                    "format": "json",
                    "size": 50,
                    "start": start,
                    "fields": "accession,protein_name,gene_names,organism_name,protein_sequence_length,function"
                }
                
                response = requests.get(url, params=params, timeout=30)
                if response.status_code == 200:
                    data = response.json()
                    results = data.get("results", [])
                    
                    if not results:
                        break
                    
                    for entry in results:
                        gene_name = entry.get("genes", [{}])[0].get("geneName", {}).get("value", "N/A") if entry.get("genes") else "N/A"
                        
                        self.add_term(
                            term=entry.get("primaryAccession", ""),
                            term_type="protein",
                            category="molecular_biology",
                            definition=entry.get("function", "")[:500] if entry.get("function") else f"Protein: {entry.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', 'N/A')}",
                            function="Protein function as described in UniProt",
                            importance=f"Organism: {entry.get('organism', {}).get('scientificName', 'N/A')}, Length: {entry.get('sequence', {}).get('length', 'N/A')} aa",
                            source="UniProt",
                            external_id=entry.get("primaryAccession"),
                            external_link=f"https://www.uniprot.org/uniprot/{entry.get('primaryAccession')}",
                            metadata=json.dumps({
                                "gene": gene_name,
                                "organism": entry.get("organism", {}).get("scientificName"),
                                "length": entry.get("sequence", {}).get("length")
                            })
                        )
                        collected += 1
                
                time.sleep(0.1)
                
        except Exception as e:
            print(f"UniProt collection error: {e}")
        
        self.update_collection_status("uniprot", collected)
        print(f"UniProt collection complete: {collected} terms")
        return collected
    
    # ============ Common Biological Terms ============
    
    def add_common_biological_terms(self):
        """Add common biological, medical, and biotechnology terms"""
        print("Adding common biological terms...")
        
        common_terms = [
            # Molecular Biology
            {"term": "DNA", "type": "molecule", "cat": "molecular_biology", "def": "Deoxyribonucleic acid - the hereditary material in humans and almost all other organisms", "func": "Carries genetic instructions for development, functioning, growth and reproduction", "imp": "Fundamental to all known life forms"},
            {"term": "RNA", "type": "molecule", "cat": "molecular_biology", "def": "Ribonucleic acid - present in all living cells", "func": "Plays crucial roles in coding, decoding, regulation and expression of genes", "imp": "Essential for protein synthesis and gene regulation"},
            {"term": "Protein", "type": "molecule", "cat": "molecular_biology", "def": "Large biomolecules consisting of amino acid chains", "func": "Perform a vast array of functions including catalyzing reactions, providing structure, and transporting molecules", "imp": "Essential for the structure, function, and regulation of the body's tissues and organs"},
            {"term": "Gene", "type": "unit", "cat": "genetics", "def": "Unit of heredity that is transferred from parent to offspring", "func": "Contains genetic information that determines characteristics", "imp": "Fundamental to inheritance and genetic diseases"},
            {"term": "Enzyme", "type": "protein", "cat": "biochemistry", "def": "Biological catalyst that speeds up chemical reactions", "func": "Catalyzes metabolic processes essential for life", "imp": "Critical for biochemical reactions in all living organisms"},
            {"term": "Kinase", "type": "enzyme", "cat": "biochemistry", "def": "Enzyme that transfers phosphate groups from ATP to other molecules", "func": "Regulates many cellular processes including metabolism, cell signaling, and apoptosis", "imp": "Important drug targets for cancer and other diseases"},
            {"term": "Receptor", "type": "protein", "cat": "molecular_biology", "def": "Protein that receives chemical signals from outside a cell", "func": "Transduces signals into cellular responses", "imp": "Major targets for drug development"},
            {"term": "Transporter", "type": "protein", "cat": "molecular_biology", "def": "Membrane protein that facilitates movement of molecules across membranes", "func": "Essential for nutrient uptake, waste removal, and ion balance", "imp": "Important for cellular homeostasis"},
            
            # Genetics
            {"term": "Chromosome", "type": "structure", "cat": "genetics", "def": "Thread-like structure of nucleic acids and protein found in the nucleus of cells", "func": "Carries genetic information in the form of genes", "imp": "Defines species and genetic traits"},
            {"term": "Allele", "type": "variant", "cat": "genetics", "def": "One of two or more alternative forms of a gene", "func": "Determines different phenotypic traits", "imp": "Important for genetic variation and inheritance"},
            {"term": "Mutation", "type": "process", "cat": "genetics", "def": "Change in the DNA sequence", "func": "Can lead to changes in protein function or regulation", "imp": "Source of genetic variation, can cause diseases"},
            {"term": "Genome", "type": "set", "cat": "genetics", "def": "Complete set of genetic information in an organism", "func": "Contains all genes and non-coding sequences", "imp": "Defines the complete hereditary information"},
            {"term": "Genotype", "type": "characteristic", "cat": "genetics", "def": "Genetic constitution of an organism", "func": "Determines inherited characteristics", "imp": "Used in genetic testing and breeding"},
            {"term": "Phenotype", "type": "characteristic", "cat": "genetics", "def": "Observable characteristics of an organism", "func": "Result of genotype and environment interaction", "imp": "Used to study genetic traits"},
            
            # Cell Biology
            {"term": "Cell", "type": "unit", "cat": "cell_biology", "def": "Basic structural and functional unit of all living organisms", "func": "Performs all metabolic activities", "imp": "Fundamental unit of life"},
            {"term": "Mitochondria", "type": "organelle", "cat": "cell_biology", "def": "Organelle found in cells that generates most of the cell's ATP", "func": "Produces energy through oxidative phosphorylation", "imp": "Powerhouse of the cell, important for metabolism"},
            {"term": "Ribosome", "type": "organelle", "cat": "cell_biology", "def": "Cellular machine that synthesizes proteins", "func": "Translates mRNA into protein sequences", "imp": "Essential for protein synthesis"},
            {"term": "Nucleus", "type": "organelle", "cat": "cell_biology", "def": "Membrane-bound organelle containing genetic material", "func": "Controls cell activities and contains DNA", "imp": "Command center of the cell"},
            {"term": "Cell Membrane", "type": "structure", "cat": "cell_biology", "def": "Biological membrane separating cell interior from exterior", "func": "Controls what enters and exits the cell", "imp": "Essential for cell survival"},
            {"term": "Cytoplasm", "type": "substance", "cat": "cell_biology", "def": "Gel-like substance inside cells excluding nucleus", "func": "Contains organelles and supports cellular processes", "imp": "Medium for cellular reactions"},
            
            # Biochemistry
            {"term": "Metabolism", "type": "process", "cat": "biochemistry", "def": "Set of life-sustaining chemical reactions in organisms", "func": "Converts food to energy and building blocks", "imp": "Essential for life"},
            {"term": "ATP", "type": "molecule", "cat": "biochemistry", "def": "Adenosine triphosphate - primary energy currency of cells", "func": "Stores and transfers chemical energy", "imp": "Essential for cellular energy"},
            {"term": "Hemoglobin", "type": "protein", "cat": "biochemistry", "def": "Protein in red blood cells that carries oxygen", "func": "Transports oxygen from lungs to tissues", "imp": "Essential for oxygen transport"},
            {"term": "Insulin", "type": "hormone", "cat": "biochemistry", "def": "Peptide hormone that regulates blood glucose", "func": "Controls glucose uptake by cells", "imp": "Critical for diabetes treatment"},
            {"term": "Glucose", "type": "molecule", "cat": "biochemistry", "def": "Simple sugar and primary source of energy", "func": "Primary fuel for cellular respiration", "imp": "Essential energy source"},
            {"term": "Lipid", "type": "molecule", "cat": "biochemistry", "def": "Group of naturally occurring molecules soluble in nonpolar solvents", "func": "Energy storage, cell membrane structure, signaling", "imp": "Essential for cell structure and energy"},
            
            # Microbiology
            {"term": "Bacteria", "type": "organism", "cat": "microbiology", "def": "Single-celled microorganisms without nucleus", "func": "Can be beneficial or pathogenic", "imp": "Important for health, industry, and disease"},
            {"term": "Virus", "type": "agent", "cat": "microbiology", "def": "Infectious agent that replicates inside living cells", "func": "Uses host cell machinery to reproduce", "imp": "Causes many diseases"},
            {"term": "Fungus", "type": "organism", "cat": "microbiology", "def": "Eukaryotic organisms including yeasts and molds", "func": "Decomposers, can be pathogenic or beneficial", "imp": "Important for ecology and industry"},
            {"term": "Pathogen", "type": "agent", "cat": "microbiology", "def": "Biological agent that causes disease", "func": "Infects and damages host organisms", "imp": "Target of medical treatments"},
            {"term": "Antibiotic", "type": "compound", "cat": "medicine", "def": "Substance that kills or inhibits bacteria", "func": "Treats bacterial infections", "imp": "Revolutionized medicine"},
            
            # Medicine
            {"term": "Disease", "type": "condition", "cat": "medicine", "def": "Abnormal condition affecting body function", "func": "Medical condition requiring treatment", "imp": "Target of healthcare"},
            {"term": "Cancer", "type": "disease", "cat": "medicine", "def": "Group of diseases involving uncontrolled cell growth", "func": "Can invade and spread to other body parts", "imp": "Major cause of death worldwide"},
            {"term": "Diabetes", "type": "disease", "cat": "medicine", "def": "Metabolic disease characterized by high blood sugar", "func": "Impaired insulin production or response", "imp": "Affects millions worldwide"},
            {"term": "Immune System", "type": "system", "cat": "medicine", "def": "Network of cells and proteins that defends the body", "func": "Protects against infection and disease", "imp": "Essential for survival"},
            {"term": "Antibody", "type": "protein", "cat": "medicine", "def": "Y-shaped protein produced by immune system", "func": "Recognizes and neutralizes pathogens", "imp": "Key to immunity and diagnostics"},
            {"term": "Vaccine", "type": "agent", "cat": "medicine", "def": "Substance that stimulates immune response", "func": "Provides immunity against diseases", "imp": "Prevents infectious diseases"},
            {"term": "Drug", "type": "compound", "cat": "medicine", "def": "Substance used to diagnose, cure, treat, or prevent disease", "func": "Modifies biological processes", "imp": "Fundamental to medical treatment"},
            
            # Biotechnology
            {"term": "CRISPR", "type": "technique", "cat": "biotechnology", "def": "Gene editing technology allowing precise DNA modification", "func": "Enables targeted genetic changes", "imp": "Revolutionary genetic engineering tool"},
            {"term": "PCR", "type": "technique", "cat": "biotechnology", "def": "Polymerase Chain Reaction - method to amplify DNA", "func": "Makes millions of copies of DNA", "imp": "Essential for molecular biology"},
            {"term": "Cloning", "type": "technique", "cat": "biotechnology", "def": "Process of creating identical copy of DNA or organism", "func": "Produces genetically identical copies", "imp": "Important for research and agriculture"},
            {"term": "Gene Therapy", "type": "technique", "cat": "biotechnology", "def": "Medical technique treating diseases by modifying genes", "func": "Corrects defective genes", "imp": "Potential cure for genetic diseases"},
            {"term": "Recombinant DNA", "type": "molecule", "cat": "biotechnology", "def": "DNA molecules formed by combining genetic material", "func": "Used to create genetically modified organisms", "imp": "Foundation of genetic engineering"},
            {"term": "Vector", "type": "molecule", "cat": "biotechnology", "def": "DNA vehicle used to transfer genetic material", "func": "Delivers DNA into cells", "imp": "Essential for gene therapy and cloning"},
            {"term": "Expression Vector", "type": "molecule", "cat": "biotechnology", "def": "Vector designed for gene expression in host cells", "func": "Produces proteins in host organisms", "imp": "Used for protein production"},
            {"term": "Bioreactor", "type": "device", "cat": "biotechnology", "def": "Device that supports biologically active environment", "func": "Used for cell growth and product formation", "imp": "Essential for bioprocessing"},
            
            # Bioinformatics
            {"term": "Sequencing", "type": "technique", "cat": "bioinformatics", "def": "Determining order of nucleotides in DNA or amino acids in proteins", "func": "Reads genetic code", "imp": "Essential for genomics"},
            {"term": "Genome Sequencing", "type": "technique", "cat": "bioinformatics", "def": "Determining complete DNA sequence of an organism", "func": "Reads entire genetic blueprint", "imp": "Foundation of genomics"},
            {"term": "Bioinformatics", "type": "field", "cat": "bioinformatics", "def": "Interdisciplinary field using computation to analyze biological data", "func": "Processes and analyzes biological data", "imp": "Essential for modern biology"},
            {"term": "Database", "type": "resource", "cat": "bioinformatics", "def": "Organized collection of biological data", "func": "Stores and retrieves biological information", "imp": "Essential for research"},
            {"term": "Alignment", "type": "technique", "cat": "bioinformatics", "def": "Arranging sequences to identify similarities", "func": "Compares genetic sequences", "imp": "Essential for evolutionary studies"},
            {"term": "Phylogeny", "type": "study", "cat": "bioinformatics", "def": "Study of evolutionary relationships among organisms", "func": "Shows evolutionary tree", "imp": "Understands evolution"},
            
            # Additional important terms
            {"term": "Apoptosis", "type": "process", "cat": "cell_biology", "def": "Programmed cell death", "func": "Removes damaged or unwanted cells", "imp": "Essential for development and cancer prevention"},
            {"term": "Metastasis", "type": "process", "cat": "medicine", "def": "Spread of cancer cells from original site", "func": "Cancer progression mechanism", "imp": "Major cause of cancer death"},
            {"term": "Stem Cell", "type": "cell", "cat": "cell_biology", "def": "Undifferentiated cell with self-renewal capacity", "func": "Can differentiate into various cell types", "imp": "Potential for regenerative medicine"},
            {"term": "Telomere", "type": "structure", "cat": "genetics", "def": "DNA-protein structures at chromosome ends", "func": "Protects chromosomes from degradation", "imp": "Related to aging and cancer"},
            {"term": "Oxidative Stress", "type": "process", "cat": "biochemistry", "def": "Cellular damage from reactive oxygen species", "func": "Associated with aging and disease", "imp": "Target of antioxidants"},
            {"term": "Signal Transduction", "type": "process", "cat": "molecular_biology", "def": "Process of cellular communication", "func": "Transmits signals from outside to inside cell", "imp": "Essential for cellular response"},
            {"term": "Transcription Factor", "type": "protein", "cat": "molecular_biology", "def": "Protein that controls gene transcription", "func": "Regulates gene expression", "imp": "Important for development and disease"},
            {"term": "Growth Factor", "type": "protein", "cat": "molecular_biology", "def": "Substance that stimulates cell growth", "func": "Promotes cell division and differentiation", "imp": "Important for development and cancer"},
            {"term": "Cytokine", "type": "molecule", "cat": "molecular_biology", "def": "Small protein used for cell signaling", "func": "Modulates immune response", "imp": "Important for immune system"},
            {"term": "Hormone", "type": "molecule", "cat": "biochemistry", "def": "Chemical messenger in multicellular organisms", "func": "Regulates physiological processes", "imp": "Essential for homeostasis"},
        ]
        
        for item in common_terms:
            self.add_term(
                term=item["term"],
                term_type=item["type"],
                category=item["cat"],
                definition=item["def"],
                function=item["func"],
                importance=item["imp"],
                source="Knowledge Base Core"
            )
        
        print(f"Added {len(common_terms)} common biological terms")
        return len(common_terms)
    
    def update_collection_status(self, source, terms_collected):
        """Update collection status"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        cursor.execute('''
            INSERT OR REPLACE INTO collection_status 
            (source, last_collection, terms_collected, status)
            VALUES (?, ?, ?, ?)
        ''', (source, datetime.now(), terms_collected, "completed"))
        
        conn.commit()
        conn.close()
    
    def collect_all(self):
        """Collect data from all sources"""
        print("\n" + "="*50)
        print("Starting Knowledge Base Collection")
        print("="*50)
        
        # Add common terms first
        self.add_common_biological_terms()
        
        # Collect from databases
        self.collect_ncbi_taxonomy(max_taxa=500)
        self.collect_ncbi_genes(search_terms=[
            "kinase", "receptor", "transporter", "enzyme", 
            "hormone", "growth factor", "immune", "neural",
            "cancer", "diabetes", "heart", "brain"
        ], max_per_term=50)
        
        self.collect_ncbi_proteins(search_terms=[
            "kinase", "receptor", "enzyme", "antibody",
            "hormone", "transport", "structural"
        ], max_per_term=30)
        
        self.collect_uniprot_proteins(max_results=1000)
        
        print("\n" + "="*50)
        print(f"Collection complete! Total terms: {self.get_total_terms()}")
        print("="*50)


# Run collection if executed directly
if __name__ == "__main__":
    collector = KnowledgeBaseCollector()
    collector.collect_all()
