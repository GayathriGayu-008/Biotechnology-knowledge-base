"""
Biotechnology Knowledge Base - Simplified Main Module
"""

import sqlite3
import json
import os
from datetime import datetime

class KnowledgeBase:
    """Main knowledge base class"""
    
    def __init__(self, db_path="biotech_kb.db"):
        self.db_path = db_path
        self.init_database()
        
    def init_database(self):
        """Initialize database"""
        conn = sqlite3.connect(self.db_path)
        c = conn.cursor()
        
        c.execute('''CREATE TABLE IF NOT EXISTS terms (
            id INTEGER PRIMARY KEY,
            term TEXT NOT NULL,
            term_type TEXT,
            category TEXT,
            definition TEXT,
            function TEXT,
            importance TEXT,
            source TEXT,
            created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )''')
        
        c.execute('CREATE INDEX IF NOT EXISTS idx_term ON terms(term)')
        c.execute('CREATE INDEX IF NOT EXISTS idx_cat ON terms(category)')
        
        conn.commit()
        conn.close()
        
    def add_term(self, term, term_type, category, definition, function, importance, source):
        """Add a term"""
        conn = sqlite3.connect(self.db_path)
        c = conn.cursor()
        c.execute('''INSERT INTO terms (term, term_type, category, definition, function, importance, source)
            VALUES (?, ?, ?, ?, ?, ?, ?)''',
            (term, term_type, category, definition, function, importance, source))
        conn.commit()
        conn.close()
        
    def search(self, query, limit=20):
        """Search terms"""
        conn = sqlite3.connect(self.db_path)
        c = conn.cursor()
        c.execute('''SELECT * FROM terms WHERE term LIKE ? OR definition LIKE ? LIMIT ?''',
            (f'%{query}%', f'%{query}%', limit))
        results = c.fetchall()
        conn.close()
        return results
    
    def get_all(self, category=None, limit=100):
        """Get all terms"""
        conn = sqlite3.connect(self.db_path)
        c = conn.cursor()
        if category:
            c.execute('SELECT * FROM terms WHERE category=? LIMIT ?', (category, limit))
        else:
            c.execute('SELECT * FROM terms LIMIT ?', (limit,))
        results = c.fetchall()
        conn.close()
        return results
    
    def get_count(self):
        """Get total count"""
        conn = sqlite3.connect(self.db_path)
        c = conn.cursor()
        c.execute('SELECT COUNT(*) FROM terms')
        count = c.fetchone()[0]
        conn.close()
        return count
    
    def populate(self):
        """Populate with comprehensive terms"""
        terms = [
            # Molecular Biology
            ("DNA", "molecule", "molecular_biology", "Deoxyribonucleic acid - hereditary material", "Carries genetic instructions", "Fundamental genetics", "Core"),
            ("RNA", "molecule", "molecular_biology", "Ribonucleic acid - single-stranded nucleic acid", "Gene expression", "Essential biology", "Core"),
            ("Protein", "molecule", "molecular_biology", "Large biomolecule of amino acids", "Catalysis, structure, transport", "Essential life", "Core"),
            ("Enzyme", "protein", "molecular_biology", "Biological catalyst", "Speeds chemical reactions", "Critical biochemistry", "Core"),
            ("Kinase", "enzyme", "molecular_biology", "Enzyme transferring phosphate groups", "Regulates signaling pathways", "Drug targets", "Core"),
            ("Receptor", "protein", "molecular_biology", "Protein binding molecules", "Signal transduction", "Drug targets", "Core"),
            ("Transporter", "protein", "molecular_biology", "Membrane protein moving molecules", "Nutrient uptake", "Cell function", "Core"),
            ("Hemoglobin", "protein", "molecular_biology", "Oxygen-carrying protein in blood", "Oxygen transport", "Blood function", "Core"),
            ("Insulin", "hormone", "molecular_biology", "Glucose-regulating hormone", "Controls blood sugar", "Diabetes treatment", "Core"),
            ("ATP", "molecule", "molecular_biology", "Cell energy currency", "Energy transfer", "Metabolism", "Core"),
            ("Gene", "unit", "genetics", "Unit of heredity", "Contains DNA instructions", "Inheritance", "Core"),
            ("Chromosome", "structure", "genetics", "DNA-protein structure", "Carries genes", "Genome", "Core"),
            ("Genome", "set", "genetics", "Complete genetic material", "All genes", "Heredity", "Core"),
            ("Mutation", "process", "genetics", "DNA sequence change", "Genetic variation", "Disease cause", "Core"),
            ("CRISPR", "technique", "biotechnology", "Gene editing technology", "Precise DNA modification", "Revolutionary biotech", "Core"),
            ("PCR", "technique", "biotechnology", "Polymerase Chain Reaction", "DNA amplification", "Diagnostics", "Core"),
            ("Cloning", "technique", "biotechnology", "Creating DNA copies", "Gene amplification", "Biotech", "Core"),
            ("Cell", "unit", "cell_biology", "Basic life unit", "All life functions", "Life", "Core"),
            ("Nucleus", "organelle", "cell_biology", "Cell's control center", "DNA storage", "Cell biology", "Core"),
            ("Mitochondria", "organelle", "cell_biology", "Energy producer", "ATP synthesis", "Powerhouse", "Core"),
            ("Ribosome", "organelle", "cell_biology", "Protein synthesizer", "Translation", "Protein factory", "Core"),
            ("Apoptosis", "process", "cell_biology", "Programmed cell death", "Cell removal", "Development", "Core"),
            ("Metabolism", "process", "biochemistry", "Chemical reactions in life", "Energy production", "Life processes", "Core"),
            ("Glycolysis", "pathway", "biochemistry", "Glucose breakdown", "Energy production", "Metabolism", "Core"),
            ("Hormone", "molecule", "biochemistry", "Chemical messenger", "Body regulation", "Homeostasis", "Core"),
            ("Virus", "agent", "microbiology", "Infectious particle", "Uses host cells", "Disease", "Core"),
            ("Bacteria", "organism", "microbiology", "Single-celled prokaryote", "Decomposition", "Industry/disease", "Core"),
            ("Antibody", "protein", "medicine", "Immune protein", "Antigen binding", "Immunity", "Core"),
            ("Vaccine", "agent", "medicine", "Immunity stimulator", "Disease prevention", "Prevention", "Core"),
            ("Cancer", "disease", "medicine", "Uncontrolled cell growth", "Can spread", "Leading death", "Core"),
            ("Diabetes", "disease", "medicine", "Blood sugar disease", "Insulin issues", "Metabolic disease", "Core"),
            ("Stem Cell", "cell", "cell_biology", "Undifferentiated cell", "Can become any cell", "Regenerative medicine", "Core"),
            ("Neuron", "cell", "cell_biology", "Nerve cell", "Signal transmission", "Nervous system", "Core"),
            ("Cytokine", "protein", "cell_biology", "Cell signaling protein", "Immune regulation", "Immunology", "Core"),
            ("Transcription Factor", "protein", "molecular_biology", "Gene expression regulator", "Controls genes", "Development", "Core"),
            ("Growth Factor", "protein", "molecular_biology", "Cell growth stimulator", "Cell division", "Development/cancer", "Core"),
            ("Phosphatase", "enzyme", "molecular_biology", "Removes phosphates", "Signal regulation", "Cell signaling", "Core"),
            ("Protease", "enzyme", "molecular_biology", "Protein breaker", "Digestion", "Drug targets", "Core"),
            ("Lipid", "molecule", "biochemistry", "Fat molecule", "Energy storage", "Membranes", "Core"),
            ("Glucose", "molecule", "biochemistry", "Simple sugar", "Energy source", "Metabolism", "Core"),
            ("Cholesterol", "molecule", "biochemistry", "Steroid alcohol", "Membrane component", "Heart disease", "Core"),
            ("Fungus", "organism", "microbiology", "Eukaryotic microbe", "Decomposition", "Industry/disease", "Core"),
            ("Pathogen", "agent", "microbiology", "Disease cause", "Infection", "Medicine", "Core"),
            ("Antibiotic", "compound", "medicine", "Bacteria killer", "Infection treatment", "Medicine", "Core"),
            ("Immune System", "system", "medicine", "Body defense", "Fights infection", "Health", "Core"),
            ("Hypertension", "condition", "medicine", "High blood pressure", "Heart risk", "Cardiovascular", "Core"),
            ("Alzheimer's", "disease", "medicine", "Neurodegenerative disease", "Memory loss", "Dementia", "Core"),
            ("Parkinson's", "disease", "medicine", "Movement disorder", "Dopamine loss", "Neurology", "Core"),
            ("HIV", "virus", "medicine", "Human Immunodeficiency Virus", "Causes AIDS", "Infectious disease", "Core"),
            ("Tumor", "mass", "medicine", "Abnormal cell growth", "Benign/malignant", "Oncology", "Core"),
            ("Metastasis", "process", "medicine", "Cancer spread", "Cancer progression", "Oncology", "Core"),
            ("Gene Therapy", "technique", "medicine", "Treating genes", "Genetic correction", "Future medicine", "Core"),
            ("Bioinformatics", "field", "bioinformatics", "Biological data analysis", "Computation for biology", "Modern biology", "Core"),
            ("Sequencing", "technique", "bioinformatics", "Reading DNA order", "Genome projects", "Genomics", "Core"),
            ("Microbiome", "community", "microbiology", "Microbe community", "Body health", "Medicine", "Core"),
            ("Biofilm", "structure", "microbiology", "Microbe community", "Infection", "Medicine", "Core"),
            ("Probiotic", "organism", "microbiology", "Beneficial bacteria", "Gut health", "Nutrition", "Core"),
            ("Autoimmune", "disease", "medicine", "Immune attacking body", "Various conditions", "Immunology", "Core"),
            ("Cardiovascular", "system", "medicine", "Heart and blood vessels", "Blood circulation", "Heart disease", "Core"),
            ("Pulmonology", "field", "medicine", "Lung study", "Respiratory disease", "Medicine", "Core"),
            ("Oncology", "field", "medicine", "Cancer study", "Cancer treatment", "Medicine", "Core"),
            ("Immunology", "field", "medicine", "Immune system study", "Immunity", "Medicine", "Core"),
            ("Neuroscience", "field", "medicine", "Nervous system study", "Brain function", "Medicine", "Core"),
            ("Genomics", "field", "genetics", "Genome study", "Complete DNA", "Modern biology", "Core"),
            ("Proteomics", "field", "molecular_biology", "Protein study", "All proteins", "Biology", "Core"),
            ("Metabolomics", "field", "biochemistry", "Metabolite study", "Metabolism", "Biochemistry", "Core"),
            ("Transcriptomics", "field", "genetics", "RNA study", "Gene expression", "Genomics", "Core"),
            ("Synthetic Biology", "field", "biotechnology", "Designing life", "New organisms", "Future biotech", "Core"),
            ("Bioreactor", "device", "biotechnology", "Biological reactor", "Cell culture", "Industry", "Core"),
            ("Fermentation", "process", "biotechnology", "Biological production", "Food/industry", "Industry", "Core"),
            ("Viral Vector", "tool", "biotechnology", "Virus delivery tool", "Gene therapy", "Medicine", "Core"),
            ("Monoclonal Antibody", "product", "biotechnology", "Single antibody type", "Cancer treatment", "Medicine", "Core"),
            ("Recombinant Protein", "product", "biotechnology", "Engineered protein", "Drugs", "Medicine", "Core"),
            ("Tissue Engineering", "field", "biotechnology", "Creating tissues", "Organ replacement", "Future medicine", "Core"),
            ("3D Bioprinting", "technique", "biotechnology", "Printing tissues", "Organ printing", "Future medicine", "Core"),
            ("Telomere", "structure", "genetics", "Chromosome end", "Aging", "Genetics", "Core"),
            ("Epigenetics", "field", "genetics", "Gene expression changes", "Non-DNA inheritance", "Modern genetics", "Core"),
            ("Cell Signaling", "process", "cell_biology", "Cell communication", "Coordination", "Biology", "Core"),
            ("Oxidative Stress", "process", "biochemistry", "Cell damage", "Aging", "Disease", "Core"),
            ("Apoptosis", "process", "cell_biology", "Cell death", "Development", "Cell biology", "Core"),
            ("Autophagy", "process", "cell_biology", "Self-digestion", "Recycling", "Cell biology", "Core"),
            ("Cytoplasm", "substance", "cell_biology", "Cell interior", "Reactions", "Cell biology", "Core"),
            ("Cell Membrane", "structure", "cell_biology", "Cell boundary", "Protection/transport", "Cell biology", "Core"),
            ("Golgi Apparatus", "organelle", "cell_biology", "Protein packaging", "Secretion", "Cell biology", "Core"),
            ("Endoplasmic Reticulum", "organelle", "cell_biology", "Protein synthesis", "Folding/modification", "Cell biology", "Core"),
            ("Lysosome", "organelle", "cell_biology", "Digestive organelle", "Cleanup", "Cell biology", "Core"),
            ("Chloroplast", "organelle", "cell_biology", "Plant energy organelle", "Photosynthesis", "Plants", "Core"),
            ("Peroxisome", "organelle", "cell_biology", "Oxidative organelle", "Detoxification", "Cell biology", "Core"),
            ("Centrosome", "organelle", "cell_biology", "Cell division organizer", "Mitosis", "Cell biology", "Core"),
            ("Fibroblast", "cell", "cell_biology", "Connective tissue cell", "Collagen production", "Healing", "Core"),
            ("Macrophage", "cell", "cell_biology", "Immune cell", "Phagocytosis", "Immunity", "Core"),
            ("T Cell", "cell", "cell_biology", "Immune cell", "Cell immunity", "Immunology", "Core"),
            ("B Cell", "cell", "cell_biology", "Antibody cell", "Humoral immunity", "Immunology", "Core"),
            ("Erythrocyte", "cell", "cell_biology", "Red blood cell", "Oxygen transport", "Blood", "Core"),
            ("Leukocyte", "cell", "cell_biology", "White blood cell", "Immunity", "Blood", "Core"),
            ("Meiosis", "process", "genetics", "Cell division for reproduction", "Gamete production", "Genetics", "Core"),
            ("Mitosis", "process", "genetics", "Cell division", "Growth/repair", "Cell biology", "Core"),
            ("Transcription", "process", "genetics", "DNA to RNA", "Gene expression", "Central dogma", "Core"),
            ("Translation", "process", "genetics", "RNA to protein", "Protein synthesis", "Central dogma", "Core"),
            ("Replication", "process", "genetics", "DNA copying", "Inheritance", "Cell division", "Core"),
            ("Promoter", "sequence", "genetics", "Gene start", "Transcription initiation", "Gene regulation", "Core"),
            ("Intron", "sequence", "genetics", "Non-coding DNA", "Splicing", "Gene structure", "Core"),
            ("Exon", "sequence", "genetics", "Coding DNA", "Protein coding", "Gene structure", "Core"),
            ("Plasmid", "molecule", "genetics", "Circular DNA", "Gene transfer", "Cloning", "Core"),
            ("SNP", "variant", "genetics", "Single base change", "Genetic variation", "Personalized medicine", "Core"),
            ("Vector", "tool", "biotechnology", "DNA vehicle", "Gene delivery", "Gene therapy", "Core"),
            ("RT-PCR", "technique", "biotechnology", "RNA to cDNA to PCR", "Gene expression", "Research", "Core"),
            ("Western Blot", "technique", "biotechnology", "Protein detection", "Protein analysis", "Research", "Core"),
            ("Flow Cytometry", "technique", "biotechnology", "Cell analysis", "Cell sorting", "Research", "Core"),
        ]
        
        for t in terms:
            self.add_term(*t)
        print(f"Added {len(terms)} terms")


if __name__ == "__main__":
    kb = KnowledgeBase()
    kb.populate()
    print(f"Total terms: {kb.get_count()}")
    print("Sample search for 'cancer':")
    for r in kb.search("cancer"):
        print(r)
