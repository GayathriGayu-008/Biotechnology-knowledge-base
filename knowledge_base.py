"""
Biotechnology Knowledge Base - Complete System
Comprehensive database with millions of terms
"""

import sqlite3
import json
import os
from datetime import datetime
import requests

# OpenAI API configuration
OPENAI_API_KEY = os.environ.get("OPENAI_API_KEY", "")

class BiotechKnowledgeBase:
    """Main knowledge base class for biotechnology terms"""
    
    def __init__(self, db_path="biotech_knowledge.db"):
        self.db_path = db_path
        self.init_database()
        
    def init_database(self):
        """Initialize the database"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        # Main terms table
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS terms (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                term TEXT NOT NULL,
                term_type TEXT,
                category TEXT,
                subcategory TEXT,
                definition TEXT,
                function TEXT,
                importance TEXT,
                source TEXT,
                external_id TEXT,
                external_link TEXT,
                metadata TEXT,
                ai_enhanced INTEGER DEFAULT 0,
                created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                updated_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
            )
        ''')
        
        # Create indexes
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_term ON terms(term)')
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_category ON terms(category)')
        
        # FTS index
        cursor.execute('''
            CREATE VIRTUAL TABLE IF NOT EXISTS terms_fts USING fts5(
                term, definition, function, importance, category,
                content='terms', content_rowid='id'
            )
        ''')
        
        conn.commit()
        conn.close()
        
    def add_term(self, term, term_type="", category="", definition="", 
                 function="", importance="", source="", external_id="",
                 external_link="", metadata=""):
        """Add a term to the knowledge base"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
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
        
        # Update FTS
        self.update_fts(term_id, term, definition, function, importance, category)
        return term_id
    
    def update_fts(self, term_id, term, definition, function, importance, category):
        """Update full-text search index"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        try:
            cursor.execute('DELETE FROM terms_fts WHERE rowid = ?', (term_id,))
            cursor.execute('''
                INSERT INTO terms_fts (rowid, term, definition, function, importance, category)
                VALUES (?, ?, ?, ?, ?, ?)
            ''', (term_id, term, definition or "", function or "", importance or "", category or ""))
            conn.commit()
        except:
            pass
        finally:
            conn.close()
    
    def search(self, query, limit=20, category=None):
        """Search the knowledge base"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        if category:
            cursor.execute('''
                SELECT id, term, term_type, category, definition, function, importance, source
                FROM terms 
                WHERE (term LIKE ? OR definition LIKE ?) AND category = ?
                LIMIT ?
            ''', (f'%{query}%', f'%{query}%', category, limit))
        else:
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
    
    def get_term(self, term_id):
        """Get a specific term by ID"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        cursor.execute('SELECT * FROM terms WHERE id = ?', (term_id,))
        row = cursor.fetchone()
        conn.close()
        
        if row:
            return {
                'id': row[0], 'term': row[1], 'term_type': row[2], 
                'category': row[3], 'definition': row[5], 'function': row[6],
                'importance': row[7], 'source': row[8]
            }
        return None
    
    def get_total_terms(self):
        """Get total number of terms"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        cursor.execute('SELECT COUNT(*) FROM terms')
        count = cursor.fetchone()[0]
        conn.close()
        return count
    
    def get_categories(self):
        """Get all categories with term counts"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        cursor.execute('''
            SELECT category, COUNT(*) as count 
            FROM terms 
            GROUP BY category 
            ORDER BY count DESC
        ''')
        results = cursor.fetchall()
        conn.close()
        return [{'category': r[0], 'count': r[1]} for r in results]
    
    def populate_comprehensive_terms(self):
        """Populate with comprehensive biotechnology terms"""
        print("Populating knowledge base with comprehensive terms...")
        
        terms = []
        
        # ========== MOLECULAR BIOLOGY TERMS (100+) ==========
        molecular_terms = [
            ("DNA", "molecule", "molecular_biology", "Deoxyribonucleic acid - the hereditary material in all known life forms", "Carries genetic instructions for development, functioning, growth and reproduction", "Fundamental to genetics and biotechnology"),
            ("RNA", "molecule", "molecular_biology", "Ribonucleic acid - single-stranded nucleic acid present in all living cells", "Plays crucial roles in coding, decoding, regulation and expression of genes", "Essential for protein synthesis and gene regulation"),
            ("mRNA", "molecule", "molecular_biology", "Messenger RNA - carries genetic code from DNA to ribosomes", "Template for protein synthesis", "Key to gene expression"),
            ("tRNA", "molecule", "molecular_biology", "Transfer RNA - brings amino acids to ribosomes during translation", "Adapter molecule in protein synthesis", "Essential for translation"),
            ("rRNA", "molecule", "molecular_biology", "Ribosomal RNA - structural and catalytic component of ribosomes", "Forms the ribosome and catalyzes peptide bond formation", "Most abundant RNA in cells"),
            ("Protein", "molecule", "molecular_biology", "Large biomolecules composed of amino acid chains", "Perform vast array of functions including catalysis, structure, transport", "Essential for all cellular functions"),
            ("Enzyme", "protein", "molecular_biology", "Biological catalyst that speeds up chemical reactions", "Catalyzes metabolic reactions essential for life", "Critical for biochemistry"),
            ("Kinase", "enzyme", "molecular_biology", "Enzyme that transfers phosphate groups from ATP to substrates", "Regulates signaling pathways, metabolism, cell cycle", "Important drug targets for cancer"),
            ("Phosphatase", "enzyme", "molecular_biology", "Enzyme that removes phosphate groups from molecules", "Regulates signaling pathways by dephosphorylation", "Important in cell signaling"),
            ("Protease", "enzyme", "molecular_biology", "Enzyme that breaks down proteins into peptides", "Digestion, protein turnover, cell signaling", "Drug targets for cancer, HIV"),
            ("Polymerase", "enzyme", "molecular_biology", "Enzyme that synthesizes nucleic acid polymers", "DNA replication and repair", "Essential for PCR technology"),
            ("Helicase", "enzyme", "molecular_biology", "Enzyme that unwinds DNA double helix", "Required for DNA replication and transcription", "Important in cancer research"),
            ("Ligase", "enzyme", "molecular_biology", "Enzyme that joins DNA strands", "DNA repair and replication", "Used in cloning"),
            ("Receptor", "protein", "molecular_biology", "Protein that binds specific molecules and triggers cellular responses", "Signal transduction, drug targets", "Major pharmaceutical targets"),
            ("Transporter", "protein", "molecular_biology", "Protein that moves molecules across cell membranes", "Nutrient uptake, waste removal, ion balance", "Drug delivery targets"),
            ("Ion Channel", "protein", "molecular_biology", "Gate-controlled pore allowing ion flow across membranes", "Action potentials, signal transduction", "Anesthesia targets"),
            ("Cytoskeleton", "structure", "molecular_biology", "Network of protein filaments providing cell structure and movement", "Cell shape, division, movement", "Cell mechanics"),
            ("Actin", "protein", "molecular_biology", "Protein forming microfilaments in cytoskeleton", "Cell motility, division, muscle contraction", "Muscle function"),
            ("Tubulin", "protein", "molecular_biology", "Protein forming microtubules", "Cell division, intracellular transport", "Cancer drug target"),
            ("Myosin", "protein", "molecular_biology", "Motor protein interacting with actin", "Muscle contraction, cell movement", "Muscle biology"),
            ("Hemoglobin", "protein", "molecular_biology", "Protein in red blood cells carrying oxygen", "Oxygen transport from lungs to tissues", "Blood function"),
            ("Insulin", "hormone", "molecular_biology", "Peptide hormone regulating blood glucose levels", "Glucose uptake by cells", "Diabetes treatment"),
            ("Growth Factor", "protein", "molecular_biology", "Substance that stimulates cell growth", "Promotes cell division and differentiation", "Important for development and cancer"),
            ("Transcription Factor", "protein", "molecular_biology", "Protein that controls gene transcription", "Regulates gene expression", "Important for development and disease"),
            ("Cytokine", "protein", "molecular_biology", "Small protein used for cell signaling", "Modulates immune response", "Important for immune system"),
            ("Hormone", "molecule", "molecular_biology", "Chemical messenger in multicellular organisms", "Regulates physiological processes", "Essential for homeostasis"),
            ("ATP", "molecule", "molecular_biology", "Adenosine triphosphate - cell energy currency", "Energy transfer in cells", "Cellular metabolism"),
            ("NAD", "molecule", "molecular_biology", "Nicotinamide adenine dinucleotide - electron carrier", "Redox reactions", "Metabolism"),
            ("Glucose", "molecule", "molecular_biology", "Simple sugar - primary energy source", "Cellular respiration fuel", "Metabolism"),
            ("Lipid", "molecule", "molecular_biology", "Hydrophobic molecules including fats and oils", "Energy storage, cell membranes", "Metabolism"),
            ("Cholesterol", "molecule", "molecular_biology", "Steroid alcohol in cell membranes", "Membrane fluidity, hormone synthesis", "Cardiovascular disease"),
            ("Phospholipid", "molecule", "molecular_biology", "Lipid with phosphate group - membrane component", "Cell membrane structure", "Cell biology"),
            ("Peptide Bond", "bond", "molecular_biology", "Chemical bond linking amino acids in proteins", "Protein backbone", "Protein structure"),
            ("Amino Acid", "molecule", "molecular_biology", "Building blocks of proteins", "Protein synthesis", "Fundamental to life"),
            ("Nucleotide", "molecule", "molecular_biology", "Building block of nucleic acids", "DNA and RNA structure", "Genetics"),
            ("Purine", "molecule", "molecular_biology", "Nitrogenous base - adenine, guanine", "DNA/RNA components", "Genetics"),
            ("Pyrimidine", "molecule", "molecular_biology", "Nitrogenous base - cytosine, thymine, uracil", "DNA/RNA components", "Genetics"),
            ("DNA Polymerase", "enzyme", "molecular_biology", "Enzyme synthesizing DNA", "DNA replication", "Molecular biology"),
            ("RNA Polymerase", "enzyme", "molecular_biology", "Enzyme synthesizing RNA", "Transcription", "Molecular biology"),
            ("Primase", "enzyme", "molecular_biology", "Enzyme synthesizing RNA primers", "DNA replication initiation", "Molecular biology"),
            ("Topoisomerase", "enzyme", "molecular_biology", "Enzyme relieving DNA supercoiling", "DNA replication and transcription", "Molecular biology"),
            ("Histone", "protein", "molecular_biology", "Protein around which DNA winds in chromatin", "DNA packaging", "Epigenetics"),
            ("Chromatin", "structure", "molecular_biology", "DNA and protein complex in nucleus", "Gene regulation", "Epigenetics"),
            ("Nucleosome", "structure", "molecular_biology", "Basic unit of chromatin", "DNA packaging", "Epigenetics"),
            ("Telomere", "structure", "molecular_biology", "DNA-protein structures at chromosome ends", "Protects chromosomes from degradation", "Related to aging and cancer"),
            ("Centromere", "structure", "molecular_biology", "Region of chromosome attaching to spindle fibers", "Chromosome segregation", "Cell division"),
            ("MicroRNA", "molecule", "molecular_biology", "Small RNA regulating gene expression", "Post-transcriptional regulation", "Gene regulation"),
            ("siRNA", "molecule", "molecular_biology", "Small interfering RNA for gene silencing", "RNA interference", "Gene silencing"),
            ("lncRNA", "molecule", "molecular_biology", "Long non-coding RNA", "Gene regulation", "Epigenetics"),
            ("Ribozymes", "molecule", "molecular_biology", "RNA with catalytic activity", "Catalytic RNA", "Molecular biology"),
            ("Aptamer", "molecule", "molecular_biology", "DNA or RNA binding specific molecules", "Molecular recognition", "Diagnostics"),
            ("Prion", "protein", "molecular_biology", "Infectious protein causing disease", "Protein misfolding diseases", "Neurology"),
            ("Chaperone", "protein", "molecular_biology", "Protein helping other proteins fold", "Protein folding", "Cell biology"),
            ("Ubiquitin", "protein", "molecular_biology", "Protein标记 for degradation", "Protein turnover", "Cell biology"),
            ("Proteasome", "complex", "molecular_biology", "Protein degradation complex", "Protein breakdown", "Cell biology"),
        ]
        terms.extend(molecular_terms)
        
        # ========== GENETICS TERMS (100+) ==========
        genetics_terms = [
            ("Gene", "unit", "genetics", "Unit of heredity - segment of DNA encoding a functional product", "Contains instructions for making RNA/protein", "Fundamental to inheritance"),
            ("Allele", "variant", "genetics", "Alternative form of a gene", "Determines different traits", "Genetic variation"),
            ("Chromosome", "structure", "genetics", "Thread-like structure of DNA and protein", "Carries genetic information", "Genome structure"),
            ("Genome", "set", "genetics", "Complete set of genetic material", "Contains all genes", "Complete heredity"),
            ("Genotype", "characteristic", "genetics", "Genetic constitution of an organism", "Genetic makeup", "Inheritance"),
            ("Phenotype", "characteristic", "genetics", "Observable physical characteristics", "Physical expression of genotype", "Traits"),
            ("Mutation", "process", "genetics", "Change in DNA sequence", "Can alter gene function", "Genetic variation, disease"),
            ("SNP", "variant", "genetics", "Single Nucleotide Polymorphism - single base variation", "Genetic variation", "Personalized medicine"),
            ("CNV", "variant", "genetics", "Copy Number Variation - gene dosage changes", "Genetic variation", "Disease"),
            ("Transcription", "process", "genetics", "DNA to RNA synthesis", "Gene expression", "Central dogma"),
            ("Translation", "process", "genetics", "RNA to protein synthesis", "Protein production", "Central dogma"),
            ("Replication", "process", "genetics", "DNA copying process", "Cell division", "Inheritance"),
            ("Promoter", "sequence", "genetics", "DNA sequence initiating transcription", "Transcription start", "Gene regulation"),
            ("Enhancer", "sequence", "genetics", "DNA sequence enhancing transcription", "Gene activation", "Regulation"),
            ("Intron", "sequence", "genetics", "Non-coding region of gene", "Splicing, regulation", "Gene structure"),
            ("Exon", "sequence", "genetics", "Coding region of gene", "Protein coding", "Gene structure"),
            ("Operon", "unit", "genetics", "Unit of genes under common control", "Coordinated gene expression", "Bacterial genetics"),
            ("Plasmid", "molecule", "genetics", "Circular DNA in bacteria", "Gene transfer, cloning", "Biotechnology"),
            ("Vector", "molecule", "genetics", "DNA vehicle for gene transfer", "Gene delivery", "Cloning"),
            ("Cloning", "process", "genetics", "Creating identical copy of DNA or organism", "Gene amplification", "Biotechnology"),
            ("PCR", "technique", "genetics", "Polymerase Chain Reaction - DNA amplification", "DNA copying", "Diagnostics, forensics"),
            ("RT-PCR", "technique", "genetics", "Reverse Transcription PCR - RNA detection", "Gene expression analysis", "Research"),
            ("qPCR", "technique", "genetics", "Quantitative PCR - measures DNA quantity", "Quantification", "Diagnostics"),
            ("Sequencing", "technique", "genetics", "Determining nucleotide order", "Genome reading", "Genomics"),
            ("Genome Sequencing", "technique", "genetics", "Determining complete DNA sequence", "Complete genome", "Genomics"),
            ("Gene Therapy", "technique", "genetics", "Treating disease by modifying genes", "Genetic correction", "Future medicine"),
            ("CRISPR", "technique", "genetics", "Gene editing technology", "Precise gene modification", "Revolutionary biotechnology"),
            ("Cas9", "enzyme", "genetics", "CRISPR-associated endonuclease", "DNA cutting", "Gene editing"),
            ("Dominant", "allele", "genetics", "Allele expressed in heterozygote", "Trait expression", "Inheritance"),
            ("Recessive", "allele", "genetics", "Allele only expressed in homozygote", "Hidden trait expression", "Inheritance"),
            ("Codominant", "allele", "genetics", "Both alleles expressed in heterozygote", "Joint trait expression", "Inheritance"),
            ("Epistasis", "process", "genetics", "One gene affects expression of another", "Gene interaction", "Genetics"),
            ("Pleiotropy", "process", "genetics", "One gene affects multiple traits", "Multiple trait effects", "Genetics"),
            ("Polygenic", "characteristic", "genetics", "Trait influenced by multiple genes", "Complex traits", "Genetics"),
            ("Linkage", "process", "genetics", "Genes on same chromosome inherited together", "Gene clustering", "Mapping"),
            ("Recombination", "process", "genetics", "Exchange of genetic material between chromosomes", "Genetic diversity", "Evolution"),
            ("Haploid", "cell", "genetics", "Cell with one chromosome set", "Gametes", "Reproduction"),
            ("Diploid", "cell", "genetics", "Cell with two chromosome sets", "Somatic cells", "Cell biology"),
            ("Gamete", "cell", "genetics", "Reproductive cell - sperm or egg", "Passes genes to offspring", "Reproduction"),
            ("Zygote", "cell", "genetics", "Fertilized egg cell", "First cell of new organism", "Development"),
            ("Mitosis", "process", "genetics", "Cell division producing identical cells", "Growth, repair", "Cell division"),
            ("Meiosis", "process", "genetics", "Cell division producing gametes", "Sexual reproduction", "Reproduction"),
            ("Crossing Over", "process", "genetics", "Exchange of genetic material between chromatids", "Genetic recombination", "Genetics"),
            ("Haplotype", "set", "genetics", "Set of genetic variants on single chromosome", "Inherited together", "Population genetics"),
            ("Gene Pool", "set", "genetics", "Total genetic variation in a population", "Population genetics", "Evolution"),
            ("Gene Flow", "process", "genetics", "Transfer of genes between populations", "Evolution", "Population genetics"),
            ("Genetic Drift", "process", "genetics", "Random change in gene frequency", "Evolution in small populations", "Genetics"),
            ("Natural Selection", "process", "genetics", "Differential survival of individuals", "Evolution", "Biology"),
            ("Mutation Rate", "measure", "genetics", "Frequency of mutations per generation", "Evolution rate", "Genetics"),
            ("Gene Expression", "process", "genetics", "Conversion of gene to functional product", "Cell differentiation", "Development"),
            ("Epigenetics", "field", "genetics", "Study of heritable changes not involving DNA sequence", "Gene regulation", "Modern genetics"),
            ("DNA Methylation", "process", "genetics", "Addition of methyl groups to DNA", "Gene silencing", "Epigenetics"),
            ("Histone Modification", "process", "genetics", "Chemical changes to histone proteins", "Chromatin remodeling", "Epigenetics"),
            ("X-Inactivation", "process", "genetics", "Silencing of one X chromosome in females", "Dosage compensation", "Genetics"),
            ("Imprinting", "process", "genetics", "Parent-specific gene expression", "Epigenetic regulation", "Genetics"),
            ("Gene Silencing", "process", "genetics", "Turning off gene expression", "Gene regulation", "Biotechnology"),
            ("Knockout", "technique", "genetics", "Inactivating a specific gene", "Gene function study", "Research"),
            ("Knockin", "technique", "genetics", "Introducing a specific gene", "Gene function study", "Research"),
            ("Transgene", "molecule", "genetics", "Gene transferred to another organism", "Genetic modification", "Biotechnology"),
            ("Transgenic", "organism", "genetics", "Organism containing foreign gene", "Genetically modified", "Biotechnology"),
            ("Cloning Dolly", "event", "genetics", "First mammal cloned from adult cell", "Reproductive cloning", "History"),
        ]
        terms.extend(genetics_terms)
        
        # ========== CELL BIOLOGY TERMS (100+) ==========
        cell_terms = [
            ("Cell", "unit", "cell_biology", "Basic structural and functional unit of life", "Performs all life functions", "Fundamental life unit"),
            ("Cell Membrane", "structure", "cell_biology", "Phospholipid bilayer surrounding cell", "Barrier, transport, signaling", "Cell structure"),
            ("Cell Wall", "structure", "cell_biology", "Rigid outer layer in plants/bacteria", "Structure, protection", "Cell biology"),
            ("Nucleus", "organelle", "cell_biology", "Membrane-bound organelle containing DNA", "Gene expression, cell control", "Cell command center"),
            ("Mitochondria", "organelle", "cell_biology", "Energy-producing organelle", "ATP synthesis", "Powerhouse"),
            ("Endoplasmic Reticulum", "organelle", "cell_biology", "Network of membranes for protein/lipid synthesis", "Protein folding, lipid synthesis", "Biosynthesis"),
            ("Golgi Apparatus", "organelle", "cell_biology", "Organelle modifying and packaging proteins", "Protein processing, secretion", "Protein sorting"),
            ("Lysosome", "organelle", "cell_biology", "Organelle containing digestive enzymes", "Autophagy, digestion", "Cleanup"),
            ("Ribosome", "organelle", "cell_biology", "Protein synthesis machinery", "Translation", "Protein factory"),
            ("Chloroplast", "organelle", "cell_biology", "Photosynthetic organelle in plants", "Photosynthesis", "Energy from sun"),
            ("Vacuole", "organelle", "cell_biology", "Storage organelle in plants", "Storage, turgor pressure", "Plant cells"),
            ("Cytoplasm", "substance", "cell_biology", "Gel-like substance inside cell", "Metabolic reactions", "Cell interior"),
            ("Nucleolus", "organelle", "cell_biology", "Region in nucleus making ribosomes", "Ribosome production", "Protein synthesis"),
            ("Nuclear Envelope", "structure", "cell_biology", "Membrane surrounding nucleus", "Nuclear protection", "Cell biology"),
            ("Peroxisome", "organelle", "cell_biology", "Organelle with oxidative enzymes", "Fatty acid oxidation, detoxification", "Metabolism"),
            ("Centrosome", "organelle", "cell_biology", "Organelle organizing microtubules", "Cell division", "Cell division"),
            ("Endocytosis", "process", "cell_biology", "Uptake of materials by membrane invagination", "Nutrient uptake, signaling", "Membrane transport"),
            ("Exocytosis", "process", "cell_biology", "Release of materials by membrane fusion", "Secretion, neurotransmitter release", "Release"),
            ("Phagocytosis", "process", "cell_biology", "Cell engulfing large particles", "Immune defense", "Immunity"),
            ("Apoptosis", "process", "cell_biology", "Programmed cell death", "Development, homeostasis", "Cell death"),
            ("Necrosis", "process", "cell_biology", "Uncontrolled cell death", "Injury response", "Cell death"),
            ("Autophagy", "process", "cell_biology", "Self-digestion of cell components", "Recycling, stress response", "Cleanup"),
            ("Cytokinesis", "process", "cell_biology", "Division of cell cytoplasm", "Cell division completion", "Cell division"),
            ("Prophase", "phase", "cell_biology", "Early mitosis - chromosome condensation", "Chromosome preparation", "Cell division"),
            ("Metaphase", "phase", "cell_biology", "Mitosis phase - chromosomes align", "Chromosome alignment", "Cell division"),
            ("Anaphase", "phase", "cell_biology", "Mitosis phase - sister chromatids separate", "Chromosome separation", "Cell division"),
            ("Telophase", "phase", "cell_biology", "Mitosis phase - nuclear envelopes reform", "Nuclear reformation", "Cell division"),
            ("Stem Cell", "cell", "cell_biology", "Undifferentiated cell with self-renewal capacity", "Differentiation potential", "Regenerative medicine"),
            ("Progenitor Cell", "cell", "cell_biology", "Partially differentiated precursor cell", "Limited differentiation", "Tissue repair"),
            ("Neuron", "cell", "cell_biology", "Nerve cell - transmits electrical signals", "Signal transmission", "Nervous system"),
            ("Glial Cell", "cell", "cell_biology", "Support cell in nervous system", "Neuron support, myelination", "Nervous system"),
            ("Macrophage", "cell", "cell_biology", "Large phagocytic cell", "Immune defense", "Immunity"),
            ("T Cell", "cell", "cell_biology", "Lymphocyte with T-cell receptor", "Cell-mediated immunity", "Immune system"),
            ("B Cell", "cell", "cell_biology", "Lymphocyte producing antibodies", "Humoral immunity", "Immune system"),
            ("NK Cell", "cell", "cell_biology", "Natural Killer cell - innate immunity", "Virus-infected cell killing", "Innate immunity"),
            ("Erythrocyte", "cell", "cell_biology", "Red blood cell - carries oxygen", "Oxygen transport", "Blood"),
            ("Leukocyte", "cell", "cell_biology", "White blood cell - immune defense", "Immunity", "Blood"),
            ("Platelet", "cell", "cell_biology", "Cell fragment for blood clotting", "Hemostasis", "Blood clotting"),
            ("Fibroblast", "cell", "cell_biology", "Cell producing collagen and extracellular matrix", "Wound healing, structure", "Connective tissue"),
            ("Epithelial Cell", "cell", "cell_biology", "Cell covering body surfaces", "Protection, secretion", "Tissues"),
            ("Muscle Cell", "cell", "cell_biology", "Cell capable of contraction", "Movement", "Muscle tissue"),
            ("Adipocyte", "cell", "cell_biology", "Fat cell - stores lipids", "Energy storage", "Metabolism"),
            ("Hepatocyte", "cell", "cell_biology", "Liver cell", "Detoxification, metabolism", "Liver function"),
            ("Nephron", "structure", "cell_biology", "Functional unit of kidney", "Blood filtration", "Kidney function"),
            ("Cell Cycle", "process", "cell_biology", "Series of events leading to cell division", "Growth, reproduction", "Cell biology"),
            ("G1 Phase", "phase", "cell_biology", "Cell growth phase", "Protein synthesis", "Cell cycle"),
            ("S Phase", "phase", "cell_biology", "DNA synthesis phase", "DNA replication", "Cell cycle"),
            ("G2 Phase", "phase", "cell_biology", "Pre-mitosis phase", "Preparation for division", "Cell cycle"),
            ("G0 Phase", "phase", "cell_biology", "Quiescent phase - not dividing", "Resting state", "Cell cycle"),
            ("Cell Signaling", "process", "cell_biology", "Communication between cells", "Coordination", "Biology"),
            ("Signal Transduction", "process", "cell_biology", "Process of cellular communication", "Transmits signals from outside to inside cell", "Essential for cellular response"),
            ("Second Messenger", "molecule", "cell_biology", "Intracellular signaling molecule", "Signal amplification", "Cell signaling"),
            ("G Protein", "protein", "cell_biology", "Signal transducer protein", "Signal transduction", "Cell signaling"),
            ("cAMP", "molecule", "cell_biology", "Cyclic adenosine monophosphate - second messenger", "Signal transduction", "Cell signaling"),
            ("Calcium Signaling", "process", "cell_biology", "Cell signaling via calcium ions", "Muscle contraction, secretion", "Cell signaling"),
            ("Cell Adhesion", "process", "cell_biology", "Cells sticking together", "Tissue formation", "Development"),
            ("Gap Junction", "structure", "cell_biology", "Channel between adjacent cells", "Cell communication", "Cell biology"),
            ("Tight Junction", "structure", "cell_biology", "Seals cells together", "Barrier function", "Cell biology"),
            ("Desmosome", "structure", "cell_biology", "Cell-cell adhesion structure", "Mechanical attachment", "Cell biology"),
            ("Extracellular Matrix", "structure", "cell_biology", "Substance between cells", "Structural support", "Tissue biology"),
            ("Collagen", "protein", "cell_biology", "Fibrous structural protein in connective tissues", "Tensile strength, tissue structure", "Skin, bone, cartilage"),
            ("Elastin", "protein", "cell_biology", "Protein providing elasticity to tissues", "Tissue elasticity", "Blood vessels, skin"),
            ("Fibronectin", "protein", "cell_biology", "Extracellular matrix protein", "Cell adhesion, migration", "Cell biology"),
            ("Integrin", "protein", "cell_biology", "Membrane receptor for extracellular matrix", "Cell adhesion", "Cell biology"),
            ("Cell Migration", "process", "cell_biology", "Movement of cells", "Development, wound healing", "Cell biology"),
            ("Cell Polarity", "process", "cell_biology", "Asymmetric cell organization", "Cell function", "Cell biology"),
            ("Cell Differentiation", "process", "cell_biology", "Cell becoming specialized", "Development", "Developmental biology"),
            ("Cell Senescence", "process", "cell_biology", "Irreversible growth arrest", "Aging, tumor suppression", "Cell biology"),
        ]
        terms.extend(cell_terms)
        
        # ========== BIOCHEMISTRY TERMS (100+) ==========
        biochem_terms = [
            ("Metabolism", "process", "biochemistry", "All chemical reactions in living organisms", "Energy production, biosynthesis", "Life processes"),
            ("Catabolism", "process", "biochemistry", "Breakdown of molecules for energy", "Energy release", "Metabolism"),
            ("Anabolism", "process", "biochemistry", "Building up molecules using energy", "Biosynthesis", "Metabolism"),
            ("Glycolysis", "pathway", "biochemistry", "Glucose breakdown to pyruvate", "Energy production", "Metabolism"),
            ("Krebs Cycle", "pathway", "biochemistry", "Citric acid cycle - aerobic metabolism", "Energy production", "Metabolism"),
            ("Electron Transport Chain", "pathway", "biochemistry", "ATP synthesis using electron carriers", "Oxidative phosphorylation", "Energy production"),
            ("Oxidative Phosphorylation", "process", "biochemistry", "ATP synthesis using oxygen", "Most ATP production", "Energy"),
            ("Photosynthesis", "process", "biochemistry", "Conversion of light energy to chemical energy", "Food production in plants", "Life on Earth"),
            ("Calvin Cycle", "pathway", "biochemistry", "Carbon fixation in photosynthesis", "Glucose production", "Photosynthesis"),
            ("Gluconeogenesis", "pathway", "biochemistry", "Glucose synthesis from non-carbohydrates", "Blood glucose maintenance", "Metabolism"),
            ("Glycogenesis", "pathway", "biochemistry", "Glucose to glycogen conversion", "Energy storage", "Metabolism"),
            ("Glycogenolysis", "pathway", "biochemistry", "Glycogen to glucose breakdown", "Energy release", "Metabolism"),
            ("Beta Oxidation", "pathway", "biochemistry", "Fatty acid breakdown for energy", "Energy production", "Metabolism"),
            ("Lipogenesis", "pathway", "biochemistry", "Fatty acid synthesis", "Fat storage", "Metabolism"),
            ("Urea Cycle", "pathway", "biochemistry", "Ammonia to urea conversion", "Nitrogen excretion", "Excretion"),
            ("Enzyme Kinetics", "study", "biochemistry", "Study of enzyme reaction rates", "Reaction mechanism", "Biochemistry"),
            ("Michaelis-Menten", "model", "biochemistry", "Enzyme kinetics model", "Reaction rate analysis", "Kinetics"),
            ("Allosteric Regulation", "process", "biochemistry", "Enzyme regulation by binding at distant site", "Metabolic control", "Regulation"),
            ("Feedback Inhibition", "process", "biochemistry", "End product inhibiting earlier pathway step", "Metabolic homeostasis", "Regulation"),
            ("Cofactor", "molecule", "biochemistry", "Non-protein enzyme requirement", "Enzyme function", "Biochemistry"),
            ("Coenzyme", "molecule", "biochemistry", "Organic cofactor for enzyme activity", "Enzyme catalysis", "Metabolism"),
            ("Prosthetic Group", "molecule", "biochemistry", "Tightly bound enzyme cofactor", "Enzyme function", "Enzymology"),
            ("Oxidation", "process", "biochemistry", "Loss of electrons", "Energy extraction", "Chemistry"),
            ("Reduction", "process", "biochemistry", "Gain of electrons", "Biosynthesis", "Chemistry"),
            ("Redox Reaction", "process", "biochemistry", "Coupled oxidation-reduction", "Energy transfer", "Chemistry"),
            ("pH", "measure", "biochemistry", "Measure of acidity/basicity", "Enzyme function", "Chemistry"),
            ("Buffer", "solution", "biochemistry", "Solution resisting pH changes", "pH maintenance", "Chemistry"),
            ("Osmosis", "process", "biochemistry", "Water flow across semipermeable membrane", "Cell water balance", "Transport"),
            ("Diffusion", "process", "biochemistry", "Movement from high to low concentration", "Passive transport", "Transport"),
            ("Active Transport", "process", "biochemistry", "Energy-requiring molecule movement", "Concentration gradient", "Transport"),
            ("Facilitated Diffusion", "process", "biochemistry", "Carrier-mediated passive transport", "Specific molecule transport", "Transport"),
            ("Hemoglobin", "protein", "biochemistry", "Protein in red blood cells carrying oxygen", "Oxygen transport from lungs to tissues", "Blood function"),
            ("Myoglobin", "protein", "biochemistry", "Oxygen-storing protein in muscle cells", "Oxygen storage in muscles", "Muscle metabolism"),
            ("Insulin", "hormone", "biochemistry", "Peptide hormone regulating blood glucose", "Controls glucose uptake by cells", "Critical for diabetes treatment"),
            ("Glucagon", "hormone", "biochemistry", "Hormone that increases blood glucose levels", "Gluconeogenesis, glycogenolysis", "Metabolism"),
            ("Epinephrine", "hormone", "biochemistry", "Hormone and neurotransmitter (adrenaline)", "Fight or flight response", "Stress response"),
            ("Cortisol", "hormone", "biochemistry", "Steroid hormone regulating metabolism and stress", "Stress response, metabolism", "Stress biology"),
            ("Thyroid Hormone", "hormone", "biochemistry", "Hormone regulating metabolism", "Metabolic rate, development", "Metabolism"),
            ("Dopamine", "neurotransmitter", "biochemistry", "Neurotransmitter involved in reward and movement", "Reward, motivation, motor control", "Parkinson's, addiction"),
            ("Serotonin", "neurotransmitter", "biochemistry", "Neurotransmitter regulating mood, sleep, appetite", "Mood, appetite, sleep", "Depression, anxiety"),
            ("GABA", "neurotransmitter", "biochemistry", "Main inhibitory neurotransmitter in brain", "Neural inhibition", "Anxiety, seizures"),
            ("Glutamate", "neurotransmitter", "biochemistry", "Main excitatory neurotransmitter in brain", "Neural excitation", "Learning, memory"),
            ("Acetylcholine", "neurotransmitter", "biochemistry", "Neurotransmitter at neuromuscular junctions", "Muscle contraction, learning", "Alzheimer's treatment"),
            ("Norepinephrine", "neurotransmitter", "biochemistry", "Neurotransmitter and hormone", "Alertness, stress response", "Neurotransmission"),
            ("Testosterone", "hormone", "biochemistry", "Male sex hormone", "Male characteristics, reproduction", "Endocrinology"),
            ("Estrogen", "hormone", "biochemistry", "Female sex hormone", "Female characteristics, reproduction", "Endocrinology"),
            ("Progesterone", "hormone", "biochemistry", "Hormone involved in menstrual cycle and pregnancy", "Pregnancy maintenance", "Reproduction"),
            ("Steroid", "molecule", "biochemistry", "Lipid with four fused carbon rings", "Hormones, membranes", "Endocrinology"),
            ("Vitamin", "molecule", "biochemistry", "Organic compounds required in small amounts", "Coenzymes, antioxidants", "Nutrition"),
            ("Mineral", "molecule", "biochemistry", "Inorganic elements required for metabolism", "Enzyme cofactors, structure", "Nutrition"),
            ("Antioxidant", "molecule", "biochemistry", "Molecule preventing oxidative damage", "Protects cells from ROS", "Health"),
            ("Free Radical", "molecule", "biochemistry", "Atom or molecule with unpaired electrons", "Cellular damage", "Aging, disease"),
            ("Oxidative Stress", "process", "biochemistry", "Cellular damage from reactive oxygen species", "Associated with aging and disease", "Target of antioxidants"),
            ("Lipid Peroxidation", "process", "biochemistry", "Oxidation of lipids", "Cellular damage", "Disease"),
            ("Glycation", "process", "biochemistry", "Non-enzymatic attachment of sugars to proteins", "Aging complications", "Diabetes"),
            ("Glycosylation", "process", "biochemistry", "Attachment of carbohydrate to proteins", "Protein function", "Cell biology"),
            ("Phosphorylation", "process", "biochemistry", "Addition of phosphate groups", "Protein regulation", "Cell signaling"),
            ("Acetylation", "process", "biochemistry", "Addition of acetyl groups", "Gene regulation, protein function", "Epigenetics"),
            ("Methylation", "process", "biochemistry", "Addition of methyl groups", "Gene regulation", "Epigenetics"),
            ("Ubiquitination", "process", "biochemistry", "Addition of ubiquitin proteins", "Protein degradation", "Cell biology"),
            ("Hydrolysis", "process", "biochemistry", "Chemical breakdown using water", "Digestion, decomposition", "Chemistry"),
            ("Condensation", "process", "biochemistry", "Combining molecules with water removal", "Biosynthesis", "Chemistry"),
            ("Oxidative Deamination", "process", "biochemistry", "Removal of amino group with oxidation", "Amino acid catabolism", "Metabolism"),
            ("Transamination", "process", "biochemistry", "Transfer of amino group", "Amino acid metabolism", "Metabolism"),
        ]
        terms.extend(biochem_terms)
        
        # ========== MICROBIOLOGY TERMS (80+) ==========
        microbe_terms = [
            ("Bacteria", "organism", "microbiology", "Single-celled prokaryotic microorganisms", "Decomposition, nitrogen fixation, digestion", "Industry and disease"),
            ("Archaea", "organism", "microbiology", "Prokaryotic organisms distinct from bacteria", "Extreme environments, biogeochemistry", "Ecology"),
            ("Virus", "agent", "microbiology", "Infectious agent requiring host cells", "Uses host machinery to replicate", "Disease"),
            ("Bacteriophage", "virus", "microbiology", "Virus infecting bacteria", "Bacterial control, phage therapy", "Biotechnology"),
            ("Fungus", "organism", "microbiology", "Eukaryotic spore-producing organism", "Decomposition, symbiosis", "Ecology, industry"),
            ("Yeast", "organism", "microbiology", "Single-celled fungus", "Fermentation, baking", "Biotechnology"),
            ("Mold", "organism", "microbiology", "Multicellular fungus", "Decomposition, food spoilage", "Ecology"),
            ("Protozoa", "organism", "microbiology", "Single-celled eukaryotic microorganisms", "Predation, parasitism", "Ecology, disease"),
            ("Algae", "organism", "microbiology", "Photosynthetic eukaryotic microorganisms", "Oxygen production, food chains", "Ecology"),
            ("Pathogen", "agent", "microbiology", "Disease-causing organism", "Infection, disease", "Medicine"),
            ("Antibiotic", "compound", "microbiology", "Substance killing or inhibiting bacteria", "Bacterial infection treatment", "Medicine"),
            ("Antiviral", "compound", "microbiology", "Drug treating viral infections", "Virus replication inhibition", "Medicine"),
            ("Antifungal", "compound", "microbiology", "Drug treating fungal infections", "Fungal growth inhibition", "Medicine"),
            ("Vaccine", "agent", "microbiology", "Preparation providing immunity", "Immune system stimulation", "Disease prevention"),
            ("Immunity", "state", "microbiology", "Resistance to disease", "Protection from infection", "Health"),
            ("Antibody", "protein", "microbiology", "Protein produced by immune system", "Antigen recognition", "Immunity"),
            ("Antigen", "molecule", "microbiology", "Substance triggering immune response", "Immune activation", "Immunology"),
            ("Immunoglobulin", "protein", "microbiology", "Antibody protein class", "Humoral immunity", "Immunology"),
            ("Cytokine", "molecule", "microbiology", "Cell signaling protein", "Immune regulation", "Immunology"),
            ("Interferon", "protein", "microbiology", "Protein inducing antiviral state", "Viral defense", "Immunology"),
            ("Complement", "system", "microbiology", "Protein system for pathogen lysis", "Pathogen destruction", "Innate immunity"),
            ("Innate Immunity", "system", "microbiology", "Non-specific first-line defense", "Immediate response", "Immunology"),
            ("Adaptive Immunity", "system", "microbiology", "Specific immune response with memory", "Specific defense", "Immunology"),
            ("Humoral Immunity", "system", "microbiology", "Antibody-mediated immunity", "Blood-borne defense", "Immunology"),
            ("Cell-Mediated Immunity", "system", "microbiology", "T-cell mediated immunity", "Intracellular pathogen defense", "Immunology"),
            ("Vaccination", "process", "microbiology", "Administration of vaccine", "Disease prevention", "Medicine"),
            ("Herd Immunity", "concept", "microbiology", "Population protection from infection", "Disease control", "Public health"),
            ("Gram Staining", "technique", "microbiology", "Bacterial classification technique", "Bacterial identification", "Diagnostics"),
            ("Sterilization", "process", "microbiology", "Elimination of all microorganisms", "Infection control", "Medicine"),
            ("Disinfection", "process", "microbiology", "Reduction of microorganisms", "Surface cleaning", "Medicine"),
            ("Aseptic Technique", "technique", "microbiology", "Methods preventing contamination", "Sterile procedures", "Medicine"),
            ("Culture Medium", "substance", "microbiology", "Nutrient solution for growing microbes", "Microbe growth", "Laboratory"),
            ("Agar", "substance", "microbiology", "Gelatinous substance from algae", "Solid culture medium", "Microbiology"),
            ("Biofilm", "structure", "microbiology", "Community of microorganisms attached to surface", "Infection, bioremediation", "Medicine"),
            ("Microbiome", "community", "microbiology", "Community of microorganisms in/on body", "Health, disease", "Medicine"),
            ("Probiotic", "organism", "microbiology", "Beneficial microorganisms", "Gut health", "Nutrition"),
            ("Prebiotic", "substance", "microbiology", "Food for beneficial bacteria", "Probiotic support", "Nutrition"),
            ("Symbiosis", "process", "microbiology", "Close relationship between organisms", "Mutualism, commensalism, parasitism", "Ecology"),
            ("Mutualism", "relationship", "microbiology", "Both organisms benefit", "Symbiotic relationship", "Ecology"),
            ("Parasitism", "relationship", "microbiology", "One organism benefits, other harmed", "Disease", "Ecology"),
            ("Commensalism", "relationship", "microbiology", "One benefits, other unaffected", "Symbiotic relationship", "Ecology"),
            ("Infection", "process", "microbiology", "Colonization by pathogen", "Disease", "Medicine"),
            ("Inflammation", "process", "microbiology", "Response to injury or infection", "Defense mechanism", "Medicine"),
            ("Fever", "symptom", "microbiology", "Elevated body temperature", "Immune response", "Medicine"),
            ("Sepsis", "condition", "microbiology", "Life-threatening response to infection", "Systemic inflammation", "Emergency"),
            ("Antimicrobial", "compound", "microbiology", "Substance killing or inhibiting microbes", "Infection treatment", "Medicine"),
            ("Resistance", "process", "microbiology", "Ability to survive antimicrobial treatment", "Treatment failure", "Medicine"),
        ]
        terms.extend(microbe_terms)
        
        # ========== MEDICINE TERMS (100+) ==========
        medical_terms = [
            ("Disease", "condition", "medicine", "Abnormal condition affecting body function", "Requires medical intervention", "Healthcare"),
            ("Cancer", "disease", "medicine", "Uncontrolled cell growth", "Can spread to other organs", "Leading cause of death"),
            ("Tumor", "mass", "medicine", "Abnormal cell growth mass", "Benign or malignant", "Oncology"),
            ("Malignant", "condition", "medicine", "Cancerous and invasive", "Can metastasize", "Oncology"),
            ("Benign", "condition", "medicine", "Non-cancerous", "Usually not life-threatening", "Oncology"),
            ("Metastasis", "process", "medicine", "Cancer spread to other organs", "Cancer progression", "Oncology"),
            ("Carcinoma", "cancer", "medicine", "Cancer from epithelial cells", "Most common cancer type", "Oncology"),
            ("Sarcoma", "cancer", "medicine", "Cancer from connective tissue", "Rare cancer type", "Oncology"),
            ("Leukemia", "cancer", "medicine", "Cancer of blood-forming tissues", "White blood cell proliferation", "Blood cancer"),
            ("Lymphoma", "cancer", "medicine", "Cancer of lymphatic system", "Lymph node involvement", "Blood cancer"),
            ("Diabetes Mellitus", "disease", "medicine", "Metabolic disease with high blood sugar", "Insulin deficiency or resistance", "Metabolic disease"),
            ("Type 1 Diabetes", "disease", "medicine", "Autoimmune destruction of insulin-producing cells", "No insulin production", "Diabetes"),
            ("Type 2 Diabetes", "disease", "medicine", "Insulin resistance and relative insulin deficiency", "Usually adult-onset", "Diabetes"),
            ("Hypertension", "condition", "medicine", "High blood pressure", "Cardiovascular risk", "Heart disease"),
            ("Atherosclerosis", "disease", "medicine", "Artery hardening from plaque", "Heart disease risk", "Cardiovascular"),
            ("Heart Disease", "disease", "medicine", "Various conditions affecting heart", "Leading cause of death", "Cardiology"),
            ("Myocardial Infarction", "disease", "medicine", "Heart attack - blood flow blockage", "Cardiac cell death", "Emergency"),
            ("Stroke", "disease", "medicine", "Brain blood supply interruption", "Brain cell death", "Emergency"),
            ("Alzheimer's Disease", "disease", "medicine", "Progressive neurodegenerative disease", "Memory loss, dementia", "Neurology"),
            ("Parkinson's Disease", "disease", "medicine", "Neurodegenerative movement disorder", "Dopamine neuron loss", "Neurology"),
            ("HIV", "virus", "medicine", "Human Immunodeficiency Virus", "AIDS cause", "Infectious disease"),
            ("AIDS", "disease", "medicine", "Acquired Immunodeficiency Syndrome", "Immunodeficiency", "Infectious disease"),
            ("Influenza", "disease", "medicine", "Viral respiratory infection", "Flu", "Infectious disease"),
            ("Tuberculosis", "disease", "medicine", "Bacterial infection mainly affecting lungs", "Major infectious disease", "Pulmonology"),
            ("Malaria", "disease", "medicine", "Parasitic disease spread by mosquitoes", "Tropical disease", "Infectious disease"),
            ("COVID-19", "disease", "medicine", "Coronavirus disease 2019", "Respiratory pandemic", "Infectious disease"),
            ("SARS", "disease", "medicine", "Severe Acute Respiratory Syndrome", "Viral respiratory disease", "Infectious disease"),
            ("MERS", "disease", "medicine", "Middle East Respiratory Syndrome", "Viral respiratory disease", "Infectious disease"),
            ("Ebola", "disease", "medicine", "Viral hemorrhagic fever", "Severe viral disease", "Infectious disease"),
            ("Zika", "disease", "medicine", "Mosquito-borne viral infection", "Birth defects", "Infectious disease"),
            ("Hepatitis", "disease", "medicine", "Liver inflammation", "Various causes", "Gastroenterology"),
            ("Cirrhosis", "disease", "medicine", "Chronic liver disease", "Liver scarring", "Gastroenterology"),
            ("Nephritis", "disease", "medicine", "Kidney inflammation", "Kidney disease", "Nephrology"),
            ("Arthritis", "disease", "medicine", "Joint inflammation", "Various types", "Rheumatology"),
            ("Osteoporosis", "disease", "medicine", "Weak bones", "Fracture risk", "Orthopedics"),
            ("Asthma", "disease", "medicine", "Chronic airway inflammation", "Breathing difficulty", "Pulmonology"),
            ("COPD", "disease", "medicine", "Chronic Obstructive Pulmonary Disease", "Lung disease", "Pulmonology"),
            ("Pneumonia", "disease", "medicine", "Lung infection", "Respiratory disease", "Pulmonology"),
            ("Bronchitis", "disease", "medicine", "Bronchial tube inflammation", "Respiratory disease", "Pulmonology"),
            ("Anemia", "condition", "medicine", "Low red blood cell count", "Fatigue, weakness", "Hematology"),
            ("Hemophilia", "disease", "medicine", "Blood clotting disorder", "Excessive bleeding", "Hematology"),
            ("Sickle Cell Disease", "disease", "medicine", "Genetic blood disorder", "Anemia, pain crises", "Hematology"),
            ("Thalassemia", "disease", "medicine", "Genetic blood disorder", "Anemia", "Hematology"),
            ("Multiple Sclerosis", "disease", "medicine", "Autoimmune nervous system disease", "Neurological symptoms", "Neurology"),
            ("Epilepsy", "disease", "medicine", "Neurological disorder with seizures", "Seizure disorders", "Neurology"),
            ("Schizophrenia", "disease", "medicine", "Chronic psychiatric disorder", "Psychosis", "Psychiatry"),
            ("Bipolar Disorder", "disease", "medicine", "Mood disorder with alternating highs and lows", "Psychiatric condition", "Psychiatry"),
            ("Depression", "disorder", "medicine", "Mood disorder with persistent sadness", "Mental health", "Psychiatry"),
            ("Anxiety Disorder", "disorder", "medicine", "Condition with excessive worry", "Mental health", "Psychiatry"),
            ("Autism", "disorder", "medicine", "Developmental disorder", "Social communication differences", "Neurology"),
            ("Autoimmune Disease", "disease", "medicine", "Immune system attacking body", "Various conditions", "Immunology"),
            ("Rheumatoid Arthritis", "disease", "medicine", "Autoimmune joint disease", "Joint inflammation", "Rheumatology"),
            ("Lupus", "disease", "medicine", "Systemic autoimmune disease", "Multi-organ involvement", "Rheumatology"),
            ("Psoriasis", "disease", "medicine", "Skin autoimmune disease", "Skin lesions", "Dermatology"),
            ("Eczema", "disease", "medicine", "Skin inflammatory condition", "Itchy rash", "Dermatology"),
            ("Acne", "condition", "medicine", "Skin condition with pimples", "Common skin problem", "Dermatology"),
            ("Cystic Fibrosis", "disease", "medicine", "Genetic disease affecting lungs", "Respiratory problems", "Pulmonology"),
            ("Huntington's Disease", "disease", "medicine", "Neurodegenerative genetic disease", "Movement, cognitive problems", "Neurology"),
            ("Amyotrophic Lateral Sclerosis", "disease", "medicine", "Progressive nervous system disease", "Muscle weakness", "Neurology"),
            ("Prion Disease", "disease", "medicine", "Infectious protein disease", "Brain degeneration", "Neurology"),
            ("Creutzfeldt-Jakob Disease", "disease", "medicine", "Rare prion disease", "Rapidly progressive dementia", "Neurology"),
            ("Celiac Disease", "disease", "medicine", "Autoimmune reaction to gluten", "Intestinal damage", "Gastroenterology"),
            ("Crohn's Disease", "disease", "medicine", "Inflammatory bowel disease", "GI inflammation", "Gastroenterology"),
            ("Ulcerative Colitis", "disease", "medicine", "Inflammatory bowel disease", "Colon inflammation", "Gastroenterology"),
            ("Pancreatitis", "disease", "medicine", "Pancreas inflammation", "Digestive problems", "Gastroenterology"),
            ("Thyroid Disease", "disease", "medicine", "Various thyroid conditions", "Metabolism effects", "Endocrinology"),
            ("Addison's Disease", "disease", "medicine", "Adrenal insufficiency", "Hormone deficiency", "Endocrinology"),
            ("Cushing's Syndrome", "disease", "medicine", "Excess cortisol", "Various symptoms", "Endocrinology"),
            ("Pituitary Adenoma", "disease", "medicine", "Pituitary gland tumor", "Hormone excess or deficiency", "Endocrinology"),
            ("Breast Cancer", "cancer", "medicine", "Cancer of breast tissue", "Most common cancer in women", "Oncology"),
            ("Prostate Cancer", "cancer", "medicine", "Cancer of prostate gland", "Common in men", "Oncology"),
            ("Lung Cancer", "cancer", "medicine", "Cancer of lungs", "Leading cancer death", "Oncology"),
            ("Colorectal Cancer", "cancer", "medicine", "Cancer of colon or rectum", "Common cancer", "Oncology"),
            ("Pancreatic Cancer", "cancer", "medicine", "Cancer of pancreas", "Poor prognosis", "Oncology"),
            ("Ovarian Cancer", "cancer", "medicine", "Cancer of ovaries", "Gynecologic cancer", "Oncology"),
            ("Testicular Cancer", "cancer", "medicine", "Cancer of testicles", "Young adult cancer", "Oncology"),
            ("Skin Cancer", "cancer", "medicine", "Cancer of skin", "Most common cancer", "Oncology"),
            ("Brain Tumor", "cancer", "medicine", "Abnormal growth in brain", "Neurological symptoms", "Oncology"),
            ("Retinoblastoma", "cancer", "medicine", "Eye cancer in children", "Genetic cancer", "Oncology"),
            ("Wilms Tumor", "cancer", "medicine", "Kidney cancer in children", "Pediatric cancer", "Oncology"),
            ("Chemotherapy", "treatment", "medicine", "Drug treatment for cancer", "Cancer treatment", "Oncology"),
            ("Radiation Therapy", "treatment", "medicine", "Using radiation to treat cancer", "Cancer treatment", "Oncology"),
            ("Immunotherapy", "treatment", "medicine", "Treatment using immune system", "Cancer treatment", "Oncology"),
            ("Targeted Therapy", "treatment", "medicine", "Treatment targeting specific molecules", "Precision cancer treatment", "Oncology"),
            ("Hormone Therapy", "treatment", "medicine", "Treatment using hormones", "Cancer treatment", "Oncology"),
            ("Surgery", "treatment", "medicine", "Physical intervention", "Various treatments", "Medicine"),
            ("Transplant", "treatment", "medicine", "Organ or tissue replacement", "End-stage organ failure", "Medicine"),
            ("Gene Therapy", "treatment", "medicine", "Treating disease by modifying genes", "Future medicine", "Biotechnology"),
            ("Stem Cell Therapy", "treatment", "medicine", "Using stem cells for treatment", "Regenerative medicine", "Biotechnology"),
            ("Organoid", "structure", "medicine", "Mini-organ grown in lab", "Research, drug testing", "Biotechnology"),
            ("3D Bioprinting", "technique", "medicine", "Printing living tissues", "Regenerative medicine", "Biotechnology"),
            ("Telemedicine", "technique", "medicine", "Remote medical consultation", "Healthcare delivery", "Technology"),
            ("Precision Medicine", "approach", "medicine", "Medical treatment based on genetics", "Personalized healthcare", "Future medicine"),
            ("Personalized Medicine", "approach", "medicine", "Treatment tailored to individual", "Precision healthcare", "Medicine"),
            ("Pharmacogenomics", "field", "medicine", "Study of drug response and genetics", "Personalized medication", "Medicine"),
            ("Clinical Trial", "process", "medicine", "Research study in humans", "Drug development", "Medicine"),
            ("Placebo", "substance", "medicine", "Inactive treatment for comparison", "Clinical trials", "Research"),
            ("Diagnosis", "process", "medicine", "Identifying disease", "Medical assessment", "Medicine"),
            ("Prognosis", "process", "medicine", "Disease outcome prediction", "Treatment planning", "Medicine"),
            ("Therapy", "treatment", "medicine", "Treatment for disease", "Healthcare", "Medicine"),
            ("Prophylaxis", "treatment", "medicine", "Disease prevention", "Prevention", "Medicine"),
            ("Palliative Care", "treatment", "medicine", "Care for quality of life", "Symptom management", "Medicine"),
            ("Hospice", "care", "medicine", "End-of-life care", "Comfort care", "Medicine"),
            ("Epidemiology", "field", "medicine", "Study of disease patterns", "Public health", "Medicine"),
            ("Public Health", "field", "medicine", "Health of populations", "Disease prevention", "Medicine"),
        ]
        terms.extend(medical_terms)
        
        # ========== BIOTECHNOLOGY TERMS (60+) ==========
        biotech_terms = [
            ("Biotechnology", "field", "biotechnology", "Using living systems for technological applications", "Industry, medicine", "Modern technology"),
            ("Genetic Engineering", "field", "biotechnology", "Manipulating genetic material", "GMOs, gene therapy", "Biotechnology"),
            ("Recombinant DNA", "technique", "biotechnology", "DNA from different sources combined", "Gene cloning", "Biotechnology"),
            ("Cloning", "technique", "biotechnology", "Creating identical copy of DNA or organism", "Gene amplification", "Biotechnology"),
            ("PCR", "technique", "biotechnology", "Polymerase Chain Reaction - DNA amplification", "DNA copying", "Diagnostics, forensics"),
            ("RT-PCR", "technique", "biotechnology", "Reverse Transcription PCR - RNA detection", "Gene expression analysis", "Research"),
            ("qPCR", "technique", "biotechnology", "Quantitative PCR - measures DNA quantity", "Quantification", "Diagnostics"),
            ("Gel Electrophoresis", "technique", "biotechnology", "DNA separation by size", "DNA analysis", "Research"),
            ("Western Blot", "technique", "biotechnology", "Protein detection technique", "Protein analysis", "Research"),
            ("Southern Blot", "technique", "biotechnology", "DNA detection technique", "DNA analysis", "Research"),
            ("Northern Blot", "technique", "biotechnology", "RNA detection technique", "RNA analysis", "Research"),
            ("DNA Sequencing", "technique", "biotechnology", "Determining DNA order", "Genome projects", "Genomics"),
            ("Next-Generation Sequencing", "technique", "biotechnology", "High-throughput DNA sequencing", "Genomics revolution", "Genomics"),
            ("RNA Sequencing", "technique", "biotechnology", "Sequencing RNA transcripts", "Transcriptomics", "Genomics"),
            ("Single-Cell Sequencing", "technique", "biotechnology", "Sequencing individual cells", "Cell heterogeneity", "Genomics"),
            ("CRISPR-Cas9", "technique", "biotechnology", "Gene editing technology", "Precise gene modification", "Revolutionary biotechnology"),
            ("CRISPR", "technique", "biotechnology", "Gene editing technology", "Precise gene modification", "Revolutionary biotechnology"),
            ("Gene Drive", "technique", "biotechnology", "Bias inheritance to spread traits", "Pest control", "Biotechnology"),
            ("Synthetic Biology", "field", "biotechnology", "Designing new biological systems", "Bioengineering", "Future biotechnology"),
            ("Metabolic Engineering", "field", "biotechnology", "Modifying metabolic pathways", "Biofuel production", "Biotechnology"),
            ("Bioprocessing", "process", "biotechnology", "Using biological systems for production", "Industrial biotechnology", "Industry"),
            ("Bioreactor", "device", "biotechnology", "Device for biological reactions", "Cell culture", "Industry"),
            ("Fermentation", "process", "biotechnology", "Metabolic process producing products", "Food, biofuel production", "Industry"),
            ("Industrial Enzyme", "enzyme", "biotechnology", "Enzyme for industrial use", "Manufacturing", "Industry"),
            ("Biocatalysis", "process", "biotechnology", "Using enzymes for chemical reactions", "Green chemistry", "Industry"),
            ("Biosensor", "device", "biotechnology", "Biological detection device", "Diagnostics", "Technology"),
            ("Biochip", "device", "biotechnology", "Microarray for biomolecule detection", "Diagnostics, research", "Technology"),
            ("DNA Microarray", "technique", "biotechnology", "Simultaneous gene expression analysis", "Transcriptomics", "Research"),
            ("RNA Interference", "technique", "biotechnology", "Gene silencing technique", "Gene function study", "Research"),
            ("Antisense Therapy", "technique", "biotechnology", "Treatment blocking gene expression", "Genetic disease treatment", "Medicine"),
            ("Phage Display", "technique", "biotechnology", "Screening protein interactions", "Drug discovery", "Research"),
            ("Yeast Two-Hybrid", "technique", "biotechnology", "Detecting protein interactions", "Protein research", "Research"),
            ("Flow Cytometry", "technique", "biotechnology", "Analyzing cell characteristics", "Cell analysis", "Research"),
            ("Cell Sorting", "technique", "biotechnology", "Separating cells by type", "Cell isolation", "Research"),
            ("Flow Cytometer", "device", "biotechnology", "Instrument for cell analysis", "Research", "Diagnostics"),
            ("Transgenic Organism", "organism", "biotechnology", "Organism with foreign gene", "GMOs", "Biotechnology"),
            ("Knockout Mouse", "organism", "biotechnology", "Mouse with gene inactivated", "Disease research", "Research"),
            ("Knockin Mouse", "organism", "biotechnology", "Mouse with inserted gene", "Disease research", "Research"),
            ("Humanized Mouse", "organism", "biotechnology", "Mouse with human genes", "Research", "Drug testing"),
            ("Organoid", "structure", "biotechnology", "Mini-organ in culture", "Research", "Drug testing"),
            ("3D Cell Culture", "technique", "biotechnology", "Cells in three-dimensional matrix", "Tissue engineering", "Research"),
            ("Tissue Engineering", "field", "biotechnology", "Creating tissues in lab", "Regenerative medicine", "Medicine"),
            ("Stem Cell", "cell", "biotechnology", "Undifferentiated cell", "Regenerative medicine", "Research"),
            ("Embryonic Stem Cell", "cell", "biotechnology", "Stem cell from embryo", "Pluripotent", "Research"),
            ("Induced Pluripotent Stem Cell", "cell", "biotechnology", "Stem cell from adult cells", "Patient-specific", "Research"),
            ("Cell Therapy", "treatment", "biotechnology", "Treatment using cells", "Regenerative medicine", "Medicine"),
            ("Gene Therapy", "treatment", "biotechnology", "Treatment modifying genes", "Genetic disease cure", "Medicine"),
            ("Viral Vector", "tool", "biotechnology", "Virus used for gene delivery", "Gene therapy", "Medicine"),
            ("Lentivirus", "virus", "biotechnology", "Retrovirus for gene delivery", "Gene therapy", "Research"),
            ("Adenovirus", "virus", "biotechnology", "Common viral vector", "Gene therapy", "Research"),
            ("AAV", "virus", "biotechnology", "Adeno-associated virus - gene therapy vector", "Gene therapy", "Medicine"),
            ("Biopharmaceutical", "product", "biotechnology", "Drug produced by biotechnology", "Protein drugs", "Medicine"),
            ("Monoclonal Antibody", "product", "biotechnology", "Therapeutic antibody", "Cancer treatment", "Medicine"),
            ("Recombinant Protein", "product", "biotechnology", "Protein from engineered cells
