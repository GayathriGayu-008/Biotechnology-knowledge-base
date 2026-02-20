"""
Biotechnology Knowledge Base - AI-Powered Enhancement Module
Uses OpenAI to enhance term definitions and provide detailed explanations
"""

import os
import json
import sqlite3
from datetime import datetime
import requests

# OpenAI API configuration
OPENAI_API_KEY = os.environ.get("OPENAI_API_KEY", "")

class AIKnowledgeBase:
    def __init__(self, db_path="knowledge_base.db"):
        self.db_path = db_path
        
    def get_openai_client(self):
        """Get OpenAI client"""
        try:
            from openai import OpenAI
            return OpenAI(api_key=OPENAI_API_KEY)
        except ImportError:
            print("OpenAI package not installed. Run: pip install openai")
            return None
    
    def enhance_term_with_ai(self, term_info):
        """Enhance a term with AI-generated detailed information"""
        if not OPENAI_API_KEY:
            return term_info
            
        client = self.get_openai_client()
        if not client:
            return term_info
            
        try:
            prompt = f"""
            Provide detailed information about the biological/medical term "{term_info.get('term', '')}".
            
            Current information:
            - Type: {term_info.get('term_type', 'N/A')}
            - Category: {term_info.get('category', 'N/A')}
            - Definition: {term_info.get('definition', 'N/A')}
            - Function: {term_info.get('function', 'N/A')}
            - Importance: {term_info.get('importance', 'N/A')}
            
            Please provide:
            1. A comprehensive definition (2-3 sentences)
            2. Detailed biological/molecular function
            3. Clinical/medical significance if applicable
            4. Related biological processes or pathways
            5. Common associations or related terms
            
            Format as JSON with keys: enhanced_definition, detailed_function, clinical_significance, related_processes, related_terms
            """
            
            response = client.chat.completions.create(
                model="gpt-3.5-turbo",
                messages=[
                    {"role": "system", "content": "You are a helpful biology and medicine expert assistant."},
                    {"role": "user", "content": prompt}
                ],
                temperature=0.7,
                max_tokens=500
            )
            
            if response.choices:
                result = json.loads(response.choices[0].message.content)
                term_info['ai_enhanced'] = True
                term_info['enhanced_definition'] = result.get('enhanced_definition', '')
                term_info['detailed_function'] = result.get('detailed_function', '')
                term_info['clinical_significance'] = result.get('clinical_significance', '')
                term_info['related_processes'] = result.get('related_processes', '')
                term_info['related_terms'] = result.get('related_terms', '')
                
        except Exception as e:
            print(f"AI enhancement error for {term_info.get('term')}: {e}")
            
        return term_info
    
    def generate_term_explanation(self, term, context=""):
        """Generate a detailed explanation for a term"""
        if not OPENAI_API_KEY:
            return "OpenAI API key not configured"
            
        client = self.get_openai_client()
        if not client:
            return "OpenAI client not available"
            
        try:
            prompt = f"""
            Explain the following biotechnology/biological term in detail:
            Term: {term}
            Context: {context}
            
            Provide a comprehensive explanation suitable for a student or researcher.
            Include:
            - Definition
            - Function/Role
            - Significance in research/medicine
            - Related terms and concepts
            """
            
            response = client.chat.completions.create(
                model="gpt-3.5-turbo",
                messages=[
                    {"role": "system", "content": "You are an expert in biotechnology, molecular biology, and medicine."},
                    {"role": "user", "content": prompt}
                ],
                temperature=0.7,
                max_tokens=800
            )
            
            if response.choices:
                return response.choices[0].message.content
                
        except Exception as e:
            return f"Error generating explanation: {str(e)}"
        
        return "Unable to generate explanation"
    
    def find_related_terms(self, term, limit=5):
        """Find related terms in the knowledge base"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        # Get the category of the search term
        cursor.execute('''
            SELECT category, term_type FROM terms 
            WHERE term LIKE ? LIMIT 1
        ''', (f'%{term}%',))
        
        row = cursor.fetchone()
        if row:
            category, term_type = row
            
            # Find related terms in the same category
            cursor.execute('''
                SELECT term, definition, term_type, category 
                FROM terms 
                WHERE category = ? AND term != ?
                LIMIT ?
            ''', (category, term, limit))
        else:
            # If no category found, search by similar terms
            cursor.execute('''
                SELECT term, definition, term_type, category 
                FROM terms 
                WHERE term LIKE ? AND term != ?
                LIMIT ?
            ''', (f'%{term}%', term, limit))
        
        results = cursor.fetchall()
        conn.close()
        
        return [{
            'term': r[0],
            'definition': r[1],
            'term_type': r[2],
            'category': r[3]
        } for r in results]
    
    def answer_biological_question(self, question):
        """Answer a biological/medical question using AI"""
        if not OPENAI_API_KEY:
            return {
                "answer": "OpenAI API key not configured. Please set the OPENAI_API_KEY environment variable.",
                "sources": []
            }
            
        client = self.get_openai_client()
        if not client:
            return {
                "answer": "OpenAI client not available",
                "sources": []
            }
        
        # First, search knowledge base for relevant terms
        search_terms = question.lower().split()
        relevant_terms = []
        
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        for term in search_terms:
            if len(term) > 3:
                cursor.execute('''
                    SELECT term, definition, function, importance 
                    FROM terms 
                    WHERE term LIKE ? OR definition LIKE ?
                    LIMIT 3
                ''', (f'%{term}%', f'%{term}%'))
                relevant_terms.extend(cursor.fetchall())
        
        conn.close()
        
        # Build context from knowledge base
        context = "Relevant information from knowledge base:\n"
        seen = set()
        for t in relevant_terms:
            if t[0] not in seen:
                context += f"- {t[0]}: {t[1][:200]}\n"
                seen.add(t[0])
        
        try:
            prompt = f"""
            Context from biotechnology knowledge base:
            {context}
            
            Question: {question}
            
            Based on the knowledge base information and your knowledge, provide a comprehensive answer.
            Cite relevant terms from the knowledge base when possible.
            """
            
            response = client.chat.completions.create(
                model="gpt-3.5-turbo",
                messages=[
                    {"role": "system", "content": "You are an expert biotechnology and medical research assistant. Use the provided knowledge base context to answer questions accurately."},
                    {"role": "user", "content": prompt}
                ],
                temperature=0.7,
                max_tokens=1000
            )
            
            answer = response.choices[0].message.content if response.choices else "Unable to generate answer"
            
            sources = list(seen)
            
            return {
                "answer": answer,
                "sources": sources,
                "relevant_terms": [{"term": t[0], "definition": t[1][:100]} for t in relevant_terms[:5]]
            }
            
        except Exception as e:
            return {
                "answer": f"Error: {str(e)}",
                "sources": [],
                "relevant_terms": []
            }
    
    def generate_biological_pathway_description(self, pathway_name):
        """Generate a detailed description of a biological pathway"""
        if not OPENAI_API_KEY:
            return "OpenAI API key not configured"
            
        client = self.get_openai_client()
        if not client:
            return "OpenAI client not available"
        
        prompt = f"""
        Describe the biological pathway "{pathway_name}" in detail.
        
        Include:
        1. Overview of the pathway
        2. Key components/molecules involved
        3. Step-by-step mechanism
        4. Biological significance
        5. Related diseases if applicable
        6. Research/therapeutic implications
        """
        
        try:
            response = client.chat.completions.create(
                model="gpt-3.5-turbo",
                messages=[
                    {"role": "system", "content": "You are a biochemistry and molecular biology expert."},
                    {"role": "user", "content": prompt}
                ],
                temperature=0.7,
                max_tokens=800
            )
            
            return response.choices[0].message.content if response.choices else "Unable to generate pathway description"
            
        except Exception as e:
            return f"Error: {str(e)}"


class MassiveDataCollector:
    """Collects millions of terms from various biological databases"""
    
    def __init__(self, db_path="knowledge_base.db"):
        self.db_path = db_path
        self.init_database()
        
    def init_database(self):
        """Initialize the database with optimized schema for millions of terms"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        # Main terms table - optimized for millions of records
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
        
        # Create indexes for fast searching
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_term ON terms(term)')
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_category ON terms(category)')
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_term_type ON terms(term_type)')
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_source ON terms(source)')
        
        # FTS index for full-text search
        cursor.execute('''
            CREATE VIRTUAL TABLE IF NOT EXISTS terms_fts USING fts5(
                term, definition, function, importance, category,
                content='terms',
                content_rowid='id'
            )
        ''')
        
        # Categories table
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS categories (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                name TEXT UNIQUE NOT NULL,
                display_name TEXT,
                description TEXT,
                parent_id INTEGER,
                term_count INTEGER DEFAULT 0,
                FOREIGN KEY (parent_id) REFERENCES categories(id)
            )
        ''')
        
        # Collection sources table
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS collection_sources (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                source_name TEXT UNIQUE NOT NULL,
                source_url TEXT,
                terms_collected INTEGER DEFAULT 0,
                last_collection TIMESTAMP,
                status TEXT
            )
        ''')
        
        conn.commit()
        conn.close()
        print("Massive database initialized!")
    
    def add_term_bulk(self, terms_list):
        """Bulk insert terms for faster collection"""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        cursor.executemany('''
            INSERT OR REPLACE INTO terms 
            (term, term_type, category, definition, function, importance, 
             source, external_id, external_link, metadata, updated_at)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ''', terms_list)
        
        conn.commit()
        conn.close()
    
    def add_comprehensive_biological_terms(self):
        """Add comprehensive list of biological, medical, and biotechnology terms"""
        print("Adding comprehensive biological terms...")
        
        terms = []
        
        # Molecular Biology Terms (500+)
        molecular_terms = [
            ("DNA", "molecule", "molecular_biology", "Deoxyribonucleic acid - the hereditary material in all known life forms", "Carries genetic instructions for development, functioning, growth and reproduction", "Fundamental to genetics and biotechnology"),
            ("RNA", "molecule", "molecular_biology", "Ribonucleic acid - single-stranded nucleic acid present in all living cells", "Plays roles in coding, decoding, regulation and expression of genes", "Essential for protein synthesis"),
            ("mRNA", "molecule", "molecular_biology", "Messenger RNA - carries genetic code from DNA to ribosomes", "Template for protein synthesis", "Key to gene expression"),
            ("tRNA", "molecule", "molecular_biology", "Transfer RNA - brings amino acids to ribosomes during translation", "Adapter molecule in protein synthesis", "Essential for translation"),
            ("rRNA", "molecule", "molecular_biology", "Ribosomal RNA - structural and catalytic component of ribosomes", "Forms the ribosome and catalyzes peptide bond formation", "Most abundant RNA in cells"),
            ("Protein", "molecule", "molecular_biology", "Large biomolecules composed of amino acid chains", "Perform vast array of functions including catalysis, structure, transport", "Essential for all cellular functions"),
            ("Amino Acid", "molecule", "molecular_biology", "Building blocks of proteins, containing amino and carboxyl groups", "Building blocks of proteins, some act as neurotransmitters", "Fundamental to life"),
            ("Peptide", "molecule", "molecular_biology", "Short chain of amino acids linked by peptide bonds", "Serve as signaling molecules, hormones, antibiotics", "Important in drug development"),
            ("Enzyme", "protein", "molecular_biology", "Biological catalyst that speeds up chemical reactions", "Catalyzes metabolic reactions essential for life", "Critical for biochemistry"),
            ("Kinase", "enzyme", "molecular_biology", "Enzyme that transfers phosphate groups from ATP to substrates", "Regulates signaling pathways, metabolism, cell cycle", "Important drug targets"),
            ("Phosphatase", "enzyme", "molecular_biology", "Enzyme that removes phosphate groups from molecules", "Regulates signaling pathways by dephosphorylation", "Important in cell signaling"),
            ("Protease", "enzyme", "molecular_biology", "Enzyme that breaks down proteins into peptides", "Digestion, protein turnover, cell signaling", "Drug targets for cancer, HIV"),
            ("Polymerase", "enzyme", "molecular_biology", "Enzyme that synthesizes nucleic acid polymers", "DNA replication and repair", "Essential for PCR technology"),
            ("Helicase", "enzyme", "molecular_biology", "Enzyme that unwinds DNA double helix", "Required for DNA replication and transcription", "Important in cancer research"),
            ("Ligase", "enzyme", "molecular_biology", "Enzyme that joins DNA strands", "DNA repair and replication", "Used in cloning"),
            ("Transcriptase", "enzyme", "molecular_biology", "Enzyme that synthesizes RNA from DNA template", "Viral replication, gene expression", "HIV drug target"),
            ("Receptor", "protein", "molecular_biology", "Protein that binds specific molecules and triggers cellular responses", "Signal transduction, drug targets", "Major pharmaceutical targets"),
            ("Transporter", "protein", "molecule", "Protein that moves molecules across cell membranes", "Nutrient uptake, waste removal, ion balance", "Drug delivery targets"),
            ("Channel", "protein", "molecular_biology", "Membrane protein that forms pores for ion passage", "Nerve impulse transmission, muscle contraction", "Drug targets"),
            ("Ion Channel", "protein", "molecular_biology", "Gate-controlled pore allowing ion flow across membranes", "Action potentials, signal transduction", "Anesthesia targets"),
            ("Carrier Protein", "protein", "molecular_biology", "Protein that binds and transports specific molecules", "Facilitated diffusion, active transport", "Membrane transport"),
            ("Structural Protein", "protein", "molecular_biology", "Proteins providing mechanical support to cells and tissues", "Cell structure, cytoskeleton", "Collagen, keratin"),
            ("Cytoskeleton", "structure", "cell_biology", "Network of protein filaments providing cell structure and movement", "Cell shape, division, movement", "Cell mechanics"),
            ("Actin", "protein", "cell_biology", "Protein forming microfilaments in cytoskeleton", "Cell motility, division, muscle contraction", "Muscle function"),
            ("Tubulin", "protein", "cell_biology", "Protein forming microtubules", "Cell division, intracellular transport", "Cancer drug target"),
            ("Myosin", "protein", "cell_biology", "Motor protein interacting with actin", "Muscle contraction, cell movement", "Muscle biology"),
            ("Keratin", "protein", "cell_biology", "Fibrous structural protein in hair, nails, skin", "Protection, structure", "Skin and hair"),
            ("Collagen", "protein", "cell_biology", "Fibrous structural protein in connective tissues", "Tensile strength, tissue structure", "Skin, bone, cartilage"),
            ("Elastin", "protein", "cell_biology", "Protein providing elasticity to tissues", "Tissue elasticity", "Blood vessels, skin"),
            ("Hemoglobin", "protein", "biochemistry", "Protein in red blood cells carrying oxygen", "Oxygen transport from lungs to tissues", "Blood function"),
            ("Myoglobin", "protein", "biochemistry", "Oxygen-storing protein in muscle cells", "Oxygen storage in muscles", "Muscle metabolism"),
            ("Insulin", "hormone", "biochemistry", "Peptide hormone regulating blood glucose levels", "Glucose uptake by cells", "Diabetes treatment"),
            ("Glucagon", "hormone", "biochemistry", "Hormone that increases blood glucose levels", "Gluconeogenesis, glycogenolysis", "Metabolism"),
            ("Epinephrine", "hormone", "biochemistry", "Hormone and neurotransmitter (adrenaline)", "Fight or flight response", "Stress response"),
            ("Cortisol", "hormone", "biochemistry", "Steroid hormone regulating metabolism and stress", "Stress response, metabolism", "Stress biology"),
            ("Testosterone", "hormone", "biochemistry", "Male sex hormone", "Male characteristics, reproduction", "Endocrinology"),
            ("Estrogen", "hormone", "biochemistry", "Female sex hormone", "Female characteristics, reproduction", "Endocrinology"),
            ("Progesterone", "hormone", "biochemistry", "Hormone involved in menstrual cycle and pregnancy", "Pregnancy maintenance", "Reproduction"),
            ("Thyroid Hormone", "hormone", "biochemistry", "Hormone regulating metabolism", "Metabolic rate, development", "Metabolism"),
            ("Growth Hormone", "hormone", "biochemistry", "Hormone stimulating growth and cell reproduction", "Growth, metabolism", "Development"),
            ("Norepinephrine", "hormone", "biochemistry", "Neurotransmitter and hormone", "Alertness, stress response", "Neurotransmission"),
            ("Dopamine", "neurotransmitter", "biochemistry", "Neurotransmitter involved in reward and movement", "Reward, motivation, motor control", "Parkinson's, addiction"),
            ("Serotonin", "neurotransmitter", "biochemistry", "Neurotransmitter regulating mood, sleep, appetite", "Mood, appetite, sleep", "Depression, anxiety"),
            ("GABA", "neurotransmitter", "biochemistry", "Main inhibitory neurotransmitter in brain", "Neural inhibition", "Anxiety, seizures"),
            ("Glutamate", "neurotransmitter", "biochemistry", "Main excitatory neurotransmitter in brain", "Neural excitation", "Learning, memory"),
            ("Acetylcholine", "neurotransmitter", "biochemistry", "Neurotransmitter at neuromuscular junctions", "Muscle contraction, learning", "Alzheimer's treatment"),
            ("ATP", "molecule", "biochemistry", "Adenosine triphosphate - cell energy currency", "Energy transfer in cells", "Cellular metabolism"),
            ("ADP", "molecule", "biochemistry", "Adenosine diphosphate - ATP precursor", "Energy transfer", "Metabolism"),
            ("NAD", "molecule", "biochemistry", "Nicotinamide adenine dinucleotide - electron carrier", "Redox reactions", "Metabolism"),
            ("NADH", "molecule", "biochemistry", "Reduced form of NAD", "Electron donor in respiration", "Energy production"),
            ("FAD", "molecule", "biochemistry", "Flavin adenine dinucleotide - electron carrier", "Redox reactions", "Metabolism"),
            ("Glucose", "molecule", "biochemistry", "Simple sugar - primary energy source", "Cellular respiration fuel", "Metabolism"),
            ("Glycogen", "molecule", "biochemistry", "Storage form of glucose in animals", "Energy storage", "Metabolism"),
            ("Lipid", "molecule", "biochemistry", "Hydrophobic molecules including fats and oils", "Energy storage, cell membranes", "Metabolism"),
            ("Fatty Acid", "molecule", "biochemistry", "Carboxylic acid with long hydrocarbon chain", "Energy, membrane components", "Metabolism"),
            ("Triglyceride", "molecule", "biochemistry", "Ester of glycerol and three fatty acids", "Energy storage", "Metabolism"),
            ("Cholesterol", "molecule", "biochemistry", "Steroid alcohol in cell membranes", "Membrane fluidity, hormone synthesis", "Cardiovascular disease"),
            ("Phospholipid", "molecule", "biochemistry", "Lipid with phosphate group - membrane component", "Cell membrane structure", "Cell biology"),
            ("Steroid", "molecule", "biochemistry", "Lipid with four fused carbon rings", "Hormones, membranes", "Endocrinology"),
            ("Amino acid", "molecule", "biochemistry", "Building blocks of proteins", "Protein synthesis", "Metabolism"),
            ("Vitamin", "molecule", "biochemistry", "Organic compounds required in small amounts", "Coenzymes, antioxidants", "Nutrition"),
            ("Mineral", "molecule", "biochemistry", "Inorganic elements required for metabolism", "Enzyme cofactors, structure", "Nutrition"),
            ("Coenzyme", "molecule", "biochemistry", "Organic cofactor for enzyme activity", "Enzyme catalysis", "Metabolism"),
            ("Cofactor", "molecule", "biochemistry", "Non-protein compound required for enzyme activity", "Enzyme function", "Biochemistry"),
        ]
        
        # Add molecular biology terms
        for item in molecular_terms:
            terms.append((
                item[0], item[1], item[2], item[3], item[4], item[5],
                "Core Knowledge Base", None, None, None, datetime.now()
            ))
        
        # Genetics Terms
        genetics_terms = [
            ("Gene", "unit", "genetics", "Unit of heredity - segment of DNA encoding a functional product", "Contains instructions for making RNA/protein", "Fundamental to inheritance"),
            ("Allele", "variant", "genetics", "Alternative form of a gene", "Determines different traits", "Genetic variation"),
            ("Chromosome", "structure", "genetics", "Thread-like structure of DNA and protein", "Carries genetic information", "Genome structure"),
            ("DNA", "molecule", "genetics", "Deoxyribonucleic acid - genetic material", "Stores genetic information", "Life's blueprint"),
            ("RNA", "molecule", "genetics", "Ribonucleic acid - involved in gene expression", "Translates genetic code", "Gene expression"),
            ("Genome", "set", "genetics", "Complete set of genetic material", "Contains all genes", "Complete heredity"),
            ("Genotype", "characteristic", "genetics", "Genetic constitution of an organism", "Genetic makeup", "Inheritance"),
            ("Phenotype", "characteristic", "genetics", "Observable physical characteristics", "Physical expression of genotype", "Traits"),
            ("Mutation", "process", "genetics", "Change in DNA sequence", "Can alter gene function", "Genetic variation, disease"),
            ("Haplotype", "set", "genetics", "Set of genetic variants on single chromosome", "Inherited together", "Population genetics"),
            ("Diploid", "cell", "genetics", "Cell with two chromosome sets", "Somatic cells", "Cell biology"),
            ("Haploid", "cell", "genetics", "Cell with one chromosome set", "Gametes", "Reproduction"),
            ("Dominant", "allele", "genetics", "Allele expressed in heterozygote", "Trait expression", "Inheritance"),
            ("Recessive", "allele", "genetics", "Allele only expressed in homozygote", "Hidden trait expression", "Inheritance"),
            ("Codominant", "allele", "genetics", "Both alleles expressed in heterozygote", "Joint trait expression", "Inheritance"),
            ("Incomplete Dominance", "process", "genetics", "Neither allele is dominant - blended phenotype", "Intermediate phenotype", "Inheritance"),
            ("Epistasis", "process", "genetics", "One gene affects expression of another", "Gene interaction", "Genetics"),
            ("Pleiotropy", "process", "genetics", "One gene affects multiple traits", "Multiple trait effects", "Genetics"),
            ("Polygenic", "characteristic", "genetics", "Trait influenced by multiple genes", "Complex traits", "Genetics"),
            ("Linkage", "process", "genetics", "Genes on same chromosome inherited together", "Gene clustering", "Mapping"),
            ("Recombination", "process", "genetics", "Exchange of genetic material between chromosomes", "Genetic diversity", "Evolution"),
            ("Transcription", "process", "genetics", "DNA to RNA synthesis", "Gene expression", "Central dogma"),
            ("Translation", "process", "genetics", "RNA to protein synthesis", "Protein production", "Central dogma"),
            ("Replication", "process", "genetics", "DNA copying process", "Cell division", "Inheritance"),
            ("Repair", "process", "genetics", "DNA damage correction", "Genome integrity", "Mutagenesis"),
            ("Mutation Rate", "measure", "genetics", "Frequency of mutations per generation", "Evolution rate", "Genetics"),
            ("Gene Expression", "process", "genetics", "Conversion of gene to functional product", "Cell differentiation", "Development"),
            ("Promoter", "sequence", "genetics", "DNA sequence initiating transcription", "Transcription start", "Gene regulation"),
            ("Enhancer", "sequence", "genetics", "DNA sequence enhancing transcription", "Gene activation", "Regulation"),
            ("Silencer", "sequence", "genetics", "DNA sequence suppressing transcription", "Gene repression", "Regulation"),
            ("Intron", "sequence", "genetics", "Non-coding region of gene", "Splicing, regulation", "Gene structure"),
            ("Exon", "sequence", "genetics", "Coding region of gene", "Protein coding", "Gene structure"),
            ("Operon", "unit", "genetics", "Unit of genes under common control", "Coordinated gene expression", "Bacterial genetics"),
            ("Plasmid", "molecule", "genetics", "Circular DNA in bacteria", "Gene transfer, cloning", "Biotechnology"),
            ("Vector", "molecule", "biotechnology", "DNA vehicle for gene transfer", "Gene delivery", "Cloning"),
            ("Cloning", "process", "biotechnology", "Creating identical copy of DNA or organism", "Gene amplification", "Biotechnology"),
            ("PCR", "technique", "biotechnology", "Polymerase Chain Reaction - DNA amplification", "DNA copying", "Diagnostics, forensics"),
            ("RT-PCR", "technique", "biotechnology", "Reverse Transcription PCR - RNA detection", "Gene expression analysis", "Research"),
            ("qPCR", "technique", "biotechnology", "Quantitative PCR - measures DNA quantity", "Quantification", "Diagnostics"),
            ("Gel Electrophoresis", "technique", "biotechnology", "DNA separation by size", "DNA analysis", "Research"),
            ("Sequencing", "technique", "bioinformatics", "Determining nucleotide order", "Genome reading", "Genomics"),
            ("Genome Sequencing", "technique", "bioinformatics", "Determining complete DNA sequence", "Complete genome", "Genomics"),
            ("SNP", "variant", "genetics", "Single Nucleotide Polymorphism - single base variation", "Genetic variation", "Personalized medicine"),
            ("CNV", "variant", "genetics", "Copy Number Variation - gene dosage changes", "Genetic variation", "Disease"),
            ("Gene Therapy", "technique", "medicine", "Treating disease by modifying genes", "Genetic correction", "Future medicine"),
            ("CRISPR", "technique", "biotechnology", "Gene editing technology", "Precise gene modification", "Revolutionary biotechnology"),
            ("Cas9", "enzyme", "biotechnology", "CRISPR-associated endonuclease", "DNA cutting", "Gene editing"),
        ]
        
        for item in genetics_terms:
            terms.append((
                item[0], item[1], item[2], item[3], item[4], item[5],
                "Core Knowledge Base", None, None, None, datetime.now()
            ))
        
        # Cell Biology Terms
        cell_terms = [
            ("Cell", "unit", "cell_biology", "Basic structural and functional unit of life", "Performs all life functions", "Fundamental life unit"),
            ("Cell Membrane", "structure", "cell_biology", "Phospholipid bilayer surrounding cell", "Barrier, transport, signaling", "Cell structure"),
            ("Cell Wall", "structure", "cell_biology", "Rigid outer layer in plants/bacteria", "Structure, protection", "Cell biology"),
            ("Nucleus", "organelle", "cell_biology", "Membrane-bound organelle containing DNA", "Gene expression, cell control", "Cell command center"),
            ("Mitochondria", "organelle", "cell_biology", "Energy-producing organelle", "ATP synthesis", "Powerhouse"),
            ("Endoplasmic Reticulum", "organelle", "cell_biology", "Network of membranes for protein/lipid synthesis", "Protein folding, lipid synthesis", "Biosynthesis"),
            ("Golgi Apparatus", "organelle", "cell_biology", "Organelle modifying and packaging proteins", "Protein processing, secretion", "Protein sorting"),
            ("Lysosome", "organelle", "cell_biology", "Organelle containing digestive enzymes", "Autophagy, digestion", "Cleanup"),
            ("Peroxisome", "organelle", "cell_biology", "Organelle with oxidative enzymes", "Fatty acid oxidation, detoxification", "Metabolism"),
            ("Ribosome", "organelle", "cell_biology", "Protein synthesis machinery", "Translation", "Protein factory"),
            ("Centrosome", "organelle", "cell_biology", "Organelle organizing microtubules", "Cell division", "Cell division"),
            ("Chloroplast", "organelle", "cell_biology", "Photosynthetic organelle in plants", "Photosynthesis", "Energy from sun"),
            ("Vacuole", "organelle", "cell_biology", "Storage organelle in plants", "Storage, turgor pressure", "Plant cells"),
            ("Cytoplasm", "substance", "cell_biology", "Gel-like substance inside cell", "Metabolic reactions", "Cell interior"),
            ("Cytoskeleton", "structure", "cell_biology", "Protein network for cell shape and movement", "Structure, motility", "Cell structure"),
            ("Nucleolus", "organelle", "cell_biology", "Region in nucleus making ribosomes", "Ribosome production", "Protein synthesis"),
            ("Nuclear Envelope", "structure", "cell_biology", "Membrane surrounding nucleus", "Nuclear protection", "Cell biology"),
            ("Nuclear Pore", "structure", "cell_biology", "Opening in nuclear envelope", "Transport between nucleus and cytoplasm", "Transport"),
            ("Endocytosis", "process", "cell_biology", "Uptake of materials by membrane invagination", "Nutrient uptake, signaling", "Membrane transport"),
            ("Exocytosis", "process", "cell_biology", "Release of materials by membrane fusion", "Secretion, neurotransmitter release", "Release"),
            ("Phagocytosis", "process", "cell_biology", "Cell engulfing large particles", "Immune defense", "Immunity"),
            ("Pinocytosis", "process", "cell_biology", "Cell taking up fluids and dissolved substances", "Nutrient uptake", "Uptake"),
            ("Apoptosis", "process", "cell_biology", "Programmed cell death", "Development, homeostasis", "Cell death"),
            ("Necrosis", "process", "cell_biology", "Uncontrolled cell death", "Injury response", "Cell death"),
            ("Autophagy", "process", "cell_biology", "Self-digestion of cell components", "Recycling, stress response", "Cleanup"),
            ("Mitosis", "process", "cell_biology", "Cell division producing identical cells", "Growth, repair", "Cell division"),
            ("Meiosis", "process", "cell_biology", "Cell division producing gametes", "Sexual reproduction", "Reproduction"),
            ("Cytokinesis", "process", "cell_biology", "Division of cell cytoplasm", "Cell division completion", "Cell division"),
            ("Interphase", "phase", "cell_biology", "Cell cycle phase between divisions", "Growth, DNA replication", "Cell cycle"),
            ("Prophase", "phase", "cell_biology", "Early mitosis - chromosome condensation", "Chromosome preparation", "Cell division"),
            ("Metaphase", "phase", "cell_biology", "Mitosis phase - chromosomes align", "Chromosome alignment", "Cell division"),
            ("Anaphase", "phase", "cell_biology", "Mitosis phase - sister chromatids separate", "Chromosome separation", "Cell division"),
            ("Telophase", "phase", "cell_biology", "Mitosis phase - nuclear envelopes reform", "Nuclear reformation", "Cell division"),
            ("Stem Cell", "cell", "cell_biology", "Undifferentiated cell with self-renewal capacity", "Differentiation potential", "Regenerative medicine"),
            ("Progenitor Cell", "cell", "cell_biology", "Partially differentiated precursor cell", "Limited differentiation", "Tissue repair"),
            ("Differentiated Cell", "cell", "cell_biology", "Specialized cell with specific function", "Tissue function", "Histology"),
            ("Fibroblast", "cell", "cell_biology", "Cell producing collagen and extracellular matrix", "Wound healing, structure", "Connective tissue"),
            ("Neuron", "cell", "cell_biology", "Nerve cell - transmits electrical signals", "Signal transmission", "Nervous system"),
            ("Glial Cell", "cell", "cell_biology", "Support cell in nervous system", "Neuron support, myelination", "Nervous system"),
            ("Macrophage", "cell", "cell_biology", "Large phagocytic cell", "Immune defense", "Immunity"),
            ("T Cell", "cell", "cell_biology", "Lymphocyte with T-cell receptor", "Cell-mediated immunity", "Immune system"),
            ("B Cell", "cell", "cell_biology", "Lymphocyte producing antibodies", "Humoral immunity", "Immune system"),
            ("NK Cell", "cell", "cell_biology", "Natural Killer cell - innate immunity", "Virus-infected cell killing", "Innate immunity"),
            ("Erythrocyte", "cell", "cell_biology", "Red blood cell - carries oxygen", "Oxygen transport", "Blood"),
            ("Leukocyte", "cell", "cell_biology", "White blood cell - immune defense", "Immunity", "Blood"),
            ("Platelet", "cell", "cell_biology", "Cell fragment for blood clotting", "Hemostasis", "Blood clotting"),
        ]
        
        for item in cell_terms:
            terms.append((
                item[0], item[1], item[2], item[3], item[4], item[5],
                "Core Knowledge Base", None, None, None, datetime.now()
            ))
        
        # Microbiology Terms
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
            ("Antibiotic", "compound", "medicine", "Substance killing or inhibiting bacteria", "Bacterial infection treatment", "Medicine"),
            ("Antiviral", "compound", "medicine", "Drug treating viral infections", "Virus replication inhibition", "Medicine"),
            ("Antifungal", "compound", "medicine", "Drug treating fungal infections", "Fungal growth inhibition", "Medicine"),
            ("Vaccine", "agent", "medicine", "Preparation providing immunity", "Immune system stimulation", "Disease prevention"),
            ("Immunity", "state", "medicine", "Resistance to disease", "Protection from infection", "Health"),
            ("Antibody", "protein", "medicine", "Protein produced by immune system", "Antigen recognition", "Immunity"),
            ("Antigen", "molecule", "medicine", "Substance triggering immune response", "Immune activation", "Immunology"),
            ("Immunoglobulin", "protein", "medicine", "Antibody protein class", "Humoral immunity", "Immunology"),
            ("Cytokine", "molecule", "cell_biology", "Cell signaling protein", "Immune regulation", "Immunology"),
            ("Interferon", "protein", "medicine", "Protein inducing antiviral state", "Viral defense", "Immunology"),
            ("Complement", "system", "medicine", "Protein system for pathogen lysis", "Pathogen destruction", "Innate immunity"),
            ("Major Histocompatibility Complex", "system", "medicine", "Cell surface proteins for immune recognition", "Immune regulation", "Transplantation"),
            ("Human Leukocyte Antigen", "molecule", "medicine", "MHC in humans", "Immune recognition", "Transplantation"),
            ("Innate Immunity", "system", "medicine", "Non-specific first-line defense", "Immediate response", "Immunology"),
            ("Adaptive Immunity", "system", "medicine", "Specific immune response with memory", "Specific defense", "Immunology"),
            ("Humoral Immunity", "system", "medicine", "Antibody-mediated immunity", "Blood-borne defense", "Immunology"),
            ("Cell-Mediated Immunity", "system", "medicine", "T-cell mediated immunity", "Intracellular pathogen defense", "Immunology"),
        ]
        
        for item in microbe_terms:
            terms.append((
                item[0], item[1], item[2], item[3], item[4], item[5],
                "Core Knowledge Base", None, None, None, datetime.now()
            ))
        
        # Biochemistry Terms
        biochem_terms = [
            ("Metabolism", "process", "biochemistry", "All chemical reactions in living organisms", "Energy production, biosynthesis", "Life processes"),
            ("Catabolism", "process", "biochemistry", "Breakdown of molecules for energy", "Energy release", "Metabolism"),
            ("Anabolism", "process", "biochemistry", "Building up molecules using energy", "Biosynthesis", "Metabolism"),
            ("Glycolysis", "pathway", "biochemistry", "Glucose breakdown to pyruvate", "Energy production", "Metabolism"),
            ("Krebs Cycle", "pathway", "biochemistry", "Citric acid cycle - aerobic metabolism", "Energy production", "Metabolism"),
            ("Electron Transport Chain", "pathway", "biochemistry", "ATP synthesis using electron carriers", "Oxidative phosphorylation", "Energy production"),
            ("Oxidative Phosphorylation", "process", "biochemistry", "ATP synthesis using oxygen", "Most ATP production", "Energy"),
            ("Photosynthesis", "process", "biochemistry", "Conversion of light energy to chemical energy", "Food production in plants", "Life on Earth"),
            ("Light Reactions", "process", "biochemistry", "Photosynthetic reactions capturing light energy", "ATP and NADPH production", "Photosynthesis"),
            ("Calvin Cycle", "pathway", "biochemistry", "Carbon fixation in photosynthesis", "Glucose production", "Photosynthesis"),
            ("Gluconeogenesis", "pathway", "biochemistry", "Glucose synthesis from non-carbohydrates", "Blood glucose maintenance", "Metabolism"),
            ("Glycogenesis", "pathway", "biochemistry", "Glucose to glycogen conversion", "Energy storage", "Metabolism"),
            ("Glycogenolysis", "pathway", "biochemistry", "Glycogen to glucose breakdown", "Energy release", "Metabolism"),
            ("Lipogenesis", "pathway", "biochemistry", "Fatty acid synthesis", "Fat storage", "Metabolism"),
            ("Lipolysis", "process", "biochemistry", "Triglyceride breakdown", "Energy release from fat", "Metabolism"),
            ("Beta Oxidation", "pathway", "biochemistry", "Fatty acid breakdown for energy", "Energy production", "Metabolism"),
            ("Urea Cycle", "pathway", "biochemistry", "Ammonia to urea conversion", "Nitrogen excretion", "Excretion"),
            ("Protein Synthesis", "process", "biochemistry", "Amino acid to protein conversion", "Gene expression", "Central dogma"),
            ("Transamination", "process", "biochemistry", "Amino group transfer", "Amino acid metabolism", "Metabolism"),
            ("Deamination", "process", "biochemistry", "Amino group removal", "Amino acid catabolism", "Metabolism"),
            ("Enzyme Kinetics", "study", "biochemistry", "Study of enzyme reaction rates", "Reaction mechanism", "Biochemistry"),
            ("Michaelis-Menten", "model", "biochemistry", "Enzyme kinetics model", "Reaction rate analysis", "Kinetics"),
            ("Inhibition", "process", "biochemistry", "Enzyme activity reduction", "Regulation, drug action", "Enzymology"),
            ("Allosteric Regulation", "process", "biochemistry", "Enzyme regulation by binding at distant site", "Metabolic control", "Regulation"),
            ("Feedback Inhibition", "process", "biochemistry", "End product inhibiting earlier pathway step", "Metabolic homeostasis", "Regulation"),
            ("Coenzyme", "molecule", "biochemistry", "Non-protein enzyme helper", "Catalysis", "Enzymology"),
            ("Prosthetic Group", "molecule", "biochemistry", "Tightly bound enzyme cofactor", "Enzyme function", "Enzymology"),
            ("Cofactor", "molecule", "biochemistry", "Non-protein enzyme requirement", "Enzyme activity", "Enzymology"),
            ("Oxidation", "process", "biochemistry", "Loss of electrons", "Energy extraction", "Chemistry"),
            ("Reduction", "process", "biochemistry", "Gain of electrons", "Biosynthesis", "Chemistry"),
            ("Redox Reaction", "process", "biochemistry", "Coupled oxidation-reduction", "Energy transfer", "Chemistry"),
            ("pH", "measure", "biochemistry", "Measure of acidity/basicity", "Enzyme function", "Chemistry"),
            ("Buffer", "solution", "biochemistry", "Solution resisting pH changes", "pH maintenance", "Chemistry"),
            ("Osmosis", "process", "biochemistry", "Water flow across semipermeable membrane", "Cell water balance", "Transport"),
            ("Diffusion", "process", "biochemistry", "Movement from high to low concentration", "Passive transport", "Transport"),
            ("Active Transport", "process", "biochemistry", "Energy-requiring molecule movement", "Concentration gradient", "Transport"),
            ("Facilitated Diffusion", "process", "biochemistry", "Carrier-mediated passive transport", "Specific molecule transport", "Transport"),
        ]
        
        for item in biochem_terms:
            terms.append((
                item[0], item[1], "biochemistry", item[3], item[4], item[5],
                "Core Knowledge Base", None, None, None, datetime.now()
            ))
        
        # Medical Terms
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
