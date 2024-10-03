import yaml

# Load the configuration file
configfile: "config.yaml"

# Access the variables from the config
GENES = config["genes"]
NUM_FOLDINGS = config["num_foldings"]
colabfoldExec = config["colabfold_exec"]

rule python_path:
    # this helps when debugging
    shell:
        "which python"
rule plot_DAG:
    #plot the dag
    input:
        "Snakefile"
    output:
        "dag.png"
    shell:
        "snakemake --rulegraph all | dot -Gdpi=300 -Tpng > {output}"


rule download_variation_file_from_ClinVar:
    #gets clinical variation data from https://www.ncbi.nlm.nih.gov/clinvar/
    output:
        "data/clinvar.txt"
    shell:
        """
        wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
        gunzip -c variant_summary.txt.gz > {output}
        """

rule preprocess_variation_file:
    #preprocess the clinical variation data into a dictionary with gene as key so future retrival is faster
    input:
        "data/clinvar.txt"
    output:
        "data/clinvar.pk"
    run:
        from pprint import pprint
        import pickle
        from tqdm import tqdm
        out={}
        num_lines = sum(1 for _ in open(input[0]))
        with open(input[0], 'r') as f:
            header = f.readline()
            header = header[1:].strip().split("\t")
            for line in tqdm(f, total=num_lines):
                line = dict(zip(header, line.strip().split("\t")))
                gene = line["GeneSymbol"]
                data = line['Name']
                ClinicalSignificance = line['ClinicalSignificance']
                PhenotypeIDS = line['PhenotypeIDS']
                PhenotypeList = line['PhenotypeList']
                Type = line['Type'] 
                HGNC_ID = line['HGNC_ID']
                if not gene in out:
                    out[gene] = []
                out[gene].append(
                    {
                        "data": data,
                        "ClinicalSignificance": ClinicalSignificance,
                        "PhenotypeIDS": PhenotypeIDS,
                        "PhenotypeList": PhenotypeList,
                        "Type": Type,
                        "HGNC_ID": HGNC_ID
                    }
                )

        with open(output[0], 'wb') as f:
            pickle.dump(out, f)
            

rule set_gene_specific_clinVar_file:
    #sets the gene specific clinvar file
    input:
        table="data/clinvar.pk"
    output:
        table="output/{gene}/clinvar_raw.csv"
    run:    
        import pandas as pd
        import pickle
        with open(input.table, 'rb') as f:
            clinvar = pickle.load(f)[wildcards.gene]
        pd.DataFrame(clinvar).to_csv(output.table, index=False)
rule add_protein_WT_sequences_to_raw_clinVar:
    #adds wildtype sequence to raw clinvar using entrez efetch
    input:
        table="output/{gene}/clinvar_raw.csv"

    output:
        table="output/{gene}/clinvar_seqWT.csv"
    threads: 1
    run:

        import requests
        import pandas as pd
        import xmltodict

        clinvar = pd.read_csv(input.table)
        seqs = []
        cSeqs = []
        import functools

        @functools.cache
        def fetch(gene):
            if gene.startswith("NC_"): return ("", "")
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            # Define the parameters for the query
            params = {
                "db": "nuccore",
                "id": gene,   # The RefSeq accession for which you want data
                "rettype": "gb",       # "gb" for GenBank format, could be changed to "fasta" or "xml"
                "retmode": "xml"      # Return as plain text
            }
            # Make the request
            response = requests.get(url, params=params)
            extra = xmltodict.parse(response.text)
            from pprint import pprint
            GBfeatures =  extra["GBSet"]["GBSeq"]["GBSeq_feature-table"]["GBFeature"]
            cSeq = extra["GBSet"]["GBSeq"]['GBSeq_sequence']
            for feature in GBfeatures:
                if feature["GBFeature_key"] != "CDS": continue
                for qual in feature["GBFeature_quals"]['GBQualifier']:
                    if not qual["GBQualifier_name"] == "translation": continue
                    pSeq = (qual['GBQualifier_value'])
                    return pSeq, cSeq
            return "", cSeq
            
        from tqdm import tqdm
        for i, row in tqdm(clinvar.iterrows(), total=len(clinvar)):
            geneCode = row["data"].split("(")[0]
            pSeq, cSeq = fetch(geneCode)

            seqs.append(pSeq)
            cSeqs.append(cSeq)
        clinvar["WTpSeq"] = seqs
        clinvar["WTcSeq"] = cSeqs
        clinvar.to_csv(output.table, index=False)
        
rule add_mutation_sequence_to_raw_clinVar:
    #parses the mutation and adds it to the sequence
    input:
        table="output/{gene}/clinvar_seqWT.csv"

    output:
        table="output/{gene}/clinvar_seqMUT.csv"
    run:

        import pandas as pd 
        from tqdm import tqdm
        import re
        from Bio.Seq import Seq

        aa_map = {
                    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C', 'Glu': 'E',
                    'Gln': 'Q', 'Gly': 'G', 'His': 'H', 'Ile': 'I', 'Leu': 'L', 'Lys': 'K',
                    'Met': 'M', 'Phe': 'F', 'Pro': 'P', 'Ser': 'S', 'Thr': 'T', 'Trp': 'W',
                    'Tyr': 'Y', 'Val': 'V',"Ter":'.'
                }
        clinvar = pd.read_csv(input.table, index_col=0)

        seqs = {}
        mutations = {}
        for i, row in tqdm(clinvar.iterrows(), total=len(clinvar)):
            data = i
            seq = row.WTpSeq
            seq2 = seq
            mutation = data.split(":")[1]
            if len(mutation.split(" "))==1:continue
            mutation = mutation.split(" ")[1].replace("(","").replace(")","").replace("p.","")
            #replace tripple amino acids name with single
            for aaa,a in aa_map.items():
                mutation = mutation.replace(aaa,a)


            if row.Type == "single nucleotide variant":
                if mutation.endswith("="): continue
                originalAA=mutation[0]
                newAA=mutation[-1]
                position = int(mutation[1:-1])-1
                seq2 = list(seq2)
                seq2[position] = newAA
                seq2 = "".join(seq2)
                mutations[data] = (mutation)
                seqs[data] = seq2
            elif row.Type == "Indel":
                location=data.split(":")[1].split("c.")[-1].split("delins")[-0].split("_")
                new = data.split(":")[1].split("c.")[-1].split("delins")[1].split()[0]
                if len(location)==1: location = [location[0],location[0]]
                start,end = [int(x)-1 for x in location]
                cSeq = row.WTcSeq
                cSeq = cSeq[:start]+new+cSeq[end+1:]
                cSeq = Seq(cSeq)
                pSeq = cSeq.translate(to_stop=True)
                mutations[data] = (mutation)
                seqs[data] = pSeq 

            elif row.Type == "Deletion":
                location=data.split(":")[1].split("c.")[-1].split("del")[-0].split("_")
                if len(location)==1: location = [location[0],location[0]]
                start,end = [int(x)-1 for x in location]
                cSeq = row.WTcSeq
                cSeq = cSeq[:start]+cSeq[end+1:]
                cSeq = Seq(cSeq)
                pSeq = cSeq.translate(to_stop=True)
                mutations[data] = (mutation)
                seqs[data] = pSeq 
            elif row.Type == "Microsatellite" or row.Type == "Duplication":
                if data.endswith("del)"):
                    if mutation.endswith("del"):
                        mutation = mutation[:-3]
                        originalAA=mutation[0]
                        location = int(mutation[1:])-1
                        pSeq = (seq2[:location-1]+seq2[location:])
                        mutations[data] = (mutation)
                        seqs[data] = pSeq 
                        continue     
                if "del (" in data:
                        cMutation = data.split(":")[1].split("c.")[-1].split("del")[-0].split("_")
                        if len(cMutation)==1: cMutStart,cMutEnd = [int(cMutation[0])-1,int(cMutation[0])-1]
                        else: cMutStart,cMutEnd = [int(x)-1 for x in cMutation]
                        cSeq = row.WTcSeq
                        cSeq = cSeq[:cMutStart]+cSeq[cMutEnd+1:]
                        cSeq = Seq(cSeq)
                        pSeq = cSeq.translate(to_stop=True)
                        mutations[data] = (mutation)
                        seqs[data] = pSeq
                        continue
                location=data.split(":")[1].split("c.")[-1].split(" ")[0].replace("dup","").split("_")
                if len(location)==1: location = [location[0],location[0]]
                start,end = [int(x)-1 for x in location]
                cSeq = row.WTcSeq
                cSeq = cSeq[:start]+ cSeq[start:end]*2+cSeq[end+1:]
                cSeq = Seq(cSeq)
                pSeq = cSeq.translate(to_stop=True)
                mutations[data] = (mutation)
                seqs[data] = pSeq 

            else:
                print(row.Type)
                print(data)
                break
            
            mutations[data] = mutations[data].split(".")[0]
        seqs = pd.Series(seqs)
        mutations = pd.Series(mutations)

        clinvar["Mutation"] = mutations
        clinvar["MUTpSeq"] = seqs

        #filter out mutations that did not change the protein sequence
        clinvar = clinvar[clinvar.MUTpSeq != clinvar.WTpSeq]
        #drop rows with MUTpSeq as nan
        clinvar = clinvar.dropna()
        #remove duplicate rows
        clinvar = clinvar.drop_duplicates()
        clinvar.to_csv(output.table)

checkpoint set_output_folders:
    # Creates a folder for a mutation output and includes an isoform fasta file with the mutations
    input:
        table="output/{gene}/clinvar_seqMUT.csv"   
    output:
        path="output/{gene}/isoforms.fa"
    threads: 1
    run:
        import pandas as pd
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord

        clinvar = pd.read_csv(input.table, index_col=0)
        records = []

        record = SeqRecord(
            Seq(clinvar['WTpSeq'].iloc[0]),           # The sequence itself
            id="canonical",            # The sequence ID
            description=""           # Optional description
            )
        records.append(record)
    

        for i, row in clinvar.iterrows():
            transcript = i.replace(" ","").replace("(","_").replace(")","").replace(":","_").replace(">","-")
            protein_seq = row["MUTpSeq"].split(".")[0]
        

            record = SeqRecord(
                Seq(protein_seq),           # The sequence itself
                id=transcript,            # The sequence ID
                description=""           # Optional description
                )
            if len(protein_seq)>1:
                records.append(record)
        

        with open(output.path, "w") as handle:
            SeqIO.write(records, handle, "fasta")

rule prepare_isoform_for_alphafold:
    #Extracts the specific transcript from the isoform fasta file
    input:
        "output/{gene}/isoforms.fa"
    output:
        "output/{gene}/{mutation}/seq.fa"
    run:
        from Bio import SeqIO
        
        # Extract the specific transcript ID from wildcards
        mutation_id = wildcards.mutation
        
        # Read the input fasta file and search for the transcript
        with open(input[0], "r") as fasta_in:
            for record in SeqIO.parse(fasta_in, "fasta"):
                if record.id == mutation_id:
                    # Write the specific transcript to the output file
                    with open(output[0], "w") as fasta_out:
                        record.id = "isoform"
                        record.description = "isoform"

                        SeqIO.write(record , fasta_out, "fasta")
                    break
            else:
                raise ValueError(f"Transcript {mutation_id} not found in {input[0]}")

def get_mutations(gene):
    # returns a list of mutations from the mutations fasta file
    mutations_file = checkpoints.set_output_folders.get(gene=gene).output[0]
    mutations=[]
    with open(mutations_file, 'r') as f: 
        for line in f:
            if line.startswith(">"):
                
                mutation = line.split()[0][1:]
                mutations.append(mutation)
    return mutations
    


rule run_colabfold:
    # Runs colabfold on the isoform
    input:
        fasta = "output/{gene}/{mutation}/seq.fa"
    output:
        relaxed   = "output/{gene}/{mutation}/isoform_relaxed_rank_001_alphafold2_ptm_model_1_seed_{seed}.pdb",
        unrelaxed = "output/{gene}/{mutation}/isoform_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_{seed}.pdb",
        scores    = "output/{gene}/{mutation}/isoform_scores_rank_001_alphafold2_ptm_model_1_seed_{seed}.json",
    log:
        stdout = "output/{gene}/{mutation}/isoform_scores_rank_001_alphafold2_ptm_model_1_seed_{seed}.log",
        stderr = "output/{gene}/{mutation}/isoform_scores_rank_001_alphafold2_ptm_model_1_seed_{seed}.err",
    shell:
        """
        {colabfoldExec} \
            --amber     \
            --num-seeds 1 \
            --random-seed {wildcards.seed} \
            --overwrite-existing-results \
            --num-models 1 \
            "{input.fasta}"  \
            "output/{wildcards.gene}/{wildcards.mutation}/" > {log.stdout} 2> {log.stderr}; \
        rm -rf \
            output/{wildcards.gene}/{wildcards.mutation}/isoform_env \
            output/{wildcards.gene}/{wildcards.mutation}/config.json \
            output/{wildcards.gene}/{wildcards.mutation}/isoform_coverage.png \
            output/{wildcards.gene}/{wildcards.mutation}/isoform_pae.png \
            output/{wildcards.gene}/{wildcards.mutation}/isoform_plddt.png \
            output/{wildcards.gene}/{wildcards.mutation}/isoform_predicted_aligned_error_v1.json \
            output/{wildcards.gene}/{wildcards.mutation}/isoform.a3m \
            output/{wildcards.gene}/{wildcards.mutation}/isoform.done.txt \
            output/{wildcards.gene}/{wildcards.mutation}/log.txt 
        """


        
rule rotate_foldings_to_match_first:
    # Rotates the foldings to match the first protein (ie seed=000)
    input:
        protein     = "output/{gene}/{mutation}/isoform_relaxed_rank_001_alphafold2_ptm_model_1_seed_{seed}.pdb",
        reference   = "output/{gene}/{mutation}/isoform_relaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb",

    output: 
        "output/{gene}/{mutation}/{seed}/isoform_relaxed_rotated.pdb",

    run:
        from biopandas.pdb import PandasPdb
        from Bio.PDB import Superimposer, PDBParser, PDBIO
        from glob import glob
        from tqdm import tqdm
        import numpy as np
        # Load the proteins
        protein_file  = input.protein
        reference_file = input.reference
        if protein_file == reference_file:
            shell(f"cp {reference_file} {output[0]}")
        else:
            # Biopython PDB parser
            parser = PDBParser(QUIET=True)

            # Select a reference structure (use the first protein as reference)
            reference_structure = parser.get_structure("reference", reference_file)

            # Extract the CA atoms from the reference structure
            ref_atoms = [atom for atom in reference_structure.get_atoms() if atom.id == 'CA']

            # Initialize Superimposer from Biopython
            sup = Superimposer()

            # Iterate over the proteins and superimpose each on the reference
            target_structure = parser.get_structure("target", protein_file)
            
            # Extract the CA atoms from the target structure
            target_atoms = [atom for atom in target_structure.get_atoms() if atom.id == 'CA']
            
            # Ensure both structures have the same number of CA atoms
            if len(ref_atoms) == len(target_atoms):
                # Superimpose the target on the reference
                sup.set_atoms(ref_atoms, target_atoms)
                sup.apply(target_structure.get_atoms())  # Apply rotation and translation to all atoms
                
                # Save the aligned structure
                io = PDBIO()
                io.set_structure(target_structure)
                io.save(output[0])
            else:
                print(f"Skipping {protein_file} due to mismatched atom counts.")

rule average_foldings:
    # Averages the foldings
    input:
        lambda wc : expand("output/{gene}/{mutation}/{seed:03d}/isoform_relaxed_rotated.pdb", 
                            gene=[wc.gene], 
                            mutation=[wc.mutation], 
                            seed=range(NUM_FOLDINGS)),
    output:
        "output/{gene}/{mutation}/isoform_relaxed_averaged_rotated.pdb",
    run:
        from Bio.PDB import PDBParser, Superimposer, PDBIO, Select
        from biopandas.pdb import PandasPdb

        import pandas as pd
        from tqdm import tqdm
        pdbs = [PandasPdb().read_pdb(i) for i in tqdm(input)]#
        new = pdbs[0].df["ATOM"].copy()

        for col in ["x_coord", "y_coord", "z_coord","b_factor"]:
            new[col] = sum([i.df['ATOM'][col] for i in pdbs])/len(pdbs)
        pdb = PandasPdb().read_pdb(input[0])
        pdb.df["ATOM"] = new
        pdb.to_pdb(output[0]) 




rule run_all_isoforms:
    # Uses checkpoints to reasses the DAG, and run all isoforms
    input:
        lambda wc : expand("output/{gene}/{mutation}/{file}", 
            gene=[wc.gene], 
            mutation=get_mutations(wc.gene),
            n=[NUM_FOLDINGS],
            file=[
                "isoform_relaxed_averaged_rotated.pdb",
                "plddt.csv",
                "plddt.html"
                ]
            )

    output:
        "output/{gene}/done.txt",

    shell:
        "touch {output}"

          
rule compare_folding:
    # Compares the folding using TMalign
    input:
        canonical = "output/{gene}/canonical/isoform_relaxed_averaged_rotated.pdb",
        mutant = "output/{gene}/{mut}/isoform_relaxed_averaged_rotated.pdb"

    output:
        "output/{gene}/{mut}/score.txt",

    shell:
        """
        TMalign {input.canonical} {input.mutant} > {output}
        """

rule parse_tmscore:
    # Parses the TMscore into JSON
    input:
        table="output/{gene}/{mut}/score.txt"
    output:
        json="output/{gene}/{mut}/score.json"
    run:
        import re
        import json
        with open(input.table, 'r') as f:
            txt = f.read()
        x = txt.split("\n\n")
        out = {}
        chains = x[1]
        scores = x[2]
        comparison = x[3]
        out["comparison"] = comparison

        out["chain1"] = f"output/{wildcards.gene}/canonical/isoform_relaxed_averaged_rotated.pdb"
        out["chain2"] = f"output/{wildcards.gene}/{wildcards.mut}/isoform_relaxed_averaged_rotated.pdb"
        out["chain1_length"] = int(re.search(r'Length of Chain_1:\s+(\d+)\s+residues',chains).group(1))
        out["chain2_length"] = int(re.search(r'Length of Chain_2:\s+(\d+)\s+residues',chains).group(1))

        search = (re.search(r"Aligned length=\s*([\d.]+),\s*RMSD=\s*([\d.]+),\s*Seq_ID=.*?=\s*([\d.]+)",scores))
        out["aligned_length"] = float(search.group(1))
        out["rmsd"] = float(search.group(2))
        out["identical/aligned"] = float(search.group(3))

        with open(output.json, 'w') as f:
            json.dump(out, f, indent=4)



rule add_scores_to_table:
    # Adds the TMalign scores to the table
    input:
        "output/{gene}/done.txt",
        table="output/{gene}/clinvar_seqMUT.csv" ,
        scores=lambda wc : expand("output/{gene}/{transcript}/score.json",gene=[wc.gene], transcript=get_mutations(wc.gene))

    output:
        table="output/{gene}/clinvar_seqMUT_scores.csv"

    run:
        import pandas as pd
        df = pd.read_csv(input.table, index_col=0)
        #transcript = i.replace(" ","").replace("(","_").replace(")","").replace(":","_").replace(">","-")
        newTable = {}
        for i in df.index:
            transcript = i.replace(" ","").replace("(","_").replace(")","").replace(":","_").replace(">","-")
            scores_json=f"output/{wildcards.gene}/{transcript}/score.json"

            with open(scores_json, 'r') as f:
                score = json.load(f)
            newTable[i]=score
        newTable = pd.DataFrame(newTable).T
        df = df.join(newTable)
        df.to_csv(output.table)

rule simplify_table:
    # Turns the table into html
    input:
        "output/{gene}/clinvar_seqMUT_scores.csv",
    output:
        "output/{gene}/clinvar_seqMUT_scores_summary.html",
    run:
        import pandas as pd
        df = pd.read_csv(input[0])
        df = df[['data','ClinicalSignificance', 'PhenotypeList','aligned_length','rmsd','identical/aligned']]
        df.columns = ['mutation','clinical significance', 'phenotype list', 'aligned length', 'rmsd', 'identical/aligned']
        df["mutation"] = df["mutation"].apply( # insert links
            lambda x: f"<a href='isoform.html?gene={wildcards.gene}&mutation={x}'>"+x.split(":")[1]+"</a>"
            )

        df.to_html(
            output[0],
            render_links=True,
            escape=False,
            index=False
        )
        #df.to_csv(output[0])


rule get_biomart:
    # gets the biomaRt table
    output:
        table = "data/hsapiens_gene_ensembl.csv"
    run:
        from pybiomart import Dataset

        dataset = Dataset(name='hsapiens_gene_ensembl',
                        host='http://www.ensembl.org')
        data =dataset.query(attributes=['ensembl_gene_id', 'external_gene_name',"description"])
        data.to_csv(output.table)


            

rule dot_plot:
    # Makes a dot plot comparing tmalign scores with clinical significance
    input:
        table="output/{gene}/clinvar_seqMUT_scores.csv"
    output:
        plot="output/{gene}/dot_plot.{metric}.png",
        table="output/{gene}/clinvar_seqMUT_scores.split_categories.{metric}.csv",
    run:
        import seaborn as sns
        import pandas as pd
        import matplotlib.pyplot as plt
        df = pd.read_csv(input.table)
        metric = wildcards.metric
        if metric == "identical_v_aligned":
            metric = "identical/aligned"

        # Define the correct order for ClinicalSignificance
        order = [
            'Benign',
            'Likely benign',
            'Uncertain significance',
            'Uncertain risk allele',
            'Conflicting classifications of pathogenicity',
            'not provided',
            'Likely risk allele',
            'Likely pathogenic',
            'Pathogenic',
        ]

        # Convert ClinicalSignificance to a categorical type with the correct order
        d2={}
        for i, row in df.iterrows():
            for significance in row["ClinicalSignificance"].split("/"):
                row2 = row.copy()
                row2["ClinicalSignificance"] = significance
                d2[i]=row2

        df = pd.DataFrame(d2).T
        df['ClinicalSignificance'] = pd.Categorical(df['ClinicalSignificance'], categories=[i  for i in order if i in df["ClinicalSignificance"].unique()], ordered=True)
        df = df.iloc[::-1]

        df.to_csv(output.table)
        # Create the figure and axis
        plt.figure(figsize=(10, 6))

        # Overlay the stripplot (scatter with jitter)
        sns.stripplot(
            y=metric, 
            x="ClinicalSignificance", 
            data=df, 
            jitter=True,  # Add jitter for horizontal wiggle
            dodge=True,   # Spread overlapping points
            alpha=0.6,    # Make points slightly transparent for better visibility
            linewidth=1,  # Add border to points
            hue="PhenotypeList"

        )

        # Add vertical lines to separate categories
        for i in range(1, len(df.ClinicalSignificance.unique())):
            plt.axvline(i - 0.5, color='gray', linestyle='--', alpha=0.6)

        # Rotate the x-axis labels vertically
        plt.xticks(rotation=90)

        # Add labels and title
        plt.xlabel('Clinical Significance')
        plt.ylabel(metric)
        plt.legend([],[], frameon=False)
        plt.title(f"{wildcards.gene} - {metric} vs. Clinical Significance")

        # Save the figure
        plt.savefig(output.plot,bbox_inches='tight')

rule plddt_plot:
    # Makes a plotly plddt plot
    input:
        jsons = lambda wc: expand("output/{gene}/{allele}/isoform_scores_rank_001_alphafold2_ptm_model_1_seed_{seed:03d}.json", 
                        gene=[wc.gene],   
                        allele=[wc.allele],
                        seed=range(NUM_FOLDINGS)),
        fasta = "output/{gene}/{allele}/seq.fa"
    output:
        table="output/{gene}/{allele}/plddt.csv",
        plot="output/{gene}/{allele}/plddt.html"
    run:
        import json
        import pandas as pd
        import plotly.graph_objects as go
        from Bio import SeqIO
        from Bio.SeqUtils import seq3  # To convert amino acid to three-letter codes

                
        # Define variables
        gene = wildcards.gene
        allele = wildcards.allele
        metric = "plddt"
        seeds = range(NUM_FOLDINGS)
        df = pd.DataFrame({i: json.loads(open(f"output/{gene}/{allele}/isoform_scores_rank_001_alphafold2_ptm_model_1_seed_{i:03d}.json").read())[metric] for i in seeds})
        df.to_csv(output.table)
        # Calculate the average across all seeds
        df['average'] = df.mean(axis=1)

        # Parse the FASTA file to get the amino acid sequence
        fasta_file = f"output/{gene}/{allele}/seq.fa"
        sequence = ""
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequence = str(record.seq)

        # Ensure the sequence length matches the DataFrame length
        if len(sequence) != len(df):
            raise ValueError("The length of the amino acid sequence and the DataFrame rows do not match!")

        # Convert the single-letter amino acid code to the three-letter code
        three_letter_sequence = [seq3(aa).capitalize() for aa in sequence]  # List of three-letter codes

        # Use regex to extract the mutation position and amino acid change if present
        mutation_match = re.search(r'p\.([A-Za-z]+)(\d+)([A-Za-z]+)', allele)

        mutation_position = None
        mutation_annotation = None
        if mutation_match:
            original_aa, position, mutated_aa = mutation_match.groups()
            mutation_position = int(position)  # Mutation position as an integer
            mutation_annotation = f"{(original_aa)} to {(mutated_aa)}"  # Formatted annotation

        # Create the Plotly plot
        fig = go.Figure()

        # Add a line for each seed (with x position in hover template)
        for i in seeds:
            fig.add_trace(go.Scatter(
                x=df.index + 1,  # Adding 1 to match amino acid positions starting from 1
                y=df[i], 
                mode='lines', 
                name=f'{i}',
                #hovertemplate="%{i}: %{x}<br>Score: %{y}<extra></extra>",  # Display x position and y-value
            ))

        # Add the average line (with amino acid and x position in hover template)
        fig.add_trace(go.Scatter(
            x=df.index + 1,  # Adding 1 to match amino acid positions starting from 1
            y=df['average'], 
            mode='lines', 
            name='Average', 
            line=dict(color='firebrick', width=4),
            hovertemplate="<b>%{text}</b> : %{y}<extra></extra>",  # Custom hover with amino acid and position
            text=three_letter_sequence  # Three-letter amino acid code
        ))

        # Add an annotation arrow for the mutation if applicable
        if mutation_position and (mutation_position-1 in df.index):
            fig.add_annotation(
                x=mutation_position,  # X-coordinate is the mutation position
                y=df['average'].iloc[mutation_position - 1],  # Y-coordinate is the average score at mutation position
                text=mutation_annotation,  # Display the amino acid change
                showarrow=True,
                arrowhead=2,
                ax=0,  # X offset for the annotation arrow
                ay=-40,  # Y offset for the annotation arrow
                #bgcolor="green",
                #font=dict(color="black"),
            )


        # Update layout to adjust hover mode and axes
        fig.update_layout(
            title=f"{metric} scores for {gene} ({allele})",
            xaxis_title="Position",
            yaxis_title=f"{metric.capitalize()} Score",
            legend_title="Seed",
            hovermode="x unified",  # Unified hover mode to show all values at a given x position
        )

        # To export the plot as an embeddable HTML file
        fig.write_html(output.plot)




rule network_calculation:
    #calculates the network by doing TMalign between all PDB pairs
    input:
        table="output/{gene}/clinvar_seqMUT_scores.csv"
    output:
        "output/{gene}/network.csv"
    run:
        from glob import glob
        from  tqdm import tqdm
        import numpy as np
        import itertools
        import pandas as pd
        pdb_files  = glob(f"output/{wildcards.gene}/*/isoform_relaxed_averaged_rotated.pdb")


        def calculate_rmsd(pdb1, pdb2):
            result = subprocess.run(["TMalign", pdb1, pdb2], capture_output=True, text=True)
            # Check if the command ran successfully
            if result.returncode != 0:
                raise RuntimeError(f"TM-align failed: {result.stderr}")
            
            # Parse the output
            output = result.stdout
            
            # Extract RMSD from the output using basic string manipulation
            for line in output.splitlines():
                if "RMSD=" in line:
                    # Split by spaces and '=' to extract the RMSD value
                    parts = line.split("RMSD=")[1].split()[0].replace(",", "")
                    rmsd = float(parts)  # The RMSD value is the 5th element in the split line
                    return rmsd
            
            # If RMSD was not found in the output, raise an error
            raise ValueError("RMSD not found in TM-align output.")

        # Initialize distance matrix
        n = len(pdb_files)
        distance_matrix = np.zeros((n, n))
        out=[]
        # Calculate pairwise RMSD for all PDB pairs
        for i, j in tqdm(itertools.combinations(range(n), 2), total=n*(n-1)//2):
            rmsd = calculate_rmsd(pdb_files[i], pdb_files[j])
            out.append((i, j, rmsd))
            
        # Save the distance matrix
        pd.DataFrame(out, columns=["Protein1","Protein2","RMSD"]).to_csv(output[0], index=False)


rule plotly_isoform_plots:
    #plots isoform dot plots with plotly
    input:
        table= "output/{gene}/clinvar_seqMUT_scores.split_categories.rmsd.csv"
    output:
        "output/{gene}/isoform_plot.html"
    run:
        import pandas as pd
        import plotly.express as px


        df = pd.read_csv(input.table)
        metric = "rmsd"
        if metric == "identical_v_aligned":
            metric = "identical/aligned"


        order = df.ClinicalSignificance.unique()
        # Create hover text that includes useful information (customize as needed)
        df['hover_text'] = df['data'] +"<br>" + df['PhenotypeList'] + '<br>' + df[metric].astype(str) + ' ' + metric
        df['url'] = f"isoform.html?gene={wildcards.gene}&mutation=" + df['data']

        # Create the interactive plot using Plotly Express
        fig = px.strip(df,
                        x="ClinicalSignificance",
                        y=metric,
                        hover_name="hover_text",  # Display hover information
                        #hover_data=["PhenotypeList", metric],  # Other hover data
                        color="PhenotypeList",  # Color points based on the phenotype
                        category_orders={"ClinicalSignificance": order},  # Set category order
                        #stripmode="overlay",  # Overlay scatter points
                        )


        # Add vertical lines to separate categories
        for i in range(1, len(order)):
            fig.add_vline(x=i - 0.5, line_width=1, line_dash="dash", line_color="gray")

        # Customize the layout of the plot
        fig.update_layout(
            xaxis_title='Clinical Significance',
            yaxis_title=metric,
            title=f"{wildcards.gene}- {metric} vs. Clinical Significance",
            xaxis={'categoryorder': 'array', 'categoryarray': order},
            template='plotly_white',
            showlegend=False  # This removes the legend
        )
        fig.layout.xaxis.fixedrange = True
        fig.layout.yaxis.fixedrange = True

        # Customize the hover text with hovertemplate and make the URL clickable
        fig.update_traces(
            marker=dict(size=10),
            #customdata=df['url'],  # Use the 'url' column from the dataframe for redirection
            hovertemplate=(
                "%{hovertext}<br>"  # Display the hover_text column
                #"Clinical Significance: %{x}<br>"  # Show the x-axis value (Clinical Significance)
                #f"{metric.capitalize()}: %{y}<br>"  # Show the y-axis value (RMSD or other metric)
                "<extra></extra>"  # Remove default extra info
            )
        )


        fig.write_html(output[0], include_plotlyjs='cdn', full_html=False, div_id='plotly_graph')
        with open(output[0], 'a') as f:
            f.write(f'''<script>
            var plot = document.querySelectorAll('.js-plotly-plot')[0];  // Get the first Plotly plot
            plot.on('plotly_click', function(data) {{
                var point = data.points[0];
                let mutation = point.hovertext.split("<br>")[0];
                mutation = mutation.replace(/\\s/g, '') 
                                .replace(/\\(/g, '_') 
                                .replace(/\\)/g, '')  
                                .replace(/:/g, '_')  
                                .replace(/>/g, '-'); 
                console.log(mutation);
                window.location.href = `isoform.html?mutation=${{mutation}}&gene={wildcards.gene}`; 
            }});
            </script>
            ''')
rule get_gene_info:
    #Get gene info from NCBI
    output:
        "output/{gene}/gene_info.json"
    run:
        import time
        import xmltodict
        from Bio import Entrez
        from random import randint
        from time import sleep
        sleep(randint(3, 6)) #to space them out if using multiple threads

        def get_entrez_gene_summary(
            gene_name, email, organism="human", max_gene_ids=100
        ):
            """Returns the 'Summary' contents for provided input
            gene from the Entrez Gene database. All gene IDs 
            returned for input gene_name will have their docsum
            summaries 'fetched'.
            
            Args:
                gene_name (string): Official (HGNC) gene name 
                (e.g., 'KAT2A')
                email (string): Required email for making requests
                organism (string, optional): defaults to human. 
                Filters results only to match organism. Set to None
                to return all organism unfiltered.
                max_gene_ids (int, optional): Sets the number of Gene
                ID results to return (absolute max allowed is 10K).
                
            Returns:
                dict: Summaries for all gene IDs associated with 
                gene_name (where: keys → [orgn][gene name],
                            values → gene summary)
            """
            Entrez.email = email

            query = (
                f"{gene_name}[Gene Name]"
                if not organism
                else f"({gene_name}[Gene Name]) AND {organism}[Organism]"
            )
            handle = Entrez.esearch(db="gene", term=query, retmax=max_gene_ids)
            record = Entrez.read(handle)
            handle.close()

            gene_summaries = {}
            gene_ids = record["IdList"]

            print(
                f"{len(gene_ids)} gene IDs returned associated with gene {gene_name}."
            )
            for gene_id in gene_ids:
                print(f"\tRetrieving summary for {gene_id}...")
                handle = Entrez.efetch(db="gene", id=gene_id, rettype="docsum")
                gene_dict = xmltodict.parse(
                    "".join([x.decode(encoding="utf-8") for x in handle.readlines()]),
                    dict_constructor=dict,
                )
                gene_docsum = gene_dict["eSummaryResult"]["DocumentSummarySet"][
                    "DocumentSummary"
                ]
                name = gene_docsum.get("Name")
                summary = gene_docsum.get("Summary")
                gene_summaries[name] = summary
                handle.close()
                time.sleep(1)  # Requests to NCBI are rate limited to 3 per second

            return gene_summaries

        email = "x"
        gene_summaries = get_entrez_gene_summary(wildcards.gene, email)
        with open(f"output/{wildcards.gene}/gene_info.json", "w") as f:
            json.dump(dict(gene_summaries), f)



rule process_gene:
    #Process all isoforms for a given gene
    input:
        "output/{gene}/isoform_plot.html",
        "output/{gene}/clinvar_seqMUT_scores.csv",
        "output/{gene}/clinvar_seqMUT_scores_summary.html",
        "output/{gene}/dot_plot.rmsd.png",
        "output/{gene}/dot_plot.identical_v_aligned.png",
        "output/{gene}/gene_info.json",
        "output/{gene}/tmalign_network.csv",
        "output/{gene}/isoform_plot.html"


    output:
        "output/{gene}/DONE_all.txt"
    shell:
        "touch {output}"



rule network_calculation_tmalign:
    input:
        table="output/{gene}/clinvar_seqMUT_scores.csv"
    output:
        "output/{gene}/tmalign_network.csv"
    run:
        from glob import glob
        from tqdm import tqdm
        import numpy as np
        import itertools
        import pandas as pd
        import subprocess

        pdb_files = glob(f"output/{wildcards.gene}/*/isoform_relaxed_averaged_rotated.pdb")

        def run_tmalign(pdb1, pdb2):
            """
            Run TMalign on two PDB files and return the full output.
            """
            result = subprocess.run(["TMalign", pdb1, pdb2], capture_output=True, text=True)
            
            # Check if the command ran successfully
            if result.returncode != 0:
                raise RuntimeError(f"TM-align failed: {result.stderr}")
            
            # Return the full output of TMalign
            return result.stdout

        # Initialize output list to store results
        out = []
        n=len(pdb_files)
        def processName(name):
            return name.split("/")[-2].replace(" ","").replace("(","_").replace(")","").replace(":","_").replace(">","-")

        # Calculate pairwise TMalign output for all PDB pairs
        for i in tqdm(pdb_files):
            for j in pdb_files:
                tmalign_output = run_tmalign(i, j)
                out.append((processName(i), processName(j), tmalign_output))
        
        # Save the TMalign outputs to a CSV
        pd.DataFrame(out, columns=["Protein1", "Protein2", "TMalign_Output"]).to_csv(output[0], index=False)


rule completedDictionary:
    #create dictionary of completed genes, to be use for the search function
    input:
        DONE = expand("output/{gene}/DONE_all.txt",gene=GENES),
        table = "data/hsapiens_gene_ensembl.csv",

    output:
        "completed.json"
    run:
        import json
        import pandas as pd
        print(input)
        genes = [i.split("/")[1] for i in input.DONE]
        table = pd.read_csv(input.table)
        GENES = {}
        for gene in genes:
            df = table[table['Gene name'] == gene]
            names =[i.split(" [")[0] for i in  df['Gene description'].unique()]
            GENES[gene] = gene
            for name in names:
                GENES[name] = gene
        with open("completed.json", "w") as f:
            json.dump(GENES, f)

rule all:
    input:
        "completed.json",
        expand("output/{gene}/isoform_plot.html",gene=GENES),
        expand("output/{gene}/clinvar_seqMUT_scores.csv",gene=GENES),
        expand("output/{gene}/clinvar_seqMUT_scores_summary.html",gene=GENES),
        expand("output/{gene}/dot_plot.rmsd.png",gene=GENES),
        expand("output/{gene}/dot_plot.identical_v_aligned.png",gene=GENES),
        expand("output/{gene}/gene_info.json",gene=GENES),
        expand("output/{gene}/tmalign_network.csv",gene=GENES),
        expand("output/{gene}/isoform_plot.html",gene=GENES),
        lambda wc: [f"output/{gene}/{allele}/plddt.html" for gene,allele in gene_mutants()],
        lambda wc: [f"output/{gene}/{allele}/plddt.csv" for gene,allele in gene_mutants()],
        "completed.json"


def gene_mutants():
    mutations_files = {gene:checkpoints.set_output_folders.get(gene=gene).output[0] for gene in GENES}
    mutations=[]
    for gene,mutations_file in mutations_files.items():
        
        with open(mutations_file, 'r') as f: 
            for line in f:
                if line.startswith(">"):
                    
                    mutation = line.split()[0][1:]
                    mutations.append((gene,mutation))
    return mutations
