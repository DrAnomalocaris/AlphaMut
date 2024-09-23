# Define all chromosomes to download
CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]

# Define the output directory
OUTPUT_DIR = "data/genome/human/"

rule python_path:
    shell:
        "which python"
rule download_variation_file_from_ClinVar:
    output:
        "data/clinvar.txt"
    shell:
        """
        wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
        gunzip -c variant_summary.txt.gz > {output}
        """

rule preprocess_variation_file:
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
    input:
        #expand(f"{OUTPUT_DIR}{{chrom}}.fa", chrom=CHROMOSOMES),
        #gtf = "data/Homo_sapiens.GRCh38.109.gtf.pkl",
        table="output/{gene}/clinvar_raw.csv"

    output:
        table="output/{gene}/clinvar_seqWT.csv"
    threads: 1
    run:

        import requests
        import pandas as pd
        import xmltodict

        clinvar = pd.read_csv(input.table)
        print(clinvar)
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
                seqs[data] = seq2 

            elif row.Type == "Deletion":
                location=data.split(":")[1].split("c.")[-1].split("del")[-0].split("_")
                if len(location)==1: location = [location[0],location[0]]
                start,end = [int(x)-1 for x in location]
                cSeq = row.WTcSeq
                cSeq = cSeq[:start]+cSeq[end+1:]
                cSeq = Seq(cSeq)
                pSeq = cSeq.translate(to_stop=True)
                mutations[data] = (mutation)
                seqs[data] = seq2 
            elif row.Type == "Microsatellite" or row.Type == "Duplication":
                if data.endswith("del)"):
                    if mutation.endswith("del"):
                        mutation = mutation[:-3]
                        originalAA=mutation[0]
                        location = int(mutation[1:])-1
                        pSeq = (seq2[:location-1]+seq2[location:])
                        mutations[data] = (mutation)
                        seqs[data] = seq2 
                        continue     
                location=data.split(":")[1].split("c.")[-1].split(" ")[0].replace("dup","").split("_")
                if len(location)==1: location = [location[0],location[0]]
                start,end = [int(x)-1 for x in location]
                cSeq = row.WTcSeq
                cSeq = cSeq[:start]+ cSeq[start:end]*2+cSeq[end+1:]
                cSeq = Seq(cSeq)
                pSeq = cSeq.translate(to_stop=True)
                mutations[data] = (mutation)
                seqs[data] = seq2 

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

rule set_output_folders:
    # This needs to be altered to get mutated variants
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
    input:
        "output/{gene}/isoforms.fa"
    output:
        "output/{gene}/{transcript}/seq.fa"
    run:
        from Bio import SeqIO
        
        # Extract the specific transcript ID from wildcards
        transcript_id = wildcards.transcript
        
        # Read the input fasta file and search for the transcript
        with open(input[0], "r") as fasta_in:
            for record in SeqIO.parse(fasta_in, "fasta"):
                if record.id == transcript_id:
                    # Write the specific transcript to the output file
                    with open(output[0], "w") as fasta_out:
                        record.id = "isoform"
                        record.description = "isoform"

                        SeqIO.write(record , fasta_out, "fasta")
                    break
            else:
                raise ValueError(f"Transcript {transcript_id} not found in {input[0]}")
rule run_colabfold:
    input:
        "output/{gene}/{transcript}/seq.fa"
    output:
        "output/{gene}/{transcript}/isoform_plddt.png",
        "output/{gene}/{transcript}/isoform.done.txt",
    shell:
        """
        /mnt/c/Users/lahat/colabfold/localcolabfold/colabfold-conda/bin/colabfold_batch --amber "{input[0]}"  "output/{wildcards.gene}/{wildcards.transcript}/" 
        """


rule run_all_isoforms:
    input:
        fasta = "output/{gene}/isoforms.fa"

    output:
        "output/{gene}/done.txt",

    run:
        import pickle
        gene = wildcards.gene
        transcripts=[]
        with open(input.fasta, 'r') as f: 
            for line in f:
                if line.startswith(">"):
                    
                    transcript = line.split()[0][1:]
                    transcripts.append(transcript)
        to_do = []
        for transcript in transcripts:
            to_do.append( f"output/{gene}/{transcript}/score.json",)
        to_do = " ".join(sorted(to_do))
        shell(
            f"snakemake {to_do} -c {threads} --rerun-incomplete -p"
            )
        shell(
            f"touch output/{gene}/done.txt"
            )
rule find_best_folding:
    input:
        "output/{gene}/{transcript}/isoform.done.txt"

    output:
        pdb="output/{gene}/{transcript}/rank1_relaxed.pdb",
        json="output/{gene}/{transcript}/rank1_relaxed.json",

    shell:
        """
        cp output/{wildcards.gene}/{wildcards.transcript}/isoform_relaxed_rank_001_alphafold2_ptm_model_*_seed_*.pdb {output.pdb};
        cp output/{wildcards.gene}/{wildcards.transcript}/isoform_scores_rank_001_alphafold2_ptm_model_*_seed_*.json {output.json};
        """
rule compare_folding:
    input:
        canonical = "output/{gene}/canonical/rank1_relaxed.pdb",
        mutant = "output/{gene}/{mut}/rank1_relaxed.pdb"

    output:
        "output/{gene}/{mut}/score.txt",

    shell:
        """
        TMalign {input.canonical} {input.mutant} > {output}
        """

rule parse_tmscore:
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

        out["chain1"] = f"output/{wildcards.gene}/canonical/rank1_relaxed.pdb"
        out["chain2"] = f"output/{wildcards.gene}/{wildcards.mut}/rank1_relaxed.pdb"
        out["chain1_length"] = int(re.search(r'Length of Chain_1:\s+(\d+)\s+residues',chains).group(1))
        out["chain2_length"] = int(re.search(r'Length of Chain_2:\s+(\d+)\s+residues',chains).group(1))

        search = (re.search(r"Aligned length=\s*([\d.]+),\s*RMSD=\s*([\d.]+),\s*Seq_ID=.*?=\s*([\d.]+)",scores))
        out["aligned_length"] = float(search.group(1))
        out["rmsd"] = float(search.group(2))
        out["identical/aligned"] = float(search.group(3))

        with open(output.json, 'w') as f:
            json.dump(out, f, indent=4)


rule add_scores_to_table:
    input:
        "output/{gene}/done.txt",
        table="output/{gene}/clinvar_seqMUT.csv"   

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
            if not os.path.exists(scores_json): 
                shell(f"snakemake {scores_json} -c {threads} --rerun-incomplete -p")

            with open(scores_json, 'r') as f:
                score = json.load(f)
            newTable[i]=score
        newTable = pd.DataFrame(newTable).T
        df = df.join(newTable)
        df.to_csv(output.table)


rule autofillList:
    input:
        "data/Homo_sapiens.GRCh38.109.gtf.pkl"
        
    output:
        "autofill.js"
    run:
        import pickle
        from glob import glob
        with open(input[0], 'rb') as f: 
            general, isoforms = pickle.load(f)
        with open(output[0], 'w') as f:
            f.write(f"const autofill = [\n")
            for i in ([i.split("/")[-1] for i in glob("output/*")]):
                f.write(f"    '{i}',\n")
                name = general[i]['name']
                f.write(f"'{name}',\n")
            f.write("];\n")
            f.write("const geneDict = {\n")
            for i in ([i.split("/")[-1] for i in glob("output/*")]):
                name = general[i]['name']
                f.write(f"    '{name}':'{i}',\n")
                f.write(f"    '{i}':'{i}',\n")
            f.write("};\n")
rule makeJsonForGene:
    input:
        "data/Homo_sapiens.GRCh38.109.gtf.pkl",
        "output/{gene}/done.txt",
        gtf = "data/Homo_sapiens.GRCh38.109.gtf",
        
    output:
        "output/{gene}/data.js"
    run:
        import pickle
        from glob import glob
        import json
        from gtfparse import read_gtf
        import polars as pl
        import os

        # Load the general and isoforms data from the pickle file
        with open(input[0], 'rb') as f:
            general, isoforms = pickle.load(f)

        # Get the data for the specific gene
        data = general[wildcards.gene]
        transcripts = data["transcripts"]

        # Initialize dictionaries for storing PDB content and scores
        d = {}
        s = {}

        for transcript in transcripts:
            # Read and store the content of the PDB file
            if not os.path.exists(f"output/{wildcards.gene}/{transcript}"):continue
            folded_path = glob(f"output/{wildcards.gene}/{transcript}/isoform_relaxed_rank_001_alphafold2_ptm_model_*_seed_*.pdb")[0]
            with open(folded_path, 'r') as pdb_file:
                pdb_content = pdb_file.read()
            d[transcript] = pdb_content  # Store the PDB content directly

            # Store the path to the scores file (or its content if needed)
            scores_path = glob(f"output/{wildcards.gene}/{transcript}/isoform_scores_rank_001_alphafold2_ptm_model_*_seed_*.json")[0]
            s[transcript] = scores_path  # Store the path to the scores file

        # Create the JavaScript object with PDB content and scores
        json_object = "const relaxed_pdb = " + json.dumps(d, indent=4) 
        json_object += ";\nconst scores = " + json.dumps(s, indent=4)
        gtf =  read_gtf(input["gtf"])
        # Filter for transcripts of the specific gene
        gtf = gtf.filter(
            (pl.col("gene_biotype") == "protein_coding") & 
            (pl.col("feature")      == "transcript") & 
            (pl.col("gene_id")      == wildcards.gene)
        )
        transcript_dict = dict(zip(gtf['transcript_id'].to_list(), gtf['transcript_name'].to_list()))
        json_object += ";\nconst transcript_names = " + json.dumps(transcript_dict, indent=4)

        # Write the JavaScript object to the output file
        with open(output[0], 'w') as f:
            f.write(json_object)

genes = [
    "ENSG00000141510", #TP53
    "ENSG00000254647", #insulin
    "ENSG00000170315", #Ubiquitin UBB
    "ENSG00000197061", #Histone H4, HIST1H4A
    "ENSG00000164825", #Beta-Defensin1 DEFB1
    "ENSG00000154620", #thymosin beta 4 Y-linked 
    "ENSG00000164128", #Neuropeptide Y
    "ENSG00000101200", #Vasopressin
    ]

rule all:
    input:
        expand("output/{gene}/data.js", gene=genes),
    shell:
        "rm autofill.js ; "
        "snakemake autofill.js -c1 --rerun-incomplete -p"

rule dot_plot:
    input:
        table="output/{gene}/clinvar_seqMUT_scores.csv"
    output:
        "output/{gene}/dot_plot.{metric}.png"
    run:
        import seaborn as sns
        import pandas as pd
        import matplotlib.pyplot as plt
        df = pd.read_csv(input.table)

        # Define the correct order for ClinicalSignificance
        order = [
            'Likely benign',
            'Uncertain significance',
            'Uncertain significance/Uncertain risk allele',
            'Conflicting classifications of pathogenicity',
            'not provided',
            'Likely risk allele',
            'Likely pathogenic/Likely risk allele',
            'Pathogenic/Likely risk allele',
            'Likely pathogenic',
            'Pathogenic/Likely pathogenic/Likely risk allele',
            'Pathogenic/Likely pathogenic',
            'Pathogenic',
        ]

        # Convert ClinicalSignificance to a categorical type with the correct order
        df['ClinicalSignificance'] = pd.Categorical(df['ClinicalSignificance'], categories=order, ordered=True)

        # Create the figure and axis
        plt.figure(figsize=(10, 6))

        # Add the boxplot to show distribution and mean
        sns.boxplot(
            y=wildcards.metric, 
            x="ClinicalSignificance", 
            data=df, 
            order=order,
            showcaps=True,  # Show caps on boxplot whiskers
            boxprops={'facecolor': 'None'},  # Transparent boxplot to show points behind it
            medianprops={'color': 'red'},  # Red line for the median
            showmeans=True,  # Show the mean
            #meanprops={"marker": "o", "markerfacecolor": "white", "markeredgecolor": "black", "markersize": "8"},  # Customize mean marker
        )

        # Overlay the stripplot (scatter with jitter)
        sns.stripplot(
            y=wildcards.metric, 
            x="ClinicalSignificance", 
            data=df, 
            jitter=True,  # Add jitter for horizontal wiggle
            dodge=True,   # Spread overlapping points
            alpha=0.6,    # Make points slightly transparent for better visibility
            linewidth=1,  # Add border to points
            order=order
        )

        # Add vertical lines to separate categories
        for i in range(1, len(order)):
            plt.axvline(i - 0.5, color='gray', linestyle='--', alpha=0.6)

        # Rotate the x-axis labels vertically
        plt.xticks(rotation=90)

        # Add labels and title
        plt.xlabel('Clinical Significance')
        plt.ylabel('RMSD')
        plt.title(f"{wildcards.gene} - {wildcards.metric} vs. Clinical Significance")

        # Save the figure
        plt.savefig(output[0])


