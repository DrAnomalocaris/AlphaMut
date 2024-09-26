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
rule find_and_rotate_best_folding:
    input:
        "output/{gene}/{transcript}/isoform.done.txt",
        "output/{gene}/canonical/isoform.done.txt"

    output:
        pdb="output/{gene}/{transcript}/rank1_relaxed.pdb",
        json="output/{gene}/{transcript}/rank1_relaxed.json",

    run:
        from Bio.PDB import PDBParser, Superimposer, PDBIO, Select
        from glob import glob
        pdbMutant    = glob(f"output/{wildcards.gene}/{wildcards.transcript}/isoform_relaxed_rank_001_alphafold2_ptm_model_*_seed_*.pdb")[0]
        pdbCanonical = glob(f"output/{wildcards.gene}/canonical/isoform_relaxed_rank_001_alphafold2_ptm_model_*_seed_*.pdb")[0]
        if pdbMutant == pdbCanonical: 
            shell(f"cp output/{wildcards.gene}/{wildcards.transcript}/isoform_relaxed_rank_001_alphafold2_ptm_model_*_seed_*.pdb {output.pdb};")   
        shell(f"cp output/{wildcards.gene}/{wildcards.transcript}/isoform_scores_rank_001_alphafold2_ptm_model_*_seed_*.json {output.json};")

        parser = PDBParser(QUIET=True)
        structure1 = parser.get_structure('canonical',pdbCanonical)  # Canonical protein PDB
        structure2 = parser.get_structure('mutant',pdbMutant)        # Mutant protein PDB
        # Select alpha-carbons (CA atoms) for alignment from both structures
        atoms1 = []
        atoms2 = []

        for model1, model2 in zip(structure1, structure2):
            for chain1, chain2 in zip(model1, model2):
                for residue1, residue2 in zip(chain1, chain2):
                    # Select alpha carbon atoms
                    if 'CA' in residue1 and 'CA' in residue2:
                        atoms1.append(residue1['CA'])
                        atoms2.append(residue2['CA'])
        

        # Initialize the superimposer
        super_imposer = Superimposer()

        # Perform the superimposition
        super_imposer.set_atoms(atoms1, atoms2)

        # Print RMSD
        print(f"RMSD after superimposition: {super_imposer.rms:.4f} Å")

        # Apply the transformation (rotate and translate mutant to align with canonical)
        super_imposer.apply(structure2.get_atoms())

        # Save the aligned mutant structure to a new PDB file
        io = PDBIO()
        io.set_structure(structure2)
        io.save(output.pdb)

          
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

rule simplify_table:
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


rule autofillList:
    output:
        "autofill.js"
    run:
        import pickle
        from glob import glob
        with open(output[0], 'w') as f:
            f.write(f"const autofill = [\n")
            for i in ([i.split("/")[1] for i in glob("output/*/isoform_plot.html")]):
                f.write(f"    '{i}',\n")
            f.write("];\n")





rule dot_plot:
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


rule network_calculation:
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
        pdb_files  = glob(f"output/{wildcards.gene}/*/rank1_relaxed.pdb")


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
            customdata=df['url'],  # Use the 'url' column from the dataframe for redirection
            hovertemplate=(
                "%{hovertext}<br>"  # Display the hover_text column
                #"Clinical Significance: %{x}<br>"  # Show the x-axis value (Clinical Significance)
                #f"{metric.capitalize()}: %{y}<br>"  # Show the y-axis value (RMSD or other metric)
                "<extra></extra>"  # Remove default extra info
            )
        )


        fig.write_html(output[0], include_plotlyjs='cdn', full_html=False, div_id='plotly_graph')
        with open(output[0], 'a') as f:
            f.write('''<script>
            var plot = document.querySelectorAll('.js-plotly-plot')[0];  // Get the first Plotly plot
            plot.on('plotly_click', function(data) {
                var point = data.points[0];
                var url = point.customdata;  // customdata contains the URL
                window.open(url, "_blank");  // Opens the link in a new tab
            });
            </script>
            ''')

rule get_gene_indo:
    output:
        "output/{gene}/gene_info.json"
    run:


        import time
        import xmltodict
        from Bio import Entrez

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
                time.sleep(0.34)  # Requests to NCBI are rate limited to 3 per second

            return gene_summaries

        email = "x"
        gene_summaries = get_entrez_gene_summary("INS", email)
        with open(f"output/{wildcards.gene}/gene_info.json", "w") as f:
            json.dump(dict(gene_summaries), f)



rule process_gene:
    input:
        "output/{gene}/isoform_plot.html",
        "output/{gene}/clinvar_seqMUT_scores.csv",
        "output/{gene}/clinvar_seqMUT_scores_summary.html",
        "output/{gene}/dot_plot.rmsd.png",
        "output/{gene}/dot_plot.identical_v_aligned.png",
        "output/{gene}/gene_info.json",
        "output/{gene}/tmalign_network.csv"

    output:
        "output/{gene}/DONE.txt"
    shell:
        "touch {output}"

rule all:
    input:
        expand("output/{gene}/DONE.txt",
            gene=[
                "INS",
                "GCG",
                "ADM",
                "SST"]
        )


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

        pdb_files = glob(f"output/{wildcards.gene}/*/rank1_relaxed.pdb")

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
