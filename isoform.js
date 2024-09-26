const urlParams = new URLSearchParams(window.location.search);
const mutation = urlParams.get('mutation'); // Replace 'mutation' with the appropriate parameter
const gene = urlParams.get('gene');
const mutation_path=mutation.replace(/\s/g, '')  // Replace spaces
                            .replace(/\(/g, '_') // Replace "(" with "_"
                            .replace(/\)/g, '')  // Remove ")"
                            .replace(/:/g, '_')  // Replace ":" with "_"
                            .replace(/>/g, '-'); // Replace ">" with "-"


document.getElementById('gene-name').innerHTML = `<a href="gene.html?gene=${gene}">${gene}</a>`;
document.getElementById('gene-mutation').innerHTML = `${mutation}`;
document.getElementById('download-p-sequence').href = `output/${gene}/${mutation_path}/seq.fa`;

fetch(`output/${gene}/gene_info.json`)
    .then(response => response.json())
    .then(data => {
        if (data && data[gene]) {
            document.getElementById('gene-description').innerHTML = `${data[gene]}`;
        }
    })
    .catch(error => {
        console.error('Error fetching the gene info:', error);

    })



const pdbPath = `output/${gene}/${mutation_path}/rank1_relaxed.pdb`; // Path to the PDB file
const plddtPath = `output/${gene}/${mutation_path}/rank1_relaxed.json`; // Path to the plddt file
const viewerContainer = document.getElementById('3d-viewer');
const viewer = $3Dmol.createViewer(viewerContainer, {
    defaultcolors: $3Dmol.rasmolElementColors,
});

let pdbUri = `output/${gene}/${mutation_path}/rank1_relaxed.pdb`;  // PDB file path
let plddtUri = `output/${gene}/${mutation_path}/plddt.json`;      // PLDDT file path (or another file)

// Use jQuery.ajax to load the PDB file
jQuery.ajax(pdbUri, {
    success: function(data) {
        let v = viewer;
        v.addModel(data, "pdb");                      // Load the PDB data
        v.setStyle({}, { cartoon: { color: 'spectrum' } }); // Apply cartoon style
        v.zoomTo();                                   // Adjust camera zoom
        v.render();                                   // Render the structure
        v.zoom(1.2, 1000);                            // Slight zoom

        // Adjust canvas position (remove absolute positioning)
        const viewerDiv = document.getElementById('3d-viewer');
        const canvasElement = viewerDiv.querySelector('canvas');
        if (canvasElement) {
            canvasElement.style.position = 'relative';  // Reset position to relative
        }

        // Set download links for the PDB and other files
        document.getElementById('download-pdb').href = pdbUri;   // Link for PDB download
        document.getElementById('download-plddt').href = plddtUri; // Link for PLDDT (or another file)
    },
    error: function(hdr, status, err) {
        console.error("Failed to load PDB " + pdbUri + ": " + err);
        viewerContainer.innerHTML = '<p>Error loading 3D structure.</p>';
    }
});

    // Define the path to the score.txt file
const scoreFilePath = `output/${gene}/${mutation_path}/score.txt`;

// Get the TMscore container
const TMscoreDiv = document.getElementById('TMscore');


const pdbPath1 = `output/${gene}/canonical/rank1_relaxed.pdb`;  // Canonical protein
const pdbPath2 = pdbPath; // Mutant protein (existing mutant PDB path)

// Get the viewer container
const compareContainer = document.getElementById('3d-compare');

// Initialize the 3Dmol.js viewer
const viewerCompare = $3Dmol.createViewer(compareContainer, {
    defaultcolors: $3Dmol.rasmolElementColors,
    backgroundColor: 'white' // Optional: Set background color
});

let protein1, protein2; // Variables to store the protein models
function renameChain(pdbData, oldChain, newChain) {
    // PDB format has chain identifiers on atom lines (ATOM and HETATM) at position 22
    const regex = new RegExp(`^(.{21})${oldChain}(.{4}.*)$`, 'gm');
    return pdbData.replace(regex, `$1${newChain}$2`);
}

// Load the canonical PDB file (Protein 1) using jQuery.ajax
jQuery.ajax( pdbPath1, { 
  success: function(data) {
    let v = viewerCompare;
    v.addModel( data, "pdb" );                       /* load data */
    v.setStyle({ chain: 'A' }, {cartoon: {color: 'blue'}});  /* style all atoms */
    v.zoomTo();                                      /* set camera */
    v.render();                                      /* render scene */
    v.zoom(1.2, 1000);                               /* slight zoom */
    const viewerDiv = document.getElementById('3d-compare');
    const canvasElement = viewerDiv.querySelector('canvas');    if (canvasElement) {
        canvasElement.style.position = 'relative';  // Reset position to relative
    }
  },
  error: function(hdr, status, err) {
    console.error( "Failed to load PDB " + pdbUri + ": " + err );
  },
});


jQuery.ajax( pdbPath2, { 
    success: function(data) {
      let modifiedPdbData = renameChain(data, 'A', 'B');

      let v = viewerCompare;
      v.addModel( modifiedPdbData, "pdb" );                       /* load data */
      v.setStyle({ chain: 'B' }, { cartoon: { color: 'red' } }); // Color mutant protein red
      v.zoomTo();                                      /* set camera */
      v.render();                                      /* render scene */
      v.zoom(1.2, 1000);                               /* slight zoom */
      const viewerDiv = document.getElementById('3d-compare');
      const canvasElement = viewerDiv.querySelector('canvas');    if (canvasElement) {
          canvasElement.style.position = 'relative';  // Reset position to relative
      }
    },
    error: function(hdr, status, err) {
      console.error( "Failed to load PDB " + pdbUri + ": " + err );
    },
  });
  const aminoAcidMap = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V",
    "Sec": "U", "Pyl": "O"  // Include special cases like Selenocysteine (Sec) and Pyrrolysine (Pyl)
};

// Populate the dropdown with mutations from the CSV file using PapaParse
jQuery.get("output/INS/clinvar_seqMUT_scores.csv", function(csvData) {
    const select = document.getElementById('mutant-select');
    // Use PapaParse to reliably parse the CSV data
    Papa.parse(csvData, {
        header: true, // Assumes the first row contains column headers
        complete: function(results) {
            results.data.forEach((row) => {
                const mutation = row['data']; // Assuming 'data' is the name of the first column
                if (mutation && mutation.trim() !== "") {
                    // Add the mutation to the dropdown
                    let option = document.createElement('option');
                    option.value = mutation.trim(); // Ensure there's no extra space in the value
                    option.text = mutation.trim();  // Trim any spaces for cleaner display
                    select.appendChild(option);
                };
                let mutName=mutation.replace(/\s/g, '')  // Replace spaces
                                    .replace(/\(/g, '_') // Replace "(" with "_"
                                    .replace(/\)/g, '')  // Remove ")"
                                    .replace(/:/g, '_')  // Replace ":" with "_"
                                    .replace(/>/g, '-'); // Replace ">" with "-"
                if (mutName === mutation_path) {
                    let dnaSequence = row['WTcSeq'].toUpperCase();
                    let proteinSequence = row['MUTpSeq'];
                    let pathogenicity = row['ClinicalSignificance'].split('|');
                    let phenotypes = row['PhenotypeList'].split('|');
                    let phenotypesCodes = row['PhenotypeIDS'].split('|');

                    let list = document.getElementById('pathogenicity');
                    pathogenicity.forEach(item => {
                        const listItem = document.createElement('li'); // Create a new list item
                        listItem.textContent = item; // Set the text content
                        list.appendChild(listItem); // Append to the li
                    });
                    
                    databseDic={
                        "MONDO":"https://www.ebi.ac.uk/ols4/ontologies/mondo/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252F",
                        "MedGen":"https://www.ncbi.nlm.nih.gov/medgen/",
                        "OMIM":"https://omim.org/entry/",
                    };
                    let list2 = document.getElementById('Phenotypes');
                    phenotypes.forEach((item,index) => {
                        let PhenotypeIDs = phenotypesCodes[index];
                        idString=""
                        
                        PhenotypeIDs.split(",").forEach(id => {
                            id = id.replace("MONDO:MONDO:","MONDO:").split(":");
                            source = id[0];
                            code = id[1];
                            if (source in databseDic) {
                                idString += `<a href="${databseDic[source]}${code}">  ${source}:${code}</a>`;
                            } else {
                                idString += `${source}:${code} `;
                                console.log(source,code);
                            }

                            
                        }

                        )

                        const listItem = document.createElement('li'); // Create a new list item
                        listItem.innerHTML  = `${item} (${idString})`; // Set the text content
                        list2.appendChild(listItem); // Append to the list
                    });
                    


                    const dnaMatch = mutation.match(/c\.(\d+)([A-Z])>([A-Z])/); // Matches c.266G>A
                    const proteinMatch = mutation.match(/p\.(\w{3})(\d+)(\w{3})/);  // Matches p.Arg89His
                    if (dnaMatch) {
                        let position = parseInt(dnaMatch[1], 10);   // Get the position (266)
                        let original = dnaMatch[2];                 // Get the original nucleotide (G)
                        let mutated = dnaMatch[3];                  // Get the mutated nucleotide (A)
                
                        // Replace the nucleotide at the given position with a span containing mutation info
                        let highlightedDNA = dnaSequence.slice(0, position - 1) +
                            `<span class="highlight" data-mutation="${position} ${original}>${mutated}">${mutated}</span>` +
                            dnaSequence.slice(position);
                
                        document.getElementById('c-sequence').innerHTML = highlightedDNA;
                    }
                
                    if (proteinMatch) {
                        let position = parseInt(proteinMatch[2], 10);  // Get the position (89)
                        let original = proteinMatch[1];                // Get the original amino acid (Arg)
                        let mutated = proteinMatch[3];                 // Get the mutated amino acid (His)
                
                        // Replace the amino acid at the given position with a span containing mutation info
                        let highlightedProtein = proteinSequence.slice(0, position - 1) +
                            `<span class="highlight" data-mutation="${position} ${original}>${mutated}">${aminoAcidMap[mutated]}</span>` +
                            proteinSequence.slice(position);
                
                        document.getElementById('p-sequence').innerHTML = highlightedProtein;
                    }
                }
            });
        }
    });
});

// Handle mutation selection from the dropdown
document.getElementById('mutant-select').addEventListener('change', function(event) {
    const selectedMutation = event.target.value;
    if (selectedMutation) {
        // Update the path for the selected mutation's PDB
        loadMutantPDB(selectedMutation); // Load the new mutation PDB
        updateTMScore(selectedMutation); // Update the TMscore card with the selected mutation
    }
});


// Function to load PDB for both canonical and mutant proteins
function loadMutantPDB(mutationName) {
    // Clear the 3D viewer (remove all models)
    viewerCompare.removeAllModels();
    let newName = mutationName.replace(/\s/g, '')  // Replace spaces
                                .replace(/\(/g, '_') // Replace "(" with "_"
                                .replace(/\)/g, '')  // Remove ")"
                                .replace(/:/g, '_')  // Replace ":" with "_"
                                .replace(/>/g, '-'); // Replace ">" with "-"
    let newMutationPath = `output/${gene}/${newName}/rank1_relaxed.pdb`; // Path to the PDB file
    // Load the canonical protein first
    jQuery.ajax(newMutationPath, {
        success: function(data) {
            let protein1 = viewerCompare.addModel(data, "pdb"); // Load canonical protein
            viewerCompare.setStyle({ chain: 'A' }, { cartoon: { color: 'blue' } }); // Color canonical protein blue
            viewerCompare.zoomTo(); // Adjust zoom
            viewerCompare.render(); // Render the viewer

            // Load the mutant protein after the canonical protein is loaded
            jQuery.ajax( pdbPath2, { 
                success: function(data) {
                  let modifiedPdbData = renameChain(data, 'A', 'B');
            
                  let v = viewerCompare;
                  v.addModel( modifiedPdbData, "pdb" );                       /* load data */
                  v.setStyle({ chain: 'B' }, { cartoon: { color: 'red' } }); // Color mutant protein red
                  v.zoomTo();                                      /* set camera */
                  v.render();                                      /* render scene */
                  v.zoom(1.2, 1000);                               /* slight zoom */
                  const viewerDiv = document.getElementById('3d-compare');
                  const canvasElement = viewerDiv.querySelector('canvas');    if (canvasElement) {
                      canvasElement.style.position = 'relative';  // Reset position to relative
                  }
                },
                error: function(hdr, status, err) {
                  console.error( "Failed to load PDB " + pdbUri + ": " + err );
                },
              });;
        },
        error: function(hdr, status, err) {
            console.error("Failed to load canonical PDB file " + pdbPath1 + ": " + err);
        }
    });
}

Papa.parse(`output/${gene}/tmalign_network.csv`, {
    download: true,
    header: true, // Treat the first row as the header
    complete: function(results) {
        // Store the parsed data
        tmalignData = results.data;
        updateTMScore("canonical");

    },
    error: function(err) {
        console.error("Error parsing TMalign CSV:", err);
    }
});
// Function to update the TMscore card based on the selected mutation
function updateTMScore(selectedMutation) {
    const tmscoreDiv = document.getElementById('TMscore');
    tmscoreDiv.innerHTML = ''; // Clear previous content
    const selectedMutationNew = selectedMutation.replace(/\s/g, '')  // Replace spaces
                                                .replace(/\(/g, '_') // Replace "(" with "_"
                                                .replace(/\)/g, '')  // Remove ")"
                                                .replace(/:/g, '_')  // Replace ":" with "_"
                                                .replace(/>/g, '-'); // Replace ">" with "-"

    // Loop through the TMalign data to find the matching rows for the canonical and selected mutation
    const match = tmalignData.find(row => {
        // Ensure Protein1 and Protein2 are defined before checking with .includes()
        const protein1 = row.Protein1 ? row.Protein1 : '';
        const protein2 = row.Protein2 ? row.Protein2 : '';
 
        return (
            (protein1 == selectedMutationNew && protein2 == mutation_path)
        );
    });
    // If a match is found, update the TMscore card
    if (match) {
        tmscoreDiv.innerHTML = `<pre>${match.TMalign_Output}</pre>`;
    } else {
        tmscoreDiv.innerHTML = 'No TM-align result available for the selected mutation.';
    }
}

// Event listener for mutation selection
document.getElementById('mutant-select').addEventListener('change', function(event) {
    const selectedMutation = event.target.value;
    if (selectedMutation) {
        updateTMScore(selectedMutation); // Update the TMscore card with the selected mutation
    }
});

