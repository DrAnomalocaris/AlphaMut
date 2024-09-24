const urlParams = new URLSearchParams(window.location.search);
const mutation = urlParams.get('mutation'); // Replace 'mutation' with the appropriate parameter
const gene = urlParams.get('gene');
const mutation_path=mutation.replace(/\s/g, '')  // Replace spaces
                            .replace(/\(/g, '_') // Replace "(" with "_"
                            .replace(/\)/g, '')  // Remove ")"
                            .replace(/:/g, '_')  // Replace ":" with "_"
                            .replace(/>/g, '-'); // Replace ">" with "-"


const titleContainer = document.getElementById('title');
titleContainer.innerHTML = `
                        <h1>${gene}</h1>
                        <h2>${mutation}</h2>
                    `;

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

// Fetch and load the contents of score.txt
fetch(scoreFilePath)
    .then(response => response.text())
    .then(scoreText => {
        // Insert the file contents into the TMscore div
        TMscoreDiv.innerHTML = `<pre>${scoreText}</pre>`; // Use <pre> to preserve formatting
    })
    .catch(error => {
        console.error('Error loading score.txt:', error);
        TMscoreDiv.innerHTML = '<p>Error loading TMscore.</p>';
    });


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
  
  