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

// Fetch and load the PDB file
fetch(pdbPath)
    .then(response => response.text())
    .then(pdbContent => {
        viewer.addModel(pdbContent, "pdb"); // Load PDB structure
        viewer.setStyle({}, { cartoon: { color: 'spectrum' } }); // Apply cartoon style
        viewer.zoomTo(); // Zoom to fit the structure
        viewer.render(); // Render the 3D structure

        // Remove the inline style attribute from the canvas element after rendering
        const canvasElement = document.querySelector('canvas#undefined');
        if (canvasElement) {
            canvasElement.removeAttribute('style'); // Remove the style attribute
        }
        document.getElementById('download-pdb').href = pdbPath; // Set the link to the PDB file
        document.getElementById('download-plddt').href = plddtPath; // Set the link to the other file

    })
    .catch(error => {
        console.error('Error loading PDB file:', error);
        viewerContainer.innerHTML = '<p>Error loading 3D structure.</p>';
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