// Extract the ENSG code directly from the URL
const urlParams = new URLSearchParams(window.location.search);
const ensgId = urlParams.get('gene');

const geneInfoDiv = document.getElementById('gene-info');
const loadingElement = document.getElementById('loading');

if (ensgId) {
    // Load the gene-specific data.js file using the ENSG code
    const script = document.createElement('script');
    script.src = `output/${ensgId}/data.js`; // Path to the gene's data.js file
    script.onload = () => {
        // Now the relaxed_pdb and scores objects are available
        if (typeof relaxed_pdb === 'undefined' || typeof scores === 'undefined') {
            // Gene data not found, redirect to 404 page with gene parameter
            window.location.href = `404.html?gene=${encodeURIComponent(ensgId)}`;
        } else {
            displayGeneInformation(ensgId); // Pass the ENSG code directly
        }
    };
    script.onerror = () => {
        // Error loading gene data, redirect to 404 page with gene parameter
        window.location.href = `404.html?gene=${encodeURIComponent(ensgId)}`;
    };
    document.head.appendChild(script);
} else {
    // No gene information available, redirect to 404 page
    window.location.href = '404.html';
}

function displayGeneInformation(ensgId) {
    // Clear loading message
    loadingElement.style.display = 'none';

    // Display basic gene information
    let geneInfoHTML = `
        <div class="card-content">
            <h2>Gene Information for ${ensgId}</h2>
            <p><strong>Ensembl ID:</strong> ${ensgId}</p>
        </div>
    `;

    // Array to keep track of Blob URLs for revocation
    const pdbBlobUrls = [];

    // Add a card for each transcript with its own 3D viewer and image
    if (typeof relaxed_pdb !== 'undefined' && typeof scores !== 'undefined') {
        for (const [transcriptId, pdbContent] of Object.entries(relaxed_pdb)) {
            // Create unique IDs for the viewer container
            const viewerContainerId = `protein-viewer-${transcriptId}`;

            // Construct the paths to the image and files
            const imagePath = `output/${ensgId}/${transcriptId}/isoform_plddt.png`;
            const scoresFilePath = `output/${ensgId}/${transcriptId}/data.json`; // Updated path to data.json

            // Create a Blob URL for the PDB content
            const pdbBlob = new Blob([pdbContent], { type: 'chemical/x-pdb' });
            const pdbBlobUrl = URL.createObjectURL(pdbBlob);
            pdbBlobUrls.push(pdbBlobUrl); // Keep track of Blob URLs to revoke later

            // Add a new card for each transcript
            geneInfoHTML += `
                <div class="card">
                    <h3>${transcript_names[transcriptId]}</h3>
                    <h4>${transcriptId}</h4>
                    <p>
                        <a href="${scoresFilePath}" download>Download Scores</a> |
                        <a href="${pdbBlobUrl}" download="${transcriptId}_model.pdb">Download PDB</a>
                    </p>
                    <div class="card-content">
                        <!-- Display the image on the left -->
                        <div class="image-container">
                            <img src="${imagePath}" alt="Isoform pLDDT Image" style="max-width: 100%; cursor: pointer;" onclick="openImageModal('${imagePath}')">
                        </div>
                        <!-- Display the 3D viewer on the right -->
                        <div id="${viewerContainerId}" class="viewer-container" style="width: 100%; height: 400px; position: relative;"></div>
                    </div>
                </div>
            `;
        }

        // Revoke Blob URLs when the page unloads to prevent memory leaks
        window.addEventListener('beforeunload', function() {
            pdbBlobUrls.forEach(function(url) {
                URL.revokeObjectURL(url);
            });
        });

    } else {
        geneInfoHTML += '<p>No additional structural data available.</p>';
    }

    // Add the modal container for the enlarged image at the end of the gene info div
    geneInfoHTML += `
        <!-- Image Modal -->
        <div id="image-modal" class="modal" onclick="closeImageModal()">
            <span class="close-button" onclick="closeImageModal()">&times;</span>
            <img class="modal-content" id="modal-image">
        </div>
    `;

    // Set the inner HTML of the gene info div
    geneInfoDiv.innerHTML = geneInfoHTML;

    // Ensure the viewer is loaded after the content is inserted into the DOM
    if (typeof relaxed_pdb !== 'undefined') {
        // Delay initialization slightly to ensure the DOM is fully rendered
        setTimeout(() => {
            for (const [transcriptId, pdbContent] of Object.entries(relaxed_pdb)) {
                // Load the PDB content into the 3D viewer for this transcript
                loadPDBInCard(transcriptId, pdbContent);
            }
        }, 0); // Delay by 0 milliseconds to allow the DOM to update
    }
}

function loadPDBInCard(transcriptId, pdbContent) {
    const viewerContainerId = `protein-viewer-${transcriptId}`;
    const viewerContainer = document.getElementById(viewerContainerId);

    if (!viewerContainer) return;

    // Set position to relative
    viewerContainer.style.position = 'relative';

    // Initialize 3Dmol.js Viewer for this card
    const viewer = $3Dmol.createViewer(viewerContainer, {
        defaultcolors: $3Dmol.rasmolElementColors
    });

    // Load the PDB content into the viewer
    viewer.addModel(pdbContent, "pdb"); // Load PDB from string
    viewer.setStyle({}, { cartoon: { color: 'spectrum' } });
    viewer.zoomTo();
    viewer.render();
}

// Function to open the image modal
function openImageModal(imageSrc) {
    const modal = document.getElementById('image-modal');
    const modalImage = document.getElementById('modal-image');

    modal.style.display = 'block';
    modalImage.src = imageSrc;
}

// Function to close the image modal
function closeImageModal() {
    const modal = document.getElementById('image-modal');
    modal.style.display = 'none';
}
