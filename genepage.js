// Extract the 'gene' parameter from the URL
const urlParams = new URLSearchParams(window.location.search);
const gene = urlParams.get('gene');

// Select the gene-info div where you want to display the information
const geneInfoDiv = document.getElementById('gene-info');
const mutantPlotDiv = document.getElementById('mutant-plot');
const mutantCardsDiv = document.getElementById('mutatnt-table');
const loadingElement = document.getElementById('loading');


// Function to load the pre-rendered mutant table
function loadMutationTable(gene) {
    const tableFilePath = `output/${gene}/clinvar_seqMUT_scores_summary.html`;

    // Fetch the pre-rendered table HTML and embed it into the mutantCardsDiv
    fetch(tableFilePath)
        .then(response => {
            if (!response.ok) {
                throw new Error('Network response was not ok');
            }
            return response.text();
        })
        .then(html => {
            // Hide the loading message
            loadingElement.style.display = 'none';

            // Insert the pre-rendered HTML table into the div
            mutantCardsDiv.innerHTML = html;
        })
        .catch(error => {
            console.error('Error fetching the mutant table:', error);
            mutantCardsDiv.innerHTML = `
                <div class="card">
                    <h2>Error</h2>
                    <p>There was an issue loading the mutant table. Please try again later.</p>
                </div>
            `;
        });
}
// Fetch the gene_info.json file for the specific gene
if (gene) {
    const filePath = `output/${gene}/gene_info.json`;
    const plotFilePath = `output/${gene}/isoform_plot.html`;

    // Fetch gene info JSON
    fetch(filePath)
        .then(response => response.json())
        .then(data => {
            if (data && data[gene]) {
                loadingElement.style.display = 'none';
                const geneInfo = data[gene];

                geneInfoDiv.innerHTML = `
                    <div class="card">
                        <h2>Gene Information for ${gene}</h2>
                        <p><strong>Description:</strong> ${geneInfo}</p>
                    </div>
                `;
            } else {
                geneInfoDiv.innerHTML = `
                    <div class="card">
                        <h2>Error: Gene not found</h2>
                        <p>The gene you are looking for is not in our database.</p>
                    </div>
                `;
            }
        })
        .catch(error => {
            console.error('Error fetching the gene info:', error);
            geneInfoDiv.innerHTML = `
                <div class="card">
                    <h2>Error</h2>
                    <p>There was an issue fetching the gene information. Please try again later.</p>
                </div>
            `;
        });

    // Fetch the isoform plot HTML and embed it into the mutantPlotDiv
    fetch(plotFilePath)
        .then(response => response.text())
        .then(html => {
            // Hide the loading message and display the embedded HTML
            document.getElementById('loading').style.display = 'none';

            // Create a temporary container to evaluate the HTML and inject scripts
            const tempDiv = document.createElement('div');
            tempDiv.innerHTML = html;

            // Append the contents of the temp div to the mutant plot div
            while (tempDiv.firstChild) {
                mutantPlotDiv.appendChild(tempDiv.firstChild);
            }

            // Manually execute any scripts in the loaded HTML
            const scripts = mutantPlotDiv.getElementsByTagName('script');
            for (let i = 0; i < scripts.length; i++) {
                const script = document.createElement('script');
                script.textContent = scripts[i].textContent;
                document.body.appendChild(script);
            }
        })
        .catch(error => {
            console.error('Error fetching the plot:', error);
            mutantPlotDiv.innerHTML = `
                <div class="card">
                    <h2>Error</h2>
                    <p>There was an issue loading the plot. Please try again later.</p>
                </div>
            `;
        });
        loadMutationTable(gene);
} else {
    // If no gene parameter is found, redirect to 404
    window.location.href = '404.html';
}
