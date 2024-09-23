document.addEventListener("DOMContentLoaded", () => {
    // Common code for all pages (Search bar functionality)
    const searchInput = document.getElementById('search-input');
    const autocompleteList = document.getElementById('autocomplete-list');
    const searchButton = document.getElementById('search-button'); 

    // Check if the search elements exist before adding event listeners
    if (searchInput && autocompleteList && searchButton) {
        // Use the autofill array directly
        const exampleGenes = autofill;

        searchInput.addEventListener('input', () => {
            const inputValue = searchInput.value;
            autocompleteList.innerHTML = '';

            if (inputValue.length === 0) {
                return;
            }

            // Filter and display suggestions
            const suggestions = exampleGenes.filter(gene => gene.toLowerCase().includes(inputValue.toLowerCase()));
            suggestions.forEach(suggestion => {
                const item = document.createElement('div');
                item.textContent = suggestion;
                item.addEventListener('click', () => {
                    searchInput.value = suggestion;
                    autocompleteList.innerHTML = ''; // Clear the autocomplete list
                });
                autocompleteList.appendChild(item);
            });
        });

        // Hide autocomplete list when clicking outside
        document.addEventListener('click', (e) => {
            if (e.target !== searchInput) {
                autocompleteList.innerHTML = '';
            }
        });

        // Handle the send button click
        searchButton.addEventListener('click', () => {
            const gene = searchInput.value.trim();

            if (gene) {
                // Check if the gene exists in geneDict
                const ensgCode = geneDict[gene];

                if (ensgCode) {
                    // Redirect to gene.html with the ENSG code as a query parameter
                    window.location.href = `gene.html?gene=${encodeURIComponent(ensgCode)}`;
                } else {
                    // Gene not found in geneDict, but proceed anyway
                    // Redirect to gene.html with the user input as the gene parameter
                    window.location.href = `gene.html?gene=${encodeURIComponent(gene)}`;
                }
            } else {
                alert('Please enter a gene to search.');
            }
        });
    }

    // Gene-specific code for gene.html
    const currentPage = window.location.pathname.split('/').pop();

    if (currentPage === 'gene.html') {
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

        // Function to display gene information using ENSG code directly
        function displayGeneInformation(ensgId) {
            // ... your existing code to display gene information ...
        }
    }
});
