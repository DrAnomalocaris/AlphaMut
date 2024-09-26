let autofill = [];

// Load the JSON file
fetch('completed.json')
    .then(response => response.json()) // Parse JSON data
    .then(data => {
        autofill = Object.keys(data); // Extract keys from the data; // Assign data to the autofill variable
        
        // Call the function to set up event listeners after data is loaded
        setupSearchFunctionality();
    })
    .catch(error => {
        console.error('Error loading the JSON file:', error);
    });

function setupSearchFunctionality() {
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
                // Gene not found in geneDict, but proceed anyway
                // Redirect to gene.html with the user input as the gene parameter
                window.location.href = `search.html?q=${encodeURIComponent(gene)}`;
            } else {
                alert('Please enter a gene to search.');
            }
        });
    }

    // Gene-specific code for gene.html
    const currentPage = window.location.pathname.split('/').pop();
}

// Ensure the search functionality is set up only after the DOM is fully loaded
document.addEventListener("DOMContentLoaded", () => {
    // The fetch logic and setupSearchFunctionality() are already handled above,
    // so no need to duplicate the code here.
});
