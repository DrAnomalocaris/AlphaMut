
// Load the JSON file
fetch('completed.json')
    .then(response => response.json()) // Parse JSON data
    .then(data => {
        displayResults(data); // Call the function to display results
    })
    .catch(error => {
        console.error('Error loading the JSON file:', error);
    });
// Levenshtein Distance Function
function levenshteinDistance(a, b) {
    const matrix = [];

    for (let i = 0; i <= b.length; i++) {
        matrix[i] = [i];
    }

    for (let j = 0; j <= a.length; j++) {
        matrix[0][j] = j;
    }

    for (let i = 1; i <= b.length; i++) {
        for (let j = 1; j <= a.length; j++) {
            if (b.charAt(i - 1) === a.charAt(j - 1)) {
                matrix[i][j] = matrix[i - 1][j - 1];
            } else {
                matrix[i][j] = Math.min(
                    matrix[i - 1][j - 1] + 1,
                    Math.min(matrix[i][j - 1] + 1, matrix[i - 1][j] + 1)
                );
            }
        }
    }

    return matrix[b.length][a.length];
}

// Sort Genes by Similarity
function sortBySimilarity(searchTerm, strings) {
    return strings.sort((a, b) => {
        const distanceA = levenshteinDistance(searchTerm.toLowerCase(), a.toLowerCase());
        const distanceB = levenshteinDistance(searchTerm.toLowerCase(), b.toLowerCase());
        return distanceA - distanceB; // Sort by distance
    });
}
const urlParams = new URLSearchParams(window.location.search)
const searchTerm = urlParams.get('q'); // Get the 'q' parameter from the URL
document.getElementById('search-query').innerHTML = searchTerm;


// Display Results Based on the Query Parameter
function displayResults(genes) {
    const urlParams = new URLSearchParams(window.location.search);
    const searchTerm = urlParams.get('q'); // Get the 'q' parameter from the URL
    const resultsList = document.getElementById('results');

    resultsList.innerHTML = ''; // Clear previous results

    if (searchTerm) {
        const sortedResults = sortBySimilarity(searchTerm, Object.keys(genes));

        sortedResults.forEach(gene => {
            const listItem = document.createElement('li');
            const link = document.createElement('a');
            link.href = `gene.html?gene=${encodeURIComponent([gene])}`; // Adjust to your gene page
            link.textContent = gene; // Display the gene name
            listItem.appendChild(link);
            resultsList.appendChild(listItem);
        });
    } else {
        resultsList.innerHTML = '<li>No search term provided.</li>';
    }
}
