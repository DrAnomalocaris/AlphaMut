// latestFoldings.js
const foldedProteinsContainer = document.getElementById('folded-proteins-container');
fetch('completed.json')
.then(response => {
    if (!response.ok) {
        throw new Error('Network response was not ok');
    }
    return response.json(); // Parse JSON data
})
.then(data => {
    // Assuming data is an array of completed genes
    data.forEach(gene => {
        const listItem = document.createElement('li');
        listItem.textContent = gene; // Change according to your JSON structure
        foldedProteinsContainer.appendChild(listItem);
    });
})
.catch(error => {
    console.error('There was a problem with the fetch operation:', error);
});