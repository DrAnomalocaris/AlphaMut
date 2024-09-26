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
    const genes = [...new Set(Object.values(data))];
    
    // Clear the container before appending new items (optional)
    foldedProteinsContainer.innerHTML = ''; 

    genes.forEach(gene => {
        const listItem = document.createElement('li'); // Create a list item
        const link = document.createElement('a'); // Create an anchor element
        
        // Set the text and href attributes for the link
        link.textContent = gene; // Set the display text
        link.href = `gene.html?gene=${encodeURIComponent(gene)}`; // Set the URL
        
        // Append the link to the list item and the item to the container
        listItem.appendChild(link); 
        foldedProteinsContainer.appendChild(listItem);
    });
})
.catch(error => {
    console.error('There was a problem with the fetch operation:', error);
});
