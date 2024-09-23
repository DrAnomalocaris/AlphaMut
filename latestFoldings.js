// latestFoldings.js

document.addEventListener("DOMContentLoaded", () => {
    // Ensure geneDict is available
    if (typeof geneDict !== 'undefined') {
        const foldedProteinsContainer = document.getElementById('folded-proteins-container');

        if (foldedProteinsContainer) {
            const uniqueEnsgIds = new Set();
            const ensgIdToGeneName = {};

            for (const [key, value] of Object.entries(geneDict)) {
                uniqueEnsgIds.add(value);

                // If the key is a gene symbol (not starting with 'ENSG')
                if (!key.startsWith('ENSG')) {
                    ensgIdToGeneName[value] = key;
                }
            }

            const foldedProteins = Array.from(uniqueEnsgIds).map(ensgId => {
                return {
                    ensgId: ensgId,
                    geneName: ensgIdToGeneName[ensgId] || ensgId,
                };
            });

            // Sort the folded proteins alphabetically by gene name
            foldedProteins.sort((a, b) => a.geneName.localeCompare(b.geneName));

            // Create and append the list to the container
            const heading = document.createElement('h2');
            heading.textContent = 'Available Folded Proteins';
            foldedProteinsContainer.appendChild(heading);

            const list = document.createElement('ul');
            list.className = 'folded-proteins-list';

            foldedProteins.forEach(protein => {
                const listItem = document.createElement('li');
                const link = document.createElement('a');
                link.href = `gene.html?gene=${encodeURIComponent(protein.ensgId)}`;
                link.textContent = `${protein.geneName} (${protein.ensgId})`;
                listItem.appendChild(link);
                list.appendChild(listItem);
            });

            foldedProteinsContainer.appendChild(list);
        }
    }
});
