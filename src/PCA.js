
// Load data from CSV file
d3.csv(`output/${gene}/PCA.csv`).then(function(data) {
    data.forEach((d, i) => {
        d.PC1 = +d.PC1;
        d.PC2 = +d.PC2;
        d.rmsd = +d.rmsd;
    });

    // Set up dimensions and margins
    const margin = {top: 20, right: 30, bottom: 40, left: 50},
        width = 800 - margin.left - margin.right,
        height = 600 - margin.top - margin.bottom;

    // Create SVG container with zoom capability
    const svg = d3.select("#PCA")
                .append("svg")
                .attr("width", width + margin.left + margin.right)
                .attr("height", height + margin.top + margin.bottom)
                .call(d3.zoom().scaleExtent([0.5, 10]).on("zoom", zoomed))
                .append("g")
                .attr("transform", `translate(${margin.left},${margin.top})`);

    // Set up scales
    const x = d3.scaleLinear()
                .domain(d3.extent(data, d => d.PC1))
                .range([0, width]);

    const y = d3.scaleLinear()
                .domain(d3.extent(data, d => d.PC2))
                .range([height, 0]);

    const color = d3.scaleOrdinal()
                    .domain([...new Set(data.map(d => d.ClinicalSignificance))])
                    .range(d3.schemeCategory10);

    // Add tooltip
    const tooltip = d3.select("body")
                    .append("div")
                    .attr("class", "tooltip")
                    .style("opacity", 0);

  // Add points
// Add points
const points = svg.selectAll(".point")
.data(data)
.enter()
.append(function(d) {
  if (d.mutation === 'canonical') {
      return document.createElementNS(d3.namespaces.svg, 'path');
  } else {
      return document.createElementNS(d3.namespaces.svg, 'circle');
  }
})
.attr("class", d => d.mutation === 'canonical' ? 'canonical-point' : 'circle-point')
.attr("d", d => d.mutation === 'canonical' ? d3.symbol().type(d3.symbolCross).size(100)() : null)
.attr("transform", d => d.mutation === 'canonical' ? `translate(${x(d.PC1)},${y(d.PC2)})` : null)
.attr("cx", d => d.mutation !== 'canonical' ? x(d.PC1) : null)
.attr("cy", d => d.mutation !== 'canonical' ? y(d.PC2) : null)
.attr("r", d => d.mutation !== 'canonical' ? 8 : null)  // Larger radius for better visibility
.style("fill", d => color(d.ClinicalSignificance))
.style("opacity", 0.7)
.on("mouseover", function(event, d) {
      tooltip.transition()
             .duration(200)
             .style("opacity", .9);
      tooltip.html(`ID: ${d.mutation}<br>Phenotype: ${d.PhenotypeList}<br>Significance: ${d.ClinicalSignificance}<br>RMSD: ${d.rmsd}`)
             .style("left", (event.pageX + 10) + "px")
             .style("top", (event.pageY - 28) + "px");

      points.style("opacity", point => point.PhenotypeList === d.PhenotypeList ? 1 : 0.1);
})
.on("mouseout", function() {
      tooltip.transition()
             .duration(500)
             .style("opacity", 0);

      points.style("opacity", 0.7);
})
.on("click", function(event, d) {
      window.location.href = `isoform.html?gene=${encodeURIComponent(gene)}&mutation=${encodeURIComponent(d.mutation)}`;
});

// Zoom function
function zoomed({transform}) {
svg.attr("transform", transform);
svg.selectAll(".canonical-point")
.attr("transform", d => `translate(${transform.applyX(x(d.PC1))},${transform.applyY(y(d.PC2))})`);
svg.selectAll(".circle-point")
.attr("cx", d => transform.applyX(x(d.PC1)))
.attr("cy", d => transform.applyY(y(d.PC2)));
}
});