var velocity = 0.1;

document.addEventListener("DOMContentLoaded", function () {
  var slider = document.getElementById("myRange");
  var output = document.getElementById("demo");
  output.innerHTML = slider.value / 10;

  slider.oninput = function () {
    output.innerHTML = this.value / 10;
  };

  velocity = output.innerHTML;
});

function generateValues() {
  // Simulating a 3-second delay
  // Get the number of stations entered by the user
  const stationCount = document.getElementById("stationCount").value;

  // Get the table element
  const table = document.getElementById("outputTable");

  // Clear existing table rows
  table.innerHTML = "";

  // Create table header
  const headerRow = table.insertRow(0);
  const xHeader = headerRow.insertCell(0);
  const yHeader = headerRow.insertCell(1);
  const zHeader = headerRow.insertCell(2);
  xHeader.innerHTML = "<b>X</b>";
  yHeader.innerHTML = "<b>Y</b>";
  zHeader.innerHTML = "<b>Z</b>";

  // Generate rows for each station
  for (let i = 1; i <= stationCount; i++) {
    const row = table.insertRow(i);
    const xCell = row.insertCell(0);
    const yCell = row.insertCell(1);
    const zCell = row.insertCell(2);

    // Create input fields for X, Y, Z
    const xInput = document.createElement("input");
    xInput.type = "text";
    xInput.placeholder = "Enter X value";
    xCell.appendChild(xInput);

    const yInput = document.createElement("input");
    yInput.type = "text";
    yInput.placeholder = "Enter Y value";
    yCell.appendChild(yInput);

    const zInput = document.createElement("input");
    zInput.type = "text";
    zInput.placeholder = "Enter Z value";
    zCell.appendChild(zInput);
  }
  document.getElementById("submitBtn").style.display = "block";
  document.getElementById("slide-container").style.display = "block";
  document.getElementById("actual-earthquake").style.display = "block";
}

var sta_loc = [
  [35, 9, 0],
  [-44, 10, 0],
  [-11, -25, 0],
  [23, -39, 0],
  [42, -27, 0],
  [-12, 50, 0],
  // [-45, 16, 0],
  // [12, -16, 0],
  // [5, -19, 0],
  // [-1, -11, 0],
  // [20, 11, 0],
  // [15, 8, 0],
  // [-45, 12, 0],
  // [-13, 45, 0],
  // [-26, 32, 0],
  // [-45, 65, 0],
  // [-23, 12, 0],
];
let x, y, z;

function submitValues() {
  // sta_loc = [];

  // const table = document.getElementById("outputTable");

  // const rows = table.querySelectorAll("tr");

  // for (let i = 1; i < rows.length; i++) {
  //   const cells = rows[i].querySelectorAll("input");

  //   const coordinates = [];
  //   cells.forEach((input) => {
  //     const value = input.value.trim();
  //     if (value !== "") {
  //       coordinates.push(parseFloat(value));
  //     }
  //   });
  //   sta_loc.push(coordinates);
  // }

  // Log updated sta_loc array
  console.log({ sta_loc: sta_loc });

  document.getElementById("graphs").style.display = "block";
  document.getElementById("Input-container").style.display = "none";
  showGraphs();
}

function showGraphs() {
  var margin = { top: 10, right: 30, bottom: 30, left: 60 },
    width = 460 - margin.left - margin.right,
    height = 400 - margin.top - margin.bottom;

  x = document.getElementById("xInput").value;
  y = document.getElementById("yInput").value;
  z = document.getElementById("zInput").value;

  const eqloc_actual = [x, y, z];
  const v0 = velocity; // km/s // Homogeneous earth velocity
  const t0 = 0; // origin time of earthquake

  console.log({ actual: eqloc_actual });

  // Station locations (x y z)
  var size = 100,
    x = new Array(size),
    y = new Array(size),
    z = new Array(size);

  for (var i = 0; i < size; i++) {
    x[i] = y[i] = -2 * Math.PI + (4 * Math.PI * i) / size;
    z[i] = new Array(size);
  }

  // X, Y, Z coordinates of stations
  const xrec = sta_loc.map((coords) => coords[0]);
  const yrec = sta_loc.map((coords) => coords[1]);
  const zrec = sta_loc.map((coords) => coords[2]);

  const minX = Math.min(...xrec);
  const maxX = Math.max(...xrec);
  const minY = Math.min(...yrec);
  const maxY = Math.max(...yrec);

  // Number of data and model
  const N = sta_loc.length; // 10 stations (= P arrival times)
  const M = 5; // xr, yr, zr, t0, v

  // Plot station-receiver configuration
  var svg = d3
    .select("#u")
    .append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

  var xScale = d3
    .scaleLinear()
    .domain([minX, maxX * 1.5])
    .range([0, width]);

  var yScale = d3
    .scaleLinear()
    .domain([minY - 20, maxY * 1.5])
    .range([height, 0]);

  var tooltip = svg.append("text").attr("class", "tooltip").style("opacity", 0);

  // Plot stations
  svg
    .selectAll(".station")
    .data(sta_loc)
    .enter()
    .append("circle")
    .attr("class", "station")
    .attr("cx", (d) => xScale(d[0]))
    .attr("cy", (d) => yScale(d[1]))
    .attr("r", 5)
    .attr("fill", "black")
    .on("mouseover", function (event, d) {
      // Show tooltip with coordinates

      const circle = d3.select(this); // Select the circle that triggered the event
      const circleRect = circle.node().getBoundingClientRect(); // Get the bounding rectangle of the circle

      tooltip
        .html("(" + d[0] + ", " + d[1] + ")")
        .attr("x", circleRect.x - 60)
        .attr("y", circleRect.y - 200) // Adjust this value as needed
        .style("opacity", 1);
    })
    .on("mouseout", function () {
      // Hide tooltip
      tooltip.style("opacity", 0);
    });

  // Plot actual location
  svg
    .selectAll(".actual-location")
    .data([eqloc_actual])
    .enter()
    .append("circle")
    .attr("class", "actual-location")
    .attr("cx", (d) => xScale(d[0]))
    .attr("cy", (d) => yScale(d[1]))
    .attr("r", 5)
    .attr("fill", "red")
    .on("mouseover", function (event, d) {
      // Show tooltip with coordinates
      const circle = d3.select(this); // Select the circle that triggered the event
      const circleRect = circle.node().getBoundingClientRect(); // Get the bounding rectangle of the circle

      tooltip
        .html("(" + d[0] + ", " + d[1] + ")")
        .attr("x", circleRect.x - 60)
        .attr("y", circleRect.y - 200) // Adjust this value as needed
        .style("opacity", 1);
    })
    .on("mouseout", function () {
      // Hide tooltip
      tooltip.style("opacity", 0);
    });

  svg
    .append("g")
    .attr("transform", "translate(0," + height + ")")
    .call(d3.axisBottom(xScale));

  svg.append("g").call(d3.axisLeft(yScale));

  // Compute travel-time
  const travel_time = (xr, yr, zr, xs, ys, zs, t, v) =>
    t + Math.sqrt((xr - xs) ** 2 + (yr - ys) ** 2 + (zr - zs) ** 2) / v;

  // Compute synthetic observed data
  let t_obs = xrec.map((xr, idx) =>
    travel_time(
      xr,
      yrec[idx],
      zrec[idx],
      eqloc_actual[0],
      eqloc_actual[1],
      eqloc_actual[2],
      t0,
      v0
    )
  );

  console.log({ t_obs: t_obs });

  // Add noise to data
  const sigma_d = 0.1;
  t_obs = t_obs.map((t) => t + Math.random() * sigma_d);

  var size = sta_loc.length;

  // Create an empty 2D array to represent the weight matrix
  var W = [];

  // Loop to populate the weight matrix with diagonal elements
  for (var i = 0; i < size; i++) {
    W[i] = [];
    for (var j = 0; j < size; j++) {
      // Set diagonal elements to 1/sigma_d, off-diagonal elements to 0
      if (i === j) {
        W[i][j] = 1 / sigma_d;
      } else {
        W[i][j] = 0;
      }
    }
  }

  const Cd = Array.from({ length: N }, (_, i) =>
    Array.from({ length: N }, (_, j) => (i === j ? sigma_d * sigma_d : 0))
  );

  const Nsamples = 10000;
  const xsamp = Array.from(
    { length: Nsamples },
    () => minX + (maxX - minX) * 2 * Math.random()
  );

  const ysamp = Array.from(
    { length: Nsamples },
    () => minY + (maxY - minY) * 2 * Math.random()
  );
  const zsamp = Array(Nsamples).fill(eqloc_actual[2]);
  const tsamp = Array(Nsamples).fill(t0);
  const vsamp = Array(Nsamples).fill(v0);

  const t_samp = [];
  const residual_samp = [];
  const residual_samp_all = [];

  for (let ii = 0; ii < Nsamples; ii++) {
    const t_samp_current = xrec.map((xr, idx) =>
      travel_time(
        xr,
        yrec[idx],
        zrec[idx],
        xsamp[ii],
        ysamp[ii],
        zsamp[ii],
        tsamp[ii],
        vsamp[ii]
      )
    );

    t_samp.push(t_samp_current);

    const residual_current = t_obs.map((t, idx) => t - t_samp_current[idx]);
    residual_samp.push(residual_current);

    const sumOfSquares = residual_current.reduce(
      (sum, val) => sum + val ** 2,
      0
    );

    residual_samp_all.push(sumOfSquares / t_obs.length);
  }

  // Define color scale
  console.log({ maxi: d3.max(residual_samp) });
  const colorScale = d3
    .scaleSequential(d3.interpolateTurbo) // or any other color scale you prefer
    .domain([0, d3.max(residual_samp_all)]);

  const minXsamp = Math.min(...xsamp);
  const maxXsamp = Math.max(...xsamp);
  const minYsamp = Math.min(...ysamp);
  const maxYsamp = Math.max(...ysamp);

  // Calculate the range to include the circles
  const xRange = [0, width];
  const yRange = [height, 0];

  // Adjust range to cover all circles
  if (minXsamp < 0) xRange[0] += Math.abs(minXsamp);
  if (maxXsamp * 1.5 > width) xRange[1] += maxXsamp * 1.5 - width;
  if (minYsamp < 0) yRange[1] += Math.abs(minYsamp);
  if (maxYsamp * 1.5 > height) yRange[0] += maxYsamp * 1.5 - height;

  var xScale = d3
    .scaleLinear()
    .domain([minXsamp, maxXsamp * 1.5])
    .range([0, width]);

  var yScale = d3
    .scaleLinear()
    .domain([minYsamp, maxYsamp * 1.5])
    .range([height, 0]);

  // Adjust view box and scales
  const viewBoxWidth = maxXsamp * 10;
  const viewBoxHeight = maxYsamp * 10;

  // Create SVG element
  const GridSvg = d3
    .select("#grid-search-results")
    .append("svg")
    .attr("width", 800)
    .attr("height", 400)
    .attr("viewBox", `0 0 ${200} ${400}`);

  GridSvg.selectAll(".grid-search-result")
    .data(xsamp)
    .enter()
    .append("circle")
    .attr("class", "grid-search-result")
    .attr("cx", (d) => xScale(d))
    .attr("cy", (d, i) => yScale(ysamp[i]))
    .attr("r", 6)
    .attr("fill", (d, i) => colorScale(Math.abs(residual_samp_all[i])))
    .attr("opacity", 0.9)
    .attr("stroke", "black")
    .attr("stroke-width", 1);

  // Create circles for stations
  GridSvg.selectAll(".station")
    .data(xrec)
    .enter()
    .append("circle")
    .attr("class", "station")
    .attr("cx", (d) => xScale(d))
    .attr("cy", (d, i) => yScale(yrec[i]))
    .attr("r", 7)
    .attr("fill", "black");

  // Create circle for actual location
  GridSvg.append("circle")
    .attr("cx", xScale(eqloc_actual[0]))
    .attr("cy", yScale(eqloc_actual[1]))
    .attr("r", 7)
    .attr("fill", "red");

  // Add title
  GridSvg.append("text")
    .attr("x", width / 1.5)
    .attr("y", 20)
    .attr("text-anchor", "middle")
    .text("Residuals for Random Samples (Forward Simulation - No Inversion)");

  // Start Inversion
  // Inital model (or starting model) - Prior
  const eqloc_initial = [10, 5, 20]; // assumed earthquake location
  const ti = 2; // assumed origin at 2 sec
  const vi = 4; // km/s // assumed velocity

  // ... (plotting initial model, similar setup to MATLAB)

  // model parameters
  let m1 = eqloc_initial[0];
  let m2 = eqloc_initial[1];
  let m3 = eqloc_initial[2];
  let m4 = ti;
  let m5 = vi;

  // Create functions for estimating partial derivative matrix (G)
  const dgdx = (xr, yr, zr, xs, ys, zs, v) =>
    (xs - xr) /
    (v * Math.sqrt((xr - xs) ** 2 + (yr - ys) ** 2 + (zr - zs) ** 2));
  const dgdy = (xr, yr, zr, xs, ys, zs, v) =>
    (ys - yr) /
    (v * Math.sqrt((xr - xs) ** 2 + (yr - ys) ** 2 + (zr - zs) ** 2));
  const dgdz = (xr, yr, zr, xs, ys, zs, v) =>
    (zs - zr) /
    (v * Math.sqrt((xr - xs) ** 2 + (yr - ys) ** 2 + (zr - zs) ** 2));
  const dgdt = 1;
  const dgdv = (xr, yr, zr, xs, ys, zs, v) =>
    (1 / v ** 2) * Math.sqrt((xr - xs) ** 2 + (yr - ys) ** 2 + (zr - zs) ** 2);

  const n_iterations = 4;

  var residuals = [];

  function matrixConditionNumber(matrix) {
    // Compute the singular value decomposition (SVD) of the matrix
    const svd = numeric.svd(matrix);

    // Get the singular values
    const singularValues = svd.S;

    // Compute the condition number using the Frobenius norm
    const maxSingularValue = Math.max(...singularValues);
    const minSingularValue = Math.min(...singularValues);

    // Condition number is the ratio of the maximum singular value to the minimum singular value
    const conditionNumber = maxSingularValue / minSingularValue;

    return conditionNumber;
  }

  for (let ii = 0; ii < n_iterations; ii++) {
    // Compute partial derivative matrix (G)
    var G = Array.from({ length: N }, (_, idx) => [
      dgdx(xrec[idx], yrec[idx], zrec[idx], m1, m2, m3, m5),
      dgdy(xrec[idx], yrec[idx], zrec[idx], m1, m2, m3, m5),
      dgdz(xrec[idx], yrec[idx], zrec[idx], m1, m2, m3, m5),
      1,
      dgdv(xrec[idx], yrec[idx], zrec[idx], m1, m2, m3, m5),
    ]);

    const ts_ii = xrec.map((xr, idx) =>
      travel_time(xr, yrec[idx], zrec[idx], m1, m2, m3, m4, m5)
    );
    console.log({ W: W });
    console.log({ G: G });

    const delta_d = numeric.sub(t_obs, ts_ii);
    const Gw = numeric.dot(W, G);
    const dw = numeric.dot(W, delta_d);

    console.log({ delta_d: delta_d });
    console.log({ Gw: Gw });
    console.log({ dw: dw });

    var GwTranspose = numeric.transpose(Gw);
    console.log({ GwTranspose: GwTranspose });
    var GwTransposeTimesGw = numeric.dot(GwTranspose, Gw);
    console.log({ GwTransposeTimesGw });
    var GwPseudoInverse = numeric.inv(GwTransposeTimesGw);
    console.log({ GwPseudoInverse });
    var dm = numeric.dot(numeric.dot(GwPseudoInverse, GwTranspose), dw);

    // G = matrixConditionNumber(G);

    console.log({ G: G });

    // Residual
    var residualSum = 0;
    for (var idx = 0; idx < t_obs.length; idx++) {
      var delta = t_obs[idx] - ts_ii[idx];
      residualSum += delta * delta;
    }
    residuals[ii] = residualSum / t_obs.length / 1000;

    // Update model
    m1 += dm[0];
    m2 += dm[1];
    m3 += dm[2];
    m4 += dm[3];
    m5 += dm[4];

    console.log({ m1: m1 });
    console.log({ m2: m2 });
  }

  const plotSvg = d3
    .select("#plot-svg")
    .append("svg")
    .attr("width", width + margin.left + margin.right + 100)
    .attr("height", height + margin.top + margin.bottom)
    .append("g")
    .attr("transform", `translate(${margin.left + 100},${margin.top})`);

  // Define scales
  var xScale = d3
    .scaleLinear()
    .domain([0, residuals.length]) // Assuming 1-based iteration numbering
    .range([0, width]);

  var yScale = d3
    .scaleLinear()
    .domain([0, d3.max(residuals)])
    .range([height, 0]);

  // Plot the line
  console.log(residuals);
  plotSvg
    .selectAll(".line")
    .data([residuals])
    .enter()
    .append("path")
    .attr("class", "line")
    .attr("fill", "none")
    .attr("stroke", "steelblue")
    .attr("stroke-width", 2)
    .attr(
      "d",
      d3
        .line()
        .x((d, i) => xScale(i + 1))
        .y((d) => yScale(d))
    );

  // Plot the markers
  plotSvg
    .selectAll(".dot")
    .data(residuals)
    .enter()
    .append("circle")
    .attr("class", "dot")
    .attr("cx", (d, i) => xScale(i + 1))
    .attr("cy", (d) => yScale(d))
    .attr("r", 5)
    .attr("fill", "red")
    .attr("stroke", "black")
    .attr("stroke-width", 1);

  // Add x axis
  plotSvg
    .append("g")
    .attr("transform", `translate(0, ${height})`)
    .call(d3.axisBottom(xScale).ticks(residuals.length));

  // Add y axis
  plotSvg.append("g").call(d3.axisLeft(yScale));
  // Add labels
  plotSvg
    .append("text")
    .attr("x", width / 2)
    .attr("y", height + margin.top + 20)
    .style("text-anchor", "middle")
    .text("Number of Iterations");

  plotSvg
    .append("text")
    .attr("transform", "rotate(-90)")
    .attr("x", 0 - height / 2)
    .attr("y", 0 - margin.left - 60)
    .attr("dy", "1em")
    .style("text-anchor", "middle")
    .text("Residual");
}
