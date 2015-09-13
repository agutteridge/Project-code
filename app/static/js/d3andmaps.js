// initialize Google Maps
function initialize() {
  var mapOptions = {
    // centres on London
    center: { lat: 0, lng: 0},
    zoom: 1,
    mapTypeControlOptions: {
      mapTypeIds: [google.maps.MapTypeId.ROADMAP, 'map_style']
    }
  };
  
  // Style for the map
  var styles = [
    {
      stylers: [
        { hue: "#00ffe6" },
        { saturation: -20 }
      ]
    },{
      featureType: "road",
      elementType: "geometry",
      stylers: [
        { lightness: 100 },
        { visibility: "simplified" }
      ]
    },{
      featureType: "road",
      elementType: "labels",
      stylers: [
        { visibility: "off" }
      ]
    },{
      featureType: "poi.park",
      elementType: "labels",
      stylers: [
        { visibility: "off" }
      ]
    },{
      featureType: "poi.attraction",
      elementType: "labels",
      stylers: [
        { visibility: "off" }
      ]
    },{
      featureType: "poi.place_of_worship",
      elementType: "labels",
      stylers: [
        { visibility: "off" }
      ]
    },{
      featureType: "poi.sports_complex",
      elementType: "labels",
      stylers: [
        { visibility: "off" }
      ]
    }
  ];

  var styledMap = new google.maps.StyledMapType(styles, {name: "Styled Map"});
  map = new google.maps.Map(document.getElementById('map-canvas'), mapOptions);
  map.mapTypes.set('map_style', styledMap);
  map.setMapTypeId('map_style');
}

function placeMarkers(data) {
  // Removes old markers
  clearMarkers();

  var marker
  var i
  var infowindow = new google.maps.InfoWindow();

  for (i = 0; i < data.length; i++) {
    var opacity = 0.2;

    if (data[i].PMID === firstpmid) {
      opacity = 1;
    };

    marker = new google.maps.Marker({
      position: new google.maps.LatLng(data[i].place.geometry.location.lat, 
        data[i].place.geometry.location.lng),
      map: map,
      pmid: data[i].PMID,
      opacity: opacity,
      address: data[i].place.name
    });

    // add marker to global array
    markers.push(marker)

    google.maps.event.addListener(marker, 'click', (function(marker, i) {
      return function() {
        infowindow.setContent(marker.address + "<br><button id=\"paperlink\">" + 
          marker.pmid + "</button>");
        infowindow.open(map, marker);
        thispaper = getPaperbyPMID(marker.pmid)
        document.getElementById("paperlink").addEventListener("click", function (e) { sendPaper(thispaper); });
      }
    })(marker, i));
  }
}

function clearMarkers() {
  for (var i = 0; i < markers.length; i++) {
    markers[i].setMap(null);
  }
  markers.length = 0;
}

function sendPaper(paper) {
  console.log('sending paper')
  console.log(paper)
  title = document.getElementById('bar-title');
  info = document.getElementById('bar-info');

  title.innerHTML = "<a href=\"http://ncbi.nlm.nih.gov/pubmed/" + 
    paper.PMID + "\" target=\"_blank\">" + 
    paper.title +
    "</a>";

  info.innerHTML = paper.authors + ", <i>" +
    paper.journal + "</i>, " +
    paper.date;

  var markerBounds = new google.maps.LatLngBounds();

  for (var i = 0; i < markers.length; i++) {
    if (markers[i].pmid === paper.PMID) {
      markers[i].setOpacity(1);
      markerBounds.extend(markers[i].position);
    } else {
      markers[i].setOpacity(0.2);
    };
  };

  google.maps.event.addListener(map, 'zoom_changed', function() {
      zoomChangeBoundsListener = 
          google.maps.event.addListener(map, 'bounds_changed', function(event) {
              if (this.getZoom() > 5 && this.initialZoom == true) {
                  console.log('zoom too close')
                  // Change max/min zoom here
                  this.setZoom(5);
                  this.initialZoom = false;
              }
          google.maps.event.removeListener(zoomChangeBoundsListener);
      });
  });
  map.initialZoom = true;

  if (!(markerBounds.isEmpty()) && (!initial)) {
    map.fitBounds(markerBounds);
  } else {
    initial = false;
    map.setCenter(new google.maps.LatLng(0, 0))
    map.setZoom(1)
  }
  changeColour(paper.PMID);
}

// global vars
var markers = [];
var map;
var papers = [];
var paper_position = 0; 
var initial = true;
var text_field = document.getElementById('text_field');
var btn_search = document.getElementById('btn_search');
var previous_btn = document.getElementById('previous_btn');
var next_btn = document.getElementById('next_btn');

function next() {
  if (paper_position > -1 && paper_position < 19) {
    sendPaper(papers[paper_position + 1]);
    paper_position += 1;
  };
}

function previous() {
  if (paper_position != 0) {
    sendPaper(papers[paper_position - 1]);
    paper_position -= 1;
  };
}

function getPaperbyPMID(pmid) {
  for (i = 0; i < papers.length; i++) {
    if (papers[i].PMID === pmid) {
      paper_position = i;
      return papers[i];
    };
  };
}

btn_search.onclick = function(){
  console.log(text_field.value);
  document.getElementById("loading").style.display = "block";
  d3.json("/data?search_term=" + text_field.value, d3_callback);
  svg.selectAll("circle,text").remove();
};

previous_btn.onclick = previous;
next_btn.onclick = next;

document.onkeydown = function(e) {
    e = e || window.event;
    switch(e.which || e.keyCode) {
      case 37: // left
        previous();
        break;

        case 39: // right
        next();
        break;

        default: return; // exit this handler for other keys
    }
    e.preventDefault(); // prevent the default action (scroll / move caret)
};

//D3: mostly taken from http://bl.ocks.org/mbostock/7607535
var margin = 0,
    diameter = 500;

var color = d3.scale.linear()
    .domain([-1, 5])
    .range(["hsl(174,100%,80%)", "hsl(180,10.3%,40%)"])
    .interpolate(d3.interpolateHcl);

var pack = d3.layout.pack()
    .padding(5)
    .size([diameter - margin, diameter - margin])
    .value(function(d) { return d.PMIDs.length * 500; });

var svg = d3.select("#svg-container").append("svg")
    .attr("width", diameter)
    .attr("height", diameter)
  .append("g")
    .attr("transform", "translate(" + diameter / 2 + "," + diameter / 2 + ")");

function changeColour (pmid) {
  svg.selectAll("circle")
    .style("fill", function(d) { return d.children ? color(d.depth) : colourNode(d, pmid); })
}

function colourNode(d, thispmid) {
  if ("PMIDs" in d) {
    if (d.PMIDs.indexOf(thispmid) > -1) {
      return "#f2756a";
    } else {
      return null
    };
  };
}

function d3_callback(error, data) {
  console.log('in d3_callback')
  if (error) throw error;

  if (data == 'No results found') {
    window.alert('No results found!\nPlease enter another search term.');
    document.getElementById("loading").style.display = "none";
  } else {
    // assigning global var papers
    papers = data.papers

    var root = data.concepts;
    // sets global var as first PMID so markers can be formatted accordingly
    firstpmid = data.papers[0].PMID;
    // passing location data to Google Maps function
    placeMarkers(data.places);
    //displaying paper information
    sendPaper(data.papers[0]);

    //show buttons
    document.getElementById('previous_btn').style.display = "inline";
    document.getElementById('next_btn').style.display = "inline";

    var focus = root,
        nodes = pack.nodes(root),
        view;

    var circle = svg.selectAll("circle")
        .data(nodes)
      .enter().append("circle")
        .attr("class", function(d) { return d.parent ? d.children ? "node" : "node node--leaf" : "node node--root"; })
        .style("fill", function(d) { return d.children ? color(d.depth) : colourNode(d, firstpmid); })
        .on("click", function(d) { if (focus !== d) zoom(d), d3.event.stopPropagation(); });

    var text = svg.selectAll("text")
        .data(nodes)
      .enter().append("text")
        .attr("class", "label")
        .style("fill-opacity", function(d) { return d.parent === root ? 1 : 0; })
        .style("display", function(d) { return d.parent === root ? null : "none"; })
        .text(function(d) { return d.name; });

    var node = svg.selectAll("circle,text");

    d3.select("body")
        //.style("background", color(-1))
        .on("click", function() { zoom(root); });

    zoomTo([root.x, root.y, root.r * 2 + margin]);
    document.getElementById("loading").style.display = "none";

  }

  function zoom(d) {
    var focus0 = focus; focus = d;

    var transition = d3.transition()
        .duration(d3.event.altKey ? 7500 : 750)
        .tween("zoom", function(d) {
          var i = d3.interpolateZoom(view, [focus.x, focus.y, focus.r * 2 + margin]);
          return function(t) { zoomTo(i(t)); };
        });

    transition.selectAll("text")
      .filter(function(d) { return d.parent === focus || this.style.display === "inline"; })
        .attr("y", -26)
        .style("fill-opacity", function(d) { return d.parent === focus ? 1 : 0; })
        .each("start", function(d) { if (d.parent === focus) this.style.display = "inline"; })
        .each("end", function(d) { if (d.parent !== focus) this.style.display = "none"; });
  }

  function zoomTo(v) {
    var k = diameter / v[2]; view = v;
    node.attr("transform", function(d) { return "translate(" + (d.x - v[0]) * k + "," + (d.y - v[1]) * k + ")"; });
    circle.attr("r", function(d) { return d.r * k; });
  }
}

// Load map and place D3 SVG
google.maps.event.addDomListener(window, 'load', initialize);
d3.select(self.frameElement)
  .style("height", diameter + "px")
  .style("stroke", "black")
  .style("fill", "none")
  .style("stroke-width", "1");
