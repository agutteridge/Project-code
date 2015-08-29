// initialize Google Maps
function initialize() {
  var mapOptions = {
    // centres on London
    center: { lat: 51.508742, lng: -0.120850},
    zoom: 2,
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
    marker = new google.maps.Marker({
      position: new google.maps.LatLng(data[i].place.geometry.location.lat, 
        data[i].place.geometry.location.lng),
      map: map,
      pmid: data[i].PMID
    });

    // add marker to global array
    markers.push(marker)

    google.maps.event.addListener(marker, 'click', (function(marker, i) {
      return function() {
        infowindow.setContent(marker.pmid);
        infowindow.open(map, marker);
      }
    })(marker, i));
  }

  var markerCluster = new MarkerClusterer(map, markers);
}

function clearMarkers() {
  for (var i = 0; i < markers.length; i++ ) {
    markers[i].setMap(null);
  }
  markers.length = 0;
}

// global vars
var markers = [];
var map;
var papers;
var text_field = document.getElementById('text_field');
var btn_search = document.getElementById('btn_search');

btn_search.onclick = function(){
  console.log(text_field.value);
  document.getElementById("loading").style.display = "block";
  d3.json("/data?search_term=" + text_field.value, d3_callback);
  svg.selectAll("circle,text").remove();
};

//D3: mostly taken from http://bl.ocks.org/mbostock/7607535
var margin = 20,
    diameter = 600;

var color = d3.scale.linear()
    .domain([-1, 5])
    .range(["hsl(152,80%,80%)", "hsl(228,30%,40%)"])
    .interpolate(d3.interpolateHcl);

var pack = d3.layout.pack()
    .padding(5)
    .size([diameter - margin, diameter - margin])
    .value(function(d) { return d.PMIDs.length * 500; })

var svg = d3.select("body").append("svg")
    .attr("width", diameter)
    .attr("height", diameter)
  .append("g")
    .attr("transform", "translate(" + diameter / 2 + "," + diameter / 2 + ")");

function d3_callback(error, data) {
  console.log('in d3_callback')
  if (error) throw error;

  var root = data.concepts;
  // passing location data to Google Maps function
  placeMarkers(data.places)
  // assigning global var papers
  papers = data.papers

  var focus = root,
      nodes = pack.nodes(root),
      view;

  var circle = svg.selectAll("circle")
      .data(nodes)
    .enter().append("circle")
      .attr("class", function(d) { return d.parent ? d.children ? "node" : "node node--leaf" : "node node--root"; })
      .style("fill", function(d) { return d.children ? color(d.depth) : null; })
      .on("click", function(d) { if (focus !== d) zoom(d), d3.event.stopPropagation(); });

// function(d) { return d.name }

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
d3.select(self.frameElement).style("height", diameter + "px");
