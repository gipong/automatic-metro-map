import * as d3 from "./js/d3.v3.min";
import * as numeric from "./js/numeric";
import * as turf from "./js/turf.min";
import queue from "./js/queue.v1.min";
import * as Coordsolver from "./js/coordsolver";
import "./js/d3-polygon.v1.min";
import "./js/d3-voronoi.v1.min";

import "./main.css"

// layout setting
var w = 800;
var h = 800;

// map projection setting
var projection = d3.geo.mercator().scale(1).translate([0, 0]);
var path = d3.geo.path().projection(projection);
//var zoom = d3.behavior.zoom()
//    .translate([0, 0])
//    .scale(1).scaleExtent([1, 8])
//    .on("zoom", function() {
//        svg.attr("transform", `translate(${d3.event.translate})scale(${d3.event.scale})`);
//    });

// svg setting
var tooltip = d3.select("body").append("div").attr("class", "tooltip").style("opacity", 0);
var svg = d3.select("body").append("svg")
    .attr({"width": w, "height": h, "class": "framed", "id": "svg"});

var g = svg.append("g");

//svg.call(zoom).call(zoom.event);

var drag = d3.behavior.drag()
    .origin(function(d) {return d;})
    .on("dragstart", dragstartFun)
    .on("drag", dragFun)
    .on("dragend", dragendFun);

function dragstartFun() {
//        d3.event.sourceEvent.stopPropagation();
}

function dragFun(d) {
    let coord = d3.mouse(this);
    d3.select(this)
        .attr("cx", coord[0])
        .attr("cy", coord[1]);

    distorted[d.properties.id][0] = coord[0];
    distorted[d.properties.id][1] = coord[1];

    redraw();
}

function dragendFun() {
    d3.select(this).classed("dragging", false);
}

var topo, nodes, topoGeo, topoTempGeo, TIN, TINi, TINic, TINinv, predistorted, edgeLength;
var origin = [], distorted = [], lines = [], linec = [], linecDistorted = [], net = [], edgeSet = [], svgList = [], svgIdx = [];

// routes open data, http://data.taipei/opendata/datalist/datasetMeta?oid=afccd2ac-75b1-4362-9099-45983e332776
// stations open data, http://data.taipei/opendata/datalist/datasetMeta?oid=758e5ae0-e6ee-448b-81f5-316eb68a5ba7
queue()
    .defer(d3.json, "geodata/TpeMRTStations_WGS84_2011.geojson") //TpeMRTStations_WGS84_2011.geojson
    .defer(d3.json, "geodata/Tpe_boundary.geojson")
    .defer(d3.json, "geodata/topology_2011.json")
    .await(metroMap);

// line 24, [[55, 18], [18, 52]]
// line 26, [[93, 98], [98, 27]]
// line 71, [[54, 50], [50, 55]]

function metroMap(error, stations, tpe, topology) {
    if(error) throw error;

    topo = topology;
    nodes = stations;

    let b = path.bounds(tpe.features[0]);
    let s = .95/Math.max((b[1][0] - b[0][0])/w, (b[1][1] - b[0][1])/h);
    let t = [(w - s*(b[1][0] + b[0][0]))/2, (h - s*(b[1][1] + b[0][1]))/2];

    projection
        .scale(s)
        .translate(t);

    topoGeo = [];
    topo.forEach(function(e){
        e.topology.forEach(function(c, ci){
            let coords = [];
            for(let i=0; i< c.length; i++) {
                coords.push(nodes.features[c[i][0]].geometry.coordinates);
            }
            coords.push(nodes.features[c.slice(-1)[0][1]].geometry.coordinates);
            topoGeo.push({
                "geometry": {"type": "LineString", "coordinates": coords}, "type": "Feature",
                "properties": {
                    "lineName": e.lineName,
                    "lineName_en": e.lineName_en,
                    "c": e.c,
                    "ci": ci
                }
            })
        });

    });

    g.selectAll("path .line").data(topoGeo).enter().append("path")
        .attr({
            "class": "line",
            "id": function(d) {return `${d.properties.lineName_en}_${d.properties.ci}`;},
            "d": path,
            "stroke": function(d) {return d.properties.c;},
            "stroke-width": 3,
            "stroke-linejoin": "round",
            "fill": "none"
        });

    g.selectAll("circle").data(nodes.features).enter().append("circle")
        .attr({
            "cx": function(d) {return projection(d.geometry.coordinates)[0];},
            "cy": function(d) {return projection(d.geometry.coordinates)[1];},
            "id": function(d) {return d.properties.id;},
            "r": 3,
            "stroke": "black",
            "stroke-width": 2,
            "fill": "white",
        })
        .call(drag)
        .on("mouseover", function(d) {
            tooltip.style("opacity",.9).text(d.properties.NAME+d.properties.id);
        })
        .on("mouseout", function() {
            tooltip.style("opacity", 0);
        });

    g.selectAll("text").data(nodes.features).enter().append("text")
        .attr({
            "x": function(d) {return projection(d.geometry.coordinates)[0]+5;},
            "y": function(d) {return projection(d.geometry.coordinates)[1]+5;},
            "font-size": 8,
            "fill": "#555",
        })
        .text(function(d) { return d.properties.id;});

    nodes.features.forEach(function(e) {
        let proCoords = projection(e.geometry.coordinates);
        origin.push(proCoords);
        distorted.push(proCoords);
    });

    edgeLength = d3.mean(
        d3.selectAll("path")[0]
            .map(function(d){return d.getTotalLength();})
    );

    edgeSet = [];
    topo.forEach(function(t) {
        edgeSet.push(t.topology.reduce(function(a, b) {return b.concat(a);}));
    });
    edgeSet = edgeSet.reduce(function(a, b) {return b.concat(a);});

    for (let i=0; i<distorted.length; i++) {
        net.push([]);
    }

    for(let id=0; id<edgeSet.length; id++) {
        let i = edgeSet[id][0];
        let j = edgeSet[id][1];

        net[i].push(j);
        net[j].push(i);
    }

    TINmodel();
    // dev
    //reset2optimal();
}

function conjugateGradient(A, b, x) {
    // Ax = b
    //let tol = (typeof tol != 'undefined')?tol:1e-5;
    let r = numeric.sub(b, numeric.dotMV(A, x));
    let p = r.slice(0);

    let rsold = numeric.dotVV(r, r);
    for(let i=0; i< r.length; i++) {
        let Ap = numeric.dotMV(A, p);
        let alpha = rsold/numeric.dotVV(p, Ap);
        x = numeric.add(x, numeric.mul(p, numeric.linspace(alpha, alpha, p.length)));
        r = numeric.sub(r, numeric.mul(Ap, numeric.linspace(alpha, alpha, Ap.length)));
        let rsoldN = numeric.dotVV(r, r);

        if(Math.sqrt(rsoldN) < 1e-10) break;
        p = numeric.add(r, numeric.mul(p, numeric.linspace(rsoldN/rsold, rsoldN/rsold, p.length)));
        rsold = rsoldN;
    }
    return x;
}

function redraw() {

    let x = [],
        y = [];

    for(let i=0; i<distorted.length; i++) {
        x.push(distorted[i][0]);
        y.push(distorted[i][1]);
    }

    let xscale = d3.scale.linear()
        .domain([d3.min(x), d3.max(x)])
        .range([100, 700]);

    let yscale = d3.scale.linear()
        .domain([d3.min(y), d3.max(y)])
        .range([100, 700]);

    for(let i=0; i<distorted.length; i++) {
        distorted[i][0] = xscale(distorted[i][0]);
        distorted[i][1] = yscale(distorted[i][1]);
    }

    g.selectAll("circle").transition()
        .attr({
            "cx": function(d, i) {return distorted[i][0];},
            "cy": function(d, i) {return distorted[i][1];}
        });

    g.selectAll("text").transition()
        .attr({
            "x": function(d, i) {return distorted[i][0]+5;},
            "y": function(d, i) {return distorted[i][1]+5;},
            "font-size": 8,
            "fill": "#555",
        })
        .text(function(d) { return d.properties.id;});

    let topoTemp = [];
    topoTempGeo = [];
    topo.forEach(function(e){
        e.topology.forEach(function(c, ci){
            let coords = [];
            for(let i=0; i< c.length; i++) {
                coords.push(distorted[c[i][0]]);
            }
            coords.push(distorted[c.slice(-1)[0][1]]);
            topoTemp.push(coords);
            coords = coords.map(function(d) {return projection.invert(d);});
            topoTempGeo.push({
                "geometry": {"type": "LineString", "coordinates": coords}, "type": "Feature",
                "properties": {
                    "lineName": e.lineName,
                    "lineName_en": e.lineName_en,
                    "c": e.c,
                    "ci": ci
                }
            });
        });

    });

    g.selectAll(".line").data(topoTemp).transition()
        .attr("d", d3.svg.line()
            .x(function(d) {return d[0];})
            .y(function(d) {return d[1];})
        );

}

// using numeric.js quadratic solver
function optimize() {
    let mat = setVarNumber(distorted.length*2);
    let Amat = [];
    let bvec = [];
    //let deg45 = Math.PI/4;

    for(let id=0; id<edgeSet.length; id++) {
        let i = edgeSet[id][0];
        let j = edgeSet[id][1];

        let vec = numeric.sub(distorted[i], distorted[j]);
        let psd_vec = [0, edgeLength];

        let theta = calAngle(psd_vec, vec);

        let rMat = [
            [Math.cos(theta), -Math.sin(theta)],
            [Math.sin(theta), Math.cos(theta)]
        ];
        //let b = numeric.dotMV(rMat, vec);

        //Coordsolver.addSoftConstraint(mat[1], mat[2],
        //        `vars[${i}]-vars[${j}]-${b[0]}`, 0.1);
        //Coordsolver.addSoftConstraint(mat[1], mat[2],
        //        `vars[${i+distorted.length}]-vars[${j+distorted.length}]-${b[1]}`, 0.1);

        Coordsolver.addSoftConstraint(mat[1], mat[2],
            `vars[${i}]-vars[${j}]-${Math.sqrt(Math.abs(distorted[i][0]-distorted[j][0]))}`, 1);
        Coordsolver.addSoftConstraint(mat[1], mat[2],
            `vars[${i+distorted.length}]-vars[${j+distorted.length}]-${Math.sqrt(Math.abs(distorted[i][1]-distorted[j][1]))}`, 1);


        //Coordsolver.addSoftConstraint(mat[1], mat[2],
        //        `vars[${i}]-vars[${j}]-${20}`);
        //Coordsolver.addSoftConstraint(mat[1], mat[2],
        //        `vars[${i+distorted.length}]-vars[${j+distorted.length}]-${20}`);
        //
        //Coordsolver.addSoftConstraint(mat[1], mat[2], `vars[${i}]-vars[${j}]+vars[${i+distorted.length}]-vars[${j+distorted.length}]-`+(b[0]-b[1]));
    }


    let passSet = [6, 15, 36, 52, 54, 55, 56, 71];
    for(let u=0; u<net.length; u++) {
        let o = 360/net[u].length;
        let leni = +(( Math.tan( ((180-o)/2)*Math.PI/180 ) -1 )/2).toFixed(3);
        let leniSign = leni>0?"+":"-";
        let leni1 = +((1+Math.tan( ((180-o)/2)*Math.PI/180 ))/2).toFixed(3);
        let leni1Sign = leni1>0?"-":"+";

        if(net[u].length == 1) {
            // boundary point
            //let p = net[net[u][0]];
            //Coordsolver.addSoftConstraint(mat[1], mat[2],
            //        `vars[${p[0]}]-vars[${net[u][0]}]`, 1);
            //Coordsolver.addSoftConstraint(mat[1], mat[2],
            //        `vars[${p[0]+distorted.length}]-vars[${net[u][0]+distorted.length}]`, 1);

        }

        for (let v=0; v<net[u].length-1; v++) {
            //let ts = (passSet.indexOf(u)>0)
            Coordsolver.addSoftConstraint(mat[1], mat[2],
                `vars[${net[u][v]}]-2vars[${u}]+vars[${net[u][v+1]}]`, 0.1);
            Coordsolver.addSoftConstraint(mat[1], mat[2],
                `vars[${net[u][v]+distorted.length}]-2vars[${u+distorted.length}]+vars[${net[u][v+1]+distorted.length}]`, 0.1);

            let di = net[u][v];
            //let dj = net[u][v+1];
            //Coordsolver.addSoftConstraint(mat[1], mat[2],
            //        `vars[${net[u][v]}]-vars[${u}]-${distorted[di][0]-distorted[u][0]}`, 0.05);
            //Coordsolver.addSoftConstraint(mat[1], mat[2],
            //        `vars[${net[u][v]+distorted.length}]-vars[${u+distorted.length}]-${distorted[di][1]-distorted[u][1]}`, 0.05);
            //
            //Coordsolver.addSoftConstraint(mat[1], mat[2],
            //        `vars[${net[u][v+1]}]-vars[${u}]-${distorted[dj][0]-distorted[u][0]}`, 0.1);
            //Coordsolver.addSoftConstraint(mat[1], mat[2],
            //        `vars[${net[u][v+1]+distorted.length}]-vars[${u+distorted.length}]-${distorted[dj][1]-distorted[u][1]}`, 0.1);

            if(net[u].length>2) {
                Coordsolver.addSoftConstraint(mat[1], mat[2],
                    `vars[${u}]${leniSign+Math.abs(leni)}vars[${net[u][v]}]${leni1Sign+Math.abs(leni1)}vars[${net[u][v+1]}]`, 0.5);
                Coordsolver.addSoftConstraint(mat[1], mat[2],
                    `vars[${u+distorted.length}]${leniSign+Math.abs(leni)}vars[${net[u][v]+distorted.length}]${leni1Sign+Math.abs(leni1)}vars[${net[u][v+1]+distorted.length}]`, 0.5);

                //Coordsolver.addSoftConstraint(mat[1], mat[2],
                //        `vars[${net[u][v]}]-`+distorted[net[u][v]][0], 0.01);
                //Coordsolver.addSoftConstraint(mat[1], mat[2],
                //        `vars[${net[u][v]}]-`+distorted[net[u][v]][1], 0.01);
                //Coordsolver.addSoftConstraint(mat[1], mat[2],
                //        `vars[${net[u][v+1]}]-`+distorted[net[u][v+1]][0], 0.01);
                //Coordsolver.addSoftConstraint(mat[1], mat[2],
                //        `vars[${net[u][v+1]}]-`+distorted[net[u][v+1]][1], 0.01);
                //
                //let si = 1;
                //Coordsolver.addSoftConstraint(mat[1], mat[2],
                //        `vars[${u}]-vars[${net[u][v]}]-${si}`);
                //Coordsolver.addSoftConstraint(mat[1], mat[2],
                //        `vars[${u}]-vars[${net[u][v+1]}]-${si}`);
                //Coordsolver.addSoftConstraint(mat[1], mat[2],
                //        `vars[${u+distorted.length}]-vars[${net[u][v]+distorted.length}]-${si}`);
                //Coordsolver.addSoftConstraint(mat[1], mat[2],
                //        `vars[${u+distorted.length}]-vars[${net[u][v+1]+distorted.length}]-${si}`);

                Coordsolver.addSoftConstraint(mat[1], mat[2],
                    `vars[${u}]-vars[${di}]-${2*Math.abs(distorted[u][0]-distorted[di][0])/edgeLength}`, 0.5);
                Coordsolver.addSoftConstraint(mat[1], mat[2],
                    `vars[${u+distorted.length}]-vars[${di+distorted.length}]-${2*Math.abs(distorted[u][1]-distorted[di][1])/edgeLength}`, 0.5);

            }
        }
    }

    for(let i=0; i<distorted.length*2; i++) {
        Coordsolver.addHardConstraint(Amat, bvec, mat[1], `vars[${i}]>0`);
    }

    //for(let i=0; i<distorted.length; i++) {
    //    for(let j=0; j<distorted.length; j++) {
    //        if(i!=j) {
    //            Coordsolver.addSoftConstraint(mat[1], mat[2], `vars[${i}]-vars[${j}]-${edgeLength}`, 0.01);
    //            Coordsolver.addSoftConstraint(mat[1], mat[2], `vars[${i+distorted.length}]-vars[${j+distorted.length}]-${edgeLength}`, 0.01);
    //        }
    //    }
    //}

    for(let i=0; i<distorted.length; i++) {
        Coordsolver.addSoftConstraint(mat[1], mat[2], `vars[${i}]-`+distorted[i][0], 0.05);
        Coordsolver.addSoftConstraint(mat[1], mat[2], `vars[${i+distorted.length}]-`+distorted[i][1], 0.05);
    }

    for(let ti=0; ti<TINinv.length; ti++) {
        let ntiarr = [...TINinv[ti]];
        //ntiarr.push(ntiarr[0]);

        let x = [], y = [];
        for(let ni=0; ni<ntiarr.length; ni++) {
            x.push(distorted[ni][0]);
            y.push(distorted[ni][1]);
        }

        let xmin = Math.min(...x);
        let xmax = Math.max(...x);
        let ymin = Math.min(...y);
        let ymax = Math.max(...y);

        //Coordsolver.addHardConstraint(Amat, bvec, mat[1], `vars[${ti}]>${xmin}`);
        //Coordsolver.addHardConstraint(Amat, bvec, mat[1], `vars[${ti}]>${-xmax}`);
        //Coordsolver.addHardConstraint(Amat, bvec, mat[1], `vars[${ti+distorted.length}]>${ymin}`);
        //Coordsolver.addHardConstraint(Amat, bvec, mat[1], `vars[${ti+distorted.length}]>${-ymax}`);

        Coordsolver.addSoftConstraint(mat[1], mat[2], `vars[${ti}]-`+(xmin+xmax)/2, 0.01);
        Coordsolver.addSoftConstraint(mat[1], mat[2], `vars[${ti+distorted.length}]-`+(ymin+ymax)/2, 0.01);

        Coordsolver.addSoftConstraint(mat[1], mat[2], `2vars[${ti}]-vars[${x[x.indexOf(xmin)]}]-vars[${x[x.indexOf(xmax)]}]`, 0.01);
        Coordsolver.addSoftConstraint(mat[1], mat[2], `2vars[${ti+distorted.length}]-vars[${y[y.indexOf(ymin)]+distorted.length}]-vars[${y[y.indexOf(ymax)+distorted.length]}]`, 0.01);
    }


    let result = numeric.solveQP(mat[1], mat[2], numeric.transpose(Amat), bvec);

    predistorted = distorted.slice(0);

    for(let i=0; i<distorted.length; i++) {
        distorted[i] = [result.solution[i], result.solution[i+distorted.length]];
    }

    redraw();

}

// using numeric.js quadratic solver
function optimizeSmooth() {
    let mat = setVarNumber(distorted.length*2);
    let Amat = [];
    let bvec = [];
    let skipNode = net.filter(d=>{return d.length>2;}).reduce((a, b)=>{return a.concat([...b])}, []);

    for(let id=0; id<skipNode.length; id++) {
        if(net[skipNode[id]].length < 2) continue;
        let i = net[skipNode[id]][0];
        let j = net[skipNode[id]][1];

        Coordsolver.addSoftConstraint(mat[1], mat[2],
            `vars[${skipNode[id]}]-${(distorted[i][0]+distorted[j][0])*0.5}`, 0.5);
        Coordsolver.addSoftConstraint(mat[1], mat[2],
            `vars[${skipNode[id]+distorted.length}]-${(distorted[i][1]+distorted[j][1])*0.5}`, 0.5);

        Coordsolver.addSoftConstraint(mat[1], mat[2],
            `2vars[${skipNode[id]}]-vars[${i}]-vars[${j}]`, 0.05);
        Coordsolver.addSoftConstraint(mat[1], mat[2],
            `2vars[${skipNode[id]+distorted.length}]-vars[${i}]-vars[${j}]`, 0.05);
    }



    for(let i=0; i<distorted.length; i++) {
        //let w = (net[i].length > 2)?5:1;
        let w = 1;
        Coordsolver.addSoftConstraint(mat[1], mat[2], `vars[${i}]-`+distorted[i][0], w);
        Coordsolver.addSoftConstraint(mat[1], mat[2], `vars[${i+distorted.length}]-`+distorted[i][1], w);
    }

    for(let i=0; i<distorted.length*2; i++) {
        Coordsolver.addHardConstraint(Amat, bvec, mat[1], `vars[${i}]>0`);
    }

    let result = numeric.solveQP(mat[1], mat[2], numeric.transpose(Amat), bvec);

    predistorted = distorted.slice(0);

    for(let i=0; i<distorted.length; i++) {
        distorted[i] = [result.solution[i], result.solution[i+distorted.length]];
    }

    redraw();

}

function setVarNumber(number) {
    let vars = [];
    let Dmat = [];
    let dvec = [];
    for(let i=0; i<number; i++) {
        vars.push("v"+i);
        Dmat.push(Array.apply(null, Array(number)).map(Number.prototype.valueOf,0));
        dvec.push(0);
    }
    return [vars, Dmat, dvec];
}

function calAngle(v1, v2) {
    return Math.atan2(
        (v1[0]*v2[1] - v1[1]*v2[0]),
        v1[0]*v2[0] + v1[1]*v2[1]
    );
}

function optimialcg(netNum) {
    netNum = (typeof netNum !== "undefined")?netNum:2;
    let weight = {
        "length": 5,
        "maxAngle": 2.5
    };
    let x = [];
    let y = [];

    for(let i=0; i<distorted.length; i++) {
        x.push(distorted[i][0]);
        y.push(distorted[i][1]);
    }

    //let deg45 = Math.PI/4;
    let A = [];
    let bx = [];
    let by = [];

    for(let id=0; id<edgeSet.length; id++) {
        let i = edgeSet[id][0];
        let j = edgeSet[id][1];

        let rMat = [
            [Math.cos(0), -Math.sin(0)],
            [Math.sin(0), Math.cos(0)]
        ];

        let vec = numeric.sub(distorted[i], distorted[j]);
        let L  = Math.sqrt(numeric.norm2Squared(numeric.sub(distorted[i], distorted[j])));

        let b = numeric.dotMV(numeric.mul(1*(edgeLength/L), rMat), vec);

        let r = new Array(distorted.length).fill(0);

        let w = weight["length"];
        r[i] = w;
        r[j] = -w;

        A.push(r);
        bx.push(w*b[0]);
        by.push(w*b[1]);

    }

    let skipNode = net.filter(d=>{return d.length>2;}).reduce((a, b)=>{return a.concat([...b])}, []);

    for(let i=0; i<distorted.length; i++) {
        if(net[i].length < netNum)
            continue;

        let theta = 2*Math.PI/net[i].length;
        let tanCoef = Math.tan((Math.PI - theta)*0.5);

        let w = (skipNode.indexOf(i) > 0 && net[i].length < 3)?1:weight["maxAngle"];

        for(let j=0; j<net[i].length; j++) {
            let r = new Array(distorted.length).fill(0);
            let k = (j+1) % net[i].length;
            let idx_j = net[i][j];
            let idx_k = net[i][k];
            //console.log(`${i}, ${idx_j}, ${idx_k}`);
            //
            //w = weight["maxAngle"];

            r[i] = 1*w;
            r[idx_j] = (-1 + 0.5 + tanCoef * 0.5)*w;
            r[idx_k] = (-0.5 - tanCoef * 0.5)*w;

            A.push(r);
            bx.push(0);
            by.push(0);
        }
    }

    let k = 0;
    let iterNum = 100;
    let Ax = A.slice(0);
    let Ay = A.slice(0);
    let engterm =
        numeric.norm2Squared(numeric.sub(numeric.dotMV(Ax,x),bx))+
        numeric.norm2Squared(numeric.sub(numeric.dotMV(Ay,y),by));
    let pre_engterm = engterm;
    let diff =  (pre_engterm-engterm)/pre_engterm;

    do {
        let ATx = numeric.transpose(Ax);
        let ATy = numeric.transpose(Ay);
        x = conjugateGradient(numeric.dotMMbig(ATx, Ax), numeric.dotMV(ATx, bx), x);
        y = conjugateGradient(numeric.dotMMbig(ATy, Ay), numeric.dotMV(ATy, by), y);

        engterm =
            numeric.norm2Squared(numeric.sub(numeric.dotMV(Ax,x),bx))+
            numeric.norm2Squared(numeric.sub(numeric.dotMV(Ay,y),by));

        k++;

        diff = (pre_engterm-engterm)/pre_engterm;
    }while (k<iterNum && diff>1e-5);

    predistorted = distorted.slice(0);

    for(let i=0; i<distorted.length; i++) {
        distorted[i] = [x[i], y[i]];
    }

    redraw();

}

function octilinearity() {
    let x = [];
    let y = [];

    for(let i=0; i<distorted.length; i++) {
        x.push(distorted[i][0]);
        y.push(distorted[i][1]);
    }

    let deg45 = Math.PI/4;
    let A = [];
    let bx = [];
    let by = [];

    for(let id=0; id<edgeSet.length; id++) {
        let i = edgeSet[id][0];
        let j = edgeSet[id][1];

        let vec = numeric.sub(distorted[i], distorted[j]);
        let psd_vec = [0, 1];

        let diffAngle = calAngle(psd_vec, vec);
        let theta = (Math.round(diffAngle/deg45)*deg45)-diffAngle;

        let rMat = [
            [Math.cos(theta), -Math.sin(theta)],
            [Math.sin(theta), Math.cos(theta)]
        ];


        let L  = Math.sqrt(numeric.norm2Squared(numeric.sub(distorted[i], distorted[j])));
        let s = (net[i].length>=2 || net[j].length>=2)?edgeLength/L:1;
        //if(edgeLength > Math.sqrt(numeric.norm2Squared(vec)))
        //    vec = numeric.mul(edgeLength/Math.sqrt(numeric.norm2Squared(vec)), vec);
        //
        //vec = numeric.mul(s, vec);
        let arr = calDensity();
        let avg = (d3.max(arr)+d3.min(arr))/2;
        if(arr[i] < avg || arr[j] < avg) {
            s = avg/arr[i];
        } else {
            s = arr[i]/avg;
        }

        s = (net[i].length>=2 || net[j].length>=2)?edgeLength/L:1;
        if(s< 0.2*edgeLength/L) s = 0.2*edgeLength/L;

        //vec = numeric.mul(s, vec);
        let b = numeric.mul(1, numeric.dotMV(rMat, vec));
        let r = new Array(distorted.length).fill(0);

        let w = 10;
        r[i] = w;
        r[j] = -w;

        //if(Math.abs(theta) < 0.1) continue;

        A.push(r);
        bx.push(w*b[0]);
        by.push(w*b[1]);
    }

    for(let i=0; i<distorted.length; i++) {
        if(net[i].length < 2)
            continue;

        let theta = 2*Math.PI/net[i].length;
        let tanCoef = Math.tan((Math.PI - theta)*0.5);


        for(let j=0; j<net[i].length; j++) {
            let r = new Array(distorted.length).fill(0);
            let k = (j+1) % net[i].length;
            let idx_j = net[i][j];
            let idx_k = net[i][k];
            //console.log(`${i}, ${idx_j}, ${idx_k}`);

            let w = 1;
            if(net[i].length > 2) w = 5;
            w = 0.5;

            r[i] = 1*w;
            r[idx_j] = (-1 + 0.5 + tanCoef * 0.5)*w;
            r[idx_k] = (-0.5 - tanCoef * 0.5)*w;

            A.push(r);
            bx.push(0);
            by.push(0);
        }
    }

    for(let ti=0; ti<TINinv.length; ti++) {
        //let ntiarr = [...TINinv[ti]];
        //ntiarr.push(ntiarr[0]);

        //let x = [], y = [];
        //for(let ni=0; ni<ntiarr.length; ni++) {
        //    x.push(distorted[ni][0]);
        //    y.push(distorted[ni][1]);
        //}
        //
        //let xmin = Math.min(...x);
        //let xmax = Math.max(...x);
        //let ymin = Math.min(...y);
        //let ymax = Math.max(...y);
        //
        //let r = new Array(distorted.length).fill(0);
        //r[ti] = 1;
        //r[x[x.indexOf(xmin)]] = -1;
        ////r[x[x.indexOf(xmax)]] = -0.5;
        //
        //console.log(ti, x.indexOf(xmin), x.indexOf(xmax));
        //A.push(r);
        //bx.push(edgeLength);
        //
        //r = new Array(distorted.length).fill(0);
        //r[ti] = 1;
        //r[y[y.indexOf(ymin)]] = -1;
        ////r[y[y.indexOf(ymax)]] = -0.5;
        //by.push(0);
    }


    for(let id=0; id<distorted.length; id++) {
        let r = new Array(distorted.length).fill(0);
        r[id] = 1;
        A.push(r);
        bx.push(distorted[id][0]);
        by.push(distorted[id][1]);
    }

    let k = 0;
    let iterNum = 100;
    let Ax = A.slice(0);
    let Ay = A.slice(0);
    let engterm =
        numeric.norm2Squared(numeric.sub(numeric.dotMV(Ax,x),bx))+
        numeric.norm2Squared(numeric.sub(numeric.dotMV(Ay,y),by));
    let pre_engterm = engterm;
    let diff =  (pre_engterm-engterm)/pre_engterm;

    do {
        let ATx = numeric.transpose(Ax);
        let ATy = numeric.transpose(Ay);
        x = conjugateGradient(numeric.dotMMbig(ATx, Ax), numeric.dotMV(ATx, bx), x);
        y = conjugateGradient(numeric.dotMMbig(ATy, Ay), numeric.dotMV(ATy, by), y);

        engterm =
            numeric.norm2Squared(numeric.sub(numeric.dotMV(Ax,x),bx))+
            numeric.norm2Squared(numeric.sub(numeric.dotMV(Ay,y),by));

        k++;

        diff = (pre_engterm-engterm)/pre_engterm;
    }while (k<iterNum && diff>1e-4);

    predistorted = distorted.slice(0);

    for(let i=0; i<distorted.length; i++) {
        distorted[i] = [x[i], y[i]];
    }

    redraw();
}

function laplacianDeformation(p, svgUrl) {
    let weight = {
        "length": 1
    };
    let x = [];
    let y = [];

    for(let i=0; i<distorted.length; i++) {
        x.push(distorted[i][0]);
        y.push(distorted[i][1]);
    }

    let preL = 0;
    let A = [];
    let bx = [];
    let by = [];

    //for(let ti=0; ti<TINinv.length; ti++) {
    //    let ntiarr = [...TINinv[ti]];
    //
    //    console.log(ntiarr);
    //}

    for(let id=0; id<edgeSet.length; id++) {
        let i = edgeSet[id][0];
        let j = edgeSet[id][1];

        let rMat = [
            [Math.cos(0), -Math.sin(0)],
            [Math.sin(0), Math.cos(0)]
        ];

        let vec = numeric.sub(distorted[i], distorted[j]);
        let L  = Math.sqrt(numeric.norm2Squared(numeric.sub(distorted[i], distorted[j])));

        let enl = 1;
        if(i == p || j == p) {
            let psd_vec = [0, 1];
            let vec;
            if(i == p) { vec = numeric.sub(distorted[i], distorted[j]); }
            if(j == p) { vec = numeric.sub(distorted[j], distorted[i]); }
            let diffAngle = calAngle(psd_vec, vec)*180/Math.PI;

            if(diffAngle > -45 && diffAngle < 135) {
                if(i == p) { preL = j; }
                if(j == p) { preL = i; }
                enl = 4.5;
            } else {
                enl = 1;
            }
        }

        let b = numeric.dotMV(numeric.mul(enl*(edgeLength/L), rMat), vec);

        let r = new Array(distorted.length).fill(0);

        let w = weight["length"];
        r[i] = w;
        r[j] = -w;

        A.push(r);
        bx.push(w*b[0]);
        by.push(w*b[1]);

    }


    //for(let i=0; i<distorted.length; i++) {
    let r = new Array(distorted.length).fill(0);
    let b = [0, 0];

    let calCotList = TINinv[p].concat(TINinv[p].slice(0, TINinv[p].length%3));
    for(let idx = 1; idx<calCotList.length-1; idx++) {
        //console.log(calCotList[idx]);
        r[calCotList[idx]] = ptDegree(
            distorted[p],
            distorted[calCotList[idx]],
            distorted[calCotList[idx-1]],
            distorted[calCotList[idx+1]]
        );
        b[0] += distorted[calCotList[idx]][0];
        b[1] += distorted[calCotList[idx]][1];
    }

    //TINinv[p].forEach(function(t) {
    //    console.log(TINinv[p]);
    //    r[t] = 1/TINinv[p].length;
    //    b[0] += distorted[t][0];
    //    b[1] += distorted[t][1];
    //    console.log(b);
    //});

    r[p] = -1;
    A.push(r);
    bx.push(b[0]/TINinv[p].length-distorted[p][0]);
    by.push(b[1]/TINinv[p].length-distorted[p][1]);
    //}

    //for(let i=0; i<distorted.length; i++) {
    //    let r = new Array(distorted.length).fill(0);
    //    r[i] = 1;
    //
    //    A.push(r);
    //    bx.push(distorted[i][0]);
    //    by.push(distorted[i][1]);
    //
    //}

    let k = 0;
    let iterNum = 100;
    let Ax = A.slice(0);
    let Ay = A.slice(0);
    let engterm =
        numeric.norm2Squared(numeric.sub(numeric.dotMV(Ax,x),bx))+
        numeric.norm2Squared(numeric.sub(numeric.dotMV(Ay,y),by));
    let pre_engterm = engterm;
    let diff =  (pre_engterm-engterm)/pre_engterm;

    do {
        let ATx = numeric.transpose(Ax);
        let ATy = numeric.transpose(Ay);
        x = conjugateGradient(numeric.dotMMbig(ATx, Ax), numeric.dotMV(ATx, bx), x);
        y = conjugateGradient(numeric.dotMMbig(ATy, Ay), numeric.dotMV(ATy, by), y);

        engterm =
            numeric.norm2Squared(numeric.sub(numeric.dotMV(Ax,x),bx))+
            numeric.norm2Squared(numeric.sub(numeric.dotMV(Ay,y),by));

        k++;

        diff = (pre_engterm-engterm)/pre_engterm;
    }while (k<iterNum && diff>1e-5);

    predistorted = distorted.slice(0);

    for(let i=0; i<distorted.length; i++) {
        distorted[i] = [x[i], y[i]];
    }

    redraw();
    putSvg(p, preL, svgUrl);
}

function ptDegree(vi, vj, lp, rp) {
    return (calAngle(numeric.sub(lp, vi), numeric.sub(lp, vj)) + calAngle(numeric.sub(rp, vi), numeric.sub(rp, vj)))*0.5;
}

function putSvg(p, preL, svgUrl) {
    svgList.push([p, preL]);
    svgIdx.push(p);
    fetch(svgUrl)
        .then(function(response) { return response.text(); })
        .then(function(svg) {
            let t = /<g>(.*)<\/g>/g.exec(svg)[1];
            let L = 0.55*Math.sqrt(numeric.norm2Squared(numeric.sub(distorted[p], distorted[preL])));
            let s = L/parseInt(d3.select("svg").attr("width"));
            d3.select("svg").append("g")
                .attr("id", `svg${svgList.length-1}`)
                .attr("transform", function() {
                    return `translate(${distorted[p][0]}, ${distorted[p][1]-L}) scale(${s}) `;
                })
                .html(t);
        });

    updateSvg();
}

function updateSvg() {
    svgList.forEach(function(svg, idx) {
        let L = 0.55*Math.sqrt(numeric.norm2Squared(numeric.sub(distorted[svg[0]], distorted[svg[1]])));
        let s = L/parseInt(d3.select("svg").attr("width"));
        d3.select(`#svg${idx}`).transition()
            .attr("transform", function() {
                return `translate(${distorted[svg[0]][0]}, ${distorted[svg[0]][1]-L}) scale(${s}) `;
            });
    });
}

function calDiff() {
    let sum = 0;
    for(let i=0; i<distorted.length; i++) {
        sum += Math.sqrt(numeric.norm2Squared(numeric.sub(predistorted[i], distorted[i])));
    }
    return sum;
}

function calDensity() {
    let arr = [];
    for(let i=0; i<distorted.length; i++) {
        let sum = 0;
        for(let j=0; j<distorted.length; j++) {
            sum += Math.sqrt(numeric.norm2Squared(numeric.sub(distorted[j], distorted[i])));
        }
        arr.push(sum);
    }

    return arr;
}

function runOptimialProcess() {
    let p = new Promise(function(resolve, reject) {
        optimialcg();
    });

    let c = calDiff();
    let diff;
    //let i = 0;
    do {
        p.then(optimialcg());
        diff = c-calDiff();
        c = calDiff();
        //console.log(`iternum ${++i}, diff ${diff}`);
    } while(diff>5*1e-1);

    //console.log(`run the octilinear process ...`);
    //runOctProcess();
}

function runOctProcess() {
    let p = new Promise(function(resolve, reject) {
        octilinearity();
    });

    let c = calDiff();
    let diff;
    do {
        p.then(octilinearity());
        diff = Math.abs(c-calDiff());
        c = calDiff();
    } while(diff>10);

}

function TINmodel() {
    TIN = turf.tin(nodes);
    topo.forEach(function(t) {
        lines.push(flatten(t.topology).filter(function(e, i, f) {
            return i == f.indexOf(e);
        }));
    });
    TINi = [];
    TINic = [];

    TIN.features.forEach(function(t) {
        let coords = [];
        let lineCheck = new Set();
        for(let i=0; i<t.geometry.coordinates[0].length-1; i++) {
            let ei = findNodeIdx(t.geometry.coordinates[0][i]);
            let ej = findNodeIdx(t.geometry.coordinates[0][i+1]);
            coords.push([ei[0], ej[0]]);

            lineCheck.add(ei[1]);
            lineCheck.add(ej[1]);
        }
        if(lineCheck.size > 1) {
            let u = numeric.sub(
                nodes.features[coords[0][0]].geometry.coordinates,
                nodes.features[coords[0][1]].geometry.coordinates
            );
            let v = numeric.sub(
                nodes.features[coords[1][1]].geometry.coordinates,
                nodes.features[coords[1][0]].geometry.coordinates
            );
            linec.push(crossVec(u.concat(0), v.concat(0)));
            TINi.push(coords);
        }
        TINic.push(coords);
    });

    TINinv = [];
    for(let i=0; i<distorted.length; i++) {
        TINinv[i] = new Set();
    }
    // reverse index

    TINic.forEach(function(e, i) {
        TINinv[e[0][0]].add(TINic[i][0][0]).add(TINic[i][1][0]).add(TINic[i][2][0]).delete(e[0][0]);
        TINinv[e[1][0]].add(TINic[i][0][0]).add(TINic[i][1][0]).add(TINic[i][2][0]).delete(e[1][0]);
        TINinv[e[2][0]].add(TINic[i][0][0]).add(TINic[i][1][0]).add(TINic[i][2][0]).delete(e[2][0]);
    });

    TINinv = TINinv.map(function(e) {
        return [...e];
    });

    TINinv.forEach(function(e, i) {
        e.sort(function(a, b) {
            return calAngle(numeric.sub(distorted[i], distorted[a]), [1, 0]) - calAngle(numeric.sub(distorted[i], distorted[b]), [1, 0])
        });
    });

    //console.log(TINinv);

    g.selectAll("path .line").data(TIN.features).enter().append("path")
        .attr({
            "class": "line",
            "d": path,
            "stroke": "#ccc",
            "stroke-width": 1,
            "fill": "none",
            "opacity":.5
        });
}

function updateLayoutByForce() {
    let force = [];
    for(let i=0; i<distorted.length; i++) {
        force[i] = [0, 0];
        //if(TINinv[i].size >= 5) {
        for(let j=0; j<TINinv[i].size; j++) {
            let delta = numeric.sub(distorted[i], distorted[[...TINinv[i]][j]]);
            let dist = numeric.norm2Squared(delta);

            force[i][0] -= (200000/dist)*delta[0];
            force[i][1] -= (200000/dist)*delta[1];

            let d = Math.sqrt(dist);
            force[i][0] += (d - 50)*delta[0];
            force[i][1] += (d - 50)*delta[1];
        }
        //}

    }

    for(let i=0; i<distorted.length; i++) {
        distorted[i][0] += force[i][0]*0.00001;
        distorted[i][1] += force[i][1]*0.00001;
    }

    redraw();
}

function checkTIN() {
    linecDistorted = [];
    TINi.forEach(function(t) {
        let u = numeric.sub(
            distorted[t[0][0]],
            distorted[t[0][1]]
        );
        let v = numeric.sub(
            distorted[t[1][1]],
            distorted[t[1][0]]
        );
        linecDistorted.push(crossVec(u.concat(0), v.concat(0)));
    });

    linecDistorted.forEach(function(e, i) {
        if(e[2]>0) {
            flatten(TINi[i]).filter(function(te, ti, tf) {
                return ti == tf.indexOf(te);
            }).forEach(function(c) {
                distorted[c] = predistorted[c];
            });
        }
    });

    redraw();
}

function findNodeIdx(arr) {
    for(let i=0; i<nodes.features.length; i++) {
        if(JSON.stringify(arr) == JSON.stringify(nodes.features[i].geometry.coordinates)) {
            for(let j=0; j<lines.length; j++) {
                if(lines[j].indexOf(i) != -1) {
                    return [i, j];
                }
            }
        }
    }
}

function flatten(arr) {
    return arr.reduce(function(flat, toFlatten) {
        return flat.concat(Array.isArray(toFlatten)?flatten(toFlatten):toFlatten);
    }, []);
}

function crossVec(a, b) {
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]];
}

function exportGeoJson(geojson) {
    return `{"type": "FeatureCollection","features": ${JSON.stringify(geojson)}}`;
}

function exportCoord() {
    let objSet = {
        "coords": distorted,
        "net": net,
        "edgeSet": edgeSet
    };

    let  a = window.document.createElement("a");
    a.href = window.URL.createObjectURL(new Blob([JSON.stringify(objSet)], {
        type: "text/json"
    }));
    a.download = "export.json";
    a.click();
}

function reset2optimal() {
    distorted = [[100,119.50901636574338],[129.16558528971547,119.51526232271698],[158.16413505216366,119.52771379448595],[187.00613527149665,119.54632900416178],[215.69364295929824,119.57103400348873],[288.96215084363246,100],[269.39372794150455,119.6385259643769],[289.5767774798936,139.98809779148746],[311.9593197490579,162.54964424039247],[334.52916299437237,185.3085314355897],[357.0812149996669,208.05741515667887],[379.5978409556842,230.7770340472578],[401.9068820336813,324.4907261986656],[401.98027601201284,289.0672046136959],[402.0237782931299,253.45447506722266],[401.83215644265846,352.9296946120418],[401.75526727857823,381.42691641728743],[401.67769062916943,416.57490326191686],[420.6940238743069,464.1953533005871],[670.3823226770598,594.0647524129689],[640.9584066552721,594.0600566444533],[611.5164498924232,594.0229675862837],[589.269093963861,571.7150123969373],[566.9969414936191,549.3381937044203],[544.7307164063332,526.9612504702201],[522.4607841910134,504.57579370177797],[500.18873239448027,482.1846710959028],[478.0218414914344,459.89678552557046],[460.3474344223747,414.7109605994157],[656.4209973311271,420.4755867326353],[634.4661516884166,442.3865579974234],[604.6943680896941,442.3708139864887],[574.7215631852896,442.3274457686754],[544.657208006971,442.28381731596687],[514.532167153297,442.2400198647714],[484.48567555389195,442.19610956889335],[460.38440559314233,442.15225042726377],[463.2521222829371,515.425426009642],[530.1529677738937,582.6942297469424],[419.37303057042845,535.1910098335246],[328.8226253491644,461.3429181160828],[358.6470428190058,461.2929655632679],[284.43857836918284,505.8972629844971],[262.2428017788326,528.1940259880712],[240.0488708533672,550.4854739153412],[195.6387394978405,595.0757099314642],[173.4100256432573,617.3848140186718],[552.2273677602312,604.9257777687297],[419.3426255493006,602.8181329286906],[306.6443318428493,483.5882478119099],[424.32354256326005,444.8595989846141],[403.259345064383,479.2360536051747],[420.70618366277336,496.80156046954824],[439.0888818967462,515.3548936138213],[443.6697865982477,425.3335420616409],[401.60039890051644,444.97879046888113],[385.43155584648974,461.2707089313244],[678.2927188027732,398.61256124228265],[505.05725860081753,369.54776257538873],[549.640424599534,324.55833750678363],[571.7492464125881,302.2913156238396],[600.5157875207997,302.2021676715216],[629.1700615614182,302.1444946582186],[657.5285666102135,302.1174430309885],[678.8598016822735,323.4783039913992],[700,344.72396903549736],[699.958482782475,376.9185499640766],[217.84963639181905,572.7776210477657],[507.85212156992606,560.2664682696319],[485.50015849867907,537.7812077039245],[244.15591866192574,119.60173113246712],[552.2294153405882,635.9351112349303],[552.2328633043129,666.2026732792384],[552.2292277401997,700],[419.34145754301346,636.6105293270055],[419.344982035625,568.951975911783],[443.6439101188199,395.1022742315248],[420.6015684493232,371.8527965436127],[349.84825448434344,352.9097561757566],[320.68521732456554,352.90246938317944],[291.479072353265,352.896995232037],[262.2321895993413,352.8933414739074],[232.93695709252944,352.8915110146727],[526.4813978358372,635.9410468032669],[482.6477685429817,392.16626639199956],[378.8416719655186,352.918845270017],[527.3927203259489,347.0067061814476]];
    //distorted = [[100,100],[129.81410802936273,108.86192269974407],[159.66112344889126,117.61217757228457],[189.5659981688959,126.1618709445709],[219.5402851748487,134.45864795775543],[291.1665510084174,119.10259107230371],[277.27854788183055,142.93348174553603],[292.4965186548079,166.3747854412718],[309.40554573996246,192.0859031589626],[326.2740093612231,218.16327477562857],[342.6544543100831,244.5819837695658],[358.305554286772,271.44075264575275],[400.47284689458934,354.5482944566937],[387.12764876260013,326.51571268027226],[373.1195557765392,298.76924048070777],[399.7246947737564,379.41561891999623],[410.6789278929716,401.78112063042886],[408.9405757806004,432.7561510089303],[412.1290588711049,472.0058137772797],[684.8530513619502,596.9810535141445],[656.8012748439135,583.5512963214802],[628.9015284051234,569.8106982763784],[601.2963791449091,555.490852829214],[574.1071389785023,540.3995214347472],[547.4199987854279,524.4393297881742],[521.2737240894859,507.60929225844427],[495.65596365055154,489.9927685625941],[470.53309496250785,471.75215882283663],[460.0695429060004,423.1812077736435],[666.1704130224782,451.71734700288704],[638.8407826421402,465.861107176872],[608.7294611788739,472.9543024198606],[577.706404628506,474.322834862975],[546.703095759282,471.5109634095031],[516.0236644200619,465.9795769584225],[485.6374784096948,459.0384624595338],[463.1986501218654,447.9894929736425],[452.81378261939125,520.7459512134943],[508.39450577666855,595.5215301933841],[410.4806827021551,543.176326067042],[323.9524584155167,473.6539256277326],[351.4523629810636,459.2468280384963],[269.9477686196554,504.4850805618483],[243.4629736571428,520.7831640420239],[217.28547938885325,537.5706522151384],[165.61963127755064,572.1960904439527],[139.9721879547228,589.7835070628818],[525.4590119601086,621.429158552978],[398.8353737701949,604.2182807820268],[296.7733154069249,488.7550936323583],[420.0095973876575,450.47864578149296],[383.9994303760716,494.83898539919653],[410.5922345619732,499.07762216272783],[425.15697012370583,519.9098322447621],[443.69709015468993,436.3678108440752],[395.32263637406555,453.6040020080164],[377.58152717716763,468.10704634099034],[687.761691003264,430.07016677451566],[499.78950677679205,375.2406273706324],[547.3604840214988,335.38594347040555],[575.0482216471629,321.49705931631155],[604.9477039419404,313.87450528489484],[635.6426924397134,314.42073114523157],[664.403413453123,324.626652164209],[687.2539180207831,344.5914747518141],[699.9771253089119,372.0199801883343],[700,402.2667801934045],[191.36179029530746,554.7477860735122],[490.56177760099376,570.0652387952249],[472.01571057738164,545.1223381992739],[249.5711658679782,142.48341560820836],[515.7181959003535,647.3077525678297],[530.2495863760298,670.8212140951456],[519.7285161607025,700],[393.2839855615473,634.8034555215247],[404.5248963093276,573.6608138346925],[445.6752101662501,408.9294421687726],[423.22110465523684,387.6117756352771],[344.8960599575605,368.9828681788447],[314.8570562691593,360.92637668720806],[284.8446377230158,352.7613990615803],[254.85286872947108,344.519909521123],[224.87201867489983,336.23869351547154],[488.3421254480944,650.2352487780678],[479.2510511333903,398.66932720084765],[374.945070193464,376.9001536923836],[522.3277240367654,353.777972800306]];
    redraw();
}

function rollback() {
    distorted = predistorted;
    redraw();
}

function getCplexAns(arr) {
    for(let i=0; i<distorted.length; i++) {
        distorted[i] = [arr[i], arr[i+distorted.length]];
    }
    redraw();
}

function labelreset() {
    let svg = d3.select("svg");
    let mv = 1;
    while(mv > 0) {
        mv = 0;
        svg.selectAll("text").each(function() {
            var that = this;
            var label = this.getBoundingClientRect();

            svg.selectAll("text").each(function() {
                if(this!=that) {
                    var anotherLabel = this.getBoundingClientRect();
                    if((Math.abs(label.left - anotherLabel.left)*2 < (label.width + anotherLabel.width)) &&
                        Math.abs(label.top - anotherLabel.top)*2 < (label.height + anotherLabel.height)) {
                        var dx = (Math.max(0, label.right - anotherLabel.left) + Math.min(0, label.left - anotherLabel.right)) * 0.01,
                            dy = (Math.max(0, label.bottom - anotherLabel.top) + Math.min(0, label.top - anotherLabel.bottom)) * 0.02,
                            tt = d3.transform(d3.select(this).attr("transform")),
                            to = d3.transform(d3.select(that).attr("transform"));
                        mv += Math.abs(dx) + Math.abs(dy);

                        to.translate = [ to.translate[0] + dx, to.translate[1] + dy ];
                        tt.translate = [ tt.translate[0] - dx, tt.translate[1] - dy ];
                        d3.select(this).attr("transform", `translate(${tt.translate})`);
                        d3.select(this).attr("fill", "red");
                        d3.select(that).attr("transform", `translate(${to.translate})`);
                        d3.select(that).attr("fill", "red");
                        label = this.getBoundingClientRect();
                    }
                }
            });
        });
    }
}

window.exportGeoJson = exportGeoJson;
window.reset2optimal = reset2optimal;
window.rollback = rollback;
window.getCplexAns = getCplexAns;
window.labelreset = labelreset;
window.runOctProcess = runOctProcess;
window.runOptimialProcess = runOptimialProcess;