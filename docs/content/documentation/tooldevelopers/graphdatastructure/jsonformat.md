---
title: "JSON Format"
menu:
  main:
    identifier: "JsonFormat"
    parent: "GraphDataStructure"
    weight: 40
---


The JSON graph format has five basic elements: nodes, ports,
labels, edges, and edge sections. Details about each element can be found below, with some
sections describing features common to multiple elements. Note that in the JSON code, _mandatory_
fields are marked with an _asterisk_.


## Nodes, Ports, Labels, Edges, and Edge Sections

All elements, except labels, must have an _id_ that uniquely identifies them.
Labels are usually not referred to from other parts of the graph,
which is why the id is optional.
The id can be a string or an integer.
All elements furthermore can have [_layout options_]({{< relref "reference/options.md" >}}).
Layout options are basically a list of key-value pairs that are used to
assign layout option values to the element.


```json
{
  id*: "ID",
  layoutOptions: { ..object with key value pairs.. }
}
```

## Nodes, Ports, and Labels

Nodes, ports, and labels have a two-dimensional location and size. Each of these elements
can also have an arbitrary number of labels to describe them properly. Yes, even labels can
have labels, although it depends on the layout algorithm whether or not it supports labeled
labels.

```json
{
  x: ...,
  y: ...,
  width: ...,
  height: ...,
  labels: [ ..array of label objects.. ]
}
```


## Nodes

Nodes can have an arbitrary number of ports. Edges can connect to a node either directly or
through one of its ports. A node can also contain an arbitrary number of child nodes. A graph
is actually nothing more than a simple node whose children are the top-level nodes of the graph.
Finally, a node can contain edges. While it is common to define those edges under a given node
that connect that node's children, in fact any edge may be defined under any node, regardless of its
end points. This allows for flexibility when defining hierarchy-crossing edges, as well as for
alternative schemes, such as defining all edges at the root level.
See [_Coordinate System_]({{< relref "documentation/tooldevelopers/graphdatastructure/coordinatesystem.md" >}})
for the rules for interpreting edge coordinates.

```json
{
  ports: [ ..array of port objects.. ],
  children: [ ..array of child node objects.. ],
  edges: [ ..array of edge objects.. ]
}
```


## Ports

Ports do not have any more interesting properties. Ports are boring.


## Labels

Labels can additionally contain text. Note that layout algorithms
generally don't perform any size estimation of the text.
Therefore you should specify a reasonable width and height.

```json
{
  text: "A magnificent text"
}
```

## Edges

There are two types of edges: primitive edges and extended edges.
Primitive edges are solely supported for legacy models to work.
Exported graphs will always be made up of extended edges. Both kind 
of edges support labels.

### Primitive Edges

Primitive edges have a source and target node and can optionally connect
to a source port and target port.
```json
{
  source*: node identifier,
  sourcePort: port identifier,
  target*: node identifier,
  targetPort: port identifier,
  sourcePoint: {x, y},
  targetPoint: {x, y},
  bendPoints: [ .. {x, y} pairs .. ],
  labels: [ ..array of label objects.. ]
}
```

### Extended Edges

Extended edges have two mandatory arrays consisting of the identifiers of nodes and ports. One
array defines the edge's source elements, the other defines its target elements. Edges may well
connect more than one source to more than one target, making them hyperedges.
Note that many layout algorithms don't support hyperedges.
If an edge has a layout, it can specify an arbitrary number of edge sections
that define said layout. A simple edge with one source and one target only needs a single section.

```json
{
  sources*: [ ..array of node and / or port identifiers.. ],
  targets*: [ ..array of node and / or port identifiers.. ],
  sections: [ ..array of edge sections.. ],
  labels: [ ..array of label objects.. ]
}
```


## Edge Sections

Edge sections are only used in conjunction with extended edges and
capture the routing of an edge through a drawing. Each section connects two
end points. An end point can be one of the end points of the section's edge (a node or a port),
or one or more other edge sections. The points where edge sections meet are _junction
points_ where one part of the edge branches off. An edge section can only have either an
incoming shape or incoming edge sections (the same is true of course for outgoing shapes and
outgoing edge sections). In the simplest case, an edge only has a single edge section which
runs from the edge's single source to its single target. In this case, it is enough to define
the section's start and end point and possibly bendpoints.
Incoming and outgoing shapes are then filled in automatically by the importer.

```json
{
  startPoint*: {x, y},
  endPoint*: {x, y},
  bendPoints: [ ..array of {x, y} pairs.. ],
  incomingShape: node and / or port identifier,
  outgoingShape: node and / or port identifier,
  incomingSections: [ ..array of edge section identifiers.. ],
  outgoingSections: [ ..array of edge section identifiers.. ]
}
```

## Junction Points

_Junction points_ are the split and merge points of hyperedges. 
Layout algorithms supporting hyperedges may compute these points 
such that a rendering framework can use them to position visual cues, 
for instance, small circles. 
In case an algorithm computes junction points, 
an edge's representation might look like the following after layout.

```json
{
  id: "edge0",
  junctionPoints: [ ..array of {x, y} pairs.. ]
}
```


# Small Example
```json
{
  id: "root",
  properties: { "elk.direction": "RIGHT" },
  children: [
    { id: "n1", width: 10, height: 10 },
    { id: "n2", width: 10, height: 10 }
  ],
  edges: [{
    id: "e1", sources: [ "n1" ], targets: [ "n2" ]
  }]
}
```


# Java API

Importing and exporting can be done using the `ElkGraphJson` utility class.
It provides a set of methods that are outlined next.

## JSON to ELK Graph

Simple import using a json string as input:
```java
ElkNode root = ElkGraphJson.forGraph(jsonString)
                           .toElk();
```

Sometimes one desires to transfer the computed layout back to
the initial json graph.
This is possible if the graph is already available
as `JsonObject`. Pass a
`Maybe` instance to the builder which is populated with
the `JsonImporter` instance that has been used to
import the json graph. Internally, the importer maintains
several maps from json elements to ELK graph elements.

```java
JsonObject jsonGraph = ...;
Maybe<JsonImporter> importerMaybe = new Maybe<>();
ElkNode root = ElkGraphJson.forGraph(jsonGraph)
                           .rememberImporter(importerMaybe)
                           .toElk();
// [ perform layout ]

importerMaybe.get().transferLayout(root);
```

## ELK Graph to JSON

```java
ElkNode root = ...;
String json = ElkGraphJson.forGraph(root)
                          .omitLayout(true)
                          .omitZeroDimension(true)
                          .omitZeroPositions(true)
                          .shortLayoutOptionKeys(false)
                          .prettyPrint(true)
                          .toJson();
```