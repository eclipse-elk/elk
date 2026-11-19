---
title: "ELK Text Format"
menu:
  main:
    identifier: "ElkTextFormat"
    parent: "GraphDataStructure"
    weight: 50
---

The ELK text format has five basic elements: nodes, ports, labels, edges and edge sections.
Details about each element can be found below.
Note that the ELK text format uses default node sizes if they are not specified to make writing elkt files more convenient.
This is different from the JSON format, where no such defaults are assumed.

## Nodes, Ports, Labels, Edges, and Edge Sections

Nodes and ports must have an _id_ that uniquely identifies them.
Labels are usually not referred to from other parts of the graph, which is why the id is optional.
The id is a string that has to begin with a letter.
All elements furthermore can have [_layout options_]({{< relref "reference/options.md" >}}).
Layout options are basically a list of key-value pairs that are used to assign layout option values to the element.
Layout options can be an enum value, an enum set, an integer, a double, or other more complex objects.
For each layout option the type is specified in the [_reference_]({{< relref "reference/options.html" >}}).

```elkt
layoutOption1: ENUM_VALUE
layoutOption2: "ENUM_VALUE_1 ENUM_VALUE_2"
layoutOption3: 10.0
```

A special case is [_padding_]({{< relref "reference/options/org-eclipse-elk-padding.html" >}}).

```elkt
padding: "[top=20.0,left=20.0,bottom=20.0,right=20.0]"
```

[_Individual spacings_]({{< relref "reference/options/org-eclipse-elk-spacing-individual.html" >}}) are essentially a map of options that use a special syntax.

```elkt
spacing.individual: "spacing.nodeNode: 10.0 :: spacing portsSurrounding: [top=40.0,left=10.0,bottom=10.0,right=30.0]"
```

## Nodes, Ports, and Labels

Nodes, ports, and labels have a two-dimensional location and size.
Each of these elements can also have an arbitrary number of labels to describe them properly.
Yes, even labels can have labels, although it depends on the layout algorithm whether or not it supports labeled labels.
The optional `layout` block may be used to specify the size and positions of these elements.

```elkt
layout [
  position: x, y
  size: w, h
] 
```

## Nodes

Nodes can have an arbitrary number of ports.
Edges can connect to a node either directly or through one of its ports.
A node can also contain an arbitrary number of child nodes.
A graph is actually nothing more than a simple node whose children are the top-level nodes of the graph.
Finally, a node can contain edges.
While it is common to define those edges under a given node that connect that node's children, in fact any edge may be defined under any node, regardless of its end points.
This allows for flexibility when defining hierarchy-crossing edges, as well as for alternative schemes, such as defining all edges at the root level.
Note that this is not recommended and in some cases poorly supported.
See [_Coordinate System_]({{< relref "documentation/tooldevelopers/graphdatastructure/coordinatesystem.md" >}}) for the rules for interpreting edge coordinates.

```elkt
node n {
  layout [
    position: x, y
    size: w, h
  ]
  // list of layout options
  
  // definitions of any contained labels, ports, edges and other nodes
}
```

## Ports

Ports do not have any more interesting properties. Ports are boring.

## Labels

Labels contain text.
Note that layout algorithms generally don't perform any size estimation of the text.
Therefore you should specify a reasonable width and height.

```elkt
label "A magnificent text"
```

## Edges

There are two types of edges: primitive edges and extended edges.
Primitive edges are solely supported for legacy models to work.
Exported graphs will always be made up of extended edges.
Both kind of edges support labels.

### Primitive Edges

Primitive edges have a source and target node and can optionally connect to a source port and target port.

```elkt
edge source_node_id.source_port_id -> target_node_id.target_port_id {
  layout [
    start: x, y
    bends: x, y | x, y | x, y
    end: x, y
  ],
  // any attached labels
}
```

### Extended Edges

Extended edges have two mandatory arrays consisting of the identifiers of nodes and ports.
One array defines the edge's source elements, the other defines its target elements.
Edges may well connect more than one source to more than one target, making them hyperedges.
It should be noted that currently none of the algorithms implemented in ELK explicitly support hyperedges.
Note that many layout algorithms don't support hyperedges.
If an edge has a layout, it can specify an arbitrary number of edge sections
that define said layout.
A simple edge with one source and one target only needs a single section.

```elkt
edge source_node_id.source_port_id, another_source -> target_node_id.target_port_id, another_target {
  layout [
    // list of edge sections
  ],
  // any attached labels
}
```

## Edge Sections

Edge sections are only used in conjunction with extended edges and capture the routing of an edge through a drawing.
Each section connects two end points.
An end point can be one of the end points of the section's edge (a node or a port), or one or more other edge sections.
The points where edge sections meet are _junction points_ where one part of the edge branches off.
An edge section can only have either an incoming shape or incoming edge sections (the same is true of course for outgoing shapes and outgoing edge sections).
In the simplest case, an edge only has a single edge section which runs from the edge's single source to its single target.
In this case, it is enough to define the section's start and end point and possibly bendpoints.
Incoming and outgoing shapes are then filled in automatically by the importer.
Specifying edge sections in the elkt format is unnecessary in all normal use cases, because the layout algorithms will create them where necessary.

```elkt
section incoming_sections -> outgoing_sections [
  start: x, y
  end: x, y
  bends: x, y | x, y | x, y
  incoming: node and / or port identifier
  outgoing: node and / or port identifier
]
```

## Examples
Many examples can be found in [_elklive_](https://rtsys.informatik.uni-kiel.de/elklive/examples.html).
