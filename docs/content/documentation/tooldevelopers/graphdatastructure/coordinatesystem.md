---
title: "Coordinate System"
menu:
  main:
    parent: "GraphDataStructure"
    weight: 10
---

When talking about the layout of a graph, we have to agree on how coordinates of graph elements are to be stored and interpreted. For the ELK graph, this is how it works:

{{< image src="layout-reference-points.png" alt="Coordinate system reference." >}}

The coordinates of most elements are relative to their parent element. There are a few exceptions:

* Edge labels are relative to the coordinate system their edge is relative to.
* Source points, bend points, and target points of edges are relative to the edge's containing node. Following the usual containment conventions (as implemented in `ElkGraphUtil`) results in what is displayed in the image above.

## JSON

One way of using ELK is to define a graph using JSON, and then
transfer the computed layout back to the initial JSON graph.
See [_JSON Format_]({{< relref "documentation/tooldevelopers/graphdatastructure/jsonformat.md" >}}).
When working this way, it is possible to control what type of coordinates you get,
using the `org.eclipse.elk.json.shapeCoords` and `org.eclipse.elk.json.edgeCoords` layout options.

### Shape Coords

This setting controls the coordinate system that will be used for all nodes, ports, and labels of nodes and ports.

There are three possible values:

* `INHERIT`: Inherit the value from the parent element.
* `PARENT`: Coordinates are relative to the parent element.
* `ROOT`: Coordinates are relative to the root node, i.e. global coordinates.

The default value is `INHERIT`, except at the root node, where the default is `PARENT`.
Thus, the default behavior is as described in the image at the top of this page.
Setting `ROOT` allows you to instead use a "flat" coordinate system.

### Edge Coords

This setting controls the coordinate system that will be used for all source, bend, and target points of edges,
as well as all edge labels.

There are four possible values:

* `INHERIT`: Inherit the value from the parent element.
* `CONTAINER`: Coordinates are relative to the edge's proper container node (see below).
* `PARENT`: Coordinates are relative to the node in which the edge was defined in the given JSON.
* `ROOT`: Coordinates are relative to the root node, i.e. global coordinates.

The default value is `INHERIT`, except at the root node, where the default is (for historical reasons) `CONTAINER`.
Thus, the default behavior is as described at the top of this page.

The definition of an edge's proper container is given at the end of the page on the
[_Graph Data Structure_]({{< relref "documentation/tooldevelopers/graphdatastructure.md" >}}).
It involves lowest common ancestors, and can be tricky to compute correctly.
In particular, since you are allowed to define an edge under any node whatsoever in the input JSON, the edge's
*parent* node (the one it is defined under) may be different from its container.

If you elect to work in `CONTAINER` mode, ELK helps you by writing a `container` property into each edge after layout,
giving the `id` of the edge's container node. This can then be used to interpret the coordinates.

If you elect to work in `PARENT` or `ROOT` mode instead, then you do not need to know about edge containers at all.
