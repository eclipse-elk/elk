Vertiflex is a special tree drawing algorithm that diverges from the classic aesthetic constraint that requires nodes to be assigned to discrete layers.
Instead each node may now come with its own defined height constraint.
Nodes that do not define their own height constraint are assigned a location based on the defined spacings and the other nodes in the graph.

![Vertiflex example](/images/vertiflex.png)

The wrapper node in the example illustrates the vertical constraints.
The node *n1* has a vertical constraint of 115 (see [_ELK's coordinate system_]({{< relref "documentation/tooldevelopers/graphdatastructure/coordinatesystem.html" >}}).
The node's *n4*, *n7* and *n8* have no vertical constraint and are placed below *n3* to satisfy the specified node node spacing.

