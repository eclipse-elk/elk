---
title: "Layered (overview)"
menu:
  main:
    identifier: "25-08-21-layered"
    parent: "2025"
    weight: 10
---

_By Sören Domrös, August 21, 2025_

The layered algorithm or Sugiyama algorithm is in its implementation in ELK divided into five phases with intermediate processing slots: Cycle Breaking, Layer Assignment, Crossing Minimization, Node Placement, and Edge Routing.

{{< image src="layered-algorithm-structure.svg" alt="ELK Layered phases." >}}

Using the occasion of the upcoming ELK 0.11.0 release, I want to shed more light into the ELK Layered algorithm and how to find the correct layout phase to configure.

In the following, I will assume that the layout direction is left-to-right, or speaking in ELK layout options `elk.direction: RIGHT` since this makes talking about coordinates easier and allows me to use horizontal and vertical as well as left and right without discussion potential alternatives.

## Cycle Breaking

The [cycle breaking](https://eclipse.dev/elk/reference/options/org-eclipse-elk-layered-cycleBreaking-strategy.html) phase (`cycleBreaking.strategy`) of the ELK layered algorithm make a graph acyclic by reversing edges such that they go from right-to-left instead of left-to-right creating *backward edges*. However, in its implementation in ELK, it does a lot more than that since this phase typically also handles constraints that somehow determine the horizontal ordering of connected nodes.

This phase is implemented by the following strategies:

- `GREEDY` (default): Aims to minimize the number of backward edges and uses random decisions on tie
- `INTERACTIVE`: Use the ordering given by the positions of nodes (e.g. from a previous layout run) to determine the direction of edges
- `MODEL_ORDER`: Use the ordering of nodes in the input model to determine the direction of edges
- `DEPTH_FIRST`: Traverse the graph depth first beginning with the first  source in the list of nodes and using the edge order as the depth-first visiting order. All edges that lead back to already visited nodes will be reversed
- `GREEDY_MODEL_ORDER`: The `GREEDY` approach but using the model order as a tie-breaker
- `SCC_CONNECTIVITY`: Determines the strongly connected components using Tarjan's algorithm while using the in-out degree to determine edges to reverse, which results in quadratic runtime compared to the other approaches
- `SCC_NODE_TYPE`: Same as `SCC_CONNECTIVITY` but reverses edges going out of or going to specifically marked nodes using group model order
- `DFS_NODE_ORDER`: Same as `DEPTH_FIRST` but uses the node model order instead of the edge order to traverse the graph depth-first.
- `BFS_NODE_ORDER`: Same as `DFS_NODE_ORDER` but breadth-first.

Moreover, the following options might constrain the cycle breaking step:

- Partitioning: [Activating partitioning](https://eclipse.dev/elk/reference/options/org-eclipse-elk-partitioning-activate.html) and [partitions](https://eclipse.dev/elk/reference/options/org-eclipse-elk-partitioning-partition.html) not only determine partitions in the graph but also determine the direction of edges
- Layer constraints: Setting [layer constraints](https://eclipse.dev/elk/reference/options/org-eclipse-elk-layered-layering-layerConstraint.html) to `FIRST` or `LAST` also influences the direction of edges to make this possible
- Direction priority: Setting [`priority.direction`](https://eclipse.dev/elk/reference/options/org-eclipse-elk-layered-priority-direction.html) makes it possible to weight the importance of edges being not reversed

## Layer Assignment

The [layer assignment](https://eclipse.dev/elk/reference/options/org-eclipse-elk-layered-layering-strategy.html) phase (`layering.strategy`) assigns each node (and also each edge label, which is essentially only a special kind of node) a vertical layer. I.e., using the acyclic graph created by cycle breaking, it aims to create a minimal width assignment of nodes and edge labels to layers, which may create desired vertical alignment.

This phase supports the following strategies:

- `NETWORK_SIMPLEX` (default): Encodes the layering problem as a network simplex problem as proposed by [Gansner et al.](https://doi.org/10.1109/32.221135)
- `LONGEST_PATH`: Start from the sinks and greedily assign nodes to layers
- `LONGEST_PATH_SOURCE`: Start from the sources and greedily assign nodes to layers
- `COFFMAN_GRAHAM`: Uses a precedence-constraints scheduling approach by [Coffman and Graham](https://doi.org/10.1007/BF00288685) to find a layer assignment bounded by the [number of nodes in a layer](https://eclipse.dev/elk/reference/options/org-eclipse-elk-layered-layering-coffmanGraham-layerBound.html)
- `INTERACTIVE`: Uses previously set positions of nodes to determine non-overlapping vertical layers
- `STRETCH_WIDTH`  (Experimental): Similar to `LONGEST_PATH` but tries to reduce the maximum number of nodes in a layer. This is based on the work of [Nikolov et al.](https://doi.org/10.1145/1064546.1180618)
- `MIN_WIDTH` (Experimental): A heuristic to solve the NP-hard minimum-width layering problem while considering dummy nodes create for routing edges. This is based on the work of [Tarassov et al.](https://doi.org/10.1007/978-3-540-24838-5_42)
- `BF_MODEL_ORDER` (Experimental): Assign nodes to layers breadth-first based on their model order. This assumes that the cycle breaking was also done by model order.
- `DF_MODEL_ORDER` (Experimental): Same as `BF_MODEL_ORDER` but depth-first. I must note here that this is not a very efficient solution and that `LONGEST_PATH_SOURCE` is typically a better alternative. 

Moreover, the following options can be used to constrain, configure, or post-process the layering:

- Partitioning: See partitioning explained above
- Layer constraints: See layer constraints explained above
- Node promotion: Uses [Nikolov et al.'s](http://doi.acm.org/10.1145/1064546.1180618) node promotion strategy to achieve a layering with fewer dummy nodes
- Layer splitting: Split a very high layer into multiple uses using [layer unzipping](https://eclipse.dev/elk/reference/options/org-eclipse-elk-layered-layerUnzipping-strategy.html).
- Edge labels: TODO I think there are options to control their placement, right?

## Crossing Minimization

The [crossing minimization](https://eclipse.dev/elk/reference/options/org-eclipse-elk-layered-crossingMinimization-strategy.html) phase (`crossingMinimization.strategy`) determines the vertical ordering of nodes and the relative routes by also ordering the dummy nodes associated with them. Therefore, this phase makes sure that the drawing has fewer edge crossings and conforms to relative vertical constraints and also determines the ordering of ports/edges.

It supports the following strategies:

- `LAYER_SWEEP` (default): This strategy uses the barycenter heuristic proposed by [Sugiyama et al.](https://doi.org/10.1109/TSMC.1981.4308636) to find a crossing minimal ordering for nodes and the ports edges connect to
- `MEDIAN_LAYER_SWEEP`: This strategy uses the median heuristic proposed by Di Battista et al. instead of the barycenter heuristic
- `INTERACTIVE`: This strategy uses previously calculated node positions to determine the ordering of nodes in a layer
- `NONE`: This does no crossing minimization but only executes desired pre- and post-processing steps. It is useful to constrain the order of nodes and edges

This strategy can be constrained as explained in detail in the [Layered: Constraining the Model blog post](https://eclipse.dev/elk/blog/posts/2023/23-01-09-constraining-the-model.html).

Additionally, one can greedily post-process the ordering to reduce edge crossings using the [greedy switch heuristic](https://eclipse.dev/elk/reference/options/org-eclipse-elk-layered-crossingMinimization-greedySwitch-type.html).

## Node Placement

The [node placement](https://eclipse.dev/elk/reference/options/org-eclipse-elk-layered-nodePlacement-strategy.html) phase (`nodePlacement.strategy`) determines the concrete `y` coordinates of nodes since vertical routing segments created during edge routing might change `x` coordinates later on. Moreover, this phase determine horizontal alignment and centering of nodes on its connections.

Node placement supports the following strategies:

- `SIMPLE`: Center all nodes and edges in their respective layer, separating them by their spacing. This strategy is not intended for a serious use case
- `INTERACTIVE`: Leaves the previously calculated `y` coordinates of the nodes untouched
- `LINEAR_SEGMENTS`: Align nodes edges using linear segments and nodes using the pendulum method proposed by [Sander](https://doi.org/10.1007/BFb0021828). This can be additionally configured by setting [the dampening](https://eclipse.dev/elk/reference/options/org-eclipse-elk-layered-nodePlacement-linearSegments-deflectionDampening.html).
- `BRANDES_KOEPF` (default): Use [Brandes & Köpf's](https://doi.org/10.1007/3-540-45848-4_3) method for node placement. Per default it orients nodes to the top-left but this can be [configured to balance nodes based on their connections](https://eclipse.dev/elk/reference/options/org-eclipse-elk-layered-nodePlacement-bk-fixedAlignment.html).
- `NETWORK_SIMPLEX`: Uses a network simplex approach by [Gansner et al.](https://doi.org/10.1109/32.221135) to determine node positions.

Notable options to configure this are the following:

- [Post compaction](https://eclipse.dev/elk/reference/options/org-eclipse-elk-layered-compaction-postCompaction-strategy.html): Squishes the whole layout together to save space
- [Network simplex node flexibility](https://eclipse.dev/elk/reference/options/org-eclipse-elk-layered-nodePlacement-networkSimplex-nodeFlexibility.html): Dynamically increase the node height to straighten edges

## Edge Routing

ELK layered's [edge routing](https://eclipse.dev/elk/reference/options/org-eclipse-elk-edgeRouting.html) phases supports different drawing styles and finally determines the `x` coordinates of nodes and the routes for edges. This can be used to create edge control points for polyline, orthogonal, or spline edge routing.

I hope this helps readers by giving an overview of the different strategies without having to dig through different layers of documentation, or at least helps by highlighting the layout phase that is responsible for a problem.

