---
title: "Node Labels"
menu:
  main:
    identifier: "25-08-22-node-labels"
    parent: "2025"
    weight: 20
---

_By Sören Domrös, August 22, 2025_

I always seem to forget how to properly add node labels in ELK such that it adjusts its size based on them. Therefore, this little guide should solve your problems with node labels and help me remember all the required options and configurations.

The full example model I will show in the following, can be found [here](https://rtsys.informatik.uni-kiel.de/elklive/examples.html?e=labels%2FnodeLabelAndParentSize).

My error in properly letting a parent node consider the size of a node label is that just setting [`nodeSize.constraints`](https://eclipse.dev/elk/reference/options/org-eclipse-elk-nodeSize-constraints.html) to an EnumSet containing `NODE_LABELS` is not enough.
Additionally, one needs to specify the [`nodeLabels.placement`](https://eclipse.dev/elk/reference/options/org-eclipse-elk-nodeLabels-placement.html).

### Node Labels and Node Size for Child Nodes

This model results in a node that does not consider the size of the node labels, which is not what we want.
<div style="display: flex; flex-direction: row; gap: 2em; align-items: flex-start;">
  <div style="flex: 2; min-width: 0;">
    <pre><code>node childWithLabelCenter {
    layout [size: 10, 10]
    nodeSize.constraints: "NODE_LABELS MINIMUM_SIZE"
    label "childWithLabelCenter"
    label "childWithLabelCenter2"
}
</code></pre>
  </div>
  <div style="flex: 1; min-width: 0; margin:auto;">
    {{< image src="nodeLabelPlacement-bad.svg" alt="The node's size without specified node label placement.">}}
  </div>
</div>

The following model, however, does work as expected and results in a node that considers the size of the node labels.

<div style="display: flex; flex-direction: row; gap: 2em; align-items: flex-start;">
  <div style="flex: 2; min-width: 0;">
    <pre><code>node childWithLabelCenter {
    layout [size: 10, 10]
    nodeLabels.placement: "INSIDE H_CENTER V_CENTER"
    nodeSize.constraints: "NODE_LABELS MINIMUM_SIZE"
    label "childWithLabelCenter"
    label "childWithLabelCenter2"
}
</code></pre>
  </div>
  <div style="flex: 1; min-width: 0; margin:auto;">
    {{< image src="nodeLabelPlacement-correct.svg" alt="The node's size with specified node label placement.">}}
  </div>
</div>

If labels are not centered, the node will additionally reserve space on the opposite side of the label to be symmetrical.
This can be disabled by setting `nodeLabels.options` to `ASYMMETRICAL`.

The following example shows the symmetrical space reservation for node labels.
<div style="display: flex; flex-direction: row; gap: 2em; align-items: flex-start;">
  <div style="flex: 2; min-width: 0;">
    <pre><code>node childWithLabelTopLeft {
    nodeLabels.placement: "INSIDE H_LEFT V_TOP"
    nodeSize.constraints: "NODE_LABELS"
    label "childWithLabelTopLeft"
    label "childWithLabelTopLeft2"
}
</code></pre>
  </div>
  <div style="flex: 1; min-width: 0; margin:auto;">
    {{< image src="nodeLabelPlacement-symmetrical.svg" alt="The node's size with symmetrical space reservation for node labels.">}}
  </div>
</div>

While this example shows how to disable it.

<div style="display: flex; flex-direction: row; gap: 2em; align-items: flex-start;">
  <div style="flex: 2; min-width: 0;">
    <pre><code>node childWithLabelTopLeftAsymmetric {
    nodeLabels.placement: "INSIDE H_LEFT V_TOP"
    nodeSize.options: ASYMMETRICAL
    nodeSize.constraints: "NODE_LABELS"
    label "childWithLabelTopLeftAsymmetric"
    label "childWithLabelTopLeftAsymmetric2"
}
</code></pre>
  </div>
  <div style="flex: 1; min-width: 0; margin:auto;">
    {{< image src="nodeLabelPlacement-asymmetrical.svg" alt="The node's size with asymmetrical space reservation for node labels.">}}
  </div>
</div>

### Node Labels and Node Size for Parent Nodes

The key difference between node labels for child and parent nodes is that the node label placement for parent nodes depends on the `elk.algorithm` set on the parent node.
In this case, only the `layered` algorithm supports all possible node label placements and will reserve space for the node labels to not intersect with potential child nodes.

This model adds a top and a left padding for the child nodes.
Note, however, that the `ASYMMETRICAL` option is not supported for parent nodes since there is no symmetrical space reservation for parent nodes.

<div style="display: flex; flex-direction: row; gap: 2em; align-items: flex-start;">
  <div style="flex: 2; min-width: 0;">
    <pre><code>node parentWithLabelSpaceTopLeft {
    nodeLabels.placement: "INSIDE H_LEFT V_TOP"
    nodeSize.constraints: "NODE_LABELS MINIMUM_SIZE"
    nodeSize.options: ASYMMETRICAL
    label "parentWithLabelSpaceTopLeft"
    label "parentWithLabelSpaceTopLeft2"
    node n1
    node n2
}
</code></pre>
  </div>
  <div style="flex: 1; min-width: 0; margin:auto;">
    {{< image src="nodeLabelPlacement-parent.svg" alt="The parent node's size with node label placement.">}}
  </div>
</div>

Also note that aligning the label top-left will reserve seemingly unnecessary space on the left.
Here, it would be enough to only reserve space on the top, which is a downside of the algorithm.
Typically, this is solved by estimating the size of the label and reserving space by setting a `padding` on the node.
This, however, has to be solved outside of ELK.
