---
title: "Layer Unzipping"
menu:
  main:
    identifier: "24-12-06-layer-unzipping"
    parent: "2024"
    weight: 10
---

_By Maximilian Kasperowski, December 6, 2024_

The coming update (ELK 0.9.2) introduces a new feature to control the positioning of nodes in a layer for further compaction.
Layer Unzipping can split up the nodes of a single layer into multiple layers.
This reduces the total height and can be helpful in certain situations.

The update introduces three new properties:
- `org.eclipse.elk.layerUnzipping.strategy`
- `org.eclipse.elk.layerUnzipping.layerSplit`
- `org.eclipse.elk.layerUnizpping.resetOnLongEdges`

This [elklive example](https://rtsys.informatik.uni-kiel.de/elklive/elkgraph.html?compressedContent=MYewdgzglgJgpgJwLIngGwPIPggdBAFwQEMC4BzATwC4ACABQCUBRAMWcYH1mARAcWYBlAFBpilRAFUwALygAHeVDDl8RUhRq0AggBkAKhwBy2-QEkjfYcLCo4tMAEZa1uDHL2ntALQA+BwBMLm4eDs5+DgDMwiGeQRFgACwx7p4ArD7+Sa6pYZkOybF5CWkpoV4JAGzWAPQ1tABmUAiEtGIIoWISCDZ2gbQA3sK0I23iUrIKSiq4CHAQcAQYYLrg5MypEHQNxGgLtHW0MFAQxABGaPbwOwCuaAS0Z3AAFsQAblAgNwi0xGAwtCIcFIbTWtFiEF+kLm5DuxB+tngEGEAF9evAooNhqMuhM5IplKpcQhBPI0FACHQAqj0Z5Ei5EekGX0wJUXMJDgtQP8xh17MTaYFnENRmNutJ8dMieMSWSKXRomjGVFnIKkqrlWA0hqWZVVWU6fl1cFcqyjXqcuUAuEsgB2E3lSI2hz2g0FZ1gV1FLUer2mvVG11q21AA) demonstrates the properties and the exact behaviour of the properties is explained in the sections below.

## Layer Unzipping Strategy
For now there is only the strategy `LayerUnzippingStrategy.ALTERNATING`.
It evenly distributes the nodes into several sub-layers.
The first node goes into the first sub-layer, the second goes into the second and so on.
The default configuration is to use two sub-layers.

## Configuring the Number of Sub-layers
The number of sub-layers is two by default and can be changed with the `layerSplit` property.
The property applies to an entire layer.
To use it, it must be set on any node of that layer.
If multiple values are set, then the lowest value is used.

## Long Edge Treatment
Under the hood long edges are implemented using invisible dummy nodes.
The default behaviour is begin the alternatinon anew after a long edge, but this behaviour can be disabled for a layer by setting `resetOnLongEdges` to `false` for any node in the layer.