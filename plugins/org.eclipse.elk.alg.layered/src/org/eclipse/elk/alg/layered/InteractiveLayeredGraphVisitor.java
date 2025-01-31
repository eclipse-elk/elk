/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0  
 *******************************************************************************/
package org.eclipse.elk.alg.layered;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.eclipse.elk.alg.layered.options.CrossingMinimizationStrategy;
import org.eclipse.elk.alg.layered.options.CycleBreakingStrategy;
import org.eclipse.elk.alg.layered.options.LayerConstraint;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.LayeringStrategy;
import org.eclipse.elk.alg.layered.options.OrderingStrategy;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.options.Direction;
import org.eclipse.elk.core.util.IGraphElementVisitor;
import org.eclipse.elk.graph.ElkConnectableShape;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkGraphElement;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.eclipse.elk.graph.properties.IProperty;

/**
 * Graph visitor which visits only the root node and recursively steps through the graph
 * to set interactive options and pseudo coordinates for the layered algorithm.
 * Interactive layout works via two consecutive layout runs.
 * The first one is a "normal" run. The second one will is configured using this graph visitor.
 * The graph visitor recursively sets the interactive strategies for all phases and assigns pseudo positions to all
 * node based on their last position and the desired layer and position in the layer expressed by the
 * {@code LayerChoiceConstraint} and {@code PositionChoiceConstraint}.
 * The pseudo position have to represent layers and the ordering in the layer.
 * <br>
 * This class is used by adding it as an additional graph visitor to configure a layout run
 * in the {@link DiagramLayoutEngine}.
 * <pre>
 * DiagramLayoutEngine.Parameters params = new DiagramLayoutEngine.Parameters();
 * params.addLayoutRun(InteractiveLayeredGraphVisitor);
 * DiagramLayoutEngine.invokeLayout(workbenchPart, diagramPart, params);
 * </pre>
 * 
 * This graph visitor is added after the normal layered layout configurator since it needs previous layers, positions,
 * sizes, and coordinates to work with.
 */
public class InteractiveLayeredGraphVisitor implements IGraphElementVisitor {

    /**
     * Since nodes should not overlap and their width and height are only the width and height of the previous
     * layout run and therefore not reliable, we have to add additional spacing between them.
     * Therefore, nodes are placed with the following spacing, which is hopefully enough to make up for size
     * changes caused by interactive layout.
     * Of course one could construct a model for which this spacing is not enough, but we choose to disregard this.
     * In all non constructed normal cases this spacing should be sufficient to ensure that nodes do not overlap.
     */
    public static final int PSEUDO_POSITION_SPACING = Integer.MAX_VALUE;
    
    /**
     * Constants used as the index of the last layer.
     * This is used since a {@code layerConstraint} property might be set that assigns a node to the {@code LAST} layer.
     * Since a {@code layerConstraint} should still be respected it is translated into a {@code layerChoiceConstraint}
     * using this constant as the value for the last layer.
     */
    public static final int LAST_LAYER_INDEX = Integer.MAX_VALUE;

    /**
     * Visits all nodes and sets interactive options for the {@code layered} algorithm.
     * This assumes that the {@code layered} algorithm is the default algorithm.
     */
    @Override
    public void visit(final ElkGraphElement element) {
        // Only apply to root of the graph
        if (element instanceof ElkNode) {
            ElkNode root = (ElkNode) element;
            setInteractiveOptionsAndPseudoPositions(root);
        }
    }

    /**
     * Sets pseudo positions and interactive strategies for the given graph.
     * 
     * @param root Root of the graph
     */
    private void setInteractiveOptionsAndPseudoPositions(final ElkNode root) {
        if (!root.getChildren().isEmpty()) {
            String algorithm = root.getProperty(CoreOptions.ALGORITHM);
            if (algorithm == null || LayeredOptions.ALGORITHM_ID.endsWith(algorithm)) {
                // Make sure that potential usages of the LayerConstraint are translated into
                // LayerChoiceConstraints
                for (ElkNode node : root.getChildren()) {
                    if (node.hasProperty(LayeredOptions.LAYERING_LAYER_CONSTRAINT)) {
                        LayerConstraint constraint = node.getProperty(LayeredOptions.LAYERING_LAYER_CONSTRAINT);
                        node.setProperty(LayeredOptions.LAYERING_LAYER_CONSTRAINT, LayerConstraint.NONE);
                        switch (constraint) {
                        case FIRST:
                            if (!node.hasProperty(LayeredOptions.LAYERING_LAYER_CHOICE_CONSTRAINT)) {
                                node.setProperty(LayeredOptions.LAYERING_LAYER_CHOICE_CONSTRAINT, 0);
                            }
                            break;
                        case LAST:
                            if (!node.hasProperty(LayeredOptions.LAYERING_LAYER_CHOICE_CONSTRAINT)) {
                                node.setProperty(LayeredOptions.LAYERING_LAYER_CHOICE_CONSTRAINT, LAST_LAYER_INDEX);
                            }
                            break;
                        }
                    }
                    setCoordinates(root);
                    setInteractiveStrategies(root);
                }
            }
        }
    }

    /**
     * Sets the coordinates of the nodes in the graph {@code root} according to the set constraints.
     * 
     * @param root The root of the graph that should be layouted.
     */
    private void setCoordinates(final ElkNode root) {
        List<List<ElkNode>> layers = calcLayerNodes(root.getChildren());
        Direction direction = root.getProperty(LayeredOptions.DIRECTION);
        setCoordinateInLayoutDirection(layers, direction);
        int layerId = 0;
        for (List<ElkNode> layer : layers) {
            if (layer.size() > 0) {
                setCoordinatesOrthogonalToLayoutDirection(layer, layerId, direction);
                layerId++;
            }
        }
        return;
    }

    /**
     * Calculates the layers the {@code nodes} belong to.
     * 
     * @param nodes The nodes of the graph for which the layers should be calculated.
     */
    private List<List<ElkNode>> calcLayerNodes(final List<ElkNode> nodes) {
        ArrayList<ElkNode> nodesWithoutC = new ArrayList<ElkNode>();
        ArrayList<ElkNode> nodesWithLayerConstraint = new ArrayList<ElkNode>();
        List<ElkNode> nodesWithRC = new ArrayList<ElkNode>();
        List<ElkNode> nodesWithRCAndLC = new ArrayList<ElkNode>();
        // Save the nodes in corresponding lists
        for (ElkNode node : nodes) {
            if (node.hasProperty(LayeredOptions.LAYERING_LAYER_CHOICE_CONSTRAINT)
                    && (node.hasProperty(LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_PRED_OF)
                            || node.hasProperty(LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_SUCC_OF))) {
                nodesWithRCAndLC.add(node);
            } else if (node.hasProperty(LayeredOptions.LAYERING_LAYER_CHOICE_CONSTRAINT)) {
                nodesWithLayerConstraint.add(node);
            } else if (node.hasProperty(LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_PRED_OF)
                    || node.hasProperty(LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_SUCC_OF)) {
                nodesWithRC.add(node);
            } else {
                nodesWithoutC.add(node);
            }
        }

        // Handle nodes with relative constraint and and layer constraint
        updateListsForLayerAssignment(nodesWithRCAndLC, nodes, nodesWithLayerConstraint, nodesWithRC, nodesWithoutC);
        
        nodesWithRC = sortRCNodes(nodesWithRC, nodes, nodesWithoutC);
        
        // Calculate layers for nodes without constraints based on their layerID,
        // which is set by the previous layout run.
        List<List<ElkNode>> layerNodes = initialLayers(nodesWithoutC);

        // Assign layers to nodes with constraints.
        int diff = assignLayersToNodesWithProperty(nodesWithLayerConstraint, layerNodes);

        // add nodes with rC to the correct layer
        assignLayersToNodesWithRC(nodesWithRC, layerNodes, nodes, diff);

        return layerNodes;
    }

    /**
     * Sets layer constraints for nodes that are referred by the nodes in {@code nodesWithRCAndLC} 
     * and adjusts the lists.
     * 
     * @param nodesWithRCAndLC
     *          nodes that have relative and layer constraints
     * @param nodes
     *          all nodes of the graph
     * @param nodesWithRC 
     *          nodes that have only relative constraints
     * @param nodesWithLayerConstraint 
     *          nodes that have only layer constraints
     * @param nodesWithoutC 
     *          nodes without any constraint
     */
    private void updateListsForLayerAssignment(final List<ElkNode> nodesWithRCAndLC, final List<ElkNode> nodes, 
            final ArrayList<ElkNode> nodesWithLayerConstraint, final List<ElkNode> nodesWithRC, 
            final ArrayList<ElkNode> nodesWithoutC) {        
        for (int i = 0; i < nodesWithRCAndLC.size(); i++) {
            ElkNode node = nodesWithRCAndLC.get(i);
            nodesWithLayerConstraint.add(node);
            List<ElkNode> targets = new ArrayList<ElkNode>();
            // get target nodes
            if (node.hasProperty(LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_PRED_OF)) {
                targets.add(getElkNode(nodes, 
                        node.getProperty(LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_PRED_OF)));
            } 
            if (node.hasProperty(LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_SUCC_OF)) {
                targets.add(getElkNode(nodes, 
                        node.getProperty(LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_SUCC_OF)));
            }
            
            for (ElkNode target : targets) {
                // set layer constraint of target node
                if (!target.hasProperty(LayeredOptions.LAYERING_LAYER_CHOICE_CONSTRAINT)) {
                    int val = node.getProperty(LayeredOptions.LAYERING_LAYER_CHOICE_CONSTRAINT);
                    target.setProperty(LayeredOptions.LAYERING_LAYER_CHOICE_CONSTRAINT, val);
                    // node should be in the correct lists
                    if (!nodesWithRC.contains(target)) {
                        nodesWithLayerConstraint.add(target);
                        nodesWithoutC.remove(target);
                    } else {
                        nodesWithRCAndLC.add(target);
                        nodesWithRC.remove(target);
                    }
                }
            }
        }
    }

    /**
     * Assign the nodes with relative constraints to layers.
     * 
     * @param nodesWithRC
     *          Nodes with relative constraints.
     * @param layering
     *          List that contains the layers with their corresponding nodes.
     * @param allNodes
     *          All nodes of the graph.
     */
    private void assignLayersToNodesWithRC(final List<ElkNode> nodesWithRC, final List<List<ElkNode>> layering, 
            final List<ElkNode> allNodes, final int diff) {
        for (int i = 0; i < nodesWithRC.size(); i++) {
            ElkNode node = nodesWithRC.get(i);
            ElkNode succNode = null;
            ElkNode predNode = null;
            String succName = node.getProperty(LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_PRED_OF);
            String predName = node.getProperty(LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_SUCC_OF);
            if (succName != null) {
                succNode = getElkNode(allNodes, succName);
            } 
            if (predName != null) {
                predNode = getElkNode(allNodes, predName);
            }
            
            IProperty<Integer> layProp = LayeredOptions.LAYERING_LAYER_CHOICE_CONSTRAINT;
            
            if (succNode != null && predNode != null) {
                // if both rel cons are set, it may be that the layer ids of 
                // some nodes of the chain must be uptdated too
                List<ElkNode> chain = getChainByAllNodes(node, allNodes);
                int succVal = succNode.hasProperty(layProp) ? succNode.getProperty(layProp) 
                        : succNode.getProperty(LayeredOptions.LAYERING_LAYER_ID);
                int predVal = predNode.hasProperty(layProp) ? predNode.getProperty(layProp) 
                        : predNode.getProperty(LayeredOptions.LAYERING_LAYER_ID);
                int val = succVal > predVal ? succVal : predVal;
                // update layer ids
                for (ElkNode n : chain) {
                    int oldVal = n.getProperty(LayeredOptions.LAYERING_LAYER_ID); 
                    if (oldVal != val) {
                        n.setProperty(LayeredOptions.LAYERING_LAYER_ID, val);
                        shiftOtherNodes(n, val, layering, true);
                        shiftOtherNodes(n, val, layering, false);
                    }
                    layering.get(oldVal).remove(n);
                    layering.get(val).add(n);
                }
            } else {
                // determine layer of referenced node
                int succVal = -1;
                if (succNode != null) {
                    succVal = succNode.hasProperty(layProp) ? succNode.getProperty(layProp) 
                            : succNode.getProperty(LayeredOptions.LAYERING_LAYER_ID);
                }
                int predVal = -1;
                if (predNode != null) {
                    predVal = predNode.hasProperty(layProp) ? predNode.getProperty(layProp) 
                            : predNode.getProperty(LayeredOptions.LAYERING_LAYER_ID);
                }
                
                // update layer id of the current node
                int val = succVal > predVal ? succVal : predVal;
                node.setProperty(LayeredOptions.LAYERING_LAYER_ID, val - diff);
                shiftOtherNodes(node, val - diff, layering, true);
                shiftOtherNodes(node, val - diff, layering, false);
                layering.get(val - diff).add(node);
            }
        }
    }

    /**
     * Adds the nodes with a layer constraint to the already assigned layered nodes based on their layer constraint.
     * 
     * @param nodesWithLayerConstraint
     *            Nodes with set layer constraint that should be added to the layers.
     * @param layering
     *            List that contains the layers with their corresponding nodes.
     */
    private int assignLayersToNodesWithProperty(final List<ElkNode> nodesWithLayerConstraint,
            final List<List<ElkNode>> layering) {
        // Sort nodes with constraint based on their layer.
        nodesWithLayerConstraint.sort((ElkNode a, ElkNode b) -> {
            return a.getProperty(LayeredOptions.LAYERING_LAYER_CHOICE_CONSTRAINT)
                    - b.getProperty(LayeredOptions.LAYERING_LAYER_CHOICE_CONSTRAINT);
        });

        // Add the nodes with constraints to the their desired layer.

        // diff keeps track of the difference between the layer the node should
        // be in and the layer the node is really in.
        // This way nodes with the same layer constraint go in the same layer
        // although it may not be the layer that is specified by the constraint.
        int diff = 0;
        for (ElkNode node : nodesWithLayerConstraint) {
            int currentLayer = node.getProperty(LayeredOptions.LAYERING_LAYER_CHOICE_CONSTRAINT) - diff;
            if (currentLayer < layering.size() && currentLayer >= 0) {
                List<ElkNode> nodesOfLayer = layering.get(currentLayer);
                // Shift nodes to remove in-layer edges.
                shiftOtherNodes(node, currentLayer, layering, true);
                shiftOtherNodes(node, currentLayer, layering, false);
                nodesOfLayer.add(node);
            } else if (currentLayer == -1) {
                layering.add(0, new ArrayList<>(Arrays.asList(node)));
            } else {
                diff = diff + currentLayer - layering.size();
                layering.add(new ArrayList<>(Arrays.asList(node)));
            }
        }
        return diff;
    }

    /**
     * Shifts nodes to the right such that edges in the same layer do not exist.
     * 
     * @param movedNode
     *            The node which connected nodes must be shifted .
     * @param layer
     *            The layer {@code moveNode} is in.
     * @param layerNodes
     *            All existing layers with the containing nodes.
     * @param incoming
     *            Determines if incoming or outgoing edges should be considered. True: incoming edges.
     */
    private void shiftOtherNodes(final ElkNode movedNode, final int layer, final List<List<ElkNode>> layerNodes,
            final boolean incoming) {
        List<ElkNode> nodesOfLayer = layerNodes.get(layer);
        // get edges
        List<ElkEdge> edges = new ArrayList<>();
        if (incoming) {
            ElkNode root = movedNode.getParent();
            for (ElkEdge edge : root.getContainedEdges()) {
                for (ElkConnectableShape target : edge.getTargets()) {
                    if (target.equals(movedNode)
                            || (target instanceof ElkPort && target.eContainer().equals(movedNode))) {
                        edges.add(edge);
                    }
                }
            }
        } else {
            ElkNode root = movedNode.getParent();
            for (ElkEdge edge : root.getContainedEdges()) {
                for (ElkConnectableShape target : edge.getSources()) {
                    if (target.equals(movedNode)
                            || (target instanceof ElkPort && target.eContainer().equals(movedNode))) {
                        edges.add(edge);
                    }
                }
            }
        }

        for (ElkEdge edge : edges) {
            // get connected node
            ElkNode node = null;
            if (incoming) {
                if (edge.getSources().get(0) instanceof ElkPort) {
                    node = (ElkNode) edge.getSources().get(0).eContainer();
                } else if (edge.getSources().get(0) instanceof ElkNode) {
                    node = (ElkNode) edge.getSources().get(0);
                }
            } else {
                if (edge.getTargets().get(0) instanceof ElkPort) {
                    node = (ElkNode) edge.getTargets().get(0).eContainer();
                } else if (edge.getTargets().get(0) instanceof ElkNode) {
                    node = (ElkNode) edge.getTargets().get(0);
                }
            }
            // shift node to the next layer
            if (nodesOfLayer.contains(node)) {
                nodesOfLayer.remove(node);
                List<ElkNode> newLayer;
                node.setProperty(LayeredOptions.LAYERING_LAYER_ID, layer + 1);
                if (layer + 1 < layerNodes.size()) {
                    newLayer = layerNodes.get(layer + 1);
                    newLayer.add(node);
                    // the connected nodes in the layer the node is shifted to must be shifted too
                    shiftOtherNodes(node, layer + 1, layerNodes, false);
                    shiftOtherNodes(node, layer + 1, layerNodes, true);
                } else {
                    layerNodes.add(new ArrayList<>(Arrays.asList(node)));
                }
            }
        }
    }

    /**
     * Sorts the {@code nodes} in layers based on their layerID that was set via the {@code ConstraintsPostProcessor}.
     * 
     * @param nodes
     *            The nodes of the graph which layers should be calculated.
     */
    private List<List<ElkNode>> initialLayers(final ArrayList<ElkNode> nodes) {
        // Sort by layerID.
        nodes.sort((ElkNode a, ElkNode b) -> {
            return a.getProperty(LayeredOptions.LAYERING_LAYER_ID) - b.getProperty(LayeredOptions.LAYERING_LAYER_ID);
        });

        List<List<ElkNode>> layerNodes = new ArrayList<List<ElkNode>>();
        List<ElkNode> nodesOfLayer = new ArrayList<ElkNode>();
        int layerId = 0;
        int currentLayer = -1;
        // Assign nodes to layers.
        for (ElkNode node : nodes) {
            int layer = node.getProperty(LayeredOptions.LAYERING_LAYER_ID);
            if (layer > currentLayer) {
                // Check if a node is added to a new layer.
                if (!nodesOfLayer.isEmpty()) {
                    layerNodes.add(nodesOfLayer);
                    layerId++;
                }
                nodesOfLayer = new ArrayList<ElkNode>();
                currentLayer = layer;
            }

            // Nodes with layer constraint should be ignored, since they are added later.
            if (!node.hasProperty(LayeredOptions.LAYERING_LAYER_CHOICE_CONSTRAINT)) {
                node.setProperty(LayeredOptions.LAYERING_LAYER_ID, layerId);
                nodesOfLayer.add(node);
            }
        }
        if (!nodesOfLayer.isEmpty()) {
            // Add the last layer nodes.
            layerNodes.add(nodesOfLayer);
        }
        return layerNodes;
    }

    /**
     * Sets the x coordinates of the nodes in {@code layers} according to their layer. The problem is that the height
     * and width of a node are determined by the height and width before the layout with constraints. They are therefore
     * potentially wrong. The placement has to be adjusted for this by a very high spacing between nodes.
     * 
     * @param layers
     *            The layers containing the associated nodes, already sorted regarding layers
     * @param direction
     *            The layout direction
     */
    private void setCoordinateInLayoutDirection(final List<List<ElkNode>> layers, final Direction direction) {
        double position = 0;
        double nextPosition = 0;
        // Assign x (RIGHT, LEFT)/y (UP, DOWN) coordinate such that all nodes get a pseudo position
        // that assigns them to the correct layer if interactive mode is used.
        for (List<ElkNode> nodesOfLayer : layers) {
            for (ElkNode node : nodesOfLayer) {
                switch (direction) {
                case UNDEFINED:
                case RIGHT:
                    node.setX(position);
                    if (position + node.getWidth() / 2 >= nextPosition) {
                        nextPosition = node.getX() + node.getWidth() + PSEUDO_POSITION_SPACING;
                    }
                    break;
                case LEFT:
                    node.setX(position);
                    if (node.getX() <= nextPosition) {
                        nextPosition = node.getX() - PSEUDO_POSITION_SPACING;
                    }
                    break;
                case DOWN:
                    node.setY(position);
                    if (position + node.getHeight() >= nextPosition) {
                        nextPosition = node.getY() + node.getHeight() + PSEUDO_POSITION_SPACING;
                    }
                    break;
                case UP:
                    node.setY(position);
                    if (node.getY() <= nextPosition) {
                        nextPosition = node.getY() - PSEUDO_POSITION_SPACING;
                    }
                    break;
                }
            }
            position = nextPosition;
        }
        return;
    }

    /**
     * Sets the positions of the nodes in their layer.
     * 
     * @param nodesOfLayer
     *            The list containing nodes that are in the same layer.
     * @param direction
     *            The layout direction.
     */
    private void setCoordinatesOrthogonalToLayoutDirection(final List<ElkNode> nodesOfLayer, final int layerId,
            final Direction direction) {
        // Separate nodes with and without constraints.
        List<ElkNode> allNodes = new ArrayList<ElkNode>();
        List<ElkNode> nodesWithoutC = new ArrayList<ElkNode>();
        List<ElkNode> nodesWithPCAndRC = new ArrayList<ElkNode>();
        List<ElkNode> nodesWithPC = new ArrayList<ElkNode>();
        List<ElkNode> nodesWithOneRC = new ArrayList<ElkNode>();
        List<ElkNode> nodesWithBothRC = new ArrayList<ElkNode>();
        
        for (ElkNode node : nodesOfLayer) {
            allNodes.add(node);
            if (node.hasProperty(LayeredOptions.CROSSING_MINIMIZATION_POSITION_CHOICE_CONSTRAINT)
                    && (node.hasProperty(LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_PRED_OF)
                            || node.hasProperty(LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_SUCC_OF))) {
                nodesWithPCAndRC.add(node);
            } else if (node.hasProperty(LayeredOptions.CROSSING_MINIMIZATION_POSITION_CHOICE_CONSTRAINT)) {
                nodesWithPC.add(node);
            } else if (node.hasProperty(LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_PRED_OF)
                    && node.hasProperty(LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_SUCC_OF)) {
                nodesWithBothRC.add(node);
            } else if (node.hasProperty(LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_PRED_OF)
                    || node.hasProperty(LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_SUCC_OF)) {
                nodesWithOneRC.add(node);
            } else {
                nodesWithoutC.add(node);
            }
        }
        
        // Sort nodes with constraint by their position id.
        nodesWithPC.sort((ElkNode a, ElkNode b) -> {
            return a.getProperty(LayeredOptions.CROSSING_MINIMIZATION_POSITION_ID)
                    - b.getProperty(LayeredOptions.CROSSING_MINIMIZATION_POSITION_ID);
        });
        
        // each node should only have one relative constraint
        for (int n = 0; n < nodesWithBothRC.size(); n++) {
            // pred cons is translated to succ cons for the following node
            ElkNode node = nodesWithBothRC.get(n);
            String t = node.getProperty(LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_PRED_OF);
            ElkNode target = getElkNode(allNodes, t);
            node.setProperty(LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_PRED_OF, null);
            if (target.hasProperty(LayeredOptions.CROSSING_MINIMIZATION_POSITION_CHOICE_CONSTRAINT)) {
                int posCosTarget = target.getProperty(LayeredOptions.CROSSING_MINIMIZATION_POSITION_CHOICE_CONSTRAINT);
                node.setProperty(LayeredOptions.CROSSING_MINIMIZATION_POSITION_CHOICE_CONSTRAINT, posCosTarget - 1);
                nodesWithPCAndRC.add(node);
            } else {
                target.setProperty(LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_SUCC_OF, node.getIdentifier());
                nodesWithOneRC.add(node);
                if (!target.hasProperty(LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_PRED_OF)) {
                    nodesWithoutC.remove(target);
                    nodesWithOneRC.add(target);
                } else {
                    nodesWithOneRC.remove(target);
                    nodesWithBothRC.add(n + 1, target);
                }
            }
        }
        
        // handle nodes with relative and position constraint
        // all nodes the are connected to the ones with the pos cons through relative constraints get also a pc
        handleNodesWithPCAndRC(allNodes, nodesWithoutC, nodesWithOneRC, nodesWithPC, nodesWithPCAndRC);

        // Determine the order of the nodes.
        nodesWithOneRC = sortRCNodes(nodesWithOneRC, allNodes, nodesWithoutC);
        sortNodesInLayer(nodesWithPC, nodesWithoutC, direction);
        
        // Add nodes with relative constraint
        for (ElkNode node : nodesWithOneRC) {
            ElkNode[] targets = new ElkNode[2];
            
            if (node.hasProperty(LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_PRED_OF)) {
                // pred is defined
                IProperty<String> nodeProp = LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_PRED_OF;
                targets[0] = getElkNode(allNodes, node.getProperty(nodeProp));
            } else {
                // succ is defined
                IProperty<String> nodeProp = LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_SUCC_OF;
                targets[1] = getElkNode(allNodes, node.getProperty(nodeProp));
            }
            
            for (int t = 0; t < targets.length; t++) {
                if (targets[t] != null) {
                    if (!targets[t].hasProperty(LayeredOptions.CROSSING_MINIMIZATION_POSITION_CHOICE_CONSTRAINT)) {
                        // add node at correct position
                        int index = t == 0 ? nodesWithoutC.indexOf(targets[t]) : nodesWithoutC.indexOf(targets[t]) + 1;
                        nodesWithoutC.add(index, node);
                    } else {
                        // set pos cons on node because target has pos cons 
                        int val =
                                targets[t].getProperty(LayeredOptions.CROSSING_MINIMIZATION_POSITION_CHOICE_CONSTRAINT);
                        int value = t == 0 ? val - 1 : val + 1;
                        IProperty<String> nodeProp = t == 0 ? LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_PRED_OF 
                                : LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_SUCC_OF;
                        node.setProperty(LayeredOptions.CROSSING_MINIMIZATION_POSITION_CHOICE_CONSTRAINT, value);
                        node.setProperty(nodeProp, null);
                        nodesWithPC.add(node);
                    }
                }
            }
        }
        
        // Sort nodes with constraint by their position constraint.
        nodesWithPC.sort((ElkNode a, ElkNode b) -> {
            return a.getProperty(LayeredOptions.CROSSING_MINIMIZATION_POSITION_CHOICE_CONSTRAINT)
                    - b.getProperty(LayeredOptions.CROSSING_MINIMIZATION_POSITION_CHOICE_CONSTRAINT);
        });
            
        // Add the nodes with position constraint at the desired position in their layer.
        for (ElkNode node : nodesWithPC) {
            int pos = node.getProperty(LayeredOptions.CROSSING_MINIMIZATION_POSITION_CHOICE_CONSTRAINT);
            if (pos < nodesWithoutC.size()) {
                if (pos > 0) {
                    swapNodes(nodesWithoutC, pos);
                }
                nodesWithoutC.add(pos, node);
            } else {
                nodesWithoutC.add(node);
            }
        }

        // Set the y/x positions according to the order of the nodes.
        switch (direction) {
        case UNDEFINED:
        case RIGHT:
        case LEFT:
            double yPos = nodesWithoutC.get(0).getY();
            for (ElkNode node : nodesWithoutC) {
                node.setProperty(LayeredOptions.POSITION, new KVector(node.getX(), yPos));
                yPos += node.getHeight() + PSEUDO_POSITION_SPACING;
            }
            break;
        case DOWN:
        case UP:
            double xPos = nodesWithoutC.get(0).getX() + 2 * layerId * PSEUDO_POSITION_SPACING;
            for (ElkNode node : nodesWithoutC) {
                node.setProperty(LayeredOptions.POSITION, new KVector(xPos, node.getY()));
                xPos += node.getWidth() + PSEUDO_POSITION_SPACING;
            }
            break;
        }
    }

    /**
     * Handles nodes with pos and rel cons by adding pos cons to all nodes that are connected to them by rel cons.
     * @param allNodes
     *          all nodes of the graph
     * @param nodes
     *          nodes without any constraint
     * @param nodesWithRC
     *          nodes with only relative constraints
     * @param nodesWithPositionConstraint
     *          nodes with only position constraints
     * @param nodesWithPCAndRC
     *          nodes with pos and rel constraints
     */
    private void handleNodesWithPCAndRC(final List<ElkNode> allNodes, final List<ElkNode> nodes, 
            final List<ElkNode> nodesWithRC, final List<ElkNode> nodesWithPositionConstraint, 
            final List<ElkNode> nodesWithPCAndRC) {
        for (int i = 0; i < nodesWithPCAndRC.size(); i++) {
            ElkNode node = nodesWithPCAndRC.get(i);
            nodesWithPositionConstraint.add(node);
            
            IProperty<String> nodeProp = null;
            List<ElkNode> targets = new ArrayList<>();
            List<Integer> vals = new ArrayList<>();
            int val = node.getProperty(LayeredOptions.CROSSING_MINIMIZATION_POSITION_CHOICE_CONSTRAINT);
            if (node.hasProperty(LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_PRED_OF)) {
                // pred is defined
                nodeProp = LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_PRED_OF;
                targets.add(getElkNode(allNodes, node.getProperty(nodeProp)));
                vals.add(val + 1);
                node.setProperty(nodeProp, null);
            } 
            if (node.hasProperty(LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_SUCC_OF)) {
                // succ is defined
                nodeProp = LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_SUCC_OF;
                targets.add(getElkNode(allNodes, node.getProperty(nodeProp)));
                vals.add(val - 1);
                node.setProperty(nodeProp, null);
            }
            
            // targets should be in the correct lists
            for (int t = 0; t < targets.size(); t++) {
                ElkNode target = targets.get(t);
                if (!nodesWithPositionConstraint.contains(target)) {
                    target.setProperty(LayeredOptions.CROSSING_MINIMIZATION_POSITION_CHOICE_CONSTRAINT, vals.get(t));
                    if (!nodesWithRC.contains(target)) {
                        nodesWithPositionConstraint.add(target);
                        nodes.remove(target);
                    } else {
                        nodesWithPCAndRC.add(target);
                        nodesWithRC.remove(target);
                    }
                }
            }
        }
    }

    /**
     * Swaps nodes such that a node that should be placed at {@code pos} does not interrupt nodes that
     * are connected with a relative constraint.
     * 
     * @param nodes Nodes in the current layer.
     * @param position Position a node should be placed at in the current layer.
     */
    private void swapNodes(final List<ElkNode> nodes, final int position) {
        // Get surrounding nodes for node put at the given position.
        ElkNode pred = nodes.get(position - 1);
        ElkNode succ = nodes.get(position);
        
        // Get all nodes that are connected by relative constraints to the predecessor.
        List<ElkNode> chain = getChainByLayerNodes(pred, nodes);
        if (chain.contains(succ)) {
            // Nodes only must be shifted if a relation between pred and succ exists

            // Number nodes that must be swapped
            int count = chain.indexOf(succ);
            // End of relative constraints chain
            int end = position + (chain.size() - count);
            
            // Determine nodes to swap
            List<ElkNode> swapNodes = new ArrayList<>();
            for (int i = 0; i < count; i++) {
                if (end >= nodes.size()) {
                    break;
                }
                ElkNode currentNode = nodes.get(end);
                chain = getChainByLayerNodes(currentNode, nodes);
                // Check whether the number of nodes linked by relative constraints
                if (chain.size() + i <= count) {
                    nodes.removeAll(chain);
                    swapNodes.addAll(chain);
                    i += chain.size() - 1;
                } else {
                    // skip chains of relative constraints
                    end += chain.size();
                }
            }
            
            // add swapNodes 
            nodes.addAll(position - count, swapNodes);
        }
    }

    /**
     * sorts {@code nodesWithRC} such that referred nodes that have also relative constraints come first.
     * 
     * @param nodesWithRC
     *          nodes with relative constraints
     * @param allNodes 
     *          all nodes of the graph
     * @param nodesWithoutC 
     *          nodes without any constraints
     * @return sorted nodes that have relative constraints
     */
    private List<ElkNode> sortRCNodes(final List<ElkNode> nodesWithRC, final List<ElkNode> allNodes, 
            final List<ElkNode> nodesWithoutC) {
        List<ElkNode> nodes = new ArrayList<>();
        while (!nodesWithRC.isEmpty()) {
            boolean circle = true;
            for (int i = 0; i < nodesWithRC.size(); i++) {
                ElkNode cur = nodesWithRC.get(i);
                if (cur.hasProperty(LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_PRED_OF)) {
                    // predecessor is defined
                    ElkNode target = getElkNode(allNodes, 
                            cur.getProperty(LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_PRED_OF));
                    if (!nodesWithRC.contains(target)) {
                        nodes.add(cur);
                        nodesWithRC.remove(i);
                        i--;
                        circle = false;
                    }
                } else {
                    // successor is defined
                    ElkNode target = getElkNode(allNodes, 
                            cur.getProperty(LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_SUCC_OF));
                    if (!nodesWithRC.contains(target)) {
                        nodes.add(cur);
                        nodesWithRC.remove(i);
                        i--;
                        circle = false;
                    }
                }
            }
            
            // if relative constraints create a circle, remove one of the constraints
            if (circle) {
                ElkNode cur = nodesWithRC.get(0);
                if (cur.hasProperty(LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_PRED_OF)) {
                    cur.setProperty(LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_PRED_OF, null);
                } else {
                    cur.setProperty(LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_SUCC_OF, null);
                }
                nodesWithRC.remove(0);
                nodesWithoutC.add(cur);
            }
        }
        return nodes;
    }

    /**
     * Sorts the {@code nodesWithPositionConstraint} according their position constraint and {@code nodes} according to
     * relevant coordinate (y for left, right, x for up, down).
     * 
     * @param nodesWithPositionConstraint
     *            The nodes with position constraint.
     * @param nodes
     *            The nodes without position constraints.
     * @param direction
     *            The layout direction.
     */
    private void sortNodesInLayer(final List<ElkNode> nodesWithPositionConstraint, final List<ElkNode> nodes,
            final Direction direction) {
        // Sort nodes with constraint by their position constraint.
        nodesWithPositionConstraint.sort((ElkNode a, ElkNode b) -> {
            return a.getProperty(LayeredOptions.CROSSING_MINIMIZATION_POSITION_CHOICE_CONSTRAINT)
                    - b.getProperty(LayeredOptions.CROSSING_MINIMIZATION_POSITION_CHOICE_CONSTRAINT);
        });
        // Sort other nodes based on their y(RIGHT/LEFT)/x(DOWN/UP) coordinate.
        nodes.sort((ElkNode a, ElkNode b) -> {
            switch (direction) {
            case UNDEFINED:
            case RIGHT:
            case LEFT:
                return (int) (a.getY() - b.getY());
            case DOWN:
            case UP:
                return (int) (a.getX() - b.getX());
            }
            return 0;
        });
    }

    /**
     * Sets the (semi) interactive strategies in the phases crossing minimization, layer assignment, cycle breaking for
     * the given parent node.
     * 
     * @param parent
     *            The graph which strategies should be set.
     */
    private void setInteractiveStrategies(final ElkNode parent) {
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_SEMI_INTERACTIVE, true);
        parent.setProperty(LayeredOptions.LAYERING_STRATEGY, LayeringStrategy.INTERACTIVE);
        parent.setProperty(LayeredOptions.CYCLE_BREAKING_STRATEGY, CycleBreakingStrategy.INTERACTIVE);
        // Disable model order for the final run, since it destroys the ordering created by the constraints.
        parent.setProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY, OrderingStrategy.NONE);
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_FORCE_NODE_MODEL_ORDER, false);
        parent.setProperty(LayeredOptions.CROSSING_MINIMIZATION_STRATEGY, CrossingMinimizationStrategy.LAYER_SWEEP);
    }
    
    /**
     * Gets the ElkNode with {@code id} as identifier.
     * 
     * @param nodes All nodes of the graph
     * @param id Id of the ElkNode that is searched
     * @return the ElkNode with the corresponding {@code id} or null if it could not be found.
     */
    private ElkNode getElkNode(final List<ElkNode> nodes, final String id) {
        for (ElkNode eN : nodes) {
            if (id.equals(eN.getIdentifier())) {
                return eN;
            }
        }
        return null;
    }

    /**
     * Determines the nodes that are connected to {@code node} by relative constraints.
     * The returned list of nodes is sorted based on the position of the nodes.
     * 
     * @param node One node of the chain
     * @param layerNodes Nodes that are in the same layer as {@code node}
     */
    private List<ElkNode> getChainByLayerNodes(final ElkNode node, final List<ElkNode> layerNodes) {
        int pos = layerNodes.indexOf(node);
        List<ElkNode> chainNodes = new ArrayList<ElkNode>();
        chainNodes.add(node);
        
        // from node to the start
        for (int i = pos - 1; i >= 0; i--) {
            if (layerNodes.get(i).hasProperty(LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_PRED_OF)
                || layerNodes.get(i + 1).hasProperty(LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_SUCC_OF)) {
                    chainNodes.add(0, layerNodes.get(i));
            } else {
                break;
            }
        }
        
        // count from node to the end
        for (int i = pos + 1; i < layerNodes.size(); i++) {
            if (layerNodes.get(i).hasProperty(LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_SUCC_OF)
                || layerNodes.get(i - 1).hasProperty(LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_PRED_OF)) {
                    chainNodes.add(layerNodes.get(i));
            } else {
                break;
            }
        }
        
        return chainNodes;
    }
    
    /**
     * Determines nodes that are connected to {@code node} by relative constraints.
     * @param node
     *          A node of the graph
     * @param allNodes
     *          All nodes of the graph
     * @return The nodes that are connected to {@code node} (unsorted).
     */
    private List<ElkNode> getChainByAllNodes(final ElkNode node, final List<ElkNode> allNodes) {
        List<ElkNode> chain = new ArrayList<>();
        chain.add(node);
        
        ElkNode succNode = null;
        ElkNode predNode = null;
        String succName = node.getProperty(LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_PRED_OF);
        String predName = node.getProperty(LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_SUCC_OF);
        while (succName != null) {
            succNode = getElkNode(allNodes, succName);
            chain.add(succNode);
            succName = succNode.getProperty(LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_PRED_OF);
        }
        while (predName != null) {
            predNode = getElkNode(allNodes, predName);
            chain.add(predNode);
            predName = predNode.getProperty(LayeredOptions.CROSSING_MINIMIZATION_IN_LAYER_PRED_OF);
        }
        
        return chain;
    }
    
}
