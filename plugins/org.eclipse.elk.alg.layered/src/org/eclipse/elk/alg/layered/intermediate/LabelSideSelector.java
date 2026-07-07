/*******************************************************************************
 * Copyright (c) 2012, 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate;

import java.util.ArrayDeque;
import java.util.Collection;
import java.util.Deque;
import java.util.List;
import java.util.Queue;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LLabel;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.EdgeLabelSideSelection;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.options.LabelSide;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.IElkProgressMonitor;

import com.google.common.collect.Lists;

/**
 * <p>Decides for each edge label whether to place it above or below its respective edge. How the decision is made
 * depends on the active edge label side selection strategy.</p>
 * 
 * <dl>
 *   <dt>Precondition:</dt>
 *     <dd>a properly layered graph with fixed port orders</dd>
 *     <dd>port order is reflected in the way port lists are sorted</dd>
 *     <dd>center edge labels are not attached to edges, but to label dummy nodes</dd>
 *   <dt>Postcondition:</dt>
 *     <dd>the placement side is chosen for each label, and each label is annotated accordingly. (in the case of center
 *       edge labels, the representing dummy node is annotated with the placement side.)</dd>
 *     <dd>label dummy nodes have their ports placed such that they extend above or below their edge, depending on the
 *       side of their label.</dd>
 *   <dt>Slots:</dt>
 *     <dd>Before phase 4.</dd>
 *   <dt>Same-slot dependencies:</dt>
 *     <dd>{@link InLayerConstraintProcessor}</dd>
 *     <dd>{@link LabelManagementProcessor}</dd>
 *     <dd>{@link LabelDummySwitcher}</dd>
 * </dl>
 * 
 * @see EdgeLabelSideSelection
 */
public final class LabelSideSelector implements ILayoutProcessor<LGraph> {
    
    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor monitor) {
        EdgeLabelSideSelection mode = layeredGraph.getProperty(LayeredOptions.EDGE_LABELS_SIDE_SELECTION);
        monitor.begin("Label side selection (" + mode + ")", 1);
        
        // Calculate all label sides depending on the given strategy
        switch (mode) {
        case ALWAYS_UP:
            sameSide(layeredGraph, LabelSide.ABOVE);
            break;
        case ALWAYS_DOWN:
            sameSide(layeredGraph, LabelSide.BELOW);
            break;
        case DIRECTION_UP:
            basedOnDirection(layeredGraph, LabelSide.ABOVE);
            break;
        case DIRECTION_DOWN:
            basedOnDirection(layeredGraph, LabelSide.BELOW);
            break;
        case SMART_UP:
            smart(layeredGraph, LabelSide.ABOVE);
            break;
        case SMART_DOWN:
            smart(layeredGraph, LabelSide.BELOW);
            break;
        }

        monitor.done();
    }
    
    
    ////////////////////////////////////////////////////////////////////////////////////////
    // Simple Placement Strategies
    
    /**
     * Configures all labels to be placed on the given side.
     */
    private void sameSide(final LGraph graph, final LabelSide labelSide) {
        for (Layer layer : graph) {
            for (LNode node : layer) {
                if (node.getType() == NodeType.LABEL) {
                    applyLabelSide(node, labelSide);
                }
                
                for (LEdge edge : node.getOutgoingEdges()) {
                    applyLabelSide(edge, labelSide);
                }
            }
        }
    }
    
    /**
     * Configures all labels to be placed according to their edge's direction. If their edge points right, the labels
     * will be placed on the side passed to this method. Otherwise, they will be placed on the opposite side.
     */
    private void basedOnDirection(final LGraph graph, final LabelSide sideForRightwardEdges) {
        for (Layer layer : graph) {
            for (LNode node : layer) {
                if (node.getType() == NodeType.LABEL) {
                    LabelSide side = doesEdgePointRight(node)
                            ? sideForRightwardEdges
                            : sideForRightwardEdges.opposite();
                    applyLabelSide(node, side);
                }
                
                for (LEdge edge : node.getOutgoingEdges()) {
                    LabelSide side = doesEdgePointRight(edge)
                            ? sideForRightwardEdges
                            : sideForRightwardEdges.opposite();
                    applyLabelSide(edge, side);
                }
            }
        }
    }
    
    
    ////////////////////////////////////////////////////////////////////////////////////////
    // Smart Placement Strategy
    
    /**
     * Chooses label sides depending on certain patterns. If in doubt, uses the given default side.
     */
    private void smart(final LGraph graph, final LabelSide defaultSide) {
        // We will collect consecutive runs of certain dummy nodes while we iterate through layers
        Deque<LNode> dummyNodeQueue = new ArrayDeque<>();
        
        for (Layer layer : graph) {
            // The first call to any method
            boolean topGroup = true;
            int labelDummiesInQueue = 0;
            
            for (LNode node : layer) {
                switch (node.getType()) {
                case LABEL:
                    labelDummiesInQueue++;
                    // Intended fall-through to add the label dummy to the queue
                case LONG_EDGE:
                    dummyNodeQueue.add(node);
                    break;
                    
                case NORMAL:
                    smartForRegularNode(node, defaultSide);
                    
                    // Intended fall-through to handle the most recent dummy node run, if any
                    
                default:
                    // Empty dummy node queue
                    if (!dummyNodeQueue.isEmpty()) {
                        smartForConsecutiveDummyNodeRun(
                                dummyNodeQueue, labelDummiesInQueue, topGroup, false, defaultSide);
                    }
                    
                    // Reset things
                    topGroup = false;
                    labelDummiesInQueue = 0;
                }
            }
            
            // Do stuff with the nodes in the queue
            if (!dummyNodeQueue.isEmpty()) {
                smartForConsecutiveDummyNodeRun(
                        dummyNodeQueue, labelDummiesInQueue, topGroup, true, defaultSide);
            }
        }
    }
    
    /**
     * Assigns label sides to all label dummies in the given queue and empties the queue afterwards. The queue is
     * expected to not be empty. We also want to know whether this group consists of the topmost or of the bottommost
     * nodes in their layer (both can be true simultaneously). Decisions will be made, and they might depend on these
     * information.
     */
    private void smartForConsecutiveDummyNodeRun(final Deque<LNode> dummyNodes, final int labelDummyCount,
            final boolean topGroup, final boolean bottomGroup, final LabelSide defaultSide) {
        
        assert !dummyNodes.isEmpty();
        
        // We distinguish a number of special cases whose rules seem rather complicated. 
        if (topGroup
                && (!bottomGroup || dummyNodes.size() > 1)
                && labelDummyCount == 1
                && dummyNodes.peek().getType() == NodeType.LABEL) {
            
            // The current run of dummy nodes is at the top of the layer, has only a single label dummy, and that label
            // dummy is at the top of the run; select the ABOVE side to ensure that its edge doesn't get too long
            applyLabelSide(dummyNodes.peek(), LabelSide.ABOVE);
            
        } else if (bottomGroup
                && (!topGroup || dummyNodes.size() > 1)
                && labelDummyCount == 1
                && dummyNodes.peekLast().getType() == NodeType.LABEL) {
            
            // Symmetric to the previous case, only at the bottom of the layer
            applyLabelSide(dummyNodes.peekLast(), LabelSide.BELOW);
            
        } else if (dummyNodes.size() == 2) {
            // There's only a run of two edges, so place the label of the first above its edge, and the label of the
            // second below its edge
            applyLabelSide(dummyNodes.poll(), LabelSide.ABOVE);
            applyLabelSide(dummyNodes.poll(), LabelSide.BELOW);
            
        } else {
            // Not one of the special cases above. Iterate over the dummy nodes and assign the default label side,
            // except if we find two consecutive label dummies that connect the same nodes (tight loops that can
            // happen in control flow diagrams)
            applyForDummyNodeRunWithSimpleLoops(dummyNodes, labelDummyCount, defaultSide);
        }
        
        // Ensure the list is cleared
        dummyNodes.clear();
    }
    
    /**
     * Takes a collection of dummy nodes that contains a given number of label dummy nodes and applies label sides to
     * label dummy nodes in the collection. Usually, that will be the default side. If, however, we find exactly two
     * consecutive label dummies that connect the same regular nodes, that's going to be a special case.
     */
    private void applyForDummyNodeRunWithSimpleLoops(final Collection<LNode> dummyNodes, final int labelDummyCount,
            final LabelSide defaultSide) {
        
        // We keep track of runs of consecutive label dummy nodes that connect the same two nodes
        List<LNode> labelDummyRun = Lists.newArrayListWithCapacity(dummyNodes.size());
        LNode prevLongEdgeSource = null;
        LNode prevLongEdgeTarget = null;
        
        for (LNode currentDummy : dummyNodes) {
            assert currentDummy.getType() == NodeType.LABEL || currentDummy.getType() == NodeType.LONG_EDGE;
            
            // Check if we are continuing a previous run
            LNode currLongEdgeSource = getLongEdgeEndNode(currentDummy, true);
            LNode currLongEdgeTarget = getLongEdgeEndNode(currentDummy, false);
            
            if (prevLongEdgeSource != currLongEdgeSource || prevLongEdgeTarget != currLongEdgeTarget) {
                // We're starting a new run
                applyLabelSidesToLabelDummyRun(labelDummyRun, defaultSide);
                
                prevLongEdgeSource = currLongEdgeSource;
                prevLongEdgeTarget = currLongEdgeTarget;
            }
            
            labelDummyRun.add(currentDummy);
        }
        
        // Assign label sides to whatever dummy nodes are left
        applyLabelSidesToLabelDummyRun(labelDummyRun, defaultSide);
    }
    
    /**
     * Returns either the long edge source or target node of the given label dummy. May be, but shouldn't be
     * {@code null}.
     */
    private LNode getLongEdgeEndNode(final LNode labelDummy, final boolean source) {
        LPort endPort = labelDummy.getProperty(source
                ? InternalProperties.LONG_EDGE_SOURCE
                : InternalProperties.LONG_EDGE_TARGET);
        
        return endPort == null
                ? null
                : endPort.getNode();
    }
    
    /**
     * Applies label sides to the given list of consecutive dummy nodes and empties that list afterwards.
     */
    private void applyLabelSidesToLabelDummyRun(final List<LNode> labelDummyRun, final LabelSide defaultSide) {
        if (!labelDummyRun.isEmpty()) {
            // If the list contains exactly two label dummies, we place labels differently
            if (labelDummyRun.size() == 2) {
                applyLabelSide(labelDummyRun.get(0), LabelSide.ABOVE);
                applyLabelSide(labelDummyRun.get(1), LabelSide.BELOW);
            } else {
                for (LNode dummyNode : labelDummyRun) {
                    applyLabelSide(dummyNode, defaultSide);
                }
            }
            
            labelDummyRun.clear();
        }
    }

    /**
     * Assigns label sides to all end labels incident to this node. The assigned label sides depend on how many ports
     * there are on any given side.
     */
    private void smartForRegularNode(final LNode node, final LabelSide defaultSide) {
        // Iterate over the node's list of ports on each side. Remember the ones that have edges connected to them
        // and make the label side decision based on how many such ports there are
        Queue<List<LLabel>> endLabelQueue = new ArrayDeque<>(node.getPorts().size());
        PortSide currentPortSide = null;
        
        // This is where we assume that the list of ports is properly sorted
        for (LPort port : node.getPorts()) {
            if (port.getSide() != currentPortSide) {
                if (!endLabelQueue.isEmpty()) {
                    smartForRegularNodePortEndLabels(endLabelQueue, currentPortSide, defaultSide);
                }
                
                endLabelQueue.clear();
                currentPortSide = port.getSide();
            }
            
            // Possibly add the port's end labels to our queue, if it has any
            List<LLabel> portEndLabels = EndLabelPreprocessor.gatherLabels(port);
            if (portEndLabels != null) {
                endLabelQueue.add(portEndLabels);
            }
        }

        // Clear remaining ports
        if (!endLabelQueue.isEmpty()) {
            smartForRegularNodePortEndLabels(endLabelQueue, currentPortSide, defaultSide);
        }
    }
    
    /**
     * Handle the end labels currently in the queue (which may get modified in the process).
     */
    private void smartForRegularNodePortEndLabels(final Queue<List<LLabel>> endLabelQueue, final PortSide portSide,
            final LabelSide defaultSide) {
        
        assert !endLabelQueue.isEmpty();
        assert portSide != null;
        
        if (endLabelQueue.size() == 2) {
            // What we're going to do depends on which port side we are traversing...
            if (portSide == PortSide.NORTH || portSide == PortSide.EAST) {
                applyLabelSide(endLabelQueue.poll(), LabelSide.ABOVE);
                applyLabelSide(endLabelQueue.poll(), LabelSide.BELOW);
            } else {
                applyLabelSide(endLabelQueue.poll(), LabelSide.BELOW);
                applyLabelSide(endLabelQueue.poll(), LabelSide.ABOVE);
            }
        } else {
            for (List<LLabel> labelList : endLabelQueue) {
                applyLabelSide(labelList, defaultSide);
            }
        }
    }
    
    
    ////////////////////////////////////////////////////////////////////////////////////////
    // Helper Methods
    
    /**
     * Applies the given label side to the given label dummy node. If necessary, its ports are moved to reserve space
     * for the label on the correct side.
     */
    private void applyLabelSide(final LNode labelDummy, final LabelSide side) {
        // This method only does things to label dummy nodes
        if (labelDummy.getType() == NodeType.LABEL) {
            LabelSide effectiveSide = labelDummy.isInlineEdgeLabel()
                    ? LabelSide.INLINE
                    : side;

            labelDummy.setProperty(InternalProperties.LABEL_SIDE, effectiveSide);
            
            // If the label is not below the edge, the ports need to be moved
            if (effectiveSide != LabelSide.BELOW) {
                LEdge originEdge = (LEdge) labelDummy.getProperty(InternalProperties.ORIGIN);
                double thickness = originEdge.getProperty(LayeredOptions.EDGE_THICKNESS);
                
                // The new port position depends on the new placement
                double portPos = 0;
                if (effectiveSide == LabelSide.ABOVE) {
                    portPos = labelDummy.getSize().y - Math.ceil(thickness / 2);
                } else if (effectiveSide == LabelSide.INLINE) {
                    // The label dummy has a superfluous label-edge spacing
                    portPos = Math.ceil((labelDummy.getSize().y
                                - labelDummy.getGraph().getProperty(LayeredOptions.SPACING_EDGE_LABEL)
                                - thickness))
                            / 2.0;
                    labelDummy.getSize().y -= labelDummy.getGraph().getProperty(LayeredOptions.SPACING_EDGE_LABEL);
                    labelDummy.getSize().y -= thickness;
                }
                
                for (LPort port : labelDummy.getPorts()) {
                    port.getPosition().y = portPos;
                }
            }
        }
    }
    
    /**
     * Applies the given label side to all labels of the given edge.
     */
    private void applyLabelSide(final LEdge edge, final LabelSide side) {
        for (LLabel label : edge.getLabels()) {
            label.setProperty(InternalProperties.LABEL_SIDE, side);
        }
    }
    
    /**
     * Applies the given label side to all labels in the list.
     */
    private void applyLabelSide(final List<LLabel> labels, final LabelSide side) {
        for (LLabel label : labels) {
            label.setProperty(InternalProperties.LABEL_SIDE, side);
        }
    }
    
    /**
     * Checks if the given particular edge will point right or left in the final drawing. The edge can
     * either be a regular edge or part of a long edge. In the latter case, this edge is assumed to
     * point left if the original long edge target is left of the original long edge source, and if this
     * edge is marked to have been reversed.
     * 
     * @param edge the edge to check.
     * @return {@code true} if the edge will point right in the final drawing.
     */
    private boolean doesEdgePointRight(final LEdge edge) {
        // cds: This method formerly used a much more complicated implementation that took the original
        // long edge source and target into account. However, this short implementation should be quite
        // enough...
        return !edge.getProperty(InternalProperties.REVERSED);
    }
    
    /**
     * Checks if the given label dummy node is part of an edge segment that will point right in the
     * final drawing. We assume that this is the case if at least one of the incident edges points
     * rightwards according to {@link #doesEdgePointRight(LEdge)}.
     * 
     * @param labelDummy the label dummy to check.
     * @return {@code true} if this label dummy node is part of an edge that will point right in the
     *         final drawing.
     */
    private boolean doesEdgePointRight(final LNode labelDummy) {
        assert labelDummy.getType() == NodeType.LABEL;
        assert labelDummy.getIncomingEdges().iterator().hasNext();
        assert labelDummy.getOutgoingEdges().iterator().hasNext();
        
        // Find incoming and outgoing edge
        LEdge incoming = labelDummy.getIncomingEdges().iterator().next();
        LEdge outgoing = labelDummy.getOutgoingEdges().iterator().next();
        
        return doesEdgePointRight(incoming) || doesEdgePointRight(outgoing);
    }

}
