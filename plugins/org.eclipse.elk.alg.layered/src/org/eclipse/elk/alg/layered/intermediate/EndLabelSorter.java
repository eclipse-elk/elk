/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.eclipse.elk.alg.common.nodespacing.cellsystem.LabelCell;
import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.options.EdgeLabelPlacement;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.core.util.adapters.GraphAdapters.LabelAdapter;

/**
 * Sorts end labels according to the order of nodes their respective edges come from or head to. This is to prevent
 * seemingly random label orders in diagrams that employ end labels and have multiple edges going into or leaving the
 * same port.
 * 
 * <dl>
 *   <dt>Precondition:</dt>
 *     <dd>There are no long edge dummies left.</dd>
 *     <dd>There are no label dummies left.</dd>
 *     <dd>There are no north / south port dummies left.</dd>
 *   <dt>Postcondition:</dt>
 *     <dd>End labels at each port are sorted according to order of source / target nodes.</dd>
 *     <dd>Center edge labels are sorted according to order of source / target nodes.</dd>
 *   <dt>Slots:</dt>
 *     <dd>After phase 5.</dd>
 *   <dt>Same-slot dependencies:</dt>
 *     <dd>{@link LongEdgeJoiner}</dd>
 *     <dd>{@link NorthSouthPortPostprocessor}</dd>
 * </dl>
 */
public final class EndLabelSorter implements ILayoutProcessor<LGraph> {
    
    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor monitor) {
        monitor.begin("Sort end labels", 1);
        
        layeredGraph.getLayers().stream()
            .flatMap(layer -> layer.getNodes().stream())
            .filter(node -> node.getType() == NodeType.NORMAL)
            .forEach(node -> processNode(node));
        
        monitor.done();
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Node Processing and Initialization
    
    private void processNode(final LNode node) {
        boolean initializeMethodCalled = false;
        
        if (node.hasProperty(InternalProperties.END_LABELS)) {
            Map<LPort, LabelCell> labelCellMap = node.getProperty(InternalProperties.END_LABELS);
            
            // Iterate over all ports and check for each port if it requires its labels to be sorted
            for (LPort port : node.getPorts()) {
                if (needsSorting(port)) {
                    // Check if we need to initialize
                    if (!initializeMethodCalled) {
                        initialize(node.getGraph());
                        initializeMethodCalled = true;
                    }
                    
                    sort(port, labelCellMap.get(port));
                }
            }
        }
    }

    /**
     * A port requires its end labels to be sorted if there are end labels of at least two edges there.
     */
    private boolean needsSorting(final LPort port) {
        int edgesWithEndLabels = 0;
        
        for (LEdge inEdge : port.getIncomingEdges()) {
            boolean headLabels = inEdge.getLabels().stream()
                .anyMatch(label -> label.getProperty(LayeredOptions.EDGE_LABELS_PLACEMENT) == EdgeLabelPlacement.HEAD);
            if (headLabels) {
                edgesWithEndLabels++;
            }
        }
        
        for (LEdge outEdge : port.getOutgoingEdges()) {
            boolean tailLabels = outEdge.getLabels().stream()
                .anyMatch(label -> label.getProperty(LayeredOptions.EDGE_LABELS_PLACEMENT) == EdgeLabelPlacement.TAIL);
            if (tailLabels) {
                edgesWithEndLabels++;
            }
        }
        
        return edgesWithEndLabels >= 2;
    }
    
    /**
     * Called once we find the first instance of labels that have to be sorted. This method initializes everything 
     */
    private void initialize(final LGraph lGraph) {
        // Give nodes and ports ascending IDs
        int nextElementID = 0;
        for (Layer layer : lGraph) {
            for (LNode node : layer) {
                node.id = nextElementID++;
                
                for (LPort port : node.getPorts()) {
                    port.id = nextElementID++;
                }
            }
        }
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Sorting
    
    /**
     * Sorts the labels of the given port which are contained in the given label cell.
     */
    private void sort(final LPort port, final LabelCell portLabelCell) {
        List<LabelGroup> labelGroups = createLabelGroups(portLabelCell);
        labelGroups.sort(LABEL_GROUP_COMPARATOR);
        
        // Re-add the label cell's labels in the proper order. We directly access the cell's label list here, which
        // won't case the cell to recompute its size
        List<LabelAdapter<?>> portLabelCellLabels = portLabelCell.getLabels();
        portLabelCellLabels.clear();
        for (LabelGroup group : labelGroups) {
            portLabelCellLabels.addAll(group.labels);
        }
    }
    
    
    /**
     * Creates a list of {@link LabelGroup}s that group labels from the same edge.
     */
    private List<LabelGroup> createLabelGroups(final LabelCell portLabelCell) {
        Map<LEdge, LabelGroup> edgeToGroupMap = new HashMap<>();
        
        // Make sure every label is contained in a label group
        for (LabelAdapter<?> label : portLabelCell.getLabels()) {
            LEdge edge = label.getProperty(InternalProperties.END_LABEL_EDGE);
            
            // Create a new label group for the edge if we see it for the first time. Note that this works even if the
            // label didn't have an END_LABEL_EDGE property set
            if (!edgeToGroupMap.containsKey(edge)) {
                edgeToGroupMap.put(edge, new LabelGroup(edge));
            }
            
            edgeToGroupMap.get(edge).labels.add(label);
        }
        
        return new ArrayList<>(edgeToGroupMap.values());
    }


    /**
     * A group of labels that have a certain order and belong to a single edge. Label groups know the ultimate source
     * and target of the edge they belong to, and can be sorted. Since this is an internal class, we don't bother with
     * accessors.
     */
    private static class LabelGroup {
        
        /** The edge the labels belong to. This can be a dummy edge if the original edge was broken by a label dummy. */
        private final LEdge edge;
        /** List of labels that belong to this group. The order of these labels will never be changed. */
        private final List<LabelAdapter<?>> labels = new ArrayList<>();
        
        
        LabelGroup(final LEdge edge) {
            this.labels.addAll(labels);
            this.edge = edge;
        }
        
    }
    
    private static final Comparator<LabelGroup> LABEL_GROUP_COMPARATOR = new Comparator<LabelGroup>() {

        @Override
        public int compare(final LabelGroup group1, final LabelGroup group2) {
            // If they are not connected to the same source port, use the difference as a basis for the comparison
            int sourcePortDiff = Integer.compare(
                    group1.edge.getSource().id,
                    group2.edge.getSource().id);
            if (sourcePortDiff != 0) {
                return sourcePortDiff;
            }
            
            // They are connected to the same source port. Sort by target node, if that makes a difference.
            int targetNodeDiff = Integer.compare(
                    group1.edge.getTarget().getNode().id,
                    group2.edge.getTarget().getNode().id);
            if (targetNodeDiff != 0) {
                return targetNodeDiff;
            }
            
            // They are connected to the same source port and to the same target node. Compare target ports, but
            // backwards. Since western ports are ordered from bottom to top, not sorting backwards would yield the
            // opposite of our desires.
            return Integer.compare(
                    group2.edge.getTarget().id,
                    group1.edge.getTarget().id);
        }
        
    };
    
}
