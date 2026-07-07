/*******************************************************************************
 * Copyright (c) 2010, 2018 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate;

import java.util.List;
import java.util.ListIterator;
import java.util.stream.Collectors;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.IElkProgressMonitor;

/**
 * Looks for long edge dummy nodes that can be joined together. The aim is to reduce the
 * number of edges by having edges originating from the same port or going into the same
 * port joined. This should be done after crossing minimization. Only those dummy nodes
 * are joined that the crossing minimizer placed right next to each other.
 * 
 * Merging the adjacent long edge dummies is allowed if one of the following holds true:
 * <ul>
 * <li>Neither of the two originating edges has {@link NodeType#LABEL} dummies <b>and</b> the two label dummies 
 * belong to the same hyperedge as computed by {@link #identifyHyperedges(LGraph)}</li>
 * <li>The two long edge dummies have the same {@link InternalProperties#LONG_EDGE_SOURCE LONG_EDGE_SOURCE} 
 *     <b>and</b> neither of the two has {@link NodeType#LABEL} dummies in a preceeding layer of the 
 *     current long edge dummy.
 *     <b>OR</b> the symmetric case for {@link InternalProperties#LONG_EDGE_TARGET} a
 *     nd label dummies in suceeding layers. 
 * </ul>
 * 
 * <dl>
 *   <dt>Precondition:</dt>
 *     <dd>a layered graph</dd>
 *     <dd>node orders are fixed</dd>
 *     <dd>for long edge dummies {@link InternalProperties#LONG_EDGE_SOURCE} and
 *         {@link InternalProperties#LONG_EDGE_TARGET} properties must be set.</dd>
 *   <dt>Postcondition:</dt>
 *     <dd>long edge dummy nodes belonging to the same hyperedge and being directly next to each other
 *         are merged subject to the condition discussed above.</dd>
 *   <dt>Slots:</dt>
 *     <dd>Before phase 4.</dd>
 *   <dt>Same-slot dependencies:</dt>
 *     <dd>{@link InLayerConstraintProcessor}</dd>
 *     <dd>{@link LabelDummySwitcher}</dd>
 * </dl>
 */
public final class HyperedgeDummyMerger implements ILayoutProcessor<LGraph> {

    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor monitor) {
        monitor.begin("Hyperedge merging", 1);
        
        identifyHyperedges(layeredGraph);
        
        // Iterate through the layers
        ListIterator<Layer> layerIter = layeredGraph.getLayers().listIterator();
        while (layerIter.hasNext()) {
            Layer layer = layerIter.next();
            List<LNode> nodes = layer.getNodes();
            
            // If there are no nodes anyway, just move on to the next layer
            if (nodes.isEmpty()) {
                continue;
            }
            
            LNode currNode = null;
            NodeType currNodeType = null;
            LNode lastNode = null;
            NodeType lastNodeType = null;
            
            // Iterate through the remaining nodes
            for (int nodeIndex = 0; nodeIndex < nodes.size(); nodeIndex++) {
                // Get the next node
                currNode = nodes.get(nodeIndex);
                currNodeType = currNode.getType();
                
                // We're only interested if the current and last nodes are long edge dummies
                if (currNodeType == NodeType.LONG_EDGE && lastNodeType == NodeType.LONG_EDGE) {
                    
                    // If the source or the target are identical and we are allowed to merge, merge the current node
                    // into the last
                    MergeState state = checkMergeAllowed(currNode, lastNode);
                    if (state.allowMerge) {
                        mergeNodes(currNode, lastNode, state.sameSource, state.sameTarget);
                        
                        // Remove the current node and make the last node the current node
                        nodes.remove(nodeIndex);
                        nodeIndex--;
                        currNode = lastNode;
                        currNodeType = lastNodeType;
                    }
                }
                
                // Remember this node for the next iteration
                lastNode = currNode;
                lastNodeType = currNodeType;
            }
        }
        
        monitor.done();
    }
    
    /**
     * @return whether merging the two long edge dummies is eligible as specified in the class's javadoc.
     */
    private MergeState checkMergeAllowed(final LNode currNode, final LNode lastNode) {

        boolean currHasLabelDummies = currNode.getProperty(InternalProperties.LONG_EDGE_HAS_LABEL_DUMMIES);
        boolean lastHasLabelDummies = lastNode.getProperty(InternalProperties.LONG_EDGE_HAS_LABEL_DUMMIES);
        
        // Get long edge source and target ports
        LPort currNodeSource = currNode.getProperty(InternalProperties.LONG_EDGE_SOURCE);
        LPort lastNodeSource = lastNode.getProperty(InternalProperties.LONG_EDGE_SOURCE);
        LPort currNodeTarget = currNode.getProperty(InternalProperties.LONG_EDGE_TARGET);
        LPort lastNodeTarget = lastNode.getProperty(InternalProperties.LONG_EDGE_TARGET);
        
        // Find out whether the two long edges come from the same source or target (non-null!)
        boolean sameSource = currNodeSource != null && currNodeSource == lastNodeSource;
        boolean sameTarget = currNodeTarget != null && currNodeTarget == lastNodeTarget;

        if (!currHasLabelDummies && !lastHasLabelDummies) {
            // assumption: long edge dummies always have two ports, both have the same id
            return new MergeState(
                    currNode.getPorts().iterator().next().id == lastNode.getPorts().iterator().next().id,
                    sameSource, sameTarget);
        }
        
        // If we can merge on grounds of the two long edges having the same source, we need to be
        // in front of the label dummy of each of the two edges (if any)
        boolean eligibleForSourceMerging =
                (!currNode.getProperty(InternalProperties.LONG_EDGE_HAS_LABEL_DUMMIES)
                        || currNode.getProperty(InternalProperties.LONG_EDGE_BEFORE_LABEL_DUMMY))
                &&
                (!lastNode.getProperty(InternalProperties.LONG_EDGE_HAS_LABEL_DUMMIES)
                        || lastNode.getProperty(InternalProperties.LONG_EDGE_BEFORE_LABEL_DUMMY));
        
        // If we can merge on grounds of the two long edges having the same target, we need to be
        // behind the label dummy of each of the two edges (if any)
        boolean eligibleForTargetMerging =
                (!currNode.getProperty(InternalProperties.LONG_EDGE_HAS_LABEL_DUMMIES)
                        || !currNode.getProperty(InternalProperties.LONG_EDGE_BEFORE_LABEL_DUMMY))
                &&
                (!lastNode.getProperty(InternalProperties.LONG_EDGE_HAS_LABEL_DUMMIES)
                        || !lastNode.getProperty(InternalProperties.LONG_EDGE_BEFORE_LABEL_DUMMY));
        
        return new MergeState(
                (sameSource && eligibleForSourceMerging) || (sameTarget && eligibleForTargetMerging), 
                sameSource, sameTarget);
    }
    
    // SUPPRESS CHECKSTYLE NEXT 10 Javadoc|VisibilityModifier
    private static class MergeState {
        final boolean allowMerge;
        final boolean sameSource;
        final boolean sameTarget;

        MergeState(final boolean allowMerge, final boolean sameSource, final boolean sameTarget) {
            this.allowMerge = allowMerge;
            this.sameSource = sameSource;
            this.sameTarget = sameTarget;
        }
    }
    
    /**
     * Merges the merge source node into the merge target node. All edges that were previously
     * connected to the merge source's ports are rerouted to the merge target. The merge target's
     * long edge source and target ports can be set to {@code null}.
     * 
     * @param mergeSource the merge source node.
     * @param mergeTarget the merge target node.
     * @param keepSourcePort if {@code false}, the long edge source property is set to {@code null}.
     *                       This should be done if the long edge sources of the two nodes to be merged
     *                       are different.
     * @param keepTargetPort if {@code false}, the long edge target property is set to {@code null}.
     *                       This should be done if the long edge targets of the two nodes to be merged
     *                       are different.
     */
    private void mergeNodes(final LNode mergeSource, final LNode mergeTarget,
            final boolean keepSourcePort, final boolean keepTargetPort) {
        
        // We assume that the input port is west, and the output port east
        LPort mergeTargetInputPort = mergeTarget.getPorts(PortSide.WEST).iterator().next();
        LPort mergeTargetOutputPort = mergeTarget.getPorts(PortSide.EAST).iterator().next();
        
        for (LPort port : mergeSource.getPorts()) {
            while (!port.getIncomingEdges().isEmpty()) {
                port.getIncomingEdges().get(0).setTarget(mergeTargetInputPort);
            }
            
            while (!port.getOutgoingEdges().isEmpty()) {
                port.getOutgoingEdges().get(0).setSource(mergeTargetOutputPort);
            }
        }
        
        // Possibly reset source and target ports
        if (!keepSourcePort) {
            mergeTarget.setProperty(InternalProperties.LONG_EDGE_SOURCE, null);
        }
        
        if (!keepTargetPort) {
            mergeTarget.setProperty(InternalProperties.LONG_EDGE_TARGET, null);
        }
    }

    private void identifyHyperedges(final LGraph lGraph) {
        List<LPort> ports = lGraph.getLayers().stream()
            .flatMap(l -> l.getNodes().stream())
            .flatMap(n -> n.getPorts().stream())
            .collect(Collectors.toList());
            
        ports.forEach(p -> p.id = -1);
        
        int index = 0;
        for (LPort p : ports) {
            if (p.id == -1) {
                dfs(p, index++);
            }
        }
    }

    private void dfs(final LPort p, final int index) {
        p.id = index;
        // follow edges connected to the same port
        for (LPort p2 : p.getConnectedPorts()) {
            if (p2.id == -1) {
                dfs(p2, index);
            }
        }
        // follow edges connected to the same long edge dummy
        if (p.getNode().getType() == NodeType.LONG_EDGE) {
            for (LPort p2 : p.getNode().getPorts()) {
                if (p2 != p && p2.id == -1) {
                    dfs(p2, index);
                }
            }
        }
    }
    
}
