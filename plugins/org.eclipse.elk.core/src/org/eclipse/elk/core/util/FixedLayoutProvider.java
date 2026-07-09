/*******************************************************************************
 * Copyright (c) 2010, 2019 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.

 * This Source Code may also be made available under the following Secondary
 * Licenses when the conditions for such availability set forth in the Eclipse
 * Public License v. 2.0 are satisfied: GPL-3.0 which is available at
 * https://www.gnu.org/licenses/gpl-3.0-standalone.html.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.util;

import java.util.ListIterator;

import org.eclipse.elk.core.AbstractLayoutProvider;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.math.KVectorChain;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.options.EdgeRouting;
import org.eclipse.elk.core.options.FixedLayouterOptions;
import org.eclipse.elk.core.options.SizeConstraint;
import org.eclipse.elk.graph.ElkBendPoint;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkEdgeSection;
import org.eclipse.elk.graph.ElkGraphFactory;
import org.eclipse.elk.graph.ElkLabel;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.eclipse.elk.graph.util.ElkGraphUtil;

/**
 * A layout provider that sets fixed positions for all elements. These positions are taken
 * from the {@link CoreOptions#POSITION} and {@link CoreOptions#BEND_POINTS} options.
 * Elements that have no position option attached just stay where they are.
 * This is useful for at least two things:
 * <ul>
 *   <li>Fix the layout of a part of the diagram so it won't be affected by automatic layout.</li>
 *   <li>Apply a layout imported from somewhere else, e.g. the original layout that was
 *     manually created in another modeling tool.</li>
 * </ul>
 * 
 * <p>
 * MIGRATE The fixed layout provider does not support hyperedges yet.
 * </p>
 */
public class FixedLayoutProvider extends AbstractLayoutProvider {

    @Override
    public void layout(final ElkNode layoutNode, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Fixed Layout", 1);
        
        EdgeRouting edgeRouting = layoutNode.getProperty(CoreOptions.EDGE_ROUTING);
        
        double maxx = 0;
        double maxy = 0;
        
        for (ElkNode node : layoutNode.getChildren()) {
            // set the fixed position of the node, or leave it as it is
            KVector pos = node.getProperty(FixedLayouterOptions.POSITION);
            if (pos != null) {
                node.setLocation(pos.x, pos.y);
                
                // set the fixed size of the node
                if (node.getProperty(FixedLayouterOptions.NODE_SIZE_CONSTRAINTS).contains(
                        SizeConstraint.MINIMUM_SIZE)) {
                    
                    KVector minSize = node.getProperty(FixedLayouterOptions.NODE_SIZE_MINIMUM);
                    if (minSize.x > 0 && minSize.y > 0) {
                        ElkUtil.resizeNode(node, minSize.x, minSize.y, true, true);
                    }
                }
            }
            maxx = Math.max(maxx, node.getX() + node.getWidth());
            maxy = Math.max(maxy, node.getY() + node.getHeight());
            
            // set the fixed position of the node labels, or leave them as they are
            for (ElkLabel label : node.getLabels()) {
                pos = label.getProperty(FixedLayouterOptions.POSITION);
                if (pos != null) {
                    label.setLocation(pos.x, pos.y);
                }
                maxx = Math.max(maxx, node.getX() + label.getX() + label.getWidth());
                maxy = Math.max(maxy, node.getY() + label.getY() + label.getHeight());
            }
                
            // set the fixed position of the ports, or leave them as they are
            for (ElkPort port : node.getPorts()) {
                pos = port.getProperty(FixedLayouterOptions.POSITION);
                if (pos != null) {
                    port.setLocation(pos.x, pos.y);
                }
                
                double portx = node.getX() + port.getX();
                double porty = node.getY() + port.getY();
                maxx = Math.max(maxx, portx + port.getWidth());
                maxy = Math.max(maxy, porty + port.getHeight());
                
                // set the fixed position of the port labels, or leave them as they are
                for (ElkLabel label : port.getLabels()) {
                    pos = label.getProperty(FixedLayouterOptions.POSITION);
                    if (pos != null) {
                        label.setLocation(pos.x, pos.y);
                    }
                    maxx = Math.max(maxx, portx + label.getX() + label.getWidth());
                    maxy = Math.max(maxy, porty + label.getY() + label.getHeight());
                }
            }
            
            // set fixed routing for the node's outgoing edges (both simple and hierarchical), or leave them as they are
            for (ElkEdge edge : ElkGraphUtil.allOutgoingEdges(node)) {
                KVector maxv = processEdge(edge, edgeRouting);
                maxx = Math.max(maxx, maxv.x);
                maxy = Math.max(maxy, maxv.y);
            }
            
            // additionally handle incoming hierarchical edges (both short and long) 
            //  note that this potentially results in hierarchical edges being handled twice
            for (ElkEdge edge : ElkGraphUtil.allIncomingEdges(node)) {
                if (ElkGraphUtil.getSourceNode(edge).getParent() != layoutNode) {
                    KVector maxv = processEdge(edge, edgeRouting);
                    maxx = Math.max(maxx, maxv.x);
                    maxy = Math.max(maxy, maxv.y);
                }
            }
        }
        
        // if orthogonal routing is selected, determine the junction points
        if (edgeRouting == EdgeRouting.ORTHOGONAL) {
            for (ElkNode node : layoutNode.getChildren()) {
                for (ElkEdge edge : ElkGraphUtil.allOutgoingEdges(node)) {
                    generateJunctionPoints(edge);
                }  
            }
        }
        
        // set size of the parent node unless its size should be fixed as well
        if (!layoutNode.getProperty(FixedLayouterOptions.NODE_SIZE_FIXED_GRAPH_SIZE)) {
            ElkPadding padding = layoutNode.getProperty(FixedLayouterOptions.PADDING);
            double newWidth = maxx + padding.getLeft() + padding.getRight();
            double newHeight = maxy + padding.getTop() + padding.getBottom();
            ElkUtil.resizeNode(layoutNode, newWidth, newHeight, true, true);
        }
        progressMonitor.done();
    }
    
    /**
     * Process an edge and its labels.
     * 
     * @param edge an edge
     * @param edgeRouting the global edge routing setting
     */
    private KVector processEdge(final ElkEdge edge, final EdgeRouting edgeRouting) {
        // MIGRATE Does this properly support hyperedges?
        ElkNode sourceParent = ElkGraphUtil.connectableShapeToNode(edge.getSources().get(0)).getParent();
        ElkNode targetParent = ElkGraphUtil.connectableShapeToNode(edge.getTargets().get(0)).getParent();
        boolean sameHierarchy = sourceParent == targetParent;
        
        KVector maxv = new KVector();
        KVectorChain bendPoints = edge.getProperty(FixedLayouterOptions.BEND_POINTS);
        
        // we need at least two bend points, since the source point and target point must be included
        if (bendPoints != null && bendPoints.size() >= 2) {
            if (edge.getSections().isEmpty()) {
                // We need an edge section to apply the bend points to
                ElkEdgeSection edgeSection = ElkGraphFactory.eINSTANCE.createElkEdgeSection();
                edge.getSections().add(edgeSection);
            } else if (edge.getSections().size() > 1) {
                // We can only apply bend points to a single edge section, so throw away all except for the last one
                ListIterator<ElkEdgeSection> sections = edge.getSections().listIterator();
                while (sections.hasNext()) {
                    sections.remove();
                }
            }
            
            ElkUtil.applyVectorChain(bendPoints, edge.getSections().get(0));
        }
        
        // determine maximal coordinates
        if (sameHierarchy) {
            for (ElkEdgeSection edgeSection : edge.getSections()) {
                for (ElkBendPoint point : edgeSection.getBendPoints()) {
                    maxv.x = Math.max(maxv.x, point.getX());
                    maxv.y = Math.max(maxv.y, point.getY());
                }
            }
        }
        
        // set the fixed position of the edge labels, or leave them as they are
        for (ElkLabel label : edge.getLabels()) {
            KVector pos = label.getProperty(FixedLayouterOptions.POSITION);
            if (pos != null) {
                label.setLocation(pos.x, pos.y);
            }
            
            if (sameHierarchy) {
                maxv.x = Math.max(maxv.x, label.getX() + label.getWidth());
                maxv.y = Math.max(maxv.y, label.getY() + label.getHeight());
            }
        }
        
        return maxv;
    }

    /**
     * Determine if given edge has junction points and add appropriate layout option if necessary.
     * 
     * @param edge the edge to determine junction points of
     */
    private void generateJunctionPoints(final ElkEdge edge) {
        // Note: if the edge coordinates are not modified, the junction points are also ignored
        KVectorChain junctionPoints = ElkUtil.determineJunctionPoints(edge);
        if (junctionPoints.isEmpty()) {
            edge.setProperty(CoreOptions.JUNCTION_POINTS, null);
        } else {
            edge.setProperty(CoreOptions.JUNCTION_POINTS, junctionPoints);
        }
    }
    
}
