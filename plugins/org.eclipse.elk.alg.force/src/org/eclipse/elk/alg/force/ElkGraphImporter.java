/*******************************************************************************
 * Copyright (c) 2011, 2015 Kiel University and others.
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
package org.eclipse.elk.alg.force;

import java.util.HashMap;
import java.util.Map;

import org.eclipse.elk.alg.force.graph.FBendpoint;
import org.eclipse.elk.alg.force.graph.FEdge;
import org.eclipse.elk.alg.force.graph.FGraph;
import org.eclipse.elk.alg.force.graph.FLabel;
import org.eclipse.elk.alg.force.graph.FNode;
import org.eclipse.elk.alg.force.options.ForceOptions;
import org.eclipse.elk.alg.force.options.InternalProperties;
import org.eclipse.elk.core.UnsupportedGraphException;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkEdgeSection;
import org.eclipse.elk.graph.ElkLabel;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;

/**
 * Manages the transformation of KGraphs to FGraphs.
 *
 * @author msp
 */
public class ElkGraphImporter implements IGraphImporter<ElkNode> {
    
    ///////////////////////////////////////////////////////////////////////////////
    // Transformation KGraph -> FGraph

    @Override
    public FGraph importGraph(final ElkNode kgraph) {
        FGraph fgraph = new FGraph();
        
        // copy the properties of the KGraph to the force graph
        fgraph.copyProperties(kgraph);
        // TODO Find another way to do this kind of bounds checking
//        fgraph.checkProperties(Properties.SPACING, Properties.ASPECT_RATIO, Properties.TEMPERATURE,
//                Properties.ITERATIONS, Properties.REPULSION);
        fgraph.setProperty(InternalProperties.ORIGIN, kgraph);
                
        // keep a list of created nodes in the force graph
        Map<ElkNode, FNode> elemMap = new HashMap<ElkNode, FNode>();
        
        // transform everything
        transformNodes(kgraph, fgraph, elemMap);
        transformEdges(kgraph, fgraph, elemMap);
                
        return fgraph;
    }
    
    /**
     * Transforms the nodes and ports defined by the given layout node.
     * 
     * @param parentNode the layout node whose edges to transform.
     * @param fgraph the force graph.
     * @param elemMap the element map that maps the original {@code KGraph} elements to the
     *                transformed {@code FGraph} elements.
     */
    private void transformNodes(final ElkNode parentNode, final FGraph fgraph, final Map<ElkNode, FNode> elemMap) {
        int index = 0;
        for (ElkNode knode : parentNode.getChildren()) {
            String label = "";
            if (!knode.getLabels().isEmpty()) {
                label = knode.getLabels().get(0).getText();
            }
            FNode newNode = new FNode(label);
            newNode.copyProperties(knode);
            newNode.setProperty(InternalProperties.ORIGIN, knode);
            
            newNode.id = index++;
            newNode.getPosition().x = knode.getX() + knode.getWidth() / 2;
            newNode.getPosition().y = knode.getY() + knode.getHeight() / 2;
            newNode.getSize().x = Math.max(knode.getWidth(), 1);
            newNode.getSize().y = Math.max(knode.getHeight(), 1);
            
            fgraph.getNodes().add(newNode);
            
            elemMap.put(knode, newNode);
            
            // port constraints cannot be undefined
            PortConstraints portConstraints = knode.getProperty(ForceOptions.PORT_CONSTRAINTS);
            if (portConstraints == PortConstraints.UNDEFINED) {
                portConstraints = PortConstraints.FREE;
            }
            
            // TODO consider ports
        }
    }
    
    /**
     * Transforms the edges defined by the given layout node.
     * 
     * @param parentNode the layout node whose edges to transform.
     * @param fgraph the force graph.
     * @param elemMap the element map that maps the original {@code KGraph} elements to the
     *                transformed {@code FGraph} elements.
     */
    private void transformEdges(final ElkNode parentNode, final FGraph fgraph, final Map<ElkNode, FNode> elemMap) {
        for (ElkNode knode : parentNode.getChildren()) {
            for (ElkEdge kedge : ElkGraphUtil.allOutgoingEdges(knode)) {
                // We don't support hyperedges
                if (kedge.isHyperedge()) {
                    throw new UnsupportedGraphException("Graph must not contain hyperedges.");
                }
                
                // exclude edges that pass hierarchy bounds as well as self-loops
                if (!kedge.isHierarchical() 
                        && knode != ElkGraphUtil.connectableShapeToNode(kedge.getTargets().get(0))) {
                    // create a force edge
                    FEdge newEdge = new FEdge();
                    newEdge.copyProperties(kedge);
                    // TODO
//                    newEdge.checkProperties(Properties.LABEL_SPACING, Properties.REPULSIVE_POWER);
                    newEdge.setProperty(InternalProperties.ORIGIN, kedge);
                    
                    newEdge.setSource(elemMap.get(knode));
                    newEdge.setTarget(elemMap.get(ElkGraphUtil.connectableShapeToNode(kedge.getTargets().get(0))));
                    
                    fgraph.getEdges().add(newEdge);
                    
                    // transform the edge's labels
                    for (ElkLabel klabel : kedge.getLabels()) {
                        FLabel newLabel = new FLabel(newEdge, klabel.getText());
                        newLabel.copyProperties(klabel);
                        newLabel.setProperty(InternalProperties.ORIGIN, klabel);
                        
                        newLabel.getSize().x = Math.max(klabel.getWidth(), 1);
                        newLabel.getSize().y = Math.max(klabel.getHeight(), 1);
                        newLabel.refreshPosition();
                        
                        fgraph.getLabels().add(newLabel);
                    }
                }
            }
        }
    }    
    
    ///////////////////////////////////////////////////////////////////////////////
    // Apply Layout Results

    @Override
    public void applyLayout(final FGraph fgraph) {
        ElkNode kgraph = (ElkNode) fgraph.getProperty(InternalProperties.ORIGIN);
        
        // calculate the offset from border spacing and node distribution
        double minXPos = Integer.MAX_VALUE;
        double minYPos = Integer.MAX_VALUE;
        double maxXPos = Integer.MIN_VALUE;
        double maxYPos = Integer.MIN_VALUE;
        
        for (FNode node : fgraph.getNodes()) {
            KVector pos = node.getPosition();
            KVector size = node.getSize();
            minXPos = Math.min(minXPos, pos.x - size.x / 2);
            minYPos = Math.min(minYPos, pos.y - size.y / 2);
            maxXPos = Math.max(maxXPos, pos.x + size.x / 2);
            maxYPos = Math.max(maxYPos, pos.y + size.y / 2);
        }
        for (FBendpoint bendpoint : fgraph.getBendpoints()) {
            KVector pos = bendpoint.getPosition();
            KVector size = bendpoint.getSize();
            minXPos = Math.min(minXPos, pos.x - size.x / 2);
            minYPos = Math.min(minYPos, pos.y - size.y / 2);
            maxXPos = Math.max(maxXPos, pos.x + size.x / 2);
            maxYPos = Math.max(maxYPos, pos.y + size.y / 2);
        }
        
        ElkPadding padding = kgraph.getProperty(ForceOptions.PADDING);
        KVector offset = new KVector(padding.getLeft() - minXPos, padding.getTop() - minYPos);

        // process the nodes
        for (FNode fnode : fgraph.getNodes()) {
            Object object = fnode.getProperty(InternalProperties.ORIGIN);
            
            if (object instanceof ElkNode) {
                // set the node position
                ElkNode knode = (ElkNode) object;
                KVector nodePos = new KVector(fnode.getPosition()).add(offset);
                knode.setLocation(nodePos.x - knode.getWidth() / 2, nodePos.y - knode.getHeight() / 2);
            }
        }
        
        // process the edges
        for (FEdge fedge : fgraph.getEdges()) {
            ElkEdge kedge = (ElkEdge) fedge.getProperty(InternalProperties.ORIGIN);
            ElkEdgeSection kedgeSection = ElkGraphUtil.firstEdgeSection(kedge, true, true);
            
            // Add bendpoints
            KVector startLocation = new KVector(fedge.getSourcePoint());
            startLocation.add(offset);
            kedgeSection.setStartLocation(startLocation.x, startLocation.y);
            fedge.getBendpoints().forEach(bp -> {
                KVector position = new KVector(bp.getPosition());
                position.add(offset);
                ElkGraphUtil.createBendPoint(kedgeSection, position.x, position.y);
            });
            
            KVector endLocation = new KVector(fedge.getTargetPoint());
            endLocation.add(offset);
            
            kedgeSection.setEndLocation(endLocation.x, endLocation.y);
        }
        
        // process the labels
        for (FLabel flabel : fgraph.getLabels()) {
            ElkLabel klabel = (ElkLabel) flabel.getProperty(InternalProperties.ORIGIN);
            KVector labelPos = new KVector(flabel.getPosition()).add(offset);
            klabel.setLocation(labelPos.x, labelPos.y);
        }
        
        // set up the parent node
        double width = (maxXPos - minXPos) + padding.getHorizontal();
        double height = (maxYPos - minYPos) + padding.getVertical();
        if (!kgraph.getProperty(CoreOptions.NODE_SIZE_FIXED_GRAPH_SIZE)) {
            ElkUtil.resizeNode(kgraph, width, height, false, true);
        }
        kgraph.setProperty(CoreOptions.CHILD_AREA_WIDTH, width - padding.getHorizontal());
        kgraph.setProperty(CoreOptions.CHILD_AREA_HEIGHT, height - padding.getVertical());
    }
    
}
