/*******************************************************************************
 * Copyright (c) 2013 - 2022 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.mrtree;

import java.util.HashMap;
import java.util.Map;

import org.eclipse.elk.alg.mrtree.graph.TEdge;
import org.eclipse.elk.alg.mrtree.graph.TGraph;
import org.eclipse.elk.alg.mrtree.graph.TNode;
import org.eclipse.elk.alg.mrtree.options.InternalProperties;
import org.eclipse.elk.alg.mrtree.options.MrTreeOptions;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.math.KVectorChain;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkEdgeSection;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;

/**
 * Manages the transformation of KGraphs to TGraphs.
 * 
 * @author sor
 * @author sgu
 * @author sdo
 */
public class ElkGraphImporter implements IGraphImporter<ElkNode> {

    // /////////////////////////////////////////////////////////////////////////////
    // Transformation KGraph -> TGraph

    /**
     * {@inheritDoc}
     */
    public TGraph importGraph(final ElkNode elkgraph) {
        TGraph tGraph = new TGraph();

        // copy the properties of the KGraph to the t-graph
        tGraph.copyProperties(elkgraph);
        // TODO Find a new concept for checking validity of bounds
//        tGraph.checkProperties(Properties.SPACING, Properties.ASPECT_RATIO);
        tGraph.setProperty(InternalProperties.ORIGIN, elkgraph);

        // keep a list of created nodes in the t-graph
        Map<ElkNode, TNode> elemMap = new HashMap<ElkNode, TNode>();

        // transform nodes and edges
        transformNodes(elkgraph, tGraph, elemMap);
        transformEdges(elkgraph, tGraph, elemMap);

        return tGraph;
    }

    /**
     * Transforms the nodes and ports defined by the given layout node.
     * 
     * @param parentNode
     *            the layout node whose edges to transform.
     * @param tGraph
     *            the t-graph.
     * @param elemMap
     *            the element map that maps the original {@code KGraph} elements to the transformed
     *            {@code TGraph} elements.
     */
    private void transformNodes(final ElkNode parentNode, final TGraph tGraph, final Map<ElkNode, TNode> elemMap) {
        int index = 0;

        for (ElkNode elknode : parentNode.getChildren()) {
            // copy label
            String label = "";
            if (!elknode.getLabels().isEmpty()) {
                label = elknode.getLabels().get(0).getText();
            }

            // create new tNode
            TNode newNode = new TNode(index++, tGraph, label);
            newNode.copyProperties(elknode);
            newNode.setProperty(InternalProperties.ORIGIN, elknode);
            
            newNode.getPosition().y = elknode.getY() + elknode.getHeight() / 2;
            newNode.getSize().x = Math.max(elknode.getWidth(), 1);
            newNode.getPosition().x = elknode.getX() + elknode.getWidth() / 2;
            newNode.getSize().y = Math.max(elknode.getHeight(), 1);
            
            tGraph.getNodes().add(newNode);

            // keep the corresponding tNode of each kNode for edge transformation
            elemMap.put(elknode, newNode);
        }
    }

    /**
     * Transforms the edges defined by the given layout node.
     * 
     * @param parentNode
     *            the layout node whose edges to transform.
     * @param tGraph
     *            the t-graph.
     * @param elemMap
     *            the element map that maps the original {@code KGraph} elements to the transformed
     *            {@code TGraph} elements.
     */
    private void transformEdges(final ElkNode parentNode, final TGraph tGraph, final Map<ElkNode, TNode> elemMap) {
        for (ElkNode elknode : parentNode.getChildren()) {
            for (ElkEdge elkedge : ElkGraphUtil.allOutgoingEdges(elknode)) {
                // exclude hyperedges, edges that pass hierarchy bounds, and self-loops
                if (!elkedge.isHierarchical() && !elkedge.isHierarchical() && !elkedge.isSelfloop()) {
                    // find the corresponding source and target tNode of edge
                    TNode source = elemMap.get(elknode);
                    TNode target = elemMap.get(ElkGraphUtil.connectableShapeToNode(elkedge.getTargets().get(0)));

                    if (source != null && target != null) {
                        // create a edge and add edge to tGraph
                        TEdge newEdge = new TEdge(source, target);
                        newEdge.setProperty(InternalProperties.ORIGIN, elkedge);

                        // TODO transform the edge's labels

                        // set properties of the new edge
                        newEdge.copyProperties(elkedge);

                        // update tNode accordingly
                        source.getOutgoingEdges().add(newEdge);
                        target.getIncomingEdges().add(newEdge);

                        tGraph.getEdges().add(newEdge);
                    }
                }
            }
        }
    }

    // /////////////////////////////////////////////////////////////////////////////
    // Apply Layout Results

    /**
     * {@inheritDoc}
     */
    public void applyLayout(final TGraph tGraph) {
        // get the corresponding kGraph
        ElkNode elkgraph = (ElkNode) tGraph.getProperty(InternalProperties.ORIGIN);
        
        // calculate the offset from border spacing and node distribution
        double minXPos = Integer.MAX_VALUE;
        double minYPos = Integer.MAX_VALUE;
        double maxXPos = Integer.MIN_VALUE;
        double maxYPos = Integer.MIN_VALUE;
        for (TNode tNode : tGraph.getNodes()) {
            KVector pos = tNode.getPosition();
            KVector size = tNode.getSize();
            minXPos = Math.min(minXPos, pos.x - size.x / 2);
            minYPos = Math.min(minYPos, pos.y - size.y / 2);
            maxXPos = Math.max(maxXPos, pos.x + size.x / 2);
            maxYPos = Math.max(maxYPos, pos.y + size.y / 2);
        }

        ElkPadding padding = elkgraph.getProperty(MrTreeOptions.PADDING);
        
        // apply tNode positions to elkNodes
        for (TNode tNode : tGraph.getNodes()) {
            Object object = tNode.getProperty(InternalProperties.ORIGIN);
            if (object instanceof ElkNode) {
                ElkNode elknode = (ElkNode) object;
                elknode.setLocation(tNode.getPosition().x, tNode.getPosition().y);
                elknode.copyProperties(tNode);
            }
        }
        
        // copy tEdge bendpoints to elkEdges
        for (TEdge tEdge : tGraph.getEdges()) {
            ElkEdge elkedge = (ElkEdge) tEdge.getProperty(InternalProperties.ORIGIN);
            if (elkedge != null) {
                KVectorChain bendPoints = tEdge.getBendPoints();
                ElkEdgeSection edgeSection = ElkGraphUtil.firstEdgeSection(elkedge, true, true);
                ElkUtil.applyVectorChain(bendPoints, edgeSection);
            }
        }

        // set up the graph
        double width = maxXPos - minXPos + padding.getHorizontal();
        double height = maxYPos - minYPos + padding.getVertical();
        if (!elkgraph.getProperty(CoreOptions.NODE_SIZE_FIXED_GRAPH_SIZE)) {
            ElkUtil.resizeNode(elkgraph, width, height, false, false);
        }
        elkgraph.setProperty(CoreOptions.CHILD_AREA_WIDTH, width - padding.getHorizontal());
        elkgraph.setProperty(CoreOptions.CHILD_AREA_HEIGHT, height - padding.getVertical());
        
    }


    /**
     * Modify the given center position to the border of the node.
     * 
     * @param center the node center position (modified by this method)
     * @param next the next point of the edge vector chain
     * @param size the node size
     */
    private static void toNodeBorder(final KVector center, final KVector next, final KVector size) {
        double wh = size.x / 2, hh = size.y / 2;
        double absx = Math.abs(next.x - center.x), absy = Math.abs(next.y - center.y);
        double xscale = 1, yscale = 1;
        if (absx > wh) {
            xscale = wh / absx;
        }
        if (absy > hh) {
            yscale = hh / absy;
        }
        double scale = Math.min(xscale, yscale);
        center.x += scale * (next.x - center.x);
        center.y += scale * (next.y - center.y);
    }
}
