/*******************************************************************************
 * Copyright (c) 2009, 2015 Kiel University and others.
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
package org.eclipse.elk.conn.gmf.layouter;

import java.util.HashMap;
import java.util.Map;

import org.eclipse.draw2d.PositionConstants;
import org.eclipse.draw2d.geometry.Dimension;
import org.eclipse.draw2d.geometry.Insets;
import org.eclipse.draw2d.geometry.Point;
import org.eclipse.draw2d.geometry.PointList;
import org.eclipse.draw2d.graph.CompoundDirectedGraphLayout;
import org.eclipse.draw2d.graph.DirectedGraph;
import org.eclipse.draw2d.graph.DirectedGraphLayout;
import org.eclipse.draw2d.graph.Edge;
import org.eclipse.draw2d.graph.Node;
import org.eclipse.elk.core.AbstractLayoutProvider;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.options.Direction;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkEdgeSection;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;

/**
 * Layout provider that uses the layout algorithm shipped with Draw2D. Either the
 * default version or the Compound version of the algorithm can be applied.
 * 
 * TODO implement compound graph layout using CompoundDirectedGraphLayout
 *
 * @author msp
 */
public class Draw2DLayoutProvider extends AbstractLayoutProvider {

    /** parameter value for the compound version of the layout algorithm. */
    public static final String PARAM_COMPOUND = "Compound";
    /** default value for minimal spacing. */
    private static final double DEF_MIN_SPACING = 16.0;
    
    /** indicates whether the compound graph version of the algorithm shall be used. */
    private boolean compoundMode = false;
    
    @Override
    public void initialize(final String parameter) {
        compoundMode = PARAM_COMPOUND.equals(parameter);
    }
    
    @Override
    public void layout(final ElkNode layoutNode, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Draw2D Directed Graph layout", 1);
        
        DirectedGraph graph = buildDraw2dGraph(layoutNode);
        DirectedGraphLayout draw2dLayout 
                = compoundMode ? new CompoundDirectedGraphLayout() : new DirectedGraphLayout();
        draw2dLayout.visit(graph);
        applyLayout(layoutNode, graph);
        
        progressMonitor.done();
    }
    
    /**
     * Builds the graph which is used internally by the Draw2D layout algorithm.
     * 
     * @param layoutNode parent layout node
     * @return a graph that contains the children of the given parent
     */
    @SuppressWarnings("unchecked")
    private DirectedGraph buildDraw2dGraph(final ElkNode layoutNode) {
        DirectedGraph graph = new DirectedGraph();
        
        // set layout options for the graph
        double minSpacing = layoutNode.getProperty(Draw2DOptions.SPACING_NODE_NODE);
        if (minSpacing < 0) {
            minSpacing = DEF_MIN_SPACING;
        }
        graph.setDefaultPadding(new Insets((int) minSpacing));
        
        ElkPadding padding = layoutNode.getProperty(Draw2DOptions.PADDING);
        graph.setMargin(new Insets((int) padding.top, (int) padding.left, (int) padding.bottom, (int) padding.right));
        Direction layoutDirection = layoutNode.getProperty(Draw2DOptions.DIRECTION);
        switch (layoutDirection) {
        case UP:
        case DOWN:
            graph.setDirection(PositionConstants.SOUTH);
            break;
        default:
            graph.setDirection(PositionConstants.EAST);
            break;
        }
        
        // add nodes to the graph
        Map<ElkNode, Node> nodeMap = new HashMap<>();
        for (ElkNode elknode : layoutNode.getChildren()) {
            Node draw2dNode = new Node(elknode);
            ElkUtil.resizeNode(elknode);
            draw2dNode.width = (int) elknode.getWidth();
            draw2dNode.height = (int) elknode.getHeight();
            nodeMap.put(elknode, draw2dNode);
            graph.nodes.add(draw2dNode);
        }
        
        // add edges to the graph
        for (ElkNode source : layoutNode.getChildren()) {
            Node draw2dSource = nodeMap.get(source);
            for (ElkEdge kedge : ElkGraphUtil.allOutgoingEdges(source)) {
                // we don't support hyperedges
                if (kedge.isHyperedge()) {
                    continue;
                }
                
                ElkNode target = ElkGraphUtil.connectableShapeToNode(kedge.getTargets().get(0));
                Node draw2dTarget = nodeMap.get(target);
                if (draw2dTarget != null && draw2dTarget != draw2dSource) {
                    Edge draw2dEdge = new Edge(kedge, draw2dSource, draw2dTarget);
                    graph.edges.add(draw2dEdge);
                }
            }
        }
        
        return graph;
    }
    
    /**
     * Applies the layout determined by Draw2D to the original graph.
     * 
     * @param parentNode the parent layout node
     * @param graph the Draw2D graph
     */
    private void applyLayout(final ElkNode parentNode, final DirectedGraph graph) {
        // apply node layouts
        // MIGRATE The original code does not apply parent node insets to the child node coordinates, but it should?
        for (int i = 0; i < graph.nodes.size(); i++) {
            Node node = graph.nodes.getNode(i);
            if (node.data instanceof ElkNode) {
                ElkNode knode = (ElkNode) node.data;
                knode.setLocation(node.x, node.y);
            }
        }
        
        // apply edge layouts
        for (int i = 0; i < graph.edges.size(); i++) {
            Edge edge = graph.edges.getEdge(i);
            if (edge.data instanceof ElkEdge) {
                ElkEdge kedge = (ElkEdge) edge.data;
                ElkEdgeSection edgeSection = ElkGraphUtil.firstEdgeSection(kedge, true, true);
                PointList pointList = edge.getPoints();
                
                // Start point
                Point startPoint = pointList.getFirstPoint();
                edgeSection.setStartLocation(startPoint.x, startPoint.y);
                
                // Bend points
                for (int j = 1; j < pointList.size() - 1; j++) {
                    Point point = pointList.getPoint(j);
                    ElkGraphUtil.createBendPoint(edgeSection, point.x, point.y);
                }
                
                // End point
                Point endPoint = pointList.getFirstPoint();
                edgeSection.setStartLocation(endPoint.x, endPoint.y);
            }
        }
        
        // apply parent node layout
        ElkPadding insets = parentNode.getProperty(Draw2DOptions.PADDING);
        Dimension layoutSize = graph.getLayoutSize();
        double width = insets.getLeft() + layoutSize.width + insets.getRight();
        double height = insets.getTop() + layoutSize.height + insets.getBottom();
        ElkUtil.resizeNode(parentNode, width, height, false, true);
    }

}
