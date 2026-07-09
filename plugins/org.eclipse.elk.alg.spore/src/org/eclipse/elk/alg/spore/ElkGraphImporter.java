/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
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
package org.eclipse.elk.alg.spore;

import java.util.Map;

import org.eclipse.elk.alg.common.ICostFunction;
import org.eclipse.elk.alg.common.nodespacing.NodeDimensionCalculation;
import org.eclipse.elk.alg.common.spore.InternalProperties;
import org.eclipse.elk.alg.common.spore.Node;
import org.eclipse.elk.alg.common.utils.Utils;
import org.eclipse.elk.alg.spore.graph.Graph;
import org.eclipse.elk.alg.spore.options.CompactionStrategy;
import org.eclipse.elk.alg.spore.options.RootSelection;
import org.eclipse.elk.alg.spore.options.SporeCompactionOptions;
import org.eclipse.elk.alg.spore.options.SpanningTreeCostFunction;
import org.eclipse.elk.alg.spore.options.TreeConstructionStrategy;
import org.eclipse.elk.core.math.ElkMargin;
import org.eclipse.elk.core.math.ElkMath;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.math.ElkRectangle;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.core.util.Pair;
import org.eclipse.elk.core.util.adapters.ElkGraphAdapters;
import org.eclipse.elk.core.util.adapters.ElkGraphAdapters.ElkGraphAdapter;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkEdgeSection;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;

import com.google.common.collect.Maps;

/**
 * Provides the methods for conversions between ElkGraphs and Graphs.
 */
public class ElkGraphImporter implements IGraphImporter<ElkNode> {
    private Graph graph;
    /** Associates the vertices used in the SPOrE process with the original ElkNodes. */
    private Map<KVector, Pair<Node, ElkNode>> nodeMap;
    /** The original graph to be resized after the execution phase. */
    private ElkNode elkGraph;
    /** Spacing between nodes. */
    private double spacingNodeNode;

    // callback functions ----------------------------------------------------------------------------------------------
    /** A cost function for the spanning tree calculation, that represents the euclidean distance between 
     * the center points of two elements. */
    private final ICostFunction centerDistance = e -> {
        return e.u.distance(e.v);
    };
    
    /** A cost function for the spanning tree calculation, that represents the minimum of the euclidean distances 
     * between the center point of each element and the center of the root element. */
    private final ICostFunction minimumRootDistance = e -> {
        return Math.min(e.u.distance(graph.preferredRoot.vertex),
                e.v.distance(graph.preferredRoot.vertex));
    };
    
    /** A cost function for the spanning tree calculation, that represents the euclidean distance between 
     * the circumcircles of two elements. */
    private final ICostFunction circleUnderlap = e -> {
        Node n1 = nodeMap.get(e.u).getFirst();
        Node n2 = nodeMap.get(e.v).getFirst();
        return e.u.distance(e.v) - e.u.distance(n1.rect.getPosition())
                - e.v.distance(n2.rect.getPosition());
    };
    
    /** A cost function for the spanning tree calculation, that represents the distance in a certain direction between 
     * the bounding boxes of two elements. The direction is given by a vector connecting the center points. */
    private final ICostFunction rectangleUnderlap = e -> {
        Node n1 = nodeMap.get(e.u).getFirst();
        Node n2 = nodeMap.get(e.v).getFirst();
        return n1.underlap(n2);
    };
    
    /** Cost function for overlap removal by Nachmanson et al. "Node Overlap Removal by Growing a Tree".*/
    private final ICostFunction invertedOverlap = e -> {
        Node n1 = nodeMap.get(e.u).getFirst();
        Node n2 = nodeMap.get(e.v).getFirst();
        ElkRectangle r1 = n1.rect;
        ElkRectangle r2 = n2.rect;
        double dist = ElkMath.shortestDistance(r1, r2);
        if (dist >= 0) {
            return dist;
        }
        
        double s = r2.getCenter().sub(r1.getCenter()).length();
        return -(Utils.overlap(r1, r2) - 1) * s;
    };
    
    @Override
    public Graph importGraph(final ElkNode inputGraph) {
        
        elkGraph = inputGraph;
        
        nodeMap = Maps.newHashMap();
        
        // calculate margins
        ElkGraphAdapter adapter = ElkGraphAdapters.adapt(elkGraph);
        NodeDimensionCalculation.calculateNodeMargins(adapter);
        
        // retrieve layout options
        String preferredRootID = elkGraph.getProperty(SporeCompactionOptions.PROCESSING_ORDER_PREFERRED_ROOT);
        SpanningTreeCostFunction costFunctionID = 
                elkGraph.getProperty(SporeCompactionOptions.PROCESSING_ORDER_SPANNING_TREE_COST_FUNCTION);
        TreeConstructionStrategy treeConstructionStrategy = 
                elkGraph.getProperty(SporeCompactionOptions.PROCESSING_ORDER_TREE_CONSTRUCTION);
        CompactionStrategy compactionStrategy =
                elkGraph.getProperty(SporeCompactionOptions.COMPACTION_COMPACTION_STRATEGY);
        RootSelection rootSelection = elkGraph.getProperty(SporeCompactionOptions.PROCESSING_ORDER_ROOT_SELECTION);
        spacingNodeNode = elkGraph.getProperty(SporeCompactionOptions.SPACING_NODE_NODE);
        
        ICostFunction costFunction = centerDistance;
        switch (costFunctionID) {
        case CENTER_DISTANCE:
            costFunction = centerDistance;
            break;
        case CIRCLE_UNDERLAP:
            costFunction = circleUnderlap;
            break;
        case RECTANGLE_UNDERLAP:
            costFunction = rectangleUnderlap;
            break;
        case INVERTED_OVERLAP:
            costFunction = invertedOverlap;
            break;
        case MINIMUM_ROOT_DISTANCE:
            costFunction = minimumRootDistance;
            break;
        default:
            throw new IllegalArgumentException(
                    "No implementation available for " + costFunctionID.toString());
        }
        
        // instantiate Graph
        graph = new Graph(costFunction, treeConstructionStrategy, compactionStrategy);
        graph.setProperty(InternalProperties.DEBUG_SVG, elkGraph.getProperty(SporeCompactionOptions.DEBUG_MODE));
        graph.orthogonalCompaction = elkGraph.getProperty(SporeCompactionOptions.COMPACTION_ORTHOGONAL);
        
        if (elkGraph.getChildren().isEmpty()) {
            // don't bother
            return graph;
        }
        
        // create Nodes representing the ElkNodes
        for (ElkNode elkNode : elkGraph.getChildren()) {
            double halfWidth = elkNode.getWidth() / 2;
            double halfHeight = elkNode.getHeight() / 2;
            KVector vertex = new KVector(elkNode.getX() + halfWidth, 
                                         elkNode.getY() + halfHeight);
            
            // randomly shift identical points a tiny bit to make them unique
            while (nodeMap.containsKey(vertex)) {
                // SUPPRESS CHECKSTYLE NEXT 1 MagicNumber
                vertex.add((Math.random() - 0.5) * 0.001, (Math.random() - 0.5) * 0.001);
                // If two positions were identical, their corresponding edge in the spanning tree would be
                // of zero length, had no direction, and couldn't be scaled by anything.
            }
            
            ElkMargin margin = elkNode.getProperty(CoreOptions.MARGINS);
            
            Node node = new Node(vertex, new ElkRectangle(vertex.x - halfWidth - spacingNodeNode / 2 - margin.left, 
                    vertex.y - halfHeight - spacingNodeNode / 2 - margin.top, 
                    elkNode.getWidth() + spacingNodeNode + margin.getHorizontal(), 
                    elkNode.getHeight() + spacingNodeNode + margin.getVertical()));
            
            graph.vertices.add(node);
            nodeMap.put(vertex, Pair.of(node, elkNode));
        }
        
        // spanning tree root selection method
        switch (rootSelection) {
        case FIXED:
            if (preferredRootID == null) {
                // get first Node in list if no ID specified
                graph.preferredRoot = graph.vertices.get(0);
            } else {
                // find Node associated with the ElkNode containing the ID
                for (Node node : graph.vertices) {
                    String id = nodeMap.get(node.originalVertex).getSecond().getIdentifier();
                    if (id != null && id.equals(preferredRootID)) {
                        graph.preferredRoot = node;
                    }
                }
            }
            break;

        case CENTER_NODE:
            // find node that is most central in the drawing
            KVector center = new KVector(elkGraph.getWidth(), elkGraph.getHeight());
            // CHECKSTYLEOFF MagicNumber
            center.scale(0.5);  
            center.add(elkGraph.getX(), elkGraph.getY());
            double closest = Double.POSITIVE_INFINITY;
            for (Node node : graph.vertices) {
                double distance = node.originalVertex.distance(center);
                if (distance < closest) {
                    closest = distance;
                    graph.preferredRoot = node;
                }
            }
            break;
            
        default:
            throw new IllegalArgumentException(
                    "No implementation available for " + rootSelection.toString());
        }
        
        return graph;
    }
    
    @Override
    public void updateGraph(final Graph g) {
        Map<KVector, Pair<Node, ElkNode>> updatedNodeMap = Maps.newHashMap();
        // reset graph
        g.tEdges = null;
        g.tree = null;
        
        // update nodes
        for (Node n : g.vertices) {
            Pair<Node, ElkNode> original = nodeMap.get(n.originalVertex);
            n.originalVertex = n.rect.getCenter();
            updatedNodeMap.put(n.originalVertex, original);
        }
        nodeMap = updatedNodeMap;
    }

    @Override
    public void applyPositions(final Graph g) {
        // set new node positions
        double minX = Double.POSITIVE_INFINITY;
        double minY = Double.POSITIVE_INFINITY;
        double maxX = Double.NEGATIVE_INFINITY;
        double maxY = Double.NEGATIVE_INFINITY;
        for (Node node : g.vertices) {
            ElkNode elkNode = nodeMap.get(node.originalVertex).getSecond();
            elkNode.setLocation(node.rect.x, node.rect.y);
            minX = Math.min(minX, elkNode.getX());
            minY = Math.min(minY, elkNode.getY());
            maxX = Math.max(maxX, elkNode.getX() + elkNode.getWidth());
            maxY = Math.max(maxY, elkNode.getY() + elkNode.getHeight());
        }
        
        // set new dimensions of parent node
        ElkPadding padding = elkGraph.getProperty(SporeCompactionOptions.PADDING);
        ElkUtil.resizeNode(elkGraph, maxX - minX + padding.getHorizontal(), 
                maxY - minY + padding.getVertical(), true, true);
        ElkUtil.translate(elkGraph, -minX + padding.left, -minY + padding.top);
        
        // update edges and route them as straight lines
        for (ElkEdge e : elkGraph.getContainedEdges()) {
            ElkEdgeSection kedgeSection = ElkGraphUtil.firstEdgeSection(e, true, true);
            
            ElkNode source = ElkGraphUtil.getSourceNode(e);
            ElkNode target = ElkGraphUtil.getTargetNode(e);
            KVector startLocation = new KVector(source.getX() + source.getWidth() / 2,
                    source.getY() + source.getHeight() / 2);
            KVector endLocation = new KVector(target.getX() + target.getWidth() / 2,
                    target.getY() + target.getHeight() / 2);
            
            KVector uv = endLocation.clone().sub(startLocation);
            ElkMath.clipVector(uv, source.getWidth(), source.getHeight());
            startLocation.add(uv);
            KVector vu = startLocation.clone().sub(endLocation);
            ElkMath.clipVector(vu, target.getWidth(), target.getHeight());
            endLocation.add(vu);
            
            kedgeSection.setStartLocation(startLocation.x, startLocation.y);
            kedgeSection.setEndLocation(endLocation.x, endLocation.y);
        }
    }
}
