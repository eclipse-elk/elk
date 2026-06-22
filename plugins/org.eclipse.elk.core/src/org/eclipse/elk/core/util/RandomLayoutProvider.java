/*******************************************************************************
 * Copyright (c) 2010, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.core.util;

import java.util.ListIterator;
import java.util.Random;

import org.eclipse.elk.core.AbstractLayoutProvider;
import org.eclipse.elk.core.math.ElkPadding;
import org.eclipse.elk.core.options.RandomLayouterOptions;
import org.eclipse.elk.graph.ElkBendPoint;
import org.eclipse.elk.graph.ElkConnectableShape;
import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkEdgeSection;
import org.eclipse.elk.graph.ElkGraphFactory;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.ElkPort;
import org.eclipse.elk.graph.util.ElkGraphUtil;

/**
 * Layout provider that computes random layouts. Can be useful to demonstrate the difference
 * between a good layout and an extremely bad one.
 * 
 * <p>
 * MIGRATE The random layout provider does not support hyperedges yet.
 * 
 * For that to work, the {@link #randomize(ElkEdge, Random, double, double)} method needs to be extended.
 * </p>
 *
 * @author msp
 */
public class RandomLayoutProvider extends AbstractLayoutProvider {
    
    @Override
    public void layout(final ElkNode parentNode, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Random Layout", 1);
        
        if (parentNode.getChildren().isEmpty()) {
            progressMonitor.done();
            return;
        }
        
        // initialize random seed
        Random random;
        Integer randomSeed = parentNode.getProperty(RandomLayouterOptions.RANDOM_SEED);
        if (randomSeed != null && randomSeed != 0) {
            random = new Random(randomSeed);
        } else {
            random = new Random();
        }
        
        // retrieve some layout option values
        float aspectRatio = parentNode.getProperty(RandomLayouterOptions.ASPECT_RATIO).floatValue();
        float spacing = parentNode.getProperty(RandomLayouterOptions.SPACING_NODE_NODE).floatValue();
        ElkPadding padding = parentNode.getProperty(RandomLayouterOptions.PADDING);
        
        // randomize the layout
        randomize(parentNode, random, aspectRatio, spacing, padding);

        progressMonitor.done();
    }
    
    /**
     * Randomize the given graph.
     * 
     * @param parent the parent node of the graph
     * @param random the random number generator
     * @param aspectRatio desired aspect ratio (must not be 0)
     * @param spacing desired object spacing
     * @param offset offset to the border
     */
    private void randomize(final ElkNode parent, final Random random, final double aspectRatio,
            final double spacing, final ElkPadding padding) {
        
        // determine width and height of the drawing and count the number of edges
        double nodesArea = 0.0f, maxWidth = 0.0f, maxHeight = 0.0f;
        int m = 1;
        for (ElkNode node : parent.getChildren()) {
            for (ElkEdge ignored : ElkGraphUtil.allOutgoingEdges(node)) {
                m++;
            }
            
            double width = node.getWidth();
            maxWidth = Math.max(maxWidth, width);
            double height = node.getHeight();
            maxHeight = Math.max(maxHeight, height);
            
            nodesArea += width * height;
        }
        int n = parent.getChildren().size();
        
        // a heuristic formula that determines an area in which nodes are randomly distributed
        double drawArea = nodesArea + 2 * spacing * spacing * m * n;
        double areaSqrt = (double) Math.sqrt(drawArea);
        double drawWidth = Math.max(areaSqrt * aspectRatio, maxWidth);
        double drawHeight = Math.max(areaSqrt / aspectRatio, maxHeight);
        
        // randomize node positions
        for (ElkNode node : parent.getChildren()) {
            double x = padding.getLeft() + random.nextDouble() * (drawWidth - node.getWidth());
            double y = padding.getLeft() + random.nextDouble() * (drawHeight - node.getHeight());
            node.setLocation(x, y);
        }
        
        // randomize edge positions
        double totalWidth = drawWidth + padding.getHorizontal();
        double totalHeight = drawHeight + padding.getVertical();
        for (ElkNode source : parent.getChildren()) {
            for (ElkEdge edge : ElkGraphUtil.allOutgoingEdges(source)) {
                if (!edge.isHierarchical()) {
                    randomize(edge, random, totalWidth, totalHeight);
                }
            }
        }
        
        totalWidth += padding.getLeft() + padding.getRight();
        totalHeight += padding.getTop() + padding.getBottom();
        ElkUtil.resizeNode(parent, totalWidth, totalHeight, false, true);
    }
    
    /** the maximal number of generated bend points for each edge. */
    private static final int MAX_BENDS = 5;
    /** a factor for the distance between source and target node, determines how much edge
     *  bend point may deviate from the straight line between those nodes. */
    private static final double RAND_FACT = 0.2f;
    
    /**
     * Randomize the given edge by adding bend points in the area between the source and target node.
     * 
     * @param edge an edge
     * @param source the source node
     * @param target the target node
     * @param random the random number generator
     * @param drawWidth the total width of the drawing
     * @param drawHeight the total Height of the drawing
     */
    private void randomize(final ElkEdge edge, final Random random, final double drawWidth, final double drawHeight) {
        // Determine source coordinate
        ElkConnectableShape sourceShape = edge.getSources().get(0);
        double sourceX = sourceShape.getX();
        double sourceY = sourceShape.getY();
        double sourceWidth = sourceShape.getWidth() / 2;
        double sourceHeight = sourceShape.getHeight() / 2;
        
        if (sourceShape instanceof ElkPort) {
            ElkPort sourcePort = (ElkPort) sourceShape;
            sourceX += sourcePort.getParent().getX();
            sourceX += sourcePort.getParent().getX();
        }
        
        sourceX += sourceWidth;
        sourceY += sourceHeight;
        
        // Determine target coordinate
        ElkConnectableShape targetShape = edge.getSources().get(0);
        double targetX = targetShape.getX();
        double targetY = targetShape.getY();
        double targetWidth = targetShape.getWidth() / 2;
        double targetHeight = targetShape.getHeight() / 2;
        
        if (targetShape instanceof ElkPort) {
            ElkPort targetPort = (ElkPort) targetShape;
            targetX += targetPort.getParent().getX();
            targetX += targetPort.getParent().getX();
        }
        
        targetX += targetWidth;
        targetY += targetHeight;
        
        // Ensure that the edge has a single edge section
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
        
        ElkEdgeSection edgeSection = edge.getSections().get(0);
        
        // set the source point onto the border of the source element
        double sourcePX = targetX;
        if (targetX > sourceX + sourceWidth) {
            sourcePX = sourceX + sourceWidth;
        } else if (targetX < sourceX - sourceWidth) {
            sourcePX = sourceX - sourceWidth;
        }
        
        double sourcePY = targetY;
        if (targetY > sourceY + sourceHeight) {
            sourcePY = sourceY + sourceHeight;
        } else if (targetY < sourceY - sourceHeight) {
            sourcePY = sourceY - sourceHeight;
        }
        
        if (sourcePX > sourceX - sourceWidth && sourcePX < sourceX + sourceWidth
                && sourcePY > sourceY - sourceHeight && sourcePY < sourceY + sourceHeight) {
            sourcePX = sourceX + sourceWidth;
        }
        
        edgeSection.setStartLocation(sourcePX, sourcePY);
        
        // set the target point onto the border of the target element
        double targetPX = sourceX;
        if (sourceX > targetX + targetWidth) {
            targetPX = targetX + targetWidth;
        } else if (sourceX < targetX - targetWidth) {
            targetPX = targetX - targetWidth;
        }
        
        double targetPY = sourceY;
        if (sourceY > targetY + targetHeight) {
            targetPY = targetY + targetHeight;
        } else if (sourceY < targetY - targetHeight) {
            targetPY = targetY - targetHeight;
        }
        
        if (targetPX > targetX - targetWidth && targetPX < targetX + targetWidth
                && targetPY > targetY - targetHeight && targetPY < targetY + targetHeight) {
            targetPY = targetY + targetHeight;
        }
        
        edgeSection.setEndLocation(targetPX, targetPY);
        
        // add a random number of bend points
        edgeSection.getBendPoints().clear();
        int bendsNum = random.nextInt(MAX_BENDS);
        if (sourceShape == targetShape) {
            bendsNum++;
        }
        double xdiff = targetPX - sourcePX;
        double ydiff = targetPY - sourcePY;
        double totalDist = (double) Math.sqrt(xdiff * xdiff + ydiff * ydiff);
        double maxRand = totalDist * RAND_FACT;
        double xincr = xdiff / (bendsNum + 1);
        double yincr = ydiff / (bendsNum + 1);
        double x = sourcePX, y = sourcePY;
        for (int i = 0; i < bendsNum; i++) {
            // determine coordinates that deviate from the straight connection by a random amount
            x += xincr;
            y += yincr;
            double randx = x + random.nextFloat() * maxRand - maxRand / 2;
            if (randx < 0) {
                randx = 1;
            } else if (randx > drawWidth) {
                randx = drawWidth - 1;
            }
            double randy = y + random.nextFloat() * maxRand - maxRand / 2;
            if (randy < 0) {
                randy = 1;
            } else if (randy > drawHeight) {
                randy = drawHeight - 1;
            }
            ElkBendPoint bendPoint = ElkGraphFactory.eINSTANCE.createElkBendPoint();
            bendPoint.set(randx, randy);
            edgeSection.getBendPoints().add(bendPoint);
        }
    }

}
