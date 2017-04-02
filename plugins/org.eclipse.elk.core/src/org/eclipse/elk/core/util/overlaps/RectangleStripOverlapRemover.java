/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * which accompanies this distribution, and is available at
 * http://www.eclipse.org/legal/epl-v10.html
 *
 * Contributors:
 *    Kiel University - initial API and implementation
 *******************************************************************************/
package org.eclipse.elk.core.util.overlaps;

import java.util.List;
import java.util.SortedSet;

import org.eclipse.elk.core.math.ElkRectangle;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

/**
 * Main class for removing overlaps between rectangles that have a fixed position along one of the two dimsions. Yes,
 * only two dimensions. Sorry.
 * 
 * <p>To give an example of the typical use case, this thing might be used for labels of northern ports. Their x
 * position is fixed, which can result in overlaps. This is resolved by choosing y positions that remove overlaps.</p>
 * 
 * <p>This class is in charge of the main overlap removal configuration and of building the overlap graph, in which
 * each label is represented by a node and nodes are connected if their labels overlap. The actual overlap removal is
 * performed by overlap removal strategies.</p>
 */
public final class RectangleStripOverlapRemover {

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Constants
    
    /** The default gap to be left between labels. */
    private static final double DEFAULT_GAP = 5;
    

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Fields

    /** The direction in which to move rectangles to remove overlaps. */
    private OverlapRemovalDirection overlapRemovalDirection;
    /** Gap to be left between labels for them to not be considered overlapping. */
    private double gap = DEFAULT_GAP;
    /** The coordinate where the first rectangles, viewed along the overlap removal direction, are placed. */
    private double startCoordinate = 0;
    /** The overlap strategy to use to actually remove overlaps. */
    private IRectangleStripOverlapRemovalStrategy overlapRemovalStrategy;
    /** Rectangle nodes sorted by the x coordinate of their rectangle. */
    private SortedSet<RectangleNode> rectangleNodes = Sets.newTreeSet(
            RectangleStripOverlapRemover::compareLeftRectangleBorders);


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Creation
    
    /**
     * Create a new remover that removes overlaps between rectangles by moving them in the given direction.
     * 
     * @param direction
     *            the direction to move rectangles in.
     * @return the new remover.
     */
    public static RectangleStripOverlapRemover createForDirection(final OverlapRemovalDirection direction) {
        RectangleStripOverlapRemover remover = new RectangleStripOverlapRemover();
        remover.overlapRemovalDirection = direction;
        return remover;
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Configuration
    
    /**
     * Configures the gap to be left between rectangles. Returns this instance to enable method chaining.
     */
    public RectangleStripOverlapRemover withGap(final double theGap) {
        gap = theGap;
        return this;
    }
    
    /**
     * Configures the start coordinate. See the class documentation for more information on what the start coordinate
     * is. Returns this instance to enable method chaining.
     */
    public RectangleStripOverlapRemover withStartCoordinate(final double coordinate) {
        startCoordinate = coordinate;
        return this;
    }
    
    /**
     * Sets the overlap removal strategy to be used to actually remove overlaps. Returns this instance to enable method
     * chaining.
     */
    public RectangleStripOverlapRemover withOverlapRemovalStrategy(
            final IRectangleStripOverlapRemovalStrategy strategy) {
        
        overlapRemovalStrategy = strategy;
        return this;
    }
    
    /**
     * Adds the given rectangle to the list of rectangles to remove overlaps between. Returns this instance to enable
     * method chaining.
     */
    public RectangleStripOverlapRemover addRectangle(final ElkRectangle rectangle) {
        rectangleNodes.add(new RectangleNode(rectangle, importRectangle(rectangle)));
        return this;
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Getters
    
    /**
     * Returns the gap to be left between rectangles.
     */
    public double getGap() {
        return gap;
    }
    
    /**
     * Returns the set of rectangle nodes in the rectangle overlap graph, sorted by their rectangle's x coordinate.
     */
    public SortedSet<RectangleNode> getRectangleNodes() {
        return rectangleNodes;
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Coordinate Transformation
    
    /**
     * Returns a rectangle which represents the given rectangle when transformed into the coordinate space we will be
     * working in.
     */
    private ElkRectangle importRectangle(final ElkRectangle rectangle) {
        switch (overlapRemovalDirection) {
        case UP:
        case DOWN:
            return rectangle;
            
        case LEFT:
        case RIGHT:
            return new ElkRectangle(rectangle.y, 0, rectangle.height, rectangle.width);
            
        default:
            assert false;
            return null;
        }
    }
    
    /**
     * Applies the movement computed for the given rectangle node's rectangle to its original rectangle, paying
     * attention to the overlap removal direction.
     * 
     * @param rectangleNode
     *            the rectangle node whose results to apply to the original rectangle.
     * @param stripSize
     *            strip size computed by the overlap removal algorithm.
     */
    private void exportRectangle(final RectangleNode rectangleNode, final double stripSize) {
        ElkRectangle rectangle = rectangleNode.rectangle;
        ElkRectangle originalRectangle = rectangleNode.originalRectangle;
        
        switch (overlapRemovalDirection) {
        case UP:
            originalRectangle.y = startCoordinate - rectangle.height - rectangle.y;
            break;
            
        case DOWN:
            originalRectangle.y += startCoordinate;
            break;
            
        case LEFT:
            originalRectangle.x = startCoordinate - rectangle.height - rectangle.y;
            break;
            
        case RIGHT:
            originalRectangle.x = startCoordinate + rectangle.y;
            break;
        }
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Actual Algorithm
    
    /**
     * Triggers the overlap removal algorithm.
     * 
     * @return the height of the resulting label strip when removing horizontal overlaps, or the width of the resulting
     *         label strip when removing vertical overlaps.
     */
    public double removeOverlaps() {
        // Check if an overlap removal strategy is used; if not, fall back on the greedy strategy
        if (overlapRemovalStrategy == null) {
            overlapRemovalStrategy = new GreedyRectangleStripOverlapRemover();
        }
        
        // Compute and remove overlaps
        computeOverlaps();
        double stripSize = overlapRemovalStrategy.removeOverlaps(this);
        
        // Apply the results
        rectangleNodes.stream()
                .forEach(node -> exportRectangle(node, stripSize));
        
        return stripSize;
    }
    
    /**
     * Populates the {@link RectangleNode#overlappingNodes} lists. This method uses a scanline algorithm that remembers
     * which rectangles currently intersect the scanline. As it encounters a new rectangle, it first removes rectangles
     * that do not intersect the scanline anymore and then adds overlaps between all currently intersecting rectangles.
     */
    private void computeOverlaps() {
        SortedSet<RectangleNode> intersectingNodes = Sets.newTreeSet(
                RectangleStripOverlapRemover::compareRightRectangleBorders);
        double scanlinePos = Double.MIN_VALUE;
        
        // We iterate over the nodes according to their x coordinate
        for (RectangleNode currNode : rectangleNodes) {
            // Move the scanline to the new node's left border
            scanlinePos = currNode.rectangle.x;
            
            // Remove intersecting node which do not intersect the scanline anymore
            while (!intersectingNodes.isEmpty()) {
                RectangleNode intersectingRectangle = intersectingNodes.first();
                
                if (intersectingRectangle.rectangle.x + intersectingRectangle.rectangle.width < scanlinePos) {
                    intersectingNodes.remove(intersectingRectangle);
                } else {
                    // Since we iterate over intersecting nodes by the coordinate of their right border, once we have
                    // found a node the scanline has not moved past, we can stop looking
                    break;
                }
            }
            
            // Add overlaps between the currently intersecting nodes and the new node
            for (RectangleNode intersectingNode : intersectingNodes) {
                intersectingNode.overlappingNodes.add(currNode);
                currNode.overlappingNodes.add(intersectingNode);
            }
            
            // The new node is now part of the set of rectangles that intersect the scanline
            intersectingNodes.add(currNode);
        }
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Utility Methods

    /**
     * Compares two rectangles by the coordinates of their left borders.
     */
    private static int compareLeftRectangleBorders(final RectangleNode rn1, final RectangleNode rn2) {
        return Double.compare(rn1.rectangle.x, rn2.rectangle.x);
    }
    
    /**
     * Compares two rectangles by the coordinates of their right borders.
     */
    private static int compareRightRectangleBorders(final RectangleNode rn1, final RectangleNode rn2) {
        return Double.compare(rn1.rectangle.x + rn1.rectangle.width, rn2.rectangle.x + rn2.rectangle.width);
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Support Classes
    
    /**
     * The direction in which to move rectangles in order to remove overlaps. The direction also defines whether to
     * remove horizontal or vertical overlaps.
     */
    public static enum OverlapRemovalDirection {
        
        /** Remove horizontal overlaps by moving rectangles upwards. */
        UP,
        /** Remove horizontal overlaps by moving rectangles downwards. */
        DOWN,
        /** Remove vertical overlaps by moving rectangles leftwards. */
        LEFT,
        /** Remove vertical overlaps by moving rectangles rightwards. */
        RIGHT;
        
    }
    
    /**
     * A node in the rectangle overlap graph.
     */
    public static final class RectangleNode {
        
        /** The original rectangle represented by this node. */
        private ElkRectangle originalRectangle;
        /**
         * The rectangle represented by this node after the coordinate transformation. If no coordinate transformation
         * was necessary, this may well be the same rectangle as {@link #originalRectangle}.
         */
        private ElkRectangle rectangle;
        /** List of nodes that this node overlaps with. */
        private List<RectangleNode> overlappingNodes = Lists.newLinkedList();
        
        
        /**
         * Creates a new instance holding the given data.
         * 
         * @param originalRectangle
         *            the original rectangle represented by this node.
         * @param rectangle
         *            the rectangle transformed into our standard coordinate space.
         */
        private RectangleNode(final ElkRectangle originalRectangle, final ElkRectangle rectangle) {
            this.originalRectangle = originalRectangle;
            this.rectangle = rectangle;
        }


        /**
         * Returns this node's rectangle.
         */
        public ElkRectangle getRectangle() {
            return rectangle;
        }

        /**
         * Returns the nodes whose rectangles overlap with this node's rectangle.
         */
        public List<RectangleNode> getOverlappingNodes() {
            return overlappingNodes;
        }
        
    }
    
}
