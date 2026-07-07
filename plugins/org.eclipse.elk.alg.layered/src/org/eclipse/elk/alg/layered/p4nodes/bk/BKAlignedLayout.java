/*******************************************************************************
 * Copyright (c) 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p4nodes.bk;

import java.util.Arrays;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.Spacings;

/**
 * Class which holds all information about a layout in one of the four direction
 * combinations.
 * 
 * @author jjc
 * @author uru
 */
public final class BKAlignedLayout {
    
    // Allow the fields of this container to be accessed from package siblings.
    // SUPPRESS CHECKSTYLE NEXT 28 VisibilityModifier
    /** The root node of each node in a block. */
    LNode[] root;
    /** The size of a block. */
    Double[] blockSize;
    /** The next node in a block, or the first if the current node is the last, forming a ring. */
    LNode[] align;
    /** The value by which a node must be shifted to stay straight inside a block. */
    Double[] innerShift;
    /** The root node of a class, mapped from block root nodes to class root nodes. */
    LNode[] sink;
    /** The value by which a block must be shifted for a more compact placement. */
    Double[] shift;
    /** The y-coordinate of every node, forming the final layout. */
    Double[] y;
    /** The vertical direction of the current layout. */
    VDirection vdir;
    /** The horizontal direction of the current layout. */
    HDirection hdir;
    /** Flags blocks, represented by their root node, that are part of a straightened edge. */
    Boolean[] su;
    /** Flags blocks, represented by their root node, that they are solely made up of dummy nodes. */
    Boolean[] od;

    /** The graph to process. */
    LGraph layeredGraph;
    /** Spacing values. */
    Spacings spacings;
    
    /**
     * Basic constructor for a layout.
     * 
     * @param layeredGraph
     *            the layered graph.
     * @param nodeCount
     *            number of nodes in this layout
     * @param vdir
     *            vertical traversal direction of the algorithm
     * @param hdir
     *            horizontal traversal direction of the algorithm
     */
    public BKAlignedLayout(final LGraph layeredGraph, final int nodeCount, 
            final VDirection vdir, final HDirection hdir) {

        this.layeredGraph = layeredGraph;
        // Initialize spacing value from layout options.
        spacings = layeredGraph.getProperty(InternalProperties.SPACINGS);
        
        root = new LNode[nodeCount];
        blockSize = new Double[nodeCount];
        align = new LNode[nodeCount];
        innerShift = new Double[nodeCount];
        sink = new LNode[nodeCount];
        shift = new Double[nodeCount];
        y = new Double[nodeCount];
        su = new Boolean[nodeCount];
        Arrays.fill(su, false);
        od = new Boolean[nodeCount];
        Arrays.fill(od, true);
        this.vdir = vdir;
        this.hdir = hdir;
    }

    /**
     * Explicitly release any allocated resources.
     */
    public void cleanup() {
        root = null;
        blockSize = null;
        align = null;
        innerShift = null;
        sink = null;
        shift = null;
        y = null;
    }
    
    /**
     * Calculate the layout size for comparison.
     * 
     * @return The size of the layout
     */
    public double layoutSize() {
        double min = Double.POSITIVE_INFINITY;
        double max = Double.NEGATIVE_INFINITY;
        // Prior to KIPRA-1426 the size of the layout was determined  
        // only based on y coordinates, neglecting any block sizes.
        // We now determine the maximal extend of the layout based on
        // the minimum y coordinate of any node and the maximum
        // y coordinate _plus_ the size of any block.
        for (Layer layer : layeredGraph.getLayers()) {
            for (LNode n : layer.getNodes()) {
                double yMin = y[n.id];
                double yMax = yMin + blockSize[root[n.id].id];
                min = Math.min(min, yMin);
                max = Math.max(max, yMax);
            }
        }
        return max - min;
    }

    
    /**
     * @param src
     *            source port of the tested edge
     * @param tgt
     *            target port of the tested edge
     * @return A delta larger than 0 if the {@code tgt} port has a larger y coordinate than
     *         {@code src} and a delta smaller than zero if {@code src} has the larger y coordinate.
     *         This means that for {@code delta > 0} the target node has to be shifted upwards to
     *         straighten the edge.
     */
    public double calculateDelta(final LPort src, final LPort tgt) {
        double srcPos = y[src.getNode().id] + innerShift[src.getNode().id] 
                + src.getPosition().y + src.getAnchor().y;
        double tgtPos = y[tgt.getNode().id] + innerShift[tgt.getNode().id] 
                + tgt.getPosition().y + tgt.getAnchor().y;
        return tgtPos - srcPos;
    }

    /**
     * Shifts the y-coordinates of all nodes of the block represented by {@code root} by the
     * specified {@code delta}.
     * 
     * @param rootNode
     *            root node of the block.
     * @param delta
     *            the delta by which the node should be move. Can be either positive or negative.
     */
    public void shiftBlock(final LNode rootNode, final double delta) {
        LNode current = rootNode;
        do {
            double newPos = y[current.id] + delta;
            y[current.id] = newPos;
            current = align[current.id];
        } while (current != rootNode);
    }
    
    /**
     * Checks whether a block with root node {@code blockRoot} can be shifted upwards by
     * {@code delta} without overlapping any of the block's nodes upper neighbors.
     * 
     * @param blockRoot
     *            root node of a block
     * @param delta
     *            a positive value
     * @return A value smaller or equal to {@code delta} indicating the maximal distance the
     *         block can be moved upward.
     */
    public double checkSpaceAbove(final LNode blockRoot, final double delta, final NeighborhoodInformation ni) {
        
        double availableSpace = delta;
        final LNode rootNode = blockRoot;
        // iterate through the block
        LNode current = rootNode;
        do {
            current = align[current.id];
            // get minimum possible position of the current node
            double minYCurrent = getMinY(current);

            LNode neighbor = getUpperNeighbor(current, ni);
            if (neighbor != null) {
                double maxYNeighbor = getMaxY(neighbor);
                // minimal position at which the current block node could validly be placed
                availableSpace =
                        Math.min(availableSpace,
                                minYCurrent 
                                - (maxYNeighbor + spacings.getVerticalSpacing(current, neighbor)));
            }
            // until we wrap around
        } while (rootNode != current);

        return availableSpace;
    }
    
    /**
     * Checks whether a block with root node {@code blockRoot} can be shifted downwards by
     * {@code delta} without overlapping any of the block's nodes lower neighbors.
     * 
     * @param blockRoot
     *            root node of a block
     * @param delta
     *            a positive value
     * @return A value smaller or equal to {@code delta} indicating the maximal distance the
     *         block can be moved upward.
     */
    public double checkSpaceBelow(final LNode blockRoot, final double delta, final NeighborhoodInformation ni) {

        double availableSpace = delta;
        final LNode rootNode = blockRoot;
        // iterate through the block
        LNode current = rootNode;
        do {
            current = align[current.id];
            // get maximum possible position of the current node
            double maxYCurrent = getMaxY(current);

            // get the lower neighbor and check its position allows shifting
            LNode neighbor = getLowerNeighbor(current, ni);
            if (neighbor != null) {
                double minYNeighbor = getMinY(neighbor);

                // minimal position at which the current block node could validly be placed
                availableSpace = Math.min(availableSpace,
                                minYNeighbor 
                                - (maxYCurrent + spacings.getVerticalSpacing(current, neighbor)));
            }
            // until we wrap around
        } while (rootNode != current);

        return availableSpace;
    }
    
    /**
     * Returns the minimum position of node {@code n} and its margins, that is,
     * {@code node.y + node.innerShift - node.margin.top}. Note that no spacing is accounted for.
     * 
     * @param n
     *            a node
     * @return the minimum position.
     */
    public double getMinY(final LNode n) {

        // node size + margins + inside shift etc
        LNode rootNode = root[n.id];
        return y[rootNode.id]
            + innerShift[n.id]
            - n.getMargin().top;
    }
    
    /**
     * Returns the maximum position of node {@code n} and its margins, that is,
     * {@code node.y + node.innerShift + node.size + node.margin.bottom}. Note that no spacing is
     * accounted for.
     * 
     * @param n
     *            a node
     * @return the minimum position.
     */
    public double getMaxY(final LNode n) {
        
        // node size + margins + inside shift etc
        LNode rootNode = root[n.id];
        return y[rootNode.id]
            + innerShift[n.id]
            + n.getSize().y
            + n.getMargin().bottom;
    }
    
    
    /**
     * @param n
     *            the node for which the neighbor is requested.
     * @param ni
     *            the {@link NeighborhoodInformation} of {@code n}.
     * @return the node with a <b>larger</b> y than {@code n} within {@code n}'s layer if it exists,
     *         otherwise {@code null}.
     */
    private LNode getLowerNeighbor(final LNode n, final NeighborhoodInformation ni) {
        Layer l = n.getLayer();
        int layerIndex = ni.nodeIndex[n.id];
        if (layerIndex < l.getNodes().size() - 1) {
            return l.getNodes().get(layerIndex + 1);
        }
        return null;
    }

    /**
     * @param n
     *            the node for which the neighbor is requested.
     * @param ni
     *            the {@link NeighborhoodInformation} of {@code n}.
     * @return the node with a <b>smaller</b> y than {@code n} within {@code n}'s layer if it
     *         exists, otherwise {@code null}.
     */
    private LNode getUpperNeighbor(final LNode n, final NeighborhoodInformation ni) {
        Layer l = n.getLayer();
        int layerIndex = ni.nodeIndex[n.id];
        if (layerIndex > 0) {
            return l.getNodes().get(layerIndex - 1); 
        }
        return null;
    }
    
    @Override
    public String toString() {
        String result = "";
        if (hdir == HDirection.RIGHT) {
            result += "RIGHT";
        } else if (hdir == HDirection.LEFT) {
            result += "LEFT";
        }
        if (vdir == VDirection.DOWN) {
            result += "DOWN";
        } else if (vdir == VDirection.UP) {
            result += "UP";
        } else {
            result += "BALANCED";
        }
        return result;
    }

    /////////////////////////////////////////////////////////////////////////////////////////////////////
    // Enumerations
    
    /**
     * Vertical direction enumeration.
     */
    public static enum VDirection {
        /** Iteration direction top-down. */
        DOWN, 
        /** Iteration direction bottom-up. */
        UP;
    }
    
    /**
     * Horizontal direction enumeration.
     */
    public static enum HDirection {
        /** Iterating from right to left. */
        RIGHT, 
        /** Iterating from left to right. */
        LEFT;
    }
}
