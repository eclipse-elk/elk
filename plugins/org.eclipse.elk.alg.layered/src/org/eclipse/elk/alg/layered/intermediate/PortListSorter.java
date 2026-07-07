/*******************************************************************************
 * Copyright (c) 2011, 2018 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate;

import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.function.Function;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.PortSortingStrategy;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.core.util.Pair;

/**
 * Sorts the port lists of nodes with fixed port orders or fixed port positions. 
 * The node's list of ports is sorted beginning at the leftmost northern port, going clockwise.
 * This order of ports may be used during crossing minimization for calculating port ranks.
 * 
 * In case of {@link PortConstraints#FIXED_SIDE FIXED_SIDE} the ports are order by side
 * according to the order north - east - south - west (and depending on the value of the {@link PortSortingStrategy} 
 * layout option)
 * In case of {@link PortConstraints#FIXED_ORDER FIXED_ORDER} the side and
 * {@link LayeredOptions#PORT_INDEX PORT_INDEX} are used if specified. Otherwise the order is inferred
 * from specified port positions. For {@link PortConstraints#FIXED_POS FIXED_POS} solely the position of
 * the ports determines the order.
 * 
 * <dl>
 *   <dt>Precondition:</dt><dd>a layered graph.</dd>
 *   <dt>Postcondition:</dt><dd>the port lists of nodes are sorted.</dd>
 *   <dt>Slots:</dt><dd>Before phase 3.</dd>
 *   <dt>Same-slot dependencies:</dt><dd>{@link NorthSouthPortPreprocessor}.</dd>
 * </dl>
 * 
 * @see LNode#getPorts()
 */
public final class PortListSorter implements ILayoutProcessor<LGraph> {

    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor monitor) {
        monitor.begin("Port order processing", 1);

        PortSortingStrategy pss = layeredGraph.getProperty(LayeredOptions.PORT_SORTING_STRATEGY);
        
        // Iterate through the nodes of all layers
        for (Layer layer : layeredGraph) {
            for (LNode node : layer) {
                
                PortConstraints portConstraints = node.getProperty(LayeredOptions.PORT_CONSTRAINTS);
                List<LPort> ports = node.getPorts();
                
                // We need to sort the port list according to port constraints
                if (portConstraints.isOrderFixed()) {
                    Collections.sort(ports, CMP_COMBINED);
                } else if (portConstraints.isSideFixed()) {
                    // by default try to preserve the input order
                    Collections.sort(ports, CMP_PORT_SIDE);
                    reverseWestAndSouthSide(ports);
                    
                    if (pss == PortSortingStrategy.PORT_DEGREE) {
                        Collections.sort(ports, CMP_PORT_DEGREE_EAST_WEST);
                    }
                }
                node.cachePortSides();
            }
        }

        monitor.done();
    }
    
    private void reverseWestAndSouthSide(final List<LPort> ports) {
        if (ports.size() <= 1) {
            return;
        }
        
        Pair<Integer, Integer> southIndices = findPortSideRange(ports, PortSide.SOUTH);
        reverse(ports, southIndices.getFirst(), southIndices.getSecond());
        
        Pair<Integer, Integer> westIndices = findPortSideRange(ports, PortSide.WEST);
        reverse(ports, westIndices.getFirst(), westIndices.getSecond());
    }
    
    /**
     * Identifies indices {@code [i, j)} for the <em>pre-sorted</em> list {@code ports} such that 
     * {@code i} is the first index of a port on the given {@code side} and {@code j} is either an index 
     * of a port of a different side or equal to the list's size. 
     * If no ports are on the given side, {@code i} and {@code j} will be equal.  
     */
    private Pair<Integer, Integer> findPortSideRange(final List<LPort> ports, final PortSide side) {
        if (ports.isEmpty()) {
            return Pair.of(0, 0);
        }
            
        PortSide currentSide = ports.get(0).getSide();
        int lowIdx = 0;
        int lb = side.ordinal();
        int hb = side.ordinal() + 1;

        while (lowIdx < ports.size() - 1 && currentSide.ordinal() < lb) {
            lowIdx++;
            currentSide = ports.get(lowIdx).getSide();
        }
        int highIdx = lowIdx;
        while (highIdx < ports.size() - 1 && currentSide.ordinal() < hb) {
            highIdx++;
            currentSide = ports.get(lowIdx).getSide();
        }
    
        return Pair.of(lowIdx, highIdx);
    }
    
    /**
     * Reverses the {@code ports} list between {@code lowerIdx} (inclusive) and {@code higherIdx} (exclusive).
     * That is, given a sequence {@code p_0, ..., p_i, ..., p_j, ..., p_n}, where {@code i} and {@code j} are the 
     * low and high index, respectively, the method alters the list to 
     * {@code p_0, ..., p_j, p_{j-1}, ..., p_{i+1}, p_i, ..., p_n}.    
     */
    private void reverse(final List<LPort> ports, final int lowIdx, final int highIdx) {
        if (highIdx <= lowIdx + 2) {
            return;
        }
        int n = (highIdx - lowIdx) / 2;
        for (int i = 0; i < n; ++i) {
            LPort tmp = ports.get(lowIdx + i);
            ports.set(lowIdx + i, ports.get(highIdx - i - 1));
            ports.set(highIdx - i - 1, tmp);
        }
    }
    
        
    private static final Function<LPort, Iterable<LEdge>> IN_EDGES = p -> p.getIncomingEdges();
    private static final Function<LPort, Iterable<LEdge>> OUT_EDGES = p -> p.getOutgoingEdges();
    
    
    /**
     * Comparator for sorting a node's ports according to the {@link PortSide}: north - east - south - west. The order 
     * for the individual sides is not affected by this comparator. It can, however, be combined with other comparators,
     * like so: 
     * <pre>
     * CMP_PORT_SIDE.thenComparing(CMP_FIXED_ORDER_AND_FIXED_POS)
     * </pre>
     */
    public static final Comparator<LPort> CMP_PORT_SIDE = (p1, p2) -> {
        // Sort by side first (if the comparison ends here, the ports were on different sides;
        // otherwise, the ports must be on the same side)
        int ordinalDifference = p1.getSide().ordinal() - p2.getSide().ordinal();
        if (ordinalDifference != 0) {
            return ordinalDifference;
        }
                
        return 0;
    };
    
    /**
     * Comparator for sorting a node's ports such that on each the {@link PortSide#EAST} side ports are further
     * separated based on their out-degree (in the original input graph that is), and on the {@link PortSide#WEST} side
     * based on their in-degree. The comparator assumes that the ports are pre-sorted according to their sides (e.g. by
     * using {@link #CMP_PORT_SIDE}).
     */
    public static final Comparator<LPort> CMP_PORT_DEGREE_EAST_WEST = (p1, p2) -> {

        int ordinalDifference = p1.getSide().ordinal() - p2.getSide().ordinal();
        if (ordinalDifference != 0) {
            return 0; // ports on different sides -- not our job
        }
        
        switch (p1.getSide()) {
        case EAST:
            return realDegree(p2, OUT_EDGES) - realDegree(p1, OUT_EDGES);
            
        case WEST: 
            return realDegree(p1, IN_EDGES) - realDegree(p2, IN_EDGES);
        }
        
        return 0;
    };

    private static int realDegree(final LPort p, final Function<LPort, Iterable<LEdge>> edgesFun) {
        int realDegree = 0; 
        for (LEdge e : edgesFun.apply(p)) {
            if (!e.getProperty(InternalProperties.REVERSED)) {
                realDegree++;
            }
        }
        return realDegree;
    }

    /**
     * Comparator that handles the order of ports that lie on the same side when {@link PortConstraints} are
     * {@link PortConstraints#FIXED_ORDER} or {@link PortConstraints#FIXED_POS}. Note: The comparator exports the
     * collection of ports to be sorted according to their sides, which can be achieved by using {@link #CMP_PORT_SIDE}.
     */
    public static final Comparator<LPort> CMP_FIXED_ORDER_AND_FIXED_POS = (p1, p2) -> {
        
        final PortConstraints portConstraints =  p1.getNode().getProperty(LayeredOptions.PORT_CONSTRAINTS);
        
        // Do nothing if not at least the port order is fixed or the ports are on different sides
        int ordinalDifference = p1.getSide().ordinal() - p2.getSide().ordinal();
        if (ordinalDifference != 0 || !portConstraints.isOrderFixed()) {
            return 0;
        }

        // If the ports are on the same side and the node has FIXED_ORDER port constraints (that is,
        // the coordinates of the ports don't necessarily make sense), we check if the port index
        // has been explicitly set
        if (portConstraints == PortConstraints.FIXED_ORDER) {
            // In case of equal sides, sort by port index property
            Integer index1 = p1.getProperty(LayeredOptions.PORT_INDEX);
            Integer index2 = p2.getProperty(LayeredOptions.PORT_INDEX);
            if (index1 != null && index2 != null) {
                int indexDifference = index1 - index2;
                if (indexDifference != 0) {
                    return indexDifference;
                }
            }
        }

        // In case of equal index (or FIXED_POS), sort by position
        switch (p1.getSide()) {
        case NORTH:
            // Compare x coordinates
            return Double.compare(p1.getPosition().x, p2.getPosition().x);

        case EAST:
            // Compare y coordinates
            return Double.compare(p1.getPosition().y, p2.getPosition().y);

        case SOUTH:
            // Compare x coordinates in reversed order
            return Double.compare(p2.getPosition().x, p1.getPosition().x);

        case WEST:
            // Compare y coordinates in reversed order
            return Double.compare(p2.getPosition().y, p1.getPosition().y);

        default:
            // Port sides should not be undefined
            throw new IllegalStateException("Port side is undefined");
        }
    };
    
    /**
     * Combination of {@link #CMP_PORT_SIDE} and {@link #CMP_FIXED_ORDER_AND_FIXED_POS}. Should be the go to comparator
     * outside of the {@link PortListSorter}.
     */
    public static final Comparator<LPort> CMP_COMBINED = 
            CMP_PORT_SIDE.thenComparing(CMP_FIXED_ORDER_AND_FIXED_POS);
}