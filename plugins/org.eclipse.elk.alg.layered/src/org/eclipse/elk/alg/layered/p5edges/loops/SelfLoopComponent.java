/*******************************************************************************
 * Copyright (c) 2018 Kiel University and others.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * which accompanies this distribution, and is available at
 * http://www.eclipse.org/legal/epl-v10.html
 *******************************************************************************/
package org.eclipse.elk.alg.layered.p5edges.loops;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.core.options.PortSide;

/**
 * A self-loop component comprises a number of ports self-loops are incident to. Each can be thought of as a connected
 * component of a graph whose nodes are the ports, and whose edges are the self-loop edges. Each self loop component
 * can have labels and dependencies to other components it crosses.
 */
public class SelfLoopComponent {

    /** The ports of the self loop component. They are all connected by at least one edge. */
    private final List<SelfLoopPort> ports;
    /** All center edge labels that occur at the component edges, if any. */
    private SelfLoopLabel selfLoopLabel;
    /** The components that the self loop component depends on. */
    private final Map<SelfLoopNodeSide, List<SelfLoopComponent>> dependencyComponents = new HashMap<>();
    /** Any edges the component depends on. */
    private final Map<PortSide, List<SelfLoopEdge>> edgeDependencies = new HashMap<>();

    
    /**
     * Create a SelfLoopComponent by passing the ports of a node which are connected with each other.
     * 
     * @param ports
     *            Ports of the connected self loop component.
     */
    public SelfLoopComponent(final List<LPort> ports) {
        this.ports = new ArrayList<>(ports.size());
        addAllPorts(ports);
        calculateConnectedEdges();
    }

    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Component Calculation

    /**
     * Creates all self loop components for the given node and adds them to the node.
     */
    public static void createSelfLoopComponents(final SelfLoopNode slNode) {
        LNode lNode = slNode.getNode();
        List<LPort> nodePorts = new ArrayList<LPort>(lNode.getPorts());
        
        while (!nodePorts.isEmpty()) {
            // find all connected ports
            Set<LPort> connectedPorts = new LinkedHashSet<>();
            findAllConnectedPorts(lNode, nodePorts.get(0), connectedPorts);
            
            // sort the ports according to their original index at the node
            List<LPort> sortedPorts = new ArrayList<>(connectedPorts);
            sortedPorts.sort(Comparator.comparing(lNode.getPorts()::indexOf));
            
            // create the actual component
            SelfLoopComponent component = new SelfLoopComponent(sortedPorts);
            slNode.getSelfLoopComponents().add(component);

            nodePorts.removeAll(sortedPorts);
        }
    }

    /**
     * Given a port this method searches via incoming and outgoing edges others port which can be reached. Pass in an
     * empty set as the last argument.
     */
    private static void findAllConnectedPorts(final LNode node, final LPort lPort,
            final Set<LPort> connectedPorts) {
        // only process nodes which where not visited yet
        if (!connectedPorts.add(lPort)) {
            return;
        }

        // follow all incoming edges to find connected ports
        for (LEdge edge : lPort.getIncomingEdges()) {
            LPort source = edge.getSource();
            if (source.getNode().equals(node)) {
                findAllConnectedPorts(node, source, connectedPorts);
            }
        }
        
        // follow all outgoing edges to find connected ports
        for (LEdge edge : lPort.getOutgoingEdges()) {
            LPort target = edge.getTarget();
            if (target.getNode().equals(node)) {
                findAllConnectedPorts(node, target, connectedPorts);
            }
        }
    }

    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Ports

    /**
     * Adds a {@link SelfLoopPort} to the component that represents the given regular port.
     */
    private void addPort(final LPort port) {
        ports.add(new SelfLoopPort(port, this));
    }

    /**
     * Calls {@link #addPort(LPort)} on all ports in the list.
     */
    private void addAllPorts(final List<LPort> newPorts) {
        for (LPort port : newPorts) {
            addPort(port);
        }
    }

    /**
     * Returns the list of ports this component represents.
     */
    public List<SelfLoopPort> getPorts() {
        return ports;
    }

    /**
     * Find the corresponding {@link SelfLoopPort} for an {@link LPort} or {@code null} if it could not be found.
     */
    public SelfLoopPort findPort(final LPort port) {
        for (SelfLoopPort selfLoopPort : ports) {
            if (selfLoopPort.getLPort() == port) {
                return selfLoopPort;
            }
        }
        return null;
    }

    /**
     * Returns the port that follows the given port in the list of ports.
     */
    public SelfLoopPort getNextPort(final SelfLoopPort port) {
        int index = ports.indexOf(port);
        if (index != ports.size() - 1) {
            return ports.get(index + 1);
        } else {
            return null;
        }
    }

    /**
     * Returns the last port in the list of ports.
     */
    public SelfLoopPort getLastPort() {
        if (ports.isEmpty()) {
            return null;
        } else {
            return ports.get(ports.size() - 1);
        }
    }

    /**
     * Returns the set of port sides this component has ports on.
     */
    public Set<PortSide> getPortSides() {
        return ports.stream()
                .map(port -> port.getPortSide())
                .collect(Collectors.toSet());
    }

    /**
     * Returns a list of of ports that are on the given port side.
     */
    public List<SelfLoopPort> getPortsOfSide(final PortSide side) {
        return ports.stream()
                .filter(port -> port.getPortSide() == side)
                .collect(Collectors.toList());
    }

    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Labels

    /**
     * Return the self loop label. May be {@code null} if there are no labels attached to this component.
     */
    public SelfLoopLabel getSelfLoopLabel() {
        return selfLoopLabel;
    }

    /**
     * Set the label.
     */
    public void setSelfLoopLabel(final SelfLoopLabel selfLoopLabel) {
        this.selfLoopLabel = selfLoopLabel;
    }

    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Edges

    /**
     * Return the edges that belong to this component.
     */
    public Set<SelfLoopEdge> getConnectedEdges() {
        return ports.stream()
                .flatMap(port -> port.getConnectedEdges().stream())
                .collect(Collectors.toSet());
    }
    
    /**
     * Calculate the edges connecting the ports of the component. For each edge a SelfLoopEdge is created.
     */
    private void calculateConnectedEdges() {
        for (SelfLoopPort port : ports) {
            for (LEdge edge : port.getLPort().getOutgoingEdges()) {
                SelfLoopPort targetPort = findPort(edge.getTarget());
                if (targetPort != null) {
                    SelfLoopEdge selfLoopEdge = new SelfLoopEdge(this, port, targetPort, edge);
                    port.getConnectedEdges().add(selfLoopEdge);
                    targetPort.getConnectedEdges().add(selfLoopEdge);
                }
            }
        }
    }

    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Dependencies

    /**
     * Returns the component dependencies for each self loop node side.
     */
    public Map<SelfLoopNodeSide, List<SelfLoopComponent>> getDependencyComponents() {
        return dependencyComponents;
    }

    /**
     * Returns the edge dependencies for each port side.
     */
    public Map<PortSide, List<SelfLoopEdge>> getEdgeDependencies() {
        return edgeDependencies;
    }

    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Component Type

    /**
     * Returns a list of all the port sides this component spans.
     */
    public List<PortSide> getComponentSpanningSides() {
        Set<PortSide> sides = new LinkedHashSet<>();
        PortSide currentPortSide = ports.get(0).getPortSide();
        PortSide lastPortSide = ports.get(ports.size() - 1).getPortSide();
        sides.add(lastPortSide);

        if (currentPortSide == lastPortSide
                && ports.get(0).getOriginalIndex() > ports.get(ports.size() - 1).getOriginalIndex()) {
            do {
                sides.add(currentPortSide);
                currentPortSide = currentPortSide.right();
            } while (currentPortSide != lastPortSide);
        } else {
            while (currentPortSide != lastPortSide) {
                sides.add(currentPortSide);
                currentPortSide = currentPortSide.right();
            }
        }
        return new ArrayList<PortSide>(sides);
    }

    /**
     * Returns {@code true} or {@code false} as this component does or does not span exactly two adjacent port sides.
     */
    public boolean isCornerComponent() {
        Set<PortSide> portSide = getPortSides();
        
        if (portSide.size() == 2) {
            Iterator<PortSide> portSideIterator = portSide.iterator();
            PortSide side1 = portSideIterator.next();
            PortSide side2 = portSideIterator.next();
            
            return side1.areAdjacent(side2);
        } else {
            return false;
        }
    }

    /**
     * Returns {@code true} or {@code false} as this component does or does not span exactly two opposing port sides.
     */
    public boolean isOpposingComponent() {
        return getPortSides().size() == 2 && !isCornerComponent();
    }

    /**
     * Determines this component's component type.
     */
    public SelfLoopType getType(final SelfLoopNode nodeRep) {
        if (ports.size() == 1) {
            return SelfLoopType.NON_LOOP;
        }
        
        SelfLoopPort source = ports.get(0);
        PortSide sourceSide = source.getPortSide();
        SelfLoopPort target = ports.get(ports.size() - 1);
        PortSide targetSide = target.getPortSide();
        
        boolean rightDir = source.getDirection() == SelfLoopRoutingDirection.RIGHT;
        boolean leftDir = source.getDirection() == SelfLoopRoutingDirection.LEFT;

        if (source.getLPort().getNode() != target.getLPort().getNode()) {
            // The ports don't belong to the same node in the first place
            return SelfLoopType.NON_LOOP;
            
        } else if (sourceSide == targetSide) {
            // We have a self-loop, either routed on the same side or once around the node
            SelfLoopNodeSide nodeRepside = nodeRep.getNodeSide(sourceSide);
            int sourceIndex = nodeRepside.getPorts().indexOf(source);
            int targetIndex = nodeRepside.getPorts().indexOf(target);
            
            if (leftDir && sourceIndex < targetIndex || rightDir && targetIndex < sourceIndex) {
                return SelfLoopType.FOUR_CORNER;
            } else {
                return SelfLoopType.SIDE;
            }
            
        } else if (sourceSide.left() == targetSide || sourceSide.right() == targetSide) {
            // We either have a corner self-loop or one routed almost completely around the node
            if (leftDir && source.getPortSide().right() == target.getPortSide()
                    || rightDir && source.getPortSide().left() == target.getPortSide()) {
                return SelfLoopType.THREE_CORNER;
            } else {
                return SelfLoopType.CORNER;
            }
        }

        return SelfLoopType.OPPOSING;
    }

    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Miscellaneous

    @Override
    public String toString() {
        StringBuilder builder = new StringBuilder();
        builder.append("(");
        for (SelfLoopPort port : ports) {
            builder.append(port);
            if (ports.indexOf(port) != ports.size() - 1) {
                builder.append(", ");
            }
        }
        builder.append(")");
        return builder.toString();
    }
}
