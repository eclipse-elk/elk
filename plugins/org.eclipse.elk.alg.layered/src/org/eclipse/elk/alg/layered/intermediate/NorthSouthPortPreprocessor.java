/*******************************************************************************
 * Copyright (c) 2011, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate;

import java.util.ArrayList;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LGraphUtil;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.OrderingStrategy;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.IElkProgressMonitor;

import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;

/**
 * Inserts dummy nodes to cope with northern and southern ports.
 * 
 * <p>
 * For each node with northern or southern ports, dummy nodes are inserted in the node's layer above
 * or below the node. (for northern and southern ports, respectively) Assume that a node has several
 * northern ports. First, the ports are assembled in left-to-right order according to their
 * position. Then, the ports are partitioned into ports with only incoming edges (in ports), ports
 * with only outgoing edges (out ports) and ports with both, incoming and outgoing edges (in/out
 * ports). The way the dummy nodes are then created can either be according to the
 * <em>old approach</em> or the <em>new approach</em>.
 * </p>
 * 
 * 
 * <h3>The Old Approach</h3>
 * 
 * <p>
 * In and out ports are matched up left to right and right to left, respectively, as long as the out
 * port is right of the in port. In the example below, ports P1 and P6 will be matched, as will be
 * ports P2 and P5. Ports P3 and P4 will not be matched since P3 is an out port and not right of P4.
 * </p>
 * 
 * <p>
 * For each pair of matched ports, a dummy node is then inserted that is shared by both ports. This
 * leads to edges of both ports sharing the same y coordinate, which is fine as long as the out port
 * is right of the in port. When this is no longer the case, ports cannot share a dummy node
 * anymore. (the horizontal edge segments would overlap and thus introduce ambiguity) Thus, P4 and
 * P3 each get their own dummy node.
 * </p>
 * 
 * <p>
 * Once in and out ports are processed, in/out ports are given their own dummy node.
 * </p>
 * 
 * <p>
 * This approach minimizes the amount of edge crossings among the edges connected to northern ports
 * locally, without paying attention to the rest of the graph. To do that, this approach fixes the
 * order of the dummy nodes in the layer: the node successor constraints are set in a way that each
 * dummy node points to the next, with the bottom-most northern dummy node pointing to the regular
 * node it was created for.
 * </p>
 * 
 * <pre>
 *                      ------------------------------
 *                      |
 * ---------------------+----
 *                      |   |
 * ------------------   |   |   ----------------------
 *                  |   |   |   |
 * --------------   |   |   |   |   ------------------
 *              |   |   |   |   |   |
 *              P1  P2  P3  P4  P5  P6
 * </pre>
 * 
 * 
 * <h3>The New Approach</h3>
 * 
 * <p>
 * The new approach is designed to not fix the order of dummy nodes at this point and to instead
 * leave the task of ordering them to the crossing minimization phase. It creates a single dummy
 * node for every northern port, without two ports sharing a dummy node. The only necessary
 * constraints are that dummy nodes created for northern ports must have their regular node as their
 * successor. (Similarly, the regular node must have all dummy nodes created for southern ports as
 * successors.)
 * 
 * 
 * <h3>Self-Loops</h3>
 * 
 * <p>
 * Self-loops are a special case that is handled partly by this processor. For this to work, the
 * {@link DummySelfLoopProcessor} must have been executed prior to this processor's execution. Then, this
 * processor correctly processes all kinds of self-loops involving northern or southern ports as
 * follows: (due to the {@link DummySelfLoopProcessor}, only the cases detailled below must be handled)
 * </p>
 * 
 * <ul>
 * <li>West-north, west-south, north-east and south-east self-loops are handled just like every
 * other connection.</li>
 * <li>North-north and south-south self-loops are handled by creating a dummy node just for that
 * self-loop.</li>
 * <li>North-south self-loops are handled by adding a northern and a southern dummy, both having
 * just one port via which they are connected. At the moment, the port is created on the eastern
 * side of the dummy node, but that could be handled more clerverly later.</li>
 * </ul>
 * 
 * <dl>
 *   <dt>Precondition:</dt>
 *     <dd>a layered graph; nodes have fixed port sides.</dd>
 *   <dt>Postcondition:</dt>
 *     <dd>dummy nodes have been inserted for edges connected to ports on north and south sides, with
 *         layout groups and node successor constraints set.</dd>
 *   <dt>Slots:</dt>
 *     <dd>Before phase 3.</dd>
 *   <dt>Same-slot dependencies:</dt>
 *     <dd>{@link PortListSorter}</dd>
 *     <dd>{@link DummySelfLoopProcessor}</dd>
 * </dl>
 * 
 * @see NorthSouthPortPostprocessor
 * @see PortListSorter
 * @see DummySelfLoopProcessor
 */
public final class NorthSouthPortPreprocessor implements ILayoutProcessor<LGraph> {

    /**
     * Whether the new or the old approach should be used. The new approach does not impose an order
     * on the dummy nodes created for northern and southern ports.
     */
    private static final boolean USE_NEW_APPROACH = true;

    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor monitor) {
        monitor.begin("Odd port side processing", 1);

        int pointer;
        List<LNode> northDummyNodes = Lists.newArrayList();
        List<LNode> southDummyNodes = Lists.newArrayList();

        // Iterate through the layers
        for (Layer layer : layeredGraph) {
            // The pointer indicates the index of the current node while northern ports are
            // processed,
            // and the index of the most recently inserted dummy while south ports are processed
            pointer = -1;

            // Iterate through the nodes (use an array to avoid concurrent modification exceptions)
            LNode[] nodeArray = LGraphUtil.toNodeArray(layer.getNodes());
            for (LNode node : nodeArray) {
                pointer++;

                // We only care about non-dummy nodes with fixed port sides
                if (!(node.getType() == NodeType.NORMAL
                        && node.getProperty(LayeredOptions.PORT_CONSTRAINTS).isSideFixed())) {

                    continue;
                }

                // Sort the port list if we have control over the port order
                if (!node.getProperty(LayeredOptions.PORT_CONSTRAINTS).isOrderFixed()
                        && node.getGraph().getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY) == OrderingStrategy.NONE) {
                    sortPortList(node);
                }

                // Nodes form their own layout unit
                node.setProperty(InternalProperties.IN_LAYER_LAYOUT_UNIT, node);

                // Clear the lists of northern and southern dummy nodes
                northDummyNodes.clear();
                southDummyNodes.clear();

                // Create a list of barycenter associates for the node
                List<LNode> barycenterAssociates = Lists.newArrayList();

                // Prepare a list of ports on the northern side, sorted from left to right (when viewed
                // in the diagram); create the appropriate dummy nodes and assign them to the layer
                LinkedList<LPort> portList = Lists.newLinkedList();
                Iterables.addAll(portList, node.getPorts(PortSide.NORTH));

                if (node.getGraph().getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY) != OrderingStrategy.NONE) {                    
                    portList = modelOrderNorthSouthInputReversing(portList, node);
                }

                createDummyNodes(layeredGraph, portList, northDummyNodes, southDummyNodes,
                        barycenterAssociates);

                /*
                 * We will now iterate over all nodes and insert successor constraints to each one
                 * of them. If we're using the old constrained approach, each dummy has the next
                 * dummy as its successor. If we use the new unconstrained approach, each dummy on
                 * the northern side has its regular node as a successor while the node has all
                 * southern dummies as its successors.
                 */

                int insertPoint = pointer;
                LNode successor = node;
                for (LNode dummy : northDummyNodes) {
                    dummy.setLayer(insertPoint, layer);
                    pointer++;

                    // The dummy nodes form a layout unit identified by the node they were created from.
                    // In addition, northern dummy nodes must appear before the regular node
                    dummy.setProperty(InternalProperties.IN_LAYER_LAYOUT_UNIT, node);
                 
                    // If originPort has port constraint ALLOW_NON_FLOW_PORTS_TO_SWITCH_SIDES,
                    // do not apply successor constraints to the dummy node dummy.
                    // Their position will be determined according to their barycenter value.
                    // If originPort does not have the port constraint ALLOW_NON_FLOW_PORTS_TO_SWITCH_SIDES,
                    // the dummy node dummy needs to appear before the regular node.

                    // Each dummy node has least one port (there may be two if an odd port
                    // has both an incoming and an outgoing edge, however the origin is the same)
                    assert dummy.getPorts().size() >= 1;
                    LPort dummyPort = dummy.getPorts().get(0);
                    // The port the dummy node was created for
                    LPort originPort = (LPort) dummyPort.getProperty(InternalProperties.ORIGIN);

                    if (!originPort.getProperty(LayeredOptions.ALLOW_NON_FLOW_PORTS_TO_SWITCH_SIDES)) {
                        dummy.getProperty(InternalProperties.IN_LAYER_SUCCESSOR_CONSTRAINTS).add(
                                successor);
                    }

                    if (!USE_NEW_APPROACH) {
                        // The old approach needs the successor to always point to the most recently
                        // created dummy node
                        successor = dummy;
                    }
                }

                // Do the same for ports on the southern side; the list of ports must be built in
                // reversed order, since ports on the southern side are listed from right to left
                portList.clear();
                for (LPort port : node.getPorts(PortSide.SOUTH)) {
                    portList.addFirst(port);
                }

                if (node.getGraph().getProperty(LayeredOptions.CONSIDER_MODEL_ORDER_STRATEGY) != OrderingStrategy.NONE) {                    
                    portList = modelOrderNorthSouthInputReversing(portList, node);
                }

                createDummyNodes(layeredGraph, portList, southDummyNodes, null,
                        barycenterAssociates);

                LNode predecessor = node;
                for (LNode dummy : southDummyNodes) {
                    dummy.setLayer(++pointer, layer);

                    // The dummy nodes form a layout unit identified by the node they were created from.
                    // In addition, southern dummy nodes must appear after the regular node
                    dummy.setProperty(InternalProperties.IN_LAYER_LAYOUT_UNIT, node);

                    // If originPort has port constraint ALLOW_NON_FLOW_PORTS_TO_SWITCH_SIDES,
                    // do not apply successor constraints to the dummy node dummy.
                    // Their position will be determined according to their barycenter value.
                    // If originPort does not have the port constraint ALLOW_NON_FLOW_PORTS_TO_SWITCH_SIDES,
                    // the dummy node dummy needs to appear before the regular node.

                    // Each dummy node has at least one port (there may be two if an odd port
                    // has both an incoming and an outgoing edge, however the origin is the same)
                    assert dummy.getPorts().size() >= 1;
                    LPort dummyPort = dummy.getPorts().get(0);
                    // The port the dummy node was created for
                    LPort originPort = (LPort) dummyPort.getProperty(InternalProperties.ORIGIN);

                    if (!originPort.getProperty(LayeredOptions.ALLOW_NON_FLOW_PORTS_TO_SWITCH_SIDES)) {
                        predecessor.getProperty(InternalProperties.IN_LAYER_SUCCESSOR_CONSTRAINTS).add(
                                dummy);
                    }

                    if (!USE_NEW_APPROACH) {
                        // The old approach needs the predecessor to always point to the most
                        // recently
                        // created dummy node
                        predecessor = dummy;
                    }
                }

                // If the list of barycenter associates contains nodes, set the appropriate property
                if (!barycenterAssociates.isEmpty()) {
                    node.setProperty(InternalProperties.BARYCENTER_ASSOCIATES, barycenterAssociates);
                }
            }
        }

        monitor.done();
    }

    // /////////////////////////////////////////////////////////////////////////////
    // PORT LIST SORTING

    /**
     * Sorts the list of northern and southern ports such that ports with only incoming edges end up
     * left, ports with only outgoing edges end up right, and ports with both end up in between.
     * Ports on the eastern and western sides are left untouched.
     * 
     * @param node
     *            the node whose ports to sort.
     */
    private void sortPortList(final LNode node) {
        int ports = node.getPorts().size();

        // Next IDs for ports with a given configuration of input and output edges. The choice of
        // initial IDs ensures that port IDs will be unique
        int inPorts = 0;
        int inOutPorts = ports;
        int outPorts = 2 * ports;

        // Iterate over the list of ports and set their IDs
        for (LPort port : node.getPorts()) {
            switch (port.getSide()) {
            case EAST:
            case WEST:
                port.id = -1;
                break;

            case NORTH:
            case SOUTH:
                int incoming = port.getIncomingEdges().size();
                int outgoing = port.getOutgoingEdges().size();

                if (incoming > 0 && outgoing > 0) {
                    port.id = inOutPorts++;
                } else if (incoming > 0) {
                    port.id = inPorts++;
                } else if (outgoing > 0) {
                    port.id = outPorts++;
                } else {
                    // Unconnected ports are placed between input ports...
                    port.id = inPorts++;
                }

                break;
            }
        }

        // With all IDs assigned, sort the port list
        Collections.sort(node.getPorts(), (port1, port2) -> {
            PortSide side1 = port1.getSide();
            PortSide side2 = port2.getSide();

            if (side1 != side2) {
                // sort according to the node side
                return side1.ordinal() - side2.ordinal();
            } else {
                if (port1.id == port2.id) {
                    // Eastern and western ports have the same ID and have to retain their order
                    return 0;
                } else {
                    if (side1 == PortSide.NORTH) {
                        return port1.id - port2.id;
                    } else {
                        return port2.id - port1.id;
                    }
                }
            }
        });
    }
    
    private LinkedList<LPort> modelOrderNorthSouthInputReversing(LinkedList<LPort> portList, LNode node) {
        // Reverse port list if edges are incoming edges.
        LinkedList<LPort> incoming = new LinkedList<LPort>();
        LinkedList<LPort> outgoing = new LinkedList<LPort>();
        for (LPort port : portList) {
            if (!port.getIncomingEdges().isEmpty()) {
                // Incoming edge
                incoming.add(port);
            } else {
                outgoing.add(port);
            }
        }
        Lists.reverse(incoming).addAll(outgoing);
        return incoming;
    }

    // /////////////////////////////////////////////////////////////////////////////
    // DUMMY NODE CREATION

    /**
     * Returns a list of dummy nodes for the given ports. The list of ports must be sorted by
     * position from left to right.
     * 
     * @param layeredGraph
     *            the layered graph.
     * @param ports
     *            the list of ports to create dummy nodes for.
     * @param dummyNodes
     *            all created dummy nodes are added to this list.
     * @param opposingSideDummyNodes
     *            all dummy nodes created due to north-south self-loops for the southern side are
     *            placed in this list. When called for southern ports, this may be {@code null}.
     * @param barycenterAssociates
     *            dummy nodes created for anything other than self-loops are put in this list to
     *            remember to include them in the barycenter calculations later.
     */
    private void createDummyNodes(final LGraph layeredGraph, final List<LPort> ports,
            final List<LNode> dummyNodes, final List<LNode> opposingSideDummyNodes,
            final List<LNode> barycenterAssociates) {

        // We'll assemble lists of ports with only incoming, ports with only outgoing
        // and ports with both, incoming and outgoing edges
        List<LPort> inPorts = new ArrayList<LPort>(ports.size());
        List<LPort> outPorts = new ArrayList<LPort>(ports.size());
        List<LPort> inOutPorts = new ArrayList<LPort>(ports.size());
        List<LEdge> sameSideSelfLoopEdges = new ArrayList<LEdge>(ports.size());
        List<LEdge> northSouthSelfLoopEdges = new ArrayList<LEdge>(ports.size());

        for (LPort port : ports) {
            // Go through the port's outgoing edges, looking for self-loops that need special
            // handling
            for (LEdge edge : port.getOutgoingEdges()) {
                // Check for self loops we'd be interested in
                if (edge.getSource().getNode() == edge.getTarget().getNode()) {
                    // Check which sides the ports are on
                    if (port.getSide() == edge.getTarget().getSide()) {
                        // Same side
                        sameSideSelfLoopEdges.add(edge);
                        continue;

                    } else if (port.getSide() == PortSide.NORTH
                            && edge.getTarget().getSide() == PortSide.SOUTH) {

                        // North->South self-loop. Due to the SelfLoopProcessor, a
                        // South->North self-loop cannot happen
                        northSouthSelfLoopEdges.add(edge);
                        continue;
                    }
                }
            }
        }

        // First, create the dummy nodes that handle north->south self-loops. For now,
        // we always route north->south self-loops east to the node. This could later
        // change, though.
        for (LEdge edge : northSouthSelfLoopEdges) {
            createNorthSouthSelfLoopDummyNodes(layeredGraph, edge, dummyNodes,
                    opposingSideDummyNodes, PortSide.EAST);
        }

        // Second, create the dummy nodes that handle same-side self-loops
        for (LEdge edge : sameSideSelfLoopEdges) {
            createSameSideSelfLoopDummyNode(layeredGraph, edge, dummyNodes);
        }

        // Now we iterate over the ports again, with certain self-loop edges already removed, and
        // check
        // if they are input ports, output ports, or both
        for (LPort port : ports) {
            // Find out if the port has incoming or outgoing edges
            boolean in = port.getIncomingEdges().size() > 0;
            boolean out = port.getOutgoingEdges().size() > 0;

            if (in && out) {
                inOutPorts.add(port);
            } else if (in) {
                inPorts.add(port);
            } else if (out) {
                outPorts.add(port);
            }
        }

        // We will now create dummy nodes for input ports and for output ports. How we create them
        // depends on whether we're using the old approach or the new one
        if (USE_NEW_APPROACH) {
            // Give the rest of input and output ports their dummy nodes
            for (LPort inPort : inPorts) {
                barycenterAssociates.add(createDummyNode(layeredGraph, inPort, null, dummyNodes));
            }

            for (LPort outPort : outPorts) {
                barycenterAssociates.add(createDummyNode(layeredGraph, null, outPort, dummyNodes));
            }
        } else {
            // Iterate through the lists of input and output ports while both lists still
            // have elements and while output ports are still located right of input ports.
            // While this is true, input and output ports may share the same dummy node
            int inPortsIndex = 0;
            int outPortsIndex = outPorts.size() - 1;

            while (inPortsIndex < inPorts.size() && outPortsIndex >= 0) {
                LPort inPort = inPorts.get(inPortsIndex);
                LPort outPort = outPorts.get(outPortsIndex);

                // If the out port is not right of the in port, they cannot share the same
                // dummy node anymore
                if (ports.indexOf(outPort) < ports.indexOf(inPort)) {
                    break;
                }

                // Otherwise, create a dummy node for them
                barycenterAssociates
                        .add(createDummyNode(layeredGraph, inPort, outPort, dummyNodes));

                inPortsIndex++;
                outPortsIndex--;
            }

            // Give the rest of input and output ports their dummy nodes
            while (inPortsIndex < inPorts.size()) {
                barycenterAssociates.add(createDummyNode(layeredGraph, inPorts.get(inPortsIndex),
                        null, dummyNodes));
                inPortsIndex++;
            }

            while (outPortsIndex >= 0) {
                barycenterAssociates.add(createDummyNode(layeredGraph, null,
                        outPorts.get(outPortsIndex), dummyNodes));
                outPortsIndex--;
            }
        }

        // in / out ports get their own dummy nodes
        for (LPort inOutPort : inOutPorts) {
            barycenterAssociates
                    .add(createDummyNode(layeredGraph, inOutPort, inOutPort, dummyNodes));
        }
    }

    /**
     * Creates a dummy node for the given ports. Edges going into {@code inPort} are rerouted to the
     * dummy node's input port. Edges leaving the {@code outPort} are rerouted to the dummy node's
     * output port. Both arguments may refer to the same port. The dummy's port have their
     * {@code ORIGIN} property set to the port whose edges have been rerouted to them.
     * 
     * @param layeredGraph
     *            the layered graph.
     * @param inPort
     *            the input port whose edges to reroute. May be {@code null}.
     * @param outPort
     *            the output port whose edges to reroute. May be {@code null}.
     * @param dummyNodes
     *            list the created dummy node should be added to.
     * @return the created dummy node.
     */
    private LNode createDummyNode(final LGraph layeredGraph, final LPort inPort,
            final LPort outPort, final List<LNode> dummyNodes) {

        LNode dummy = new LNode(layeredGraph);
        dummy.setType(NodeType.NORTH_SOUTH_PORT);
        dummy.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_POS);

        int crossingHint = 0;

        // Input port
        if (inPort != null) {
            // The port is expected to have edges connected to it
            assert !(inPort.getIncomingEdges().isEmpty() && inPort.getOutgoingEdges().isEmpty());

            LPort dummyInputPort = new LPort();
            dummyInputPort.setProperty(InternalProperties.ORIGIN, inPort);
            dummy.setProperty(InternalProperties.ORIGIN, inPort.getNode());
            dummyInputPort.setSide(PortSide.WEST);
            dummyInputPort.setNode(dummy);

            // Reroute edges
            LEdge[] edgeArray = LGraphUtil.toEdgeArray(inPort.getIncomingEdges());
            for (LEdge edge : edgeArray) {
                edge.setTarget(dummyInputPort);
            }

            // Make sure the inPort knows about the dummy node
            inPort.setProperty(InternalProperties.PORT_DUMMY, dummy);

            crossingHint++;
        }

        // Output port
        if (outPort != null) {
            // The port is expected to have edges connected to it
            assert !(outPort.getIncomingEdges().isEmpty() && outPort.getOutgoingEdges().isEmpty());

            LPort dummyOutputPort = new LPort();
            dummy.setProperty(InternalProperties.ORIGIN, outPort.getNode());
            dummyOutputPort.setProperty(InternalProperties.ORIGIN, outPort);
            dummyOutputPort.setSide(PortSide.EAST);
            dummyOutputPort.setNode(dummy);

            // Reroute edges
            LEdge[] edgeArray = LGraphUtil.toEdgeArray(outPort.getOutgoingEdges());
            for (LEdge edge : edgeArray) {
                edge.setSource(dummyOutputPort);
            }

            // Make sure the outPort knows about the dummy node
            outPort.setProperty(InternalProperties.PORT_DUMMY, dummy);

            crossingHint++;
        }

        // Set the crossing hint used for cross counting later
        dummy.setProperty(InternalProperties.CROSSING_HINT, crossingHint);

        dummyNodes.add(dummy);

        return dummy;
    }

    /**
     * Creates a dummy node for the given non-north-south self-loop edge. The dummy node's
     * {@code ORIGIN} property is set to the edge. The dummy node has two ports, one for each port
     * the node was connected to. Their {@code ORIGIN} property is set accordingly.
     * 
     * @param layeredGraph
     *            the layered graph.
     * @param selfLoop
     *            the self-loop edge.
     * @param dummyNodes
     *            list the created dummy node should be added to.
     */
    private void createSameSideSelfLoopDummyNode(final LGraph layeredGraph, final LEdge selfLoop,
            final List<LNode> dummyNodes) {

        LNode dummy = new LNode(layeredGraph);
        dummy.setType(NodeType.NORTH_SOUTH_PORT);
        dummy.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_POS);
        dummy.setProperty(InternalProperties.ORIGIN, selfLoop);

        // Input port
        LPort dummyInputPort = new LPort();
        dummyInputPort.setProperty(InternalProperties.ORIGIN, selfLoop.getTarget());
        dummyInputPort.setSide(PortSide.WEST);
        dummyInputPort.setNode(dummy);

        // Output port
        LPort dummyOutputPort = new LPort();
        dummyOutputPort.setProperty(InternalProperties.ORIGIN, selfLoop.getSource());
        dummyOutputPort.setSide(PortSide.EAST);
        dummyOutputPort.setNode(dummy);

        // Make sure the ports know about the dummy node
        selfLoop.getSource().setProperty(InternalProperties.PORT_DUMMY, dummy);
        selfLoop.getTarget().setProperty(InternalProperties.PORT_DUMMY, dummy);

        // Disconnect the edge
        selfLoop.setSource(null);
        selfLoop.setTarget(null);

        dummyNodes.add(dummy);

        // Set the crossing hint used for cross counting later
        dummy.setProperty(InternalProperties.CROSSING_HINT, 2);
    }

    /**
     * Creates two dummy nodes for the given north-south self-loop edge. Each dummy node has only
     * one port, on the specified side of the node. Their {@code ORIGIN} property is set
     * accordingly.
     * 
     * @param layeredGraph
     *            the layered graph.
     * @param selfLoop
     *            the self-loop edge.
     * @param northDummyNodes
     *            list the created northern dummy node should be added to.
     * @param southDummyNodes
     *            list the created southern dummy node should be added to.
     * @param portSide
     *            on which sides on the dummy nodes the self-loop ports should be placed.
     */
    private void createNorthSouthSelfLoopDummyNodes(final LGraph layeredGraph,
            final LEdge selfLoop, final List<LNode> northDummyNodes,
            final List<LNode> southDummyNodes, final PortSide portSide) {

        // North dummy
        LNode northDummy = new LNode(layeredGraph);
        northDummy.setType(NodeType.NORTH_SOUTH_PORT);
        northDummy.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_POS);
        northDummy.setProperty(InternalProperties.ORIGIN, selfLoop.getSource().getNode());

        LPort northDummyOutputPort = new LPort();
        northDummyOutputPort.setProperty(InternalProperties.ORIGIN, selfLoop.getSource());
        northDummyOutputPort.setSide(portSide);
        northDummyOutputPort.setNode(northDummy);

        selfLoop.getSource().setProperty(InternalProperties.PORT_DUMMY, northDummy);

        // South dummy
        LNode southDummy = new LNode(layeredGraph);
        southDummy.setType(NodeType.NORTH_SOUTH_PORT);
        southDummy.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_POS);
        southDummy.setProperty(InternalProperties.ORIGIN, selfLoop.getTarget().getNode());

        LPort southDummyInputPort = new LPort();
        southDummyInputPort.setProperty(InternalProperties.ORIGIN, selfLoop.getTarget());
        southDummyInputPort.setSide(portSide);
        southDummyInputPort.setNode(southDummy);

        selfLoop.getTarget().setProperty(InternalProperties.PORT_DUMMY, southDummy);

        // Reroute the edge
        selfLoop.setSource(northDummyOutputPort);
        selfLoop.setTarget(southDummyInputPort);

        northDummyNodes.add(0, northDummy);
        southDummyNodes.add(southDummy);

        // Set the crossing hints used for cross counting later
        northDummy.setProperty(InternalProperties.CROSSING_HINT, 1);
        southDummy.setProperty(InternalProperties.CROSSING_HINT, 1);
    }
}
