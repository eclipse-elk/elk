/*******************************************************************************
 * Copyright (c) 2013, 2020 Kiel University and others.
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
package org.eclipse.elk.alg.layered.compound;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LGraphUtil;
import org.eclipse.elk.alg.layered.graph.LLabel;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.options.GraphProperties;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.PortType;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.Direction;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.options.PortLabelPlacement;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.ElkUtil;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.graph.properties.IPropertyHolder;
import org.eclipse.elk.graph.properties.MapPropertyHolder;

/**
 * Preprocess a compound graph by splitting cross-hierarchy edges. The result is stored in
 * {@link LayeredOptions#CROSS_HIERARCHY_MAP}, which is attached to the top-level graph.
 * 
 * <p>This processor assumes that external port dummy nodes only occur on the uppermost level of
 * hierarchy in the input graph. In all deeper levels, it is the job of this processor to create
 * external ports and associated dummy nodes.</p>
 * 
 * <strong>Implementation Notes</strong>
 * 
 * <p>Basically, the algorithm replaces cross-hierarchy edges by hierarchy-local edge segments. It
 * distinguishes between two types of segments: <em>outer segments</em> and <em>inner segments</em>.
 * Outer segments are those two segments that connect to the original source or target port of a
 * hierarchical edge. Inner segments are the remaining segments between the two outer segments.</p>
 * 
 * <p>To split cross-hierarchy edges, the algorithm dives depth-first into the graph's hierarchy tree
 * and begins working its way from the deepest levels of hierarchy upwards. For each contained graph,
 * it looks for cross-hierarchy edges beginning or ending there and starts by splitting those up and
 * thus creating the first outer segments. This will result in external ports being created and added
 * to the graph's parent node. Those are published to the upper level of hierarchy.</p>
 * 
 * </p>All external ports thus published by child nodes are then processed further. With all outer
 * segments created, the algorithm then creates required inner segments. To this end, it continues
 * adding external ports and publishes them to the upper level of hierarchy.</p>
 * 
 * <p>Remember when graph layout was easy? Pepperidge Farm remembers...</p>
 * 
 * <dl>
 *   <dt>Precondition:</dt>
 *     <dd>a compound graph with no layers.</dd>
 *     <dd>no external port dummy nodes except on the uppermost hierarchy level.</dd>
 *   <dt>Postcondition:</dt>
 *     <dd>a compound graph with no layers and no cross-hierarchy edges, but with external ports.</dd>
 * </dl>
 *
 * @author msp
 * @author cds
 */
public class CompoundGraphPreprocessor implements ILayoutProcessor<LGraph> {
    
    ////////////////////////////////////////////////////////////////////////////////////////////
    // Variables
    
    /** map of original edges to generated cross-hierarchy edges. */
    private Map<LEdge, Set<CrossHierarchyEdge>> crossHierarchyMap;
    /** map of ports to their assigned dummy nodes in the nested graphs. */
    private final Map<LPort, LNode> dummyNodeMap = new HashMap<>();

    
    ////////////////////////////////////////////////////////////////////////////////////////////
    // Processing
    
    @Override
    public void process(final LGraph graph, final IElkProgressMonitor monitor) {
        monitor.begin("Compound graph preprocessor", 1);
        
        // This can be a HashMultimap since the values are sorted prior to working with them
        crossHierarchyMap = new HashMap<>();
        
        // create new dummy edges at hierarchy bounds and move the labels around accordingly
        transformHierarchyEdges(graph, null);
        moveLabelsAndRemoveOriginalEdges(graph);

        setSidesOfPortsToSidesOfDummyNodes();

        // Attach cross hierarchy map to the graph and cleanup
        graph.setProperty(InternalProperties.CROSS_HIERARCHY_MAP, crossHierarchyMap);
        crossHierarchyMap = null;
        dummyNodeMap.clear();
        
        monitor.done();
    }
    
    /**
     * Ensures that for each dummy node the external port and vice versa is set. Since the side of the dummy node has
     * also already been calculated, we set this fixed here. Therefore the compound node has a fixed side constraint,
     * with some sides still set to UNDEFINED. This must be dealt with later.
     */
   private void setSidesOfPortsToSidesOfDummyNodes() {
       for (Entry<LPort, LNode> e : dummyNodeMap.entrySet()) {
           LPort externalPort = e.getKey();
           LNode dummyNode = e.getValue();
           dummyNode.setProperty(InternalProperties.ORIGIN, externalPort);
           externalPort.setProperty(InternalProperties.PORT_DUMMY, dummyNode);
           externalPort.setProperty(InternalProperties.INSIDE_CONNECTIONS, true);
           externalPort.setSide(dummyNode.getProperty(InternalProperties.EXT_PORT_SIDE)); 
           dummyNode.getProperty(InternalProperties.EXT_PORT_SIDE);
           externalPort.getNode().setProperty(LayeredOptions.PORT_CONSTRAINTS,
                   PortConstraints.FIXED_SIDE);
           externalPort.getNode().getGraph().getProperty(InternalProperties.GRAPH_PROPERTIES)
                   .add(GraphProperties.NON_FREE_PORTS); 
       }
   }

   /**
     * Recursively transform cross-hierarchy edges into sequences of dummy ports and dummy edges.
     * 
     * @param graph
     *            the layered graph to process
     * @param parentNode
     *            the node that represents the graph in the upper hierarchy level, or {@code null}
     *            if it already is on top-level
     * @return the external ports created to split edges that cross the boundary of the parent node
     */
    private List<ExternalPort> transformHierarchyEdges(final LGraph graph, final LNode parentNode) {
        // process all children and recurse down to gather their external ports
        List<ExternalPort> containedExternalPorts = new ArrayList<>();
        
        for (LNode node : graph.getLayerlessNodes()) {
            LGraph nestedGraph = node.getNestedGraph();
            if (nestedGraph != null) {
                // recursively process the child graph
                List<ExternalPort> childPorts = transformHierarchyEdges(nestedGraph, node);
                containedExternalPorts.addAll(childPorts);
                
                // process inside self loops
                processInsideSelfLoops(nestedGraph, node);
                
                // make sure that all hierarchical ports have had dummy nodes created for them (some
                // will already have been created, but perhaps not all)
                if (nestedGraph.getProperty(InternalProperties.GRAPH_PROPERTIES).contains(
                        GraphProperties.EXTERNAL_PORTS)) {
                    
                    // We need the port constraints to position external ports (Issues KIPRA-1528, ELK-48) 
                    PortConstraints portConstraints = node.getProperty(LayeredOptions.PORT_CONSTRAINTS);
                    boolean insidePortLabels =
                            node.getProperty(LayeredOptions.PORT_LABELS_PLACEMENT).contains(PortLabelPlacement.INSIDE);

                    for (LPort port : node.getPorts()) {
                        // Make sure that every port has a dummy node created for it
                        LNode dummyNode = dummyNodeMap.get(port);
                        if (dummyNode == null) {
                            dummyNode = LGraphUtil.createExternalPortDummy(port,
                                    portConstraints, port.getSide(), -port.getNetFlow(),
                                    null, new KVector(), port.getSize(),
                                    nestedGraph.getProperty(LayeredOptions.DIRECTION), nestedGraph);
                            dummyNode.setProperty(InternalProperties.ORIGIN, port);
                            dummyNodeMap.put(port, dummyNode);
                            nestedGraph.getLayerlessNodes().add(dummyNode);
                        }
                        
                        // We need to reserve space for external port labels by adding those labels to our dummy nodes
                        LPort dummyNodePort = dummyNode.getPorts().get(0);
                        
                        for (LLabel extPortLabel : port.getLabels()) {
                            LLabel dummyPortLabel = new LLabel();
                            dummyPortLabel.getSize().x = extPortLabel.getSize().x;
                            dummyPortLabel.getSize().y = extPortLabel.getSize().y;
                            dummyNodePort.getLabels().add(dummyPortLabel);
                            
                            // If port labels are placed outside, modify the size. At this point, the port's side
                            // may not be known yet if port constraints are free. If they are, however, we know that
                            // the port will end up on either the east or west side. (see #596).
                            // But if the port labels are fixed, we should consider the part that is inside the node.
                            if (!insidePortLabels) {
                                PortSide side = port.getSide();
                                double insidePart = 0;
                                if (PortLabelPlacement
                                        .isFixed(node.getProperty(LayeredOptions.PORT_LABELS_PLACEMENT))) {
                                    // We use 0 as port border offset here, as we only want the label part that is
                                    // inside the node "after" the port.
                                    insidePart = ElkUtil.computeInsidePart(extPortLabel.getPosition(),
                                            extPortLabel.getSize(), port.getSize(), 0, side);
                                }
                                if (portConstraints == PortConstraints.FREE
                                        || PortSide.SIDES_EAST_WEST.contains(side)) {
                                    dummyPortLabel.getSize().x = insidePart;
                                } else {
                                    dummyPortLabel.getSize().y = insidePart;
                                }
                            }
                        }
                    }
                }
            }
        }
        
        // this will be the list of external ports we will export
        List<ExternalPort> exportedExternalPorts = new ArrayList<>();
        
        // process the cross-hierarchy edges connected to the inside of the child nodes
        processInnerHierarchicalEdgeSegments(graph, parentNode, containedExternalPorts,
                exportedExternalPorts);
        
        // process the cross-hierarchy edges connected to the outside of the parent node
        if (parentNode != null) {
            processOuterHierarchicalEdgeSegments(graph, parentNode, exportedExternalPorts);
        }
        
        return exportedExternalPorts;
    }

    /**
     * Moves all labels of the original edges to the appropriate dummy edges and removes the
     * original edges from the graph.
     * 
     * @param graph
     *            the top-level graph.
     */
    private void moveLabelsAndRemoveOriginalEdges(final LGraph graph) {
        // move all labels of the original edges to the appropriate dummy edges and remove the original
        // edges from the graph
        for (LEdge origEdge : crossHierarchyMap.keySet()) {
            // if the original edge had any labels, we need to move them to the newly introduced edge
            // segments
            if (origEdge.getLabels().size() > 0) {
                // retrieve and sort the edge segments introduced for the original edge
                List<CrossHierarchyEdge> edgeSegments = new ArrayList<CrossHierarchyEdge>(
                        crossHierarchyMap.get(origEdge));
                Collections.sort(edgeSegments, new CrossHierarchyEdgeComparator(graph));
                
                // iterate over the labels and move them to the edge segments
                Iterator<LLabel> labelIterator = origEdge.getLabels().listIterator();
                while (labelIterator.hasNext()) {
                    LLabel currLabel = labelIterator.next();
                    
                    // find the index of the dummy edge we will move the label to
                    int targetDummyEdgeIndex = -1;
                    switch (currLabel.getProperty(LayeredOptions.EDGE_LABELS_PLACEMENT)) {
                    case HEAD:
                        targetDummyEdgeIndex = edgeSegments.size() - 1;
                        break;
                    
                    case CENTER:
                        targetDummyEdgeIndex = getShallowestEdgeSegment(edgeSegments);
                        break;
                        
                    case TAIL:
                        targetDummyEdgeIndex = 0;
                        break;
                        
                    default:
                        // we have no idea what to do with the label, so ignore it    
                    }
                    
                    // move the label if we were lucky enough to find a new home for it
                    if (targetDummyEdgeIndex != -1) {
                        CrossHierarchyEdge targetSegment = edgeSegments.get(targetDummyEdgeIndex);
                        targetSegment.getEdge().getLabels().add(currLabel);
                        targetSegment.getEdge().getSource().getNode().getGraph().getProperty(
                                InternalProperties.GRAPH_PROPERTIES).add(GraphProperties.END_LABELS);
                        targetSegment.getEdge().getSource().getNode().getGraph().getProperty(
                                InternalProperties.GRAPH_PROPERTIES).add(GraphProperties.CENTER_LABELS);
                        
                        labelIterator.remove();
                        currLabel.setProperty(InternalProperties.ORIGINAL_LABEL_EDGE, origEdge);
                    }
                }
            }
            
            // remove original edge
            origEdge.setSource(null);
            origEdge.setTarget(null);
        }
    }
    
    /**
     * Determines the index of the shallowest edge segment in a {@link CrossHierarchyEdge}. The
     * shallowest segment is the last edge segement with type set to {@link PortType.OUTPUT}. If no
     * such element exists, the first segment is the shallowest segment.
     *
     * @param edgeSegments
     *            The sorted list of all edge segments in the {@link CrossHierarchyEdge}
     * @return The index of the shallowest edge segment
     */
    private int getShallowestEdgeSegment(final List<CrossHierarchyEdge> edgeSegments) {
        int result = -1;
        int index = 0;

        for (CrossHierarchyEdge crossHierarchyEdge : edgeSegments) {
            if (crossHierarchyEdge.getType().equals(PortType.INPUT)) {
                // We have an edge pointing downwards here.
                // If this is the first segment we take this element.
                // If there was a segement before, we take that one
                result = index == 0 ? 0 : index - 1;
                break;
            } else if (index == edgeSegments.size() - 1) {
                // If this is the last segment and we didn't find a descending edge,
                // This edge is always ascending. Place the label at the last
                // segment here.
                result = index;
            }
            index += 1;
        }

        return result;
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////
    // Inner Hierarchical Edge Segment Processing
    
    /**
     * Deals with the inner segments of hierarchical edges by breaking them between external ports.
     * For each external port contained in child nodes, this method adds appropriate new external
     * ports and / or dummy edges.
     * 
     * @param graph
     *            the graph whose child nodes have exposed external ports.
     * @param parentNode
     *            the node that represents the graph in the upper hierarchy level, or {@code null}
     *            if the graph already is at the top level.
     * @param containedExternalPorts
     *            list of external ports exposed by the graph's child nodes.
     * @param exportedExternalPorts
     *            list that will be filled with the external ports this method creates.
     */
    private void processInnerHierarchicalEdgeSegments(final LGraph graph, final LNode parentNode,
            final List<ExternalPort> containedExternalPorts,
            final List<ExternalPort> exportedExternalPorts) {
        
        // we remember the ports and the dummy nodes we create to add them to the graph afterwards
        // (this is not strictly necessary, but allows us to reuse methods we also use for outer
        // hierarchy edge segments)
        List<ExternalPort> createdExternalPorts = new ArrayList<>();
        
        // iterate over the list of contained external ports
        for (ExternalPort externalPort : containedExternalPorts) {
            ExternalPort currentExternalPort = null;
            
            // we process output ports and input ports a bit differently
            if (externalPort.type == PortType.OUTPUT) {
                // iterate over the port's original edges
                for (LEdge outEdge : externalPort.origEdges) {
                    /* at this point, we distinguish three cases:
                     *  1. The edge connects to a direct child of the parent node.
                     *  2. connects two external ports of direct children of the parent node.
                     *  3. The edge comes from a child node and either connects directly to the parent
                     *     node or leaves it.
                     */
                    LNode targetNode = outEdge.getTarget().getNode();
                    if (targetNode.getGraph() == graph) {
                        // case 1: edge connects to a direct chlid
                        connectChild(graph, externalPort, outEdge, externalPort.dummyPort,
                                outEdge.getTarget());
                    } else if (parentNode == null || LGraphUtil.isDescendant(targetNode, parentNode)) {
                        // case 2: edge connects two direct children
                        connectSiblings(graph, externalPort, containedExternalPorts, outEdge);
                    } else {
                        // case 3: edge connects to parent node or to the outside world
                        ExternalPort newExternalPort = introduceHierarchicalEdgeSegment(
                                graph,
                                parentNode,
                                outEdge,
                                externalPort.dummyPort,
                                PortType.OUTPUT,
                                currentExternalPort);
                        if (newExternalPort != currentExternalPort) {
                            createdExternalPorts.add(newExternalPort);
                        }
                        
                        // the port is our new current external port if it is exported
                        if (newExternalPort.exported) {
                            currentExternalPort = newExternalPort;
                        }
                    }
                }
            } else {
                // iterate over the port's original edges
                for (LEdge inEdge : externalPort.origEdges) {
                    /* at this point, we distinguish three cases:
                     *  1. The edge comes from a direct child of the parent node.
                     *  2. connects two external ports of direct children of the parent node. (we don't
                     *     deal with that case here; the code that handles output ports above does that)
                     *  3. The edge comes from the parent node or from the outside.
                     */
                    LNode sourceNode = inEdge.getSource().getNode();
                    if (sourceNode.getGraph() == graph) {
                        // case 1: edge comes from a direct child
                        connectChild(graph, externalPort, inEdge, inEdge.getSource(),
                                externalPort.dummyPort);
                    } else if (parentNode == null || LGraphUtil.isDescendant(sourceNode, parentNode)) {
                        // case 2: edge connects two direct children; this case is handled in the code
                        //         for output ports above, so there's nothing to do here
                        continue;
                    } else {
                        // case 3: edge comes from the parent node or from the outside
                        ExternalPort newExternalPort = introduceHierarchicalEdgeSegment(
                                graph,
                                parentNode,
                                inEdge,
                                externalPort.dummyPort,
                                PortType.INPUT,
                                currentExternalPort);
                        if (newExternalPort != currentExternalPort) {
                            createdExternalPorts.add(newExternalPort);
                        }
                        
                        // the port is our new current external port if it is exported
                        if (newExternalPort.exported) {
                            currentExternalPort = newExternalPort;
                        }
                    }
                }
            }
        }
        
        // add dummy nodes and exported external ports
        for (ExternalPort externalPort : createdExternalPorts) {
            if (!graph.getLayerlessNodes().contains(externalPort.dummyNode)) {
                graph.getLayerlessNodes().add(externalPort.dummyNode);
            }
            
            if (externalPort.exported) {
                exportedExternalPorts.add(externalPort);
            }
        }
    }

    /**
     * Connects an external port with a child node of the given graph. To this end, a new dummy edge
     * is inserted and associated with the original hierarchy-crossing edge in the cross hierarchy
     * map.
     * 
     * @param graph
     *            the graph whose child to connect.
     * @param externalPort
     *            the external port that provides the other end of the connection.
     * @param origEdge
     *            the original hierarchy-crossing edge.
     * @param sourcePort
     *            the source port the edge shall be connected to.
     * @param targetPort
     *            the target port the edge shall be connected to.
     */
    private void connectChild(final LGraph graph, final ExternalPort externalPort, final LEdge origEdge,
            final LPort sourcePort, final LPort targetPort) {
        
        // add new dummy edge and connect properly
        LEdge dummyEdge = createDummyEdge(graph, origEdge);
        dummyEdge.setSource(sourcePort);
        dummyEdge.setTarget(targetPort);
        
        crossHierarchyMap.computeIfAbsent(origEdge, k -> new HashSet<>())
                .add(new CrossHierarchyEdge(dummyEdge, graph, externalPort.type));
    }

    /**
     * Connects external ports of two child nodes of the given graph. To this end, the provided list
     * of external ports is searched for the counterpart of the provided external output port, and a
     * new dummy edge is created to connect the two. The dummy edge is associated with the original
     * hierarchy-crossing edge in the cross hierarchy map.
     * 
     * @param graph
     *            the graph whose child nodes to connect.
     * @param externalOutputPort
     *            the external output port.
     * @param containedExternalPorts
     *            list of external ports exposed by children of the graph. This list is searched for
     *            the external target port.
     * @param origEdge
     *            the original edge that is being broken.
     */
    private void connectSiblings(final LGraph graph, final ExternalPort externalOutputPort,
            final List<ExternalPort> containedExternalPorts, final LEdge origEdge) {
        
        // find the opposite external port
        ExternalPort targetExternalPort = null;
        for (ExternalPort externalPort2 : containedExternalPorts) {
            if (externalPort2 != externalOutputPort && externalPort2.origEdges.contains(origEdge)) {
                targetExternalPort = externalPort2;
                break;
            }
        }
        assert targetExternalPort.type == PortType.INPUT;
        
        // add new dummy edge and connect properly
        LEdge dummyEdge = createDummyEdge(graph, origEdge);
        dummyEdge.setSource(externalOutputPort.dummyPort);
        dummyEdge.setTarget(targetExternalPort.dummyPort);
        
        crossHierarchyMap.computeIfAbsent(origEdge, k -> new HashSet<>())
                .add(new CrossHierarchyEdge(dummyEdge, graph, externalOutputPort.type));
    }
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////
    // Outer Hierarchical Edge Segment Processing
    
    /**
     * Deals with the outer segments of hierarchical edges by breaking them at their source or
     * target. For each hierarchical edge that starts or ends at one of the graph's children, this
     * method adds appropriate new external ports and / or dummy edges to the graph.
     * 
     * @param graph
     *            the graph whose child nodes have exposed external ports.
     * @param parentNode
     *            the node that represents the graph in the upper hierarchy level, or {@code null}
     *            if the graph already is at the top level.
     * @param exportedExternalPorts
     *            list that will be filled with the external ports this method creates.
     */
    private void processOuterHierarchicalEdgeSegments(final LGraph graph, final LNode parentNode,
            final List<ExternalPort> exportedExternalPorts) {
        
        // we need to remember the ports and the dummy nodes we create to add them to the graph
        // afterwards (to avoid concurrent modification exceptions)
        List<ExternalPort> createdExternalPorts = new ArrayList<>();
        
        // iterate over all ports of the graph's child nodes
        for (LNode childNode : graph.getLayerlessNodes()) {
            for (LPort childPort : childNode.getPorts()) {
                // we treat outgoing and incoming edges separately
                ExternalPort currentExternalOutputPort = null;
                for (LEdge outEdge : LGraphUtil.toEdgeArray(childPort.getOutgoingEdges())) {
                    if (!LGraphUtil.isDescendant(outEdge.getTarget().getNode(), parentNode)) {
                        // the edge goes to the outside or to the parent node itself, so create an
                        // external port if necessary and introduce a new dummy edge
                        ExternalPort newExternalPort = introduceHierarchicalEdgeSegment(
                                graph,
                                parentNode,
                                outEdge,
                                outEdge.getSource(),
                                PortType.OUTPUT,
                                currentExternalOutputPort);
                        if (newExternalPort != currentExternalOutputPort) {
                            createdExternalPorts.add(newExternalPort);
                        }
                        
                        // the port is our new current external port if it is exported
                        if (newExternalPort.exported) {
                            currentExternalOutputPort = newExternalPort;
                        }
                    }
                }

                ExternalPort currentExternalInputPort = null;
                for (LEdge inEdge : LGraphUtil.toEdgeArray(childPort.getIncomingEdges())) {
                    if (!LGraphUtil.isDescendant(inEdge.getSource().getNode(), parentNode)) {
                        // the edge comes from the outside or from the parent node itself, so create an
                        // external port if necessary and introduce a new dummy edge
                        ExternalPort newExternalPort = introduceHierarchicalEdgeSegment(
                                graph,
                                parentNode,
                                inEdge,
                                inEdge.getTarget(),
                                PortType.INPUT,
                                currentExternalInputPort);
                        if (newExternalPort != currentExternalInputPort) {
                            createdExternalPorts.add(newExternalPort);
                        }
                        
                        // the port is our new current external port if it is exported
                        if (newExternalPort.exported) {
                            currentExternalInputPort = newExternalPort;
                        }
                    }
                }
            }
        }
        
        // add dummy nodes and exported external ports
        for (ExternalPort externalPort : createdExternalPorts) {
            if (!graph.getLayerlessNodes().contains(externalPort.dummyNode)) {
                graph.getLayerlessNodes().add(externalPort.dummyNode);
            }
            
            if (externalPort.exported) {
                exportedExternalPorts.add(externalPort);
            }
        }
    }
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////
    // Inside Self Loop Processing
    
    private void processInsideSelfLoops(final LGraph nestedGraph, final LNode node) {
        // Check if inside self loops are enabled for the node
        if (!node.getProperty(LayeredOptions.INSIDE_SELF_LOOPS_ACTIVATE)) {
            return;
        }
        
        // Iterate over the edges and look for an inside self loop
        for (LPort lport : node.getPorts()) {
            // Avoid ConcurrentModificationExceptions
            LEdge[] outEdges = LGraphUtil.toEdgeArray(lport.getOutgoingEdges());
            
            for (LEdge outEdge : outEdges) {
                boolean isSelfLoop = outEdge.getTarget().getNode() == node;
                boolean isInsideSelfLoop = isSelfLoop
                        && outEdge.getProperty(LayeredOptions.INSIDE_SELF_LOOPS_YO);
                
                if (isInsideSelfLoop) {
                    // Check if the ports have already been transformed into external port dummies
                    LPort sourcePort = outEdge.getSource();
                    LNode sourceExtPortDummy = dummyNodeMap.get(sourcePort);
                    if (sourceExtPortDummy == null) {
                        sourceExtPortDummy = LGraphUtil.createExternalPortDummy(
                                sourcePort,
                                PortConstraints.FREE,
                                sourcePort.getSide(),
                                -1,
                                null,
                                null,
                                sourcePort.getSize(),
                                nestedGraph.getProperty(LayeredOptions.DIRECTION),
                                nestedGraph);
                        sourceExtPortDummy.setProperty(InternalProperties.ORIGIN, sourcePort);
                        dummyNodeMap.put(sourcePort, sourceExtPortDummy);
                        nestedGraph.getLayerlessNodes().add(sourceExtPortDummy);
                    }
                    
                    LPort targetPort = outEdge.getTarget();
                    LNode targetExtPortDummy = dummyNodeMap.get(targetPort);
                    if (targetExtPortDummy == null) {
                        targetExtPortDummy = LGraphUtil.createExternalPortDummy(
                                targetPort,
                                PortConstraints.FREE,
                                targetPort.getSide(),
                                1,
                                null,
                                null,
                                targetPort.getSize(),
                                nestedGraph.getProperty(LayeredOptions.DIRECTION),
                                nestedGraph);
                        targetExtPortDummy.setProperty(InternalProperties.ORIGIN, targetPort);
                        dummyNodeMap.put(targetPort, targetExtPortDummy);
                        nestedGraph.getLayerlessNodes().add(targetExtPortDummy);
                    }
                    
                    // Create a new dummy edge
                    LEdge dummyEdge = createDummyEdge(nestedGraph, outEdge);
                    dummyEdge.setSource(sourceExtPortDummy.getPorts().get(0));
                    dummyEdge.setTarget(targetExtPortDummy.getPorts().get(0));
                    
                    // Remember the new edge
                    crossHierarchyMap.computeIfAbsent(outEdge, k -> new HashSet<>())
                            .add(new CrossHierarchyEdge(dummyEdge, nestedGraph, PortType.OUTPUT));
                    
                    nestedGraph.getProperty(InternalProperties.GRAPH_PROPERTIES).add(
                            GraphProperties.EXTERNAL_PORTS);
                }
            }
        }
    }

    
    ////////////////////////////////////////////////////////////////////////////////////////////
    // General Hierarchical Edge Segment Processing
    
    /**
     * Does the actual work of creating a new hierarchical edge segment between an external port and
     * a given opposite port. The external port used for the segment is returned. This method does
     * not put any created edges into the cross hierarchy map!
     * 
     * <p>
     * The method first decides on an external port to use for the segment. If the default external
     * port passed to the method is not {@code null} and if external ports are to be merged in the
     * current graph, the default external port is reused. An exception are segments that start or
     * end in the parent node; each such segments gets its own external port.
     * </p>
     * 
     * <p>
     * If a new external port is created, the method also creates a dummy node for it as well as an
     * actual port on the parent node, if no such port already exists, as well as a dummy edge for
     * the connection. Thus, the newly created external port has everything it needs to be properly
     * represented and initialized.
     * </p>
     * 
     * <p>
     * The original edge is added to the list of original edges in the external port used for the
     * segment. The dummy edge is associated with the original hierarchy-crossing edge in the cross
     * hierarchy map.
     * </p>
     * 
     * @param graph
     *            the layered graph.
     * @param parentNode
     *            the graph's parent node, or {@code null} if the graph is at the top level.
     * @param origEdge
     *            the hierarchy-crossing edge that is being broken.
     * @param oppositePort
     *            the port that will be one of the two end points of the new segment.
     * @param portType
     *            the type of the port to create as one of the segment's edge points.
     * @param defaultExternalPort
     *            a default external port we can reuse if external ports should be merged. If this
     *            is {@code null}, a new external port is always created. If this port is reused, it
     *            is returned by this method.
     * @return the (created or reused) external port used as one endpoint of the edge segment.
     */
    private ExternalPort introduceHierarchicalEdgeSegment(final LGraph graph, final LNode parentNode,
            final LEdge origEdge, final LPort oppositePort, final PortType portType,
            final ExternalPort defaultExternalPort) {
        
        // check if external ports are to be merged
        boolean mergeExternalPorts = graph.getProperty(LayeredOptions.MERGE_HIERARCHY_EDGES);
        
        // check if the edge connects to the parent node instead of to the outside world; if so, the
        // parentEndPort will be non-null
        LPort parentEndPort = null;
        if (portType == PortType.INPUT && origEdge.getSource().getNode() == parentNode) {
            parentEndPort = origEdge.getSource();
        } else if (portType == PortType.OUTPUT && origEdge.getTarget().getNode() == parentNode) {
            parentEndPort = origEdge.getTarget();
        }
        
        // only create a new external port if the current one is null or if ports are not to be merged
        // or if the connection actually ends at the parent node
        ExternalPort externalPort = defaultExternalPort;
        if (externalPort == null || !mergeExternalPorts || parentEndPort != null) {
            // create a dummy node that will represent the external port
            PortSide externalPortSide = PortSide.UNDEFINED;
            if (parentEndPort != null) {
                externalPortSide = parentEndPort.getSide();
            } else {
                // We try to infer the port side from the port type if its node has its port constraints
                // set to at least FIXED_SIDE; this may produce strange effects, so the safest thing is
                // for people to set compound node port constraints to FREE
                if (parentNode.getProperty(LayeredOptions.PORT_CONSTRAINTS).isSideFixed()) {
                    externalPortSide = portType == PortType.INPUT ? PortSide.WEST : PortSide.EAST;
                }
            }
            LNode dummyNode = createExternalPortDummy(
                    graph, parentNode, portType, externalPortSide, origEdge);
            
            // create a dummy edge to be connected to the port
            LEdge dummyEdge = createDummyEdge(parentNode.getGraph(), origEdge);

            if (portType == PortType.INPUT) {
                // if the external port is an input port, the source of the edge must be connected to
                // the new dummy node
                dummyEdge.setSource(dummyNode.getPorts().get(0));
                dummyEdge.setTarget(oppositePort);
            } else {
                // if the external port is an output port, the target of the edge must be connected to
                // the new dummy node
                dummyEdge.setSource(oppositePort);
                dummyEdge.setTarget(dummyNode.getPorts().get(0));
            }
            
            // create the external port (the port is to be exported if the connection is not just to the
            // parent node)
            externalPort = new ExternalPort(origEdge, dummyEdge, dummyNode,
                    (LPort) dummyNode.getProperty(InternalProperties.ORIGIN), portType,
                    parentEndPort == null);
        } else {
            // we use an existing external port, so simply add the original edge to its list of
            // original edges
            externalPort.origEdges.add(origEdge);
            
            // merge the properties of the original edges
            double thickness = Math.max(externalPort.newEdge.getProperty(LayeredOptions.EDGE_THICKNESS),
                    origEdge.getProperty(LayeredOptions.EDGE_THICKNESS));
            externalPort.newEdge.setProperty(LayeredOptions.EDGE_THICKNESS, thickness);
        }

        crossHierarchyMap.computeIfAbsent(origEdge, k -> new HashSet<>())
                .add(new CrossHierarchyEdge(externalPort.newEdge, graph, portType));
        
        return externalPort;
    }
    
    /**
     * Creates and initializes a new dummy edge for the given original hierarchy-crossing edge. All
     * that remains to be done afterwards is to properly connect the edge. Nice!
     * 
     * @param graph
     *            the graph the edge will be placed in.
     * @param origEdge
     *            the original hierarchy-crossing edge.
     * @return a new dummy edge.
     */
    private LEdge createDummyEdge(final LGraph graph, final LEdge origEdge) {
        LEdge dummyEdge = new LEdge();
        dummyEdge.copyProperties(origEdge);
        dummyEdge.setProperty(LayeredOptions.JUNCTION_POINTS, null);
        return dummyEdge;
    }
    
    /**
     * Adapted from {@code ElkGraphImporter}.
     * Count how many edges want the port to be an output port of the parent and how many want it to
     * be an input port. An edge coming into the port from the inside votes for the port to be an
     * output port of the parent, as does an edge leaving the port for the outside. The result returned
     * by this method is the so-called net flow as fed into {@code createExternalPortDummy(..)}.
     * 
     * @return the port's net flow.
     */
    private int calculateNetFlow(final LPort port) {
        final LNode node = port.getNode();
        final boolean insideSelfLoopsEnabled = node.getProperty(LayeredOptions.INSIDE_SELF_LOOPS_ACTIVATE);

        int outputPortVote = 0, inputPortVote = 0;
        
        // Iterate over outgoing edges
        for (LEdge outgoingEdge : port.getOutgoingEdges()) {
            final boolean isSelfLoop = outgoingEdge.isSelfLoop();
            final boolean isInsideSelfLoop = isSelfLoop && insideSelfLoopsEnabled
                    && outgoingEdge.getProperty(LayeredOptions.INSIDE_SELF_LOOPS_YO);
            final LNode targetNode = outgoingEdge.getTarget().getNode();

            if (isSelfLoop && isInsideSelfLoop) {
                inputPortVote++;
            } else if (isSelfLoop && !isInsideSelfLoop) {
                outputPortVote++;
            } else if (targetNode.getGraph().getParentNode()  == node) {
                inputPortVote++;
            } else {
                outputPortVote++;
            }
        }
        
        // Iterate over incoming edges
        for (LEdge incomingEdge : port.getIncomingEdges()) {
            final boolean isSelfLoop = incomingEdge.isSelfLoop();
            final boolean isInsideSelfLoop = isSelfLoop && insideSelfLoopsEnabled
                    && incomingEdge.getProperty(LayeredOptions.INSIDE_SELF_LOOPS_YO);
            final LNode sourceNode = incomingEdge.getSource().getNode();

            if (isSelfLoop && isInsideSelfLoop) {
                outputPortVote++;
            } else if (isSelfLoop && !isInsideSelfLoop) {
                inputPortVote++;
            } else if (sourceNode.getGraph().getParentNode() == node) {
                outputPortVote++;
            } else {
                inputPortVote++;
            }
        }
        
        return outputPortVote - inputPortVote;
    }
    
    /**
     * Retrieves a dummy node to be used to represent a new external port of the parent node and to
     * connect a new segment of the given hierarchical edge to. A proper dummy node might already
     * have been created; if so, that one is returned.
     * 
     * @param graph
     *            the graph.
     * @param parentNode
     *            the graph's parent node.
     * @param portType
     *            the type of the new external port.
     * @param edge
     *            the edge that will be connected to the external port.
     * @return an appropriate external port dummy.
     */
    private LNode createExternalPortDummy(final LGraph graph, final LNode parentNode,
            final PortType portType, final PortSide portSide, final LEdge edge) {
        
        LNode dummyNode = null;
        
        // find the port on the outside of its parent node that the edge connects to
        LPort outsidePort = portType == PortType.INPUT ? edge.getSource() : edge.getTarget();
        Direction layoutDirection = LGraphUtil.getDirection(graph);
        
        // check if the edge connects to the parent node or to something way outside...
        if (outsidePort.getNode() == parentNode) {
            // we need to check if a dummy node has already been created for the port
            dummyNode = dummyNodeMap.get(outsidePort);
            if (dummyNode == null) {
                // Ticket #160 explains why we need to pass on a position here. While that works fine to determine port
                // orders, the exact port positions might yet need to be adjusted for some inside paddings
                dummyNode = LGraphUtil.createExternalPortDummy(
                        outsidePort,
                        parentNode.getProperty(LayeredOptions.PORT_CONSTRAINTS),
                        portSide,
                        calculateNetFlow(outsidePort),
                        null,
                        outsidePort.getPosition(),
                        outsidePort.getSize(),
                        layoutDirection,
                        graph
                );
                dummyNode.setProperty(InternalProperties.ORIGIN, outsidePort);
                dummyNodeMap.put(outsidePort, dummyNode);
            }
        } else {
            // We create a new dummy node in any case, and since there is no port yet we have to create one as well. We
            // do need to pass a position vector to keep an NPE from being thrown (#160), but it is questionable
            // whether this part of the code should actually be used with anything other than fixed port constraints.
            dummyNode = LGraphUtil.createExternalPortDummy(
                    createExternalPortProperties(graph),
                    parentNode.getProperty(LayeredOptions.PORT_CONSTRAINTS),
                    portSide,
                    portType == PortType.INPUT ? -1 : 1,
                    null,
                    new KVector(),
                    new KVector(0, 0), // A dummy port is not considered to have a size (#766)
                    layoutDirection,
                    graph
            );
            LPort dummyPort = createPortForDummy(dummyNode, parentNode, portType);
            dummyNode.setProperty(InternalProperties.ORIGIN, dummyPort);
            dummyNodeMap.put(dummyPort, dummyNode);
        }
        
        // set a few graph properties
        graph.getProperty(InternalProperties.GRAPH_PROPERTIES).add(GraphProperties.EXTERNAL_PORTS);
        if (graph.getProperty(LayeredOptions.PORT_CONSTRAINTS).isSideFixed()) {
            graph.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
        } else {
            graph.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FREE);
        }
        
        return dummyNode;
    }
    
    /**
     * Create suitable port properties for dummy external ports.
     * 
     * @param graph
     *            the graph for which the dummy external port is created
     * @return properties to apply to the dummy port
     */
    private static IPropertyHolder createExternalPortProperties(final LGraph graph) {
        IPropertyHolder propertyHolder = new MapPropertyHolder();
        // FIXME No idea why ...
        double offset = graph.getProperty(LayeredOptions.SPACING_EDGE_EDGE) / 2;
        propertyHolder.setProperty(LayeredOptions.PORT_BORDER_OFFSET, offset);
        return propertyHolder;
    }
    
    /**
     * Create a port for an existing external port dummy node.
     * 
     * @param dummyNode
     *            the dummy node
     * @param parentNode
     *            the parent node to which it is attached
     * @param type
     *            the port type
     * @return a new port
     */
    private LPort createPortForDummy(final LNode dummyNode, final LNode parentNode,
            final PortType type) {
        
        LGraph graph = parentNode.getGraph();
        Direction layoutDirection = LGraphUtil.getDirection(graph);
        LPort port = new LPort();
        port.setNode(parentNode);
        switch (type) {
        case INPUT:
            port.setSide(PortSide.fromDirection(layoutDirection).opposed());
            break;
        case OUTPUT:
            port.setSide(PortSide.fromDirection(layoutDirection));
            break;
        }
        port.setProperty(LayeredOptions.PORT_BORDER_OFFSET, dummyNode.getProperty(LayeredOptions.PORT_BORDER_OFFSET));
        return port;
    }
    
    
    ////////////////////////////////////////////////////////////////////////////////////////////
    // Class ExternalPort

    /**
     * An internal representation for external ports. This class is used to pass information
     * gathered on one hierarchy level to the containing hierarchy level. Instances are created
     * whenever a cross-hierarchy edge crosses the hierarchy bounds of a parent node; the instance
     * represents the split point of the edge.
     */
    private static class ExternalPort {
        /** the list of original edges for which the port is created. */
        private List<LEdge> origEdges = new ArrayList<>();
        /** the new edge by which the original edge is replaced. */
        private LEdge newEdge;
        /** the dummy node used by the algorithm as representative for the external port. */
        private LNode dummyNode;
        /** the dummy port used by the algorithm as representative for the external port. */
        private LPort dummyPort;
        /** the flow direction: input or output. */
        private PortType type = PortType.UNDEFINED;
        /**
         * whether the external port will be exported to the outside or not. (it will not be exported
         * if the port was introduced for connections from an inside node to its parent)
         */
        private boolean exported;
        
        /**
         * Create an external port.
         * 
         * @param origEdge
         *            the original edge for which the port is created
         * @param newEdge
         *            the new edge by which the original edge is replaced
         * @param dummyNode
         *            the dummy node used by the algorithm as representative for the external port
         * @param portType
         *            the flow direction: input or output
         * @param exported
         *            whether the external port is to be exported by its parent.
         */
        ExternalPort(final LEdge origEdge, final LEdge newEdge, final LNode dummyNode,
                final LPort dummyPort, final PortType portType, final boolean exported) {
            
            this.origEdges.add(origEdge);
            this.newEdge = newEdge;
            this.dummyNode = dummyNode;
            this.dummyPort = dummyPort;
            this.type = portType;
            this.exported = exported;
        }
    }

}
