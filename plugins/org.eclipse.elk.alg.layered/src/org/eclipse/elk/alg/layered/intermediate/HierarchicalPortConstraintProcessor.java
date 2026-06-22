/*******************************************************************************
 * Copyright (c) 2010, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LGraphUtil;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.options.Alignment;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.IElkProgressMonitor;


/**
 * Processes constraints imposed on hierarchical node dummies.
 * 
 * <p>Eastern and western ports cannot be ordered arbitrarily by the crossing minimizer if
 * the port order is fixed. Thus, this processor inserts appropriate in-layer successor
 * constraints to restrict the node ordering.</p>
 * 
 * <p>Northern and southern external ports can sadly not be processed in the usual way with
 * port constraints at least at {@code FIXED_ORDER}. Instead, we need to replace them by new
 * external port dummies. For each node connected to a northern or southern hierarchical port
 * dummy, we need to place a new dummy in the next layer, rerouting the edges appropriately.
 * The original dummies are removed, to be reinserted later by
 * {@link HierarchicalPortOrthogonalEdgeRouter}. For simplification, this is also done in all
 * other port constraint cases. This saves us the trouble of having to differentiate between
 * the different port constraints later on.</p>
 * 
 * <dl>
 *   <dt>Precondition:</dt>
 *     <dd>a layered graph</dd>
 *     <dd>long edge dummies have not yet been inserted</dd>
 *     <dd>layer constraints have not yet been applied.</dd>
 *   <dt>Postcondition:</dt>
 *     <dd>Hierarchical port dummies for northern and southern ports are replaced by multiple dummies if
 *       the port constraints are at least {@code FIXED_ORDER}.</dd>
 *   <dt>Slots:</dt>
 *     <dd>Before phase 3.</dd>
 *   <dt>Same-slot dependencies:</dt>
 *     <dd>{@link LayerConstraintProcessor}</dd>
 * </dl>
 * 
 * @see HierarchicalPortDummySizeProcessor
 * @see HierarchicalPortOrthogonalEdgeRouter
 * @see HierarchicalPortPositionProcessor
 * @author cds
 */
public final class HierarchicalPortConstraintProcessor implements ILayoutProcessor<LGraph> {
    
    /**
     * Comparator to compare nodes by their position values in ascending order. Nodes
     * that are not hierarchical port dummies are sorted to the bottom of a list.
     * 
     * @author cds
     */
    private static class NodeComparator implements Comparator<LNode> {
        @Override
        public int compare(final LNode node1, final LNode node2) {
            NodeType nodeType1 = node1.getType();
            double nodePos1 = node1.getProperty(InternalProperties.PORT_RATIO_OR_POSITION);
            NodeType nodeType2 = node2.getType();
            double nodePos2 = node2.getProperty(InternalProperties.PORT_RATIO_OR_POSITION);
            
            if (nodeType2 != NodeType.EXTERNAL_PORT) {
                return -1;
            } else if (nodeType1 != NodeType.EXTERNAL_PORT) {
                return 1;
            } else {
                if (nodePos1 == nodePos2) {
                    return 0;
                } else {
                    return nodePos1 < nodePos2 ? -1 : 1;
                }
            }
        }
    }
    
    
    // CONSTANTS
    /**
     * Index of the input port in the list of ports of newly created north / south port dummy nodes.
     */
    private static final int DUMMY_INPUT_PORT = 0;

    /**
     * Index of the output port in the list of ports of newly created north / south port dummy nodes.
     */
    private static final int DUMMY_OUTPUT_PORT = 1;
    

    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor monitor) {
        monitor.begin("Hierarchical port constraint processing", 1);
        
        processEasternAndWesternPortDummies(layeredGraph);
        processNorthernAndSouthernPortDummies(layeredGraph);
        
        monitor.done();
    }


    ///////////////////////////////////////////////////////////////////////////////
    // East / West Hierarchical Port Dummies
    
    /**
     * Process eastern and western hierarchical port dummies.
     * 
     * @param layeredGraph the layered graph
     */
    private void processEasternAndWesternPortDummies(final LGraph layeredGraph) {
        // If the port constraints are not at least FIXED_ORDER, there's nothing to be done here
        if (!layeredGraph.getProperty(LayeredOptions.PORT_CONSTRAINTS).isOrderFixed()) {
            return;
        }
        
        List<Layer> layers = layeredGraph.getLayers();
        
        // This affects the first and last layer
        processEasternAndWesternPortDummies(layers.get(0));
        processEasternAndWesternPortDummies(layers.get(layers.size() - 1));
    }
    
    /**
     * Process the eastern and western hierarchical port dummies present in the given layer.
     * 
     * @param layer the layer.
     */
    private void processEasternAndWesternPortDummies(final Layer layer) {
        // Put the nodes into an array
        LNode[] nodes = LGraphUtil.toNodeArray(layer.getNodes());
        
        // Sort the array; hierarchical port dummies are at the top, sorted by
        // position or ratio in ascending order
        Arrays.sort(nodes, new NodeComparator());
        
        // Insert in-layer successor constraints where appropriate
        LNode lastHierarchicalDummy = null;
        
        for (LNode node : nodes) {
            if (node.getType() != NodeType.EXTERNAL_PORT) {
                // No hierarchical port dummy nodes any more
                break;
            }
            
            // Only process dummies created for eastern or western external ports
            PortSide externalPortSide = node.getProperty(InternalProperties.EXT_PORT_SIDE);
            if (externalPortSide != PortSide.WEST && externalPortSide != PortSide.EAST) {
                continue;
            }
            
            if (lastHierarchicalDummy != null) {
                lastHierarchicalDummy.getProperty(InternalProperties.IN_LAYER_SUCCESSOR_CONSTRAINTS).add(node);
            }
            
            lastHierarchicalDummy = node;
        }
    }


    ///////////////////////////////////////////////////////////////////////////////
    // North / South Hierarchical Port Dummies
    
    /**
     * Process northern and southern hierarchical port dummies.
     * 
     * @param layeredGraph the layered graph.
     */
    private void processNorthernAndSouthernPortDummies(final LGraph layeredGraph) {
        // If the port constraints are not at least FIXED_SIDE, there's nothing to do here
        PortConstraints portConstraints = layeredGraph.getProperty(LayeredOptions.PORT_CONSTRAINTS);
        if (!portConstraints.isSideFixed()) {
            return;
        }
        
        List<Layer> layers = layeredGraph.getLayers();
        int layerCount = layers.size();
        
        // For each layer, we keep a map of dummy nodes created for a given original external port
        // dummy. This lets us remember which new dummy node was created for which original external
        // port dummy in which layer. We also keep a list of dummy nodes to be added to each layer
        // since iterating over the new dummy nodes in the map may mess up their order since we don't
        // have hashCode() implemented. Once we do, this code can be reworked to just use the maps
        // again.
        // We keep enough space to hold layerCount + 2 instances of each data structure because we
        // might have to add a new first and last layer, even though that shouldn't normally happen
        List<Map<Object, LNode>> extPortToDummyNodeMap = new ArrayList<>(layerCount + 2);
        List<List<LNode>> newDummyNodes = new ArrayList<>(layerCount + 2);
        
        // Add maps and lists for a new first layer that might have to be created as well as for the
        // current first layer. A map for the next layer is added on each iteration of the for loop
        extPortToDummyNodeMap.add(new HashMap<Object, LNode>());
        extPortToDummyNodeMap.add(new HashMap<Object, LNode>());
        newDummyNodes.add(new ArrayList<LNode>());
        newDummyNodes.add(new ArrayList<LNode>());
        
        // We remember each original external port dummy we encounter (they must be removed from the layers later)
        List<LNode> originalExternalPortDummies = new ArrayList<>();
        
        // Iterate through each layer
        for (int currLayerIdx = 0; currLayerIdx < layerCount; currLayerIdx++) {
            Layer currentLayer = layers.get(currLayerIdx);
            
            // Dummy node maps and lists for the next and previous layer
            Map<Object, LNode> prevExtPortToDummyNodesMap = extPortToDummyNodeMap.get(currLayerIdx);
            Map<Object, LNode> nextExtPortToDummyNodesMap = new HashMap<>();
            extPortToDummyNodeMap.add(nextExtPortToDummyNodesMap);
            
            List<LNode> prevNewDummyNodes = newDummyNodes.get(currLayerIdx);
            List<LNode> nextNewDummyNodes = new ArrayList<>();
            newDummyNodes.add(nextNewDummyNodes);
            
            // Iterate through the layer's nodes, looking for normal nodes connected to
            // northern / southern hierarchical port dummies
            for (LNode currentNode : currentLayer) {
                if (isNorthernOrSouthernDummy(currentNode)) {
                    // It's a northern or southern external port dummy. Schedule for removal and move on to next node
                    originalExternalPortDummies.add(currentNode);
                    continue;
                }
                
                // Iterate over the node's incoming edges
                for (LEdge edge : currentNode.getIncomingEdges()) {
                    LNode sourceNode = edge.getSource().getNode();
                    
                    // Check if it's a northern / southern dummy node
                    if (!isNorthernOrSouthernDummy(sourceNode)) {
                        continue;
                    }
                    
                    // See if a dummy has already been created for the previous layer
                    LNode prevLayerDummy = prevExtPortToDummyNodesMap.get(
                            sourceNode.getProperty(InternalProperties.ORIGIN));
                    
                    if (prevLayerDummy == null) {
                        // No. Create one.
                        prevLayerDummy = createDummy(layeredGraph, sourceNode);
                        prevExtPortToDummyNodesMap.put(
                                sourceNode.getProperty(InternalProperties.ORIGIN), prevLayerDummy);
                        prevNewDummyNodes.add(prevLayerDummy);
                    }
                    
                    // Reroute the edge
                    edge.setSource(prevLayerDummy.getPorts().get(DUMMY_OUTPUT_PORT));
                }
                
                // Iterate over the node's outgoing edges
                for (LEdge edge : currentNode.getOutgoingEdges()) {
                    LNode targetNode = edge.getTarget().getNode();
                    
                    // Check if it's a northern / southern dummy node
                    if (!isNorthernOrSouthernDummy(targetNode)) {
                        continue;
                    }
                    
                    // See if a dummy has already been created for the next layer
                    LNode nextLayerDummy = nextExtPortToDummyNodesMap.get(
                            targetNode.getProperty(InternalProperties.ORIGIN));
                    
                    if (nextLayerDummy == null) {
                        // No. Create one.
                        nextLayerDummy = createDummy(layeredGraph, targetNode);
                        nextExtPortToDummyNodesMap.put(
                                targetNode.getProperty(InternalProperties.ORIGIN), nextLayerDummy);
                        nextNewDummyNodes.add(nextLayerDummy);
                    }
                    
                    // Reroute the edge
                    edge.setTarget(nextLayerDummy.getPorts().get(DUMMY_INPUT_PORT));
                }
            }
        }
        
        // Add the newly created dummy nodes
        for (int i = 0; i < newDummyNodes.size(); i++) {
            List<LNode> nodeList = newDummyNodes.get(i);
            if (nodeList.isEmpty()) {
                // No dummy nodes, so just move on
                continue;
            }
            
            // Find the layer the dummy nodes should be added to
            Layer layer = null;
            if (i == 0) {
                // A new first layer must be created
                layer = new Layer(layeredGraph);
                layers.add(0, layer);
            } else if (i == extPortToDummyNodeMap.size() - 1) {
                // A new layer layer must be created
                layer = new Layer(layeredGraph);
                layers.add(layer);
            } else {
                layer = layers.get(i - 1);
            }
            
            for (LNode dummy : nodeList) {
                dummy.setLayer(layer);
            }
        }
        
        // Iterate through the hierarchical port dummies and remove them
        for (LNode originalDummy : originalExternalPortDummies) {
            // Remove the original dummy; new dummy nodes have already been created for it
            originalDummy.setLayer(null);
        }
        
        // Remember the original external port dummies in the graph
        layeredGraph.setProperty(InternalProperties.EXT_PORT_REPLACED_DUMMIES, originalExternalPortDummies);
    }
    
    /**
     * Checks if the node represents a northern or southern external port.
     * 
     * @param node the node to check.
     * @return {@code true} if the node represents a northern or southern external port,
     *         {@code false} otherwise.
     */
    private boolean isNorthernOrSouthernDummy(final LNode node) {
        NodeType nodeType = node.getType();
        
        if (nodeType == NodeType.EXTERNAL_PORT) {
            PortSide portSide = node.getProperty(InternalProperties.EXT_PORT_SIDE);
            
            return portSide == PortSide.NORTH || portSide == PortSide.SOUTH;
        }
        
        return false;
    }
    
    /**
     * Creates a dummy for the given original dummy. The dummy's {@code ORIGIN} property
     * will point to the original dummy. This way, the original dummy can later be restored,
     * and the newly created dummy can be told apart from the original.
     * 
     * @param layeredGraph the layered graph.
     * @param originalDummy the original dummy.
     * @return the newly created dummy node.
     */
    private LNode createDummy(final LGraph layeredGraph, final LNode originalDummy) {
        LNode newDummy = new LNode(layeredGraph);
        newDummy.copyProperties(originalDummy);
        newDummy.setProperty(InternalProperties.EXT_PORT_REPLACED_DUMMY, originalDummy);
        newDummy.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_POS);
        newDummy.setProperty(LayeredOptions.ALIGNMENT, Alignment.CENTER);
        newDummy.setType(NodeType.EXTERNAL_PORT);
        
        LPort inputPort = new LPort();
        inputPort.setNode(newDummy);
        inputPort.setSide(PortSide.WEST);
        
        LPort outputPort = new LPort();
        outputPort.setNode(newDummy);
        outputPort.setSide(PortSide.EAST);
        
        return newDummy;
    }
    
}
