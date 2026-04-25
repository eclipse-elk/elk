/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate.wrapping;

import java.util.List;
import java.util.ListIterator;
import java.util.function.Predicate;
import java.util.ArrayList;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LGraphElement;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.intermediate.LongEdgeJoiner;
import org.eclipse.elk.alg.layered.intermediate.wrapping.BreakingPointInserter.BPInfo;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.util.IElkProgressMonitor;
import org.eclipse.elk.alg.layered.graph.LGraphUtil;


/**
 * The second of the triumvirate of breaking point processors. While the other two
 * processors are responsible for inserting and removing the breaking point dummies, 
 * this processors performs the actual 'wrapping' of the graph.
 * 
 * As such it locates breaking point start dummies and moves any subsequent 
 * layers to the previous (already existing) layers. Care has to be taken
 * to rearrange everything in a proper order within the layers 
 * as the crossing minimization phase has already finished. 
 *
 *
 * During the (optional) improvement steps, internally the {@link LGraphElement#id} field of nodes is used 
 * to remember a node's index within its layer. It is important that this index is updated when nodes are moved.
 *
 * <dl>
 *   <dt>Preconditions:</dt>
 *     <dd>crossings minimization has been performed, breaking point dummies were inserted</dd>
 *   <dt>Postconditions:</dt>
 *     <dd>the graph's layering has been adjusted as desired by the inserted breaking points</dd>
 *   <dt>Slots:</dt>
 *     <dd>After phase 3.</dd>
 *   <dt>Same-slot dependencies:</dt>
 *     <dd>none</dd>
 * </dl>
 * 
 * @see BreakingPointInserter
 * @see BreakingPointRemover
 */
public class BreakingPointProcessor implements ILayoutProcessor<LGraph> {

    @Override
    public void process(final LGraph graph, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Breaking Point Processor", 1);
        
        // #1 perform initial wrapping
        performWrapping(graph);
        
        // #2 if desired, improve edge lengths
        if (graph.getProperty(LayeredOptions.WRAPPING_MULTI_EDGE_IMPROVE_WRAPPED_EDGES)) {
            
            // assign indexes to nodes, s.t. the node's index
            // equals the nodes position within its layer
            for (Layer l : graph.getLayers()) {
                int index = 0;
                for (LNode n : l.getNodes()) {
                    n.id = index++;
                }
            }
            
            improveMultiCutIndexEdges(graph);
            
            improveUnneccesarilyLongEdges(graph, true);
            improveUnneccesarilyLongEdges(graph, false);
        }
        
        progressMonitor.done();
    }
    
    /**
     * Main wrapping procedure, placing the determined chunks one below the other. 
     */
    private void performWrapping(final LGraph graph) {
        
        final List<Layer> layers = graph.getLayers();
        final ListIterator<Layer> layerIt = layers.listIterator();
    
        // add initial empty layer to account for break point start dummies 
        layerIt.add(new Layer(graph));
        
        // iterate the layers from left to right
        //  as soon as a layer with a BREAKING POINT start dummy is encountered,
        //  the nodes of the next layer should be moved to the very first layer
        //  consequently 'wrapping' the graph
        boolean reverse = false;
        int idx = 1;
        
        while (layerIt.hasNext()) {
            
            Layer layer = layerIt.next();
            Layer newLayer = layers.get(idx);
            List<LNode> nodesToMove = new ArrayList<>(layer.getNodes());
           
            // remember an offset used for adding in-layer dummies later
            int offset = nodesToMove.size();
            
            // move the nodes to their new layer
            for (LNode n : nodesToMove) {
                n.setLayer(newLayer);
            }
            
            if (reverse) {
                // important to introduce the chains of long edge dummies in reversed order
                for (LNode n : LGraphUtil.reversed(nodesToMove)) {
                    
                    List<LEdge> incEdges = new ArrayList<>();
                    n.getIncomingEdges().forEach(incEdges::add);
                    for (LEdge e : incEdges) {
                        
                        // reverse the edge
                        e.reverse(graph, true);
                        graph.setProperty(InternalProperties.CYCLIC, true);
    
                        // insert proper dummy nodes for the newly created long edge
                        List<LEdge> dummyEdges = CuttingUtils.insertDummies(graph, e, offset);
                        
                        // ameliorate breaking point info
                        BPInfo bpi = n.getProperty(InternalProperties.BREAKING_POINT_INFO);
                        LEdge startInLayerEdge = dummyEdges.get(dummyEdges.size() - 1);
                        bpi.startInLayerDummy = startInLayerEdge.getSource().getNode();
                        bpi.startInLayerEdge = startInLayerEdge;
                        bpi.endInLayerDummy = e.getTarget().getNode();
                        bpi.endInLayerEdge = e;
                    }
                }
                reverse = false;
                
            } else {   
                if (!nodesToMove.isEmpty()) {
                    LNode aNode = nodesToMove.get(0);
                    if (aNode.getType() == NodeType.BREAKING_POINT) {
                        // next layer should be moved (it contains the breaking point end dummies)
                        reverse = true;
                        // start moving nodes to the very first layer
                        idx = -1;
                    }
                }
            }
            
            idx++;
        }
    
        // remove old layers that are now empty
        ListIterator<Layer> it = graph.getLayers().listIterator();
        while (it.hasNext()) {
            Layer l = it.next();
            if (l.getNodes().isEmpty()) {
                it.remove();
            }
        }
    }

    /**
     * Improvement step, reducing the overall edge length of edges 
     * that span multiple cut indexes. 
     */
    private void improveMultiCutIndexEdges(final LGraph graph) {

        for (Layer l : graph.getLayers()) {
            for (LNode n : new ArrayList<>(l.getNodes())) {
                if (BPInfo.isStart(n)) {
                    BPInfo info = n.getProperty(InternalProperties.BREAKING_POINT_INFO);

                    if (info.prev == null && info.next != null) {

                        BPInfo current = info;
                        BPInfo next = info.next;

                        while (next != null) {

                            // drop the dummy chain 
                            dropDummies(next.start, next.startInLayerDummy, false, true);

                            // update the in-layer indexes of subsequent nodes
                            updateIndexesAfter(current.end);
                            updateIndexesAfter(next.start);
                            updateIndexesAfter(next.startInLayerDummy);
                            updateIndexesAfter(next.endInLayerDummy);
                            
                            // reconnect the edge
                            next.endInLayerEdge.setTarget(current.endInLayerEdge.getTarget());
                            current.endInLayerEdge.setTarget(null);

                            // throw out all the unnecessary stuff
                            current.end.setLayer(null);
                            next.start.setLayer(null);
                            next.startInLayerDummy.setLayer(null);
                            next.endInLayerDummy.setLayer(null);

                            // update BPInfo
                            BPInfo newInfo = new BPInfo(current.start, next.end, current.nodeStartEdge,
                                    next.startEndEdge, next.originalEdge);
                            newInfo.startInLayerDummy = current.startInLayerDummy;
                            newInfo.startInLayerEdge = current.startInLayerEdge;
                            newInfo.endInLayerDummy = current.endInLayerDummy;
                            newInfo.endInLayerEdge = next.endInLayerEdge;
                            newInfo.prev = current.prev;
                            newInfo.next = next.next;

                            current.start.setProperty(InternalProperties.BREAKING_POINT_INFO, newInfo);
                            next.end.setProperty(InternalProperties.BREAKING_POINT_INFO, newInfo);

                            next = next.next;
                            // continue with the newly assembled one (_not_ 'next')
                            current = newInfo;
                        }
                    }
                }
            }
        }

    }

    /**
     * Improvement step, reducing the length of edges that were wrapped back.
     */
    private void improveUnneccesarilyLongEdges(final LGraph graph, final boolean forwards) {

        final Predicate<LNode> check = forwards ? BPInfo::isEnd : BPInfo::isStart;

        boolean didsome = false;
        do {
            didsome = false;

            List<Layer> layers = forwards ? LGraphUtil.reversed(graph.getLayers()) : graph.getLayers();
            for (Layer layer : layers) {

                List<LNode> nodes = new ArrayList<>(layer.getNodes());
                if (!forwards) {
                    LGraphUtil.reversed(nodes);
                }

                for (LNode n : nodes) {
                    if (check.test(n)) {
                        LNode bpNode = n;
                        BPInfo bpInfo = n.getProperty(InternalProperties.BREAKING_POINT_INFO);
                        LNode dummy = forwards ? bpInfo.endInLayerDummy : bpInfo.startInLayerDummy;

                        didsome = dropDummies(bpNode, dummy, forwards, false);
                    }

                }
            }
        } while (didsome);

    }

    /**
     * The method drops dummy nodes in scenarios like the following:
     * 
     * <pre>
     *  ... -- d1 -- d2 -- startBP
     *                        |
     *  ... -- d4 -- d3 -- inLayer
     * </pre>
     * 
     * For this case the layers are iterated from right to left and after identifying the {@code startBP} node and the
     * {@code inLayer} dummy, subsequent pairs of long edge dummy nodes {@code (d2, d3)} and {@code (d1, d4)} are
     * removed under certain circumstances.
     * Either the {@code force} flag must be {@code true}, or the dummy node pair must be adjacent in their layer. 
     * 
     * The symmetric case is handled as well if {@code forwards} is true.
     * 
     * @return {@code true} if any dummies have been removed.
     */
    private boolean dropDummies(final LNode bpNode, final LNode inLayerDummy, final boolean forwards,
            final boolean force) {

        // get the adjacent long edge dummies (if any)
        LNode predOne = nextLongEdgeDummy(bpNode, forwards);
        LNode predTwo = nextLongEdgeDummy(inLayerDummy, forwards);

        // predOne and predTwo are guaranteed to be LONG_EDGE 
        boolean didsome = false;
        while (predOne != null && predTwo != null) {

            if (force || isAdjacentOrSeparatedByBreakingpoints(predOne, predTwo, forwards)) {

                LNode nextOne = nextLongEdgeDummy(predOne, forwards);
                LNode nextTwo = nextLongEdgeDummy(predTwo, forwards);
               
                updateIndexesAfter(inLayerDummy);
                updateIndexesAfter(bpNode);
                
                // the two dummies were in same layer
                Layer newLayer = predOne.getLayer();

                LongEdgeJoiner.joinAt(predOne, false);
                LongEdgeJoiner.joinAt(predTwo, false);

                if (forwards) {
                    inLayerDummy.setLayer(predTwo.id, newLayer);
                    inLayerDummy.id = predTwo.id;
                    bpNode.setLayer(predOne.id + 1, newLayer);
                    bpNode.id = predOne.id;
                    
                } else {
                    bpNode.setLayer(predOne.id, newLayer);
                    bpNode.id = predOne.id;
                    inLayerDummy.setLayer(predTwo.id + 1, newLayer);
                    inLayerDummy.id = predTwo.id;
                }

                predOne.setLayer(null);
                predTwo.setLayer(null);

                predOne = nextOne;
                predTwo = nextTwo;

                didsome = true;
            } else {
                break;
            }
        }
        return didsome;
    }

    private boolean isAdjacentOrSeparatedByBreakingpoints(final LNode dummy1, final LNode dummy2,
            final boolean forwards) {

        // d1 is the one in line with the bp dummy
        Layer layer = dummy1.getLayer();
        LNode start = forwards ? dummy2 : dummy1;
        LNode end = forwards ? dummy1 : dummy2;
        
        for (int i = start.id + 1; i < end.id; ++i) {
            // we know that nodes is an array list, thus we can easily query using index
            LNode node = layer.getNodes().get(i);
            if (!(node.getType() == NodeType.BREAKING_POINT || isInLayerDummy(node))) {
                return false;
            }
        }
        
        return true;
    }
    
    private LNode nextLongEdgeDummy(final LNode start, final boolean forwards) {
        
        Iterable<LEdge> edges = forwards ? start.getOutgoingEdges() : start.getIncomingEdges();

        // Usually long edge dummies are connected by exactly one edge. 
        // However, we allow the 'start' node to be an in-layer dummy.
        // Thus, it may have two incoming edges, one of which is the in-layer edge. 
        for (LEdge e : edges) {
            LNode other = e.getOther(start);
            
            if (other.getType() == NodeType.LONG_EDGE
                    && other.getLayer() != start.getLayer()) {
                return other;
            } 
        }
        return null;
    }
    
    private boolean isInLayerDummy(final LNode node) {
        if (node.getType() == NodeType.LONG_EDGE) {
            for (LEdge e : node.getConnectedEdges()) {
                if (!e.isSelfLoop() 
                        && node.getLayer() == e.getOther(node).getLayer()) {
                    return true;
                }
            }
        }
        return false;
    }
    
    private void updateIndexesAfter(final LNode node) {
        for (int i = node.id + 1; i < node.getLayer().getNodes().size(); ++i) {
            node.getLayer().getNodes().get(i).id--;
        }
    }

}
