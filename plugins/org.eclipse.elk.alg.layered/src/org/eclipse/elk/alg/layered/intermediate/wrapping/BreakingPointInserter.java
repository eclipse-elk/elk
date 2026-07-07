/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.layered.intermediate.wrapping;

import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;
import java.util.Set;

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.graph.LPort;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.intermediate.LongEdgeSplitter;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.options.PortConstraints;
import org.eclipse.elk.core.options.PortSide;
import org.eclipse.elk.core.util.IElkProgressMonitor;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

/**
 * Calculates cut points and inserts breaking points into the calculated layering. 
 * The goal is to improve very wide and narrow drawings, i.e. many layers but 
 * few nodes per layer. 
 * 
 * Cutting is to be performed between pairs of layers. Edges that span the cut points are split.  
 * Splitting an edge {@code e=(u,v)} looks like this:
 * <pre>
 *  u --nodeStartEdge--> startMarker --startEndEdge--> endMarker --e--> v
 * </pre>
 * Note that the original edge is pushed "to the right", potentially allowing it to be split it again.
 * 
 * During a later processing step ({@link BreakingPointProcessor}), the {@code startEndEdge} will be 
 * rerouted such that it runs from the very right of the drawing back to the very left. 
 * It gets replaced by a chain of {@link NodeType#LONG_EDGE} dummies.  
 * 
 * <dl>
 *   <dt>Preconditions:</dt>
 *     <dd>a layered graph without {@link NodeType#LONG_EDGE} dummies.</dd>
 *   <dt>Postconditions:</dt>
 *     <dd>{@link NodeType#BREAKING_POINT}s have been introduced.</dd>
 *   <dt>Slots:</dt>
 *     <dd>After phase 2.</dd>
 *   <dt>Same-slot dependencies:</dt>
 *     <dd>Must run before {@link LongEdgeSplitter}</dd>
 * </dl>
 * 
 * @see BreakingPointProcessor
 * @see BreakingPointRemover
 */
public class BreakingPointInserter implements ILayoutProcessor<LGraph> {

    @Override
    public void process(final LGraph graph, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Breaking Point Insertion", 1);

        GraphStats gs = new GraphStats(graph);
        
        // #1 determine the cut points 
        ICutIndexCalculator icic;
        switch (graph.getProperty(LayeredOptions.WRAPPING_CUTTING_STRATEGY)) {
            case MANUAL:
                icic = new ICutIndexCalculator.ManualCutIndexCalculator();
                break;
            case ARD:
                icic = new ARDCutIndexHeuristic();
                break;
            default:
                icic = new MSDCutIndexHeuristic();
                break;
        }
        List<Integer> cuts = icic.getCutIndexes(graph, gs);

        // #2 improve cuts 
        if (graph.getProperty(LayeredOptions.WRAPPING_MULTI_EDGE_IMPROVE_CUTS)) {
            cuts = improveCuts(graph, cuts);
        }
        
        // #3 if they are not guaranteed to be valid, make them valid
        if (!icic.guaranteeValid() && graph.hasProperty(LayeredOptions.WRAPPING_VALIDIFY_STRATEGY)) {
            switch (graph.getProperty(LayeredOptions.WRAPPING_VALIDIFY_STRATEGY)) {
                case LOOK_BACK:
                    cuts = SingleEdgeGraphWrapper.validifyIndexesLookingBack(gs, cuts);
                    break;
                case GREEDY:
                    cuts = SingleEdgeGraphWrapper.validifyIndexesGreedily(gs, cuts);
                    break;
                default:
                    // leave the cuts as they are
            }
        }
        
        if (cuts.isEmpty()) {
            progressMonitor.done();
            return;
        }
        
        // #4 insert the breaking points
        @SuppressWarnings("unused")
        int splits = applyCuts(graph, cuts);
        
        progressMonitor.done();
    }

    private Integer applyCuts(final LGraph graph, final List<Integer> cuts) {

        // initialize iterators
        ListIterator<Layer> layerIt = graph.getLayers().listIterator();
        Iterator<Integer> cutIt = cuts.iterator();
        
        int idx = 0;
        int cut = cutIt.next();
        int noSplitEdges = 0;
        
        Set<LEdge> alreadySplit = Sets.newHashSet();
        
        // keep track of edges that have to be split upon cut insertion
        Set<LEdge> openEdges = Sets.newLinkedHashSet();
        // run
        while (layerIt.hasNext()) {
            
            Layer layer = layerIt.next();
            
            // book keeping of 'open' edges
            for (LNode n : layer) {
                for (LEdge e : n.getOutgoingEdges()) {
                    openEdges.add(e);
                }
                for (LEdge e : n.getIncomingEdges()) {
                    openEdges.remove(e);
                }
            }
            
            // at each cut index, insert two new layers plus breaking point dummies
            if (idx + 1 == cut) {
                
                // insert two new layers
                //  note that since a list iterator is used and 'idx' is decoupled
                //  from the actual position within the list of layers, 
                //  the original 'cut' values are correct and do not have to be offset 
                
                Layer bpLayer1 = new Layer(graph);
                layerIt.add(bpLayer1);
                Layer bpLayer2 = new Layer(graph);
                layerIt.add(bpLayer2);
                
                // splitting an edge 'e' looks like this:
                //      node --nodeStartEdge--> startMarker --startEndEdge--> endMarker --e--> node
                // note that the original edge is pushed "to the right"
                //  that way it is possible to split it again
                
                for (LEdge originalEdge : openEdges) {
                    
                    if (!alreadySplit.contains(originalEdge)) {
                        noSplitEdges++;
                        alreadySplit.add(originalEdge);
                    }
                    
                    // start dummy
                    LNode bpStartMarker = new LNode(graph);
                    bpStartMarker.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
                    bpStartMarker.setLayer(bpLayer1);
                    bpStartMarker.setType(NodeType.BREAKING_POINT);
                    LPort inPortBp1 = new LPort();
                    inPortBp1.setNode(bpStartMarker);
                    inPortBp1.setSide(PortSide.WEST);
                    LPort outPortBp1 = new LPort();
                    outPortBp1.setNode(bpStartMarker);
                    outPortBp1.setSide(PortSide.EAST);
                    
                    // end dummy
                    LNode bpEndMarker = new LNode(graph);
                    bpEndMarker.setProperty(LayeredOptions.PORT_CONSTRAINTS, PortConstraints.FIXED_SIDE);
                    bpEndMarker.setLayer(bpLayer2);
                    bpEndMarker.setType(NodeType.BREAKING_POINT);
                    LPort inPortBp2 = new LPort();
                    inPortBp2.setNode(bpEndMarker);
                    inPortBp2.setSide(PortSide.WEST);
                    LPort outPortBp2 = new LPort();
                    outPortBp2.setNode(bpEndMarker);
                    outPortBp2.setSide(PortSide.EAST);
                    
                    // first dummy edge
                    LEdge nodeStartEdge = new LEdge();
                    nodeStartEdge.setSource(originalEdge.getSource());
                    nodeStartEdge.setTarget(inPortBp1);
                    nodeStartEdge.setProperty(InternalProperties.MODEL_ORDER,
                            originalEdge.getProperty(InternalProperties.MODEL_ORDER));
                    
                    // second dummy edge
                    LEdge startEndEdge = new LEdge();
                    startEndEdge.setSource(outPortBp1);
                    startEndEdge.setTarget(inPortBp2);
                    startEndEdge.setProperty(InternalProperties.MODEL_ORDER,
                            originalEdge.getProperty(InternalProperties.MODEL_ORDER));
                    
                    // reroute the original edge
                    originalEdge.setSource(outPortBp2);
                    
                    // attach information to the dummies
                    BPInfo bpi = new BPInfo(bpStartMarker, bpEndMarker, nodeStartEdge, startEndEdge, originalEdge);
                    bpStartMarker.setProperty(InternalProperties.BREAKING_POINT_INFO, bpi);
                    bpEndMarker.setProperty(InternalProperties.BREAKING_POINT_INFO, bpi);

                    // possibly link a chain of breaking point dummies 
                    LNode prevNode = nodeStartEdge.getSource().getNode();
                    if (prevNode.getType() == NodeType.BREAKING_POINT) {
                        BPInfo bpiPrev = prevNode.getProperty(InternalProperties.BREAKING_POINT_INFO);
                        bpiPrev.next = bpi;
                        bpi.prev = bpiPrev;
                    }
                }
                
                if (cutIt.hasNext()) {
                    cut = cutIt.next();
                } else { 
                    // done
                    break;
                }
            }
            
            idx++;
        }
        
        return noSplitEdges;
    }
    
    @SuppressWarnings("unused")
    private List<Integer> improveCuts(final LGraph graph, final List<Integer> cuts) {
        
        List<Integer> improvedCuts = Lists.newArrayList();
       
        // convert to cut object
        List<Cut> ccuts = Lists.newArrayList();
        Cut lastCut = null;
        for (Integer cutIdx : cuts) {
            Cut cut = new Cut(cutIdx);
            ccuts.add(cut);
            if (lastCut != null) {
                cut.prev = lastCut;
                lastCut.suc = cut;
            }
            lastCut = cut;
        }
        
        // compute edge spans
        int[] spans = computeEdgeSpans(graph);
        
        // iteratively find the best cuts
        for (int i = 0; i < ccuts.size(); ++i) {
        
            Cut lCut = null;
            Cut rCut = ccuts.get(0).selfOrNext();
            
            Cut bestCut = null;
            double bestScore = Double.POSITIVE_INFINITY;
            
            for (int idx = 1; idx < graph.getLayers().size(); ++idx) {
            
                Cut hit;
                int rDist = rCut != null ? Math.abs(rCut.index - idx) : Math.abs(idx - lCut.index) + 1;
                int lDist = lCut != null ? Math.abs(idx - lCut.index) : rDist + 1;
                int dist;
                if (lDist < rDist) {
                    hit = lCut;
                    dist = lDist;
                } else {
                    hit = rCut;
                    dist = rDist;
                }
                
                double score = computeScore(graph, idx, spans[idx], dist);
                if (score < bestScore) {
                    bestScore = score;
                    bestCut = hit;
                    bestCut.newIndex = idx;
                }
                
                if (rCut != null && idx == rCut.index) {
                    lCut = rCut;
                    rCut = rCut.next();
                }
            }
        
            // record
            if (bestCut != null) {
                improvedCuts.add(bestCut.newIndex);
                bestCut.assigned = true;
            
                // propagate offsets as far as possible
                bestCut.offset();
            }
        }
        
        
        Collections.sort(improvedCuts);
        
        return improvedCuts;
    }

    /**
     * spans + dist : weight on dist too small, cuts tend to clump up for highly connected graphs.
     * @param index
     * @param spans
     * @param dist
     * @return
     */
    private double computeScore(final LGraph graph, final int index, final int spans, final int dist) {
        double distancePenalty = graph.getProperty(LayeredOptions.WRAPPING_MULTI_EDGE_DISTANCE_PENALTY);
        return spans + Math.pow(dist, distancePenalty);
    }
    
    /**
     * @param graph
     *            a layered graph.
     * @return an array {@code spans} with the number of edges spanning a certain index, i.e. {@code spans[i]} is the
     *         number of edge spans between layer {@code L_i-1} and {@code L_i}.
     */
    private int[] computeEdgeSpans(final LGraph graph) {
        int[] spans = new int[graph.getLayers().size() + 1];
        
        Set<LEdge> open = Sets.newHashSet();
        int i = 0;
        for (Layer l : graph.getLayers()) {
            
            spans[i++] = open.size();

            for (LNode n : l.getNodes()) {
                for (LEdge e : n.getOutgoingEdges()) {
                    open.add(e);
                }
            }
            
            for (LNode n : l.getNodes()) {
                for (LEdge e : n.getIncomingEdges()) {
                    open.remove(e);
                }
            }
        }
        
        return spans;
    }
    
    /**
     * Holds information about a single breaking point.
     */
    public static class BPInfo {
        
        // SUPPRESS CHECKSTYLE NEXT 16 Javadoc|VisibilityModifier
        public final LNode start;
        public final LNode end;
        
        public final LEdge nodeStartEdge;
        public final LEdge startEndEdge;
        public final LEdge originalEdge;
        
        public LNode startInLayerDummy;
        public LEdge startInLayerEdge;
        public LNode endInLayerDummy;
        public LEdge endInLayerEdge;
        
        public BPInfo prev;
        public BPInfo next;
        
        public BPInfo(final LNode startDummy, final LNode endDummy, 
                final LEdge nodeStartEdge, final LEdge startEndEdge, 
                final LEdge originalEdge) {
            this.start = startDummy;
            this.end = endDummy;
            this.nodeStartEdge = nodeStartEdge;
            this.startEndEdge = startEndEdge;
            this.originalEdge = originalEdge;
        }

        /**
         * @return {@code true} if {@code n} equals the start dummy of this breaking point.
         */
        public static boolean isStart(final LNode n) {
            BPInfo bpi = n.getProperty(InternalProperties.BREAKING_POINT_INFO);
            if (bpi != null) {
                return bpi.start == n;
            }
            return false;
        }

        /**
         * @return {@code true} if {@code n} equals the end dummy of this breaking point.
         */
        public static boolean isEnd(final LNode n) {
            BPInfo bpi = n.getProperty(InternalProperties.BREAKING_POINT_INFO);
            if (bpi != null) {
                return bpi.end == n;
            }
            return false;
        }
        
        @Override
        public String toString() {
            StringBuilder sb = new StringBuilder();
            sb.append("BPInfo[");
            sb.append("\n\tstart=");
            sb.append(start);
            sb.append("\n\tend=");
            sb.append(end);
            sb.append("\n\tnodeStartEdge=");
            sb.append(nodeStartEdge);
            sb.append("\n\tstartEndEdge=");
            sb.append(startEndEdge);
            sb.append("\n\toriginalEdge=");
            sb.append(originalEdge);
            sb.append("\n\tstartInLayerDummy=");
            sb.append(startInLayerDummy);
            sb.append("\n\tstartInLayerEdge=");
            sb.append(startInLayerEdge);
            sb.append("\n\tendInLayerDummy=");
            sb.append(endInLayerDummy);
            sb.append("\n\tendInLayerEdge=");
            sb.append(endInLayerEdge);
            return sb.toString();
        }
    }
    
    /**
     * Internal representation of a 'cut', for convenience.
     */
    private static class Cut {
        
        /** The 'original' cut. */
        private int index;
        
        /** An suggestion (intermediate) for a better cut. */
        private int newIndex;
        
        /** The numerically closest, smaller cut. */
        private Cut prev;
        /** The numerically closest, larger cut. */
        private Cut suc; 
        
        /** Flag indicating whether {@link #newIndex} is final. */
        private boolean assigned = false;
        
        Cut(final int index) {
            this.index = index;
        }
        
        public Cut selfOrNext() {
            if (!assigned) {
                return this;
            } else if (suc != null) {
                return suc.selfOrNext();
            }
            return null;
        }
        
        public Cut next() {
            if (suc != null) {
                return suc.selfOrNext();
            }
            return null;
        }
        
        public void offset() {
            if (!assigned) {
                throw new IllegalStateException("Cannot offset an unassigned cut.");
            }
            int offset = newIndex - index;
            offset(offset);
            offsetPrev(offset);
            offsetSuc(offset);
        }
        
        private void offset(final int offset) {
            index += offset;
        }
        
        private void offsetPrev(final int offset) {
            if (prev != null && !prev.assigned) {
                prev.offset(offset);
                prev.offsetPrev(offset);
            }
        }
        
        private void offsetSuc(final int offset) {
            if (suc != null && !suc.assigned) {
                suc.offset(offset);
                suc.offsetSuc(offset);
            }
        }
        
    }

}
