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

import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.Layer;
import org.eclipse.elk.alg.layered.intermediate.wrapping.ICutIndexCalculator.ManualCutIndexCalculator;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.util.IElkProgressMonitor;

import com.google.common.collect.Lists;

/**
 * Intermediate processor that splits path-like graphs into 
 * multiple sub-paths to improve the overall aspect ratio of the drawing.
 * 
 * <p>
 * This is a rather problem-specific intermediate processor 
 * that works only on a very small subset of graphs, namely 
 * paths that have very few branches.
 * Also, it only works in conjunction with the 
 * {@link org.eclipse.elk.alg.layered.p2layers.LongestPathLayerer LongestPathLayerer}.
 * </p>
 * 
 * <p>
 * What we do is the following. After layering and crossing minimization
 * we determine into how many "rows" we want to split the path.
 * Let n denote the number of nodes of each row. We start moving 
 * the node in layer n+1 into layer 1, node n+2 into layer 2 and so on.
 * While doing this we have to assure that the order within each 
 * layer is correct, have to insert dummy nodes for the edge that now 
 * spans from layer n+1 to layer 1, and have to treat inverted ports
 * specifically.
 * </p>
 * 
 * <dl>
 *   <dt>Precondition:</dt><dd>A layering has been calculated by the longest
 *   path layerer and the crossing minimization has finished.</dd>
 *   <dt>Postcondition:</dt><dd>Long paths are split such that the 
 *   aspect ratio is closer to the currently available area.</dd>
 *   <dt>Slots:</dt><dd>Before phase 4.</dd>
 * </dl>
 */
public class SingleEdgeGraphWrapper implements ILayoutProcessor<LGraph> {

    @Override
    public void process(final LGraph graph, final IElkProgressMonitor progressMonitor) {
        progressMonitor.begin("Path-Like Graph Wrapping", 1);

        if (graph.getLayers().isEmpty()) {
            progressMonitor.done();
            return;
        }
        
        // initialization 
        GraphStats gs = new GraphStats(graph);
        
        // stop early if there's nothing we can do
        double sumWidth = gs.getMaxWidth() * gs.longestPath;
        double currentAR = sumWidth / gs.getMaxWidth();
        if (gs.dar > currentAR) {
            // nothing to do for us
            progressMonitor.done();
            return;
        }

        // get indexes for cutting
        ICutIndexCalculator icic;
        switch (graph.getProperty(LayeredOptions.WRAPPING_CUTTING_STRATEGY)) {
            case MANUAL: 
                icic = new ManualCutIndexCalculator();
                break;
            case ARD:
                icic = new ARDCutIndexHeuristic();
                break;
            default:
                icic = new MSDCutIndexHeuristic();
        }
        List<Integer> cuts = icic.getCutIndexes(graph, gs);
        
        // if they are not guaranteed to be valid, make them valid
        if (!icic.guaranteeValid()) {
            switch (graph.getProperty(LayeredOptions.WRAPPING_VALIDIFY_STRATEGY)) {
                case LOOK_BACK:
                    cuts = validifyIndexesLookingBack(gs, cuts);
                    break;
                case GREEDY:
                    cuts = validifyIndexesGreedily(gs, cuts);
                    break;
                default:
                    // leave the cuts as they are
            }
        }


        // re-build the layering
        performCuts(graph, gs, cuts);

        progressMonitor.done();
    }
    
    /**
     * Validify {@code desiredCuts} by increasing every forbidden index (and its successors) until they are valid.
     */
    public static List<Integer> validifyIndexesGreedily(final GraphStats gs, final List<Integer> cuts) {
        List<Integer> validCuts = Lists.newArrayList();
        int offset = 0;

        Iterator<Integer> cutIt = cuts.iterator();
        while (cutIt.hasNext()) {
            Integer cut = cutIt.next() + offset;
            while (cut < gs.longestPath && !gs.isCutAllowed(cut)) {
                cut++;
                offset++;
            }
            if (cut >= gs.longestPath) {
                // cut cannot be performed anymore, discard it
                break;
            }
            validCuts.add(cut);
        }
        
        return validCuts;
    }
    
    /**
     * Validify {@code desiredCuts} by finding the closest valid cuts. 
     */
    public static List<Integer> validifyIndexesLookingBack(final GraphStats gs, final List<Integer> desiredCuts) {
        if (desiredCuts.isEmpty()) {
            return Collections.emptyList();
        }
        
        // initialize array with valid cuts
        List<Integer> validCuts = Lists.newArrayList();
        validCuts.add(Integer.MIN_VALUE); // -inf
        for (int i = 1; i < gs.longestPath; ++i) {
            if (gs.isCutAllowed(i)) {
                validCuts.add(i);
            }
        }
        if (validCuts.size() == 1) {
            return Collections.emptyList();
        }
        validCuts.add(Integer.MAX_VALUE); // +inf
        
        return validifyIndexesLookingBack(desiredCuts, validCuts);
    }
    
    /**
     * Validify the {@code desiredCuts} s.t. the resulting cuts are all out of {@code validCuts}. {@code validCuts}
     * is expected to be a list that contains -inf as initial element and +inf as last element.
     */
    private static List<Integer> validifyIndexesLookingBack(final List<Integer> desiredCuts,
            final List<Integer> validCuts) {
        assert validCuts.get(0) == Integer.MIN_VALUE;
        assert validCuts.get(validCuts.size() - 1) == Integer.MAX_VALUE;

        // turn cuts into valid cuts
        List<Integer> finalCuts = Lists.newArrayList();
        int iIdx = 0;
        int cIdx = 0;
        int offset = 0;

        while (iIdx < validCuts.size() - 1 && cIdx < desiredCuts.size()) {
            int current = desiredCuts.get(cIdx) + offset;

            while (validCuts.get(iIdx + 1) < current) {
                iIdx++;
            }

            assert validCuts.get(iIdx) <= current;
            assert current <= validCuts.get(iIdx + 1);

            int select = 0;
            int distLower = current - validCuts.get(iIdx);
            int distHigher = validCuts.get(iIdx + 1) - current;
            if (distLower > distHigher) {
                select++;
            }

            finalCuts.add(validCuts.get(iIdx + select));
            offset += validCuts.get(iIdx + select) - current;
            cIdx++;
            while (cIdx < desiredCuts.size() && desiredCuts.get(cIdx) + offset <= validCuts.get(iIdx + select)) {
                cIdx++;
            }
            iIdx += 1 + select;
        }

        return finalCuts;
    }

    /**
     * Preconditions: <ul>
     *  <li>all cuts are allowed
     *  <li>no cut larger than number of layers
     *  <li>cuts are ordered increasingly.
     * </ul>
     */
    private void performCuts(final LGraph graph, final GraphStats gs, final List<Integer> cuts) {
        
        if (cuts.isEmpty()) {
            return;
        }
        
        // index in the "old" layering
        int index = 0;
        // index in the "new" layering
        int newIndex = 0;
        
        Iterator<Integer> cutIt = cuts.iterator();
        int nextCut = cutIt.next(); // guaranteed to exist
        
        // iterate all original layers
        while (index < gs.longestPath) {

            if (index == nextCut) {
                // start new "row"
                newIndex = 0;
                if (cutIt.hasNext()) {
                    nextCut = cutIt.next();
                } else {
                    // no further cut to apply
                    nextCut = gs.longestPath + 1;
                }
            }
            
            // once we started cutting, we have to move nodes
            if (index != newIndex) {

                Layer oldLayer = graph.getLayers().get(index);
                Layer newLayer = graph.getLayers().get(newIndex);

                List<LNode> nodesToMove = Lists.newArrayList(oldLayer.getNodes()); 
                for (LNode n : nodesToMove) {
                    
                    // first move the original node 
                    // to allow creating dummy nodes for the reversed edge
                    n.setLayer(newLayer.getNodes().size(), newLayer);                    

                    // at the cut points the edges have to be re-routed
                    if (newIndex == 0) {
                        List<LEdge> incEdges = Lists.newArrayList(n.getIncomingEdges());
                        for (LEdge e : incEdges) {
                            e.reverse(graph, true);
                            graph.setProperty(InternalProperties.CYCLIC, true);

                            // insert proper dummy nodes for the newly created long edge
                            CuttingUtils.insertDummies(graph, e, 1);
                        }
                    }
                }
            }

            // move to the next layer
            newIndex++;
            index++;
        }

       // remove old layers
       ListIterator<Layer> it = graph.getLayers().listIterator();
       while (it.hasNext()) {
           Layer l = it.next();
           if (l.getNodes().isEmpty()) {
               it.remove();
           }
       }
    }
    
}
