/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
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
package org.eclipse.elk.alg.layered.intermediate.compaction;

import java.util.Collections;

import org.eclipse.elk.alg.common.compaction.oned.CNode;
import org.eclipse.elk.alg.common.compaction.oned.ICompactionAlgorithm;
import org.eclipse.elk.alg.common.compaction.oned.ISpacingsHandler;
import org.eclipse.elk.alg.common.compaction.oned.OneDimensionalCompactor;
import org.eclipse.elk.alg.layered.graph.LEdge;
import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LNode;
import org.eclipse.elk.alg.layered.graph.LNode.NodeType;
import org.eclipse.elk.alg.layered.options.GraphCompactionStrategy;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.alg.layered.options.LayeredOptions;
import org.eclipse.elk.alg.layered.options.Spacings;
import org.eclipse.elk.core.alg.ILayoutProcessor;
import org.eclipse.elk.core.options.Direction;
import org.eclipse.elk.core.util.IElkProgressMonitor;

/**
 * This processor applies additional compaction to an already routed graph and can be executed after
 * {@link org.eclipse.elk.alg.layered.p5edges.OrthogonalEdgeRouter OrthogonalEdgeRouter}.
 * Therefore nodes and vertical segments of edges are repositioned in the specified direction where
 * the position is minimal considering the desired spacing between elements.
 * 
 * <p>
 *  Since the used <em>locking</em> functionality relies on the direction of
 *  incoming and outgoing edges, this processor is required to be executed before the
 *  {@link org.eclipse.elk.alg.layered.intermediate.ReversedEdgeRestorer ReversedEdgeRestorer}.
 * </p>
 * 
 * <h3>Splines</h3>
 * The processor is able to cope with edges that are routed using spline segments. 
 * North/south ports have not been tested in this scenario, however, and are likely to yield undesired results.  
 * 
 * <dl>
 *  <dt>Precondition:</dt>
 *   <dd>The edges have been routed</dd>
 *  <dt>Postcondition:</dt>
 *   <dd>Nodes and edges are positioned compact without colliding.</dd>
 *  <dt>Slots:</dt>
 *   <dd>After phase 5.</dd>
 *  <dt>Same-slot dependencies:</dt>
 *   <dd>Before {@link org.eclipse.elk.alg.layered.intermediate.LabelDummyRemover LabelDummyRemover}.
 *   Otherwise labels would be ignored during compaction.
 *   </dd>
 *   <dd>Before {@link org.eclipse.elk.alg.layered.intermediate.ReversedEdgeRestorer
 *       ReversedEdgeRestorer}</dd>
 * </dl>
 */
public class HorizontalGraphCompactor implements ILayoutProcessor<LGraph> {

    /**
     * Compaction algorithm based on the network simplex algorithm presented by Gansner et al. in
     * the context of layering.
     * 
     * @see org.eclipse.elk.alg.common.networksimplex.NetworkSimplex
     */
    public static final ICompactionAlgorithm NETWORK_SIMPLEX_COMPACTION =
            new NetworkSimplexCompaction();
    
    private LGraph lGraph;
    
    @Override
    public void process(final LGraph layeredGraph, final IElkProgressMonitor progressMonitor) {

        GraphCompactionStrategy strategy = layeredGraph.getProperty(LayeredOptions.COMPACTION_POST_COMPACTION_STRATEGY);
        if (strategy == GraphCompactionStrategy.NONE) {
            return;
        }
        progressMonitor.begin("Horizontal Compaction", 1);
        
        this.lGraph = layeredGraph;
        
        // the layered graph is transformed into a CGraph that is passed to OneDimensionalCompactor
        LGraphToCGraphTransformer transformer = new LGraphToCGraphTransformer();
        OneDimensionalCompactor odc =
                new OneDimensionalCompactor(transformer.transform(layeredGraph));
        
        // consider special spacing requirements of the lgraph's nodes and edges
        odc.setSpacingsHandler(specialSpacingsHandler);
        
        // ---
        // select constraint algorithm
        // - 
        switch (layeredGraph.getProperty(LayeredOptions.COMPACTION_POST_COMPACTION_CONSTRAINTS)) {
            case SCANLINE:
                odc.setConstraintAlgorithm(new EdgeAwareScanlineConstraintCalculation(lGraph));
                break;
            default:
                odc.setConstraintAlgorithm(OneDimensionalCompactor.QUADRATIC_CONSTRAINTS);
        }
        
        // ---
        // select compaction strategy
        // - 
        switch (strategy) {
        case LEFT:
            odc.compact();
            break;
            
        case RIGHT:
            odc.changeDirection(Direction.RIGHT)
               .compact();
            break;
            
        case LEFT_RIGHT_CONSTRAINT_LOCKING:
            // lock CNodes if they are not constrained
            odc.compact()
               .changeDirection(Direction.RIGHT)
               .setLockFunction((node, dir) -> node.cGroup.outDegreeReal == 0)
               .compact();
            break;
            
        case LEFT_RIGHT_CONNECTION_LOCKING:
            // compacting left, locking all CNodes that have fewer connections to the right,
            // then compacting right to shorten unnecessary long edges
            odc.compact()
               .changeDirection(Direction.RIGHT)
               .setLockFunction((node, dir) -> transformer.getLockMap().get(node).get(dir))
               .compact();
            break;
         
        case EDGE_LENGTH:
            odc.setCompactionAlgorithm(NETWORK_SIMPLEX_COMPACTION)
               .compact();
            break;

        default:
            // nobody should get here
            break;
        }
        
        // since changeDirection may transform hitboxes, the final direction has to be LEFT again
        odc.finish();

        // applying the compacted positions to the LGraph and updating its size and offset
        transformer.applyLayout();
        
        progressMonitor.done();
    }
    
    /**
     * @return true if both passed nodes originate from the same (set) of
     *         {@link org.eclipse.elk.alg.layered.graph.LEdge LEdge}.
     */
    public static boolean isVerticalSegmentsOfSameEdge(final CNode cNode1, final CNode cNode2) {
        VerticalSegment v1 = getVerticalSegmentOrNull(cNode1);
        VerticalSegment v2 = getVerticalSegmentOrNull(cNode2);
        return
            // if we only want north/south segments we could use the following
            (v1 != null && v2 != null)
            // this might seem quite expensive but in most cases the sets
            // contain only one element
            && !Collections.disjoint(v1.representedLEdges, v2.representedLEdges);
    }
    
    /**
     * @return the represented {@link LNode} (if set as {@link CNode#origin}, or {@code null}.
     */
    public static LNode getLNodeOrNull(final CNode cNode) {
        if (cNode.origin instanceof LNode) {
            return (LNode) cNode.origin;
        }
        return null;
    }

    /**
     * @return the represented {@link VerticalSegment} (if set as {@link CNode#origin}, or {@code null}.
     */
    public static VerticalSegment getVerticalSegmentOrNull(final CNode cNode) {
        if (cNode.origin instanceof VerticalSegment) {
            return (VerticalSegment) cNode.origin;
        }
        return null;
    }
    
    /**
     * An implementation of a {@link ISpacingsHandler} that is able to cope with the special
     * requirements of {@link LGraph}s. For instance, there are special cases for the spacing
     * between {@link LEdge}s as opposed to {@link LNode}s.
     */
    private final ISpacingsHandler specialSpacingsHandler = new ISpacingsHandler() {
        
        @Override
        public double getHorizontalSpacing(final CNode cNode1, final CNode cNode2) {
                        
            // joining north/south segments that belong to the same edge 
            // by setting their spacing to 0
            if (isVerticalSegmentsOfSameEdge(cNode1, cNode2)) {
                return 0;
            }
            
            // get the underlying LNodes 
            LNode node1 = getLNodeOrNull(cNode1);
            LNode node2 = getLNodeOrNull(cNode2);

            // if either of the two involved nodes represents an external port,
            //  it's ok to move the port as close as possible since it will be 
            //  positioned correctly by another processor later on
            if ((node1 != null && node1.getType() == NodeType.EXTERNAL_PORT)
                    || (node2 != null && node2.getType() == NodeType.EXTERNAL_PORT)) {
                return 0;
            }
            
            // default behavior, query the Spacings object
            Spacings spacings = lGraph.getProperty(InternalProperties.SPACINGS);
            
            // either get the spacing for the node's type or if there is no 
            //  corresponding element in the original graph, we use the 
            //  long edge dummy spacing.
            return spacings.getHorizontalSpacing(
                    node1 != null ? node1.getType() : NodeType.LONG_EDGE, 
                    node2 != null ? node2.getType() : NodeType.LONG_EDGE);
        }
        
        @Override
        public double getVerticalSpacing(final CNode cNode1, final CNode cNode2) {

            // Note that this method is only called when the quadratic constraint calculation is used
            //  (not if the edge aware scanline procedure is used).
            
            // joining north/south segments that belong to the same edge 
            // by setting their vertical spacing to 1, 
            // i.e. they overlap in the y dimension which results in a constraint
            if (isVerticalSegmentsOfSameEdge(cNode1, cNode2)) {
                return 1;
            }

            // get the underlying LNodes
            LNode node1 = getLNodeOrNull(cNode1);
            LNode node2 = getLNodeOrNull(cNode2);

            // default behavior, query the Spacings object
            Spacings spacings = lGraph.getProperty(InternalProperties.SPACINGS);

            return spacings.getVerticalSpacing(
                    node1 != null ? node1.getType() : NodeType.LONG_EDGE, 
                    node2 != null ? node2.getType() : NodeType.LONG_EDGE);
        }
    };

}
