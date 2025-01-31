/*******************************************************************************
 * Copyright (c) 2022 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0  *******************************************************************************/
package org.eclipse.elk.alg.mrtree.p4route;

import java.util.LinkedList;

import org.eclipse.elk.alg.mrtree.graph.TEdge;
import org.eclipse.elk.alg.mrtree.graph.TGraph;
import org.eclipse.elk.alg.mrtree.graph.TNode;
import org.eclipse.elk.alg.mrtree.options.InternalProperties;
import org.eclipse.elk.alg.mrtree.options.MrTreeOptions;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.options.Direction;
import org.eclipse.elk.core.util.Triple;

/**
 * This class keeps track of all the multi level edges bendPoints that come through this node node gap and updates the 
 * registered bendPoints if new edges are added.
 * 
 * @author jnc
 * @author sdo
 */
public class MultiLevelEdgeNodeNodeGap {
    private TNode neighborOne, neighborTwo;
    private LinkedList<Triple<KVector, KVector, TEdge>> bendPoints;
    private TGraph graph;
    private Direction d;
    private double nodeBendpointPadding;
    private boolean onFirstNodeSide, onLastNodeSide;
    
    /**
     * @return the neighborOne
     */
    public TNode getNeighborOne() {
        return neighborOne;
    }

    /**
     * @return the neighborTwo
     */
    public TNode getNeighborTwo() {
        return neighborTwo;
    }

    /**
     * @return the onFirstNodeSide
     */
    public boolean isOnFirstNodeSide() {
        return onFirstNodeSide;
    }

    /**
     * @return the onLastNodeSide
     */
    public boolean isOnLastNodeSide() {
        return onLastNodeSide;
    }

    /**
     * Creates a new MultiLevelEdgeNodeNodeGap object.
     * 
     * @param neighborOne the neighbor with the lower index
     * @param neighborTwo the neighbor with the higher index
     * @param bendTriple the bendPoints between the two nodes
     * @param graph the graph the nodes are in
     */
    public MultiLevelEdgeNodeNodeGap(final TNode neighborOne, final TNode neighborTwo,
            final Triple<KVector, KVector, TEdge> bendTriple, final TGraph graph) {
        this.neighborOne = neighborOne;
        this.neighborTwo = neighborTwo;
        this.graph = graph;
        
        LinkedList<Triple<KVector, KVector, TEdge>> bends = new LinkedList<Triple<KVector, KVector, TEdge>>();
        bends.add(bendTriple);
        this.bendPoints = bends;
        
        d = graph.getProperty(MrTreeOptions.DIRECTION);
        nodeBendpointPadding = graph.getProperty(MrTreeOptions.SPACING_EDGE_NODE);
        
        updateBendPoints();
    }
    
    /**
     * Adds new bendPoints to the gap.
     * 
     * @param newBends the two new bendPoints
     */
    public void addBendPoints(final Triple<KVector, KVector, TEdge> newBends) {
        bendPoints.add(newBends);
        
        if (d.isHorizontal()) {
            bendPoints.sort((x, y) -> Double.compare(x.getThird().getTarget().getPosition().y, 
                    y.getThird().getTarget().getPosition().y));
        } else {
            bendPoints.sort((x, y) -> Double.compare(x.getThird().getTarget().getPosition().x, 
                    y.getThird().getTarget().getPosition().x));
        }
        
        updateBendPoints();
    }
    
    /**
     * Updates the BendPoints according to the new amount of edges that come through this gap.
     */
    private void updateBendPoints() {
        int i = 0, count = bendPoints.size();
        for (Triple<KVector, KVector, TEdge> p : bendPoints) {
            KVector bend1, bend2;
            double bendTmp, interpolation = (i + 1) / (double) (count + 1);
            if (neighborOne == null && neighborTwo == null) {
                return;
            } else if (neighborOne != null && neighborTwo == null) { // Bendpoint should be placed next to last node
                onLastNodeSide = true;
                if (d == Direction.LEFT) {
                    bendTmp = neighborOne.getPosition().y + neighborOne.getSize().y + nodeBendpointPadding * (i + 1);
                    bend1 = new KVector(neighborOne.getProperty(InternalProperties.LEVELMAX) + nodeBendpointPadding, 
                            bendTmp);
                    bend2 = new KVector(neighborOne.getProperty(InternalProperties.LEVELMIN) - nodeBendpointPadding, 
                            bendTmp);
                } else if (d == Direction.RIGHT) {
                    bendTmp = neighborOne.getPosition().y + neighborOne.getSize().y + nodeBendpointPadding * (i + 1);
                    bend1 = new KVector(neighborOne.getProperty(InternalProperties.LEVELMIN) - nodeBendpointPadding, 
                            bendTmp);
                    bend2 = new KVector(neighborOne.getProperty(InternalProperties.LEVELMAX) + nodeBendpointPadding, 
                            bendTmp);
                } else if (d == Direction.UP) {
                    bendTmp = neighborOne.getPosition().x + neighborOne.getSize().x + nodeBendpointPadding * (i + 1);
                    bend1 = new KVector(bendTmp, neighborOne.getProperty(InternalProperties.LEVELMAX) 
                            + nodeBendpointPadding);
                    bend2 = new KVector(bendTmp, neighborOne.getProperty(InternalProperties.LEVELMIN) 
                            - nodeBendpointPadding);
                } else {
                    bendTmp = neighborOne.getPosition().x + neighborOne.getSize().x + nodeBendpointPadding * (i + 1);
                    bend1 = new KVector(bendTmp, neighborOne.getProperty(InternalProperties.LEVELMIN) 
                            - nodeBendpointPadding);
                    bend2 = new KVector(bendTmp, neighborOne.getProperty(InternalProperties.LEVELMAX) 
                            + nodeBendpointPadding);
                }
            } else if (neighborOne != null && neighborTwo != null) { // Bendpoint should be placed in between two nodes
                if (d == Direction.LEFT) {
                    bendTmp = neighborTwo.getPosition().y * interpolation + (neighborOne.getPosition().y 
                            + neighborOne.getSize().y) * (1 - interpolation);
                    bend1 = new KVector(neighborOne.getProperty(InternalProperties.LEVELMAX) 
                            + nodeBendpointPadding, bendTmp);
                    bend2 = new KVector(neighborOne.getProperty(InternalProperties.LEVELMIN) 
                            - nodeBendpointPadding, bendTmp);
                } else if (d == Direction.RIGHT) {
                    bendTmp = neighborTwo.getPosition().y * interpolation + (neighborOne.getPosition().y 
                            + neighborOne.getSize().y) * (1 - interpolation);
                    bend1 = new KVector(neighborOne.getProperty(InternalProperties.LEVELMIN) 
                            - nodeBendpointPadding, bendTmp);
                    bend2 = new KVector(neighborOne.getProperty(InternalProperties.LEVELMAX) 
                            + nodeBendpointPadding, bendTmp);
                } else if (d == Direction.UP) {
                    bendTmp = neighborTwo.getPosition().x * interpolation + (neighborOne.getPosition().x 
                            + neighborOne.getSize().x) * (1 - interpolation);
                    bend1 = new KVector(bendTmp, neighborOne.getProperty(InternalProperties.LEVELMAX) 
                            + nodeBendpointPadding);
                    bend2 = new KVector(bendTmp, neighborOne.getProperty(InternalProperties.LEVELMIN) 
                            - nodeBendpointPadding);
                } else {
                    bendTmp = neighborTwo.getPosition().x * interpolation + (neighborOne.getPosition().x 
                            + neighborOne.getSize().x) * (1 - interpolation);
                    bend1 = new KVector(bendTmp, neighborOne.getProperty(InternalProperties.LEVELMIN) 
                            - nodeBendpointPadding);
                    bend2 = new KVector(bendTmp, neighborOne.getProperty(InternalProperties.LEVELMAX) 
                            + nodeBendpointPadding);
                }
            } else { // Bendpoint should be placed next to the first node
                onFirstNodeSide = true;
                if (d == Direction.LEFT) {
                    bendTmp = neighborTwo.getPosition().y - nodeBendpointPadding * (i + 1);
                    bend1 = new KVector(neighborTwo.getProperty(InternalProperties.LEVELMAX) 
                            + nodeBendpointPadding, bendTmp);
                    bend2 = new KVector(neighborTwo.getProperty(InternalProperties.LEVELMIN) 
                            - nodeBendpointPadding, bendTmp);
                } else if (d == Direction.RIGHT) {
                    bendTmp = neighborTwo.getPosition().y - nodeBendpointPadding * (i + 1);
                    bend1 = new KVector(neighborTwo.getProperty(InternalProperties.LEVELMIN) 
                            - nodeBendpointPadding, bendTmp);
                    bend2 = new KVector(neighborTwo.getProperty(InternalProperties.LEVELMAX) 
                            + nodeBendpointPadding, bendTmp);
                } else if (d == Direction.UP) {
                    bendTmp = neighborTwo.getPosition().x - nodeBendpointPadding * (i + 1);
                    bend1 = new KVector(bendTmp, neighborTwo.getProperty(InternalProperties.LEVELMAX) 
                            + nodeBendpointPadding);
                    bend2 = new KVector(bendTmp, neighborTwo.getProperty(InternalProperties.LEVELMIN) 
                            - nodeBendpointPadding);
                } else {
                    bendTmp = neighborTwo.getPosition().x - nodeBendpointPadding * (i + 1);
                    bend1 = new KVector(bendTmp, neighborTwo.getProperty(InternalProperties.LEVELMIN) 
                            - nodeBendpointPadding);
                    bend2 = new KVector(bendTmp, neighborTwo.getProperty(InternalProperties.LEVELMAX) 
                            + nodeBendpointPadding);
                }
            }
            
            // Commit new values
            p.getFirst().x = bend1.x;
            p.getFirst().y = bend1.y;
            p.getSecond().x = bend2.x;
            p.getSecond().y = bend2.y;
            
            i++;
        }
    }
}
