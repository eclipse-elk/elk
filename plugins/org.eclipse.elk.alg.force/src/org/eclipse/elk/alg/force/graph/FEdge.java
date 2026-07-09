/*******************************************************************************
 * Copyright (c) 2011, 2015 Kiel University and others.
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
package org.eclipse.elk.alg.force.graph;

import java.util.ArrayList;
import java.util.List;

import org.eclipse.elk.core.math.ElkMath;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.core.math.KVectorChain;
import org.eclipse.elk.graph.properties.MapPropertyHolder;

/**
 * A physico-virtual representation of an edge, including a list of associated bend points.
 * 
 * @author owo
 */
public final class FEdge extends MapPropertyHolder {

    /** the serial version UID. */
    private static final long serialVersionUID = 4387555754824186467L;
    
    /** the bend points of the edge. */
    private List<FBendpoint> bendpoints = new ArrayList<>();
    /** the labels of the edge. */
    private List<FLabel> labels = new ArrayList<>();
    /** the source node of the edge. */
    private FNode source;
    /** the target node of the edge. */
    private FNode target;

    @Override
    public String toString() {
        if (source != null && target != null) {
            return source.toString() + "->" + target.toString();
        } else {
            return "e_" + hashCode();
        }
    }
    
    /**
     * Returns the source node.
     * 
     * @return the source node
     */
    public FNode getSource() {
        return source;
    }

    /**
     * Returns the target node.
     * 
     * @return the target node
     */
    public FNode getTarget() {
        return target;
    }

    /**
     * Returns the list of bend points associated with this edge.
     * 
     * @return list of bend points
     */
    public List<FBendpoint> getBendpoints() {
        return bendpoints;
    }
    
    /**
     * Returns the list of labels associated with this edge.
     * 
     * @return list of labels
     */
    public List<FLabel> getLabels() {
        return labels;
    }

    /**
     * Returns the docking point at the source node.
     * 
     * @return the source docking point
     */
    public KVector getSourcePoint() {
        KVector v = target.getPosition().clone().sub(source.getPosition());
        ElkMath.clipVector(v, source.getSize().x, source.getSize().y);
        return v.add(source.getPosition());
    }

    /**
     * Returns the docking point at the target node.
     * 
     * @return the target docking point
     */
    public KVector getTargetPoint() {
        KVector v = source.getPosition().clone().sub(target.getPosition());
        ElkMath.clipVector(v, target.getSize().x, target.getSize().y);
        return v.add(target.getPosition());
    }
    
    /**
     * Sets the source vertex.
     * 
     * @param theSource
     *            the source vertex set to
     */
    public void setSource(final FNode theSource) {
        source = theSource;
    }
    
    /**
     * Sets the target vertex.
     * 
     * @param theTarget
     *            the target vertex
     */
    public void setTarget(final FNode theTarget) {
        target = theTarget;
    }
    
    /**
     * Returns a vector chain with all bend points and source and target point.
     * 
     * @return a vector chain for the edge
     */
    public KVectorChain toVectorChain() {
        KVectorChain vectorChain = new KVectorChain();
        vectorChain.add(getSourcePoint());
        for (FBendpoint bendPoint : bendpoints) {
            vectorChain.add(bendPoint.getPosition());
        }
        vectorChain.add(getTargetPoint());
        return vectorChain;
    }
    
    /**
     * Distribute the bend points evenly on the edge.
     */
    public void distributeBendpoints() {
        int count = bendpoints.size();
        if (count > 0) {
            KVector sourcePos = source.getPosition();
            KVector targetPos = target.getPosition();
            KVector incr = targetPos.clone().sub(sourcePos).scale(1 / (double) (count + 1));
            KVector pos = sourcePos.clone();
            for (FBendpoint bendPoint : bendpoints) {
                bendPoint.getPosition().x = pos.x + incr.x;
                bendPoint.getPosition().y = pos.y + incr.y;
                pos.add(incr);
            }
        }
    }

}
