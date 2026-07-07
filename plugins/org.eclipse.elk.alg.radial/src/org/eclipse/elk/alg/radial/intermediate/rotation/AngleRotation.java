/*******************************************************************************
 * Copyright (c) 2023 Kiel University and others.
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
package org.eclipse.elk.alg.radial.intermediate.rotation;

import org.eclipse.elk.alg.radial.InternalProperties;
import org.eclipse.elk.alg.radial.options.RadialOptions;
import org.eclipse.elk.core.math.KVector;
import org.eclipse.elk.graph.ElkNode;

/**
 * Rotates the entire layout around the origin to a set target angle.
 *
 */
public class AngleRotation implements IRadialRotator {

    @Override
    public void rotate(ElkNode graph) {
        double targetAngle = graph.getProperty(RadialOptions.ROTATION_TARGET_ANGLE);
        
        if (graph.getProperty(RadialOptions.ROTATION_COMPUTE_ADDITIONAL_WEDGE_SPACE)) {
            // Using the target angle as our base alignment we want to further rotate the layout such that a line 
            // following the target angle runs directly through the middle of the wedge between the first and last node.
            ElkNode root = graph.getProperty(InternalProperties.ROOT_NODE);
  
            ElkNode lastNode = (ElkNode) root.getOutgoingEdges().get(root.getOutgoingEdges().size() - 1).getTargets().get(0);
            ElkNode firstNode = (ElkNode) root.getOutgoingEdges().get(0).getTargets().get(0);
            KVector lastVector = new KVector(lastNode.getX() + lastNode.getWidth() / 2, lastNode.getY() + lastNode.getHeight() / 2);
            KVector firstVector = new KVector(firstNode.getX() + firstNode.getWidth() / 2, firstNode.getY() + firstNode.getHeight() / 2);

            // we shift all angles into the range (0,pi] to avoid dealing with negative angles.
            double alpha = targetAngle;
            if (alpha <= 0) {
                alpha += 2*Math.PI;
            }
            
            double wedgeAngle = lastVector.angle(firstVector);
            if (wedgeAngle <= 0) {
                wedgeAngle += 2*Math.PI;
            }
            
            double alignmentAngle = Math.atan2(lastVector.y, lastVector.x);
            if (alignmentAngle <= 0) {
                alignmentAngle += 2*Math.PI;
            }
            
            // alpha (originally targetAngle) is the angle of the incoming edge that we wish to align ourselves with.
            // wedgeAngle is the angle between the first and last nodes of our own layout. For the case of a single 
            // node this is 360 degrees.
            // alignmentAngle is the angle of the vector pointing to the last node i.e. the end part of the segment
            // we rotate the entire layout by subtracting the incoming angle alpha and we add half the wedge angle back
            // to make the alignment go through the center of the wedge. Finally, we need to do a transformation to 
            // make all this work in our downward facing coordinate system. So everything is inverted and we need subtract
            // the result from 180 degrees.
            targetAngle = Math.PI - (alignmentAngle - alpha + wedgeAngle / 2);

        }

        // rotate all nodes around the origin, because the root node is positioned at the origin
        // nodes are positioned with their center on the radius so use that for the rotation
        for (ElkNode node : graph.getChildren()) {
            KVector pos = new KVector(node.getX() + node.getWidth() / 2, node.getY() + node.getHeight() / 2);
            pos.rotate(targetAngle);
            node.setLocation(pos.x - node.getWidth() / 2, pos.y - node.getHeight() / 2);
        }

    }

}
