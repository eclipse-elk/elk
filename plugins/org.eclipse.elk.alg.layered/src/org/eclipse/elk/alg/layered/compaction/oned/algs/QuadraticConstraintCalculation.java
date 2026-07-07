/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
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
package org.eclipse.elk.alg.layered.compaction.oned.algs;

import org.eclipse.elk.alg.layered.compaction.oned.CNode;
import org.eclipse.elk.alg.layered.compaction.oned.CompareFuzzy;
import org.eclipse.elk.alg.layered.compaction.oned.OneDimensionalCompactor;

/**
 * Creates a constraint between CNodes A and B if B collides with the right shadow of A considering
 * vertical spacing.
 */
public class QuadraticConstraintCalculation implements IConstraintCalculationAlgorithm {

    @Override
    public void calculateConstraints(final OneDimensionalCompactor compactor) {

        // resetting constraints
        for (CNode cNode : compactor.cGraph.cNodes) {
            cNode.constraints.clear();
        }
        
        // inferring constraints from hitbox intersections
        for (CNode cNode1 : compactor.cGraph.cNodes) {
            for (CNode cNode2 : compactor.cGraph.cNodes) {
                // no self constraints
                if (cNode1 == cNode2) {
                    continue;
                }
                // no constraints between nodes of the same group
                if (cNode1.cGroup != null && cNode1.cGroup == cNode2.cGroup) {
                    continue;
                }
                
                double spacing;
                if (compactor.direction.isHorizontal()) {
                    //spacing = Math.min(cNode1.getVerticalSpacing(), cNode2.getVerticalSpacing());
                    spacing = compactor.spacingsHandler.getVerticalSpacing(cNode1, cNode2);
                } else {
                    //spacing = Math.min(cNode1.getHorizontalSpacing(), cNode2.getHorizontalSpacing());
                    spacing = compactor.spacingsHandler.getHorizontalSpacing(cNode1, cNode2);
                }
                
                // add constraint if cNode2 is to the right of cNode1 and could collide if moved
                // horizontally
                // exclude parentNodes because they don't constrain their north/south segments
                if (cNode1 != cNode2.parentNode
                        // '>' avoids simultaneous constraints A->B and B->A
                        && (cNode2.hitbox.x > cNode1.hitbox.x 
                                // 
                                || (cNode1.hitbox.x == cNode2.hitbox.x 
                                && cNode1.hitbox.width < cNode2.hitbox.width))
                                
                                && CompareFuzzy.gt(cNode2.hitbox.y + cNode2.hitbox.height + spacing,
                                        cNode1.hitbox.y)
                                        
                                        && CompareFuzzy.lt(cNode2.hitbox.y, 
                                                cNode1.hitbox.y + cNode1.hitbox.height + spacing)) {
                    
                    cNode1.constraints.add(cNode2);
                }
            }
        }
    }

}
