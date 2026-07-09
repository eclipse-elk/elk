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
package org.eclipse.elk.alg.common.compaction.oned;

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
                    // don't be confused: when compacting horizontally, the vertical spacing must 
                    // be considered when computing constraints
                    spacing = compactor.spacingsHandler.getVerticalSpacing(cNode1, cNode2);
                } else {
                    spacing = compactor.spacingsHandler.getHorizontalSpacing(cNode1, cNode2);
                }
                
                // add constraint if cNode2 is to the right of cNode1 and could collide if moved
                // horizontally
                if (// '>' avoids simultaneous constraints A->B and B->A
                        (cNode2.hitbox.x > cNode1.hitbox.x 
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
