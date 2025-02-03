/*******************************************************************************
 * Copyright (c) 2023 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 
 *******************************************************************************/
package org.eclipse.elk.alg.vertiflex.p2relative;

import org.eclipse.elk.alg.vertiflex.InternalProperties;
import org.eclipse.elk.graph.ElkNode;

/**
 * Comparator for sorting ElkNodes according to their y positions and model order.
 *
 */
class NodeComparator implements java.util.Comparator<ElkNode> {

    private boolean invert = false;

    /** Default constructor */
    public NodeComparator() {}
    
    /** Constructor for inverted comparator. */
    public NodeComparator(boolean invert) {
        this.invert = invert;
    }
    
    @Override
    public int compare(final ElkNode a, final ElkNode b) {

        int sortYresult;
        if (!invert) {
            sortYresult =  Double.compare(a.getY(), b.getY());
        } else {
            sortYresult = Double.compare(b.getY(), a.getY());
        }
        if (sortYresult == 0) {
            int intSortresult =  Integer.compare(a.getProperty(InternalProperties.NODE_MODEL_ORDER), 
                    b.getProperty(InternalProperties.NODE_MODEL_ORDER));
            return intSortresult;
        } else {
            return sortYresult;
        }
    }
}
