/*******************************************************************************
 * Copyright (c) 2023 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 
 *******************************************************************************/
package org.eclipse.elk.alg.vertiflex;

import org.eclipse.elk.graph.ElkEdge;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.elk.graph.util.ElkGraphUtil;

/** A class for smaller, independent calculation units. */
public final class VertiFlexUtil {
    
    /**
     * Private constructor to prevent initialization.
     */
    private VertiFlexUtil() { }
    
    /**
    * Computes the root node of a graph.
    * @param graph
    * @return Root node of graph.
    */
   public static ElkNode findRoot(final ElkNode graph) {
       for (ElkNode child : graph.getChildren()) {
           Iterable<ElkEdge> incomingEdges = ElkGraphUtil.allIncomingEdges(child);
           if (!incomingEdges.iterator().hasNext()) {
               return child;
           }
       }
       return null;
   }
}
