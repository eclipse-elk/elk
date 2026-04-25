/*******************************************************************************
 * Copyright (c) 2010, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0
 *******************************************************************************/
package org.eclipse.elk.alg.layered.graph.transform;

import org.eclipse.elk.alg.layered.graph.LGraph;
import org.eclipse.elk.alg.layered.graph.LGraphElement;
import org.eclipse.elk.alg.layered.options.InternalProperties;
import org.eclipse.elk.graph.ElkGraphElement;
import org.eclipse.elk.graph.ElkNode;
import org.eclipse.emf.ecore.EObject;


/**
 * Manages the transformation of ELK Graphs to LayeredGraphs. Sets the
 * {@link org.eclipse.elk.alg.layered.options.LayeredOptions#GRAPH_PROPERTIES GRAPH_PROPERTIES}
 * property on imported graphs.
 * 
 * @author msp
 * @author cds
 * @see ElkGraphImporter
 * @see ElkGraphLayoutTransferrer
 */
public class ElkGraphTransformer implements IGraphTransformer<ElkNode> {
    
    @Override
    public LGraph importGraph(final ElkNode graph) {
        return new ElkGraphImporter().importGraph(graph);
    }
    
    @Override
    public void applyLayout(final LGraph layeredGraph) {
        new ElkGraphLayoutTransferrer().applyLayout(layeredGraph);
    }
    

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Utility
    
    /**
     * Returns an identifier string for the original ElkGraph element.
     * 
     * @param element an LGraph element
     * @return the original identifier, or {@code null} if none is defined
     */
    public static String getOriginIdentifier(final LGraphElement element) {
        Object origin = element.getProperty(InternalProperties.ORIGIN);
        if (origin instanceof ElkGraphElement) {
            return getIdentifier((ElkGraphElement) origin);
        }
        return null;
    }
    
    private static String getIdentifier(final ElkGraphElement element) {
        String id = element.getIdentifier();
        if (!(id == null || id.isEmpty())) {
            EObject container = element.eContainer();
            if (container instanceof ElkGraphElement) {
                String parentId = getIdentifier((ElkGraphElement) container);
                if (parentId != null) {
                    return parentId + '.' + id;
                }
            }
            return id;
        }
        return null;
    }
    
}
