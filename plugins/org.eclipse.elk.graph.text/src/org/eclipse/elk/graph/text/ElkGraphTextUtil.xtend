/*******************************************************************************
 * Copyright (c) 2017 TypeFox GmbH (http://www.typefox.io) and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.graph.text

import org.eclipse.elk.core.data.LayoutMetaDataService
import org.eclipse.elk.core.options.CoreOptions
import org.eclipse.elk.graph.ElkEdge
import org.eclipse.elk.graph.ElkGraphElement
import org.eclipse.elk.graph.ElkLabel
import org.eclipse.elk.graph.ElkNode
import org.eclipse.elk.graph.ElkPort

/**
 * Utilities for the ELK Graph Text format.
 */
final class ElkGraphTextUtil {
    
    private new() {}
    
    /**
     * Return the algorithm that is responsible for the layout of the given element. Note that this
     * might be ambiguous: e.g. a port of a composite node can be handled both by the algorithm that
     * arranges that node and the algorithm that arranges its container node.
     */
    static def getAlgorithm(ElkGraphElement element) {
        val node = switch element {
            ElkNode:
                element
            ElkPort:
                element.parent?.parent ?: element.parent
            ElkEdge:
                element.containingNode
            ElkLabel: {
                var parent = element.parent
                while (parent instanceof ElkLabel) {
                    parent = parent.parent
                }
                switch parent {
                    ElkNode: parent.parent ?: parent
                    ElkPort: parent.parent?.parent ?: parent.parent
                    ElkEdge: parent.containingNode
                }
            }
        }
        if (node !== null) {
            val algorithmId = node.getProperty(CoreOptions.ALGORITHM)
            if (!algorithmId.nullOrEmpty)
                return LayoutMetaDataService.instance.getAlgorithmDataBySuffix(algorithmId)
        }
    }
    
}