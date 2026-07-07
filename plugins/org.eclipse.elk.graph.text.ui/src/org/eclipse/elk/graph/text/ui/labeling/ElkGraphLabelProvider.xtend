/*******************************************************************************
 * Copyright (c) 2016 TypeFox GmbH (http://www.typefox.io) and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.graph.text.ui.labeling

import com.google.inject.Inject
import org.eclipse.elk.graph.ElkEdge
import org.eclipse.elk.graph.ElkLabel
import org.eclipse.elk.graph.ElkNode
import org.eclipse.elk.graph.ElkPort
import org.eclipse.emf.edit.ui.provider.AdapterFactoryLabelProvider
import org.eclipse.xtext.naming.IQualifiedNameProvider
import org.eclipse.xtext.ui.label.DefaultEObjectLabelProvider

/**
 * Provides labels for ElkGraph classes.
 */
class ElkGraphLabelProvider extends DefaultEObjectLabelProvider {
    
    @Inject extension IQualifiedNameProvider

    @Inject
    new(AdapterFactoryLabelProvider delegate) {
        super(delegate)
    }

    def image(ElkNode node) {
        if (!node.eIsProxy && node.parent === null) {
            'elkgraph.gif'
        } else {
            'elknode.gif'
        }
    }
    
    def image(ElkPort port) {
        'elkport.gif'
    }
    
    def text(ElkEdge edge) {
        val result = new StringBuilder
        if (!edge.identifier.nullOrEmpty)
            result.append(edge.identifier).append(': ')
        edge.sources.forEach[ s, i |
            if (i > 0)
                result.append(', ')
            result.append(s.fullyQualifiedName?.toString ?: '?')
        ]
        result.append(' \u2192 ')
        edge.targets.forEach[ t, i |
            if (i > 0)
                result.append(', ')
            result.append(t.fullyQualifiedName?.toString ?: '?')
        ]
        return result.toString
    }
    
    def image(ElkEdge edge) {
        'elkedge.gif'
    }
    
    def text(ElkLabel label) {
        val result = new StringBuilder
        if (!label.identifier.nullOrEmpty)
            result.append(label.identifier).append(': ')
        result.append('"').append(label.text ?: '').append('"')
        return result.toString
    }
    
    def image(ElkLabel label) {
        'elklabel.gif'
    }
    
}
