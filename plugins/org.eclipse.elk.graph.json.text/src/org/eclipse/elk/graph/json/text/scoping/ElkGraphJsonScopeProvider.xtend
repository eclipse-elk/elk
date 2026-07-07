/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
 *
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.graph.json.text.scoping

import org.eclipse.elk.graph.ElkEdgeSection
import org.eclipse.emf.ecore.EObject
import org.eclipse.emf.ecore.EReference
import org.eclipse.xtext.naming.QualifiedName
import org.eclipse.xtext.scoping.IScope

import static org.eclipse.elk.graph.ElkGraphPackage.Literals.*
import static org.eclipse.xtext.scoping.Scopes.*

/**
 * This class contains custom scoping description.
 * 
 * See https://www.eclipse.org/Xtext/documentation/303_runtime_concepts.html#scoping
 * on how and when to use it.
 */
class ElkGraphJsonScopeProvider extends AbstractElkGraphJsonScopeProvider {
    
    override getScope(EObject context, EReference reference) {
        switch reference {
            case ELK_EDGE_SECTION__OUTGOING_SECTIONS:
                getScopeForEdgeSection(context as ElkEdgeSection)
            default:
                super.getScope(context, reference)
        }
    }
    
    private def getScopeForEdgeSection(ElkEdgeSection section) {
        val edge = section.parent
        return scopeFor(edge.sections, [QualifiedName.create(identifier)], IScope.NULLSCOPE)
    }

}
