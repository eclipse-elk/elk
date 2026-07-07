/*******************************************************************************
 * Copyright (c) 2017 TypeFox GmbH (http://www.typefox.io) and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.graph.json.text.serializer

import org.eclipse.emf.ecore.EObject
import org.eclipse.emf.ecore.EStructuralFeature
import org.eclipse.xtext.serializer.sequencer.LegacyTransientValueService

import static org.eclipse.elk.graph.ElkGraphPackage.Literals.*

class ElkGraphJsonTransientValueService extends LegacyTransientValueService {
    
    override isListTransient(EObject semanticObject, EStructuralFeature feature) {
        switch feature {
            case ELK_CONNECTABLE_SHAPE__INCOMING_EDGES,
            case ELK_CONNECTABLE_SHAPE__OUTGOING_EDGES,
            case ELK_EDGE_SECTION__INCOMING_SECTIONS:
                return ListTransient.YES
            default:
                return super.isListTransient(semanticObject, feature)
        }
    }
    
}