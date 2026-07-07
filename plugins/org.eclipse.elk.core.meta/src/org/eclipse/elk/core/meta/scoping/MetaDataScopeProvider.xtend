/*******************************************************************************
 * Copyright (c) 2016 TypeFox GmbH and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.meta.scoping

import com.google.inject.Inject
import org.eclipse.emf.ecore.EObject
import org.eclipse.emf.ecore.EReference
import org.eclipse.xtext.scoping.impl.ImportedNamespaceAwareLocalScopeProvider

import static org.eclipse.elk.core.meta.metaData.MetaDataPackage.Literals.*

/**
 * Cross-references that do not target JVM elements must be scoped with the default Xtext means
 * instead of the Xbase scoping.
 */
class MetaDataScopeProvider extends AbstractMetaDataScopeProvider {
    
    @Inject ImportedNamespaceAwareLocalScopeProvider delegate
    
    override getScope(EObject context, EReference reference) {
        switch reference {
            case MD_OPTION_DEPENDENCY__TARGET,
            case MD_OPTION_SUPPORT__OPTION:
                delegate.getScope(context, reference)
            default:
                super.getScope(context, reference)
        }
    }
    
}