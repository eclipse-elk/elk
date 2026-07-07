/*******************************************************************************
 * Copyright (c) 2024 Kiel University and others and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * This Source Code may also be made available under the following Secondary
 * Licenses when the conditions for such availability set forth in the Eclipse
 * Public License v. 2.0 are satisfied: GPL-3.0 which is available at
 * https://www.gnu.org/licenses/gpl-3.0-standalone.html.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later 
 *******************************************************************************/
package org.eclipse.elk.core.meta.validation

import com.google.inject.Inject
import org.eclipse.xtext.common.types.JvmDeclaredType
import org.eclipse.xtext.common.types.TypesPackage
import org.eclipse.xtext.naming.QualifiedName
import org.eclipse.xtext.xbase.validation.UniqueClassNameValidator
import org.eclipse.xtext.resource.impl.ResourceDescriptionsProvider

/**
 * Overwrites the UniqueClassNameValidator to only check for melk files in src folders. This probably
 * prevents some special use cases where ELK is used as a library and a new melk file is created.
 *
 */
class MelkUniqueClassNameValidator extends UniqueClassNameValidator {
    
    @Inject
    var ResourceDescriptionsProvider resourceDescriptionsProvider;
    
    override doCheckUniqueName(QualifiedName name, JvmDeclaredType type) {
        val index = resourceDescriptionsProvider.getResourceDescriptions(type.eResource())
        val others = index.getExportedObjects(TypesPackage.Literals.JVM_DECLARED_TYPE, name, false)
        return checkUniqueInIndex(type, others.filter[it.EObjectURI.segments.exists[it.contains("src")]])
    }
}
