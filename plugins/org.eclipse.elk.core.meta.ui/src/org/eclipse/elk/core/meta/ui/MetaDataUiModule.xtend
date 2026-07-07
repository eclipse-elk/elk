/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.meta.ui

import org.eclipse.elk.core.meta.validation.MelkUniqueClassNameValidator
import org.eclipse.xtend.lib.annotations.FinalFieldsConstructor
import org.eclipse.xtext.service.SingletonBinding

/**
 * Use this class to register components to be used within the Eclipse IDE.
 */
@FinalFieldsConstructor
class MetaDataUiModule extends AbstractMetaDataUiModule {
    
// If the Metacompiler has hiccups this can be used to always generate the Java source files, 
// regardles of any errors.
//    override bindIShouldGenerate() {
//        return Always
//    }
    
    @SingletonBinding(eager = true)
    override bindUniqueClassNameValidator() {
        return MelkUniqueClassNameValidator
    }
}
