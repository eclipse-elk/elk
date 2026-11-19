/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
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
package org.eclipse.elk.core.meta

import com.google.inject.Binder
import java.util.Set
import org.eclipse.elk.core.meta.jvmmodel.MelkDocumentationGenerator
import org.eclipse.xtext.generator.IOutputConfigurationProvider
import org.eclipse.xtext.generator.OutputConfiguration
import org.eclipse.xtext.generator.OutputConfigurationProvider
import org.eclipse.elk.core.meta.validation.MelkUniqueClassNameValidator
import org.eclipse.xtext.service.SingletonBinding

/**
 * Use this class to register components to be used at runtime / without the Equinox extension registry.
 */
class MetaDataRuntimeModule extends AbstractMetaDataRuntimeModule {
    
    override bindIGenerator() {
        return MelkDocumentationGenerator
    }
    
    @SingletonBinding(eager = true)
    override bindUniqueClassNameValidator() {
        return MelkUniqueClassNameValidator
    }
    
    // register MelkOutputConfigurationProvider that inserts a configuration used to read files in the project
    override configure(Binder binder) {
        super.configure(binder)
        binder.bind(IOutputConfigurationProvider).to(MelkOutputConfigurationProvider)
    }
    
    static class MelkOutputConfigurationProvider extends OutputConfigurationProvider {
        public static String AD_INPUT = "AD_INPUT"

        override Set<OutputConfiguration> getOutputConfigurations() {
            val configurations = super.outputConfigurations
            val input = new OutputConfiguration(AD_INPUT);
            input.setDescription("Additional Documentation Input Folder");
            input.setOutputDirectory("./");
            configurations.add(input)
            configurations
        }

    }
}
