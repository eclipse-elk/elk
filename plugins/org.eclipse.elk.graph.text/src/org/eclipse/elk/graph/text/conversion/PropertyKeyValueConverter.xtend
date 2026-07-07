/*******************************************************************************
 * Copyright (c) 2016 TypeFox GmbH (http://www.typefox.io) and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.graph.text.conversion

import com.google.inject.Inject
import com.google.inject.Provider
import org.eclipse.elk.core.data.LayoutMetaDataService
import org.eclipse.elk.graph.properties.IProperty
import org.eclipse.elk.graph.text.services.ElkGraphGrammarAccess
import org.eclipse.xtext.conversion.ValueConverterException
import org.eclipse.xtext.conversion.impl.AbstractValueConverter
import org.eclipse.xtext.conversion.impl.IDValueConverter
import org.eclipse.xtext.nodemodel.INode
import org.eclipse.xtext.util.Strings

class PropertyKeyValueConverter extends AbstractValueConverter<IProperty<?>> {
    
    IDValueConverter idValueConverter
    
    @Inject
    def void initialize(Provider<IDValueConverter> idValueConverterProvider, ElkGraphGrammarAccess grammarAccess) {
        this.idValueConverter = idValueConverterProvider.get => [
            rule = grammarAccess.IDRule
        ]
    }
    
    override toString(IProperty<?> value) throws ValueConverterException {
        if (value === null)
            throw new ValueConverterException("IProperty value may not be null.", null, null)
        val metaDataService = LayoutMetaDataService.instance
        val split = Strings.split(value.id, '.')
        var String suffix
        var i = split.size - 1
        while (i >= 0) {
            if (suffix === null)
                suffix = idValueConverter.toString(split.get(i--))
            else
                suffix = idValueConverter.toString(split.get(i--)) + '.' + suffix
                
            val knownOption = metaDataService.getOptionDataBySuffix(suffix)
            // make sure to always add the option's group
            if (knownOption !== null && (knownOption.group.nullOrEmpty || suffix.contains(knownOption.group)))
                return suffix 
        }
        return split.map[idValueConverter.toString(it)].join('.')
    }
    
    override toValue(String string, INode node) throws ValueConverterException {
        if (string.nullOrEmpty)
            throw new ValueConverterException("Cannot convert empty string to a property identifier.", node, null)
        val idSuffix = Strings.split(string, '.').map[idValueConverter.toValue(it, node)].join('.')
        val optionData = LayoutMetaDataService.instance.getOptionDataBySuffix(idSuffix)
        if (optionData === null)
            // If the currently given id is not a known layout option, a new property instance is created and 
            // returned that shall not be registered with the LayoutMetaDataService.
            // The thinking behind this is to allow custom "key: value" pairs within an elkt document, 
            // where the purpose of the custom elements is unknown to the ELK (but hopefully known to the user). 
            // Note that the specified values will yield 'LayoutOptionProxy' instances.
            return new org.eclipse.elk.graph.properties.Property<Object>(string)
          
        return optionData
    }
    
}