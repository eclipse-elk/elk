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
import org.eclipse.elk.graph.properties.IProperty
import org.eclipse.xtext.common.services.DefaultTerminalConverters
import org.eclipse.xtext.conversion.IValueConverter
import org.eclipse.xtext.conversion.ValueConverter

class ElkGraphValueConverterService extends DefaultTerminalConverters {

    @Inject
    NumberValueConverter numberValueConverter

    @ValueConverter(rule="Number")
    def IValueConverter<Double> Number() {
        numberValueConverter
    }

    @Inject
    PropertyKeyValueConverter propertyKeyValueConverter

    @ValueConverter(rule="PropertyKey")
    def IValueConverter<IProperty<?>> PropertyKey() {
        propertyKeyValueConverter
    }
    
    @Inject
    PropertyValueValueConverter propertyValueValueConverter

    @ValueConverter(rule="StringValue")
    def IValueConverter<Object> StringValue() {
        propertyValueValueConverter
    }

    @ValueConverter(rule="QualifiedIdValue")
    def IValueConverter<Object> QualifiedIdValue() {
        propertyValueValueConverter
    }

    @ValueConverter(rule="NumberValue")
    def IValueConverter<Object> NumberValue() {
        propertyValueValueConverter
    }

    @ValueConverter(rule="BooleanValue")
    def IValueConverter<Object> BooleanValue() {
        propertyValueValueConverter
    }

}
