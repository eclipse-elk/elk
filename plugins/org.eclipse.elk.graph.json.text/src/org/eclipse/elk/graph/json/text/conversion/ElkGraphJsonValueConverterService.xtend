/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
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
package org.eclipse.elk.graph.json.text.conversion

import com.google.inject.Inject
import org.eclipse.elk.graph.properties.IProperty
import org.eclipse.xtext.common.services.DefaultTerminalConverters
import org.eclipse.xtext.conversion.IValueConverter
import org.eclipse.xtext.conversion.ValueConverter

class ElkGraphJsonValueConverterService extends DefaultTerminalConverters {

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

    @ValueConverter(rule="NumberValue")
    def IValueConverter<Object> NumberValue() {
        propertyValueValueConverter
    }

    @ValueConverter(rule="BooleanValue")
    def IValueConverter<Object> BooleanValue() {
        propertyValueValueConverter
    }

}
