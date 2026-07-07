/*******************************************************************************
 * Copyright (c) 2020 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.graph.json.text.conversion

import org.eclipse.elk.core.data.LayoutMetaDataService
import org.eclipse.elk.core.data.LayoutOptionData
import org.eclipse.elk.core.util.internal.LayoutOptionProxy
import org.eclipse.elk.graph.impl.ElkPropertyToValueMapEntryImpl
import org.eclipse.elk.graph.properties.IProperty
import org.eclipse.elk.graph.properties.IPropertyValueProxy
import org.eclipse.xtext.conversion.ValueConverterException
import org.eclipse.xtext.conversion.impl.AbstractValueConverter
import org.eclipse.xtext.nodemodel.INode

/**
 * Note that this converter is derived from ElkText's. The main difference 
 * is that all values but numbers and bools are in double-quotes.
 */
class PropertyValueValueConverter extends AbstractValueConverter<Object> {

    override toString(Object value) throws ValueConverterException {
        // Here we assume that the toString() implementation of the proxy yields a meaningful value
        val v = if(value instanceof IPropertyValueProxy) value.toString else value

        if (v === null)
            return 'null'
        else if (v instanceof Double && Math.floor(v as Double) == v as Double)
            return Integer.toString((v as Double).intValue)
        else if (v instanceof Boolean || v instanceof Number)
            return v.toString
        else if (v instanceof String) {
            if (v == 'true' || v == 'false' || v == 'null')
                return v
            try {
                Double.parseDouble(v)
                return v
            } catch (NumberFormatException e) {
                return '"' + v + '"'
            }
        } else
            return '"' + v.toString + '"'
    }

    override toValue(String string, INode node) throws ValueConverterException {
        if (!string.nullOrEmpty && string != 'null') {
            val unquotedString = string.unquoteIfNecessary
            val semanticElem = node?.semanticElement
            if (semanticElem instanceof ElkPropertyToValueMapEntryImpl) {
                val option = semanticElem.key.toLayoutOption
                val value = option?.parseValue(unquotedString)
                if (value !== null)
                    return value
            }
            return new LayoutOptionProxy(unquotedString)
        }
    }

    private def unquoteIfNecessary(String s) {
        if (s.length >= 2 && (s.startsWith('"') && s.endsWith('"') || s.startsWith("'") && s.endsWith("'")))
            return s.substring(1, s.length - 1)
        else
            return s
    }

    private def toLayoutOption(IProperty<?> property) {
        if (property instanceof LayoutOptionData)
            return property
        else if (property !== null)
            return LayoutMetaDataService.instance.getOptionData(property.id)
    }

}
