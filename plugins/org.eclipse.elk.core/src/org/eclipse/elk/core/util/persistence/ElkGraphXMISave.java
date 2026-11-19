/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.

 * This Source Code may also be made available under the following Secondary
 * Licenses when the conditions for such availability set forth in the Eclipse
 * Public License v. 2.0 are satisfied: GPL-3.0 which is available at
 * https://www.gnu.org/licenses/gpl-3.0-standalone.html.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.util.persistence;

import org.eclipse.elk.core.data.LayoutMetaDataService;
import org.eclipse.elk.core.data.LayoutOptionData;
import org.eclipse.elk.graph.ElkGraphPackage;
import org.eclipse.elk.graph.impl.ElkPropertyToValueMapEntryImpl;
import org.eclipse.elk.graph.properties.IProperty;
import org.eclipse.emf.ecore.EObject;
import org.eclipse.emf.ecore.EStructuralFeature;
import org.eclipse.emf.ecore.xmi.XMLHelper;
import org.eclipse.emf.ecore.xmi.impl.XMISaveImpl;

/**
 * Specializes serialization such that properties are not serialized if we have no chance of parsing them in again
 * later. See {@link ElkGraphResource} for details.
 * 
 * @see ElkGraphResource
 */
final class ElkGraphXMISave extends XMISaveImpl {

    /**
     * Creates a new instance that uses the given helper.
     */
    ElkGraphXMISave(final XMLHelper helper) {
        super(helper);
    }
    
    
    @Override
    protected void saveElement(final EObject o, final EStructuralFeature f) {
        // We have special handling for property to value map entries
        if (f.equals(ElkGraphPackage.Literals.EMAP_PROPERTY_HOLDER__PROPERTIES)
                && o instanceof ElkPropertyToValueMapEntryImpl) {
            
            // Check if the property's data type is deserializable
            IProperty<?> property = ((ElkPropertyToValueMapEntryImpl) o).getKey();
            LayoutOptionData optionData = LayoutMetaDataService.getInstance().getOptionData(
                    property != null ? property.getId() : null);
            
            // Check if the value we want to serialize here will be parsable later so we can deserialize it
            // again properly; if this is not the case, there is no point serializing it in the first place
            if (optionData == null || (optionData != null && !optionData.canParseValue())) {
                return;
            }
        }
        
        super.saveElement(o, f);
    }

}
