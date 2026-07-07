/*******************************************************************************
 * Copyright (c) 2017 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.util.persistence;

import org.eclipse.elk.core.data.LayoutMetaDataService;
import org.eclipse.elk.core.util.IDataObject;
import org.eclipse.elk.core.util.internal.LayoutOptionProxy;
import org.eclipse.elk.graph.ElkGraphPackage;
import org.eclipse.elk.graph.properties.IProperty;
import org.eclipse.elk.graph.properties.Property;
import org.eclipse.emf.ecore.EDataType;
import org.eclipse.emf.ecore.EFactory;
import org.eclipse.emf.ecore.xmi.XMLResource;
import org.eclipse.emf.ecore.xmi.impl.XMIHelperImpl;

/**
 * An XMI helper which has special knowledge on how to load {@link IProperty} and {@link IDataObject} instances
 * while loading ELK graphs. A fresh instance must be used for every loading operation. See {@link ElkGraphResource}
 * for details.
 * 
 * @see ElkGraphResource
 */
final class ElkGraphXMIHelper extends XMIHelperImpl {
    
    /**
     * Creates a new XMI helper for the given resource.
     * 
     * @param resource the resource to help load or save.
     */
    ElkGraphXMIHelper(final XMLResource resource) {
        super(resource);
    }
    
    
    @Override
    protected Object createFromString(final EFactory eFactory, final EDataType eDataType, final String value) {
        // We know how to obtain IProperty and parse option values
        if (eDataType.equals(ElkGraphPackage.Literals.IPROPERTY)) {
            return createIPropertyFromString(value);
        } else if (eDataType.equals(ElkGraphPackage.Literals.PROPERTY_VALUE)) {
            return createPropertyValueFromString(value);
        } else {
            return super.createFromString(eFactory, eDataType, value);
        }
    }
    
    /**
     * Asks the layout meta data service to load the property which has the given ID.
     */
    private IProperty<?> createIPropertyFromString(final String id) {
        IProperty<?> property = LayoutMetaDataService.getInstance().getOptionData(id);
        if (property == null) {
            // ELK doesn't know about the property (yet?)
            property = new Property<Object>(id);
        }
        
        return property;
    }
    
    /**
     * Creates a {@link LayoutOptionProxy} for the given property value.
     */
    private Object createPropertyValueFromString(final String value) {
        return new LayoutOptionProxy(value);
    }
    
}
